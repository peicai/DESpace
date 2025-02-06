# Make gene expression plots
#
# @param x fill; A vector of values to use as fill for each spot.
# @param y fill.name; expression legend name.
# @param spe SpatialExperiment or SingleCellExperiment with row/col in colData.
# @param platform "Visium" or "ST", used to determine spot layout.
# @param cluster_col Column name of clusters in \code{colData(spe)}.
# @param cluster Cluster names used for drawing a boundary around a group of points (belong to the specify cluster) to drive attention.
# Can be NULL, "all"/"ALL", and a vector of cluster names.
# @param label TRUE of FALSE. Adding a label and an arrow pointing to a group.
# @param title TRUE or FALSE. If true, the title name of each (subplot) is the gene name.
#
# @return Returns a list containing a (list of) ggplot object and a vertices (table of (x.pos, y.pos, spot, fill, Cluster)).
#
# keywords internal
.geneExprsPlot <- function( x, y, spe, 
                            diverging,
                            low, high, mid, legend_exprs,
                            color,
                            platform, 
                            cluster_col, cluster, 
                            label, title, title_size=10,
                            linewidth = 0.75, linecolor = NULL,
                            sf_dim = 200, concave_hull = TRUE,
                            point_size = 0.5,
                            ...){
    fill <- as.vector(x)
    fill.name <- y
    ## fill.name -> expression legend name
    ## title.name -> subtitle name
    ## if fill.name == 'Expression', don't provide title (title <- FALSE)
    ## if fill.name is gene.name: 
    #### if title == T, title name -> gene name, and the legend name (fill.name) is NULL
    #### if title != T, title name -> NULL, and fill.name is gene name
    if(fill.name == 'Expression'){
        title <-  FALSE
    } 
    if(title){
        title.name <-  fill.name
        fill.name <- NULL
    }else{
        title.name <- NULL
    }
    is_spot <- platform %in% c("Visium", "ST")
    if(is_spot){
        vertices <- .make_vertices(spe, fill, platform)
        vertices_spot <- vertices[["spot"]]
        concave_hull <- TRUE
    }else{
        vertices <- colData(spe)[,c('row', 'col')] %>% data.frame()
        colnames(vertices) <- c("x.vertex", "y.vertex")
        vertices[["fill"]] <- fill
        vertices_spot <- rownames(vertices)
    }
    if (diverging) {
        low <- ifelse(is.null(low), "#F0F0F0", low)
        mid <- NULL
        high <- ifelse(is.null(high), "#832424", high)
    } else {
        low <- ifelse(is.null(low), "#3A3A98", low)
        mid <- ifelse(is.null(mid), "#F0F0F0", mid)
        high <- ifelse(is.null(high), "#832424", high)
    }
    
    # Create the plot based on the platform type
    # Define the base plot
    splot <- vertices %>%
      ggplot(aes(x = x.vertex, y = y.vertex)) + 
      coord_equal() + 
      theme_void() + 
      theme(legend.position = "bottom", 
            legend.key.size = unit(2, 'cm'))
    
    # Add the conditional layers based on `is_spot`
    if (is_spot) {
      splot <- splot + 
        geom_polygon(aes(group = .data[["spot"]], 
                         fill = .data[["fill"]]), 
                     color = color) +
        scale_fill_gradient2(low = low, mid = mid, high = high) +
        labs(fill = fill.name) 
    } else {
      threshold <- quantile(vertices[["fill"]], 0.9999)
      splot <- splot + 
        geom_point(size = point_size, aes(color = pmin(fill, threshold))) +
        scale_colour_gradientn(colors = c(low, mid, high), breaks = c(0, threshold), labels = c("low", "high")) +
        labs(color = fill.name)
    }
    
    ## if 'cluster_col' and 'cluster' are specified, draw the shape of the outline of the group
    if(!is.null(cluster_col) && !is.null(cluster)){
        cdata <- data.frame(colData(spe))
        Cluster <- as.character(cdata[[cluster_col]])
        if(!any(cluster %in% c("all", "ALL"))){
            Cluster[Cluster %notin% (cluster)] <- 'Others'
            cluster_use <- cluster
        }else{
            cluster_use <- unique(Cluster)
        }
        Cluster <- as.factor(Cluster)
        vertices <- cbind(vertices, Cluster)
        #vertices <- vertices %>% filter(!is.na(Cluster))
        ## annotate clusters
        if(label){
            label <- Cluster
        }else{
            label <- NULL
        }
    
        if(concave_hull){
          ## Add filter
          Cluster_use <- vertices[["Cluster"]]
          splot <- splot +  
            new_scale_fill() + new_scale_color()  + 
            suppressWarnings(geom_mark_hull( aes(x=x.vertex, 
                                                 y=y.vertex,
                                                 color = Cluster_use, fill = Cluster_use, 
                                                 linewidth = I(linewidth),
                                                 label = label, filter = Cluster_use %in% (cluster_use)),
                                             alpha=0, expand = unit(0.1, "mm"),
                                             radius = unit(0.4, "mm"),
                                             show.legend = FALSE) ) +
            scale_color_manual(values = linecolor)
        }else{
          # Add layers using lapply
          geoms <- lapply(seq_along(cluster_use), function(clus) {
            sf_poly <- reconstructShapeDensityImage(
              spe,
              marks = cluster_col,
              mark_select = cluster_use[clus],
              dim = sf_dim,
              image_id = NULL,
              image_col = NULL
            )
            
            geom_sf(
              data = sf_poly,
              fill = NA,
              color = linecolor[clus],
              inherit.aes = FALSE,
              linewidth = linewidth
            )
          })
          # Combine the layers
          splot <- splot + geoms
        }
            
    splot <- splot + scale_alpha(guide="none") + theme_void() + 
      theme(legend.position = "bottom",
            legend.key.size = unit(0.5, 'cm'))
    } 
    
    if(is.na(legend_exprs)){
        legend_exprs <-  FALSE
    }
    if(legend_exprs == FALSE){
        splot <- splot + theme(legend.position = "none")
    }
    if(title){
        splot <- splot + labs(title=title.name)+ theme(plot.title = element_text(size=title_size))
    }
    vertices[["spot"]] <- vertices_spot
    return(list(plot = splot, vertices = vertices))
}
# Make vertices outlining spots/subspots for geom_polygon()
# 
# @param spe SpatialExperiment or SingleCellExperiment with row/col in colData
# @param fill Name of a column in \code{colData(spe)} or a vector of values to
#   use as fill for each spot
# @param platform "Visium" or "ST", used to determine spot layout
  
# @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#   vertices outlining the spot's border
# 
# keywords internal
.make_vertices <- function(spe, fill, platform) {
    cdata <- data.frame(colData(spe))
    if (platform == "Visium") {
        vertices <- .make_hex_spots(cdata, fill)
    } else if (platform == "ST") {
        vertices <- .make_square_spots(cdata, fill)
    } else {
        stop("Unsupported platform: \"", platform, "\". Cannot create spot layout.")
    }
    vertices
}
# Helper to extract x, y, fill ID from colData
# 
# @return Dataframe of (x.pos, y.pos, fill) for each spot
# 
# keywords internal
# @importFrom assertthat assert_that
.select_spot_positions <- function(cdata, x="col", y="row", fill="spatial.cluster") {
    ## Provide either a column name or vector of labels/values
    assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
    ## I think this is the best way to check if something is a string
    if (is.character(fill) && length(fill) == 1) {
        spot_positions <- cdata[, c(x, y, fill)]
        colnames(spot_positions) <- c("x.pos", "y.pos", "fill")    
    } else if (is.vector(fill) || is.factor(fill)) {
        assert_that(nrow(cdata) == length(fill))
        spot_positions <- cdata[, c(x, y)]
        colnames(spot_positions) <- c("x.pos", "y.pos")    
        spot_positions$fill <- fill
    }
    spot_positions$spot <- rownames(spot_positions)
    spot_positions
}
# Compute vertex coordinates for each spot in frame of plot
#
# @param spot_positions Center for hex, top left for square
# @param vertex_offsets Data frame of (x, y) offsets wrt spot position for each
#   vertex of spot
# 
# @return Cartesian product of positions and offsets, with coordinates
#   computed as (pos + offset)
#
# keywords internal
.make_spot_vertices <- function(spot_positions, vertex_offsets) {
    spot_vertices <- merge(spot_positions, vertex_offsets)
    spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
    spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
    as.data.frame(spot_vertices)
}
# Make vertices for each hex spot
# 
# @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#   vertices outlining the spot's border
# 
# keywords internal
.make_hex_spots <- function(cdata, fill) {
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r
    spot_positions <- .select_spot_positions(cdata, fill=fill)
    spot_positions <- .adjust_hex_centers(spot_positions)
    ## vertices of each hex (with respect to center coordinates)
    ## start at top center, loop clockwise
    vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
                                y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))
    spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
    ## Flip to match image orientation
    spot_vertices$y.vertex <- -spot_vertices$y.vertex
    spot_vertices
}
# Adjust hex spot positions so hexagons are adjacent to each other in plot
#
# Spots are regular hexagons with one unit of horizontal distance
# between centers
# 
# @return Shifted spot centers
# 
# keywords internal
.adjust_hex_centers <- function(spot_positions) {
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r
    ## Start at (1-indexed origin)
    spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
    spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1
    ## Shift centers up so rows are adjacent
    spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
    ## Spot columns are offset by row
    ## (i.e. odd rows have odd numbered columns, even rows have even)
    ## Shift centers to the left so columns are adjacent (but hexes stay offset)
    spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2
    spot_positions
}
# Make vertices for each square spot
#
# Squares are simple, just mae a unit square at each array coordinate
# 
# @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#   vertices outlining the spot's border
# 
# keywords internal
.make_square_spots <- function(cdata, fill="spatial.cluster", scale.factor=1) {
    spot_positions <- .select_spot_positions(cdata, fill=fill)
    vertex_offsets <- data.frame(x.offset=c(0, 1, 1, 0),
                                y.offset=c(0, 0, 1, 1))
    vertex_offsets <- vertex_offsets * scale.factor
    .make_spot_vertices(spot_positions, vertex_offsets)
}
# Helper to pull out subspot position columns
# Probably redundant with select_spot_positions above, but we need subspot.idx
# 
# @return Dataframe of (x.pos, y.pos, fill) for each spot
# 
# keywords internal
.select_subspot_positions <- function(cdata, x="spot.col", y="spot.row", fill="spatial.cluster") {
    ## Provide either a column name or vector of labels/values
    assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
    if (is.character(fill) && length(fill) == 1) {
        spot_positions <- cdata[, c(x, y, "subspot.idx", fill)]
        colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx", "fill")
    } else if (is.vector(fill) || is.factor(fill)) {
        assert_that(nrow(cdata) == length(fill))
        spot_positions <- cdata[, c(x, y, "subspot.idx")]
        colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx")    
        spot_positions$fill <- fill
    }
    spot_positions$spot <- rownames(spot_positions)
    spot_positions
}