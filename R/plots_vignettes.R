# Make gene expression plots
#
# @param x fill; A vector of values to use as fill for each spot.
# @param y fill.name; expression legend name.
# @param spe SpatialExperiment or SingleCellExperiment with row/col in colData.
# @param platform "Visium" or "ST", used to determine spot layout.
# @param is.enhanced If true, \code{spe} contains enhanced subspot data instead
#   of spot-level expression. Used to determine spot layout.
# @param spatial_cluster Column name of clusters in \code{colData(spe)}.
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
                            low, high, mid,legend_exprs,
                            color,
                            platform, is.enhanced=FALSE,
                            spatial_cluster, cluster, 
                            label, title,title_size=10,linewidth,
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
    vertices <- .make_vertices(spe, fill, platform, is.enhanced)
    if (diverging) {
        low <- ifelse(is.null(low), "#F0F0F0", low)
        mid <- NULL
        high <- ifelse(is.null(high), "#832424", high)
    } else {
        low <- ifelse(is.null(low), "#3A3A98", low)
        mid <- ifelse(is.null(mid), "#F0F0F0", mid)
        high <- ifelse(is.null(high), "#832424", high)
    }
    ## if 'spatial_cluster' and 'cluster' are specified, draw the shape of the outline of the group
    if(!is.null(spatial_cluster) && !is.null(cluster)){
        cdata <- data.frame(colData(spe))
        Cluster <- as.character(cdata[[spatial_cluster]])
    if(cluster %notin% c("all", "ALL")){
        Cluster[Cluster %notin% (cluster)] <- 'Others'}
    Cluster <- as.factor(Cluster)
    vertices <- cbind(vertices, Cluster)
    vertices <- vertices %>% filter(!is.na(Cluster))
    ## annotate clusters
    if(label){
        label <- Cluster
    }else{
        label <- NULL
    }

    if(cluster %notin% c("all", "ALL")){
        ## Add filter
        splot <- vertices %>% 
            ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
            geom_polygon(aes(group=spot, fill=fill), color=color) +
            labs(fill=fill.name) + coord_equal() +
            theme_void() + scale_fill_gradient2(low=low, mid=mid, high=high)+ 
            new_scale_fill() + new_scale_color()  + 
            suppressWarnings(geom_mark_hull( aes(x=x.vertex, y=y.vertex, 
                                color = Cluster, fill=Cluster,linewidth = I(linewidth),
                                label = label, filter = Cluster %in% (cluster)),
                                alpha=0, expand = unit(0.1, "mm"), 
                                radius = unit(0.4, "mm"),
                                show.legend = FALSE) )
    }else{
        (splot <- vertices %>% 
            ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
            geom_polygon(aes(group=spot, fill=fill), color=color) +
            labs(fill=fill.name) + coord_equal() +
            theme_void() + scale_fill_gradient2(low=low, mid=mid, high=high)+ 
            new_scale_fill() + new_scale_color()  + 
            suppressWarnings(geom_mark_hull( aes(x=x.vertex, y=y.vertex,
                                    color = Cluster, fill=Cluster,linewidth = I(linewidth),
                                    label = label),
                                    alpha=0, expand = unit(0.05, "mm"),
                                    radius = unit(0.2, "mm"),
                                    show.legend = FALSE) ))
        
    }
    splot <- splot + scale_alpha(guide="none") + theme_void() + 
        scale_fill_discrete() +
        scale_color_discrete() + 
        theme(legend.position = "bottom",
                legend.key.size = unit(0.5, 'cm'))
    
    }else{
        splot <-  vertices %>% 
        ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
        geom_polygon(aes(group=spot, fill=fill), color=color) +
        labs(fill=fill.name) +
        coord_equal() +
        theme_void() + scale_fill_gradient2(low=low, mid=mid, high=high) +
        theme(legend.position = "bottom",
            legend.key.size = unit(2, 'cm'))
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
    return(list(plot = splot, vertices = vertices))
}
# Make vertices outlining spots/subspots for geom_polygon()
# 
# @param spe SpatialExperiment or SingleCellExperiment with row/col in colData
# @param fill Name of a column in \code{colData(spe)} or a vector of values to
#   use as fill for each spot
# @param platform "Visium" or "ST", used to determine spot layout
# @param is.enhanced If true, \code{spe} contains enhanced subspot data instead
#   of spot-level expression. Used to determine spot layout.
#   
# @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#   vertices outlining the spot's border
# 
# keywords internal
.make_vertices <- function(spe, fill, platform, is.enhanced) {
    cdata <- data.frame(colData(spe))
    if (platform == "Visium") {
        if (is.enhanced) {
            vertices <- .make_triangle_subspots(cdata, fill)
        } else {
            vertices <- .make_hex_spots(cdata, fill)
        }
    } else if (platform == "ST") {
        if (is.enhanced) {
            vertices <- .make_square_spots(cdata, fill, scale.factor=(1/3))
        } else {
            vertices <- .make_square_spots(cdata, fill)
        }
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
# Make vertices for each triangle subspot of a hex
# 
# @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#   vertices outlining the spot's border
#
# keywords internal
.make_triangle_subspots <- function(cdata, fill="spatial.cluster") {
    spot_positions <- .select_subspot_positions(cdata, x="spot.col", y="spot.row", fill=fill)
    spot_positions <- .adjust_hex_centers(spot_positions)
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r
    ## Make lists of triangle vertices (with respect to hex center)
    ## subspot.idx is same ordering as `shift` in spatialEnhance
    ## that is, beginning in top right and proceeding clockwise, (1, 5, 3, 4, 6, 2)
    ## NOTE: however, we need to reflect over x-axis to match global negation of y-coordinate
    vertex_offsets <- do.call(rbind, list(
        data.frame(x.offset=c(0, 0, r), y.offset=c(0, -R, -R/2), subspot.idx=3),
        data.frame(x.offset=c(0, r, r), y.offset=c(0, -R/2, R/2), subspot.idx=5),
        data.frame(x.offset=c(0, r, 0), y.offset=c(0, R/2, R), subspot.idx=1),
        data.frame(x.offset=c(0, 0, -r), y.offset=c(0, R, R/2), subspot.idx=2),
        data.frame(x.offset=c(0, -r, -r), y.offset=c(0, R/2, -R/2), subspot.idx=6),
        data.frame(x.offset=c(0, -r, 0), y.offset=c(0, -R/2, -R), subspot.idx=4)
    ))
    ## note that instead of cartesian product, `merge()` does an outer join
    ## on subspot.idx here
    spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
    spot_vertices$y.vertex <- -spot_vertices$y.vertex
    spot_vertices
}