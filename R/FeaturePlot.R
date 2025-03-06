#' FeaturePlot
#' 
#' Plot spatial gene expression.
#' This function is a modified version of the \code{\link{FeaturePlot}} function from BayesSpace R package.
#' In comparison to the original BayesSpace function, this function allows plotting multiple genes simultaneously 
#' and drawing an outline around a specified cluster.
#' 
#' @param spe SpatialExperiment or SingleCellExperiment. If \code{feature} is specified and is a 
#'   string, it must exist as a row in the specified assay of \code{spe}.
#' @param feature Feature vector used to color each cell. May be the name of a
#'   gene/row in an assay of \code{spe}, or a vector of continuous values.
#' @param coordinates Column names for the spatial coordinates of cells stored in \code{colData(spe)}.
#'    If specified, these coordinates will be used. If not, the function defaults to using 'row' and 'col' 
#'    in \code{colData(spe)} if they exist. Otherwise, it will use \code{spatialCoords(spe)} if 'spe' is a 
#'    SpatialExperiment object and \code{spatialCoords(spe)} is not NULL.
#' @param concave_hull A logical value (TRUE or FALSE).
#'    If TRUE, the function uses `ggforce::geom_mark_hull()` to outline cluster boundaries 
#'    (recommended for non-discontinuous clusters).
#'    If FALSE, `sosta::reconstructShapeDensityImage()` is used for complex cluster shapes.
#'    For Visium or ST platforms, `concave_hull` is automatically set to TRUE.
#' @param sf_dim A numeric value for the x-dimension of the reconstruction (default is 200). 
#'    A lower value speeds up computation but reduces accuracy. 
#'    Used only when `concave_hull` is FALSE.
#' @param assay.type String indicating which assay in \code{spe} the expression
#'    vector should be taken from.
#' @param annotation_cluster A logical value (TRUE or FALSE). 
#'    If TRUE, annotated spatial clusters are plotted alongside expression plots. 
#'    If FALSE, clusters are not displayed.
#' @param annotation_title A character string for the title of the annotated spatial clusters. 
#'    Applied only when `annotation_cluster` is TRUE.
#' @param platform A character string specifying the spatial sequencing platform. 
#'   If "Visium" or "ST", a hexagonal spot layout will be used. 
#'   Otherwise, points will be plotted.
#' @param cluster_col Column name of spatial clusters in \code{colData(spe)}.
#' @param cluster Names of the spatial clusters used for drawing a boundary 
#'    around a group of points that belong to the specify cluster.
#'    It can be NULL, "all"/"ALL", or a vector of cluster names.
#' @param legend_cluster A logical value. TRUE of FALSE, 
#'    indicating whether to plot the legend for the shaped clusters (TRUE), or not (FALSE).
#'    Only used when 'cluster_col' and 'cluster' are specified.
#' @param legend_exprs A logical value. 
#'    TRUE of FALSE, indicating whether to plot the legend for 
#'    the expression level (TRUE), or not (FALSE).
#' @param diverging A logical value. 
#'    If TRUE, uses a diverging color gradient in \code{\link{FeaturePlot}} (e.g., for fold change). 
#'    If FALSE, uses a sequential gradient (e.g., for expression).
#' @param low,mid,high Optional hex codes for low, mid, and high values of the
#'   color gradient used for continuous cell values.
#' @param color Optional hex code to set color of borders around cells. Set to
#'   \code{NA} to remove borders.
#' @param linewidth The width of the boundary line around the cluster.
#' The default ('0.4') size of the boundary line is one.
#' @param linecolor The colors of the boundary lines around the cluster. 
#' If unspecified, the default color scheme is used.
#' @param label A logical. TRUE of FALSE. Adding a label and an arrow pointing to a group.
#' @param ncol The dimensions of the grid to create. 
#' By default, 1, if the length of feature equals to 1, and 3, otherwise.
#' @param title A logical. TRUE or FALSE. If true, the title name of each (subplot) is the gene name.
#' @param title_size Title font size.
#' @param point_size Point size.
#' @return Returns a ggplot object.
#' @examples
#' # load the input data:
#' data("LIBD_subset", package = "DESpace")
#' 
#' # load pre-computed results (obtained via `svg_test`)
#' data("results_svg_test", package = "DESpace")
#' 
#' # Visualize the gene expression of the top three genes
#' feature = results_svg_test$gene_results$gene_id[seq_len(3)]
#' FeaturePlot(LIBD_subset, feature, coordinates = c("array_row", "array_col"),
#'             ncol = 3, title = TRUE)
#' 
#' @family spatial plotting functions
#'
#' @importFrom ggplot2 ggplot aes geom_polygon scale_fill_gradient scale_fill_gradient2 coord_equal labs theme_void element_text geom_sf
#' @importFrom scales muted
#' @importFrom assertthat assert_that
#' 
#' @seealso \code{\link{svg_test}}, \code{\link{individual_svg}}, \code{\link{top_results}}, \code{\link{dsp_test}}, \code{\link{individual_dsp}}
#' 
#' @export
FeaturePlot <- function(spe, feature, coordinates = NULL,
                        concave_hull = FALSE, sf_dim = 200, 
                        assay.type="logcounts", annotation_cluster = FALSE,
                        annotation_title = NULL,
                        platform = "Visium", 
                        cluster_col = NULL, cluster = NULL, 
                        legend_cluster = FALSE, legend_exprs = FALSE,
                        diverging = FALSE,
                        low = NULL, high = NULL, mid = NULL,
                        color = NULL, linewidth = 0.4, linecolor = NULL,
                        label = FALSE, ncol = 3, title = FALSE,
                        title_size = 10, point_size = 0.5) {
    if (!is(spe, "SpatialExperiment") & !is(spe, "SingleCellExperiment")){
        message("'spe' must be a 'SpatialExperiment' object or a 'SingleCellExperiment' object")
        return(NULL)
    }
  
    if (is(spe, "SpatialExperiment")){
        coord <- is.null(spatialCoords(spe))
    }
  
    if ('row' %notin% colnames(colData(spe)) | 'col' %notin% colnames(colData(spe))){
        if (is.null(coordinates) & (coord | platform %in% c("Visium", "ST"))){
            message("Spatial coordinates of cells are missing in colData(spe) and spatrialCoords(spe). ")
            message("Please provide column names of spatial coordinates in 'coordinates'. ")
            message("Default setting are 'row' and 'col'. ")
            message("Alternatively, add spatial coordinates to spatrialCoords(spe) if 'spe' is a 'SpatialExperiment' object.")
        return(NULL)
    }else if (!is.null(coordinates)){
        coordinates_col <- colData(spe)[coordinates]
        colnames(coordinates_col) <- c('row', 'col')
        colData(spe) <- cbind(colData(spe), coordinates_col)
        spatialCoords(spe) <- as.matrix(coordinates_col)
    }else {
        coordinates_col <- spatialCoords(spe)
        colnames(coordinates_col) <- c('row', 'col')
        colData(spe) <- cbind(colData(spe), coordinates_col)
    }
    }
    if (assay.type %notin%  assayNames(spe)){
        message("assay.type is missing in assayNames(spe). ")
        return(NULL)
    }
    if (is.character(feature)) {
        assert_that(all(feature %in% rownames(spe)), msg="Feature not in spe.")
        fill <- assay(spe, assay.type)[feature, ]
        fill.name <-  feature
    } else {
        fill <- feature
        x <- fill
        ## this could be an argument, but it's easily overwritten with labs()
        ## and we should encourage composing ggplot functions instead
        fill.name <- "Expression"
    }
    ## No borders around subspots by default
    if (is.null(color)) {
        color <- NA ##"#d8dcd6"
    }
  
    if(any(cluster %in% c("all", "ALL"))){
        unique_values <- colData(spe)[[cluster_col]]
        n_cluster <- length(unique(unique_values))
    }else{
        n_cluster <- length(cluster)
    }
  
    if(is.null(linecolor)){
        linecolor <- c('#000000', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                '#911eb4', '#46f0f0', '#800000', '#fabebe', '#008080', 
                '#e6beff', '#000075', '#fffac8', '#aaffc3', '#808000',
                '#ffd8b1', '#808080', '#ffffff' 
        )[seq_len(n_cluster)]
    }else{
        assert_that(
        length(linecolor) == n_cluster, 
        msg = "The number of outline colors provided must match the length of 'cluster'."
    ) }
  
    MoreArgs <- list(spe, diverging, low, high, mid, legend_exprs, color,
                platform, cluster_col, cluster,
                label, title, title_size, linewidth, linecolor,
                sf_dim, concave_hull, point_size)
    ## if feature is a vector of gene names
    if(length(feature) > 1){
        x <- asplit(as.data.frame(as.matrix(fill)), 1)
        plot.list <- mapply(.geneExprsPlot, x, fill.name, MoreArgs = MoreArgs, SIMPLIFY = FALSE)
        expression.plots <- lapply(plot.list, `[[`, 1) 
        vertices <- lapply(plot.list, `[[`, 2)[[1]]
        widths <- NULL
    }else{
        x <- fill
        plot.list <- .geneExprsPlot(x, fill.name, spe, diverging, low, high, mid, legend_exprs, color,
                        platform, cluster_col, cluster,
                        label, title, title_size, linewidth, linecolor,
                        sf_dim, concave_hull, point_size)
        expression.plots <- list(plot.list$plot)
        vertices <- plot.list$vertices
        ncol <- 1; widths <- c(5, 1)
    }
    ## if 'cluster_col' and 'cluster' are specified, and legend_cluster == TRUE,  
    ## draw the shape of the outline of the group and add the legend of clusters
    if(!is.null(cluster_col) && !is.null(cluster) && legend_cluster){
        ## create the legend of cluster
        if(!any(cluster %in% c("all", "ALL"))){
            filter <- vertices$Cluster %in% (cluster)
        }else{
            filter <- NULL
        }
    my_hist <- vertices %>% 
        ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
        geom_mark_hull(aes(x = x.vertex, y = y.vertex, 
                        color = Cluster, 
                        fill = Cluster, 
                        filter = filter),
                        expand = unit(0.02, "cm"), # Adjusts the spacing
                        inherit.aes = FALSE) +
        scale_color_manual(values = linecolor) +  # Map custom border colors
        scale_fill_manual(values = linecolor) +   # Map fill colors if necessary
        theme(legend.position = "right")
    legend <- get_plot_component(my_hist, 'guide-box-top', return_all = TRUE)
    plots <- append(expression.plots, list(as_ggplot(legend)) )
    plots <- wrap_plots(plots, ncol=min(length(expression.plots),ncol)+1,widths = widths)
    }else{
        plots <- wrap_plots(expression.plots, ncol=min(length(expression.plots),ncol))
    }
    ## if 'annotation_cluster' is specified as TRUE, plot the original spatial clusters next to expression plots
    if(annotation_cluster){
        if (is.null(cluster_col))
            stop("Column names of spatial clusters are not specified.")
        df1 <- data.frame(cluster_col = colData(spe)[[cluster_col]],
                spot = rownames(colData(spe)))
        df1 <- merge(df1, vertices, by = "spot", all.x = TRUE)
        p1 <- df1 %>% 
            ggplot(mapping = aes(x=x.vertex, y=y.vertex,
                           color=factor(cluster_col))) +
            geom_point(size = point_size) + theme_void() +  coord_equal() +
            theme(legend.position="right") + labs(color = "", title = "") +
            guides(color = guide_legend(override.aes = list(size = 2))) +
            ggtitle(annotation_title)
        plots <- append(list(p1), list(plots))
        plots <- wrap_plots(plots, ncol=2, 
                widths = c(1,min(length(expression.plots),ncol)))
    }
    plots
}
