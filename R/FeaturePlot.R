#' FeaturePlot
#' 
#' Plot spatial gene expression.
#' This function is a modified version of the \code{\link{FeaturePlot}} function from BayesSpace R package.
#' In comparison to the original BayesSpace function, this function allows plotting multiple genes simultaneously and drawing an outline around a specified cluster.
#' @param spe SpatialExperiment or SingleCellExperiment. If \code{feature} is specified and is a 
#'   string, it must exist as a row in the specified assay of \code{spe}.
#' @param feature Feature vector used to color each spot. May be the name of a
#'   gene/row in an assay of \code{spe}, or a vector of continuous values.
#' @param coordinates Column names of spatial coordinates of spots stored in \code{colData(spe)}.
#' @param assay.type String indicating which assay in \code{spe} the expression
#'   vector should be taken from.
#' @param low,mid,high Optional hex codes for low, mid, and high values of the
#'   color gradient used for continuous spot values.
#' @param diverging A logical. If true, use a diverging color gradient in
#'   \code{\link{FeaturePlot}} (e.g. when plotting a fold change) instead of a
#'   sequential gradient (e.g. when plotting expression).
#' @param spatial_cluster Column name of spatial clusters in \code{colData(spe)}.
#' @param cluster Names of the spatial clusters used for drawing a boundary 
#' around a group of points that belong to the specify cluster.
#' It can be NULL, "all"/"ALL", or a vector of cluster names.
#' @param Annotated_cluster A logical. TRUE or FALSE, indicating whether to plot the annotated spatial clusters next to expression plots.
#' @param linewidth The width of the boundary line around the cluster.
#' The default ('0.4') size of the boundary line is one.
#' @param legend_cluster A logical. TRUE of FALSE, indicating whether to plot the legend for the shaped clusters (TRUE), or not (FALSE).
#' Only used when 'spatial_cluster' and 'cluster' are specified.
#' @param legend_exprs A logical. TRUE of FALSE, indicating whether to plot the legend for the expression level (TRUE), or not (FALSE).
#' @param label A logical. TRUE of FALSE. Adding a label and an arrow pointing to a group.
#' @param platform Spatial sequencing platform. If "Visium", the hex spot layout
#'   will be used, otherwise square spots will be plotted.\cr
#'   NOTE: specifying this argument is only necessary if \code{spe} was not
#'   created by \code{BayesSpace::patialCluster()} or \code{BayesSpace::spatialEnhance()}.
#' @param ncol The dimensions of the grid to create. 
#' By default, 1, if the length of feature equals to 1, and 3, otherwise.
#' @param title A logical. TRUE or FALSE. If true, the title name of each (subplot) is the gene name.
#' @param color Optional hex code to set color of borders around spots. Set to
#'   \code{NA} to remove borders.
#' @param title_size Text size.
#' @return Returns a ggplot object.
#' @examples
#' # load the input data:
#' data("LIBD_subset", package = "DESpace")
#' 
#' # load pre-computed results (obtained via `DESpace_test`)
#' data("results_DESpace_test", package = "DESpace")
#' 
#' # Visualize the gene expression of the top three genes
#' feature = results_DESpace_test$gene_results$gene_id[seq_len(3)]
#' FeaturePlot(LIBD_subset, feature, coordinates = c("array_row", "array_col"),
#'             ncol = 3, title = TRUE)
#' 
#' @family spatial plotting functions
#'
#' @importFrom ggplot2 ggplot aes geom_polygon scale_fill_gradient scale_fill_gradient2 coord_equal labs theme_void element_text
#' @importFrom scales muted
#' @importFrom assertthat assert_that
#' 
#' @seealso \code{\link{DESpace_test}}, \code{\link{individual_test}}, \code{\link{top_results}}
#' 
#' @export
FeaturePlot <- function(spe, feature, coordinates = NULL,
                        assay.type="logcounts", Annotated_cluster = FALSE,
                        diverging=FALSE,
                        low=NULL, high=NULL, mid=NULL,
                        color=NULL,
                        platform= "Visium", #is.enhanced=FALSE,
                        spatial_cluster = NULL, cluster = NULL,legend_cluster = FALSE,
                        label = FALSE, ncol = 3, title = FALSE,linewidth = 0.4,
                        legend_exprs = FALSE,title_size = 10) {
    if (!is(spe, "SpatialExperiment") & !is(spe, "SingleCellExperiment")){
        message("'spe' must be a 'SpatialExperiment' object or a 'SingleCellExperiment' object")
    return(NULL)
    }
    if('row' %notin% colnames(colData(spe)) | 'col' %notin% colnames(colData(spe))){
        if(is.null(coordinates)){
            message("Spatial coordinates of spots are missing in colData(spe). ")
            message("Please provide column names of spatial coordinates in 'coordinates'. ")
            message("Default setting are 'row' and 'col'. ")
        return(NULL)
    }else{
        coordinates_col <- colData(spe)[coordinates]
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
    MoreArgs <- list(spe, diverging, low, high, mid, legend_exprs, color,
                    platform, is.enhanced = FALSE, spatial_cluster, cluster,
                    label, title, title_size, linewidth)
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
                                platform, is.enhanced = FALSE, spatial_cluster, cluster, label, title,
                                title_size, linewidth)
        expression.plots <- list(plot.list$plot)
        vertices <- plot.list$vertices
        ncol <- 1; widths <- c(5, 1)
    }
    ## if 'spatial_cluster' and 'cluster' are specified, and legend_cluster == TRUE,  
    ## draw the shape of the outline of the group and add the legend of clusters
    if(!is.null(spatial_cluster) && !is.null(cluster) && legend_cluster){
        ## create the legend of cluster
        if(cluster %notin% c("all", "ALL")){
            filter <- vertices$Cluster %in% (cluster)
        }else{
            filter <- NULL
            }
        my_hist <- vertices %>% 
        ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
        geom_mark_hull(aes(x=x.vertex, y=y.vertex, 
                                color = Cluster, fill=Cluster, filter = filter))  +
        theme(legend.position = "right")
        legend <- get_legend(my_hist)
        plots <- append(expression.plots, list(as_ggplot(legend)) )
        plots <- wrap_plots(plots, ncol=min(length(expression.plots),ncol)+1,widths = widths)
        }else{
            plots <- wrap_plots(expression.plots, ncol=min(length(expression.plots),ncol))
        }
    ## if 'Annotated_cluster' is specified as TRUE, plot the original spatial clusters next to expression plots
    if(Annotated_cluster){
        if (is.null(spatial_cluster))
            stop("Column names of spatial clusters are not specified.")
        df1 <- data.frame(spatial_cluster = colData(spe)[[spatial_cluster]],
                        spot = rownames(colData(spe)))
        df1 <- merge(df1, vertices, by = "spot", all.x = TRUE)
        p1 <- df1 %>% 
            ggplot(mapping = aes(x=x.vertex, y=y.vertex,
                    color=factor(spatial_cluster))) +
            geom_point() + theme_void() +  coord_equal() +
            theme(legend.position="none") + labs(color = "", title = "")
        plots <- append(list(p1), list(plots))
        plots <- wrap_plots(plots, ncol=2, 
                                        widths = c(1,min(length(expression.plots),ncol)))
        }
    plots
}
