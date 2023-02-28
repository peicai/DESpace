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
#' @param legend_layer A logical. TRUE of FALSE, indicating whether to plot the legend for the shaped layers (TRUE), or not (FALSE).
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
#' # Connect to ExperimentHub
#' ehub <- ExperimentHub::ExperimentHub()
#' # Download the example spe data
#' spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)
#' spe_all
#' 
#' # Only use one sample:
#' library(SpatialExperiment)
#' spe3 <- spe_all[, colData(spe_all)$sample_id == '151673']
#' rm(spe_all)
#' # Select small set of random genes for faster runtime in this example
#' set.seed(123)
#' sel_genes <- sample(dim(spe3)[1],500)
#' spe3 <- spe3[sel_genes,]
#' 
#' # load pre-computed results (obtained via `DESpace_test`)
#' data("results_DESpace_test", package = "DESpace")
#' 
#' # DESpace_test returns of a list of 2 objects:
#' # "gene_results": a dataframe contains main edgeR test results;
#' # "estimated_y": a DGEList object contains the estimated common dispersion, 
#' #  which can later be used to speed-up calculation when testing individual clusters.
#' 
#' # We visualize differential results:
#' head(results_DESpace_test$gene_results, 3)
#' 
#' # Visualize the gene expression of the top three genes
#' feature = results_DESpace_test$gene_results$gene_id[seq_len(3)]
#' feature
#' FeaturePlot(spe3, feature, ncol = 3, title = TRUE)
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
                        spatial_cluster = NULL, cluster = NULL,legend_layer = FALSE,
                        label = FALSE, ncol = 3, title = FALSE,linewidth = 0.4,
                        legend_exprs = FALSE,title_size = 5) {
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

  #if (is.null(platform))
  #platform <- .bsData(spe, "platform", "Visium")
  #if (is.null(is.enhanced))
  #is.enhanced <- .bsData(spe, "is.enhanced", FALSE)
  
  ## extract expression from logcounts if a gene name is passed.
  ## otherwise, assume a vector of counts was passed and let
  ## .make_vertices helpers check validity
  if (is.character(feature)) {
    assert_that(all(feature %in% rownames(spe)),
                #msg=sprintf("Feature %s not in spe.", feature))
                msg="Feature not in spe.")
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

  MoreArgs = list(spe, diverging, low, high, mid, legend_exprs, color,
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
    ncol=1
    widths <- c(5, 1)
  }
  
  
  ## if 'spatial_cluster' and 'cluster' are specified, and legend_layer == TRUE,  draw the shape of the outline of the group and add the legend of layers
  if(!is.null(spatial_cluster) && !is.null(cluster) && legend_layer){
    ## create the legend of layer
    if(cluster %notin% c("all", "ALL")){
      filter <- vertices$Layer %in% (cluster)
    }else{
      filter <- NULL
    }
    
    my_hist <- vertices %>% 
      ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
      ggforce::geom_mark_hull( aes(x=x.vertex, y=y.vertex, 
                                   color = Layer, fill=Layer, filter = filter))  +
      theme(legend.position = "right") 
    
    legend <- cowplot::get_legend(my_hist)
    plots <- append(expression.plots, list(ggpubr::as_ggplot(legend)) )
    plots <- patchwork::wrap_plots(plots, ncol=min(length(expression.plots),ncol)+1,widths = widths) #& 
      #theme(legend.position = "bottom") 
  }else{
    plots <- patchwork::wrap_plots(expression.plots, ncol=min(length(expression.plots),ncol))

  }
  
  ## if 'Annotated_cluster' is specified as TRUE, plot the original spatial clusters next to expression plots
  if(Annotated_cluster){
    if (is.null(spatial_cluster))
      stop("Column names of spatial clusters not specified.")
      # Manual annotation
      # p1 <- ggplot(as.data.frame(colData(spe)), aes(x=col, y=row, color=factor(get(spatial_cluster)))) + 
      #   geom_point() + theme_void() + scale_y_reverse() + coord_equal() +
      #   scale_x_reverse() + 
      #   theme(legend.position="none") + labs(color = "", title = "")
      
      df1 = data.frame(spatial_cluster = colData(spe)[[spatial_cluster]],
                       spot = rownames(colData(spe)))
      df1 = merge(df1, vertices, by = "spot", all.x = TRUE)
      
      p1 <- df1 %>% 
        ggplot(mapping = aes(x=x.vertex, y=y.vertex,
                             color=factor(spatial_cluster))) +
        #ggplot(as.data.frame(colData(spe)), aes(x=col, y=row, color=factor(get(spatial_cluster)))) + 
        geom_point() + theme_void() +  coord_equal() +
        #scale_x_reverse() + scale_y_reverse() +
        theme(legend.position="none") + labs(color = "", title = 
                                               #paste0("Sample: ")
                                             "")
      plots <- append(list(p1), list(plots))
      plots <- patchwork::wrap_plots(plots, ncol=2, 
                                         widths = c(1,min(length(expression.plots),ncol))
                                         )
      }
  
  plots
  }
