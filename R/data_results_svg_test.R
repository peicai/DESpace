#' Results from \code{\link{svg_test}} function
#' 
#' @rdname results_svg_test
#' @name results_svg_test
#' @aliases results_svg_test
#' 
#' @param results_svg_test contains a \code{\linkS4class{list}} object,
#' with the results obtained applying \code{\link{svg_test}} function to an external dataset from the spatialLIBD package.
#' Below the code used to obtain 'results_svg_test'.
#
#' @examples
#' # load the input data:
#' # data("LIBD_subset", package = "DESpace")
#' # LIBD_subset
#' # 
#' # Fit the model via `svg_test` function. 
#' # Parameter `spe` specifies the input `SpatialExperiment` or `SingleCellExperiment` object, 
#' # while `cluster_col` defines the column names of `colData(spe)` containing spatial clusters. 
#' # To obtain all statistics, set `verbose` to `TRUE`.
#' # svg_test returns of a list of 4 objects:
#' # "gene_results": a dataframe contains main edgeR test results;
#' # "estimated_y": a DGEList object contains the estimated common dispersion, 
#' #  which can later be used to speed-up calculation when testing individual clusters.
#' # "glmFit": a DGEGLM object contains full statistics from "edgeR::glmFit"
#' # "glmLRT" a DGELRT object contains full statistics from "edgeR::glmLRT"
#' # 
#' # set.seed(123)
#' # results_svg_test <- svg_test(spe = LIBD_subset,
#' #                                          cluster_col = "layer_guess_reordered",
#' #                                          verbose = FALSE)
#' # 
#' # save(results_svg_test, file = "./DESpace/data/results_svg_test.RData")
#' 
#' @author Peiying Cai \email{peiying.cai@uzh.ch}, Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{svg_test}}
NULL