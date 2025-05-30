#' Results from \code{\link{individual_svg}} function
#' 
#' @rdname results_individual_svg
#' @name results_individual_svg
#' @aliases results_individual_svg
#' 
#' @param results_individual_svg contains a \code{\linkS4class{list}} object,
#' with the results obtained applying \code{\link{individual_svg}} function to an external dataset from the spatialLIBD package.
#' Below the code used to obtain 'results_individual_svg'.
#
#' @examples
#' # load the input data:
#' # data("LIBD_subset", package = "DESpace")
#' # LIBD_subset
#' # load pre-computed results (obtained via `svg_test`)
#' # data("results_svg_test", package = "DESpace")
#' # results_svg_test
#' 
#' # Function `individual_svg()` can be used to identify SVGs for each individual cluster.
#' # Parameter `spatial_cluster` indicates the column names of `colData(spe)` 
#' # containing spatial clusters.
#' # set.seed(123)
#' # results_individual_svg <- individual_svg(LIBD_subset,
#' #                                            edgeR_y = results_svg_test$estimated_y,
#' #                                            spatial_cluster = "layer_guess_reordered")
#' # save(results_individual_svg, file = "./DESpace/data/results_individual_svg.RData")
#' @return A List of 7 elements - one element for each spatial cluster
#' @author Peiying Cai \email{peiying.cai@uzh.ch}, Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{svg_test}}, \code{\link{individual_svg}}
NULL