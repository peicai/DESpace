#' Results from \code{\link{individual_test}} function
#' 
#' @rdname results_individual_test
#' @name results_individual_test
#' @aliases results_individual_test
#' 
#' @param results_individual_test contains a \code{\linkS4class{list}} object,
#' with the results obtained applying \code{\link{individual_test}} function to an external dataset from the spatialLIBD package.
#' Below the code used to obtain 'results_individual_test'.
#
#' @examples
#' # load the input data:
#' # data("LIBD_subset", package = "DESpace")
#' # LIBD_subset
#' # load pre-computed results (obtained via `DESpace_test`)
#' # data("results_DESpace_test", package = "DESpace")
#' # results_DESpace_test
#' 
#' # Function `individual_test()` can be used to identify SVGs for each individual cluster.
#' # Parameter `spatial_cluster` indicates the column names of `colData(spe)` 
#' # containing spatial clusters.
#' # set.seed(123)
#' # results_individual_test <- individual_test(LIBD_subset,
#' #                                            edgeR_y = results_DESpace_test$estimated_y,
#' #                                            spatial_cluster = "layer_guess_reordered")
#' # save(results_individual_test, file = "./DESpace/data/results_individual_test.RData")
#' 
#' @author Peiying Cai \email{peiying.cai@uzh.ch}, Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{DESpace_test}}, \code{\link{individual_test}}
NULL