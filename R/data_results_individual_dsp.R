#' Results from \code{\link{individual_dsp}} function
#' 
#' @rdname results_individual_dsp
#' @name results_individual_dsp
#' @aliases results_individual_dsp
#' 
#' @param results_individual_dsp contains a \code{\linkS4class{list}} object,
#' with the results obtained applying \code{\link{individual_dsp}} function to an external dataset from the spatialLIBD package.
#' Below the code used to obtain 'results_individual_dsp'.
#
#' @examples
#' # load the input data:
#' # data("LIBD_multi", package = "DESpace")
#' # LIBD_multi
#' 
#' # Function `individual_dsp()` can be used to identify SVGs for each individual cluster.
#' # Parameter `spatial_cluster` indicates the column names of `colData(spe)` 
#' # containing spatial clusters.
#' # set.seed(123)
#' # results_individual_dsp <- individual_dsp(LIBD_multi,
#' #                                           cluster_col = "layer_guess_reordered",
#' #                                           sample_col = "sample_id",
#' #                                           condition_col = "condition")
#' @return A List of 7 elements - one element for each spatial cluster
#' @author Peiying Cai \email{peiying.cai@uzh.ch}, Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{individual_dsp}}, \code{\link{dsp_test}}
NULL