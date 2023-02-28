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
#' # Connect to ExperimentHub
#' # ehub <- ExperimentHub::ExperimentHub()
#' # Download the example spe data
#' # spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)
#' # spe_all
#' # 
#' # Only use one sample:
#' # spe3 <- spe_all[, colData(spe_all)$sample_id == '151673']
#' # Select small set of random genes for faster runtime in this example
#' # set.seed(123)
#' # sel_genes <- sample(dim(spe3)[1],500)
#' # spe3 <- spe3[sel_genes,]
#' # 
#' # Fit the model via `DESpace_test` function. 
#' # set.seed(123)
#' # results_DESpace_test <- DESpace_test(spe = spe3,
#' #                                          spatial_cluster = "layer_guess_reordered",
#' #                                          verbose = TRUE)
#' #
#' # Function `individual_test()` can be used to identify SVGs for each individual cluster.
#' # Parameter `spatial_cluster` indicates the column names of `colData(spe)` 
#' # containing spatial clusters.
#' # set.seed(123)
#' # results_individual_test <- individual_test(spe3,
#' #                                            edgeR_y = results_DESpace_test$estimated_y,
#' #                                            spatial_cluster = "layer_guess_reordered")
#' # save(results_individual_test, file = "./DESpace/data/results_individual_test.RData")
#' 
#' @author Peiying Cai \email{peiying.cai@uzh.ch}, Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{DESpace_test}}, \code{\link{individual_test}}
NULL