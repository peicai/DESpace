#' Results from \code{\link{DESpace_test}} function
#' 
#' @rdname results_DESpace_test
#' @name results_DESpace_test
#' @aliases results_DESpace_test
#' 
#' @param results_DESpace_test contains a \code{\linkS4class{list}} object,
#' with the results obtained applying \code{\link{DESpace_test}} function to an external dataset from the spatialLIBD package.
#' Below the code used to obtain 'results_DESpace_test'.
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
#' # Parameter `spe` specifies the input `SpatialExperiment` or `SingleCellExperiment` object, 
#' # while `spatial_cluster` defines the column names of `colData(spe)` containing spatial clusters. 
#' # To obtain all statistics, set `verbose` to `TRUE`.
#' # DESpace_test returns of a list of 4 objects:
#' # "gene_results": a dataframe contains main edgeR test results;
#' # "estimated_y": a DGEList object contains the estimated common dispersion, 
#' #  which can later be used to speed-up calculation when testing individual clusters.
#' # "glmFit": a DGEGLM object contains full statistics from "edgeR::glmFit"
#' # "glmLRT" a DGELRT object contains full statistics from "edgeR::glmLRT"
#' # 
#' # set.seed(123)
#' # results_DESpace_test <- DESpace_test(spe = spe3,
#' #                                          spatial_cluster = "layer_guess_reordered",
#' #                                          verbose = FALSE)
#' # 
#' # save(results_DESpace_test, file = "./DESpace/data/results_DESpace_test.RData")
#' 
#' @author Peiying Cai \email{peiying.cai@uzh.ch}, Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{DESpace_test}}
NULL