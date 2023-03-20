#' Subset from the human DLPFC 10x Genomics Visium dataset of the \code{spatialLIBD} package
#' 
#' @rdname LIBD_subset
#' @name LIBD_subset
#' @aliases LIBD_subset
#' 
#' @param LIBD_subset contains a \code{\linkS4class{SpatialExperiment}} object,
#' representing a subset of the sample 151673 from the full real data of the \code{spatialLIBD} package.
#' Below the code used to subset the original dataset.
#
#' @examples
#' # Connect to ExperimentHub
#' # ehub <- ExperimentHub::ExperimentHub()
#' # Download the example spe data
#' # spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)
#' # Select one sample only:
#' # LIBD_subset <- spe_all[, colData(spe_all)$sample_id == '151673']
#' # Select small set of random genes for faster runtime 
#' # set.seed(123)
#' # sel_genes <- sample(dim(LIBD_subset)[1],500)
#' # LIBD_subset <- LIBD_subset[sel_genes,]
#' # keep_col <- c("array_row","array_col","layer_guess_reordered")
#' # library(SingleCellExperiment)
#' # LIBD_subset <- SpatialExperiment(assay = list(counts = assay(LIBD_subset),
#' #                                               logcounts = logcounts(LIBD_subset)), 
#' #                                  colData = colData(LIBD_subset)[keep_col])
#' # save(LIBD_subset, file = "./DESpace/data/LIBD_subset.RData")
#' 
#' @author Peiying Cai \email{peiying.cai@uzh.ch}, Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{DESpace_test}}, \code{\link{individual_test}}
NULL