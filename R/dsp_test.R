#' dsp_test
#' 
#' 'dsp_test' identifies differential spatial pattern (DSP) genes between conditions from spatially-resolved transcriptomics data,
#' provided spatial clusters are available.
#'
#' @param spe SpatialExperiment or SingleCellExperiment.
#' @param design Matrix or array. Numeric design matrix for a regression-like model created by `model.matrix` function. 
#' @param cluster_col Character. Column name of spatial clusters in \code{colData(spe)}. 
#' @param sample_col Character. Column name of sample ids in \code{colData(spe)}.
#' Sample ids must be either a factor or character.
#' @param condition_col Character. Column name of condition ids in \code{colData(spe)}.
#' @param min_counts Numeric. Minimum number of counts per sample (across all spots) for a gene to be analyzed.
#' @param min_non_zero_spots Numeric. Minimum number of non-zero spots per sample, for a gene to be analyzed.
#' @param min_pct_cells Numeric. Minimum percentage of cells required for each cluster to be 
#' included in the analysis across the specified conditions. 
#' Default value is 0.5 (i.e., 0.5\% of total cells per cluster per condition). 
#' @param filter_gene Logical. If TRUE, \code{\link{dsp_test}} filters genes by requiring them to be expressed 
#' in at least 'min_non_zero_spots' cells and have at least 'min_counts' counts per sample across all locations. 
#' @param filter_cluster Logical. When set to TRUE, \code{\link{dsp_test}} excludes clusters that are 
#' insufficiently represented in the dataset. Only clusters meeting the 'min_pct_cells' threshold 
#' (i.e., containing at least the specified percentage of cells across all conditions) will be retained for analysis.
#' @param verbose Logical.
#'   If TRUE, \code{\link{svg_test}} returns two more results: 
#'   'DGEGLM' and 'DGELRT' objects contain full statistics from 'edgeR::glmFit' and 'edgeR::glmLRT'.
#' @return A list of results:
#' 
#' - "gene_results": a dataframe contains main edgeR test results;
#' 
#' - "estimated_y": a DGEList object contains the estimated common dispersion, 
#' which can later be used to speed-up calculation when testing individual clusters.
#' 
#' - "glmFit" (only if \code{verbose = TRUE}): a DGEGLM object contains full statistics from "edgeR::glmFit".
#' 
#' - "glmLRT" (only if \code{verbose = TRUE}): a DGELRT object contains full statistics from "edgeR::glmLRT".
#'
#' @examples
#' ## Load the example multi-sample multi-group spe object
#' ## from the muSpaData package
#' muSpaData::Wei22_example()
#' # Fit the model via \code{\link{dsp_test}} function.
#' set.seed(123)
#' results_dsp <- dsp_test(spe = spe_example,
#'                            cluster_col = "Banksy_smooth",
#'                            sample_col = "sample_id",
#'                            condition_col = "condition",
#'                            verbose = FALSE)
#' 
#' # dsp_test returns of an object:
#' # "gene_results": a dataframe contains main edgeR test results.
#' 
#' # We visualize differential results:
#' head(results_dsp, 3)
#' 
#' @seealso \code{\link{svg_test}}, \code{\link{individual_svg}}, \code{\link{individual_dsp}}, \code{\link{FeaturePlot}}, \code{\link{top_results}}
#' 
#' @export

dsp_test <-  function(spe,
                      design = NULL,
                      cluster_col,
                      sample_col,
                      condition_col,
                      min_counts = 20,
                      min_non_zero_spots = 10,
                      min_pct_cells = 0.5,
                      filter_gene = FALSE,
                      filter_cluster = TRUE,
                      verbose = FALSE) {
    message("Using 'dsp_test' for spatial variable pattern genes detection.\n")
    # check
    uniqueness_check <- .check_columns(spe, cluster_col, sample_col, condition_col)
    spe <- uniqueness_check[["updated_spe"]]
    # filter genes
    if(filter_gene){
      message("Filter low quality genes: \n")
      # message("min_counts = ", min_counts, "; min_non_zero_spots = ", min_non_zero_spots, ".\n")
      sample_id <- factor(colData(spe)[[sample_col]])
      sample_names <- levels(sample_id)
      sel_matrix <- vapply(sample_names, function(id){
        select = sample_id == id
        spe <- spe[,select]
        sel_1 <- rowSums(counts(spe)) >= min_counts
        sel_2 <- rowSums(counts(spe) > 0) >= min_non_zero_spots 
        sel_1 & sel_2
      }, FUN.VALUE = logical(nrow(spe)))
      # sel_matrix: column names are sample names; row names are genes
      # only keep genes that pass all filters across samples
      sel <- rowMeans(sel_matrix) == 1
      spe <- spe[sel, ]
      message("The number of genes that pass filtering is ", dim(spe)[1], ".\n")
    }
    if(filter_cluster) {
      message("Filter low quality clusters: ")
      # message("The minimum percentage of cells required for each cluster: ", min_pct_cells, "%.\n")
      spe <- .filter_clusters_by_pct(spe, condition_col, cluster_col, min_pct_cells, uniqueness_check[["n_condition"]])
    }
    # aggregate data to pseudobulk via muscat  
    if (!"counts" %in% assayNames(spe)) {
      stop("The 'counts' assay slot, required as input data, is missing.")
    }
    .muscat_layer_test(spe, design = design, cluster_col, sample_col, condition_col, verbose)
}