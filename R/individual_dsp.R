#' individual_dsp
#'
#' DESpace can also be used to reveal the specific areas of the tissue affected by DSP genes; i.e., spatial clusters that are particularly over/under abundant compared to the average signal across conditions.
#' This function can be used to identify SVGs among conditions for each individual cluster.
#'
#' @param spe SpatialExperiment or SingleCellExperiment.
#' @param cluster_col Character. Column name of spatial clusters in \code{colData(spe)}. 
#' @param sample_col Character. Column name of sample ids in \code{colData(spe)}.
#' Sample ids must be either a factor or character.
#' @param condition_col Character. Column name of condition ids in \code{colData(spe)}.
#' @param min_counts Numeric. Minimum number of counts per sample (across all spots) for a gene to be analyzed.
#' @param min_non_zero_spots Numeric. Minimum number of non-zero spots per sample, for a gene to be analyzed.
#' @param min_pct_cells Numeric. Minimum percentage of cells required for each cluster to be 
#' included in the analysis across the specified conditions. 
#' Default value is 0.5 (i.e., 0.5\% of total cells per cluster per condition). 
#' @param filter_gene Logical. If TRUE, 
#' \code{\link{dsp_test}} filters genes:
#' genes have to be expressed in at least 'min_non_zero_spots' spots, 
#' and a gene requires at least 'min_counts' counts per sample (across all locations). 
#' @param filter_cluster Logical. When set to TRUE, \code{\link{dsp_test}} excludes clusters that are 
#' insufficiently represented in the dataset. Only clusters meeting the 'min_pct_cells' threshold 
#' (i.e., containing at least the specified percentage of cells across all conditions) will be retained for analysis.
#' @return A list of results, with one result per spatial cluster in each element.
#' Specifically, each item in the list is a "gene_results" dataframe which contains main edgeR test results.
#' @examples
#' # load the input data:
#' eh <- ExperimentHub()
#' spe_example <- eh[["EH9613"]]
#' set.seed(123)
#' results_individual_dsp <- individual_dsp(spe_example,
#'                                           cluster_col = "Banksy_smooth",
#'                                           sample_col = "sample_id",
#'                                           condition_col = "condition")
#'                                            
#' # We visualize results for the cluster '3'
#' results <- results_individual_dsp[['3']]
#' head(results,3)
#' 
#' @seealso \code{\link{top_results}}, \code{\link{svg_test}}, \code{\link{dsp_test}}, \code{\link{FeaturePlot}}
#' 
#' @export
individual_dsp <- function(spe,
                            cluster_col,
                            sample_col,
                            condition_col,
                            min_counts = 20,
                            min_non_zero_spots = 10,
                            min_pct_cells = 0.5,
                            filter_gene = TRUE,
                            filter_cluster = TRUE){
    # check
    uniqueness_check <- .check_columns(spe, cluster_col, sample_col, condition_col)
    spe <- uniqueness_check[["updated_spe"]]
    # if rownames(spe) is null, column 'genes' would missing from the edgeR results
    # 'res_edgeR[[1]][c("genes", "LR", "logCPM", "PValue", "FDR")]' would return an error
    if(is.null(rownames(spe))){
      message("Gene names are missing in rownames(spe)")
      return(NULL)
    }
    # filtering
    if(filter_gene == TRUE){
      message("Filter low quality genes: \n")
      message("min_counts = ", min_counts, "; min_non_zero_spots = ", min_non_zero_spots, ".\n")
      # for the multi-sample case:
      # sample ids:
      sample_names <- levels(factor(colData(spe)[[sample_col]]))
      sel_matrix <- vapply(sample_names, function(id){
        spe <- subset(spe,,get(sample_col) == id)
        sel_1 <- rowSums(counts(spe)) >= min_counts
        sel_2 <- rowSums(counts(spe) > 0) >= min_non_zero_spots 
        sel_1 & sel_2
      }, FUN.VALUE = logical(nrow(spe)))
      # sel_matrix: column names are sample names; row names are genes
      # only keep genes that pass all filters across samples
      sel <- rowMeans(sel_matrix) == 1
      spe <- spe[sel, ]
      message("The number of genes that pass filtering is ", 
              dim(spe)[1], ".\n")
      # } # end for !is.null(num_sample)
    }# end for filter_gene == TRUE
    # aggregate data to pseudobulk via muscat  
    if (!"counts" %in% assayNames(spe)) {
      stop("The 'counts' assay slot, required as input data, is missing.")
    }
    layer <- factor(colData(spe)[[cluster_col]])
    # Re-label the cluster (e.g., cluster 1 vs. rest)
    results_res <- lapply(levels(layer), function(each_layer) {
      message(sprintf("Conducting tests for layer '%s' against all other layers.\n", each_layer))
      sce_one <- spe
      layer <- .cluster_label(sce_one, cluster_list = each_layer, cluster_col = cluster_col)
      sce_one[["individual_layer"]] <- layer
      res_edgeR <- .muscat_layer_test(sce_one,
                         design = NULL,
                         cluster_col = "individual_layer",
                         sample_col,
                         condition_col,
                         verbose = FALSE)
      return(res_edgeR)
    })
    names(results_res) <- levels(layer)
    message("Returning results")
    return(cluster_results = results_res)
}