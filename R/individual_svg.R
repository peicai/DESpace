#' individual_svg
#'
#' DESpace can also be used to reveal the specific areas of the tissue affected by SVGs; i.e., spatial clusters that are particularly over/under abundant compared to the average signal.
#' This function can be used to identify SVGs for each individual cluster.
#'
#' For every spatial cluster we test, \code{edgeR} would normally re-compute the dispersion estimates based on the specific design of the test.
#' However, this calculation represents the majority of the overall computing time.
#' Therefore, to speed-up calculations, we propose to use the dispersion estimates which were previously computed for the gene-level tests.
#' This introduces a minor approximation which, in our benchmarks, does not lead to decreased accuracy.
#' If you want to use pre-computed gene-level dispersion estimates, set \code{edgeR_y} to 'estimated_y'.
#' Alternatively, if you want to re-compute dispersion estimates (significantly slower, but marginally more accurate option), leave edgeR_y empty.
#'
#' @param spe SpatialExperiment or SingleCellExperiment.
#' @param cluster_col Column name of spatial clusters in \code{colData(spe)}.
#' @param sample_col Column name of sample ids in \code{colData(spe)}.
#' @param edgeR_y Pre-estimated dispersion; if it's null, compute dispersion.
#' @param min_counts Minimum number of counts per sample (across all spots) for a gene to be analyzed.
#' @param min_non_zero_spots Minimum number of non-zero spots per sample, for a gene to be analyzed.
#' @param filter_gene A logical. If TRUE, 
#' \code{\link{svg_test}} filters genes:
#' genes have to be expressed in at least 'min_non_zero_spots' spots, 
#' and a gene requires at least 'min counts' counts per sample (across all locations). 
#' @param replicates Single sample or multi-sample test.
#' @param BPPARAM An optional parameter passed internally to bplapply.
#' We suggest using as many cores as the number of spatial clusters.
#' If unspecified, the script does not run in parallel.
#' Note that parallel coding performs better only when 
#' dispersion estimations are not provided beforehand.
#' Moreover, parallelizing the script will increase the memory requirement;
#' if memory is an issue, leave 'BPPARAM' unspecified and, hence, avoid parallelization.
#' @return A list of results, with one result per spatial cluster in each element.
#' Specifically, each item in the list is a "gene_results" dataframe which contains main edgeR test results.
#' @examples
#' # load the input data:
#' data("LIBD_subset", package = "DESpace")
#' LIBD_subset
#' 
#' # load pre-computed results (obtaines via `svg_test`)
#' data("results_svg_test", package = "DESpace")
#' 
#' # svg_test returns of a list of 2 objects:
#' # "gene_results": a dataframe contains main edgeR test results;
#' # "estimated_y": a DGEList object contains the estimated common dispersion, 
#' #  which can later be used to speed-up calculation when testing individual clusters.
#' 
#' # We visualize differential results:
#' head(results_svg_test$gene_results, 3)
#' 
#' # Individual cluster test: identify SVGs for each individual cluster
#' # set parallel computing; we suggest using as many cores as the number of spatial clusters.
#' # Note that parallelizing the script will increase the memory requirement;
#' # if memory is an issue, leave 'BPPARAM' unspecified and, hence, avoid parallelization.
#' set.seed(123)
#' results_individual_svg <- individual_svg(LIBD_subset, 
#'                                            edgeR_y = results_svg_test$estimated_y, 
#'                                            cluster_col = "layer_guess_reordered")
#'                                            
#' # We visualize results for the cluster 'WM'
#' results_WM <- results_individual_svg[[7]]
#' head(results_WM,3)
#' 
#' @seealso \code{\link{top_results}}, \code{\link{svg_test}}, \code{\link{FeaturePlot}}
#' 
#' @export
individual_svg <- function(spe,
                            cluster_col,
                            sample_col = "sample_id",
                            edgeR_y=NULL,
                            min_counts = 20,
                            min_non_zero_spots = 10,
                            filter_gene = TRUE,
                            replicates = FALSE,
                            BPPARAM = NULL){
    if( cluster_col %notin% colnames(colData(spe))){
        sprintf("'cluster_col' %s  not in colData(spe)", cluster_col)
        return(NULL)
    }
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
    message("The number of genes that pass filtering is", 
        dim(spe)[1], ".\n")
    # } # end for !is.null(num_sample)
    }# end for filter_gene == TRUE
    layer <- factor(colData(spe)[[cluster_col]])
    layer  <- droplevels(layer)
    spe <- spe[, !is.na(layer)]
    if(!is.null(edgeR_y) && !identical(dim(spe), dim(edgeR_y))){
        message("The input SpatialExperiment and pre-computed estimation do not have same dimension !")
        stop("In order to use gene-level dispersion estimates, 
            the filter has to be identical.")
    }
    # Re-label the cluster (e.g., cluster 1 vs. rest)
    message("Pre-processing")
    layer_list <- lapply(seq_len(nlevels(as.factor(layer))), function(x){
        cluster <- levels(as.factor(layer))[x]
        .cluster_label(spe, cluster_list = cluster, cluster_col = cluster_col)
    })
    message("Start modeling")
    # if edgeR_y != NULL -> use dispersion taken from gene-level test (with all clusters)
    if(!is.null(edgeR_y)){
        if(replicates == TRUE){
            if( sample_col %notin% colnames(colData(spe))){
                sprintf("'sample_col' %s  not in colData(spe)", sample_col)
                return(NULL)
            }
        sample_id <-  factor(colData(spe)[[sample_col]])
        if(is.null(BPPARAM)){
            result_list <- lapply(layer_list, .layer_test1_multi,
                                    y = edgeR_y, sample_id = sample_id)
        }else{
            result_list <- bplapply(layer_list, .layer_test1_multi,
                                    y = edgeR_y, sample_id = sample_id, BPPARAM = BPPARAM)
        }
    }else{
        if(is.null(BPPARAM)){
            result_list <- lapply(layer_list, .layer_test1, y = edgeR_y)
        }else{
            result_list <- bplapply(layer_list, .layer_test1, y = edgeR_y, BPPARAM = BPPARAM)
            }
        }
    }else if(is.null(edgeR_y)){
    # if edgeR_y == NULL -> dispersion computed each time
    message("Input data (estimated dispersion) is missing.")
    message("Re-compute dispersion.")
    edgeR_y <- DGEList(counts=assays(spe)$counts,
                        genes=rownames(assays(spe)$counts))
    edgeR_y$samples$lib.size <- colSums(edgeR_y$counts)
    edgeR_y <- calcNormFactors(edgeR_y)
    if(replicates == TRUE){
        sample_id  <-  factor(colData(spe)[[sample_col]])
        if(is.null(BPPARAM)){
            result_list <- lapply(layer_list, .layer_test2_multi,
                                    y = edgeR_y, sample_id = sample_id)
        }else{
            result_list <- bplapply(layer_list, .layer_test2_multi,
                                    y = edgeR_y, sample_id = sample_id, BPPARAM = BPPARAM)
            }
    }else{
        if(is.null(BPPARAM)){
            result_list <- lapply(layer_list, .layer_test2, y = edgeR_y)
        }else{
            result_list <- bplapply(layer_list, .layer_test2, y = edgeR_y, BPPARAM = BPPARAM)
            }
        }
    }
    names(result_list) <- levels(as.factor(layer))
    message("Returning results")
    return(cluster_results = result_list)
}
