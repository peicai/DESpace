#' DESpace_test
#' 
#' 'DESpace_test' identifies spatially variable genes (SVGs) from spatially-resolved transcriptomics data,
#' provided spatial clusters are available.
#' 
#' If 'sample_col' is not specified and 'replicates == FALSE',
#' \code{\link{DESpace_test}} assumed that data comes from an individual sample, 
#' and performs SV testing on it.
#' 
#' If 'sample_col' is provided and 'replicates == FALSE',
#' \code{\link{DESpace_test}} tests each sample individually and returns a list of results for each sample.
#' 
#' If 'sample_col' is provided and 'replicates == TRUE',
#' \code{\link{DESpace_test}} performs a joint multi-sample test.
#'  
#'
#' @param spe SpatialExperiment or SingleCellExperiment.
#' @param spatial_cluster Column name of spatial clusters in \code{colData(spe)}.
#' @param sample_col Column name of sample ids in \code{colData(spe)}.
#' @param replicates A logical, indicating whether biological replicates are provided (TRUE) or not (FALSE).
#' If biological replicates are provided, \code{\link{DESpace_test}} performs a joint test across all replicates, 
#' searching for SVGs with consistent spatial patterns across samples.
#' @param min_counts Minimum number of counts per sample (across all spots) for a gene to be analyzed.
#' @param min_non_zero_spots Minimum number of non-zero spots per sample, for a gene to be analyzed.
#' @param filter_gene A logical. If TRUE, 
#' \code{\link{DESpace_test}} filters genes:
#' genes have to be expressed in at least 'min_non_zero_spots' spots, 
#' and a gene requires at least 'min counts' counts per sample (across all locations). 
#' @param verbose A logical.
#'   If TRUE, \code{\link{DESpace_test}} returns two more results: 
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
#' # load the input data:
#' data("LIBD_subset", package = "DESpace")
#' LIBD_subset
#' 
#' # Fit the model via \code{\link{DESpace_test}} function.
#' set.seed(123)
#' results_DESpace <- DESpace_test(spe = LIBD_subset,
#'                            spatial_cluster = "layer_guess_reordered",
#'                            verbose = FALSE)
#' 
#' # DESpace_test returns of a list of 2 objects:
#' # "gene_results": a dataframe contains main edgeR test results;
#' # "estimated_y": a DGEList object contains the estimated common dispersion, 
#' #  which can later be used to speed-up calculation when testing individual clusters.
#' 
#' # We visualize differential results:
#' head(results_DESpace$gene_results, 3)
#' 
#' @seealso \code{\link{top_results}}, \code{\link{individual_test}}, \code{\link{FeaturePlot}}
#' 
#' @export
DESpace_test = function(spe,
                        spatial_cluster,
                        sample_col = NULL,
                        replicates = FALSE,
                        min_counts = 20,
                        min_non_zero_spots = 10,
                        filter_gene = TRUE,
                        verbose = FALSE) {
    genes = FDR = PValue = NULL
    message("using 'DESpace_test' for spatial gene/pattern detection.")
    # if rownames(spe) is null, column 'genes' would missing from the edgeR results
    # 'res_edgeR[[1]][c("genes", "LR", "logCPM", "PValue", "FDR")]' would return an error
    if(is.null(rownames(spe))){
        message("Gene names are missing in rownames(spe)")
    return(NULL)
    }
    # if replicates == TRUE, run a multi-sample edgeR test
    # results: a list; test results (DT_results_multi) and estimated dispersion (estimated_y_multi)
    if(!is.null(sample_col)){
        num_sample <- length(unique(colData(spe)$sample_id))
        num_sample <- if(num_sample > 1) num_sample else NULL
    if(sample_col != "sample_id" && "sample_id" %in% colnames(colData(spe))){
        message("The sample ids provided in the 'sample_col' column would replace 
                the 'sample_id' column in the colData(spe).")
        sample_id = colData(spe)[[sample_col]]
        colData(spe)$sample_id = sample_id}
    }else(num_sample = NULL)
    # filtering
    if(filter_gene == TRUE){
        message("Filter low quality genes: \n")
        message("min_counts = ", min_counts, "; min_non_zero_spots = ", min_non_zero_spots, ".\n")
    if(replicates == FALSE && is.null(num_sample)){
        sel_1 <- rowSums(counts(spe)) >= min_counts
        sel_2 <- rowSums(counts(spe) > 0) >= min_non_zero_spots
        spe = spe[sel_1 & sel_2, ]
        message("The number of genes that pass filtering is", 
                dim(spe)[1], ".\n")
    }else{
        # for the multi-sample case:
        # sample ids:
        sample_names = levels(factor(colData(spe)$sample_id))
        sel_matrix = vapply(sample_names, function(id){
                spe = subset(spe,,sample_id == id)
                sel_1 <- rowSums(counts(spe)) >= min_counts
                sel_2 <- rowSums(counts(spe) > 0) >= min_non_zero_spots 
                sel_1 & sel_2
        }, FUN.VALUE = logical(nrow(spe)))
        # sel_matrix: column names are sample names; row names are genes
        # only keep genes that pass all filters across samples
        sel = rowMeans(sel_matrix) == 1
        spe = spe[sel, ]
        message("The number of genes that pass filtering is", 
                dim(spe)[1], ".\n")
        } # end for !is.null(num_sample)
    }# end for filter_gene == TRUE
    if(replicates == TRUE){
        if( sample_col %notin% colnames(colData(spe))){
            sprintf("'sample_col' %s  not in colData(spe)", sample_col)
        return(NULL)}
    stopifnot(num_sample > 1)
    #nlevels(as.factor(with(colData(spe), get(sample_col)))) > 1
    if( spatial_cluster %notin% colnames(colData(spe))){
        sprintf("'spatial_cluster' %s  not in colData(spe)", spatial_cluster)
        return(NULL)
    }
    message("multi-sample test")
    if(verbose){
        list[gene_results, estimated_y, glmLRT, glmFit] = .multi_edgeR_test(spe = spe,
                                                                            spatial_cluster = spatial_cluster,
                                                                            sample_col = sample_col,
                                                                            verbose = verbose)
    }else{
        list[gene_results, estimated_y] = .multi_edgeR_test(spe = spe,
                                                            spatial_cluster = spatial_cluster,
                                                            sample_col = sample_col,
                                                            verbose = verbose)
    }
    DT_results = gene_results
    data.table::setorder(DT_results, FDR, PValue)
    DT_results_multi = DT_results
    estimated_y_multi = estimated_y
    }
    # if replicates == FALSE, run single-sample edgeR tests
    # if spe is a list of spe objects, run single-sample edgeR test for each sample via the function ".multi_single_edgeR_test"
    # results: a list; test results (DT_results_single) and estimated dispersion (estimated_y_single)
    if(replicates == FALSE && !is.null(num_sample) ){
    if(spatial_cluster %notin% colnames(colData(spe))){
        sprintf("'spatial_cluster' %s  not in colData(spe)", spatial_cluster)
        return(NULL)
    }
    if(sample_col %notin% colnames(colData(spe))){
        sprintf("'sample_col' %s  not in colData(spe)", sample_col)
        return(NULL)
    }
    message("using 'single' for spatial gene/pattern detection. ")
    sample_names = unique(colData(spe)$sample_id)
    n_sample = length(sample_names)
    results = data.frame()
    results = lapply(seq_len(n_sample), .multi_single_edgeR_test,
                    spe = spe,
                    sample_names = sample_names,
                    sample_col = sample_col,
                    spatial_cluster = spatial_cluster,
                    verbose = verbose
    )
    names(results) <- sample_names
    return(results)
    # results_all = do.call("rbind", results)
    # DT_results_single = results_all$DT_results1
    # estimated_y_single = results_all$estimated_y
    }
    # if replicates == FALSE, run a single-sample edgeR test
    # if spe is a spe object, run a single-sample edgeR test
    # results: a list; test results (DT_results_single) and estimated dispersion (estimated_y_single)
    if(replicates == FALSE && is.null(num_sample)){
        if( spatial_cluster %notin% colnames(colData(spe))){
        sprintf("'spatial_cluster' %s  not in colData(spe)", spatial_cluster)
        return(NULL)
        }
    message("single sample test")
    # if "logcounts" not included in assayNames(spe), calculate logcounts
    if(length(assayNames(spe)) == 1){
        counts <- assays(spe)$counts
        libsizes <- colSums(counts)
        size.factors <- libsizes/mean(libsizes)
        assay(spe, "logcounts") <- log2(t(t(counts)/size.factors) + 1)}
    if(verbose){
        list[gene_results, estimated_y, glmLRT, glmFit] = .single_edgeR_test(spe = spe, 
                                                                                spatial_cluster = spatial_cluster, 
                                                                                verbose = verbose)
    }else{
        list[gene_results, estimated_y] = .single_edgeR_test(spe = spe,
                                                                spatial_cluster = spatial_cluster,
                                                                verbose = verbose)
        }
    DT_results_single = gene_results
    estimated_y_single = estimated_y
    }
    # if replicates == TRUE, gene_results estimated_y -> from multi-sample results
    if(replicates == TRUE){
        gene_results = DT_results_multi
        estimated_y = estimated_y_multi
    }
    # if replicates == FALSE, gene_results, estimated_y -> from single-sample results
    if(replicates == FALSE){
        gene_results = DT_results_single
        estimated_y = estimated_y_single
    }
    if(verbose){
        results_list <- list(gene_results, estimated_y, glmLRT, glmFit)
        names(results_list) <- c("gene_results", "estimated_y", "glmLrt", "glmFit")
    }else{
        results_list <- list(gene_results, estimated_y)
        names(results_list) <- c("gene_results", "estimated_y")
    }
    return(results_list)
}