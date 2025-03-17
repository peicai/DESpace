#' top_results
#' 
#' Filter significant results.
#' \code{\link{top_results}} returns the significant results obtained via \code{\link{svg_test}} and \code{\link{individual_svg}}.
#' It can also be used to merge gene- and cluster-level results into a single object.
#' @param gene_results  Results returned from \code{\link{svg_test}}.
#' @param cluster_results Results returned from \code{\link{individual_svg}}.
#' @param cluster A character indicating the cluster(s) whose results have to be returned.
#' Results from all clusters are returned by default ("NULL").
#' @param select A character indicating what results should be returned ("FDR", "logFC", or "both").
#' Only used if "cluster_results" are provided.
#' By default ("both"), both FDR and logFC are returned.
#' @param high_low 
#' A character indicating whether to filter results or not.
#' Only used if "cluster_results" are provided, and one cluster is specified in "cluster" parameter.
#' By default (NULL), all results are returned in a single data.frame.
#' If "high" or "HIGH", we only return SVGs with average abundace in "cluster" higher than in the rest of the tissue (i.e., logFC > 0).
#' If "low" or "LOW", we only return SVGs with average abundace in "cluster" lower than in the rest of the tissue (i.e., logFC < 0).
#' If "both" or "BOTH", then both "high" and "low" results are returned, but in two separate data.frames.
#' @return A \code{\linkS4class{data.frame}} object or a list of \code{\linkS4class{data.frame}} with results.
#' 
#' - When only “cluster_results” is provided, results are reported as a \code{\linkS4class{data.frame}} with columns for
#'  gene names (gene_id), spatial clusters affected by SV (Cluster), cluster-specific likelihood ratio test statistics (LR),
#'  cluster-specific average (across spots) log-2 counts per million (logCPM), cluster-specific log2-fold changes (logFC),
#'  cluster-specific raw p-values (PValue), and Benjamini-Hochberg adjusted p-values (FDR) for each spatial cluster.
#'
#' - When “gene_results” and “cluster_results” are given, results are reported as a \code{\linkS4class{data.frame}} that merges gene-
#'  and cluster-level results.
#'
#' - If “cluster” is specified, the function returns a subset \code{\linkS4class{data.frame}} for the given cluster, which contains cluster name,
#'  gene name, LR, logCPM, logFC, PValue and FDR, ordered by FDR for the specified cluster.
#'
#' - If “high_low” is set, the function returns a list of \code{\linkS4class{data.frame}} that contains subsets of results for genes with
#'  higher and/or lower expression in the given cluster compared to the rest of the tissue.

#' @examples
#' # load pre-computed results (obtained via `svg_test`)
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
#' # load pre-computed results (obtained via `individual_svg`)
#' data("results_individual_svg", package = "DESpace")
#' # Function `individual_svg()` can be used to identify SVGs for each individual cluster.
#' # `individual_svg()` returns a list containing the results of individual clusters.
#' # For each cluster, results are reported as a data.frame, 
#' # where columns For each cluster, results are reported as a data.frame, 
#' # where columns contain gene names (`genes`), likelihood ratio (`LR`), 
#' # log2-fold changes (`logFC`) and adjusted p-value (`FDR`).
#' # 
#' # Combine gene-and cluster-level results
#' merge_res = top_results(results_svg_test$gene_results, 
#'                         results_individual_svg)
#' head(merge_res,3)
#' # 'select = "FDR"' can be used to visualize adjusted p-values for each spatial cluster.
#' merge_res = top_results(results_svg_test$gene_results, 
#'                         results_individual_svg, select = "FDR")
#' head(merge_res,3)
#' # Specify the cluster of interest and check top genes detected by svg_test.
#' results_WM_both = top_results(cluster_results = results_individual_svg, 
#'                               cluster = "WM", high_low = "both")
#' head(results_WM_both$high_genes, 3)
#' head(results_WM_both$low_genes, 3)
#' 
#' @seealso \code{\link{svg_test}}, \code{\link{individual_svg}}, \code{\link{FeaturePlot}}, \code{\link{dsp_test}}, \code{\link{individual_dsp}}
#' 
#' @export
top_results <-  function(gene_results = NULL,
                        cluster_results,
                        cluster = NULL,
                        select = "both",
                        high_low = NULL){
    stopifnot(
        is.list(cluster_results), is.character(select), length(select) == 1L
    )
    cluster_results <- lapply(cluster_results, as.data.frame)
    # If there is gene level results (stored in gene_results), return merged results (merge gene and cluster level results)
    if(!is.null(gene_results)){
        gene_results <- as.data.frame(gene_results)
        merge_results <- .merge_results(gene_results, cluster_results, select)
        merge_results <- as.data.frame(merge_results) %>% arrange(gene_FDR)
        return(res = merge_results)
    }
    if(!is.null(high_low)){
        stopifnot(is.character(high_low), length(high_low) == 1L)
        if( !(high_low %in% c("both", "BOTH", "high", "HIGH", "low", "LOW")) ){
            message("'high_low' should be one of: 'both', 'BOTH', 'high', 'HIGH', 'low', 'LOW'.")
            return(NULL)
        }
        if( is.null(cluster) ){
            message("The spatial cluster is missing in `cluster` argument.")
            return(NULL)
        }
    }
    if(!is.null(cluster)){
        stopifnot(is.character(cluster), length(cluster) == 1L)
    }
    # Extract and combine results of the selected column ('select') from each sample
    com_results <- .getValueRes(cluster_results = cluster_results, select = 'FDR')
    # Top genes -> for each gene/row, min p-values/FDR; for each gene, return a vector of layer names
    sel_layer_min <- apply(com_results,1,which.min)
    rank_results<-as.data.frame(colnames(com_results)[sel_layer_min])
    colnames(rank_results) <- "Cluster"
    # For each gene, extract LR, PValue and FDR based on the layer name (i.e., the layer has min FDR)
    rank_results <- cbind(rank_results,
                        lapply(seq_len(nrow(com_results)), function(x) {
                            cluster_results[[sel_layer_min[x]]] %>% 
                            filter(gene_id==rownames(com_results)[x]) }) %>%
                        bind_rows())
    # Rename columns; 3:7 -> "LR", "logCPM", "logFC", "PValue", "FDR"
    colnames(rank_results)[seq(3,7)] <- paste0("Cluster_", colnames(rank_results)[seq(3,7)])
    colnames(com_results) <- paste0(names(cluster_results), "_FDR")
    # Add the "gene id" column
    com_results['gene_id'] <- rownames(com_results)
    rownames(rank_results) <- rank_results$gene_id
    # Add the results of identified layers (i.e., rank_results) into com_results
    com_results <- merge(rank_results, com_results, by = "gene_id")
    # Check genes: p-values/FDR for all layers are same
    rownames(com_results) <- com_results$gene_id
    # Sort results based on FDR
    com_results <- as.data.frame(com_results) %>% arrange(Cluster_FDR)
    # Sort results based on the specific cluster i; return top genes for cluster i
    if(!is.null(cluster)){
        gene_id <- com_results[order(com_results[, paste0((cluster), "_FDR")], decreasing = FALSE), 'gene_id']
        # Reorder rows of rank_results based on FDR of the specific cluster i
        top_results_cluster <- rank_results[match(gene_id, rank_results$gene_id),]
        # Subset Cluster == cluster i
        top_genes_cluster <- subset(top_results_cluster, Cluster == (cluster))
    }
    # Sort results based on the specific cluster i; return top genes for cluster i
    if(!is.null(cluster) & !is.null(high_low)){
        sel_FC <- .getValueRes(cluster_results = cluster_results, select = 'logFC')
        if(high_low %in% c("HIGH","high")){
            FC_high_genes <- rownames(sel_FC[sel_FC[,(cluster)] >= 0,])
            high_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_id %in% FC_high_genes,])
            low_genes <- NULL
        }
        if(high_low %in% c("LOW","low")){
            FC_low_genes <- rownames(sel_FC[sel_FC[,(cluster)] < 0,])
            low_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_id %in% FC_low_genes,])
            high_genes <- NULL
        }
        if(high_low %in% c("BOTH","both")){
            FC_high_genes <- rownames(sel_FC[sel_FC[,(cluster)] >= 0,])
            high_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_id %in% FC_high_genes,])
            FC_low_genes <- rownames(sel_FC[sel_FC[,(cluster)] < 0,])
            low_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_id %in% FC_low_genes,])
        }
    }
    if(is.null(high_low)){
        if(is.null(cluster)){
            return(res = com_results)
        }else{
            return(top_genes = top_genes_cluster)
        }
    }else{
        return(list(
        high_genes = high_genes,
        low_genes = low_genes
        ))
    }
}