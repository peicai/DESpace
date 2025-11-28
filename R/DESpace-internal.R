.cluster_label <-
    function(spe, 
            cluster_list, 
            cluster_col = "layer_guess_reordered"){
    metadata <- as.data.frame(colData(spe))
    one_layer <- factor(
        ifelse(metadata[[cluster_col]] %in% cluster_list, 
               as.character(metadata[[cluster_col]]), 
               "Other"), 
        levels = c("Other", cluster_list)
    )
    return(one_layer)
    }

.getValueRes <-
    function(cluster_results = single_cluster_results, select = "FDR"){
    ll <- lapply(seq_len(length(cluster_results)), function(i){
            data.frame(gene_id = cluster_results[[i]]$gene_id, 
            value = get(select, cluster_results[[i]]))})
    names(ll) <- names(cluster_results)
    
    ll2 <- lapply(ll, function(LL) lapply(LL, `[`, order(LL$gene_id)) )
    genesID <- ll2[[1]]$gene_id
    ll2 <-lapply(ll2, "[",  "value")
    
    com_results <- as.data.frame(do.call(cbind, lapply(ll2, as.data.frame)))
    rownames(com_results) <- genesID
    colnames(com_results) <- names(cluster_results)
    
    return(com_results)
}
.layer_test1_multi <-
    function(layers, 
            sample_id, 
            y){
    one_layer <- layers
    design_model <- model.matrix(~one_layer+sample_id)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results <- as.data.frame(results[[1]][c("genes", "LR", 
                "logCPM", "logFC", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    return(results)
}
.layer_test1 <-
    function(layers, 
            y){
    one_layer <- layers
    design_model <- model.matrix(~one_layer)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results <- as.data.frame(results[[1]][c("genes","LR", 
                "logCPM", "logFC", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    return(results)
}
.layer_test2 <-
    function(layers, 
            y){
    one_layer <- layers
    design_model <- model.matrix(~one_layer)
    is.fullrank(design_model)
    rownames(design_model) <- colnames(y)
    y <- estimateDisp(y, robust=TRUE, design = design_model)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results <- as.data.frame(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    return(results)
}
.layer_test2_multi <-
    function(layers, 
            sample_id, 
            y){
    one_layer <- layers
    design_model <- model.matrix(~one_layer)
    is.fullrank(design_model)
    rownames(design_model) <- colnames(y)
    y <- estimateDisp(y, robust=TRUE, design = design_model)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results <- as.data.frame(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    return(results)
    }
.merge_results <-
    function(gene_results,
            cluster_results, 
            select = "both"){
    stopifnot(
        is.list(cluster_results),
        is.data.frame(gene_results),
        is.character(select), length(select) == 1L
    )
    if( !(select %in% c("p_adj", "FDR", "logFC","both")) ){
        message("'select' should be one of: 'p_adj', 'FDR', 'logFC','both'")
        return(NULL)
    }
    colnames(gene_results) <- c("gene_id", "gene_LR", "gene_logCPM", "gene_Pvalue","gene_FDR")
    if(select != 'both'){
        com_results <- .getValueRes(cluster_results = cluster_results, select = select)
        colnames(com_results) <- paste0(colnames(com_results), "_", select)
        com_results['gene_id'] <- rownames(com_results)
        # Merge gene-level and cluster-level results
        merge_results <- merge(gene_results, com_results, by = "gene_id")
    }
    if(select == 'both'){
        com_results1 <- .getValueRes(cluster_results = cluster_results, select = c('FDR'))
        colnames(com_results1) <- paste0(colnames(com_results1), "_", 'FDR')
        com_results2 <- .getValueRes(cluster_results = cluster_results, select = c('logFC'))
        colnames(com_results2) <- paste0(colnames(com_results2), "_", 'logFC')
        com_results <- cbind(com_results1, com_results2)
        com_results['gene_id'] <- rownames(com_results)
        # Merge gene-level and cluster-level results
        merge_results <- merge(gene_results, com_results, by = "gene_id")
    }
    return(merge_results)
}
.multi_edgeR_test <-
    function(spe, 
            cluster_col, 
            sample_col, 
            # num_core,
            verbose = TRUE){
    layer <- as.factor(colData(spe)[[cluster_col]] )
    layer  <- droplevels(layer)
    spe <- spe[, !is.na(layer)]
    design <-  data.frame(condition = factor(colData(spe)[[cluster_col]]),
                          sample_id = factor(colData(spe)[[sample_col]]))
    design$condition <- droplevels(design$condition)
    y <- DGEList(counts=assays(spe)$counts, 
                genes=rownames(assays(spe)$counts))
    y$samples$lib.size <- colSums(y$counts)
    y <- calcNormFactors(y)
    design_model <- model.matrix(~design$condition + design$sample_id)
    rownames(design_model) <- colnames(y)
    y <- estimateDisp(y, design_model, robust=TRUE)
    fit <- glmFit(y, design_model)
    q <- nlevels(factor(colData(spe)[[cluster_col]]))
    lrt <- glmLRT(fit, coef = seq(2,q))
    res_edgeR <- topTags(lrt, n = Inf)
    results <- as.data.frame(res_edgeR[[1]][c("genes", "LR", "logCPM", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    if(verbose){
        return(list(gene_results = results,
                    estimated_y = y,
                    glmLRT = lrt,
                    glmFit = fit))
    }else{
        return(list(gene_results = results,
                    estimated_y = y))
    }   
}
.multi_single_edgeR_test <-
    function(i,
            spe,
            sample_names,
            sample_col,
            cluster_col,
            verbose = TRUE){
    spe <- subset(spe,,sample_id == sample_names[i])
    if(verbose){
        list[gene_results, estimated_y, lrt, fit]  <-  .single_edgeR_test(spe, cluster_col, verbose)
    }else{
        list[gene_results, estimated_y]  <-  .single_edgeR_test(spe, cluster_col, verbose)
}
    DT_results1 <- as.data.frame(gene_results)
    DT_results1$sample <- sample_names[i]
    setorder(DT_results1, FDR, PValue)
    if(verbose){
        return(list(DT_results1 = DT_results1,
                    estimated_y = estimated_y,
                    glmLRT = lrt,
                    glmFit = fit))
    }else{
        return(list(DT_results1 = DT_results1,
                    estimated_y = estimated_y))
    }
}
.single_edgeR_test <-
    function(spe, 
             cluster_col, 
            verbose = TRUE){
    layer <- as.factor(colData(spe)[[cluster_col]] )
    spe <- spe[, !is.na(layer)]
    design <- data.frame(condition = layer)
    design$condition <- droplevels(design$condition)
    y <- DGEList(counts=assays(spe)$counts, 
                genes=rownames(assays(spe)$counts))
    y$samples$lib.size <- colSums(y$counts)
    y <- calcNormFactors(y)
    design_model <- model.matrix(~design$condition)
    rownames(design_model) <- colnames(y)
    y <- estimateDisp(y, design_model, robust=TRUE)
    fit <- glmFit(y, design_model)
    q <- nlevels(factor(colData(spe)[[cluster_col]]))
    lrt <- glmLRT(fit, coef = seq(2,q))
    res_edgeR <- topTags(lrt, n = Inf)
    results <- as.data.frame(res_edgeR[[1]][c("genes", "LR", "logCPM", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    if(verbose){
        return(list(gene_results = results,
                    estimated_y = y,
                    glmLRT = lrt,
                    glmFit = fit))
    }else{
        return(list(gene_results = results,
                    estimated_y = y))
    }
    }
.check_columns <- function(spe, cluster_col, sample_col, condition_col) {
  # Columns to check for existence in colData(spe)
  columns_to_check <- c(cluster_col, sample_col, condition_col)
  # Check for missing columns
  missing_columns <- columns_to_check[!columns_to_check %in% colnames(colData(spe))]
  if (length(missing_columns) > 0) {
    # If there are missing columns, stop and show an error
    message <- sprintf("The following columns are not in colData(spe): %s", paste(missing_columns, collapse = ", "))
    stop(message)
  }
  updated_spe <- spe
  
  # Remove rows with NA values
  rows_with_na <- apply(as.data.frame(colData(updated_spe)), 1, anyNA)
  updated_spe <- updated_spe[ , !rows_with_na]

  results_list <- lapply(columns_to_check, function(col) {
    .get_nlevels(updated_spe, col)  # Call .get_nlevels with the current spe
  })
  # Extract updated_spe from the last iteration
  updated_spe <- results_list[[length(results_list)]]$updated_spe
  n_levels <- vapply(results_list, function(x) x$n_levels, numeric(1))
  # Get the number of levels for sample, condition, and cluster columns
  n_cluster <- n_levels[1]
  n_sample <- n_levels[2]
  n_condition <- n_levels[3]
  # Ensure all specified columns have more than one level
  if (any(c(n_sample, n_condition, n_cluster) <= 1)) {
    stop("All specified columns must have more than one level.")
  }
  return(list(n_cluster = n_cluster, n_sample = n_sample, n_condition = n_condition, updated_spe = updated_spe))
}
.get_nlevels <- function(spe, col) {
  col_data <- colData(spe)[[col]]
  if (!is.factor(col_data) && !is.character(col_data)) {
    stop(sprintf("Column '%s' must be either a factor or character.", col))
  }
  colData(spe)[[col]] <- droplevels(factor(col_data))
  n_levels <- nlevels(colData(spe)[[col]])
  return(list(n_levels = n_levels, updated_spe = spe))
}
.filter_clusters_by_pct <- function(spe, condition_col, cluster_col, min_pct_cells, n_condition) {
  # Create a table of cell counts per condition and cluster
  cluster_table <- table(colData(spe)[[condition_col]], colData(spe)[[cluster_col]])
  # Calculate total cells per condition
  total_cells_per_condition <- rowSums(cluster_table)
  # Calculate the percentage of cells per cluster per condition
  cluster_pct <- sweep(cluster_table, 1, total_cells_per_condition, "/") * 100
  # Filter out clusters with less than the specified percentage in either condition
  filtered_clusters <- colnames(cluster_pct)[colSums(cluster_pct >= min_pct_cells) == n_condition]
  message("Cluster levels to keep: ", paste(filtered_clusters, collapse = ", "))
  # Find the indices of the clusters to keep
  cluster_to_keep <- which(colData(spe)[[cluster_col]] %in% filtered_clusters)
  # Subset the SingleCellExperiment object
  spe_filtered <- spe[, cluster_to_keep]
  colData(spe_filtered)[[cluster_col]] <- droplevels(colData(spe_filtered)[[cluster_col]])
  # Return the filtered SingleCellExperiment object
  return(spe_filtered)
}
.layer_test <- function(spe,
                               design = NULL,
                               cluster_col,
                               sample_col,
                               condition_col,
                               test,
                               verbose){
  pb <- .aggregateData(spe, assay = "counts", fun = "sum",
                      by = c(cluster_col, sample_col))
  X <- do.call(cbind, assays(pb))
  # Create a grid of all combinations
  combinations <- expand.grid(colnames(assays(pb)[[1]]), names(assays(pb)))
  pasted_combinations <- apply(combinations, 1, function(x) paste0(x[1], "_", x[2]))
  colnames(X) <- pasted_combinations
  n_samples <- ncol(pb)
  n_clusters <- length(assays(pb))
  X_df <- as.data.frame(X)
  cols_with_zeros <- vapply(X_df, function(col) all(col == 0), logical(1))
  X_df <- X_df[, !cols_with_zeros]
  if(is.null(design)){
      design <- data.frame(condition = rep(pb[[condition_col]], n_clusters),
                          cluster_id = rep(names(assays(pb)), each = n_samples))
      design <- design[!cols_with_zeros,]
      # Relevel cluster_id to "Other" as baseline, if "Other" is present
      design[["cluster_id"]] <- factor(design[["cluster_id"]])
      if("Other" %in% design$cluster_id){
          design[["cluster_id"]] <- relevel(design[["cluster_id"]], ref = "Other")
      }
      design_model <- model.matrix(~ condition * cluster_id, data = design)
      rownames(design_model) <- colnames(X_df)
      message("Design model: row names represent sample names, followed by underscores and cluster names.\n")
      show(head(design_model, n = 2))
      col_with_colon <- grep("condition.*:*.cluster_id", colnames(design_model), value = TRUE)
  }else{
      design <- as.matrix(design)
      design_model <- design
      # Provide information about the pseudo-bulk count matrix and design validation
      message(
        "Pseudo-bulk count matrix dimensions: ", paste(dim(X_df), collapse = " x "), ".\n",
        "Columns are combinations of sample names and cluster names, excluding clusters not present in a sample.\n",
        "Column names:\n"
      )
      show(colnames(X_df))
      # Check if the 'design' matrix matches the pseudo-bulk count matrix
      message(
        "Please ensure that the rows in the 'design' matrix match the columns of the pseudo-bulk count matrix, both in number and order of names.\n",
        "Additionally, ensure that the ':' symbol only appears in interaction terms, not in main effect terms.\n"
      )
      # Validate dimensions
      if (dim(design)[1] != dim(X_df)[2]) {
        stop(
          "The 'design' matrix must have the same number of rows as ",
          "the pseudo-bulk count matrix has columns.\n"
        )
      }
      colnames(X_df) <- rownames(design_model)
      # extract interaction terms
      col_with_colon <- grep(":", colnames(design_model), value = TRUE)
  }

  if (!is.fullrank(design_model)) {
    stop("The design matrix is not full rank. ",
         "Please check for multicollinearity or redundant columns in your model.")
  }
  y <- DGEList(counts= X_df, genes=rownames(X_df))
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design_model, robust=TRUE)
  if(test == "QLF"){
    fit <- glmQLFit(y, design_model)
    # include all interaction terms
    lrt <- glmQLFTest(fit, coef = col_with_colon)
  }else{ # LRT
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = col_with_colon)
  }
  res_edgeR <- topTags(lrt, n = Inf)
  results <- as.data.frame(res_edgeR[[1]])
  colnames(results)[1] <- "gene_id"
  if(verbose){
      if(test == "QLF"){
          return(list(gene_results = results,
                      estimated_y = y,
                      glmQLFTest = lrt,
                      glmQLFit = fit))
      } else{ # LRT
          return(list(gene_results = results,
                      estimated_y = y,
                      glmLRT = lrt,
                      glmFit = fit))
      }
  }else{
    return(results)
  }   
}
# Function derived from https://github.com/HelenaLC/muscat/blob/main/R/aggregateData.R
.aggregateData <- function(x, 
                          assay = NULL, by = c("cluster_id", "sample_id"), 
                          fun = c("sum", "mean", "median", "prop.detected", "num.detected"), 
                          scale = FALSE, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)) {
  
  # check validity of input arguments
  fun <- match.arg(fun)
  if (is.null(assay)) 
    assay <- assayNames(x)[1] 
  stopifnot(is.character(assay), length(assay) == 1, assay %in% assayNames(x))
  if (sum(assayNames(x) == assay) > 1)
      stop("Argument 'assay' was matched to multiple times.\n ", 
           " Please assure that the input SCE has unique 'assayNames'.")
  .check_args_aggData(as.list(environment()))
  stopifnot(is(BPPARAM, "BiocParallelParam"))
  
  # assure 'by' colData columns are factors
  # so that missing combinations aren't dropped
  for (i in by) 
    if (!is.factor(x[[i]])) 
      x[[i]] <- factor(x[[i]])
  
  # compute pseudo-bulks
  pb <- .pb(x, by, assay, fun, BPPARAM)
  if (scale & length(by) == 2) {
    # compute library sizes
    cs <- if (assay == "counts" && fun == "sum")
      pb else .pb(x, by, "counts", "sum", BPPARAM)
    ls <- lapply(cs, colSums)
    # scale pseudobulks by CPM
    pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
    names(pb) <- names(ls)
  }
  
  # construct SCE
  md <- metadata(x)
  md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
  pb <- SingleCellExperiment(pb, rowData = rowData(x), metadata = md)
  
  # tabulate number of cells
  cd <- data.frame(colData(x)[, by])
  for (i in names(cd))
    if (is.factor(cd[[i]]))
      cd[[i]] <- droplevels(cd[[i]])
  ns <- table(cd)
  if (length(by) == 2) {
    ns <- asplit(ns, 2)
    ns <- lapply(ns, function(x) c(unclass(x)))
  } else ns <- c(unclass(ns))
  int_colData(pb)$n_cells <- ns
  
  # propagate 'colData' columns that are unique across 2nd 'by'
  if (length(by) == 2) {
    cd <- colData(x)
    ids <- colnames(pb)
    counts <- vapply(ids, function(u) {
      m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(colData(x))))
    cd_keep <- apply(counts, 1, function(u) all(u == 1))
    cd_keep <- setdiff(names(which(cd_keep)), by)
    if (length(cd_keep) != 0) {
      m <- match(ids, cd[, by[2]], nomatch = 0)
      cd <- cd[m, cd_keep, drop = FALSE]
      rownames(cd) <- ids
      colData(pb) <- cd
    }
  }
  return(pb)
}
.check_args_aggData <- function(u) {
  stopifnot(is.character(u$by), length(u$by) <= 2, 
            u$by %in% colnames(colData(u$x)))
  stopifnot(is.logical(u$scale), length(u$scale) == 1)
  if (u$scale & (!u$assay %in% c("cpm", "CPM") | u$fun != "sum"))
    stop("Option 'scale = TRUE' only valid for", 
         " 'assay = \"cpm/CPM\"' and 'fun = \"sum\"'.")
}
.pb <- function(x, by, assay, fun, BPPARAM = SerialParam()) {
  # compute pseudobulks
  suppressWarnings( 
    # temporarily suppressing warnings b/c 'median' 
    # warns about unspecified 'useNames' argument
    y <- summarizeAssayByGroup(x,
                               assay.type = assay, 
                               ids = (ids <- colData(x)[by]),
                               statistics = fun,
                               BPPARAM = BPPARAM))
  colnames(y) <- y[[by[length(by)]]]
  
  if (length(by) == 1) 
    return(assay(y))
  
  # reformat into one assay per 'by[1]'
  if (is.factor(ids <- y[[by[1]]]))
    ids <- droplevels(ids)
  is <- split(seq_len(ncol(y)), ids)
  ys <- lapply(is, function(idx) assay(y)[, idx, drop = FALSE])
  
  # fill in missing combinations
  unique_by2 <- unique(y[[by[2]]])
  ys <- lapply(ys, function(mat) {
      missing_cols <- setdiff(unique_by2, colnames(mat))
      if (length(missing_cols) > 0) {
          fill_matrix <- matrix(0, nrow = nrow(x), ncol = length(missing_cols),
                            dimnames = list(NULL, missing_cols))
          mat <- cbind(mat, fill_matrix)
          mat <- mat[, sort(colnames(mat)), drop = FALSE]
      }
      mat
  })
  return(ys)
}

`%notin%` <- Negate(`%in%`)
list <- structure(NA,class="result")
#' @export
"[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(seq_len(2),length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
        a <- args[[i]]
        if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
}
