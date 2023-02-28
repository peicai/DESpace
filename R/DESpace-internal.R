.cluster_label <-
  function(spe, 
           cluster_list, 
           spatial_cluster = "layer_guess_reordered"){
    #print(cluster)
    metadata = as.data.frame(colData(spe))
    layer_rename = c()
    layer_rename = as.character(metadata[[spatial_cluster]])
    layer_rename[layer_rename != cluster_list] <- 'Other'
    layer_rename = as.factor(layer_rename)
    one_layer = layer_rename
    return(one_layer)
  }
.getValueRes <-
  function(cluster_results = single_cluster_results,
           select = "FDR"){
    ll <- lapply(seq_len(length(cluster_results)), function(i){
      data.frame(gene_id = cluster_results[[i]]$gene_id, 
                 value  = get(select, cluster_results[[i]]))})
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
    one_layer = layers
    design_model = model.matrix(~one_layer+sample_id)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results = data.table::as.data.table(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    return(results)
  }
.layer_test1 <-
  function(layers, 
           y){
    one_layer = layers
    design_model = model.matrix(~one_layer)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results = data.table::as.data.table(results[[1]][c("genes","LR", "logCPM", "logFC", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    return(results)
  }
.layer_test2 <-
  function(layers, 
           y){
    # layers = factor(one_layer, levels  = c('Other', 'Layer3'))
    one_layer = layers #== layer
    design_model = model.matrix(~one_layer)
    limma::is.fullrank(design_model)
    rownames(design_model) <- colnames(y)
    y <- estimateDisp(y, robust=TRUE, design = design_model)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results = data.table::as.data.table(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
    colnames(results)[1] <- "gene_id"
    return(results)
  }
.layer_test2_multi <-
  function(layers, 
           sample_id, 
           y){
    # layers = factor(one_layer, levels  = c('Other', 'Layer3'))
    one_layer = layers #== layer
    design_model = model.matrix(~one_layer)
    limma::is.fullrank(design_model)
    rownames(design_model) <- colnames(y)
    y <- estimateDisp(y, robust=TRUE, design = design_model)
    fit <- glmFit(y, design_model)
    lrt <- glmLRT(fit, coef = 2)
    results <- topTags(lrt, n = Inf)
    results = data.table::as.data.table(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
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
    #if(dim(gene_results)[2] == 5){
    #   colnames(gene_results) <- c("gene_id", "gene_id", "gene_LR","gene_Pvalue","gene_FDR")
    #}else(colnames(gene_results) <- c("gene_id", "gene_LR","gene_Pvalue","gene_FDR"))
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
           spatial_cluster, 
           sample_col, 
           # num_core,
           verbose = TRUE){
    layer = as.factor(colData(spe)[[spatial_cluster]] )
    layer  <- droplevels(layer)
    spe = spe[, !is.na(layer)]
    design = data.frame(condition = factor(colData(spe)[[spatial_cluster]]),
                        sample_id = factor(colData(spe)[[sample_col]]))
    
    design$condition <- droplevels(design$condition)
    y <- DGEList(counts=assays(spe)$counts, 
                 genes=rownames(assays(spe)$counts))
    y$samples$lib.size <- colSums(y$counts)
    y <- calcNormFactors(y)
    design_model <- model.matrix(~design$condition + design$sample_id)
    rownames(design_model) <- colnames(y)
    # BPPARAM = MulticoreParam(num_core)
    y <- estimateDisp(y, design_model, robust=TRUE)
    fit <- glmFit(y, design_model)
    q = nlevels(factor(colData(spe)[[spatial_cluster]]))
    
    lrt <- glmLRT(fit, coef = seq(2,q))
    res_edgeR = topTags(lrt, n = Inf)
    results = data.table::as.data.table(res_edgeR[[1]][c("genes", "LR", "logCPM", "PValue", "FDR")])
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
           spatial_cluster,
           #num_core = 4,
           verbose = TRUE){
    spe = subset(spe,,sample_id == sample_names[i])
    if(verbose){
      list[gene_results, estimated_y, lrt, fit] = .single_edgeR_test(spe, spatial_cluster, verbose)
    }else{
      list[gene_results, estimated_y] = .single_edgeR_test(spe, spatial_cluster, verbose)
    }
    
    DT_results1 = data.table::as.data.table(gene_results)
    DT_results1[, sample := sample_names[i]]
    data.table::setorder(DT_results1, FDR, PValue)
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
           spatial_cluster, 
           # num_core,
           verbose = TRUE){
    layer = as.factor(colData(spe)[[spatial_cluster]] )
    spe = spe[, !is.na(layer)]
    design = data.frame(condition = layer)
    
    design$condition <- droplevels(design$condition)
    y <- DGEList(counts=assays(spe)$counts, 
                 genes=rownames(assays(spe)$counts))
    y$samples$lib.size <- colSums(y$counts)
    y <- calcNormFactors(y)
    design_model <- model.matrix(~design$condition)
    rownames(design_model) <- colnames(y)
    # BPPARAM = MulticoreParam(num_core)
    y <- estimateDisp(y, design_model, robust=TRUE)
    fit <- glmFit(y, design_model)
    q = nlevels(factor(colData(spe)[[spatial_cluster]]))
    
    lrt <- glmLRT(fit, coef = seq(2,q))
    res_edgeR = topTags(lrt, n = Inf)
    results = data.table::as.data.table(res_edgeR[[1]][c("genes", "LR", "logCPM", "PValue", "FDR")])
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

`%notin%` <- Negate(`%in%`)
list <- structure(NA,class="result")
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
