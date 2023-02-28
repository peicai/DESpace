test_that("individual_test() works faultlessly.", {
  # Connect to ExperimentHub
  ehub <- ExperimentHub::ExperimentHub()
  # Download the example spe data
  spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)
  spe_all
  
  # Only use the sample 151507:
  spe3 <- spe_all[, colData(spe_all)$sample_id == '151673']
  # Select small set of random genes for faster runtime in this example
  set.seed(123)
  sel_genes <- sample(dim(spe3)[1],500)
  spe3 <- spe3[sel_genes,]
  
  # Fit the model via DESpace_test function.
  set.seed(123)
  results_DESpace <- DESpace_test(spe = spe3,
                                      spatial_cluster = "layer_guess_reordered",
                                      verbose = TRUE)
  
  # Individual cluster test: identify SVGs for each individual cluster
  set.seed(123)
  cluster_results <- individual_test(spe3,
                                     edgeR_y = results_DESpace$estimated_y,
                                     spatial_cluster = "layer_guess_reordered")
  
  expect_is(cluster_results, "list")
  expect_is(cluster_results[[1]], "data.frame")
  expect_is(cluster_results[[2]], "data.frame")
  expect_is(cluster_results[[3]], "data.frame")
  expect_is(cluster_results[[4]], "data.frame")
  expect_is(cluster_results[[5]], "data.frame")
  expect_is(cluster_results[[6]], "data.frame")
  expect_is(cluster_results[[7]], "data.frame")
  expect_true( length(cluster_results) == 7 ) # 7 clusters in colData(spe3)$"layer_guess_reordered"
})

