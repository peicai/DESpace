test_that("DESpace() works faultlessly.", {
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
  
  expect_is(results_DESpace, "list")
  expect_is(results_DESpace[[1]], "data.frame")
  expect_is(results_DESpace[[2]], "list")
  expect_is(results_DESpace[[3]], "list")
  expect_is(results_DESpace[[4]], "list")
})