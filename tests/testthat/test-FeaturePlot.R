test_that("FeaturePlot() works faultlessly.", {
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
  
  # load pre-computed results (obtained via `DESpace_test`)
  data("results_DESpace_test", package = "DESpace")
  feature = results_DESpace_test$gene_results$gene_id[1:3]
  library(SummarizedExperiment)
  plot_results <- FeaturePlot(spe3, feature, coordinates = c("array_row", "array_col"),
                              ncol = 3, title = TRUE)
  
  expect_is(plot_results, "ggplot")
})