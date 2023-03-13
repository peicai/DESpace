test_that("FeaturePlot() works faultlessly.", {
  # load the input data:
  data("LIBD_subset", package = "DESpace")
  # load pre-computed results (obtained via `DESpace_test`)
  data("results_DESpace_test", package = "DESpace")
  
  feature = results_DESpace_test$gene_results$gene_id[1:3]
  library(SummarizedExperiment)
  plot_results <- FeaturePlot(LIBD_subset, feature, 
                              coordinates = c("array_row", "array_col"),
                              assay.type="counts",
                              ncol = 3, title = TRUE)
  
  expect_is(plot_results, "ggplot")
})