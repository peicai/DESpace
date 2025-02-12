test_that("top_results() works faultlessly.", {
  # load pre-computed results (obtained via `svg_test`)
  data("results_svg_test", package = "DESpace")
  # load pre-computed results (obtained via `individual_svg`)
  data("results_individual_svg", package = "DESpace")
  
  # Combine gene-and cluster-level results
  merge_res <- top_results(results_svg_test$gene_results, results_individual_svg)
  expect_is(merge_res, "data.frame")
  
  merge_res <- top_results(results_svg_test$gene_results, results_individual_svg, select = "FDR")
  expect_is(merge_res, "data.frame")
  
  results_WM_both <- top_results(cluster_results = results_individual_svg, cluster = "WM", high_low = "both")
  expect_is(results_WM_both, "list")
  expect_true( length(results_WM_both) == 2 ) 
})