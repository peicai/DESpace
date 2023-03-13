test_that("individual_test() works faultlessly.", {
  # load the input data:
  data("LIBD_subset", package = "DESpace")
  
  # load pre-computed results (obtained via `DESpace_test`)
  data("results_DESpace_test", package = "DESpace")
  edgeR_y = results_DESpace_test$estimated_y
  
  # Individual cluster test: identify SVGs for each individual cluster
  set.seed(123)
  cluster_results <- individual_test(LIBD_subset,
                                     edgeR_y = edgeR_y,
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

