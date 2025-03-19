# load SPE object
eh <- ExperimentHub::ExperimentHub()
spe <- eh[["EH9613"]]

test_that("individual_dsp() works faultlessly.", {
  # Individual cluster test: identify SVGs for each individual cluster
  set.seed(123)
  cluster_results <- individual_dsp(spe,
                                    cluster_col = "Banksy_smooth",
                                    sample_col = "sample_id",
                                    condition_col = "condition")

  expect_is(cluster_results, "list")
  expect_is(cluster_results[[1]], "data.frame")
  expect_is(cluster_results[[2]], "data.frame")
  expect_is(cluster_results[[3]], "data.frame")
  expect_is(cluster_results[[4]], "data.frame")
  expect_is(cluster_results[[5]], "data.frame")
  expect_true( length(cluster_results) == 5 ) # 5 clusters in spe[["Banksy_smooth"]]
})