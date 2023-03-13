test_that("DESpace() works faultlessly.", {
  # load the input data:
  data("LIBD_subset", package = "DESpace")
  
  # Fit the model via DESpace_test function.
  set.seed(123)
  results_DESpace <- DESpace_test(spe = LIBD_subset,
                                      spatial_cluster = "layer_guess_reordered",
                                      verbose = TRUE)
  
  expect_is(results_DESpace, "list")
  expect_is(results_DESpace[[1]], "data.frame")
  expect_is(results_DESpace[[2]], "list")
  expect_is(results_DESpace[[3]], "list")
  expect_is(results_DESpace[[4]], "list")
})