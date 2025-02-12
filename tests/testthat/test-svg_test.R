test_that("svg_test() works faultlessly.", {
  # load the input data:
  data("LIBD_subset", package = "DESpace")
  
  # Fit the model via svg_test function.
  set.seed(123)
  results_svg <- svg_test(spe = LIBD_subset,
                                    cluster_col = "layer_guess_reordered",
                                    verbose = TRUE)
  
  expect_is(results_svg, "list")
  expect_is(results_svg[[1]], "data.frame")
  expect_is(results_svg[[2]], "list")
  expect_is(results_svg[[3]], "list")
  expect_is(results_svg[[4]], "list")
})