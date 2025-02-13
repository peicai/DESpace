# # load SPE object
# eh <- ExperimentHub::ExperimentHub()
# spe <- eh[["EH9613"]]
# 
# test_that("dsp_test() works faultlessly.", {
#   # Fit the model via svg_test function.
#   set.seed(123)
#   results_DESpace <- dsp_test(spe = spe,
#                               cluster_col = "Banksy_smooth",
#                               sample_col = "sample_id",
#                               condition_col = "condition",
#                               verbose = TRUE)
#   expect_is(results_DESpace, "list")
#   expect_is(results_DESpace[[1]], "data.frame")
#   expect_is(results_DESpace[[2]], "list")
#   expect_is(results_DESpace[[3]], "list")
#   expect_is(results_DESpace[[4]], "list")
# })