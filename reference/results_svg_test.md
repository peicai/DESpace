# Results from [`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md) function

Results from
[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md)
function

## Arguments

- results_svg_test:

  contains a `list` object, with the results obtained applying
  [`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md)
  function to an external dataset from the spatialLIBD package. Below
  the code used to obtain 'results_svg_test'.

## Value

Large List of 2 elements:

\- "gene_results": a dataframe contains main edgeR test results;

\- "estimated_y": a DGEList object contains the estimated common
dispersion,

## See also

[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md)

## Author

Peiying Cai <peiying.cai@uzh.ch>, Simone Tiberi <simone.tiberi@unibo.it>

## Examples

``` r
# load the input data:
# data("LIBD_subset", package = "DESpace")
# LIBD_subset
# 
# Fit the model via `svg_test` function. 
# Parameter `spe` specifies the input `SpatialExperiment` or `SingleCellExperiment` object, 
# while `cluster_col` defines the column names of `colData(spe)` containing spatial clusters. 
# To obtain all statistics, set `verbose` to `TRUE`.
# 
# set.seed(123)
# results_svg_test <- svg_test(spe = LIBD_subset,
#                                          cluster_col = "layer_guess_reordered",
#                                          verbose = FALSE)
# 
# save(results_svg_test, file = "./DESpace/data/results_svg_test.RData")
```
