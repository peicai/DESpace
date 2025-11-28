# Results from [`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md) function

Results from
[`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md)
function

## Arguments

- results_individual_svg:

  contains a `list` object, with the results obtained applying
  [`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md)
  function to an external dataset from the spatialLIBD package. Below
  the code used to obtain 'results_individual_svg'.

## Value

A List of 7 elements - one element for each spatial cluster

## See also

[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md),
[`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md)

## Author

Peiying Cai <peiying.cai@uzh.ch>, Simone Tiberi <simone.tiberi@unibo.it>

## Examples

``` r
# load the input data:
# data("LIBD_subset", package = "DESpace")
# LIBD_subset
# load pre-computed results (obtained via `svg_test`)
# data("results_svg_test", package = "DESpace")
# results_svg_test

# Function `individual_svg()` can be used to identify SVGs for each individual cluster.
# Parameter `spatial_cluster` indicates the column names of `colData(spe)` 
# containing spatial clusters.
# set.seed(123)
# results_individual_svg <- individual_svg(LIBD_subset,
#                                            edgeR_y = results_svg_test$estimated_y,
#                                            spatial_cluster = "layer_guess_reordered")
# save(results_individual_svg, file = "./DESpace/data/results_individual_svg.RData")
```
