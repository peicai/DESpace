# Subset from the human DLPFC 10x Genomics Visium dataset of the `spatialLIBD` package

Subset from the human DLPFC 10x Genomics Visium dataset of the
`spatialLIBD` package

## Arguments

- LIBD_subset:

  contains a
  [`SpatialExperiment-class`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  object, representing a subset of the sample 151673 from the full real
  data of the `spatialLIBD` package. Below the code used to subset the
  original dataset.

## Value

A spatial experiment object

## See also

[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md),
[`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md)

## Author

Peiying Cai <peiying.cai@uzh.ch>, Simone Tiberi <simone.tiberi@unibo.it>

## Examples

``` r
# Connect to ExperimentHub
# ehub <- ExperimentHub::ExperimentHub()
# Download the example spe data
# spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)
# Select one sample only:
# LIBD_subset <- spe_all[, colData(spe_all)$sample_id == '151673']
# Select small set of random genes for faster runtime 
# set.seed(123)
# sel_genes <- sample(dim(LIBD_subset)[1],500)
# LIBD_subset <- LIBD_subset[sel_genes,]
# keep_col <- c("array_row","array_col","layer_guess_reordered")
# library(SingleCellExperiment)
# LIBD_subset <- SpatialExperiment(assay = list(counts = assay(LIBD_subset),
#                                               logcounts = logcounts(LIBD_subset)), 
#                                  colData = colData(LIBD_subset)[keep_col])
# save(LIBD_subset, file = "./DESpace/data/LIBD_subset.RData")
```
