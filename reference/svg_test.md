# svg_test

'svg_test' identifies spatially variable genes (SVGs) from
spatially-resolved transcriptomics data, provided spatial clusters are
available.

## Usage

``` r
svg_test(
  spe,
  cluster_col,
  sample_col = NULL,
  replicates = FALSE,
  min_counts = 20,
  min_non_zero_spots = 10,
  filter_gene = TRUE,
  verbose = FALSE
)
```

## Arguments

- spe:

  SpatialExperiment or SingleCellExperiment.

- cluster_col:

  Column name of spatial clusters in `colData(spe)`.

- sample_col:

  Column name of sample ids in `colData(spe)`.

- replicates:

  A logical, indicating whether biological replicates are provided
  (TRUE) or not (FALSE). If biological replicates are provided,
  `svg_test` performs a joint test across all replicates, searching for
  SVGs with consistent spatial patterns across samples.

- min_counts:

  Minimum number of counts per sample (across all spots) for a gene to
  be analyzed.

- min_non_zero_spots:

  Minimum number of non-zero spots per sample, for a gene to be
  analyzed.

- filter_gene:

  A logical. If TRUE, `svg_test` filters genes: genes have to be
  expressed in at least 'min_non_zero_spots' spots, and a gene requires
  at least 'min counts' counts per sample (across all locations).

- verbose:

  A logical. If TRUE, `svg_test` returns two more results: 'DGEGLM' and
  'DGELRT' objects contain full statistics from 'edgeR::glmFit' and
  'edgeR::glmLRT'.

## Value

A list of results:

\- "gene_results": a dataframe contains main edgeR test results;

\- "estimated_y": a DGEList object contains the estimated common
dispersion, which can later be used to speed-up calculation when testing
individual clusters.

\- "glmFit" (only if `verbose = TRUE`): a DGEGLM object contains full
statistics from "edgeR::glmFit".

\- "glmLRT" (only if `verbose = TRUE`): a DGELRT object contains full
statistics from "edgeR::glmLRT".

## Details

If 'sample_col' is not specified and 'replicates == FALSE', `svg_test`
assumed that data comes from an individual sample, and performs SV
testing on it.

If 'sample_col' is provided and 'replicates == FALSE', `svg_test` tests
each sample individually and returns a list of results for each sample.

If 'sample_col' is provided and 'replicates == TRUE', `svg_test`
performs a joint multi-sample test.

## See also

[`top_results`](https://peicai.github.io/DESpace/reference/top_results.md),
[`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md),
[`FeaturePlot`](https://peicai.github.io/DESpace/reference/FeaturePlot.md),
[`dsp_test`](https://peicai.github.io/DESpace/reference/dsp_test.md),
[`individual_dsp`](https://peicai.github.io/DESpace/reference/individual_dsp.md)

## Examples

``` r
# load the input data:
data("LIBD_subset", package = "DESpace")
LIBD_subset
#> class: SpatialExperiment 
#> dim: 500 3639 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(500): ENSG00000282097 ENSG00000225129 ... ENSG00000267199
#>   ENSG00000273129
#> rowData names(0):
#> colnames(3639): AAACAAGTATCTCCCA-1 AAACAATCTACTAGCA-1 ...
#>   TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
#> colData names(4): array_row array_col layer_guess_reordered sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(0) :
#> imgData names(0):

# Fit the model via \code{\link{svg_test}} function.
set.seed(123)
results_svg_test <- svg_test(spe = LIBD_subset,
                             cluster_col = "layer_guess_reordered",
                             verbose = FALSE)
#> using 'svg_test' for spatial gene/pattern detection.
#> Filter low quality genes: 
#> min_counts = 20; min_non_zero_spots = 10.
#> The number of genes that pass filtering is 237.
#> single sample test

# svg_test returns of a list of 2 objects:
# "gene_results": a dataframe contains main edgeR test results;
# "estimated_y": a DGEList object contains the estimated common dispersion, 
#  which can later be used to speed-up calculation when testing individual clusters.

# We visualize differential results:
head(results_svg_test$gene_results, 3)
#>                         gene_id        LR   logCPM        PValue           FDR
#> ENSG00000115756 ENSG00000115756 1064.0552 15.89723 1.246863e-226 2.955065e-224
#> ENSG00000054690 ENSG00000054690 1053.5156 15.48088 2.375929e-224 2.815476e-222
#> ENSG00000111716 ENSG00000111716  752.3436 16.58993 3.038878e-159 2.400714e-157
```
