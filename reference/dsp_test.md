# dsp_test

'dsp_test' identifies differential spatial pattern (DSP) genes between
conditions from spatially-resolved transcriptomics data, provided
spatial clusters are available.

## Usage

``` r
dsp_test(
  spe,
  design = NULL,
  cluster_col,
  sample_col,
  condition_col,
  test = "QLF",
  min_counts = 20,
  min_non_zero_spots = 10,
  min_pct_cells = 0.5,
  filter_gene = FALSE,
  filter_cluster = TRUE,
  verbose = FALSE
)
```

## Arguments

- spe:

  SpatialExperiment or SingleCellExperiment.

- design:

  Matrix or array. Numeric design matrix for a regression-like model
  created by \`model.matrix\` function.

- cluster_col:

  Character. Column name of spatial clusters in `colData(spe)`.

- sample_col:

  Character. Column name of sample ids in `colData(spe)`. Sample ids
  must be either a factor or character.

- condition_col:

  Character. Column name of condition ids in `colData(spe)`.

- test:

  Character. Either 'QLF' or 'LRT'. Default is 'QLF'. Specifies whether
  to perform quasi-likelihood F-tests (\`edgeR::glmQLFTest\`) or
  likelihood ratio tests (\`edgeR::glmLRT\`). Quasi-likelihood tests
  apply empirical Bayes moderation to genewise dispersions, account for
  intra-sample variability, and provide more conservative error control.

- min_counts:

  Numeric. Minimum number of counts per sample (across all spots) for a
  gene to be analyzed.

- min_non_zero_spots:

  Numeric. Minimum number of non-zero spots per sample, for a gene to be
  analyzed.

- min_pct_cells:

  Numeric. Minimum percentage of cells required for each cluster to be
  included in the analysis across the specified conditions. Default
  value is 0.5 (i.e., 0.5% of total cells per cluster per condition).

- filter_gene:

  Logical. If TRUE, `dsp_test` filters genes by requiring them to be
  expressed in at least 'min_non_zero_spots' cells and have at least
  'min_counts' counts per sample across all locations.

- filter_cluster:

  Logical. When set to TRUE, `dsp_test` excludes clusters that are
  insufficiently represented in the dataset. Only clusters meeting the
  'min_pct_cells' threshold (i.e., containing at least the specified
  percentage of cells across all conditions) will be retained for
  analysis.

- verbose:

  Logical. If TRUE,
  [`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md)
  returns additional components containing the underlying edgeR model
  object. By default (`test = "QLF"`), these are the 'edgeR::glmQLFTest'
  and 'edgeR::glmQLFit' results. When `test = "LRT"`, the returned
  objects correspond to 'edgeR::glmFit' and 'edgeR::glmLRT'.

## Value

A list containing:

\- "gene_results": a data frame with the main edgeR test results.

\- "estimated_y": a DGEList object containing the estimated common
dispersion, which can be reused to speed up computation when testing
individual clusters.

\- "glmQLFit" (only if `verbose = TRUE`): a DGEGLM object containing the
full statistics from \`edgeR::glmQLFit\` (when `test = "QLF"`).

\- "glmQLFTest" (only if `verbose = TRUE`): a DGELRT-like object
containing the full statistics from \`edgeR::glmQLFTest\` (when
`test = "QLF"`).

\- "glmFit" (only if `verbose = TRUE`): a DGEGLM object containing the
full statistics from \`edgeR::glmFit\` (when `test = "LRT"`).

\- "glmLRT" (only if `verbose = TRUE`): a DGELRT object containing the
full statistics from \`edgeR::glmLRT\` (when `test = "LRT"`).

## See also

[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md),
[`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md),
[`individual_dsp`](https://peicai.github.io/DESpace/reference/individual_dsp.md),
[`FeaturePlot`](https://peicai.github.io/DESpace/reference/FeaturePlot.md),
[`top_results`](https://peicai.github.io/DESpace/reference/top_results.md)

## Examples

``` r
## Load the example multi-sample multi-group spe object
spe <- muSpaData::Wei22_example()
#> 
#> see ?muSpaData and browseVignettes('muSpaData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> require(“SpatialExperiment”)
# Fit the model via \code{\link{dsp_test}} function.
set.seed(123)
results_dsp <- dsp_test(spe = spe,
                        cluster_col = "Banksy_smooth",
                        sample_col = "sample_id",
                        condition_col = "condition",
                        verbose = FALSE)
#> Using 'dsp_test' for spatial variable pattern genes detection.
#> Filter low quality clusters: 
#> Cluster levels to keep: 0, 1, 2, 3, 4
#> Design model: row names represent sample names, followed by underscores and cluster names.
#>          (Intercept) condition20DPI condition2DPI cluster_id1 cluster_id2
#> 2DPI_1_0           1              0             1           0           0
#> 2DPI_2_0           1              0             1           0           0
#>          cluster_id3 cluster_id4 condition20DPI:cluster_id1
#> 2DPI_1_0           0           0                          0
#> 2DPI_2_0           0           0                          0
#>          condition2DPI:cluster_id1 condition20DPI:cluster_id2
#> 2DPI_1_0                         0                          0
#> 2DPI_2_0                         0                          0
#>          condition2DPI:cluster_id2 condition20DPI:cluster_id3
#> 2DPI_1_0                         0                          0
#> 2DPI_2_0                         0                          0
#>          condition2DPI:cluster_id3 condition20DPI:cluster_id4
#> 2DPI_1_0                         0                          0
#> 2DPI_2_0                         0                          0
#>          condition2DPI:cluster_id4
#> 2DPI_1_0                         0
#> 2DPI_2_0                         0
#> calcNormFactors has been renamed to normLibSizes

# dsp_test returns of an object:
# "gene_results": a dataframe contains main edgeR test results.

# We visualize differential results:
head(results_dsp, 3)
#>                       gene_id logFC.condition20DPI.cluster_id1
#> AMEX60DD014721 AMEX60DD014721                      -0.06935701
#> AMEX60DD045083 AMEX60DD045083                       0.09262813
#> AMEX60DD011151 AMEX60DD011151                      -0.96000279
#>                logFC.condition2DPI.cluster_id1 logFC.condition20DPI.cluster_id2
#> AMEX60DD014721                      -0.3892637                       -0.8484362
#> AMEX60DD045083                       0.5474369                       -1.3319239
#> AMEX60DD011151                      -2.0955778                        0.1933357
#>                logFC.condition2DPI.cluster_id2 logFC.condition20DPI.cluster_id3
#> AMEX60DD014721                       0.8522171                       -0.8305364
#> AMEX60DD045083                       1.1218265                        0.2418323
#> AMEX60DD011151                       0.4691729                       -0.3910272
#>                logFC.condition2DPI.cluster_id3 logFC.condition20DPI.cluster_id4
#> AMEX60DD014721                      -0.9450399                      0.224922592
#> AMEX60DD045083                       1.2749710                     -0.990433737
#> AMEX60DD011151                       0.4123802                     -0.007782088
#>                logFC.condition2DPI.cluster_id4   logCPM        F       PValue
#> AMEX60DD014721                       0.2094252 9.344907 17.05839 1.514597e-08
#> AMEX60DD045083                      -0.8071333 7.505402 13.49192 1.676387e-07
#> AMEX60DD011151                       0.1872899 7.959492 11.67885 6.888665e-07
#>                         FDR
#> AMEX60DD014721 7.572984e-05
#> AMEX60DD045083 4.190966e-04
#> AMEX60DD011151 1.097488e-03
```
