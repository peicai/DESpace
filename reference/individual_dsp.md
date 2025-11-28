# individual_dsp

DESpace can also be used to reveal the specific areas of the tissue
affected by DSP genes; i.e., spatial clusters that are particularly
over/under abundant compared to the average signal across conditions.
This function can be used to identify SVGs among conditions for each
individual cluster.

## Usage

``` r
individual_dsp(
  spe,
  cluster_col,
  sample_col,
  condition_col,
  test = "QLF",
  min_counts = 20,
  min_non_zero_spots = 10,
  min_pct_cells = 0.5,
  filter_gene = TRUE,
  filter_cluster = TRUE
)
```

## Arguments

- spe:

  SpatialExperiment or SingleCellExperiment.

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

  Logical. If TRUE,
  [`dsp_test`](https://peicai.github.io/DESpace/reference/dsp_test.md)
  filters genes: genes have to be expressed in at least
  'min_non_zero_spots' spots, and a gene requires at least 'min_counts'
  counts per sample (across all locations).

- filter_cluster:

  Logical. When set to TRUE,
  [`dsp_test`](https://peicai.github.io/DESpace/reference/dsp_test.md)
  excludes clusters that are insufficiently represented in the dataset.
  Only clusters meeting the 'min_pct_cells' threshold (i.e., containing
  at least the specified percentage of cells across all conditions) will
  be retained for analysis.

## Value

A list of results, with one result per spatial cluster in each element.
Specifically, each item in the list is a "gene_results" dataframe which
contains main edgeR test results.

## See also

[`top_results`](https://peicai.github.io/DESpace/reference/top_results.md),
[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md),
[`dsp_test`](https://peicai.github.io/DESpace/reference/dsp_test.md),
[`FeaturePlot`](https://peicai.github.io/DESpace/reference/FeaturePlot.md)

## Examples

``` r
# load the input data:
spe <- muSpaData::Wei22_example()
#> see ?muSpaData and browseVignettes('muSpaData') for documentation
#> loading from cache
set.seed(123)
results_individual_dsp <- individual_dsp(spe,
                                          cluster_col = "Banksy_smooth",
                                          sample_col = "sample_id",
                                          condition_col = "condition")
#> Filter low quality genes: 
#> min_counts = 20; min_non_zero_spots = 10.
#> The number of genes that pass filtering is 5000.
#> Conducting tests for layer '0' against all other layers.
#> Design model: row names represent sample names, followed by underscores and cluster names.
#>              (Intercept) condition20DPI condition2DPI cluster_id0
#> 2DPI_1_Other           1              0             1           0
#> 2DPI_2_Other           1              0             1           0
#>              condition20DPI:cluster_id0 condition2DPI:cluster_id0
#> 2DPI_1_Other                          0                         0
#> 2DPI_2_Other                          0                         0
#> Conducting tests for layer '1' against all other layers.
#> Design model: row names represent sample names, followed by underscores and cluster names.
#>              (Intercept) condition20DPI condition2DPI cluster_id1
#> 2DPI_1_Other           1              0             1           0
#> 2DPI_2_Other           1              0             1           0
#>              condition20DPI:cluster_id1 condition2DPI:cluster_id1
#> 2DPI_1_Other                          0                         0
#> 2DPI_2_Other                          0                         0
#> Conducting tests for layer '2' against all other layers.
#> Design model: row names represent sample names, followed by underscores and cluster names.
#>              (Intercept) condition20DPI condition2DPI cluster_id2
#> 2DPI_1_Other           1              0             1           0
#> 2DPI_2_Other           1              0             1           0
#>              condition20DPI:cluster_id2 condition2DPI:cluster_id2
#> 2DPI_1_Other                          0                         0
#> 2DPI_2_Other                          0                         0
#> Conducting tests for layer '3' against all other layers.
#> Design model: row names represent sample names, followed by underscores and cluster names.
#>              (Intercept) condition20DPI condition2DPI cluster_id3
#> 2DPI_1_Other           1              0             1           0
#> 2DPI_2_Other           1              0             1           0
#>              condition20DPI:cluster_id3 condition2DPI:cluster_id3
#> 2DPI_1_Other                          0                         0
#> 2DPI_2_Other                          0                         0
#> Conducting tests for layer '4' against all other layers.
#> Design model: row names represent sample names, followed by underscores and cluster names.
#>              (Intercept) condition20DPI condition2DPI cluster_id4
#> 2DPI_1_Other           1              0             1           0
#> 2DPI_2_Other           1              0             1           0
#>              condition20DPI:cluster_id4 condition2DPI:cluster_id4
#> 2DPI_1_Other                          0                         0
#> 2DPI_2_Other                          0                         0
#> Returning results
                                           
# We visualize results for the cluster '3'
results <- results_individual_dsp[['3']]
head(results,3)
#>                       gene_id logFC.condition20DPI.cluster_id3
#> AMEX60DD046788 AMEX60DD046788                       -1.1700137
#> AMEX60DD002984 AMEX60DD002984                        1.6467572
#> AMEX60DD005921 AMEX60DD005921                       -0.4105442
#>                logFC.condition2DPI.cluster_id3   logCPM        F       PValue
#> AMEX60DD046788                       0.2544327 5.623569 23.11317 5.933891e-05
#> AMEX60DD002984                       1.2951010 7.473809 22.94955 6.144564e-05
#> AMEX60DD005921                       1.7236189 5.262614 22.35051 7.045976e-05
#>                       FDR
#> AMEX60DD046788 0.09002682
#> AMEX60DD002984 0.09002682
#> AMEX60DD005921 0.09002682
```
