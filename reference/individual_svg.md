# individual_svg

DESpace can also be used to reveal the specific areas of the tissue
affected by SVGs; i.e., spatial clusters that are particularly
over/under abundant compared to the average signal. This function can be
used to identify SVGs for each individual cluster.

## Usage

``` r
individual_svg(
  spe,
  cluster_col,
  sample_col = "sample_id",
  edgeR_y = NULL,
  min_counts = 20,
  min_non_zero_spots = 10,
  filter_gene = TRUE,
  replicates = FALSE,
  BPPARAM = NULL
)
```

## Arguments

- spe:

  SpatialExperiment or SingleCellExperiment.

- cluster_col:

  Column name of spatial clusters in `colData(spe)`.

- sample_col:

  Column name of sample ids in `colData(spe)`.

- edgeR_y:

  Pre-estimated dispersion; if it's null, compute dispersion.

- min_counts:

  Minimum number of counts per sample (across all spots) for a gene to
  be analyzed.

- min_non_zero_spots:

  Minimum number of non-zero spots per sample, for a gene to be
  analyzed.

- filter_gene:

  A logical. If TRUE,
  [`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md)
  filters genes: genes have to be expressed in at least
  'min_non_zero_spots' spots, and a gene requires at least 'min counts'
  counts per sample (across all locations).

- replicates:

  Single sample or multi-sample test.

- BPPARAM:

  An optional parameter passed internally to bplapply. We suggest using
  as many cores as the number of spatial clusters. If unspecified, the
  script does not run in parallel. Note that parallel coding performs
  better only when dispersion estimations are not provided beforehand.
  Moreover, parallelizing the script will increase the memory
  requirement; if memory is an issue, leave 'BPPARAM' unspecified and,
  hence, avoid parallelization.

## Value

A list of results, with one result per spatial cluster in each element.
Specifically, each item in the list is a "gene_results" dataframe which
contains main edgeR test results.

## Details

For every spatial cluster we test, `edgeR` would normally re-compute the
dispersion estimates based on the specific design of the test. However,
this calculation represents the majority of the overall computing time.
Therefore, to speed-up calculations, we propose to use the dispersion
estimates which were previously computed for the gene-level tests. This
introduces a minor approximation which, in our benchmarks, does not lead
to decreased accuracy. If you want to use pre-computed gene-level
dispersion estimates, set `edgeR_y` to 'estimated_y'. Alternatively, if
you want to re-compute dispersion estimates (significantly slower, but
marginally more accurate option), leave edgeR_y empty.

## See also

[`top_results`](https://peicai.github.io/DESpace/reference/top_results.md),
[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md),
[`FeaturePlot`](https://peicai.github.io/DESpace/reference/FeaturePlot.md)

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

# load pre-computed results (obtaines via `svg_test`)
data("results_svg_test", package = "DESpace")

# svg_test returns of a list of 2 objects:
# "gene_results": a dataframe contains main edgeR test results;
# "estimated_y": a DGEList object contains the estimated common dispersion, 
#  which can later be used to speed-up calculation when testing individual clusters.

# We visualize differential results:
head(results_svg_test$gene_results, 3)
#>                         gene_id        LR   logCPM        PValue           FDR
#> ENSG00000115756 ENSG00000115756 1073.3689 15.89723 1.204871e-228 2.855545e-226
#> ENSG00000054690 ENSG00000054690 1056.7744 15.48088 4.686692e-225 5.553730e-223
#> ENSG00000111716 ENSG00000111716  752.4418 16.58993 2.894030e-159 2.286283e-157

# Individual cluster test: identify SVGs for each individual cluster
# set parallel computing; we suggest using as many cores as the number of spatial clusters.
# Note that parallelizing the script will increase the memory requirement;
# if memory is an issue, leave 'BPPARAM' unspecified and, hence, avoid parallelization.
set.seed(123)
results_individual_svg <- individual_svg(LIBD_subset, 
                                         edgeR_y = results_svg_test$estimated_y, 
                                         cluster_col = "layer_guess_reordered")
#> Filter low quality genes: 
#> min_counts = 20; min_non_zero_spots = 10.
#> The number of genes that pass filtering is237.
#> Pre-processing
#> Start modeling
#> Returning results
                                           
# We visualize results for the cluster 'WM'
results_WM <- results_individual_svg[[7]]
head(results_WM,3)
#>                         gene_id       LR   logCPM     logFC        PValue
#> ENSG00000054690 ENSG00000054690 978.6555 15.48088  2.346996 7.831449e-215
#> ENSG00000166086 ENSG00000166086 572.4992 15.47657  1.845724 1.605722e-126
#> ENSG00000111716 ENSG00000111716 409.7061 16.58993 -1.079196  4.247092e-91
#>                           FDR
#> ENSG00000054690 1.856053e-212
#> ENSG00000166086 1.902780e-124
#> ENSG00000111716  3.355203e-89
```
