# top_results

Filter significant results. `top_results` returns the significant
results obtained via
[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md) and
[`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md).
It can also be used to merge gene- and cluster-level results into a
single object.

## Usage

``` r
top_results(
  gene_results = NULL,
  cluster_results,
  cluster = NULL,
  select = "both",
  high_low = NULL
)
```

## Arguments

- gene_results:

  Results returned from
  [`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md).

- cluster_results:

  Results returned from
  [`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md).

- cluster:

  A character indicating the cluster(s) whose results have to be
  returned. Results from all clusters are returned by default ("NULL").

- select:

  A character indicating what results should be returned ("FDR",
  "logFC", or "both"). Only used if "cluster_results" are provided. By
  default ("both"), both FDR and logFC are returned.

- high_low:

  A character indicating whether to filter results or not. Only used if
  "cluster_results" are provided, and one cluster is specified in
  "cluster" parameter. By default (NULL), all results are returned in a
  single data.frame. If "high" or "HIGH", we only return SVGs with
  average abundace in "cluster" higher than in the rest of the tissue
  (i.e., logFC \> 0). If "low" or "LOW", we only return SVGs with
  average abundace in "cluster" lower than in the rest of the tissue
  (i.e., logFC \< 0). If "both" or "BOTH", then both "high" and "low"
  results are returned, but in two separate data.frames.

## Value

A `data.frame` object or a list of `data.frame` with results.

\- When only “cluster_results” is provided, results are reported as a
`data.frame` with columns for gene names (gene_id), spatial clusters
affected by SV (Cluster), cluster-specific likelihood ratio test
statistics (LR), cluster-specific average (across spots) log-2 counts
per million (logCPM), cluster-specific log2-fold changes (logFC),
cluster-specific raw p-values (PValue), and Benjamini-Hochberg adjusted
p-values (FDR) for each spatial cluster.

\- When “gene_results” and “cluster_results” are given, results are
reported as a `data.frame` that merges gene- and cluster-level results.

\- If “cluster” is specified, the function returns a subset `data.frame`
for the given cluster, which contains cluster name, gene name, LR,
logCPM, logFC, PValue and FDR, ordered by FDR for the specified cluster.

\- If “high_low” is set, the function returns a list of `data.frame`
that contains subsets of results for genes with higher and/or lower
expression in the given cluster compared to the rest of the tissue.

## See also

[`svg_test`](https://peicai.github.io/DESpace/reference/svg_test.md),
[`individual_svg`](https://peicai.github.io/DESpace/reference/individual_svg.md),
[`FeaturePlot`](https://peicai.github.io/DESpace/reference/FeaturePlot.md),
[`dsp_test`](https://peicai.github.io/DESpace/reference/dsp_test.md),
[`individual_dsp`](https://peicai.github.io/DESpace/reference/individual_dsp.md)

## Examples

``` r
# load pre-computed results (obtained via `svg_test`)
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

# load pre-computed results (obtained via `individual_svg`)
data("results_individual_svg", package = "DESpace")
# Function `individual_svg()` can be used to identify SVGs for each individual cluster.
# `individual_svg()` returns a list containing the results of individual clusters.
# For each cluster, results are reported as a data.frame, 
# where columns For each cluster, results are reported as a data.frame, 
# where columns contain gene names (`genes`), likelihood ratio (`LR`), 
# log2-fold changes (`logFC`) and adjusted p-value (`FDR`).
# 
# Combine gene-and cluster-level results
merge_res = top_results(results_svg_test$gene_results, 
                        results_individual_svg)
head(merge_res,3)
#>           gene_id   gene_LR gene_logCPM   gene_Pvalue      gene_FDR
#> 1 ENSG00000115756 1073.3689    15.89723 1.204871e-228 2.855545e-226
#> 2 ENSG00000054690 1056.7744    15.48088 4.686692e-225 5.553730e-223
#> 3 ENSG00000111716  752.4418    16.58993 2.894030e-159 2.286283e-157
#>     Layer1_FDR    Layer2_FDR   Layer3_FDR   Layer4_FDR   Layer5_FDR
#> 1 8.177323e-07 9.097439e-117 1.298693e-42 2.954756e-16 1.287636e-48
#> 2 2.458083e-15  2.905948e-06 6.466892e-44 8.444596e-07 9.986696e-17
#> 3 4.260477e-47  2.227331e-01 4.739143e-26 8.444596e-07 5.775841e-15
#>     Layer6_FDR        WM_FDR Layer1_logFC Layer2_logFC Layer3_logFC
#> 1 1.163469e-07  1.087300e-37    0.6095577   -1.6910482   -0.7174741
#> 2 1.960661e-01 1.856053e-212    1.4127189    0.7519805    1.1098266
#> 3 1.750846e-02  3.355203e-89    1.2238140   -0.1102450   -0.3652317
#>   Layer4_logFC Layer5_logFC Layer6_logFC  WM_logFC
#> 1    1.0572344    1.0578709    0.3909357 -1.067597
#> 2    0.8599344    0.7981348    0.1633057  2.346996
#> 3   -0.3473850   -0.3200502    0.1242379 -1.079196
# 'select = "FDR"' can be used to visualize adjusted p-values for each spatial cluster.
merge_res = top_results(results_svg_test$gene_results, 
                        results_individual_svg, select = "FDR")
head(merge_res,3)
#>           gene_id   gene_LR gene_logCPM   gene_Pvalue      gene_FDR
#> 1 ENSG00000115756 1073.3689    15.89723 1.204871e-228 2.855545e-226
#> 2 ENSG00000054690 1056.7744    15.48088 4.686692e-225 5.553730e-223
#> 3 ENSG00000111716  752.4418    16.58993 2.894030e-159 2.286283e-157
#>     Layer1_FDR    Layer2_FDR   Layer3_FDR   Layer4_FDR   Layer5_FDR
#> 1 8.177323e-07 9.097439e-117 1.298693e-42 2.954756e-16 1.287636e-48
#> 2 2.458083e-15  2.905948e-06 6.466892e-44 8.444596e-07 9.986696e-17
#> 3 4.260477e-47  2.227331e-01 4.739143e-26 8.444596e-07 5.775841e-15
#>     Layer6_FDR        WM_FDR
#> 1 1.163469e-07  1.087300e-37
#> 2 1.960661e-01 1.856053e-212
#> 3 1.750846e-02  3.355203e-89
# Specify the cluster of interest and check top genes detected by svg_test.
results_WM_both = top_results(cluster_results = results_individual_svg, 
                              cluster = "WM", high_low = "both")
head(results_WM_both$high_genes, 3)
#>                 Cluster         gene_id Cluster_LR Cluster_logCPM Cluster_logFC
#> ENSG00000054690      WM ENSG00000054690   978.6555       15.48088      2.346996
#> ENSG00000166086      WM ENSG00000166086   572.4992       15.47657      1.845724
#> ENSG00000117266      WM ENSG00000117266   376.6276       15.46008      1.588379
#>                 Cluster_PValue   Cluster_FDR
#> ENSG00000054690  7.831449e-215 1.856053e-212
#> ENSG00000166086  1.605722e-126 1.902780e-124
#> ENSG00000117266   6.748199e-84  3.998308e-82
head(results_WM_both$low_genes, 3)
#>                 Cluster         gene_id Cluster_LR Cluster_logCPM Cluster_logFC
#> ENSG00000111716      WM ENSG00000111716   409.7061       16.58993     -1.079196
#> ENSG00000063180      WM ENSG00000063180   258.4652       16.10214     -1.064478
#> ENSG00000135940      WM ENSG00000135940   242.9055       16.88763     -0.693015
#>                 Cluster_PValue  Cluster_FDR
#> ENSG00000111716   4.247092e-91 3.355203e-89
#> ENSG00000063180   3.707311e-58 1.757266e-56
#> ENSG00000135940   9.145322e-55 3.612402e-53
```
