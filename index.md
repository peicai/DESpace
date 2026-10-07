# DESpace: a framework to discover spatially variable genes and differential spatial patterns across conditions

![DESpace logo](reference/figures/DESpace.png)

[![Bioc
release](https://bioconductor.org/shields/build/release/bioc/DESpace.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/DESpace)
[![Bioc
devel](https://bioconductor.org/shields/build/devel/bioc/DESpace.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/DESpace)
[![Bioc
years](https://bioconductor.org/shields/years-in-bioc/DESpace.svg)](https://bioconductor.org/packages/DESpace)
[![Bioc
downloads](https://bioconductor.org/shields/downloads/release/DESpace.svg)](https://bioconductor.org/packages/stats/bioc/DESpace)

`DESpace` is a framework for identifying spatially variable genes
(SVGs), a common task in spatial transcriptomics analyses, and
differential spatial variable pattern (DSP) genes, which identify
differences in spatial gene expression patterns across experimental
conditions.

By leveraging pre-annotated spatial clusters as summarized spatial
information, `DESpace` models gene expression with a negative binomial
(NB), via
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html),
with spatial clusters as covariates. SV genes are then identified by
testing the significance of spatial clusters.

For multi-sample, multi-condition datasets, again we fit a NB model via
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html),
but this time we use spatial clusters, conditions and their interactions
as covariates. DSP genes are then identified by testing the interaction
between spatial clusters and conditions.

Check the
[vignettes](https://peicai.github.io/DESpace/articles/SVG.html) for a
description of the main conceptual and mathematical aspects, as well as
usage guidelines.

## Citation

If you use `DESpace`, please cite:

> Peiying Cai, Mark D. Robinson, and Simone Tiberi (2024). DESpace:
> spatially variable gene detection via differential expression testing
> of spatial clusters. *Bioinformatics*, 40(2), btae027.
> [doi:10.1093/bioinformatics/btae027](https://doi.org/10.1093/bioinformatics/btae027)

> Peiying Cai, Mark D. Robinson, and Simone Tiberi (2026). DESpace2:
> detection of differential spatial patterns in spatial omics data.
> *Bioinformatics*, 42(Supplement_2), btag450.
> [doi:10.1093/bioinformatics/btag450](https://doi.org/10.1093/bioinformatics/btag450)

The citation information is also available from R via
`citation("DESpace")`.

## Bioconductor installation

`DESpace` is available on
[Bioconductor](https://bioconductor.org/packages/DESpace) and can be
installed with the command:

``` r

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("DESpace")
```

## Vignette

The vignette illustrating how to use the package can be accessed on
[Bioconductor](https://bioconductor.org/packages/DESpace) or from R via:

``` r

vignette("DESpace")
```

or

``` r

browseVignettes("DESpace")
```
