# DESpace: a framework to discover spatially variable genes

<img src="inst/extdata/DESpace.png" width="200" align="right"/>

`DESpace` is an intuitive framework for identifying spatially variable genes (SVGs) via [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), one of the most common methods for performing differential expression analyses. Based on pre-annotated spatial clusters as summarized spatial information, `DESpace` models gene expression using a negative binomial (NB), via [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), with spatial clusters as covariates.
SVGs are then identified by testing the significance of spatial clusters.

Check the vignettes for a description of the main conceptual and mathematical aspects, as well as usage guidelines.

> Peiying Cai, Mark D. Robinson, and Simone Tiberi (2024).
>
> DESpace: spatially variable gene detection via differential expression testing of spatial clusters.
>
> Bioinformatics.
> Available [here](https://doi.org/10.1093/bioinformatics/btae027)

## Bioconductor installation 
`DESpace` is available on [Bioconductor](https://bioconductor.org/packages/DESpace) and can be installed with the command:
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("DESpace")
```

## Vignette
The vignette illustrating how to use the package can be accessed on 
[Bioconductor](https://bioconductor.org/packages/DESpace)
or from R via:
``` r
vignette("DESpace")
```
or
``` r
browseVignettes("DESpace")
```
