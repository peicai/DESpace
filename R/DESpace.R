#' DESpace: A package for identifying spatially variable genes
#' 
#' An intuitive framework for identifying spatially variable genes (SVGs) and differential spatial pattern (DSP) genes via edgeR, 
#' one of the most common methods for performing differential expression analyses.
#'
#' Based on pre-annotated spatial clusters as summarized spatial information,
#' \code{DESpace} models gene expression using a negative binomial (NB), via \code{edgeR}, with spatial clusters as covariates.
#' SVGs are then identified by testing the significance of spatial clusters,
#' whereas DSP genes are identified by testing the significance of the interaction terms between spatial clusters and conditions (e.g., treatment conditions or time phases).
#' 
#' Our approach assumes that the spatial structure can be summarized by spatial clusters, 
#' which should reproduce the key features of the tissue (e.g., white matter and layers in brain cortex).
#' These spatial clusters are therefore taken as proxy for the actual spatial distribution; 
#' a significant test of these covariates indicates that space influences gene expression, 
#' hence identifying spatially variable genes.
#'
#' Our model is flexible and robust, and is significantly faster than the most SV methods.
#' Furthermore, to the best of our knowledge, it is the only SV approach that allows: 
#' - performing a SV test on each individual spatial cluster, 
#' hence identifying the key regions affected by spatial variability; 
#' - jointly fitting multiple samples, 
#' targeting genes with consistent spatial patterns across replicates.
#' 
#' @details
#' For an overview of the functionality provided by the package, please see the
#' vignette:
#' \code{vignette("DESpace", package="DESpace")}
#'
#' @author 
#' Peiying Cai \email{peiying.cai@uzh.ch}, 
#' Simone Tiberi \email{simone.tiberi@unibo.it}
#' @seealso \code{\link{svg_test}}, \code{\link{individual_svg}}, \code{\link{top_results}}, \code{\link{FeaturePlot}}, \code{\link{dsp_test}}, \code{\link{individual_dsp}}
#' 
#' @name DESpace
#' 
#' @keywords internal
"_PACKAGE"
