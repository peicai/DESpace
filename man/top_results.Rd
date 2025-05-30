% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/top_results.R
\name{top_results}
\alias{top_results}
\title{top_results}
\usage{
top_results(
  gene_results = NULL,
  cluster_results,
  cluster = NULL,
  select = "both",
  high_low = NULL
)
}
\arguments{
\item{gene_results}{Results returned from \code{\link{svg_test}}.}

\item{cluster_results}{Results returned from \code{\link{individual_svg}}.}

\item{cluster}{A character indicating the cluster(s) whose results have to be returned.
Results from all clusters are returned by default ("NULL").}

\item{select}{A character indicating what results should be returned ("FDR", "logFC", or "both").
Only used if "cluster_results" are provided.
By default ("both"), both FDR and logFC are returned.}

\item{high_low}{A character indicating whether to filter results or not.
Only used if "cluster_results" are provided, and one cluster is specified in "cluster" parameter.
By default (NULL), all results are returned in a single data.frame.
If "high" or "HIGH", we only return SVGs with average abundace in "cluster" higher than in the rest of the tissue (i.e., logFC > 0).
If "low" or "LOW", we only return SVGs with average abundace in "cluster" lower than in the rest of the tissue (i.e., logFC < 0).
If "both" or "BOTH", then both "high" and "low" results are returned, but in two separate data.frames.}
}
\value{
A \code{\linkS4class{data.frame}} object or a list of \code{\linkS4class{data.frame}} with results.

- When only “cluster_results” is provided, results are reported as a \code{\linkS4class{data.frame}} with columns for
 gene names (gene_id), spatial clusters affected by SV (Cluster), cluster-specific likelihood ratio test statistics (LR),
 cluster-specific average (across spots) log-2 counts per million (logCPM), cluster-specific log2-fold changes (logFC),
 cluster-specific raw p-values (PValue), and Benjamini-Hochberg adjusted p-values (FDR) for each spatial cluster.

- When “gene_results” and “cluster_results” are given, results are reported as a \code{\linkS4class{data.frame}} that merges gene-
 and cluster-level results.

- If “cluster” is specified, the function returns a subset \code{\linkS4class{data.frame}} for the given cluster, which contains cluster name,
 gene name, LR, logCPM, logFC, PValue and FDR, ordered by FDR for the specified cluster.

- If “high_low” is set, the function returns a list of \code{\linkS4class{data.frame}} that contains subsets of results for genes with
 higher and/or lower expression in the given cluster compared to the rest of the tissue.
}
\description{
Filter significant results.
\code{\link{top_results}} returns the significant results obtained via \code{\link{svg_test}} and \code{\link{individual_svg}}.
It can also be used to merge gene- and cluster-level results into a single object.
}
\examples{
# load pre-computed results (obtained via `svg_test`)
data("results_svg_test", package = "DESpace")

# svg_test returns of a list of 2 objects:
# "gene_results": a dataframe contains main edgeR test results;
# "estimated_y": a DGEList object contains the estimated common dispersion, 
#  which can later be used to speed-up calculation when testing individual clusters.

# We visualize differential results:
head(results_svg_test$gene_results, 3)

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
# 'select = "FDR"' can be used to visualize adjusted p-values for each spatial cluster.
merge_res = top_results(results_svg_test$gene_results, 
                        results_individual_svg, select = "FDR")
head(merge_res,3)
# Specify the cluster of interest and check top genes detected by svg_test.
results_WM_both = top_results(cluster_results = results_individual_svg, 
                              cluster = "WM", high_low = "both")
head(results_WM_both$high_genes, 3)
head(results_WM_both$low_genes, 3)

}
\seealso{
\code{\link{svg_test}}, \code{\link{individual_svg}}, \code{\link{FeaturePlot}}, \code{\link{dsp_test}}, \code{\link{individual_dsp}}
}
