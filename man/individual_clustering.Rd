% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/individual_clustering.R
\name{individual_clustering}
\alias{individual_clustering}
\title{Single-cell Aggregated clustering From Ensemble (SAFE)}
\usage{
individual_clustering(inputTags, datatype = "count", mt_filter = FALSE,
  low.mt = -Inf, high.mt = 0.05, nGene_filter = FALSE, low.genes = 200,
  high.genes = 2500, SC3 = TRUE, gene_filter = FALSE,
  svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE,
  nPC.seurat = NULL, resolution = 0.9, seurat_min_cell = 200,
  resolution_min = 1.2, tSNE = TRUE, saver = FALSE, dimensions = 3,
  perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10,
  var_genes = NULL, SEED = 1)
}
\arguments{
\item{inputTags}{a G*N matrix with G genes and N cells.}

\item{datatype}{defines the type of data, which could be "count", "CPM", "RPKM" and "FPKM".
Default is "count".}

\item{mt_filter}{is a boolean variable that defines whether to filter outlier cells according to mitochondrial gene percentage.
Default is "FALSE".}

\item{low.mt}{defines a low cutoff of mitochondrial percentage (Default is -Inf) that cells having lower percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}.}

\item{high.mt}{defines a high cutoff of mitochondrial percentage (Default is 0.05) that cells having higher percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}.}

\item{nGene_filter}{is a boolean variable that defines whether to filter outlier cells according to unique gene count.
Default is "FALSE".}

\item{low.genes}{defines a low cutoff of unique gene counts (Default is 200) that cells having less than 200 genes are filtered out, when \code{nGene_filter = TRUE}.}

\item{high.genes}{defines a high cutoff of unique gene counts (Default is 2500) that cells having more than 2500 genes are filtered out, when \code{nGene_filter = TRUE}.}

\item{SC3}{a boolean variable that defines whether to cluster cells using SC3 method.
Default is "TRUE".}

\item{gene_filter}{a boolean variable defines whether to perform gene filtering
before SC3 clustering, when \code{SC3 = TRUE}.}

\item{svm_num_cells, }{if \code{SC3 = TRUE}, then defines the mimimum number of cells above which SVM will be run.}

\item{CIDR}{a boolean parameter that defines whether to cluster cells using CIDR method.
Default is "TRUE".}

\item{nPC.cidr}{defines the number of principal coordinates used in CIDR clustering, when \code{CIDR = TRUE}.
Default value is esimated by \code{nPC} of \code{CIDR}.}

\item{Seurat}{is a boolean variable that defines whether to cluster cells using Seurat method.
Default is "TRUE".}

\item{nPC.seurat}{defines the number of principal components used in Seurat clustering, when \code{Seurat = TRUE}.
Default is \code{nPC.seurat = nPC.cidr}.}

\item{resolution}{defines the value of resolution used in Seurat clustering, when \code{Seurat = TRUE}.}

\item{seurat_min_cell}{defines the mimimum number of cells in input dataset below which
\code{resolution} is set to 1.2, when \code{Seurat = TRUE}.}

\item{resolution_min}{defines the resolution used in Seurat clustering for small dataset,
when \code{Seurat == TRUE} and cell number of input file < \code{seurat_min_cell}.}

\item{tSNE}{is a boolean variable that defines whether to cluster cells using t-SNE method.
Default is "TRUE".}

\item{saver}{is a boolean variable that defines whether to revise the gene expression profile in noisy and sparse single-cell RNA-seq data for downstream tSNE analysis using SAVER method.
Default is "FALSE".}

\item{dimensions}{sets the number of dimensions wanted to be retained in t-SNE step. Default is 3.}

\item{perplexity}{sets the perplexity parameter for t-SNE dimension reduction. Default is 30 when number of cells \code{>=200}.}

\item{tsne_min_cells}{defines the number of cells in input dataset below which
\code{tsne_min_perplexity=10} would be employed for t-SNE step. Default is 200.}

\item{tsne_min_perplexity}{sets the perplexity parameter of t-SNE step for small datasets (number of cells \code{<200}).}

\item{var_genes}{defines the number of variable genes used by t-SNE analysis, when \code{tSNE = TRUE}.}

\item{SEED}{sets the seed of the random number generator. Setting the seed to a fixed value can
produce reproducible clustering results.}
}
\value{
a matrix of indiviudal clustering results is output, where each row represents the cluster results of each method.
}
\description{
This function performs single-cell clustering using four state-of-the-art methods,
SC3, CIDR, Seurat and tSNE+kmeans.
}
\examples{
# Load the example data data_SAFE
data("data_SAFE")

# Zheng dataset
# Run individual_clustering
cluster.result <- individual_clustering(inputTags=data_SAFE$Zheng.expr, SEED=123)

# Biase dataset
# Run individual_clustering
cluster.result <- individual_clustering(inputTags = data_SAFE$Biase.expr, datatype = "FPKM", seurat_min_cell = 200, resolution_min = 1.2, tsne_min_cells = 200, tsne_min_perplexity = 10, SEED=123)
}
\references{
Yuchen Yang, Ruth Huh, Houston Culpepper, Yuan Lin, Michael Love, Yun Li. SAFE (Single-cell Aggregated clustering From Ensemble): Cluster ensemble for single-cell RNA-seq data. 2017
}
\author{
Yuchen Yang <yangyuchensysu@gmail.com>, Ruth Huh <rhuh@live.unc.edu>, Yun Li <yunli@med.unc.edu>
}