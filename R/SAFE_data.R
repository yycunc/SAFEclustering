#!/usr/bin/env Rscript

#' An expression_matrix.csv (32738 genes and 500 cells of three cell types,
#' cd56_NK cells, b cells and regulatory T cells from 10X Genomices
#' (Zheng et al. 2016)), and celltype.info.csv (List of cell type and cluster
#' information for the 500 cells).
#'
#' @docType data
#' @usage data("data_SAFE")
#' @examples
#' # Load the example data data_SAFE
#' data("data_SAFE")
#' working.dir = "~/Documents/single_cell_clustering"
#'
#' # Zheng dataset
#' # Run individual_clustering
#' cluster.result <- individual_clustering(inputTags=data_SAFE$Zheng.expr, SEED=123)
#'
#' # Cluster ensemble using SAFE clustering:
#' cluster.ensemble <- SAFE(cluster_results=cluster.result, program.dir = working.dir, SEED=123)
#'
#' # Biase dataset
#' # Run individual_clustering
#' cluster.result <- individual_clustering(inputTags = data_SAFE$Biase.expr, datatype = "FPKM", seurat_min_cell = 200, resolution_min = 1.2, tsne_min_cells = 200, tsne_min_perplexity = 10, SEED=123)
#'
#' # Cluster ensemble using SAFE clustering:
#' cluster.ensemble <- SAFE(cluster_results=cluster.result, program.dir = working.dir, SEED=123)
"data_SAFE"
