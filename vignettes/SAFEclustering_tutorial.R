## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE---------------
library(knitr)
opts_chunk$set(echo = TRUE)

## ----setup for Zheng dataset, warning=FALSE-------------------------
# Setup the input expression matrix
library("SAFEclustering")
data("data_SAFE")

## ----message=FALSE---------------------------------------------------------
dim(data_SAFE$Zheng.expr)
data_SAFE$Zheng.expr[1:5, 1:5]

## ----results='hide', fig.show="hide", warning=FALSE------------------------
# Perform individual clustering
cluster.results <- individual_clustering(inputTags = data_SAFE$Zheng.expr, datatype = "count", mt_filter = FALSE, nGene_filter = FALSE, SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, tSNE = TRUE, dimensions = 3, perplexity = 30, SEED = 123)
## ----message=FALSE---------------------------------------------------------
cluster.results[1:4, 1:10]

## ----results='hide'--------------------------------------------------------
# cluster ensemble
cluster.ensemble <- SAFE(cluster_results = cluster.results, program.dir = "~/Documents/single_cell_clustering", MCLA = TRUE, CSPA = TRUE, HGPA = TRUE, SEED = 123)

## --------------------------------------------------------------------------
## [1] "HGPA partitioning at K = 2: 2 clusters at ANMI = 0.00329903476904425"
## [1] "HGPA partitioning at K = 3: 3 clusters at ANMI = 0.278691668779803"
## [1] "HGPA partitioning at K = 4: 4 clusters at ANMI = 0.00392992505505839"
## [1] "HGPA partitioning at K = 5: 5 clusters at ANMI = 0.552234460801785"
## [1] "MCLA partitioning at K = 2: 2 clusters at ANMI = 0.568294023177534"
## [1] "MCLA partitioning at K = 3: 3 clusters at ANMI = 0.929094923585274"
## [1] "MCLA partitioning at K = 4: 4 clusters at ANMI = 0.872601957447147"
## [1] "MCLA partitioning at K = 5: 4 clusters at ANMI = 0.923346490477427"
## [1] "CSPA partitioning at K = 2: 2 clusters at ANMI = 0.53144399728197"
## [1] "CSPA partitioning at K = 3: 3 clusters at ANMI = 0.850151780486274"
## [1] "CSPA partitioning at K = 4: 4 clusters at ANMI = 0.665510270422344"
## [1] "CSPA partitioning at K = 5: 5 clusters at ANMI = 0.666022118059772"
## [1] "Optimal number of clusters is 3 with ANMI = 0.929094923585274"

## ----message=FALSE---------------------------------------------------------
# ensemble results,
cluster_ensembles$Summary
cluster_ensembles$MCLA[1:10]
cluster_ensembles$MCLA_optimal_k

## ----ARI calculation-------------------------------------------------------
library(cidr)

# Cell labels of ground truth
head(data_SAFE$Zheng.celltype)

# Calculating ARI for cluster ensemble
adjustedRandIndex(cluster.ensemble$optimal_clustering, data_SAFE$Zheng.celltype)
