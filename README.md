# SAFEclustering
SAFE (Single-cell Aggregated clustering From Ensemble): Cluster ensemble for single-cell RNA-seq data

Although several methods have been recently developed for clustering cell types using single-cell RNA-seq (scRNA-Seq) data, they utilize different characteristics of data and yield varying results in terms of both the number of clusters and actual cluster assignments. Here, we present SAFE-clustering, Single-cell Aggregated (From Ensemble) clustering, a flexible, accurate and robust method for clustering scRNA-Seq data. SAFE-clustering takes as input, results from multiple clustering methods, to build one consensus solution. SAFE-clustering currently embeds four state-of-the-art methods, SC3, CIDR, Seurat and t-SNE + *k*-means; and ensembles solutions from these four methods using three hypergraph-based partitioning algorithms.

SAFEclustering is maintained by Yuchen Yang [yyuchen@email.unc.edu] and Yun Li [yunli@med.unc.edu].

## News and Updates
July 24, 2018
* Version 0.99.0 released
  + First offical release
  + Now it can only work on Mac and Linux platform

Dec 5, 2018
* Version 1.00.1 released
  + Fixing an error in Seurat clustering to allow more than 20 PCs computed
  + Fixing an error in tSNE + *k*-means clustering when specifying the maximum value of the pool of cluster numbers
  
## Installation
You can install SAFEclustering from github with:
```{r install}
install.packages("devtools")

devtools::install_github("yycunc/SAFEclustering")
```
Note that hypergraph partitioning algorithm (HGPA) is performed using the *shmetis* program (from the hMETIS package v. 1.5 (Karypis *et al.*, IEEE Transactions on Very Large Scale Integration (VLSI) Systems, 1999)), and meta-clustering algorithm (MCLA) and cluster-based similarity partitioning algorithm (CSPA) are performed using *gpmetis* program (from METIS v. 5.1.0 (Karypis and Kumar, SIAM Journal on Scientific Computing, 1998)). Please download the two programs corresponding to the operating systems you are using and put them in the working directory or provide the directory where these two programs are.

## SAFEclustering Examples
Here we provide examples using two datasets: one from Zheng *et al.*, (Nature Communications, 2016) and the other from Biase *et al.*, (Genome Research, 2014). Zheng dataset contains 500 human peripheral blood mononuclear cells (PBMCs) sequenced using GemCode platform, which consists of three cell types, CD56+ natural killer cells, CD19+ B cells and CD4+/CD25+ regulatory T cells. The original data can be downloaded from [10X GENOMICS website](https://support.10xgenomics.com/single-cell-gene-expression/datasets). The Biase dataset has 49 mouse embryo cells, which were sequenced by SMART-Seq and can be found at [NCBI GEO:GSE57249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57249).

### Load the data
```{r setup for Zheng dataset}
library("SAFEclustering")
data("data_SAFE")
```

### Zheng dataset
#### Setup the input expression matrix
```{r setup for Zheng dataset}
dim(data_SAFE$Zheng.expr)

data_SAFE$Zheng.expr[1:5, 1:5]
```

#### Perform individual clustering
Here we perform single-cell clustering using four popular methods, SC3, CIDR, Seurat and t-SNE + *k*-means, without filtering any genes or cells.

```{r individual clustering for Baron_human4 dataset, results='hide', fig.show="hide", warning=FALSE}
cluster.results <- individual_clustering(inputTags = data_SAFE$Zheng.expr, datatype = "count", 
mt_filter = FALSE, nGene_filter = FALSE, SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, 
nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, tSNE = TRUE, dimensions = 3, 
perplexity = 30, SEED = 123)
```

The function *indiviual_clustering* will output a matrix, where each row represents the cluster results of each method, and each colunm represents a cell. User can also extend SAFE-clustering to other scRNA-seq clustering methods, by putting all clustering results into a M * N matrix with M clustering methods and N cells.

```{r, message=FALSE}
cluster.results[1:4, 1:10]
```

#### Cluster ensemble

Using the clustering results generated in last step, we perform cluster ensemble using three partitioning algorithms meta-clustering algorithm (MCLA), hypergraph partitioning algorithm (HGPA) and cluster-based similarity partitioning algorithm (CSPA) (Strehl and Ghosh, Proceedings of AAAI 2002, Edmonto, Canada, 2002).

Here, the programs required, *shmetis* and *gpmetis*, are in the local working directory "~/Documents/single_cell_clustering".

```{r cluster ensemble for Baron_human4 dataset, results='hide'}
cluster.ensemble <- SAFE(cluster_results = cluster.results, program.dir = "~/Documents/single_cell_clustering", 
MCLA = TRUE, CSPA = TRUE, HGPA = TRUE, SEED = 123)
```

Here is the list of ANMI results for esemble solution of each K and each partitioning algorithm.

```{r}
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
```

Function *SAFE* will output a list for Average Normalized Mutual Information (ANMI) metric (Strehl and Ghosh Proceedings of AAAI 2002, Edmonto, Canada, 2002) between each ensemble solution and the individual solutions. The optimal clustering ensemble is selected from the ensemble solution with the highest ANMI value. 

```{r ensemble results for Baron_human4 dataset, message=FALSE}
cluster.ensemble$Summary

cluster.ensemble$MCLA[1:10]

cluster.ensemble$MCLA_optimal_k
```

We can compare the clustering results to the true labels using the Adjusted Rand Index (ARI)

```{r ARI calculation for Baron_human4 dataset}
library(cidr)

# Cell labels of ground truth
head(data_SAFE$Zheng.celltype)

# Calculating ARI for cluster ensemble
adjustedRandIndex(cluster.ensemble$optimal_clustering, data_SAFE$Zheng.celltype)
```

### Biase dataset

#### Setup the input expression matrix
```{r setup for Biase dataset}
dim(data_SAFE$Biase.expr.expr)

data_SAFE$Biase.expr[1:5, 1:5]
```

#### Perform individual clustering

Here we perform single-cell clustering using four popular methods, SC3, CIDR, Seurat and t-SNE + *k*-means, without filtering any genes or cells. Since there are only 49 cells in Biase dataset, the resolution parameter is set to 1.2 according to our benchmarking results.

```{r individual clustering for Biase dataset, results='hide', fig.show="hide", warning=FALSE}
cluster.results <- individual_clustering(inputTags = data_SAFE$Biase.expr, datatype = "FPKM",  
mt_filter = FALSE, nGene_filter = FALSE, SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, 
nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, seurat_min_cell = 200, resolution_min = 1.2, 
tSNE = TRUE, dimensions = 3, tsne_min_cells = 200, tsne_min_perplexity = 10, SEED = 123)
```

#### Cluster ensemble

Using the clustering results, we perform cluster ensemble using all the three partitioning algorithms MCLA, HGPA and CSPA.

```{r cluster ensemble for Biase dataset, results='hide', message=FALSE}
cluster.ensemble <- SAFE(cluster_results = cluster.results, program.dir = "~/Documents/single_cell_clustering", 
MCLA = TRUE, CSPA = TRUE, HGPA = TRUE, SEED = 123)
```

Here is the list of ANMI results for esemble solution of each K and each partitioning algorithm.

```{r}
## [1] "HGPA partitioning at K = 2: 2 clusters at ANMI = 0.156896209024547"
## [1] "HGPA partitioning at K = 3: 3 clusters at ANMI = 0.59768416631598"
## [1] "HGPA partitioning at K = 4: 4 clusters at ANMI = 0.614176459706577"
## [1] "MCLA partitioning at K = 2: 2 clusters at ANMI = 0.784102309002763"
## [1] "MCLA partitioning at K = 3: 3 clusters at ANMI = 0.970539568368452"
## [1] "MCLA partitioning at K = 4: 4 clusters at ANMI = 0.971666531448806"
## [1] "CSPA partitioning at K = 2: 2 clusters at ANMI = 0.601004834004939"
## [1] "CSPA partitioning at K = 3: 3 clusters at ANMI = 0.622097187347639"
## [1] "CSPA partitioning at K = 4: 4 clusters at ANMI = 0.590251500201678"
## [1] "Optimal number of clusters is 4 with ANMI = 0.971666531448806"
```

```{r ensemble results for Biase dataset, message=FALSE}
cluster.ensemble$Summary

cluster.ensemble$MCLA[1:10]

cluster.ensemble$MCLA_optimal_k
```

## Citation
Yang, Y., Huh, R., Culpepper, H., Lin, Y., Love, M., Li, Y. (2019) SAFE-clustering: Single-cell Aggregated (From Ensemble) Clustering for Single-cell RNA-seq Data. *Bioinformatics*, 35: 1269-1277. [PMID: 30202935].

## Credits
Algorithms of MCLA, HGPA and CSPA are from Strehl and Ghosh (2002) . Some codes are paraphrased from the Matlab package of ClusterEnsemble-V2.0 (http://strehl.com/soft.html).
