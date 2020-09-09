# !/usr/bin/env Rscript
# 11/10/2017

#' @import SC3
#' @importFrom e1071 svm
#' @importFrom SingleCellExperiment SingleCellExperiment normcounts normcounts<- logcounts logcounts<-
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
sc3_SAFE <- function(inputTags, gene_filter, svm_num_cells, SEED){
    exp_cell_exprs <- NULL
    sc3OUTPUT <- NULL

    # cell expression
    ### For count data, it would be normalized by the total cound number and then log2 transformed
    exp_cell_exprs <- SingleCellExperiment(assays = list(counts = inputTags))
    normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags))*1000000
    logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)

    rowData(exp_cell_exprs)$feature_symbol <- rownames(exp_cell_exprs)
    exp_cell_exprs <- exp_cell_exprs[!duplicated(rowData(exp_cell_exprs)$feature_symbol), ]

    ### Estimating optimal number of clustering
    exp_cell_exprs <- sc3_estimate_k(exp_cell_exprs)
    optimal_K <- metadata(exp_cell_exprs)$sc3$k_estimation

    ### Clustering by SC3 at the optimal K
    if (ncol(inputTags) < svm_num_cells){
        #print(optimal_K)
        exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter, n_cores = 1, rand_seed = SEED)
    } else if (ncol(inputTags) >= svm_num_cells){
        ### Runing SVM
        exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter,
                            svm_max = svm_num_cells, svm_num_cells = svm_num_cells, n_cores = 1, rand_seed = SEED)
        exp_cell_exprs <- sc3_run_svm(exp_cell_exprs, ks = optimal_K)
    }

    ### Exporting SC3 results
    p_Data <- colData(exp_cell_exprs)
    col_name <- paste("sc3_", optimal_K, "_clusters", sep = '')
    sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]
    return(sc3OUTPUT)
}

#' @import cidr
cidr_SAFE <- function(inputTags, nPC.cidr, SEED){
    set.seed(SEED)

    cidrOUTPUT <- NULL

    cidrOUTPUT <- scDataConstructor(inputTags, tagType = "raw")
    cidrOUTPUT <- determineDropoutCandidates(cidrOUTPUT)
    cidrOUTPUT <- wThreshold(cidrOUTPUT)
    cidrOUTPUT <- scDissim(cidrOUTPUT)
    cidrOUTPUT <- scPCA(cidrOUTPUT)
    cidrOUTPUT <- nPC(cidrOUTPUT)

    ### Define nPC
    if(!is.null(nPC.cidr)) {
        cidrOUTPUT@nPC <- nPC.cidr
    } else {
        nPC.cidr <- cidrOUTPUT@nPC
    }

    ### Clustering by CIDR.
    # The optimal clustering number is determined automatically
    cidrOUTPUT <- scCluster(cidrOUTPUT, nPC = nPC.cidr)

    return(cidrOUTPUT)
}

#' @import Seurat
#' @import dplyr
#' @import Matrix
#' @importFrom methods .hasSlot
seurat_SAFE <- function(inputTags, nGene_filter = TRUE, low.genes, high.genes, nPC.seurat, resolution, SEED){
    seuratOUTPUT <- NULL

    # Initialize the Seurat object with the raw data (non-normalized data)
    # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
    seuratOUTPUT <- CreateSeuratObject(counts = inputTags, min.cells = 0, min.features = low.genes, project = "single-cell clustering")
    
    # Detection of variable genes across the single cells
    if (nGene_filter == TRUE){
        seuratOUTPUT <- subset(object = seuratOUTPUT, subset = nFeature_RNA > low.genes & nFeature_RNA < high.genes)
    }
    
    # Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
    seuratOUTPUT = NormalizeData(object = seuratOUTPUT, normalization.method = "LogNormalize", scale.factor = 10000)

    # Detection of variable genes across the single cells
    seuratOUTPUT <- FindVariableFeatures(object = seuratOUTPUT, selection.method = "vst", nfeatures = 2000)

    # Scale data
    all.genes <- rownames(seuratOUTPUT)
    seuratOUTPUT <- ScaleData(seuratOUTPUT, features = all.genes)

    ### Perform linear dimensional reduction
    if (nPC.seurat <= 20){
        seuratOUTPUT <- RunPCA(object = seuratOUTPUT, features = VariableFeatures(object = seuratOUTPUT), npcs = 20, seed.use = SEED, verbose = F)
        seuratOUTPUT <- FindNeighbors(seuratOUTPUT, dims = 1:20, verbose = F)
    } else {
        seuratOUTPUT <- RunPCA(object = seuratOUTPUT, features = VariableFeatures(object = seuratOUTPUT), npcs = nPC.seurat, seed.use = SEED, verbose = F)
        seuratOUTPUT <- FindNeighbors(seuratOUTPUT, dims = 1:nPC.seurat, verbose = F)
    }
    
    seuratOUTPUT <- FindClusters(seuratOUTPUT, resolution = resolution, verbose = F)
    
    seurat_output <- t(as.matrix(as.numeric(seuratOUTPUT@active.ident)))

    return(seurat_output)
}

#' @import Rtsne
#' @import ADPclust
#' @import SAVER
#' @import doParallel
#' @importFrom S4Vectors var
#' @importFrom stats kmeans
tSNE_kmeans_SAFE <- function(inputTags, saver, dimensions, perplexity, k.min, k.max, var_genes, SEED){
    input_lcpm <- NULL
    tsne_input <- NULL
    tsne_output <- NULL
    tsne_kmeansOUTPUT <- NULL
    adpOUTPUT <- NULL

    ### Data tranformation
    if (saver == TRUE){
        cl <- makeCluster(12, outfile = "")
        registerDoParallel(cl)
        inputTags_corrected <- saver(inputTags)
        stopCluster(cl)
        tsne_input <- log2(inputTags_corrected$estimate + 1)
    } else{
        ### If the input data is original count data or CPM, it would be tranformed to CPM
        input_lcpm <- log2(t(t(inputTags)/colSums(inputTags))*1000000+1)
        tsne_input <- input_lcpm
    }

    if (is.null(var_genes)){
        set.seed(SEED)

        tsne_output <- Rtsne(t(tsne_input), dims = dimensions, perplexity = perplexity, check_duplicates = FALSE)
    } else{
        se_genes = rep(NA, nrow(tsne_input))
        for (i in 1:nrow(tsne_input)){
            se_genes[i] = sqrt(var(tsne_input[i,])/length(tsne_input[i,]))
        }
        decreasing_rank = order(se_genes, decreasing = TRUE)

        set.seed(SEED)

        tsne_output <- Rtsne(t(tsne_input[decreasing_rank[1:var_genes],]), dims = dimensions, perplexity = perplexity)
    }

    ### Determining the optimal cluster number (k) and centroid by ADPclust
    adpOUTPUT <- adpclust(tsne_output$Y, htype = "amise",centroids="auto", nclust = k.min:k.max)

    ### Clustering the cells by kmeans
    tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, tsne_output$Y[adpOUTPUT$centers,], adpOUTPUT$nclust)

    return(tsne_kmeansOUTPUT)
}

#' @title Single-cell Aggregated clustering From Ensemble (SAFE)
#'
#' @description This function performs single-cell clustering using four state-of-the-art methods,
#' SC3, CIDR, Seurat and tSNE+kmeans.
#'
#' @param inputTags a G*N matrix with G genes and N cells.
#' @param mt_filter is a boolean variable that defines whether to filter outlier cells according to mitochondrial gene percentage.
#' Default is "TRUE".
#' @param mt.pattern defines the pattern of mitochondrial gene names in the data, for example, \code{mt.pattern = "^MT-"} for human and \code{mt.pattern = "^mt-"} for mouse.
#' Default is \code{mt.pattern = "^MT-"}.
#' @param mt.cutoff defines a high cutoff of mitochondrial percentage (Default is 10%) that cells having lower percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}.
#' @param SC3 a boolean variable that defines whether to cluster cells using SC3 method.
#' Default is "TRUE".
#' @param gene_filter a boolean variable defines whether to perform gene filtering before SC3 clustering, when \code{SC3 = TRUE}.
#' @param svm_num_cells, if \code{SC3 = TRUE}, then defines the mimimum number of cells above which SVM will be run.
#' @param CIDR a boolean parameter that defines whether to cluster cells using CIDR method.
#' Default is "TRUE".
#' @param nPC.cidr defines the number of principal coordinates used in CIDR clustering, when \code{CIDR = TRUE}.
#' Default value is esimated by \code{nPC} of \code{CIDR}.
#' @param Seurat is a boolean variable that defines whether to cluster cells using Seurat method.
#' Default is "TRUE".
#' @param nGene_filter is a boolean variable that defines whether to filter outlier cells according to unique gene count before Seurat clustering.
#' Default is "TRUE".
#' @param low.genes defines a low cutoff of unique gene counts (Default is 200) that cells having less than 200 genes are filtered out, when \code{nGene_filter = TRUE}.
#' @param high.genes defines a high cutoff of unique gene counts (Default is 8000) that cells having more than 8000 genes are filtered out, when \code{nGene_filter = TRUE}.
#' @param nPC.seurat defines the number of principal components used in Seurat clustering, when \code{Seurat = TRUE}.
#' Default is \code{nPC.seurat = nPC.cidr}.
#' @param resolution defines the value of resolution used in Seurat clustering, when \code{Seurat = TRUE}. Default is \code{resolution = 0.7}.
#' @param tSNE is a boolean variable that defines whether to cluster cells using t-SNE method.
#' Default is "TRUE".
#' @param saver is a boolean variable that defines whether to revise the gene expression profile in noisy and sparse single-cell RNA-seq data for downstream tSNE analysis using SAVER method.
#' Default is "FALSE".
#' @param dimensions sets the number of dimensions wanted to be retained in t-SNE step. Default is 3.
#' @param perplexity sets the perplexity parameter for t-SNE dimension reduction. Default is 30 when number of cells \code{>=200}.
#' @param tsne_min_cells defines the number of cells in input dataset below which
#' \code{tsne_min_perplexity=10} would be employed for t-SNE step. Default is 200.
#' @param tsne_min_perplexity sets the perplexity parameter of t-SNE step for small datasets (number of cells \code{<200}).
#' @param var_genes defines the number of variable genes used by t-SNE analysis, when \code{tSNE = TRUE}.
#' @param SEED sets the seed of the random number generator. Setting the seed to a fixed value can
#' produce reproducible clustering results.
#'
#' @return a matrix of indiviudal clustering results is output, where each row represents the cluster results of each method.
#'
#' @author Yuchen Yang <yangyuchensysu@gmail.com>, Ruth Huh <rhuh@live.unc.edu>, Yun Li <yunli@med.unc.edu>
#' @references Yuchen Yang, Ruth Huh, Houston Culpepper, Yuan Lin, Michael Love, Yun Li. SAFE (Single-cell Aggregated clustering From Ensemble): Cluster ensemble for single-cell RNA-seq data. 2017
#' @examples
#' # Load the example data data_SAFE
#' data("data_SAFE")
#'
#' # Zheng dataset
#' # Run individual_clustering
#' cluster.result <- individual_clustering(inputTags=data_SAFE$Zheng.expr, SEED=123)
#' @export
individual_clustering <- function(inputTags, mt_filter = TRUE, mt.pattern = "^MT-", mt.cutoff = 0.1, 
                                SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL,
                                Seurat = TRUE, nGene_filter = TRUE, low.genes = 200, high.genes = 8000, nPC.seurat = NULL, resolution = 0.7, 
                                tSNE = TRUE, saver = FALSE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL,
                                SEED = 1){

    cluster_number <- NULL
    cluster_results <- NULL
    inputTags = as.matrix(inputTags)

    # Filter out cells that have mitochondrial genes percentage over 5%
    if (mt_filter == TRUE){
        mito.genes <- grep(pattern = mt.pattern, x = rownames(x = inputTags), value = TRUE)
        percent.mito <- Matrix::colSums(inputTags[mito.genes, ])/Matrix::colSums(inputTags)
        inputTags <- inputTags[,which(percent.mito <= mt.cutoff)]
    }

    ##### SC3
    if(SC3 == TRUE){
        message("Performing SC3 clustering...")

        sc3OUTPUT <- sc3_SAFE(inputTags = inputTags, gene_filter = gene_filter,
                            svm_num_cells = svm_num_cells, SEED = SEED)
        cluster_results <- rbind(cluster_results, matrix(c(sc3OUTPUT), nrow = 1, byrow = TRUE))
        cluster_number <- c(cluster_number, max(c(sc3OUTPUT)))
    }


    ##### CIDR
    if(CIDR == TRUE){
        message("Performing CIDR clustering...")

        cidrOUTPUT <- cidr_SAFE(inputTags = inputTags, nPC.cidr = nPC.cidr, SEED = SEED)

        if(is.null(nPC.cidr)) {
            nPC.cidr <- cidrOUTPUT@nPC
        }

        cluster_results <- rbind(cluster_results, matrix(c(cidrOUTPUT@clusters), nrow = 1, byrow = TRUE))
        cluster_number <- c(cluster_number,  cidrOUTPUT@nCluster)
    }


    ##### Seurat
    if (Seurat == TRUE){
        message("Performing Seurat clustering...")

        if(is.null(nPC.seurat)) {
            nPC.seurat <- nPC.cidr
        }

        seurat_output <- seurat_SAFE(inputTags = inputTags, nGene_filter = nGene_filter, low.genes = low.genes, high.genes = high.genes, 
                                     nPC.seurat = nPC.seurat, resolution = resolution, SEED = SEED)
        cluster_results <- rbind(cluster_results, matrix(c(seurat_output), nrow = 1, byrow = TRUE))
        cluster_number <- c(cluster_number, max(!is.na(seurat_output)))
    }


    ##### tSNE+kmeans
    if(tSNE == TRUE){
        message("Performing tSNE + k-means clustering...")

        ### Dimensionality reduction by Rtsne
        if(length(inputTags[1,]) < tsne_min_cells) {
            perplexity = tsne_min_perplexity
        }

        tsne_kmeansOUTPUT <- tSNE_kmeans_SAFE(inputTags = inputTags, saver = saver, dimensions = dimensions,
                                            perplexity = perplexity, k.min = 2, k.max = max(cluster_number), var_genes = var_genes, SEED = SEED)
        cluster_results = rbind(cluster_results, matrix(c(tsne_kmeansOUTPUT$cluster), nrow = 1, byrow = TRUE))
    }

    return(cluster_results)
}
