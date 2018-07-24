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
sc3_SAFE <- function(inputTags, datatype, gene_filter, svm_num_cells, SEED){
    exp_cell_exprs <- NULL
    sc3OUTPUT <- NULL

    # cell expression
    if (datatype == "count") {
        ### For count data, it would be normalized by the total cound number and then log2 transformed
        exp_cell_exprs <- SingleCellExperiment(assays = list(counts = inputTags))
        normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags))*1000000
        logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
        ### For CPM, FPKM, RPKM or TPM data, it would be log2 transformed
        exp_cell_exprs <- SingleCellExperiment(assays = list(normcounts = inputTags))
        logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
    }

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
cidr_SAFE <- function(inputTags, datatype, nPC.cidr, SEED){
    set.seed(SEED)

    cidrOUTPUT <- NULL

    if (datatype == "count"){
        cidrOUTPUT <- scDataConstructor(inputTags, tagType = "raw")
    }
    else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM"){
        cidrOUTPUT <- scDataConstructor(inputTags, tagType = "cpm")
    }
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
seurat_SAFE <- function(inputTags, datatype, nPC.seurat, resolution, seurat_min_cell, resolution_min, SEED){
    seuratOUTPUT <- NULL

    # Initialize the Seurat object with the raw data (non-normalized data)
    # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
    seuratOUTPUT <- CreateSeuratObject(raw.data = inputTags, min.cells = 3, min.genes = 200, project = "single-cell clustering")

    # Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
    if (datatype == "count"){
        seuratOUTPUT = NormalizeData(object = seuratOUTPUT, normalization.method = "LogNormalize", scale.factor = 10000)
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM"){
        raw.data <- GetAssayData(object = seuratOUTPUT, slot = "raw.data")
        normalized.data <- log(raw.data+1)
        colnames(x = normalized.data) <- colnames(x = raw.data)
        rownames(x = normalized.data) <- rownames(x = raw.data)
        seuratOUTPUT <- SetAssayData(object = seuratOUTPUT, assay.type = "RNA",slot = "data", new.data = normalized.data)
    }

    # Detection of variable genes across the single cells
    seuratOUTPUT = FindVariableGenes(object = seuratOUTPUT, mean.function = ExpMean, dispersion.function = LogVMR,
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

    # Regress out unwanted sources of variation
    seuratOUTPUT <- ScaleData(object = seuratOUTPUT, vars.to.regress = c("nUMI"))

    ### Perform linear dimensional reduction
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, pc.genes = seuratOUTPUT@var.genes, do.print = FALSE)

    if (length(inputTags[1,]) >= seurat_min_cell){
        ### Determine statistically significant principal components
        # NOTE: This process can take a long time for big datasets, comment out for expediency.
        # More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
        # Here we chooes the same number of PCs used in CIDR

        ### Clustering the cells by Seurat
        seuratOUTPUT <- FindClusters(object = seuratOUTPUT, reduction.type = "pca", dims.use = 1:nPC.seurat, algorithm = 3,
                                    resolution = resolution, print.output = FALSE, random.seed = SEED)
    } else {
        resolution <- resolution_min
        seuratOUTPUT <- FindClusters(object = seuratOUTPUT, reduction.type = "pca", dims.use = 1:nPC.seurat, algorithm = 3,
                                    resolution = resolution_min, print.output = FALSE, random.seed = SEED)
    }

    ### Complementing the missing data
    cells_dropout <- NULL
    num_genes <- colSums(inputTags > 0)
    cells_dropout <- names(num_genes[which(num_genes <= 200)])
    if (length(cells_dropout != 0)){
        seurat_output <- matrix(NA, ncol = ncol(inputTags), byrow = TRUE)
        colnames(seurat_output) <- colnames(inputTags)
        seurat_retained <- t(as.matrix(as.numeric(seuratOUTPUT@ident)))
        colnames(seurat_retained) <- colnames(seuratOUTPUT@data)
        for (i in 1:ncol(seurat_retained)){
            seurat_output[1,colnames(seurat_retained)[i]] <- seurat_retained[1,colnames(seurat_retained)[i]]
        }
    } else {
        seurat_output <- t(as.matrix(as.numeric(seuratOUTPUT@ident)))
    }

    return(seurat_output)
}

#' @import Rtsne
#' @import ADPclust
#' @import SAVER
#' @import doParallel
#' @importFrom S4Vectors var
#' @importFrom stats kmeans
tSNE_kmeans_SAFE <- function(inputTags, datatype, saver, dimensions, perplexity, k.min, k.max, var_genes, SEED){
    input_lcpm <- NULL
    tsne_input <- NULL
    tsne_output <- NULL
    tsne_kmeansOUTPUT <- NULL
    adpOUTPUT <- NULL

    ### Data tranformation
    if (datatype == "count") {
        if(saver == TRUE){
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
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
        ### If the input data is FPKM or RPKM, we use the transformed TPM data generated before as the input
        tsne_input <- log2(inputTags + 1)
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
#' @param datatype defines the type of data, which could be "count", "CPM", "RPKM" and "FPKM".
#' Default is "count".
#' @param mt_filter is a boolean variable that defines whether to filter outlier cells according to mitochondrial gene percentage.
#' Default is "FALSE".
#' @param low.mt defines a low cutoff of mitochondrial percentage (Default is -Inf) that cells having lower percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}.
#' @param high.mt defines a high cutoff of mitochondrial percentage (Default is 0.05) that cells having higher percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}.
#' @param nGene_filter is a boolean variable that defines whether to filter outlier cells according to unique gene count.
#' Default is "FALSE".
#' @param low.genes defines a low cutoff of unique gene counts (Default is 200) that cells having less than 200 genes are filtered out, when \code{nGene_filter = TRUE}.
#' @param high.genes defines a high cutoff of unique gene counts (Default is 2500) that cells having more than 2500 genes are filtered out, when \code{nGene_filter = TRUE}.
#' @param SC3 a boolean variable that defines whether to cluster cells using SC3 method.
#' Default is "TRUE".
#' @param gene_filter a boolean variable defines whether to perform gene filtering
#' before SC3 clustering, when \code{SC3 = TRUE}.
#' @param svm_num_cells, if \code{SC3 = TRUE}, then defines the mimimum number of cells above which SVM will be run.
#' @param CIDR a boolean parameter that defines whether to cluster cells using CIDR method.
#' Default is "TRUE".
#' @param nPC.cidr defines the number of principal coordinates used in CIDR clustering, when \code{CIDR = TRUE}.
#' Default value is esimated by \code{nPC} of \code{CIDR}.
#' @param Seurat is a boolean variable that defines whether to cluster cells using Seurat method.
#' Default is "TRUE".
#' @param nPC.seurat defines the number of principal components used in Seurat clustering, when \code{Seurat = TRUE}.
#' Default is \code{nPC.seurat = nPC.cidr}.
#' @param resolution defines the value of resolution used in Seurat clustering, when \code{Seurat = TRUE}.
#' @param seurat_min_cell defines the mimimum number of cells in input dataset below which
#' \code{resolution} is set to 1.2, when \code{Seurat = TRUE}.
#' @param resolution_min defines the resolution used in Seurat clustering for small dataset,
#' when \code{Seurat == TRUE} and cell number of input file < \code{seurat_min_cell}.
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
#'
#' # Biase dataset
#' # Run individual_clustering
#' cluster.result <- individual_clustering(inputTags = data_SAFE$Biase.expr, datatype = "FPKM", seurat_min_cell = 200, resolution_min = 1.2, tsne_min_cells = 200, tsne_min_perplexity = 10, SEED=123)
#' @export
individual_clustering <- function(inputTags, datatype = "count", mt_filter = FALSE, low.mt = -Inf, high.mt = 0.05, nGene_filter = FALSE, low.genes = 200, high.genes = 2500,
                                SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL,
                                Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, seurat_min_cell = 200, resolution_min = 1.2,
                                tSNE = TRUE, saver = FALSE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL,
                                SEED = 1){

    cluster_number <- NULL
    cluster_results <- NULL
    inputTags = as.matrix(inputTags)

    # Filter out cells that have mitochondrial genes percentage over 5%
    if (mt_filter == TRUE){
        mito.genes <- grep(pattern = "MT-", x = rownames(x = inputTags), value = TRUE)
        percent.mito <- Matrix::colSums(inputTags[mito.genes, ])/Matrix::colSums(inputTags)
        inputTags <- inputTags[,which(percent.mito >= low.mt & percent.mito <= high.mt)]
    }

    # Filter out cells that have unique gene counts over 2,500 or less than 200
    if (nGene_filter == TRUE){
        nGene <- colSums(inputTags > 0)
        inputTags <- inputTags[,which(nGene >= low.genes & nGene <= high.genes)]
    }

    ##### SC3
    if(SC3 == TRUE){
        message("Performing SC3 clustering...")

        sc3OUTPUT <- sc3_SAFE(inputTags = inputTags, datatype = datatype, gene_filter = gene_filter,
                            svm_num_cells = svm_num_cells, SEED = SEED)
        cluster_results <- rbind(cluster_results, matrix(c(sc3OUTPUT), nrow = 1, byrow = TRUE))
        cluster_number <- c(cluster_number, max(c(sc3OUTPUT)))
    }


    ##### CIDR
    if(CIDR == TRUE){
        message("Performing CIDR clustering...")

        cidrOUTPUT <- cidr_SAFE(inputTags = inputTags, datatype = datatype, nPC.cidr = nPC.cidr, SEED = SEED)

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

        seurat_output <- seurat_SAFE(inputTags = inputTags, datatype = datatype, nPC.seurat = nPC.seurat, resolution = resolution,
                                    seurat_min_cell = seurat_min_cell, resolution_min = resolution_min, SEED = SEED)
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

        tsne_kmeansOUTPUT <- tSNE_kmeans_SAFE(inputTags = inputTags, datatype = datatype, saver = saver, dimensions = dimensions,
                                            perplexity = perplexity, k.min = 2, k.max = cluster_number, var_genes = var_genes, SEED = SEED)
        cluster_results = rbind(cluster_results, matrix(c(tsne_kmeansOUTPUT$cluster), nrow = 1, byrow = TRUE))
    }

    return(cluster_results)
}
