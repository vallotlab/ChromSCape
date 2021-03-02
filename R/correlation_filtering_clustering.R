## Authors : Pacôme Prompsy, Celine Vallot
## Title : Wrappers & functions to filter and cluster
## single cell data based on correlation between cells Description : Wrappers &
## functions to filter and cluser single cell data based on correlation between
## cells

#' Correlation and hierarchical clustering
#'
#' Calculates cell to cell correlation matrix based on the PCA feature space and
#' runs hierarchical clustering taking 1 - correlation scores as distance.
#'
#' This functions takes as input a SingleCellExperiment object that must have
#' PCA calculated and outputs a SingleCellExperiment object with correlation
#' matrix and hierarchical clustering.
#'
#' @param scExp A SingleCellExperiment object, containing 'PCA' in reducedDims.
#' @param correlation A correlation method to use. See \link[stats]{hclust}.
#'   ('pearson')
#' @param hc_linkage A linkage method for hierarchical clustering. See
#'   \link[stats]{cor}. ('ward.D')
#'
#' @return Return a SingleCellExperiment object with correlation matrix &
#'   hiearchical clustering.
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom Matrix t
#' @importFrom stats cor hclust as.dist
#'
#' @examples
#' data("scExp")
#' scExp_cf = correlation_and_hierarchical_clust_scExp(scExp)
#'
correlation_and_hierarchical_clust_scExp <- function(
    scExp,correlation = "pearson", hc_linkage = "ward.D"){
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(correlation),
            is.character(hc_linkage))
    
    if (is.null(SingleCellExperiment::reducedDim(scExp, "PCA")))
        stop(
            "ChromSCape::correlation_and_hierarchical_clust_scExp -
                Run PCA on the object before correlation")
    
    pca = SingleCellExperiment::reducedDim(scExp, "PCA")
    pca_t <- Matrix::t(pca)
    cor_mat = stats::cor(pca_t, method = correlation)
    
    hc_cor = stats::hclust(stats::as.dist(1 - cor_mat), method = hc_linkage)
    hc_cor$labels = rep("", length(hc_cor$labels))
    
    scExp@metadata$hc_cor = hc_cor
    SingleCellExperiment::reducedDim(scExp, "Cor") = as.matrix(cor_mat)
    
    return(scExp)
}

#' Filter lowly correlated cells
#'
#' Remove cells that have a correlation score lower than what would be expected
#' by chance with other cells.
#'
#' This functions takes as input a SingleCellExperiment object that must have
#' correlation matrix calculated and outputs a SingleCellExperiment object
#' without lowly correlated cells. TSNE is recalculated.
#'
#' @param scExp A SingleCellExperiment object containing 'Cor', a correlation
#'   matrix, in reducedDims.
#' @param random_iter Number of random matrices to create to calculate random
#'   correlation scores. (50)
#' @param corr_threshold Quantile of random correlation score above which a cell
#'   is considered to be 'correlated' with another cell. (99)
#' @param percent_correlation Percentage of the cells that any cell must be
#'   'correlated' to in order to not be filtered. (1)
#' @param run_tsne Re-run tsne ? (FALSE)
#' @param verbose (TRUE)
#'
#' @return Returns a SingleCellExperiment object without lowly correlated cells.
#'   The calculated correlation score limit threshold is saved in metadata.
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom Matrix t
#' @importFrom Rtsne Rtsne
#' @importFrom stats cor hclust as.dist
#'
#' @examples
#' data("scExp")
#' dim(scExp)
#' scExp_cf = filter_correlated_cell_scExp(scExp,
#' corr_threshold = 99, percent_correlation = 1)
#' dim(scExp_cf)
filter_correlated_cell_scExp <- function(scExp, random_iter = 50,
    corr_threshold = 99, percent_correlation = 1, run_tsne =FALSE,
    verbose = FALSE){
    warning_filter_correlated_cell_scExp(
        scExp, random_iter,corr_threshold, percent_correlation, run_tsne,
        verbose)
    scExp@metadata$Unfiltered <- scExp
    pca_t = Matrix::t(SingleCellExperiment::reducedDim(scExp, "PCA"))
    correlation_values <- vector(length = random_iter)
    corChIP <- SingleCellExperiment::reducedDim(scExp, "Cor")
    limitC <- 0
    for (i in seq_len(random_iter)) {
        random_mat <- matrix(sample(pca_t), nrow = dim(pca_t)[1])
        threshold <-
            quantile(stats::cor(random_mat), probs = seq(0, 1, 0.01))
        limitC <- threshold[corr_threshold + 1]
        correlation_values[i] = limitC}
    limitC_mean = mean(correlation_values, na.rm = TRUE)
    selection_cor_filtered <- (apply(corChIP, 1, function(x)
            length(which(x > limitC_mean)))) >
        (percent_correlation * 0.01) * dim(corChIP)[1]
    scExp <- scExp[, selection_cor_filtered]
    tab = as.data.frame(SingleCellExperiment::reducedDim(
        scExp,"Cor"))[, selection_cor_filtered]
    SingleCellExperiment::reducedDim(scExp, "Cor") = tab
    if(run_tsne){
        scExp <- run_tsne_scExp(scExp, verbose)
    }
    config = umap::umap.defaults
    config$metric = "cosine"
    umap = umap::umap(SingleCellExperiment::reducedDim(scExp, "PCA"),
        config = config, method = "naive")
    umap = as.data.frame(umap$layout)
    colnames(umap) = c("Component_1", "Component_2")
    SingleCellExperiment::reducedDim(scExp, "UMAP") = umap
    hc_cor_cor_filtered <- stats::hclust(
        stats::as.dist(1 - SingleCellExperiment::reducedDim(
            scExp,"Cor")),method = "ward.D")
    hc_cor_cor_filtered$labels = rep("", length(hc_cor_cor_filtered$labels))
    scExp@metadata$hc_cor = hc_cor_cor_filtered
    scExp@metadata$limitC = limitC_mean #specific to filtered scExp
    return(scExp)
    }

#' warning_filter_correlated_cell_scExp
#'
#' @param scExp A SingleCellExperiment object containing 'Cor', a correlation
#'   matrix, in reducedDims.
#' @param random_iter Number of random matrices to create to calculate random
#'   correlation scores. (50)
#' @param corr_threshold Quantile of random correlation score above which a cell
#'   is considered to be 'correlated' with another cell. (99)
#' @param percent_correlation Percentage of the cells that any cell must be
#'   'correlated' to in order to not be filtered. (1)
#' @param run_tsne Re-run tsne ? (FALSE)
#' @param verbose (TRUE)
#'
#' @return Warnings or Errors if the input are not correct
warning_filter_correlated_cell_scExp <- function(
    scExp, random_iter,corr_threshold, percent_correlation, run_tsne,
    verbose){
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(random_iter),
            is.numeric(corr_threshold), is.numeric(percent_correlation))
    if (is.null(SingleCellExperiment::reducedDim(scExp, "Cor")))
        stop("ChromSCape::filter_correlated_cell_scExp - 
                No correlation, run correlation_and_hierarchical_clust_scExp")
    if (is.null(SingleCellExperiment::reducedDim(scExp, "PCA")))
        stop("ChromSCape::filter_correlated_cell_scExp - No PCA, 
                run reduced_dim before filtering.")
}

#' Run tsne on single cell experiment
#'
#' @param scExp A SingleCellExperiment Object
#' @param verbose Print ?
#'
#' @return A colored kable with the number of cells per sample for display
#'
#' @export
#'
#' @importFrom Rtsne Rtsne
#' @importFrom SingleCellExperiment reducedDim
#'
#
run_tsne_scExp <- function(scExp, verbose = FALSE){
    stopifnot(is(scExp,"SingleCellExperiment"))
    
    perp = choose_perplexity(SingleCellExperiment::reducedDim(scExp,"PCA"))

    tsne = Rtsne::Rtsne(SingleCellExperiment::reducedDim(scExp, "PCA"),
                        dims = 2, max_iter = 1000, pca = FALSE, theta = 0,
                        verbose = verbose, perplexity = perp)
    tsne = as.data.frame(tsne$Y)
    
    colnames(tsne) = c("Component_1", "Component_2")
    SingleCellExperiment::reducedDim(scExp, "TSNE") = tsne
    
    return(scExp)
}

#' Table of number of cells before correlation filtering
#'
#' @param scExp A SingleCellExperiment Object
#'
#' @return A colored kable with the number of cells per sample for display
#'
#' @export
#'
#' @importFrom dplyr bind_rows tibble left_join mutate
#' @importFrom kableExtra kable kable_styling group_rows
#' @importFrom SingleCellExperiment colData
#'
#' @examples
#' data("scExp")
#' \dontrun{num_cell_before_cor_filt_scExp(scExp)}
#'
num_cell_before_cor_filt_scExp <- function(scExp)
{
    stopifnot(is(scExp, "SingleCellExperiment"))
    table <-
        as.data.frame(table(
            as.data.frame(SummarizedExperiment::colData(scExp))$sample_id
        ))
    colnames(table) = c("Sample", "#Cells")
    rownames(table) = NULL
    
    # Retrieve sample colors from user specified colors & add to table
    colors = unique(
        as.data.frame(
            SingleCellExperiment::colData(
                scExp))[, c("sample_id","sample_id_color")])
    colors = as.vector(as.character(dplyr::left_join(
        table, colors, by = c("Sample" = "sample_id")
    )[,"sample_id_color"]))
    colors = c(col2hex(colors), "")
    
    table[, 1] = as.character(table[, 1])
    table = dplyr::bind_rows(table, dplyr::tibble("Sample" = "",
                                            "#Cells" = sum(table[,-1])))
    table %>% dplyr::mutate(
    "Sample" = cell_spec("Sample", color = "white", bold = TRUE,
                        background = colors)) %>%
        kableExtra::kable(escape = FALSE, align = "c") %>%
        kableExtra::kable_styling(
            c("striped", "condensed"), full_width = TRUE) %>% 
        kableExtra::group_rows("Total cell count", dim(table)[1], dim(table)[1])
}


#' Calculate intra correlation between cluster or samples
#'
#' @param scExp_cf A SingleCellExperiment
#' @param by On which feature to calculate correlation ("sample_id" or 
#' "cell_cluster")
#'
#' @return A data.frame of cell average intra-correlation
#' @export
#'
#' @examples
#' intra_correlation_scExp(scExp, by = "sample_id")
#' intra_correlation_scExp(scExp, by = "cell_cluster")
intra_correlation_scExp <- function(scExp_cf, by = c("sample_id",
                                                     "cell_cluster")[1]){
    stopifnot(is(scExp_cf, "SingleCellExperiment"), is.character(by))
    if (is.null(SingleCellExperiment::reducedDim(scExp, "Cor")))
        stop("ChromSCape::intra_correlation_scExp - 
                No correlation, run correlation_and_hierarchical_clust_scExp")
    if (is.null(SingleCellExperiment::reducedDim(scExp, "PCA")))
        stop("ChromSCape::intra_correlation_scExp - No PCA, 
                run reduced_dim before filtering.")
    
    # By sample
    annot = SingleCellExperiment::colData(scExp_cf)
    
    pca_t = SingleCellExperiment::reducedDim(scExp_cf,"PCA")
    cor_mat = reducedDim(scExp_cf,"Cor")
    
    intra_corr=data.frame()
    for(i in unique(annot[,by])){
        cells = as.character(annot$cell_id[which(annot[,by]==i)])
        tmp = cor_mat[cells,cells]
        tab = data.frame(tmp = rep(i,ncol(tmp)), "intra_corr" = colMeans(tmp))
        colnames(tab)[1] = by
        intra_corr=rbind(intra_corr,tab)
    }
    
    return(intra_corr)
}

#' Calculate inter correlation between cluster or samples
#'
#' @param scExp_cf A SingleCellExperiment
#' @param by On which feature to calculate correlation ("sample_id" or 
#' "cell_cluster")
#' @param reference_group Reference group to calculate correlation with. 
#' Must be in accordance with "by".
#' @param other_groups Groups on which to calculate correlation (can contain
#' multiple groups, and also reference_group). Must be in accordance with "by".
#'
#' @return A data.frame of average inter-correlation of cells in other_groups
#' with cells in reference_group
#' @export
#'
#' @examples
#' intra_correlation_scExp
inter_correlation_scExp <- function(
    scExp_cf, by = c("sample_id","cell_cluster")[1],
    reference_group = unique(scExp_cf[[by]])[1],
    other_groups = unique(scExp_cf[[by]])){
    stopifnot(is(scExp_cf, "SingleCellExperiment"), is.character(by))
    if (is.null(SingleCellExperiment::reducedDim(scExp, "Cor")))
        stop("ChromSCape::inter_correlation_scExp - 
                No correlation, run correlation_and_hierarchical_clust_scExp")
    if (is.null(SingleCellExperiment::reducedDim(scExp, "PCA")))
        stop("ChromSCape::inter_correlation_scExp - No PCA, 
                run reduced_dim before filtering.")
    
    if (! (reference_group %in% unique(scExp_cf[[by]])) )
        stop("ChromSCape::inter_correlation_scExp - Wrong reference_group.")
    if (any(!(other_groups %in% unique(scExp_cf[[by]]))))
        stop("ChromSCape::inter_correlation_scExp - Wrong reference_group.")
    
    # Select groups :
    sel = (scExp_cf[[by]] %in% c(reference_group, other_groups))
    scExp_cf = scExp_cf[,sel]
    
    # By sample
    annot = SingleCellExperiment::colData(scExp_cf)
    pca_t = SingleCellExperiment::reducedDim(scExp_cf,"PCA")
    cor_mat = reducedDim(scExp_cf,"Cor")
    
    inter_corr = data.frame()
    for(i in unique(annot[,by])){
        cells_i = as.character(annot$cell_id[which(annot[,by]==i)])
        for(j in unique(annot[,by])){
            cells_j = as.character(annot$cell_id[which(annot[,by]==j)])
            tmp = cor_mat[cells_i,cells_j]
            tab = data.frame("cells_i" = cells_i,
                             "by_i" = rep(i,nrow(tmp)),
                             "by_j" = rep(j,nrow(tmp)),
                             "inter_corr" = rowMeans(tmp))
            inter_corr=rbind(inter_corr,tab)
            
        }
    }
    colnames(inter_corr)[2:3] = c(paste0(by,"_i"), paste0(by,"_j"))

    inter_corr = inter_corr[which(
        inter_corr[,paste0(by,"_i")] %in% other_groups & 
            inter_corr[,paste0(by,"_j")] %in% reference_group),]
    rownames(inter_corr) = inter_corr$cells_i
    inter_corr$association = forcats::fct_inorder(
        paste0(inter_corr[,paste0(by,"_i")],"_",inter_corr[,paste0(by,"_j")]))
    
    return(inter_corr)
}

#' Number of cells before & after correlation filtering
#'
#' @param scExp SingleCellExperiment object before correlation filtering.
#' @param scExp_cf SingleCellExperiment object atfer correlation filtering.
#'
#' @return A colored kable with the number of cells per sample before and after
#' filtering for display
#' @export
#' @importFrom SingleCellExperiment colData
#' @importFrom kableExtra kable kable_styling group_rows cell_spec
#' @importFrom dplyr bind_rows tibble left_join mutate
#'
#' @examples
#' data("scExp")
#' scExp_cf = correlation_and_hierarchical_clust_scExp(scExp)
#' scExp_cf = filter_correlated_cell_scExp(scExp_cf,
#' corr_threshold = 99, percent_correlation = 1)
#' \dontrun{num_cell_after_cor_filt_scExp(scExp,scExp_cf)}
#'
num_cell_after_cor_filt_scExp <- function(scExp, scExp_cf)
{
    stopifnot(is(scExp, "SingleCellExperiment"),
            is(scExp_cf, "SingleCellExperiment"))
    table <-
        as.data.frame(table(
            as.data.frame(SummarizedExperiment::colData(scExp))$sample_id))
    table_filtered <-
        as.data.frame(table(
            as.data.frame(SummarizedExperiment::colData(scExp_cf))$sample_id))
    colnames(table) = c("Sample", "#Cells Before Filtering")
    rownames(table) = NULL
    colnames(table_filtered) = c("Sample", "#Cells After Filtering")
    rownames(table_filtered) = NULL
    
    colors = unique(as.data.frame(
        SummarizedExperiment::colData(scExp))[, c("sample_id",
                                                "sample_id_color")])
    colors = as.vector(as.character(dplyr::left_join(
        table, colors, by = c("Sample" = "sample_id"))[,
    "sample_id_color"]))
    colors = c(col2hex(colors), "")
    table_both = dplyr::left_join(table, table_filtered, by = c("Sample"))
    table_both[, 1] = as.character(table_both[, 1])
    table_both = dplyr::bind_rows(table_both,tibble("Sample" = "",
        "#Cells Before Filtering" = sum(table_both[,2]),
            "#Cells After Filtering" = sum(table_both[, 3])))
    table_both %>% dplyr::mutate("Sample" = kableExtra::cell_spec(
        "Sample",
        color = "white",
        bold = TRUE,
        background = colors
    )) %>% kableExtra::kable(escape = FALSE, align = "c") %>%
        kableExtra::kable_styling(
            c("striped", "condensed"), full_width = TRUE) %>%
        kableExtra::group_rows("Total cell count", dim(table_both)[1],
                            dim(table_both)[1])
}

#' Wrapper to apply ConsensusClusterPlus to scExp object
#'
#' Runs consensus hierarchical clustering on PCA feature space of scExp object.
#' Plot consensus scores for each number of clusters. See
#' \link[ConsensusClusterPlus]{ConsensusClusterPlus} - Wilkerson, M.D., Hayes,
#' D.N. (2010). ConsensusClusterPlus: a class discovery tool with confidence
#' assessments and item tracking. Bioinformatics, 2010 Jun 15;26(12):1572-3.
#'
#' This functions takes as input a SingleCellExperiment object that must have
#' 'PCA' in reducedDims and outputs a SingleCellExperiment object containing
#' consclust list calculated cluster consensus and item consensus scores in
#' metadata.
#'
#' @param scExp A SingleCellExperiment object containing 'PCA' in reducedDims.
#' @param prefix character value for output directory. Directory is created only
#'   if plot_consclust is not NULL. This title can be an abosulte or relative
#'   path.
#' @param maxK integer value. maximum cluster number to evaluate. (10)
#' @param reps integer value. number of subsamples. (100)
#' @param pItem numerical value. proportion of items to sample. (0.8)
#' @param pFeature numerical value. proportion of features to sample. (1)
#' @param distance character value. 'pearson': (1 - Pearson correlation),
#'   'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum',
#'   'canberra', 'minkowski' or custom distance function. ('pearson')
#' @param clusterAlg character value. cluster algorithm. 'hc' heirarchical
#'   (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data
#'   matrix, 'kmdist' ('hc') for k-means upon distance matrices (former km
#'   option), or a function that returns a clustering. ('hc')
#' @param innerLinkage hierarchical linkage method for subsampling. ('ward.D')
#' @param finalLinkage hierarchical linkage method for consensus matrix.
#'   ('ward.D')
#' @param plot_consclust character value. NULL - print to screen, 'pdf', 'png',
#'   'pngBMP' for bitmap png, helpful for large datasets. ('pdf')
#' @param plot_icl same as above for item consensus plot. ('png')
#'
#' @return Returns a SingleCellExperiment object containing consclust list,
#'   calculated cluster consensus and item consensus scores in metadata.
#' @export
#'
#' @references ConsensusClusterPlus package by Wilkerson, M.D., Hayes, D.N.
#'   (2010). ConsensusClusterPlus: a class discovery tool with confidence
#'   assessments and item tracking. Bioinformatics, 2010 Jun 15;26(12):1572-3.
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom Matrix t
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus calcICL
#'
#' @examples
#' data("scExp")
#' scExp_cf = correlation_and_hierarchical_clust_scExp(scExp)
#' scExp_cf = consensus_clustering_scExp(scExp)
#' 
consensus_clustering_scExp <- function(scExp, prefix = NULL, maxK = 10,
    reps = 100, pItem = 0.8, pFeature = 1, distance = "pearson",
    clusterAlg = "hc", innerLinkage = "ward.D", finalLinkage = "ward.D",
    plot_consclust = "pdf", plot_icl = "png")
{
    stopifnot(is(scExp, "SingleCellExperiment"))
    if (is.null(SingleCellExperiment::reducedDim(scExp, "PCA")))
        stop(paste0("ChromSCape::consensus_clustering_scExp - No PCA,",
            "run reduced_dim before filtering."))
    if (is.null(prefix))
    {
        plot_consclust = NULL
        plot_icl = NULL
        prefix = ""
    }
    pca_t = Matrix::t(SingleCellExperiment::reducedDim(scExp, "PCA"))
    
    consclust <- ConsensusClusterPlus::ConsensusClusterPlus(
        pca_t, maxK = maxK,  reps = reps, pItem = pItem,  
        pFeature = pFeature,  title = prefix, clusterAlg = clusterAlg,
        distance = distance, innerLinkage = innerLinkage,
        finalLinkage = finalLinkage, plot = plot_consclust)
    
    icl <- ConsensusClusterPlus::calcICL(
        consclust, plot = plot_icl, title = prefix)
    for (i in 2:maxK) {
        consclust[[i]]$consensusMatrix = NULL
        consclust[[i]]$consensusTree = NULL
        consclust[[i]]$ml = NULL
        consclust[[i]]$clrs = NULL
    }
    gc()
    scExp@metadata$consclust = consclust
    scExp@metadata$icl = icl
    return(scExp)
}

#' Choose a number of clusters
#'
#' This functions takes as input a SingleCellExperiment object with consclust
#' and a number of cluster to select. It outputs a SingleCellExperiment object
#' with each cell assigned to a correlation cluster in colData. Also calculates
#' a hierarchical clustering of the consensus associations calculated by
#' ConsensusClusterPlus.
#'
#' @param scExp A SingleCellExperiment object containing consclust in metadata.
#' @param hc_linkage A linkage method for hierarchical clustering. See
#'   \link[stats]{cor}. ('ward.D')
#' @param nclust Number of cluster to pick (3)
#' @param consensus Use consensus clustering results instead of simple
#'   hierarchical clustering ? (TRUE)
#'
#' @return Returns a SingleCellExperiment object with each cell assigned to a
#'   correlation cluster in colData.
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SummarizedExperiment colData
#' @importFrom Matrix t
#' @importFrom stats hclust as.dist
#'
#' @examples
#' data("scExp")
#' scExp_cf = correlation_and_hierarchical_clust_scExp(scExp)
#' scExp_cf = choose_cluster_scExp(scExp_cf,nclust=3,consensus=FALSE)
#' table(scExp_cf$cell_cluster)
#'
#' scExp_cf = consensus_clustering_scExp(scExp)
#' scExp_cf_consensus = choose_cluster_scExp(scExp_cf,nclust=3,consensus=TRUE)
#' table(scExp_cf_consensus$cell_cluster)
#' 
choose_cluster_scExp <- function(scExp, nclust = 3, consensus = TRUE, 
                                hc_linkage = "ward.D")
{
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(nclust),
            is.logical(consensus), is.character(hc_linkage))
    if (is.null(SingleCellExperiment::reducedDim(scExp, "PCA")))
        stop(paste0("ChromSCape::choose_cluster_scExp - No PCA, run",
                    " reduced_dim before filtering."))
    if (consensus & !"consclust" %in% names(scExp@metadata))
        stop(paste0("ChromSCape::choose_cluster_scExp - No consclust, run ",
                    "consensus_clustering_scExp before choosing cluster."))
    if (consensus & !"icl" %in% names(scExp@metadata))
        stop(paste0("ChromSCape::choose_cluster_scExp - No icl, run",
                    " consensus_clustering_scExp before choosing cluster."))
    pca_t = as.data.frame(
        Matrix::t(SingleCellExperiment::reducedDim(scExp, "PCA")))
    pca_t_ordered = pca_t[, scExp@metadata$hc_cor$order]
    if (consensus) {
        sel <- as.character(scExp$cell_id)
        cell_clusters <- 
            scExp@metadata$consclust[[nclust]]$consensusClass[sel]
    } else {
        cell_clusters = stats::cutree(scExp@metadata$hc_cor, k = nclust)
        names(cell_clusters) = SummarizedExperiment::colData(scExp)$cell_id
    }
    SummarizedExperiment::colData(scExp)[, "cell_cluster"] =
        paste("C", cell_clusters, sep = "")
    cell_clusters_list <-
        lapply(unique(cell_clusters), function(z)
            names(which(cell_clusters ==z)))
    scExp = colors_scExp(scExp = scExp, annotCol = "cell_cluster",
                        color_by = "cell_cluster", color_df = NULL)
    return(scExp)
}

#' Number of cells in each cluster
#'
#' @param scExp A SingleCellExperiment object containing chromatin groups.
#'
#' @return A formatted kable of cell assignation to each cluster.
#' @export
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom Matrix t rowSums
#' @importFrom stats chisq.test
#' @importFrom kableExtra kable kable_styling group_rows cell_spec
#'
#' @examples
#' data("scExp")
#' scExp_cf = correlation_and_hierarchical_clust_scExp(scExp)
#' scExp_cf = choose_cluster_scExp(scExp_cf,nclust=3,consensus=FALSE)
#' \dontrun{num_cell_in_cluster_scExp(scExp_cf)}
#'
num_cell_in_cluster_scExp <- function(scExp){
    stopifnot(is(scExp, "SingleCellExperiment"))
    table_raw <- as.data.frame.matrix(t(table(as.data.frame(
        SingleCellExperiment::colData(scExp))[, c("cell_cluster", "sample_id")
                                            ])))
    ord =  as.character(unique(SingleCellExperiment::colData(
        scExp)[, "sample_id"]))
    table_raw = table_raw[match(ord, rownames(table_raw)), , drop = FALSE]
    chi_pvalues = c()
    if(dim(as.matrix(table_raw))[2] > 1){

    for (i in seq_len((dim(as.matrix(table_raw))[1]))){
        contingency_tab = rbind(table_raw[i,], colSums(table_raw))
        print(contingency_tab)
        chi <- suppressWarnings(stats::chisq.test(
            x = contingency_tab, correct = FALSE))
        chi_pvalues[i] = chi$p.value}
    } else {
        chi_pvalues =  rep(1, dim(as.matrix(table_raw))[2])
    }
    tab <- table_raw
    chi_pvalues = round(chi_pvalues, 5)
    chi_pvalues[which(chi_pvalues == 0)] <- "<0.00001"
    chi_pvalues = c(chi_pvalues, "")
    colors_chromatin_group = col2hex(
        unique(SingleCellExperiment::colData(scExp)[, "cell_cluster_color"]))
    colors_sample_id = col2hex(
        unique(SingleCellExperiment::colData(scExp)[, "sample_id_color"]))
    tab <- cbind(tab, rowSums(table_raw))
    tab <- rbind(tab, `#Cells` = c(Matrix::colSums(table_raw),
                                sum(Matrix::colSums(table_raw))))
    tab <- rbind(tab, `p-value` = chi_pvalues)
    tab <- as.data.frame(
        rbind(Cluster = c(colnames(tab)[seq_along(colnames(tab))-1], ""), tab))
    tab["Cluster",] = as.character(tab["Cluster",])
    tab = cbind(rep("", nrow(tab)), tab[,(seq_len(ncol(tab)))])
    samples = kableExtra::cell_spec(
        rownames(tab)[2:(length(colors_sample_id) + 3)], color = "white",
        bold = TRUE, background = c(colors_sample_id,"black","black"))
    tab[2:(length(colors_sample_id) + 3), 1] = samples
    colnames(tab) <- NULL
    tab["Cluster",length(colors_chromatin_group) + 2] = "#Cells"
    tab["Cluster",2:(
        length(colors_chromatin_group) + 2)] = kableExtra::cell_spec(
            tab["Cluster", 2:(length(colors_chromatin_group) +2)],
            color = "white", bold = TRUE, background = c(
                colors_chromatin_group,"black"))
    rownames(tab) = NULL
    tab %>% kableExtra::kable(escape = FALSE, align = "c") %>% 
        kableExtra::kable_styling(c("striped","condensed"),
                                full_width = FALSE) %>% kableExtra::scroll_box()
    }
