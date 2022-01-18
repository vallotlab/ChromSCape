## Authors : Pacôme Prompsy, Celine Vallot
##Title : Wrappers & function to create variety of plot
## to uncover heterogeneity in single cell Dataset

#' Plotting distribution of signal
#'
#' @param scExp A SingleCellExperiment Object
#' @param raw Use raw counts ?
#' @param log10 Transform using log10 ?
#' @param pseudo_counts Pseudo-count to add if using log10 
#' @param bins Number of bins in the histogram
#'
#' @return A ggplot histogram representing the distribution of count per cell
#' @export
#'
#' 
#' @importFrom  ggplot2 ggplot
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment assayNames
#' 
#' @examples
#' data("scExp")
#' plot_distribution_scExp(scExp)
plot_distribution_scExp <- function(
    scExp, raw = TRUE, log10 = FALSE, pseudo_counts = 1, bins = 150)
{
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(pseudo_counts))
    if (!raw %in% c(TRUE, FALSE) | !log10 %in% c(TRUE, FALSE)) 
        stop(paste0("ChromSCape::plot_distribution_scExp - raw and log10 must ",
                    "be true or false."))
    if (raw == FALSE && !("normcounts" %in%
                          SummarizedExperiment::assayNames(scExp))) 
        stop(paste0("ChromSCape::plot_distribution_scExp - If raw is false, ",
                    "normcounts must not be empty - run normalize_scExp first."))
    
    if (raw) 
        cell_cov_df = data.frame(
            "coverageByCell" = Matrix::colSums(counts(scExp)))
    else cell_cov_df = data.frame(
        "coverageByCell" = Matrix::colSums(normcounts(scExp)))
    
    if (log10) 
        cell_cov_df$coverageByCell = log10(
            cell_cov_df$coverageByCell + pseudo_counts)
    
    ggplot(cell_cov_df, aes(x = .data$coverageByCell)) + 
        geom_histogram(color = "black", fill = "steelblue", bins = bins) +
        labs(x = "read coverageByCell per cell") + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = 
                  element_blank(), panel.background = element_blank(),
              axis.line = element_line(colour = "black"), 
              panel.border = element_rect(colour = "black", fill = NA))
}

#' Adding colors to cells & features
#'
#' @param scExp A SingleCellExperiment Object
#' @param annotCol Column names to color
#' @param color_by If specifying color_df, column names to color
#' @param color_df Color data.frame to specify which color for which condition
#'
#' @return A SingleCellExperiment with additionnal "color" columns in colData
#' @export
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment colData
#' 
#' @examples
#' data("scExp")
#' scExp = colors_scExp(scExp,annotCol = c("sample_id",
#' "total_counts"),
#'  color_by =  c("sample_id","total_counts"))
#' 
#' #Specific colors using a manually created data.frame :
#' color_df = data.frame(sample_id=unique(scExp$sample_id),
#'  sample_id_color=c("red","blue","green","yellow"))
#' scExp = colors_scExp(scExp,annotCol="sample_id",
#' color_by="sample_id",color_df=color_df)
#' 
colors_scExp <- function(
    scExp, annotCol = "sample_id", color_by = "sample_id", color_df = NULL)
{
    stopifnot(is(scExp, "SingleCellExperiment"),
              is.character(annotCol), is.character(color_by))
    annot = as.data.frame(SingleCellExperiment::colData(scExp))
    anocol <- annotToCol2(annotS = annot[, annotCol, drop = FALSE],
                          annotT = annot, 
                          plotLegend = FALSE, categCol = NULL)
    SummarizedExperiment::colData(
        scExp)[, paste0(annotCol, "_color")] = as.data.frame(anocol, 
                                                             stringsAsFactors = FALSE)[, annotCol]  # factor or not ?
    if (!is.null(color_df))
    {
        # add custom colors
        if (!color_by %in% colnames(color_df)) 
            stop(paste0("ChromSCape::colors_scExp - color_by must be present ",
                        "in colnames of color_df is not null."))
        
        SummarizedExperiment::colData(scExp)[, paste0(annotCol, "_color")] = 
            color_df[match(SingleCellExperiment::colData(scExp)[, color_by],
                           color_df[, color_by]), paste0(color_by, "_color"),
                     drop = FALSE]
    }
    
    # same colors for cell_cluster & cluster + 'feature'
    if(any(grepl("cluster_", annotCol))){
        SummarizedExperiment::colData(scExp)[, paste0("cell_cluster_color")] =
            SummarizedExperiment::colData(scExp)[, paste0("cluster_",mainExpName(scExp),"_color")]
    }
    if(any(grepl("cell_cluster", annotCol))){
        SummarizedExperiment::colData(scExp)[, paste0("cluster_",mainExpName(scExp),"_color")] =
            SummarizedExperiment::colData(scExp)[, "cell_cluster_color"]
    }
    
    
    return(scExp)
}

#' Get color dataframe from shiny::colorInput
#'
#' @param input Shiny input object
#' @param levels_selected Names of the features
#' @param color_by Which feature color to retrieve
#' @param input_id_prefix Prefix in front of the feature names
#'
#' @return A data.frame with the feature levels and the colors of each level of
#'   this feature.
#'
#' @importFrom tibble rownames_to_column
#'   
get_color_dataframe_from_input <- function(
    input, levels_selected, color_by = c("sample_id", "total_counts"),
    input_id_prefix = "color_")
{
    stopifnot(!is.null(input), is.character(levels_selected),
              is.character(color_by))
    
    color_list <- paste0(
        "list(", paste0(levels_selected, " = input$", input_id_prefix, 
                        levels_selected, collapse = ", "), ")")
    
    color_list <- eval(parse(text = color_list))
    # Transform into dataframe with right column names
    color_df = as.matrix(color_list) %>% as.data.frame(
        stringsAsFactors = FALSE) %>% tibble::rownames_to_column(color_by)
    color_df[, 2] = as.character(color_df[, 2])
    colnames(color_df)[2] = paste0(color_by, "_color")
    
    return(color_df)
}

# Wrapper for plotting PCA & TSNE & UMAP
#' Plot reduced dimensions (PCA, TSNE, UMAP)
#'
#' @param scExp A SingleCellExperiment Object
#' @param color_by Character of eature used for coloration. Can be cell 
#' metadata ('total_counts', 'sample_id', ...) or a gene name.
#' @param reduced_dim Reduced Dimension used for plotting
#' @param select_x Which variable to select for x axis
#' @param select_y Which variable to select for y axis
#' @param downsample Number of cells to downsample
#' @param transparency Alpha parameter, between 0 and 1
#' @param max_distanceToTSS The maximum distance to TSS to consider a gene 
#' linked to a region. Used only if "color_by" is a gene name.
#' @param size Size of the points.
#' @param annotate_clusters A logical indicating if clusters should be labelled.
#' The 'cell_cluster' column should be present in metadata.
#' @param min_quantile The lower threshold to remove outlier cells, 
#' as quantile of cell embeddings (between 0 and 0.5).
#' @param max_quantile The upper threshold to remove outlier cells, 
#' as quantile of cell embeddings (between 0.5 and 1).
#'  
#' @return A ggplot geom_point plot of reduced dimension 2D reprensentation 
#' @export
#'
#' @importFrom SingleCellExperiment colData reducedDim normcounts
#' @importFrom ggplot2 ggplot geom_point labs theme element_blank element_line
#' element_rect scale_color_gradientn geom_label aes
#' @importFrom dplyr group_by filter slice_min mutate 
#' @importFrom tidyr separate_rows
#' @importFrom tibble rownames_to_column
#' @importFrom colorRamps matlab.like
#' 
#' @examples
#' data("scExp")
#' plot_reduced_dim_scExp(scExp, color_by = "sample_id")
#' plot_reduced_dim_scExp(scExp, color_by = "total_counts")
#' plot_reduced_dim_scExp(scExp, reduced_dim = "UMAP")
#' plot_reduced_dim_scExp(scExp, color_by = "CD52",  reduced_dim = "UMAP")
#' 
plot_reduced_dim_scExp <- function(
    scExp, color_by = "sample_id", reduced_dim = c("PCA", "TSNE", "UMAP"),
    select_x = NULL, select_y = NULL, downsample = 5000, 
    transparency = 0.6,  size = 1, max_distanceToTSS = 1000,
    annotate_clusters = "cell_cluster" %in% colnames(colData(scExp)),
    min_quantile = 0.01, max_quantile = 0.99)
{
    warning_plot_reduced_dim_scExp(scExp, color_by , reduced_dim,
                                   downsample, transparency,
                                   size, max_distanceToTSS,
                                   annotate_clusters,
                                   min_quantile, max_quantile)
    
    if(ncol(scExp) > downsample) scExp = scExp[,sample(ncol(scExp),
                                                       downsample,
                                                       replace = FALSE)]
    annot = SingleCellExperiment::colData(scExp)
    if(!color_by %in% colnames(annot)) {
        annot_feature = as.data.frame(SummarizedExperiment::rowRanges(scExp))
        annot_feature = annot_feature %>% 
            tidyr::separate_rows(.data[["Gene"]], sep = ", ") %>% 
            dplyr::group_by(.data[["Gene"]]) %>%
            dplyr::slice_min(.data[["distanceToTSS"]])
        annot_feature = annot_feature %>%
            dplyr::filter( .data[["distanceToTSS"]] < max_distanceToTSS)
        
        if(!color_by %in% annot_feature$Gene) stop(
            "ChromSCape::plot_reduced_dim_scExp - The gene chosen, ",
            color_by, ", is not closer than ", max_distanceToTSS, " to any", 
            "loci. Consider increasing max_distanceToTSS to take in account ",
            "this gene.")
        
        counts = SingleCellExperiment::normcounts(scExp[annot_feature$ID[which(
            annot_feature$Gene == color_by)],])
        if(!is.null(counts) & nrow(counts) > 1 ){
            counts = Matrix::rowSums(counts)
        }
        annot[,color_by] = as.numeric(counts)

    } else if(!paste0(color_by,"_color") %in% colnames(annot)){
        scExp = colors_scExp(scExp, annotCol = color_by)
        annot = SingleCellExperiment::colData(scExp)
    }
    
    plot_df = as.data.frame(
        cbind(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]), 
              annot))
    if(is.null(select_x)) select_x = colnames(plot_df)[1]
    if(is.null(select_y)) select_y = colnames(plot_df)[2]
    
    plot_df = plot_df %>% dplyr::mutate("transparency" = transparency) %>%
        dplyr::mutate("transparency" = base::I(.data[["transparency"]]))
    plot_df = plot_df %>% 
        dplyr::filter(.data[[select_x]] > quantile(plot_df[,select_x], min_quantile),
                      .data[[select_x]] < quantile(plot_df[,select_x], max_quantile)) 
    plot_df = plot_df %>% 
        dplyr::filter(.data[[select_y]] > quantile(plot_df[,select_y], min_quantile),
                      .data[[select_y]] < quantile(plot_df[,select_y], max_quantile)) 
    annot = annot[match(rownames(plot_df), rownames(annot)),]

    
    if(!color_by %in% colnames(SingleCellExperiment::colData(scExp))){
        plot_df = plot_df %>% 
            dplyr::mutate("transparency"= ifelse(
                annot[,color_by] == min(annot[,color_by]),0.15,1)) %>%
            dplyr::mutate("transparency" = I(.data[["transparency"]]))
    }
        
    p <- ggplot(plot_df, aes_string(x = select_x, y = select_y)) + 
        geom_point(size = size, aes(alpha = plot_df[,"transparency"],
                                    color = annot[, color_by])) +
        labs(color = color_by) + 
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA))
    
    if (is.numeric(annot[,color_by]))
    {
        p <- p + scale_color_gradientn(
            colours = c("grey75",rev(viridis::inferno(100)))) #+
            #guides(alpha=FALSE)
    } else
    {
        cols = unique(as.character(annot[,paste0(color_by, "_color")]))
        names(cols) = unique(as.character(annot[,color_by]))
        p <- p + scale_color_manual(values = cols) #+ guides(alpha=FALSE)
    }
    
    if(annotate_clusters){
        centroids = as.data.frame(t(sapply(unique(annot$cell_cluster), function(clust){
            cells = annot$cell_id[which(annot$cell_cluster == clust)]
            centroid = c(mean(plot_df[cells, select_x]),
                         mean(plot_df[cells, select_y]))
            return(centroid)
        }))) %>% tibble::rownames_to_column("cluster")
        p <- p + geom_label(data=centroids, aes( x=V1, y=V2, label=cluster),
                       size=4, label.padding = unit(0.15, "lines"),
                       label.size = 0.2)  
    }
    return(p)
}

#' A warning helper for plot_reduced_dim_scExp
#'
#' @param scExp A SingleCellExperiment Object
#' @param color_by Feature used for coloration
#' @param reduced_dim Reduced Dimension used for plotting
#' @param downsample Number of cells to downsample
#' @param transparency Alpha parameter, between 0 and 1
#' @param size Size of the points.
#' @param max_distanceToTSS Numeric. Maximum distance to a gene's TSS to consider
#' a region linked to a gene. 
#' @param annotate_clusters A logical indicating if clusters should be labelled.
#' The 'cell_cluster' column should be present in metadata.
#' @param min_quantile The lower threshold to remove outlier cells, 
#' as quantile of cell embeddings (between 0 and 0.5).
#' @param max_quantile The upper threshold to remove outlier cells, 
#' as quantile of cell embeddings (between 0.5 and 1).
#' 
#' @return Warning or errors if the inputs are not correct
#' 
warning_plot_reduced_dim_scExp <- function(scExp, color_by , reduced_dim,
                                           downsample,
                                           transparency, size, max_distanceToTSS,
                                           annotate_clusters,
                                           min_quantile, max_quantile){
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(color_by),
              is.character(reduced_dim), is.numeric(downsample),
              is.numeric(transparency), is.numeric(size),
              is.logical(annotate_clusters), is.numeric(max_distanceToTSS),
              is.numeric(min_quantile), is.numeric(max_quantile))
    if (!reduced_dim[1] %in% SingleCellExperiment::reducedDimNames(scExp)) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - ", reduced_dim[1],
                    " is not present in object, please run normalize_scExp ",
                    "first."))
    if (!color_by %in% colnames(SingleCellExperiment::colData(scExp))) {
        genes = unique(
            unlist(strsplit(SummarizedExperiment::rowRanges(scExp)$Gene,
                                split = ", ", fixed=TRUE)))
        if(!color_by %in% genes){
            stop("ChromSCape::plot_reduced_dim_scExp - color_by must be a ",
        "character vector present in colnames of colData(scExp) or be a gene ",
                        "name present in the rowRanges(scExp) after calling ",
                        "feature_annotation_scExp() function.")
        } 
    }
}

#' Plot Top/Bottom most contributing features to PCA
#'
#' @details 
#' If a gene TSS is within 10,000bp of the region, the name of the gene(s) will
#' be displayed instead of the region
#' 
#' @param scExp A SingleCellExperiment containing "PCA" in reducedDims and gene
#' annotation in rowRanges
#' @param component The name of the component of interest 
#' @param n_top_bot An integer number of top and bot regions to plot 
#'
#' @return A barplot of top and bottom features with the largest absolute
#' value in the component of interest
#' 
#' @importFrom SingleCellExperiment reducedDim normcounts
#' @importFrom SummarizedExperiment rowRanges
#'  
#' @export
#'
#' @examples
#' data(scExp)
#' plot_most_contributing_features(scExp, component = "Component_1")
plot_most_contributing_features <- function(scExp, component = "Component_1",
                                            n_top_bot = 10){
    
    top_bot  = retrieve_top_bot_features_pca(
        SingleCellExperiment::reducedDim(scExp,"PCA"),
        SingleCellExperiment::normcounts(scExp), component, n_top_bot)
    
    gene_info = SummarizedExperiment::rowRanges(scExp)
    gene_info = gene_info[match(rownames(top_bot), gene_info$ID)]
    Gene = substr(gene_info$Gene,start = 1, stop = 13)
    Gene = ifelse(nchar(Gene)>=13, paste0(Gene,"..."), Gene)
    top_bot$genes = ifelse(gene_info$distanceToTSS < 10000,
                           Gene,  gene_info$ID)
   
    top_bot$genes = factor(top_bot$genes, levels = unique(top_bot$genes))
    
    ggplot(top_bot) + geom_bar(aes(x = genes, y = .data[[component]]),
                               fill = c(rep("#52C461DF",n_top_bot),
                                        rep("red",n_top_bot)),
                               alpha = 0.7, stat="identity") + theme_classic() +
        theme(axis.text.x = element_text(angle=75, hjust = 1)) +
        xlab("") + ylab(paste0("Contribution of features to ", component))
}

#' Pie chart of top contribution of chromosomes in the 100 most contributing 
#' features to PCA
#'#' 
#' @param scExp A SingleCellExperiment containing "PCA" in reducedDims and gene
#' annotation in rowRanges
#' @param component The name of the component of interest 
#' @param n_top_bot An integer number of top and bot regions to plot (100)
#'
#' @return A pie chart showing the distribution of chromosomes in the top 
#' features with the largest absolute value in the component of interest
#' 
#' @importFrom SingleCellExperiment reducedDim normcounts
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom graphics pie
#'  
#' @export
#'
#' @examples
#' data(scExp)
#'  plot_pie_most_contributing_chr(scExp, component = "Component_1")
plot_pie_most_contributing_chr <- function(scExp, component = "Component_1",
                                           n_top_bot = 100){
    
    top_bot  = retrieve_top_bot_features_pca(
        SingleCellExperiment::reducedDim(scExp,"PCA"),
        SingleCellExperiment::normcounts(scExp), component, n_top_bot,
        absolute = TRUE)
    
    chr_info = as.data.frame(SummarizedExperiment::rowRanges(scExp))
    chr_info = chr_info[match(rownames(top_bot), chr_info$ID),]
    chr_info$absolute_value = abs(top_bot[,component])
    distrib = chr_info %>% dplyr::group_by(seqnames) %>% 
        summarise(contribution = sum(absolute_value))
    distrib$contribution = 100 * distrib$contribution/sum(distrib$contribution)
    distrib = setNames(distrib$contribution,distrib$seqnames)
    distrib = sort(distrib, decreasing = TRUE)
    distrib[6] = sum(distrib[6:length(distrib)])
    names(distrib)[6] = "Others"
    distrib = distrib[1:6]
    distrib_pc = round(distrib, 1)
    names(distrib) = paste0(names(distrib), " - ",distrib_pc, " %")
    graphics::pie(distrib,clockwise = TRUE, radius = 1,
        col = c("#5DA5DA","#FAA43A","#60BD68","#F17CB0","#DECF3F","#4D4D4D"))
}

#' Retrieve Top and Bot most contributing features of PCA
#'
#' @param pca A matrix/data.frame of rotated data
#' @param counts the normalized counts used for PCA
#' @param component the componenent of interest
#' @param n_top_bot the number of top & bot features to take
#' @param absolute If TRUE, return the top features in absolute values instead.
#'
#' @return a data.frame of top bot contributing features in PCA 
#'
retrieve_top_bot_features_pca <- function(pca, counts, component, n_top_bot,
                                          absolute = FALSE){
    rotations = counts %*% as.matrix(pca)

    rotations = rotations[,component,drop=FALSE]
    top = as.data.frame(as.matrix(
        utils::head(rotations[order(rotations, decreasing = TRUE),, drop=FALSE],
             n_top_bot)))
    bot = as.data.frame(as.matrix(
        utils::tail(rotations[order(rotations, decreasing = TRUE),, drop=FALSE],
             n_top_bot)))
    top_bot = as.data.frame(rbind(top, bot))
    if(absolute) {
        rotations = abs(rotations)
        top_bot = as.data.frame(as.matrix(
            utils::head(rotations[order(rotations, decreasing = TRUE),, drop=FALSE],
             n_top_bot)))
    }
    return(top_bot)
}

#' Plot cell correlation heatmap with annotations
#'
#' @param scExp A SingleCellExperiment Object
#' @param name_hc Name of the hclust contained in the SingleCellExperiment
#'   object
#' @param corColors A palette of colors for the heatmap
#' @param color_by Which features to add as additional bands on top of plot
#' @param downsample Number of cells to downsample
#' @param hc_linkage A linkage method for hierarchical clustering. See
#'   \link[stats]{cor}. ('ward.D')
#'
#' @return A heatmap of cell to cell correlation, grouping cells by hierarchical
#'   clustering.
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim reducedDimNames colData
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' data("scExp")
#' plot_heatmap_scExp(scExp)
#' 
plot_heatmap_scExp <- function(scExp, name_hc = "hc_cor", corColors = (
    grDevices::colorRampPalette(c("royalblue", "white", "indianred1")))(256),
    color_by = NULL, downsample = 1000, hc_linkage = "ward.D")
{
    
    stopifnot(is(scExp, "SingleCellExperiment"))
    if (!"Cor" %in% SingleCellExperiment::reducedDimNames(scExp)) 
        stop(paste0("ChromSCape::plot_heatmap_scExp - No correlation, run ",
                    "correlation_and_hierarchical_clust_scExp before filtering."))
    
    if (!name_hc %in% names(scExp@metadata)) 
        stop(paste0("ChromSCape::plot_heatmap_scExp - No dendrogram, run ",
                    "correlation_and_hierarchical_clust_scExp before filtering."))
    
    if (length(scExp@metadata[[name_hc]]$order) != ncol(scExp)) 
        stop(paste0("ChromSCape::plot_heatmap_scExp - Dendrogram has different",
                    " number of cells than dataset."))
    if(ncol(scExp) > downsample) {
        samp = sample(ncol(scExp),downsample, replace = FALSE)
        scExp = scExp[,samp]
        SingleCellExperiment::reducedDim(scExp, "Cor") = 
            SingleCellExperiment::reducedDim(scExp, "Cor")[,samp]
        cor_mat = as.matrix(SingleCellExperiment::reducedDim(scExp, "Cor"))
        hc_cor = stats::hclust(as_dist(1 - cor_mat), method = hc_linkage)
        hc_cor$labels = rep("",ncol(scExp))
        scExp@metadata[[name_hc]] = hc_cor
    }
    
    anocol = as.matrix(
        SingleCellExperiment::colData(
            scExp)[, grep("_color",
                          colnames(SingleCellExperiment::colData(scExp))), 
                   drop = FALSE])
    colnames(anocol) = gsub("_color","",colnames(anocol))
    
    if(!is.null(color_by)) {
        if(length(intersect(color_by,colnames(anocol)))>0) 
            anocol = anocol[,color_by,drop=FALSE]
    }
    
    return(
        hclustAnnotHeatmapPlot(
            x = as.matrix(SingleCellExperiment::reducedDim(scExp, "Cor"))[
                scExp@metadata[[name_hc]]$order, 
                scExp@metadata[[name_hc]]$order], hc = scExp@metadata[[name_hc]],
            hmColors = corColors,
            anocol = anocol[scExp@metadata[[name_hc]]$order, ,drop=FALSE],
            xpos = c(0.15, 0.9, 0.164, 0.885),
            ypos = c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95), dendro.cex = 0.04, 
            xlab.cex = 0.8, hmRowNames = FALSE)
    )
}


#' Violin plot of intra-correlation distribution
#'
#' @param scExp_cf A SingleCellExperiment
#' @param by Color by sample_id or cell_cluster
#' @param jitter_by Add jitter points of another layer
#'  (cell_cluster or sample_id)
#' @param downsample Downsample for plotting
#'
#' @return A violin plot of intra-correlation
#' @export
#' @importFrom forcats fct_inorder
#' @examples 
#' data(scExp)
#' plot_intra_correlation_scExp(scExp)
plot_intra_correlation_scExp <- function(
    scExp_cf, by = c("sample_id", "cell_cluster")[1], jitter_by = NULL,
    downsample = 5000){
    
    if(ncol(scExp_cf) > downsample) scExp_cf = scExp_cf[,sample(ncol(scExp_cf),
                                                       downsample,
                                                       replace = FALSE)]
    
    samp = intra_correlation_scExp(scExp_cf, by)
    annot = SingleCellExperiment::colData(scExp_cf)
    if(!is.null(jitter_by)){
        if(jitter_by != by) samp[,jitter_by] = 
                annot[,jitter_by][match(rownames(samp),annot$cell_id)]
    }
    p <- ggplot(samp) + geom_violin(aes(x=.data[[by]], y=.data$intra_corr,
                                        fill=.data[[by]]), alpha=0.8) +
        theme_classic() +
        theme(axis.text.x = element_text(angle=90)) + 
        ylab(paste0("Intra-",by," correlation")) + xlab("")
    
    cols = unique(as.character(annot[,paste0(by,"_color")][
        match(rownames(samp),annot$cell_id)]))
    names(cols) = unique(as.character(annot[,by][match(
        rownames(samp),annot$cell_id)]))
    p <- p + scale_fill_manual(values = cols)
    
    if(!is.null(jitter_by)){
        if(grepl("counts", jitter_by)){
            p <- p + geom_jitter(
                aes(x=.data[[by]], y=.data$intra_corr,
                    color = .data[[jitter_by]]),
                alpha = 0.45) +
                scale_color_gradientn(colours = matlab.like(100))
            
        } else{
            p <- p + geom_jitter(
                aes(x=.data[[by]], y=.data$intra_corr,
                    color = forcats::fct_inorder(.data[[jitter_by]])),
                alpha = 0.45) +
                scale_color_manual(
                    values = unique(annot[,paste0(jitter_by,"_color")][
                        match(rownames(samp),annot$cell_id)]))
        }
    }
    return(p + theme(legend.title = element_text("")))
}

#' Violin plot of inter-correlation distribution between one or multiple groups
#' and one reference group
#'
#' @param scExp_cf A SingleCellExperiment
#' @param by Color by sample_id or cell_cluster
#' @param jitter_by Add jitter points of another layer
#'  (cell_cluster or sample_id)
#' @param reference_group Character containing the reference group name to 
#' calculate correlation from.
#' @param other_groups Character vector of the other groups for which to 
#' calculate correlation with the reference group.
#' @param downsample Downsample for plotting
#'
#' @return A violin plot of inter-correlation
#' @export
#'
#' @importFrom forcats fct_inorder
#' 
#' @examples
#' data(scExp)
#' plot_intra_correlation_scExp(scExp)
plot_inter_correlation_scExp <- function(
    scExp_cf, by = c("sample_id", "cell_cluster")[1], jitter_by = NULL,
    reference_group = unique(scExp_cf[[by]])[1],
    other_groups = unique(scExp_cf[[by]]),
    downsample = 5000){
    
    if(ncol(scExp_cf) > downsample) scExp_cf =
            scExp_cf[,sample(ncol(scExp_cf),downsample,replace = FALSE)]
    
    samp = inter_correlation_scExp(scExp_cf, by, reference_group, other_groups)
    annot = SingleCellExperiment::colData(scExp_cf)
    if(!is.null(jitter_by)){
        samp[,jitter_by] = 
                annot[,jitter_by][match(rownames(samp),annot$cell_id)]
    }
    p <- ggplot(samp) + geom_violin(aes(
        x=.data[[paste0(by,"_i")]], y=.data[["inter_corr"]],
                    fill=.data[[paste0(by,"_i")]]), alpha=0.8) +
        theme_classic() +
        theme(axis.text.x = element_text(angle=90)) + 
        ylab(paste0("Inter-",by," correlation")) + xlab("")
    
    cols = unique(as.character(annot[,paste0(by,"_color")][
        match(rownames(samp),annot$cell_id)]))
    names(cols) = unique(as.character(annot[,by][match(
        rownames(samp),annot$cell_id)]))
    
    p <- p + scale_fill_manual(values = cols)
    
    if(!is.null(jitter_by)){
        if(grepl("counts_", jitter_by)){
            p <- p + geom_jitter(
                aes(x=.data[[paste0(by,"_i")]], y=.data$inter_corr,
                    color = .data[[jitter_by]]),
                alpha = 0.45) +
                scale_color_gradientn(colours = matlab.like(100))
        } else{
            p <- p + geom_jitter(
                aes(x=.data[[paste0(by,"_i")]], y=.data$inter_corr,
                    color = forcats::fct_inorder(
                        .data[[jitter_by]])),
                alpha = 0.45) +
                scale_color_manual(
                    values = unique(annot[,paste0(jitter_by,"_color")][
                        match(rownames(samp),annot$cell_id)]))
        }
    }
    return(p + theme(legend.title = element_text("")))
}

#' Coverage plot 
#'
#' @param coverages A list containing sample coverage as GenomicRanges
#' @param label_color_list List of colors, list names are labels
#' @param peaks A GRanges object containing peaks location to plot (optional)
#' @param chrom Chromosome
#' @param start Start
#' @param end End
#' @param ref Genomic Reference
#'
#' @return A coverage plot annotated with genes
#' 
#' @importFrom ggrepel geom_text_repel
#' @importFrom gggenes geom_gene_arrow
#' @importFrom dplyr mutate filter
#' @importFrom gridExtra grid.arrange
#' @importFrom S4Vectors subjectHits 
#' @importFrom GenomicRanges GRanges findOverlaps 
#' @export
#'
#' @examples
#' data(scExp)
#' 
plot_coverage_BigWig <- function(
    coverages, label_color_list, peaks = NULL,
    chrom, start, end, ref = "hg38"){
    
    eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
    genebed = eval(parse(text = paste0("", ref, ".GeneTSS")))
    genebed$strand2 = genebed$strand
    genebed$strand = "."
    colnames(genebed)[5:6] = c("score", "strand")
    
    #Select region of interest
    roi = GenomicRanges::GRanges(
        chrom, ranges = IRanges::IRanges(start, end))
    
    coverage_list = list()
    if(class(coverages) != "list") {
        coverage_list[[1]] = coverages[S4Vectors::subjectHits(
            GenomicRanges::findOverlaps(roi,coverages)),]
    } else{
        for(i in seq_along(coverages)){
            coverage_list[[i]] = coverages[[i]][S4Vectors::subjectHits(
                GenomicRanges::findOverlaps(roi,coverages[[i]])),]
        }
    }
    
    if(sum(sapply(coverage_list, length)) >0){
        
        max = round(max(sapply(coverage_list, function(tab) max(tab$score))),3)
        
        n = length(coverage_list)
        layout.matrix <- matrix(
            c(sort(rep(seq_len(n),4)), n+1,n+1), ncol = 1)
        peaks_roi=data.frame()
        if(!is.null(peaks)){
            peaks_roi = peaks[S4Vectors::subjectHits(
                GenomicRanges::findOverlaps(roi,peaks)),]
            peaks_roi = as.data.frame(peaks_roi)
            if(nrow(peaks_roi)>0) layout.matrix <- matrix(
                c(sort(rep(seq_len(n),3)), n+1,n+1,n+2,n+2), ncol = 1)
        }
        
        list_plot = list()    
        for(i in seq_along(coverage_list)){
            coverages_df_tmp = as.data.frame(coverage_list[[i]])
            list_plot[[i]] = coverages_df_tmp %>%
                ggplot(mapping = aes(x = start, y = score)) +
                geom_area(stat = "identity", fill = label_color_list[i]) + 
                geom_hline(yintercept = 0, size = 0.1) +
                ylab(label = names(label_color_list)[i]) + ylim(c(0, max)) + 
                theme_classic() + 
                theme(panel.spacing.y = unit(x = 0, units = "line"),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.line.x = element_blank(),
                      axis.line.y = element_blank()) 
            
        }
        
        if(nrow(peaks_roi)>0) {
            
            peaks_roi$molecule = 
                rep(c(1.2,1.8),
                    length(peaks_roi$start))[1:length(peaks_roi$start)]
            
            list_plot[[length(list_plot) + 1 ]] = peaks_roi %>% 
                ggplot(aes(xmin = start, xmax = end, y = molecule)) +
                gggenes::geom_gene_arrow(colour ="#B82144" , fill = "#B82144",
                                         arrowhead_width = grid::unit(0, "mm"),
                                         arrowhead_height = grid::unit(0, "mm"),
                                         arrow_body_height = grid::unit(2, "mm"),
                                         show.legend = F) +
                theme_classic() + ylim(c(0,2)) + 
                theme(panel.spacing.y = unit(x = 0.5, units = "line"),
                      axis.title.x = element_text(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.line.y = element_blank()) +
                ylab("Peaks")
        }
        
        genebed_tmp = genebed[which(genebed$chr==chrom &
                                        genebed$start > start  & 
                                        genebed$end < end),]
        
        genebed_tmp$orientation = ifelse(genebed_tmp$strand=="+", 1, -1)
                genebed_tmp$molecule = 
            rep(c(1,2.35),length(genebed_tmp$start))[1:length(genebed_tmp$start)]
        if(nrow(genebed_tmp) > 0) { 
        list_plot[[length(list_plot) + 1 ]] = genebed_tmp %>% 
            ggplot(aes(xmin = start, xmax = end, y = molecule, forward = orientation)) +
            gggenes::geom_gene_arrow(fill = "black", 
                                     arrowhead_width = grid::unit(2, "mm"),
                                     arrowhead_height = grid::unit(4, "mm"),
                                     arrow_body_height = grid::unit(1, "mm"),
                                     show.legend = F) +
            ggrepel::geom_text_repel(data = genebed_tmp %>%
                                         dplyr::mutate(start = (.data[["start"]] + .data[["end"]])/2),
                                     aes(x = start, y = molecule, label = Gene),
                                     size = 4,   inherit.aes = F,
                                     nudge_y = -0.1) + 
            theme_classic() + ylim(c(0,2)) + 
            theme(panel.spacing.y = unit(x = 0.5, units = "line"),
                  axis.title.x = element_text(), axis.title = element_blank(),
                  axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  axis.line.y = element_blank()) + xlab(chrom)
        } else{
            list_plot[[length(list_plot) + 1 ]] == NULL
        }
        print(gridExtra::grid.arrange(grobs = list_plot, layout_matrix = layout.matrix))
        
    } else{
        message("ChromSCape::plot_coverage_BigWig - No reads found within the required ranges.")
    }
}


#' Differential summary barplot
#'
#' @param scExp_cf A SingleCellExperiment object
#'
#' @return A barplot summary of differential analysis
#' @export
#'
#' @examples
#' data("scExp")
#' plot_differential_summary_scExp(scExp)
plot_differential_summary_scExp <- function(scExp_cf)
{
    
    stopifnot(is(scExp_cf, "SingleCellExperiment"))
    if (is.null(scExp_cf@metadata$diff)) 
        stop(paste0("ChromSCape::differential_barplot_scExp - ",
                    "No DA, please run differential_analysis_scExp first."))
    
    summary = scExp_cf@metadata$diff$summary
    myylim <- range(c(summary["over", ], -summary["under", ]))
    barplot(summary["over", ], col = "red", las = 1, ylim = myylim,
            main = "Number of differentially enriched regions", 
            ylab = "Number of regions", axes = FALSE)
    barplot(-summary["under", ], col = "forestgreen", ylim = myylim,
            add = TRUE, axes = FALSE, names.arg = "")
    z <- axis(2, pos = -10)
    axis(2, at = z, labels = abs(z), las = 1)
}

#' Differential H1 distribution plot
#'
#' @param scExp_cf A SingleCellExperiment object
#' @param cell_cluster Which cluster to plot
#'
#' @return A barplot of H1 distribution
#' @export
#'
#' @importFrom graphics hist barplot axis plot abline
#' @examples
#' data("scExp")
#' plot_differential_H1_scExp(scExp)
plot_differential_H1_scExp <- function(scExp_cf, cell_cluster = "C1")
{
    
    stopifnot(is(scExp_cf, "SingleCellExperiment"), is.character(cell_cluster))
    if (is.null(scExp_cf@metadata$diff)) 
        stop(paste0("ChromSCape::differential_H1_plot_scExp - No DA, please ",
                    "run differential_analysis_scExp first."))
    
    if (!cell_cluster %in% scExp_cf@metadata$diff$groups) 
        stop(paste0("ChromSCape::differential_H1_plot_scExp - Chromatin group ",
                    "specified doesn't correspond to differential analysis, please rerun ",
                    "run_differential_analysis_scExp first with correct parameters."))
    
    res = scExp_cf@metadata$diff$res
    
    tmp <- H1proportion(res[, paste("pval", cell_cluster, sep = ".")])
    hist(res[, paste("pval", cell_cluster, sep = ".")],
         breaks = seq(0, 1, by = 0.05), xlab = "P-value", ylab = "Frequency",
         main = paste(cell_cluster, "vs the rest", "\n", "H1 proportion:",
                      round(tmp, 3)))
}

#' Volcano plot of differential features
#'
#' @param scExp_cf A SingleCellExperiment object
#' @param cell_cluster Which cluster to plot
#' @param cdiff.th Fold change threshold
#' @param qval.th Adjusted p.value threshold
#'
#' @return A volcano plot of differential analysis of a specific cluster
#' @export
#'
#' @examples
#' data("scExp")
#' plot_differential_volcano_scExp(scExp,"C1")
plot_differential_volcano_scExp <- function(
    scExp_cf, cell_cluster = "C1", cdiff.th = 1, qval.th = 0.01)
{
    
    stopifnot(is(scExp_cf, "SingleCellExperiment"), is.character(cell_cluster), 
              is.numeric(qval.th), is.numeric(cdiff.th))
    
    if (is.null(scExp_cf@metadata$diff)) 
        stop(paste0("ChromSCape::differential_volcano_plot_scExp - No DA, ",
                    "please run differential_analysis_scExp first."))
    
    if (!cell_cluster %in% scExp_cf@metadata$diff$groups) 
        stop(paste0("ChromSCape::differential_volcano_plot_scExp - Chromatin ",
                    "group specified doesn't correspond to differential analysis, please ",
                    "rerun differential_analysis_scExp first with correct parameters."))
    
    res = scExp_cf@metadata$diff$res
    summary = scExp_cf@metadata$diff$summary
    
    mycol <- rep("black", nrow(res))
    mycol[which(res[, paste("qval", cell_cluster, sep = ".")]
                <= qval.th & res[, 
                                 paste("cdiff", cell_cluster, sep = ".")] > cdiff.th)] <- "red"
    mycol[which(res[, paste("qval", cell_cluster, sep = ".")]
                <= qval.th & res[, 
                                 paste("cdiff", cell_cluster, sep = ".")] < -cdiff.th)] <- "forestgreen"
    
    idx = which(scExp_cf@metadata$diff$groups == cell_cluster)
    
    plot(res[, paste("cdiff", cell_cluster, sep = ".")],
         -log10(res[, paste("qval", cell_cluster, sep = ".")]), col = mycol,
         cex = 0.7, pch = 16, xlab = "count difference", 
         ylab = "-log10(adjusted p-value)", las = 1,
         main = paste(cell_cluster, "vs", scExp_cf@metadata$diff$refs[idx] , "\n",
                      summary["over", cell_cluster],
                      "enriched,", summary["under", cell_cluster], "depleted"))
    abline(v = cdiff.th, lty = 2)
    abline(h = -log10(qval.th), lty = 2)
    abline(v = -cdiff.th, lty = 2)
    
}

#' gg_fill_hue
#'
#' @param n num hues
#'
#' @importFrom grDevices hcl
#' @return A color in HEX format 
gg_fill_hue <- function(n)
{
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}


#' Col2Hex
#'
#' Transform character color to hexadecimal color code.
#'
#' @param cname Color name
#'
#' @return The HEX color code of a particular color
#' @importFrom grDevices col2rgb rgb
#'
col2hex <- function(cname)
{
    colMat <- grDevices::col2rgb(cname)
    grDevices::rgb(
        red = colMat[1, ]/255, green = colMat[2, ]/255, blue = colMat[3,]/255)
}

#' Plot cluster consensus
#'
#' Plot cluster consensus score for each k as a bargraph.
#'
#' @param scExp A SingleCellExperiment
#'
#' @return The consensus score for each cluster for each k as a barplot
#' @importFrom dplyr as_tibble
#' @importFrom ggplot2 ggplot geom_bar facet_grid scale_fill_manual 
#' element_blank theme_minimal theme position_dodge ylab
#' 
#' @export
#'
#' @examples
#' data("scExp")
#' plot_cluster_consensus_scExp(scExp)
plot_cluster_consensus_scExp <- function(scExp)
{
    stopifnot(is(scExp,"SingleCellExperiment"))
    if(!"icl" %in% names(scExp@metadata))
        stop(paste0("plot_cluster_consensus_scExp - please run ",
                    "consensus_clustering_scExp first."))
    
    if(!"consclust" %in% names(scExp@metadata))
        stop(paste0("plot_cluster_consensus_scExp - please run ",
                    "consensus_clustering_scExp first."))
    cc = dplyr::as_tibble(scExp@metadata$icl$clusterConsensus)
    colors = unique(as.character(scExp@metadata$consclust[[1]]))
    colors = c(colors,"#B55274")
    cc$k = factor(paste0("k=",cc$k), levels=unique(paste0("k=",cc$k)))
    cc$cluster =  factor(
        paste0("C",cc$cluster), levels = unique(paste0("C",cc$cluster)))
    p = cc %>% ggplot(aes(x=.data$cluster, y=.data$clusterConsensus,
                          fill=.data$cluster)) +
        geom_bar(stat = "identity", position=position_dodge(width=0.9)) +
        facet_grid(.~as.factor(k), scales="free_x", space="free") +
        theme_minimal() + theme(panel.grid = element_blank(),
                                axis.text.x = element_blank()) +
        scale_fill_manual(values=unique(colors)) + ylab("Consensus Score")
    return(p)
}
