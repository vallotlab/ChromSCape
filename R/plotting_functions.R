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
#' @import ggplot2
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
#' @param color_by Feature used for coloration
#' @param reduced_dim Reduced Dimension used for plotting
#' @param select_x Which variable to select for x axis
#' @param select_y Which variable to select for y axis
#' @param downsample Number of cells to downsample
#' @param transparency Alpha parameter, between 0 and 1
#' 
#' @return A ggplot geom_point plot of reduced dimension 2D reprensentation 
#' @export
#'
#' @importFrom SingleCellExperiment reducedDim reducedDimNames colData
#' @import ggplot2
#' @importFrom colorRamps matlab.like
#' 
#' @examples
#' data("scExp")
#' plot_reduced_dim_scExp(scExp, color_by = "sample_id")
#' plot_reduced_dim_scExp(scExp, color_by = "total_counts")
#' plot_reduced_dim_scExp(scExp, reduced_dim = "UMAP")
#' 
plot_reduced_dim_scExp <- function(
    scExp, color_by = "sample_id", reduced_dim = c("PCA", "TSNE", "UMAP"),
    select_x = "Component_1", select_y = "Component_2", downsample = 5000, 
    transparency = 0.6,  size = 1)
{
    warning_plot_reduced_dim_scExp(scExp, color_by , reduced_dim,
                                   select_x, select_y, downsample, transparency)
    if(ncol(scExp) > downsample) scExp = scExp[,sample(ncol(scExp),
                                                       downsample,
                                                       replace = FALSE)]
    plot_df = as.data.frame(
        cbind(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]), 
              SingleCellExperiment::colData(scExp)))
    
    p <- ggplot(plot_df, aes_string(x = select_x, y = select_y)) + 
        geom_point(alpha = transparency, size = size, aes(
            color = SingleCellExperiment::colData(scExp)[, color_by])) +
        labs(color = color_by) + 
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA))
    
    if (is.numeric(SingleCellExperiment::colData(scExp)[,color_by]))
    {
        p <- p + scale_color_gradientn(colours = matlab.like(100))
    } else
    {
        
        cols = unique(as.character(
            SingleCellExperiment::colData(scExp)[,paste0(color_by, "_color")]))
        names(cols) = unique(as.character(
            SingleCellExperiment::colData(scExp)[,color_by]))
        p <- p + scale_color_manual(values = cols)
    }
    return(p)
}

#' A warning helper for plot_reduced_dim_scExp
#'
#' @param scExp A SingleCellExperiment Object
#' @param color_by Feature used for coloration
#' @param reduced_dim Reduced Dimension used for plotting
#' @param select_x Which variable to select for x axis
#' @param select_y Which variable to select for y axis
#' @param downsample Number of cells to downsample
#' @param transparency Alpha parameter, between 0 and 1
#'
#' @return Warning or errors if the inputs are not correct
#' 
warning_plot_reduced_dim_scExp <- function(scExp, color_by , reduced_dim,
                                           select_x, select_y, downsample,
                                           transparency){
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(color_by),
              is.character(reduced_dim), is.character(select_x),
              is.character(select_y), is.numeric(downsample),
              is.numeric(transparency))
    if (!reduced_dim[1] %in% SingleCellExperiment::reducedDimNames(scExp)) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - ", reduced_dim[1],
                    " is not present in object, please run normalize_scExp ",
                    "first."))
    
    if (!color_by %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - color_by must be ",
                    "present in colnames of colData(scExp)."))
    
    if (!paste0(color_by, "_color") %in%
        colnames(SingleCellExperiment::colData(scExp))) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - color_by's color ",
                    "column must be present in colnames of colData(scExp). ",
                    "Please run colors_scExp first."))
    
    if (!select_x %in% 
        colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - select_x must be ",
                    "present in colnames of PCA of scExp."))
    
    if (!select_y %in% 
        colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - select_y must be",
                    " present in colnames of PCA of scExp."))
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
#'  plot_most_contributing_features(scExp, component = "Component_1")
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
        cor_mat = SingleCellExperiment::reducedDim(scExp, "Cor")
        hc_cor = stats::hclust(stats::as.dist(1 - cor_mat), method = hc_linkage)
        hc_cor$labels = rep("",length(hc_cor$labels))
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
            x = SingleCellExperiment::reducedDim(scExp, "Cor")[
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
#' @examples plot_intra_correlation_scExp(scExp)
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
        if(jitter_by == "total_counts"){
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
#' @param downsample Downsample for plotting
#'
#' @return A violin plot of inter-correlation
#' @export
#'
#' @importFrom forcats fct_inorder
#' 
#' @examples plot_intra_correlation_scExp(scExp)
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
        if(jitter_by == "total_counts"){
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

#' Coverage plot using Sushi
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
#' @importFrom Sushi plotBedgraph plotGenes
#' @importFrom S4Vectors subjectHits 
#' @importFrom GenomicRanges GRanges findOverlaps 
#' @export
#'
#' @examples plot_intra_correlation_scExp(scExp)
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

    gl.sub <- genebed[which(genebed[,"chr"] == chrom),]

    for(i in seq_along(coverages)){
        coverages[[i]] = coverages[[i]][S4Vectors::subjectHits(
            GenomicRanges::findOverlaps(roi,coverages[[i]])),]
    }

    if(sum(sapply(coverages, length)) >0){

        max = round(max(sapply(coverages, function(tab) max(tab$score))),3)
    
        n = length(coverages)
        layout.matrix <- matrix(
            c(sort(rep(seq_len(n),3)), n+1,n+1,n+2), ncol = 1)
        peaks_roi=data.frame()
        if(!is.null(peaks)){
            peaks_roi = peaks[S4Vectors::subjectHits(
                GenomicRanges::findOverlaps(roi,peaks)),]
            peaks_roi = as.data.frame(peaks_roi)
            if(nrow(peaks_roi)>0) layout.matrix <- matrix(
                c(sort(rep(seq_len(n),3)), n+1,n+1,n+2,n+2,n+3), ncol = 1)
        }

        graphics::layout(mat = layout.matrix,
                         heights = c(1), # Heights of the two rows
                         widths = c(1)) # Widths of the two columns
        
        par(cex=0.5)
        par(mar = c(0.75, 6, 0, 1), oma = c(1, 1, 1, 1))
        for(i in seq_along(coverages)){

                Sushi::plotBedgraph(as.data.frame(coverages[[i]])[,c(1,2,3,6)],
                                chrom, start, end,
                                range = c(0, max),
                                addscale = TRUE, 
                                ylab= names(label_color_list)[i],
                                color= label_color_list[[i]],
                                cex.lab=2,cex.main=2.1)

        }
        par(mar = c(1, 6, 1, 1),xpd=NA)
        if(nrow(peaks_roi)>0) Sushi::plotBed(peaks_roi, chrom, start,
                                             end,row = "supplied",
                                             rowlabels = "Peaks",
                                             rowlabelcex = 1.5)
        Sushi::plotGenes(gl.sub, chrom, start, end,
                  bentline=FALSE, plotgenetype = "arrow",
                  labeltext = TRUE, labelat = "start", fontsize=1,
                  labeloffset = 0.4, bheight=0.07)
        
        par(mar = c(1, 6, 1, 1))
        Sushi::labelgenome(chrom, start, end, n=4, scale="Mb", cex.axis=1.5)
        
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
    
    plot(res[, paste("cdiff", cell_cluster, sep = ".")],
         -log10(res[, paste("qval", cell_cluster, sep = ".")]), col = mycol,
         cex = 0.7, pch = 16, xlab = "count difference", 
         ylab = "-log10(adjusted p-value)", las = 1,
         main = paste(cell_cluster, "vs the rest", "\n",
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
#' @import ggplot2
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
