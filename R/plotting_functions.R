## Authors : Pacôme Prompsy Title : Wrappers & function to create variety of plot
## to uncover heterogeneity in single cell Dataset


#' Plotting distribution of signal
#'
#' @param scExp 
#' @param raw 
#' @param log10 
#' @param pseudo_counts 
#' @param bins 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @import ggplot2
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment assayNames
plot_distribution_scExp <- function(scExp, raw = T, log10 = F, pseudo_counts = 1, 
    bins = 150)
    {
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(pseudo_counts))
    if (!raw %in% c(T, F) | !log10 %in% c(T, F)) 
        stop("ChromSCape::plot_distribution_scExp - raw and log10 must be true or false.")
    if (raw == F && !("normcounts" %in% SummarizedExperiment::assayNames(scExp))) 
        stop("ChromSCape::plot_distribution_scExp - If raw is false, normcounts must not be empty - run normalize_scExp first.")
    
    if (raw) 
        cell_cov_df = data.frame(coverage = Matrix::colSums(counts(scExp))) else cell_cov_df = data.frame(coverage = Matrix::colSums(normcounts(scExp)))
    
    if (log10) 
        cell_cov_df$coverage = log10(cell_cov_df$coverage + pseudo_counts)
    
    ggplot(cell_cov_df, aes(x = coverage)) + geom_histogram(color = "black", fill = "steelblue", 
        bins = bins) + labs(x = "read coverage per cell") + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA))
}

#' Adding colors to cells & features
#'
#' @param scExp 
#' @param annotCol 
#' @param color_by 
#' @param color_df 
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment colData
colors_scExp <- function(scExp, annotCol = "sample_id", color_by = "sample_id", color_df = NULL)
{
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(annotCol), is.character(color_by))
    
    # initialization
    annot = as.data.frame(SingleCellExperiment::colData(scExp))
    anocol <- geco.annotToCol2(annotS = annot[, annotCol, drop = F], annotT = annot, 
        plotLegend = F, categCol = NULL)
    SummarizedExperiment::colData(scExp)[, paste0(annotCol, "_color")] = as.data.frame(anocol, 
        stringsAsFactors = F)[, annotCol]  # factor or not ?
    
    if (!is.null(color_df))
    {
        # add custom colors
        if (!color_by %in% colnames(color_df)) 
            stop("ChromSCape::color_scExp - color_by must be present in colnames of color_df is not null.")
        
        SummarizedExperiment::colData(scExp)[, paste0(annotCol, "_color")] = color_df[match(SingleCellExperiment::colData(scExp)[, 
            color_by], color_df[, color_by]), paste0(color_by, "_color"), drop = F]
    }
    return(scExp)
}

#' Get color dataframe from shiny::colorInput 
#'
#' @param input 
#' @param levels_selected 
#' @param color_by 
#' @param input_id_prefix 
#'
#' @return
#'
#' @export
#' @importFrom tibble rownames_to_column
get_color_dataframe_from_input <- function(input, levels_selected, color_by = c("sample_id", 
    "total_counts"), input_id_prefix = "color_")
    {
    stopifnot(!is.null(input), is.character(levels_selected), is.character(color_by))
    
    # Get colors
    color_list <- paste0("list(", paste0(levels_selected, " = input$", input_id_prefix, 
        levels_selected, collapse = ", "), ")")

    color_list <- eval(parse(text = color_list))
    # Transform into dataframe with right column names
    color_df = as.matrix(color_list) %>% as.data.frame(stringsAsFactors = F) %>% 
        tibble::rownames_to_column(color_by)
    color_df[, 2] = as.character(color_df[, 2])
    colnames(color_df)[2] = paste0(color_by, "_color")
    
    return(color_df)
}

# Wrapper for plotting PCA & TSNE & UMAP
#' Plot reduced dimensions (PCA, TSNE, UMAP)
#'
#' @param scExp 
#' @param color_by 
#' @param reduced_dim 
#' @param select_x 
#' @param select_y 
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom SingleCellExperiment reducedDim reducedDimNames colData
#' @import ggplot2
#' @importFrom colorRamps matlab.like

plot_reduced_dim_scExp <- function(scExp, color_by = "sample_id", reduced_dim = c("PCA", 
    "TSNE", "UMAP"), select_x = "Component_1", select_y = "Component_2")
    {
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(color_by), is.character(reduced_dim), 
        is.character(select_x), is.character(select_y))
    
    if (!reduced_dim[1] %in% SingleCellExperiment::reducedDimNames(scExp)) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - ", reduced_dim[1], " is not present in object, please run normalize_scExp first."))
    
    if (!color_by %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::plot_reduced_dim_scExp - color_by must be present in colnames of colData(scExp).")
    
    if (!paste0(color_by, "_color") %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::plot_reduced_dim_scExp - color_by's color column must be present 
         in colnames of colData(scExp). Please run color_scExp first.")
    
    if (!select_x %in% colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop("ChromSCape::plot_reduced_dim_scExp - select_x must be present in colnames of PCA of scExp.")
    
    if (!select_y %in% colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop("ChromSCape::plot_reduced_dim_scExp - select_y must be present in colnames of PCA of scExp.")
    
    plot_df = as.data.frame(cbind(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]), 
        SingleCellExperiment::colData(scExp)))
    
    p <- ggplot(plot_df, aes_string(x = select_x, y = select_y)) + geom_point(alpha = 0.6, 
        aes(color = SingleCellExperiment::colData(scExp)[, color_by])) + labs(color = color_by) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            panel.border = element_rect(colour = "black", fill = NA))
    
    if (color_by == "total_counts")
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

#' Plot cell correlation heatmap with annotations
#'
#' @param scExp 
#' @param name_hc 
#' @param corColors 
#' @param color_by 
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom SingleCellExperiment reducedDim reducedDimNames colData
#' @importFrom grDevices colorRampPalette
plot_heatmap_scExp <- function(scExp, name_hc = "hc_cor", corColors = (grDevices::colorRampPalette(c("royalblue", 
    "white", "indianred1")))(256),
    color_by = NULL)
    {
    
    # make a variety of sanity check
    stopifnot(is(scExp, "SingleCellExperiment"))
    if (!"Cor" %in% SingleCellExperiment::reducedDimNames(scExp)) 
        stop("ChromSCape::plot_heatmap_scExp - No correlation, run correlation_and_hierarchical_clust_scExp before filtering.")
    
    if (!name_hc %in% names(scExp@metadata)) 
        stop("ChromSCape::plot_heatmap_scExp - No dendrogram, run correlation_and_hierarchical_clust_scExp before filtering.")
    
    if (length(scExp@metadata[[name_hc]]$order) != ncol(scExp)) 
        stop("ChromSCape::plot_heatmap_scExp - Dendrogram has different number of cells than dataset.")

    anocol = as.matrix(SingleCellExperiment::colData(scExp)[, grep("_color", colnames(SingleCellExperiment::colData(scExp))), 
        drop = F])
    colnames(anocol) = gsub("_color","",colnames(anocol))
    
    if(!is.null(color_by)) {
        if(length(intersect(color_by,colnames(anocol)))>0) 
            anocol = anocol[,color_by,drop=F]
    }
    
    return(
        geco.hclustAnnotHeatmapPlot(x = SingleCellExperiment::reducedDim(scExp, "Cor")[scExp@metadata[[name_hc]]$order, 
        scExp@metadata[[name_hc]]$order], hc = scExp@metadata[[name_hc]], hmColors = corColors, 
        anocol = as.matrix(anocol[scExp@metadata[[name_hc]]$order, ]), xpos = c(0.15, 
            0.9, 0.164, 0.885), ypos = c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95), dendro.cex = 0.04, 
        xlab.cex = 0.8, hmRowNames = FALSE)
    )
}


#' Differential summary barplot
#'
#' @param scExp_cf 
#'
#' @return
#' @export
#'
#' @examples
plot_differential_summary_scExp <- function(scExp_cf)
{
    
    # make a variety of sanity check
    stopifnot(is(scExp_cf, "SingleCellExperiment"))
    if (is.null(scExp_cf@metadata$diff)) 
        stop("ChromSCape::differential_barplot_scExp - No DA, please run differential_analysis_scExp first.")
    
    summary = scExp_cf@metadata$diff$summary
    myylim <- range(c(summary["over", ], -summary["under", ]))
    barplot(summary["over", ], col = "red", las = 1, ylim = myylim, main = "Differentially bound regions", 
        ylab = "Number of regions", axes = F)
    barplot(-summary["under", ], col = "forestgreen", ylim = myylim, add = T, axes = F, 
        names.arg = "")
    z <- axis(2, pos = -10)
    axis(2, at = z, labels = abs(z), las = 1)
}

#' Differential H1 distribution plot
#'
#' @param scExp_cf 
#' @param cell_cluster 
#'
#' @return
#' @export
#'
#' @examples
plot_differential_H1_scExp <- function(scExp_cf, cell_cluster = "C1")
{
    
    # make a variety of sanity check
    stopifnot(is(scExp_cf, "SingleCellExperiment"), is.character(cell_cluster))
    if (is.null(scExp_cf@metadata$diff)) 
        stop("ChromSCape::differential_H1_plot_scExp - No DA, please run differential_analysis_scExp first.")
    
    if (!cell_cluster %in% scExp_cf@metadata$diff$groups) 
        stop("ChromSCape::differential_H1_plot_scExp - Chromatin group specified doesn't correspond to differential 
         analysis, please rerun run_differential_analysis_scExp first with correct parameters.")
    
    res = scExp_cf@metadata$diff$res
    
    tmp <- geco.H1proportion(res[, paste("pval", cell_cluster, sep = ".")])
    hist(res[, paste("pval", cell_cluster, sep = ".")], breaks = seq(0, 1, by = 0.05), 
        xlab = "P-value", ylab = "Frequency", main = paste(cell_cluster, "vs the rest", 
            "\n", "H1 proportion:", round(tmp, 3)))
}

#' Volcano plot of differential features
#'
#' @param scExp_cf 
#' @param cell_cluster 
#' @param cdiff.th 
#' @param qval.th 
#'
#' @return
#' @export
#'
#' @examples
plot_differential_volcano_scExp <- function(scExp_cf, cell_cluster = "C1", cdiff.th = 1, 
    qval.th = 0.01)
    {
    
    # make a variety of sanity check
    stopifnot(is(scExp_cf, "SingleCellExperiment"), is.character(cell_cluster), 
        is.numeric(qval.th), is.numeric(cdiff.th))
    
    if (is.null(scExp_cf@metadata$diff)) 
        stop("ChromSCape::differential_volcano_plot_scExp - No DA, please run differential_analysis_scExp first.")
    
    if (!cell_cluster %in% scExp_cf@metadata$diff$groups) 
        stop("ChromSCape::differential_volcano_plot_scExp - Chromatin group specified doesn't correspond to differential 
         analysis, please rerun differential_analysis_scExp first with correct parameters.")
    
    res = scExp_cf@metadata$diff$res
    summary = scExp_cf@metadata$diff$summary
    
    mycol <- rep("black", nrow(res))
    mycol[which(res[, paste("qval", cell_cluster, sep = ".")] <= qval.th & res[, 
        paste("cdiff", cell_cluster, sep = ".")] > cdiff.th)] <- "red"
    mycol[which(res[, paste("qval", cell_cluster, sep = ".")] <= qval.th & res[, 
        paste("cdiff", cell_cluster, sep = ".")] < -cdiff.th)] <- "forestgreen"
    plot(res[, paste("cdiff", cell_cluster, sep = ".")], -log10(res[, paste("qval", 
        cell_cluster, sep = ".")]), col = mycol, cex = 0.7, pch = 16, xlab = "count difference", 
        ylab = "-log10(adjusted p-value)", las = 1, main = paste(cell_cluster, 
            "vs the rest", "\n", summary["over", cell_cluster], "enriched,", summary["under", 
                cell_cluster], "depleted"))
    abline(v = cdiff.th, lty = 2)
    abline(h = -log10(qval.th), lty = 2)
    abline(v = -cdiff.th, lty = 2)
    
}

#' gg_fill_hue
#'
#' @param n num hues
#'
#' @importFrom grDevices hcl
#' @export
gg_fill_hue <- function(n)
{
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Col2Hex
#'
#' Transform character color to hexadecimal color code.
#'
#' @param cname 
#'
#' @return
#' @importFrom grDevices col2rgb rgb
#'
#' @examples
col2hex <- function(cname)
{
    colMat <- grDevices::col2rgb(cname)
    grDevices::rgb(red = colMat[1, ]/255, green = colMat[2, ]/255, blue = colMat[3, 
        ]/255)
}

#' Plot cluster consensus
#'
#' Plot cluster consensus score for each k as a bargraph.
#'
#' @param cname 
#'
#' @return
#' @importFrom dplyr as_tibble
#' @import ggplot2
#' 
#' @export
#'
#' @examples
plot_cluster_consensus_scExp <- function(scExp)
{
    stopifnot(is(scExp,"SingleCellExperiment"))
    if(!"icl" %in% names(scExp@metadata))
        stop("plot_cluster_consensus_scExp - please run consensus_clustering_scExp 
             first.")
    
    if(!"consclust" %in% names(scExp@metadata))
        stop("plot_cluster_consensus_scExp - please run consensus_clustering_scExp 
             first.")
    cc = dplyr::as_tibble(scExp@metadata$icl$clusterConsensus)
    colors = unique(as.character(scExp@metadata$consclust[[1]]))
    colors = c(colors,"#B55274")
    cc$k = factor(paste0("k=",cc$k), levels=unique(paste0("k=",cc$k)))
    cc$cluster =  factor(paste0("C",cc$cluster), levels = unique(paste0("C",cc$cluster)))
    p = cc %>% ggplot(aes(x=cluster, y=clusterConsensus, fill=cluster)) +
        geom_bar(stat = "identity", position=position_dodge(width=0.9)) +
        facet_grid(.~as.factor(k), scales="free_x", space="free") +
        theme_minimal() + theme(panel.grid = element_blank(),
                                axis.text.x = element_blank()) +
        scale_fill_manual(values=unique(colors)) + ylab("Consensus Score")

    return(p)
}
