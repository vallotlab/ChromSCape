
#' @title Map features onto cytobands  
#'   
#' @description Map the features of a SingleCellExperiment onto the cytobands 
#' of a given genome. Some features might not be mapped to any cytobands (e.g.
#'  if they are not in the canconical chromosomes), and are removed from the
#'  returned object.
#'
#' @details 
#' The cytobands are an arbitrary cutting of the genome that dates back to 
#' staining metaphase chromosomes with Giemsa.
#' 
#' @param scExp A SingleCellExperiment with genomic coordinate as features 
#' (peaks or bins) 
#' @param ref_genome Reference genome ('hg38' or 'mm10') 
#'
#' @return A data.frame of the SCE features with their corresponding cytoband
#' name
#' 
#' @importFrom SummarizedExperiment rowRanges seqnames
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicRanges findOverlaps
#' @export
#' @examples 
#' data("scExp")
#' matching_cyto = get_cyto_features(scExp, ref_genome="hg38")
#'  
get_cyto_features <- function(scExp, ref_genome = c("hg38", "mm10")[1] ){
    
    eval(parse(text = paste0("data(",ref_genome, ".chromosomes)")))
    canonical_chr <- eval(parse(text = paste0(ref_genome, ".chromosomes")))
    canonical_chr$start = 1
    canonical_chr <- as(canonical_chr, "GRanges")
    
    eval(parse(text = paste0("data(",ref_genome, ".cytoBand)")))
    cyto <- eval(parse(text = paste0(ref_genome, ".cytoBand")))
        
    cyto = as(cyto,"GRanges")
    
    # Get genomic bins
    bins = SummarizedExperiment::rowRanges(scExp)
    bins = bins[which(as.character(SummarizedExperiment::seqnames(bins)) %in%
                          as.character(SummarizedExperiment::seqnames(canonical_chr)))]
    
    matching_cyto = data.frame(ID = bins$ID, hits = "")
    
    # Get which bins goes into which cyto band
    for(i in unique(cyto$cytoBand)){
        cyto. = cyto[which(cyto$cytoBand == i)]
        matching_cyto$hits[S4Vectors::queryHits(GenomicRanges::findOverlaps(bins,cyto.))] = i
    }
    
    matching_cyto = matching_cyto[which(matching_cyto$hits!=""),]
    return(matching_cyto)
}

#' @title  Calculate Fraction of reads in each cytobands
#' 
#' @description
#' Re-Count binned reads onto cytobands and calculate the fraction of reads in 
#' each of the cytoband in each cell. For each cell, the fraction of reads in
#'  any given cytoband is calculated. 
#' Cytobands are considered large enough in order that a variation at the 
#' cytoband level is not considered as an epigenetic event but as a genetic 
#' event, e.g. Copy Number Alterations.
#' 
#' @param scExp A SingleCellExperiment with genomic coordinate as features 
#' (peaks or bins) 
#' @param ref_genome Reference genome ('hg38' or 'mm10') 
#'
#' @return The SCE with the fraction of reads in each cytobands in each cells
#' (of dimension cell x cytoband ) in the  reducedDim slot "cytoBand".
#' 
#' @importFrom dplyr group_by summarise everything across
#' @importFrom Matrix colSums
#' @importFrom SingleCellExperiment reducedDim
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_cyto_mat(scExp, ref_genome="hg38")
#' SingleCellExperiment::reducedDim(scExp, "cytoBand") 
#' 
calculate_cyto_mat <- function(scExp, ref_genome = c("hg38","mm10")[1]){
    matching_cyto = get_cyto_features(scExp,ref_genome)
    mat = counts(scExp)[match(matching_cyto$ID, rownames(scExp)),]
    gc()
    
    # Group counts by cytobands
    cyto_mat = as.data.frame(as.matrix(mat)) %>% dplyr::group_by(matching_cyto$hits) %>%
        dplyr::summarise(dplyr::across(.cols = dplyr::everything(),sum)) %>% 
        as.data.frame()
    gc()

    # Normalize by library size
    cyto_mat. = cyto_mat[,-1]
    rownames(cyto_mat) = cyto_mat$`matching_cyto$hits`
    cyto_mat[,-1] = 100* base::sweep(cyto_mat.,2,Matrix::colSums(cyto_mat.),`/`)
    cyto_mat = cyto_mat[,-1]
    
    SingleCellExperiment::reducedDim(scExp, "cytoBand") = t(cyto_mat)
    return(scExp)
}


#' @title Retrieve the cytobands with the most variable fraction of reads
#' 
#' @description
#' Given a SingleCellExperiment object with the slot "cytoBand" containing the
#' fraction of reads in each cytoband, calculates the variance of each cytoband
#' and returns a data.frame with the top variables cytobands. Most cytobands are
#' expected to be unchanged between normal and tumor samples, therefore focusing
#' on the top variable cytobands enable to focus on the most interseting 
#' regions. 
#' 
#' @param scExp A SingleCellExperiment with "cytoBand" reducedDim slot filled.
#' @param top Number of cytobands to return (50).
#'
#' @return A data.frame of the top variable cytoBands and their variance
#' 
#' @importFrom dplyr group_by summarise desc arrange
#' @importFrom tidyr gather
#' @importFrom SingleCellExperiment reducedDim
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_cyto_mat(scExp, ref_genome="hg38")
#' get_most_variable_cyto(scExp, top=50)
#'  
get_most_variable_cyto <- function(scExp, top = 50){
    # Retrieve top variable cytobands
    cyto_df = as.data.frame(
        t(SingleCellExperiment::reducedDim(scExp, "cytoBand")))
    cyto_df$cytoBand = rownames(cyto_df)
    
    cyto_df = cyto_df %>% tidyr::gather("cell_id","Fri_cyto",-ncol(cyto_df))
    top_var = cyto_df %>% dplyr::group_by(cytoBand) %>%
        dplyr::summarise("var" = stats::var(Fri_cyto)) %>% 
        dplyr::arrange(dplyr::desc(.data[["var"]])) %>% utils::head(n = top)
    
    return(top_var)
}


#' @title Calculate the log2-ratio of tumor vs normal fraction of reads in 
#' cytobands 
#' 
#' @description
#' Given a SingleCellExperiment object with the slot "cytoBand" containing the
#' fraction of reads in each cytoband, calculates the log2-ratio of tumor vs
#' normal fraction of reads in cytobands, cell by cell.  
#' If the average signal in normal sample in a cytoband is 0, set this value to
#' 1 so that the ratio won't affect the fraction of read value.
#' 
#' @param scExp A SingleCellExperiment with "cytoBand" reducedDim slot filled.
#' - see  \code{\link{calculate_cyto_mat}}
#' @param controls Sample IDs or Cell IDs of the normal sample to take as 
#' reference.
#'
#' @return The SCE with the log2-ratio of fraction of reads in each cytobands 
#' in each cells (of dimension cell x cytoband ) in the  reducedDim 
#' slot "logRatio_cytoBand".
#' 
#' @importFrom Matrix colMeans
#' @importFrom SingleCellExperiment reducedDim
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_cyto_mat(scExp, ref_genome="hg38")
#' scExp = calculate_logRatio_CNA(scExp, controls=unique(scExp$sample_id)[1])
#' SingleCellExperiment::reducedDim(scExp, "logRatio_cytoBand")
#' 
calculate_logRatio_CNA <- function(scExp, controls){
    cyto_mat = SingleCellExperiment::reducedDim(scExp, "cytoBand")
    
    is_cell_id = ifelse(all(controls %in% scExp$cell_id), TRUE, FALSE)
    if(is_cell_id){
        normal_signal_per_chr = setNames(Matrix::colMeans(cyto_mat[controls,]),
                                         colnames(cyto_mat[controls,]))
    } else{
        controls_cell = as.character(
            scExp$cell_id[which(scExp$sample_id %in% controls)])
        normal_signal_per_chr = setNames(Matrix::colMeans(cyto_mat[controls_cell,]),
                                         colnames(cyto_mat[controls_cell,]))
    }
    
    if(any(normal_signal_per_chr==0)){
        warning("One of the cytoBand has no reads in control. ",
             "Empty control cytobands logRatios are set to 1.")
        normal_signal_per_chr[which(normal_signal_per_chr==0)] = 1
    }
    
    # Calculate log2 ratio of signal / normal
    logRatio_mat = sapply(colnames(cyto_mat), function(i){
            log2(cyto_mat[,i] / normal_signal_per_chr[i])
    }) 
    if(any(is.infinite(logRatio_mat))) logRatio_mat[which(
        is.infinite(logRatio_mat))] = 0
    
    SingleCellExperiment::reducedDim(scExp, "logRatio_cytoBand") = logRatio_mat
    return(scExp)
} 

#' @title Estimate the copy gains/loss of tumor vs normal based on log2-ratio of
#' fraction of reads
#' 
#' @description
#' Given a SingleCellExperiment object with the slot "logRatio_cytoBand" containing the
#' log2-ratio of the fraction of reads in each cytoband, estimate if the 
#' cytoband was lost or acquired a gain in a non-quantitative way. To do so, 
#' the quantiles distribution of the normal cells are calculated, and any 
#' cytoband below or above will be considered as a loss/gain. The False 
#' Discovery Rate is directly proportional to the quantiles.  
#' 
#' @param scExp A SingleCellExperiment with "logRatio_cytoBand" reducedDim slot 
#' filled. See  \code{\link{calculate_logRatio_CNA}}
#' @param controls Sample IDs or Cell IDs of the normal sample to take as 
#' reference.
#' @param quantiles Quantiles of normal log2-ratio distribution below/above 
#' which cytoband is considered to be a loss/gain. (c(0.05,0.95))
#'
#' @return The SCE with the gain or loss in each cytobands 
#' in each cells (of dimension cell x cytoband ) in the  reducedDim 
#' slot "gainOrLoss_cytoBand".
#' 
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_cyto_mat(scExp, ref_genome="hg38")
#' scExp = calculate_logRatio_CNA(scExp, controls=unique(scExp$sample_id)[1])
#' scExp = calculate_gain_or_loss(scExp, controls=unique(scExp$sample_id)[1])
#' SingleCellExperiment::reducedDim(scExp, "gainOrLoss_cytoBand")
#' 
calculate_gain_or_loss <- function(scExp, controls, quantiles = c(0.05, 0.95)){
    logRatio_mat = SingleCellExperiment::reducedDim(scExp, "logRatio_cytoBand")
    
    is_cell_id = ifelse(all(controls %in% scExp$cell_id), TRUE, FALSE)
    if(!is_cell_id){
        controls = as.character(
            scExp$cell_id[which(scExp$sample_id %in% controls)]
        )
    }
    
    gain_or_loss_cytoBand = logRatio_mat
    gain_or_loss_cytoBand[,] = 0
    for(cyto in colnames(logRatio_mat) ){
        tab_tmp = logRatio_mat[,cyto]
        loss_th = quantile(tab_tmp[controls],quantiles[1])
        gain_th = quantile(tab_tmp[controls],quantiles[2])
        
        gain_or_loss_cytoBand[which(logRatio_mat[, cyto] < loss_th), cyto] = -1 
        gain_or_loss_cytoBand[which(logRatio_mat[, cyto] > gain_th), cyto] = 1
    }
    
    SingleCellExperiment::reducedDim(scExp, "gainOrLoss_cytoBand") =
        gain_or_loss_cytoBand
    return(scExp)
}


#' @title Estimate copy number alterations in cytobands
#' 
#' @description
#' Cytobands are considered large enough in order that a variation at the 
#' cytoband level is not considered as an epigenetic event but as a genetic 
#' event, e.g. Copy Number Alterations. The function successively : 
#' - Calculates the fraction of reads in each cytoband (FrCyto). See \code{\link{calculate_cyto_mat}}
#' - Calculates the log2-ratio FrCyto of each cell by the average FrCyto in normal cells. See  \code{\link{calculate_logRatio_CNA}}
#' - Estimates if there was a gain or a loss of copy in each cyto band. See  \code{\link{calculate_gain_or_loss}}
#'   
#'  The corresponding matrices are accessibles in the reducedDim slots
#'   "cytoBands", "logRatio_cytoBands" and "gainOrLoss_cytoBands" respectively.
#'
#' @param scExp A SingleCellExperiment with "logRatio_cytoBand" reducedDim slot 
#' filled. See  \code{\link{calculate_logRatio_CNA}}
#' @param control_samples Sample IDs of the normal sample to take as 
#' reference.
#' @param ref_genome Reference genome ('hg38' or 'mm10') 
#' @param quantiles_to_define_gol Quantiles of normal log2-ratio distribution
#'  below/above which cytoband is considered to be a loss/gain. (c(0.05,0.95)). 
#' See \code{\link{calculate_gain_or_loss}}
#'
#' @return The SCE with the fraction of reads, log2-ratio and gain or loss
#'  in each cytobands in each cells (of dimension cell x cytoband) 
#'  in the  reducedDim slots.
#' 
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_CNA(scExp,  control_samples = unique(scExp$sample_id)[1],
#' ref_genome="hg38", quantiles_to_define_gol = c(0.05,0.95))
#' SingleCellExperiment::reducedDim(scExp, "cytoBand")
#' SingleCellExperiment::reducedDim(scExp, "logRatio_cytoBand")
#' SingleCellExperiment::reducedDim(scExp, "gainOrLoss_cytoBand")
#'
calculate_CNA <- function(scExp, control_samples = unique(scExp$sample_id)[1],
                          ref_genome = c("hg38","mm10")[1],
                          quantiles_to_define_gol = c(0.05,0.95)
                          ){
    scExp = calculate_cyto_mat(scExp, ref_genome = ref_genome)
    scExp = calculate_logRatio_CNA(scExp, controls = control_samples)
    scExp = calculate_gain_or_loss(scExp, controls = control_samples,
                                   quantiles =  quantiles_to_define_gol)
    return(scExp)
}

#' @title Plot Gain or Loss of cytobands of the most variables cytobands
#' 
#' @param scExp A SingleCellExperiment with "gainOrLoss_cytoBand" reducedDim slot 
#' filled. See  \code{\link{calculate_gain_or_loss}}
#' @param cells Cell IDs of the tumor samples to 
#' @param top Number of most variables cytobands to plot
#'
#' @return Plot the gains/lost in the selected cells of interest as multiple
#' barplots
#' 
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_CNA(scExp,  control_samples = unique(scExp$sample_id)[1],
#' ref_genome="hg38", quantiles_to_define_gol = c(0.05,0.95))
#' plot_gain_or_loss_barplots(scExp, cells = scExp$cell_id[which(
#' scExp$sample_id %in% unique(scExp$sample_id)[2])])
#'
plot_gain_or_loss_barplots <- function(scExp, cells = NULL, top = 20){
    gain_or_loss_cytoBand = reducedDim(scExp, "gainOrLoss_cytoBand")
    if(!is.null(cells)){
        gain_or_loss_cytoBand = gain_or_loss_cytoBand[cells,]
    }
    top_variable_cyto = get_most_variable_cyto(scExp, top = top)
    gain_or_loss_cytoBand = as.data.frame(
        t(gain_or_loss_cytoBand[,top_variable_cyto$cytoBand]))
    gain_or_loss_cytoBand = tibble::rownames_to_column(gain_or_loss_cytoBand, var = "cytoBand")
    tab = gain_or_loss_cytoBand  %>% tidyr::gather( "Cell", "Gain_or_Loss", -1)
    tab$Gain_or_Loss = factor(tab$Gain_or_Loss, levels = c(-1,0,1))
    tab$cytoBand = factor(tab$cytoBand, levels = top_variable_cyto$cytoBand)
    p <- ggplot(tab %>% dplyr::group_by(.data[["cytoBand"]],
                                        .data[["Gain_or_Loss"]]) %>%
                    summarise(ncells = n()),
                aes(x = Gain_or_Loss, y = ncells)) + facet_wrap(vars(cytoBand)) +
        geom_bar(aes(fill = Gain_or_Loss), stat ="identity") +
        theme_classic()
    return(p + scale_fill_manual(values = c("#3429FF","#EBEBEB","#F03939"))) 
    
}



#' @title Plot Gain or Loss of cytobands of the most variables cytobands 
#' 
#' @param scExp A SingleCellExperiment with "logRatio_cytoBand" reducedDim slot 
#' filled. See  \code{\link{calculate_logRatio_CNA}}
#' @param cells Cell IDs of the tumor samples to 
#' @param top Number of most variables cytobands to plot
#'
#' @return Plot the gains/lost in the selected cells of interest as multiple
#' barplots
#' 
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_CNA(scExp,  control_samples = unique(scExp$sample_id)[1],
#' ref_genome="hg38", quantiles_to_define_gol = c(0.05,0.95))
#' plot_gain_or_loss_barplots(scExp, cells = scExp$cell_id[which(
#' scExp$sample_id %in% unique(scExp$sample_id)[2])])
#'
plot_gain_or_loss_barplots <- function(scExp, cells = NULL, top = 20){
    gain_or_loss_cytoBand = reducedDim(scExp, "gainOrLoss_cytoBand")
    if(!is.null(cells)){
        gain_or_loss_cytoBand = gain_or_loss_cytoBand[cells,]
    }
    top_variable_cyto = get_most_variable_cyto(scExp, top = top)
    gain_or_loss_cytoBand = as.data.frame(
        t(gain_or_loss_cytoBand[,top_variable_cyto$cytoBand]))
    gain_or_loss_cytoBand = tibble::rownames_to_column(gain_or_loss_cytoBand, var = "cytoBand")
    tab = gain_or_loss_cytoBand  %>% tidyr::gather( "Cell", "Gain_or_Loss", -1)
    tab$Gain_or_Loss = factor(tab$Gain_or_Loss, levels = c(-1,0,1))
    tab$cytoBand = factor(tab$cytoBand, levels = top_variable_cyto$cytoBand)
    p <- ggplot(tab %>% dplyr::group_by(.data[["cytoBand"]],
                                        .data[["Gain_or_Loss"]]) %>%
                    summarise(ncells = n()),
                aes(x = Gain_or_Loss, y = ncells)) + facet_wrap(vars(cytoBand)) +
        geom_bar(aes(fill = Gain_or_Loss), stat ="identity") +
        theme_classic()
    return(p + scale_fill_manual(values = c("#3429FF","#EBEBEB","#F03939"))) 
    
}


#' @title Plot UMAP colored by Gain or Loss of cytobands
#' 
#' @param scExp A SingleCellExperiment with "gainOrLoss_cytoBand" reducedDim 
#' slot filled. See  \code{\link{calculate_gain_or_loss}}
#' @param cytoBand Which cytoBand to color cells by
#'
#' @return Plot the gains/lost of the cytoband overlayed on the epigenetic UMAP.
#' 
#' @export
#' @examples 
#' 
#' data("scExp")
#' scExp = calculate_CNA(scExp,  control_samples = unique(scExp$sample_id)[1],
#' ref_genome="hg38", quantiles_to_define_gol = c(0.05,0.95))
#' plot_reduced_dim_scExp_CNA(scExp, get_most_variable_cyto(scExp)$cytoBand[1])
#'
plot_reduced_dim_scExp_CNA <- function(scExp, cytoBand){
    gain_or_loss_cytoBand = reducedDim(scExp, "gainOrLoss_cytoBand")
    SummarizedExperiment::colData(scExp)[,paste0("Gain_or_Loss_",cytoBand)] = 
        gain_or_loss_cytoBand[,cytoBand]
    scExp = colors_scExp(scExp, annotCol = paste0("Gain_or_Loss_",cytoBand))
    plot_df = as.data.frame(cbind(reducedDim(scExp, "UMAP"), colData(scExp)))
    
    p <- ggplot(plot_df, aes_string(x = "Component_1", y = "Component_2")) +
        geom_point(alpha = 0.35, size = 2, aes(
            color = colData(scExp)[, paste0("Gain_or_Loss_",cytoBand)])) +
        labs(color = paste0("Gain_or_Loss_",cytoBand,"_color")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill = NA)) +
        labs(color = "Gain or Loss") + ggtitle(cytoBand)
    abs_max = max(abs(colData(scExp)[, paste0("Gain_or_Loss_",cytoBand)]))
    return(p + scale_color_gradientn(colours = c("#3429FF","#C5C3E3","#EBEBEB","#DEBFBF", "#F03939"), limits=c(-abs_max, abs_max)))
}