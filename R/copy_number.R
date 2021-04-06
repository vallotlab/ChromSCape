
get_cyto_features <- function(scExp, ref_genome = c("hg38", "mm10")[1] ){
    
    canonical_chr <- eval(parse(text = paste0("ChromSCape::",
                                              ref_genome, ".chromosomes")))
    canonical_chr$start = 1
    canonical_chr <- as(canonical_chr, "GRanges")
    
    cyto <- eval(parse(text = paste0("ChromSCape::",
                                     ref_genome, ".cytoBand")))
        
    cyto = as(cyto,"GRanges")
    
    # Get genomic bins
    bins = rowRanges(scExp)
    bins = bins[which(as.character(seqnames(bins)) %in% as.character(seqnames(canonical_chr)))]
    
    matching_cyto = data.frame(ID = bins$ID, hits = "")
    
    # Get which bins goes into which cyto band
    for(i in unique(cyto$cytoBand)){
        cyto. = cyto[which(cyto$cytoBand == i)]
        matching_cyto$hits[queryHits(findOverlaps(bins,cyto.))] = i
    }
    
    matching_cyto = matching_cyto[which(matching_cyto$hits!=""),]
    return(matching_cyto)
}

calculate_cyto_mat <- function(scExp, ref_genome = c("hg38","mm10")){
    matching_cyto = get_cyto_features(scExp,ref_genome)
    mat = counts(scExp)[match(matching_cyto$ID, rownames(scExp)),]
    gc()
    
    # Group counts by cytobands
    cyto_mat = as.data.frame(as.matrix(mat)) %>% dplyr::group_by(matching_cyto$hits) %>%
        dplyr::summarise(across(.cols = dplyr::everything(),sum)) %>% 
        as.data.frame()
    gc()

    # Normalize by library size
    cyto_mat. = cyto_mat[,-1]
    rownames(cyto_mat) = cyto_mat$`matching_cyto$hits`
    cyto_mat[,-1] = 100* sweep(cyto_mat.,2,colSums(cyto_mat.),`/`)
    cyto_mat = cyto_mat[,-1]
    
    SingleCellExperiment::reducedDim(scExp, "cytoBand") = t(cyto_mat)
    return(scExp)
}

get_most_variable_cyto <- function(scExp, top = 50){
    # Retrieve top variable cytobands
    cyto_df = as.data.frame(t(SingleCellExperiment::reducedDim(scExp, "cytoBand")))
    cyto_df$cytoBand = rownames(cyto_df)
    
    cyto_df = cyto_df %>% tidyr::gather("cell_id","Fri_cyto",-ncol(cyto_df))
    top_var = cyto_df %>% dplyr::group_by(cytoBand) %>% dplyr::summarise(var = var(Fri_cyto)) %>% 
        dplyr::arrange(dplyr::desc(var)) %>% head(n = top)
    
    return(top_var)
}

calculate_logRatio_CNA <- function(scExp, controls){
    cyto_mat = SingleCellExperiment::reducedDim(scExp, "cytoBand")
    
    is_cell_id = ifelse(all(controls %in% scExp$cell_id), TRUE, FALSE)
    if(is_cell_id){
        normal_signal_per_chr = setNames(colMeans(cyto_mat[controls,]),
                                         colnames(cyto_mat[controls,]))
    } else{
        controls_cell = scExp$cell_id[which(scExp$sample_id %in% controls)]
        normal_signal_per_chr = setNames(colMeans(cyto_mat[controls_cell,]),
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

calculate_gain_or_loss <- function(scExp, controls, quantiles = c(0.01, 0.99)){
    logRatio_mat = SingleCellExperiment::reducedDim(scExp, "logRatio_cytoBand")
    
    is_cell_id = ifelse(all(controls %in% scExp$cell_id), TRUE, FALSE)
    if(!is_cell_id){
        controls = scExp$cell_id[which(scExp$sample_id %in% controls)]
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


calculate_CNA <- function(scExp, control_samples = unique(scExp$sample_id)[1],
                          ref_genome = c("hg38","mm10")[1]){
    scExp = calculate_cyto_mat(scExp, ref_genome = ref_genome)
    scExp = calculate_logRatio_CNA(scExp, controls = control_samples)
    scExp = calculate_gain_or_loss(scExp, controls = control_samples)
    return(scExp)
}

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