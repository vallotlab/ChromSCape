## Authors : Pacôme Prompsy, Celine Vallot
## Title : Wrappers & functions to run differential
## analysis on clusters of single-cell based on epigenomic profiling and Gene
## Set Enrichment Analysis on the genes associated with differentially enriched
## features

#' Runs differential analysis between cell clusters
#'
#' Based on clusters of cell defined previously, runs non-parametric Wilcoxon
#' Rank Sum test to find significantly depleted or enriched features, in
#' 'one_vs_rest' mode or 'pairwise' mode. In pairwise mode, each cluster is
#' compared to all other cluster individually, and then pairwise comparisons
#' between clusters are combined to find overall differential features using
#' combineMarkers function from scran.
#'
#' This functions takes as input a SingleCellExperiment object with consclust,
#' the type of comparison, either 'one_vs_rest' or 'pairwise', the adjusted
#' p-value threshold (qval.th) and the fold-change threshold (logFC.th). It
#' outputs a SingleCellExperiment object containing a differential list.
#'
#' @param scExp A SingleCellExperiment object containing consclust with selected
#'   number of cluster.
#' @param by = A character specifying the column of the object containing
#' the groups of cells to compare. Exclusive with de_type == custom
#' @param de_type Type of comparisons. Either 'one_vs_rest', to compare each
#'   cluster against all others, or 'pairwise' to make 1 to 1 comparisons.
#'   ('one_vs_rest')
#' @param method Differential testing method, either 'wilcox' for Wilcoxon non-
#' parametric testing or 'neg.binomial' for edgerGLM based testing. ('wilcox')
#' @param block Use batches as blocking factors ? If TRUE, block will be taken
#' as the column "batch_id" from the SCE. Cells will be compared only within
#' samples belonging to the same batch. 
#' @param group If de_type = "custom", the sample / cluster of interest as a 
#' one- column data.frame. The name of the column is the group name and the
#'  values are character either cluster ("C1", "C2", ...) or sample_id.
#' @param ref If de_type = "custom", the sample / cluster of reference as a one-
#' column data.frame. The name of the column is the group name and the values
#' are character either cluster ("C1", "C2", ...) or sample_id.
#' @param prioritize_genes First filter by loci being close to genes ? E.g. for
#' differential analysis, it is more relevant to keep features close to genes
#' @param max_distanceToTSS If prioritize_genes is TRUE, the maximum distance to 
#' consider a feature close to a gene.
#' @param progress A shiny Progress instance to display progress bar. 
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return Returns a SingleCellExperiment object containing a differential list.
#' @export
#'
#' @importFrom Matrix rowMeans
#' @importFrom scran combineMarkers
#' @importFrom SingleCellExperiment colData normcounts rowData
#' @importFrom rlist list.append
#'
#' @examples
#' data("scExp")
#' scExp_cf = differential_analysis_scExp(scExp)
#' 
differential_analysis_scExp = function(
    scExp, de_type = c("one_vs_rest_fast", "one_vs_rest", "pairwise")[1],
    by = "cell_cluster",
    method = "wilcox", block = NULL, group = NULL, ref = NULL,
    prioritize_genes = nrow(scExp) > 20000, max_distanceToTSS = 1000, 
    progress = NULL, BPPARAM = BiocParallel::bpparam())
{
    warning_DA(scExp, by, de_type, method, block, group, ref)
    if (!is.null(progress)) progress$set(detail = "Retrieving counts...", value = 0.15)
    if (isFALSE(block)) block = NULL
    if(prioritize_genes){
        message("ChromSCape::differential_analysis_scExp - Large number of loci",
                " detected, selecting features closest to ", max_distanceToTSS,
                "bp from genes TSS only.")
        scExp = find_top_features(
            scExp,
            n =  nrow(scExp),
            keep_others = TRUE,
            prioritize_genes = TRUE,
            max_distanceToTSS = max_distanceToTSS)
    }
    if(de_type == "one_vs_rest_fast"){
      res = differential_activation(scExp, group_by = by,
                                    progress = progress)
    } else {

      if(method == "wilcox"){
        counts = SingleCellExperiment::normcounts(
          scExp[SingleCellExperiment::rowData(scExp)$top_feature,])
      } else{
        counts = SingleCellExperiment::counts(
          scExp[SingleCellExperiment::rowData(scExp)$top_feature,])
      }
      feature <- as.data.frame(SummarizedExperiment::rowRanges(
        scExp[SingleCellExperiment::rowData(scExp)$top_feature,]))
      feature = data.frame(ID = feature[, "ID"], chr = feature[, "seqnames"],
                           start = feature[, "start"], end = feature[, "end"])
      affectation = as.data.frame(SingleCellExperiment::colData(scExp))
      
      if (de_type == "one_vs_rest"){
        res <- DA_one_vs_rest(affectation, by, counts,
                              method, feature, block, progress,
                              BPPARAM = BPPARAM)
      } else if(de_type == "pairwise"){
        res <- DA_pairwise(affectation, by, counts, method, feature, block,
                           progress, BPPARAM = BPPARAM)
      } else if(de_type == "custom"){
        res <- DA_custom(affectation, counts, method, feature, block,
                         ref, group, progress, BPPARAM = BPPARAM)
      }
    }
    if (!is.null(progress)) progress$inc(
        detail = paste0("Generating differential table..."), amount = 0.1)
  
     rowD <- SummarizedExperiment::rowData(scExp)
     rowD = rowD[,grep("logFC.|group_activation.|reference_activation.|qval.",colnames(rowD), invert = TRUE)]
     rowD =  dplyr::left_join(as.data.frame(rowD),
                       res[,-match(c("chr", "start", "end"),colnames(res))],
                       by = c("ID")) 
     SummarizedExperiment::rowData(scExp) <- rowD
     scExp@metadata$DA_parameters$de_type = de_type
    return(scExp)
}

#' Warning for differential_analysis_scExp
#'
#' @param scExp A SingleCellExperiment object containing consclust with selected
#'   number of cluster.
#' @param by = A character specifying the column of the object containing
#' the groups of cells to compare. Exclusive with de_type == custom
#' @param de_type Type of comparisons. Either 'one_vs_rest', to compare each
#'   cluster against all others, or 'pairwise' to make 1 to 1 comparisons.
#'   ('one_vs_rest')
#' @param method Wilcoxon or edgerGLM
#' @param block Use batches as blocking factors ?
#' @param group If de_type is custom, the group to compare (data.frame), must
#' be a one-column data.frame with cell_clusters or sample_id as character in 
#' rows
#' @param ref If de_type is custom, the reference to compare (data.frame), must
#' be a one-column data.frame with cell_clusters or sample_id as character in 
#' rows
#'
#' @return Warnings or Errors if the input are not correct
warning_DA <- function(scExp, by, de_type, method, block, group, ref){
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(de_type))
  
    if (!de_type %in% c("one_vs_rest_fast","one_vs_rest", "pairwise", "custom")) 
        stop("ChromSCape::run_differential_analysis_scExp - de_type ",
                    "must be 'one_vs_rest', 'pairwise' 'custom'.")
    if (!method %in% c("wilcox", "neg.binomial")) 
        stop("ChromSCape::run_differential_analysis_scExp - method must",
                    "be 'wilcox' or 'neg.binomial'.")
    if (!by %in% colnames(SingleCellExperiment::colData(scExp))
        & de_type != "custom") 
        stop("ChromSCape::run_differential_analysis_scExp - scExp ",
                    "object must have selected number of clusters.")
    if (FALSE %in% (c("ID", "Gene", "distanceToTSS") %in% 
                    colnames(SingleCellExperiment::rowData(scExp)))) 
        stop("ChromSCape::run_differential_analysis_scExp - Please run ",
                    "feature_annotation_scExp first.")
    if(de_type == "custom"){
        stopifnot(is.data.frame(group), is.data.frame(ref))
    }
}

#' Summary of the differential analysis
#'
#' @param scExp A SingleCellExperiment object containing consclust with selected
#'   number of cluster. 
#' @param qval.th Adjusted p-value threshold. (0.01)
#' @param logFC.th Fold change threshold. (1)
#' @param min.percent Minimum fraction of cells having the feature active to
#' consider it as significantly differential. (0.01)
#'
#' @return A table summary of the differential analysis
#' @export
#'
#' @examples
#' data('scExp')
#' summary_DA(scExp)
summary_DA <- function(scExp, qval.th = 0.01, 
                       logFC.th = 1, min.percent = 0.01){
  stopifnot(!is.null(scExp))
  
    res = as.data.frame(SingleCellExperiment::rowData(scExp))
  groups = gsub(".*\\.","",grep("qval",
                                colnames(SingleCellExperiment::rowData(scExp)),
                                value = TRUE))
  summary = matrix(nrow = 3, ncol = length(groups), dimnames = list(
    c("differential", "over", "under"), groups))
  for (k in seq_along(groups)){
    gpsamp = groups[k]
    qval.col <- paste("qval", gpsamp, sep = ".")
    logFC.col <- paste("logFC", gpsamp, sep = ".")
    
    
    summary["over", gpsamp] = sum(res[,qval.col] <= qval.th & 
                                    res[,logFC.col] > logFC.th &
                                    res[,paste("group_activation", gpsamp, sep = ".")] > min.percent, na.rm = TRUE)
    summary["under", gpsamp] = sum(res[,qval.col] <= qval.th & 
                                     res[,logFC.col] < -logFC.th &
                                     res[,paste("reference_activation", gpsamp, sep = ".")] > min.percent, na.rm = TRUE)
    summary["differential", gpsamp] = summary["under", gpsamp] +
      summary["over", gpsamp]
  }
  
  return(summary)
}


#' Differential Analysis in 'One vs Rest' mode
#'
#' @param affectation An annotation data.frame with cell_id and cell_cluster
#'   columns
#' @param by = A character specifying the column of the object containing
#' the groups of cells to compare.
#' @param counts Count matrix
#' @param method DA method : Wilcoxon or EdgeR
#' @param feature Feature tables
#' @param block Blocking feature
#' @param progress A shiny Progress instance to display progress bar. 
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return A list of results, groups compared and references
#'   
DA_one_vs_rest <- function(affectation, by, counts, method, feature,
                            block, progress = NULL,
                            BPPARAM = BiocParallel::bpparam()){
    stopifnot(is.data.frame(affectation),
                is(counts,"dgCMatrix")|is.matrix(counts), is.character(method),
                is.data.frame(feature))
    groups = unique(affectation[,by])
    # compare each cluster to all the rest
    mygps = lapply(groups, function(i)
    {
        affectation[which(affectation[,by] == i), "cell_id"]
    })
    names(mygps) = groups
    
    myrefs =  lapply(groups, function(i)
    {
        affectation[which(affectation[,by] != i), "cell_id"]
    })
    
    names(myrefs) = paste0("not_", groups)
    refs = names(myrefs)
    
    if (!is.null(progress)) progress$inc(
        detail = paste0("Comparing Group vs Refs"), amount = 0.2)
    if(method == "wilcox"){ 
      res = CompareWilcox(
        dataMat = counts, annot = affectation, ref_group = myrefs, 
        groups = mygps, featureTab = feature, block = block,
        BPPARAM = BPPARAM)
    } else {
        res = CompareedgeRGLM( dataMat = counts, 
                            annot = affectation, ref_group = myrefs,
                            groups = mygps, featureTab = feature)
    }
    
    return(res)
}

#' Run differential analysis in Pairwise mode
#'
#' @param affectation An annotation data.frame with cell_cluster and cell_id 
#' columns
#' @param by = A character specifying the column of the object containing
#' the groups of cells to compare.
#' @param counts Count matrix
#' @param feature Feature data.frame
#' @param method DA method, Wilcoxon or edgeR
#' @param block Blocking feature
#' @param progress A shiny Progress instance to display progress bar. 
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return A list of results, groups compared and references
#'
DA_pairwise <- function(affectation, by, counts,
                        method, feature, block, progress = NULL,
                        BPPARAM = BiocParallel::bpparam()){
    stopifnot(is.data.frame(affectation),
                is(counts,"dgCMatrix")|is.matrix(counts), is.character(method),
                is.data.frame(feature))
    res = feature
    out <- run_pairwise_tests(affectation, by, counts, feature, method,
                              progress = progress, BPPARAM = BPPARAM)
    if (!is.null(progress)) progress$inc(detail = "Merging results...",
                                         amount = 0.1)
    count_save <- out$count_save
    pairs <- out$pairs
    single_results <- out$single_results
    all_groups = unique(affectation[, by])
    tmp.gp = list(affectation[which(
        affectation[, by] == all_groups[length(all_groups)]), "cell_id"])
    count_save[, all_groups[length(all_groups)]] = apply(
        counts, 1, function(x) mean(x[as.character(tmp.gp[[1]])]))
    
    names = c("ID",  "chr", "start", "end", "logFC",
              "qval", "group_activation", "reference_activation")
    single_results = lapply(single_results, setNames, names)

    combinedTests = scran::combineMarkers(
        de.lists = single_results, pairs = pairs, pval.field = "qval",
        effect.field = "logFC", pval.type = "any", log.p.in = FALSE, 
        log.p.out = FALSE, output.field = "stats", full.stats = TRUE,
        sorted = FALSE)
    
    for (i in all_groups)
    {
        logFCs = vapply(seq_len(all_groups[-length(all_groups)]), function(k){
            combinedTests[[i]][,k + 4]$logFC
        }, FUN.VALUE = rep(0,nrow(feature)))
        group_activation = vapply(seq_len(all_groups[-length(all_groups)]), function(k){
          combinedTests[[i]][,k + 4]$group_activation
        }, FUN.VALUE = rep(0,nrow(feature)))
        reference_activation = vapply(seq_len(all_groups[-length(all_groups)]), function(k){
          combinedTests[[i]][,k + 4]$reference_activation
        }, FUN.VALUE = rep(0,nrow(feature)))
        res[, paste0("logFC.", i)] = Matrix::rowMeans(logFCs)
        res[, paste0("qval.", i)] =  combinedTests[[i]]$p.value
        res[, paste0("group_activation.", i)] = Matrix::rowMeans(group_activation)
        res[, paste0("reference_activation.", i)] = Matrix::rowMeans(reference_activation)
        
    }
    return(res)
}

#' Differential Analysis in 'One vs Rest' mode
#'
#' @param affectation An annotation data.frame with cell_id and cell_cluster
#'   columns
#' @param counts Count matrix
#' @param method DA method : Wilcoxon or EdgeR
#' @param feature Feature tables
#' @param block Blocking feature
#' @param group If de_type is custom, the group to compare (data.frame), must
#' be a one-column data.frame with cell_clusters or sample_id as character in 
#' rows
#' @param ref If de_type is custom, the reference to compare (data.frame), must
#' be a one-column data.frame with cell_clusters or sample_id as character in 
#' rows
#' @param progress A shiny Progress instance to display progress bar. 
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return A list of results, groups compared and references
#'   
DA_custom <- function(affectation, counts, method, feature,
                      block, ref, group, progress = NULL,
                      BPPARAM = BiocParallel::bpparam()){
    stopifnot(is.data.frame(affectation),
              is(counts,"dgCMatrix")|is.matrix(counts), is.character(method),
              is.data.frame(feature), is.data.frame(ref), is.data.frame(group))
    if (!is.null(progress)) progress$inc(detail = "Retrieving group and ref...",
                                         amount = 0.1)
    # compare each cluster to all the rest
    groups = colnames(group)
    mygps = unique(c(affectation[which(affectation$cell_cluster %in% group[,1]
                        | affectation$sample_id %in% group[,1]), "cell_id"]))
    mygps = list( "a" = mygps)
    names(mygps) = groups
    refs = colnames(ref)
    myrefs = unique(c(affectation[which(affectation$cell_cluster %in% ref[,1]
                        | affectation$sample_id %in% ref[,1]), "cell_id"]))
    myrefs = list("a" = myrefs)
    names(myrefs) = refs

    if (!is.null(progress)) progress$inc(
        detail = paste0("Comparing ", groups, " vs ", refs), amount = 0.2)
    if(method == "wilcox"){
        res = CompareWilcox(
        dataMat = counts, annot = affectation, ref_group = myrefs, 
        groups = mygps, featureTab = feature, block = block, BPPARAM=BPPARAM)
    } else {
        res = CompareedgeRGLM( dataMat = counts, 
                               annot = affectation, ref_group = myrefs,
                               groups = mygps, featureTab = feature)
    }
    return(res)
}

#' Run pairwise tests
#'
#' @param affectation An annotation data.frame with cell_cluster and cell_id 
#' columns
#' @param by = A character specifying the column of the object containing
#' the groups of cells to compare.
#' @param counts Count matrix
#' @param feature Feature data.frame
#' @param method DA method, Wilcoxon or edgeR
#' @param progress A shiny Progress instance to display progress bar. 
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return A list containing objects for DA function
#'
run_pairwise_tests <- function(affectation, by, counts, 
                            feature, method, progress = NULL,
                            BPPARAM = BiocParallel::bpparam()){
    stopifnot(is.data.frame(affectation),
                is(counts,"dgCMatrix")|is.matrix(counts), is.character(method))
    count_save = data.frame(ID = feature$ID)
    single_results = list()
    pairs = setNames(data.frame(matrix(ncol = 2, nrow = 0)),
                    c("group", "ref"))
    all_groups = unique(affectation[,by])
    for (i in all_groups){
        mygps = list(affectation[which(affectation[,by] == i), "cell_id"])
        names(mygps) = i
        groups = names(mygps)
        for (j in setdiff(all_groups, groups)){
            myrefs = list(affectation[which(
                affectation[,by] == j), "cell_id"])
            names(myrefs) = j
            refs = names(myrefs)
            if (!is.null(progress)) progress$set(
                detail = paste0(groups, " vs ",refs),
                value = 0.15 + (0.85 / ((length(all_groups) - 1) * (length(all_groups) - 1) )))
            if(method == "wilcox"){ tmp_result = CompareWilcox(
                dataMat = counts, annot = affectation, ref_group = myrefs, 
                groups = mygps, featureTab = feature, BPPARAM =BPPARAM)
            } else {
                tmp_result = CompareedgeRGLM(
                    dataMat = counts, annot = affectation, ref_group = myrefs,
                    groups = mygps, featureTab = feature)
                }
            tmp_result = tmp_result[match(count_save$ID, tmp_result$ID), ]
            rownames(tmp_result) = tmp_result$ID
          
            single_results = rlist::list.append(single_results, tmp_result)
            pairs[nrow(pairs) + 1, ] = list(groups[1], refs[1])
            tmp_mirror = tmp_result
            tmp_mirror[, paste0("logFC.", groups)] = tmp_result[,paste0("logFC.", groups)] * (-1)
            tmp_mirror[, paste0("group_activation.", groups)] = tmp_result[, paste0("reference_activation.", groups)]
            tmp_mirror[, paste0("reference_activation.", groups)] =  tmp_result[, paste0("group_activation.", groups)]
            single_results = rlist::list.append(single_results, tmp_mirror)
            pairs[nrow(pairs) + 1, ] = list(refs[1], groups[1])}}
    out <- list("single_results" = single_results, "pairs" = pairs,
                "count_save" = count_save)
    return(out)}

#' Find Differentialy Activated Features (One vs All)
#' 
#' @description 
#' 
#' Based on the statement that single-cell epigenomic dataset are very sparse,
#' specifically when analysis small bins or peaks, we can define each feature as
#' being 'active' or not simply by the presence or the absence of reads in this 
#' feature. This is the equivalent of binarize the data. When trying to find 
#' differences in signal for a feature between multiple cell groups, this 
#' function simply compare the percentage of cells 'activating' the feature 
#' in each of the group. The p.values are then calculated using a Pearson's 
#' Chi-squared Test for Count Data (comparing the number of active cells in one 
#' group vs the other) and corrected using Benjamini-Hochberg correction for 
#' multiple testing. 
#' 
#' @param scExp  A SingleCellExperiment object containing consclust with selected
#'   number of cluster.
#' @param by Which grouping to run the marker enrichment ?
#' @param progress A shiny Progress instance to display progress bar. 
#' @param verbose Print ? 
#'  
#' @details 
#' To calculate the logFC, the percentage of activation of the features are 
#' corrected for total number of reads to correct for library size bias.  
#' For each cluster ('group') the function consider the rest of the cells as
#' the reference.
#' @seealso
#' For Pearson's Chi-squared Test for Count Data  \link[stats]{chisq.test}.
#' For other differential analysis see \link[ChromSCape]{differential_analysis_scExp}.
#' 
#' @return Returns a dataframe of differential activation results that contains
#' the rowData of the SingleCellExperiment with additional logFC, q.value, 
#' group activation (fraction of cells active for each feature in the group cells),
#' reference activation (fraction of cells active for each feature in the 
#' reference cells).
#' @export
#'
#' @importFrom Matrix rowSums colSums
#' @importFrom SingleCellExperiment colData counts rowData
#' @importFrom rlist list.append
#'
#' @examples
#' data("scExp")
#' res = differential_activation(scExp, by = "cell_cluster")
#' res = differential_activation(scExp, by = "sample_id")
differential_activation <- function(scExp, by = c("cell_cluster","sample_id")[1],
                                    verbose = TRUE, progress = NULL){
  .par_chisq = function(row) { return(chisq.test(
    matrix(c(row[1], row[2], row[3], row[4]), ncol = 2),
    simulate.p.value = FALSE)$p.value)}
  
  groups = unique(colData(scExp)[, by])
  list_res = list()
  
  mat = Matrix::Matrix(SingleCellExperiment::counts(
    scExp[SingleCellExperiment::rowData(scExp)$top_feature,]) > 0 + 0, sparse = TRUE)
  
  feature <- as.data.frame(SummarizedExperiment::rowRanges(
    scExp[SingleCellExperiment::rowData(scExp)$top_feature,]))
  feature = data.frame(ID = feature[, "ID"], chr = feature[, "seqnames"],
                       start = feature[, "start"], end = feature[, "end"])
  
  
  for(group in groups){
    if (!is.null(progress)) progress$inc(
      detail = paste0("Calculating differential activation - ", group, "..."), amount = 0.9/length(groups))
    if(verbose) cat("ChromSCape::differential_activation - Calculating differential activation for", group,".\n")
    cluster_bin_mat = mat[,which(SingleCellExperiment::colData(scExp)[, by] %in% group)]
    reference_bin_mat = mat[,which(!SingleCellExperiment::colData(scExp)[, by] %in% group)]
    
    rectifier = mean(Matrix::colSums(cluster_bin_mat)) / mean(Matrix::colSums(reference_bin_mat))
    group_sum = Matrix::rowSums(cluster_bin_mat)
    group_activation = group_sum / ncol(cluster_bin_mat)
    group_corrected_activation =  group_activation / rectifier
    
    reference_sum = Matrix::rowSums(reference_bin_mat) 
    reference_activation = reference_sum / ncol(reference_bin_mat)
    
    n_cell_cluster = ncol(cluster_bin_mat)
    n_cell_reference = ncol(reference_bin_mat)
    other_group =  n_cell_cluster - group_sum
    other_ref = n_cell_reference - reference_sum
    
    chisq_mat = cbind(group_sum, other_group, reference_sum, other_ref)

    suppressWarnings({pvalues = apply(chisq_mat, 1, .par_chisq)})
    
    q.values = p.adjust(pvalues, method = "BH")
    logFCs = log2(group_corrected_activation/reference_activation)
    if(any(is.nan(logFCs))) logFCs[which(is.nan(logFCs))] = 0
    if(any(logFCs == Inf)) logFCs[which(logFCs == Inf)] = max(
      logFCs[which(!is.infinite(logFCs))])
    if(any(logFCs == -Inf)) logFCs[which(logFCs == -Inf)] = min(
      logFCs[which(!is.infinite(logFCs))])
    
    res = data.frame(
      logFC.gpsamp = logFCs,
      qval.gpsamp = q.values,
      group_activation.gpsamp = group_activation,
      reference_activation.gpsamp = reference_activation
    )
    colnames(res) = gsub("gpsamp",  group, colnames(res))
    list_res[[group]] = res
  }
  gc()
  
  names(list_res) = NULL
  res =  cbind(feature, do.call("cbind", list_res))
  
  return(res)
}

#' Runs Gene Set Enrichment Analysis on genes associated with differential
#' features
#'
#' This function takes previously calculated differential features and runs
#' hypergeometric test to look for enriched gene sets in the genes associated
#' with differential features, for each cell cluster. This functions takes as
#' input a SingleCellExperiment object with consclust, the type of comparison,
#' either 'one_vs_rest' or 'pairwise', the adjusted p-value threshold (qval.th)
#' and the fold-change threshold (logFC.th). It outputs a SingleCellExperiment
#' object containing a differential list.
#'
#' @param scExp A SingleCellExperiment object containing list of differential
#'   features.
#' @param ref A reference annotation, either 'hg38' or 'mm10'. ('hg38')
#' @param enrichment_qval Adjusted p-value threshold for gene set enrichment.
#'   (0.1)
#' @param GeneSets A named list of gene sets. If NULL will automatically load
#'   MSigDB list of gene sets for specified reference genome. (NULL)
#' @param GeneSetsDf A dataframe containing gene sets & class of gene sets. If
#'   NULL will automatically load MSigDB dataframe of gene sets for specified
#'   reference genome. (NULL)
#' @param GenePool The pool of genes to run enrichment in. If NULL will
#'   automatically load Gencode list of genes fro specified reference genome.
#'   (NULL)
#' @param qval.th Adjusted p-value threshold to define differential features.
#'   (0.01)
#' @param logFC.th Fold change threshold to define differential features. (1)
#' @param min.percent Minimum fraction of cells having the feature active to
#' consider it as significantly differential. (0.01)
#' @param peak_distance Maximum distanceToTSS of feature to gene TSS to consider
#'   associated, in bp. (1000)
#' @param use_peaks Use peak calling method (must be calculated beforehand).
#'   (FALSE)
#' @param GeneSetClasses Which classes of MSIGdb to look for.  
#' @param progress A shiny Progress instance to display progress bar. 
#' 
#' @return Returns a SingleCellExperiment object containing list of enriched
#'   Gene Sets for each cluster, either in depleted features, enriched features
#'   or simply differential features (both).
#'
#' @export
#' @importFrom SingleCellExperiment colData normcounts rowData
#' @importFrom msigdbr msigdbr
#'
#' @examples
#' data("scExp")
#' 
#' #Usually recommanding qval.th = 0.01 & logFC.th = 1 or 2
#' \dontrun{scExp_cf = gene_set_enrichment_analysis_scExp(scExp,
#'  qval.th = 0.4, logFC.th = 0.3)}
#' 
gene_set_enrichment_analysis_scExp = function(
    scExp, enrichment_qval = 0.1, ref = "hg38", GeneSets = NULL,
    GeneSetsDf = NULL, GenePool = NULL, qval.th = 0.01, logFC.th = 1,
    min.percent = 0.01, peak_distance = 1000, use_peaks = FALSE,
    GeneSetClasses = c(
        "c1_positional","c2_curated","c3_motif","c4_computational",
        "c5_GO","c6_oncogenic","c7_immunologic","hallmark"), progress=NULL)
{
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(ref),
            is.numeric(enrichment_qval), is.numeric(peak_distance),
            is.numeric(qval.th), is.numeric(logFC.th), is.numeric(min.percent),
            is.character(GeneSetClasses))
    if (!is.null(progress)) progress$inc(
        detail = "Loading reference gene sets...", amount = 0.1)

    if (!any(grepl("qval", colnames(SingleCellExperiment::rowData(scExp))))) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - No DA, ",
        "please run run_differential_analysis_scExp first.")
    if (is.null(SingleCellExperiment::rowData(scExp)$Gene)) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - No Gene ",
                    "annotation, please annotate features with genes using ",
                    "feature_annotation_scExp first.")
    if (is.null(use_peaks)) use_peaks = FALSE
    
    if (use_peaks & (!"refined_annotation" %in% names(scExp@metadata))) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - ",
        "When use_peaks is TRUE, metadata must contain refined_annotation ",
        "object.")
    
    refined_annotation = NULL
    if(use_peaks) refined_annotation = scExp@metadata$refined_annotation
    
    database_MSIG <- load_MSIGdb(ref,GeneSetClasses)
    GeneSets = database_MSIG$GeneSets
    GeneSetsDf = database_MSIG$GeneSetsDf
    GenePool = database_MSIG$GenePool
    GeneSets <- lapply(GeneSets, function(x) unique(intersect(x, GenePool)))
    
    annotFeat_long = as.data.frame(tidyr::separate_rows(
        as.data.frame(SummarizedExperiment::rowRanges(scExp)), 
        .data$Gene, sep = ", "))

    enr <- combine_enrichmentTests(
        diff = as.data.frame(SingleCellExperiment::rowData(scExp)),
        enrichment_qval = enrichment_qval, qval.th = qval.th,
        logFC.th = logFC.th, min.percent = min.percent,
        annotFeat_long = annotFeat_long, peak_distance = peak_distance,
        refined_annotation = refined_annotation, GeneSets = GeneSets,
        GeneSetsDf = GeneSetsDf, GenePool = GenePool, progress = progress)
    scExp@metadata$enr <- enr
    return(scExp)
}


#' Load and format MSIGdb pathways using msigdbr package
#'
#' @param ref Reference genome, either mm10 or hg38
#' @param GeneSetClasses Which classes of MSIGdb to load
#'
#' @return A list containing the GeneSet (list), GeneSetDf (data.frame) and 
#' GenePool character vector of all possible genes
#'
load_MSIGdb <- function(ref, GeneSetClasses){
    if ((!ref %in% c("hg38", "mm10")) ) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - ",
                    "Reference genome (ref) must be ",
                "'hg38' or 'mm10' if gene sets not specified.")
    stopifnot(is.character(GeneSetClasses))
    message(
        paste0(
            "ChromSCape::gene_set_enrichment_analysis_scExp - Loading ",
            ref,
            " MSigDB gene sets."
        )
    )
    columns = c("gs_name", "gs_cat", "gene_symbol")
    if (ref == "hg38")
        GeneSetsDf = msigdbr::msigdbr("Homo sapiens")[, columns]
    if (ref == "mm10")
        GeneSetsDf = msigdbr::msigdbr("Mus musculus")[, columns]
    colnames(GeneSetsDf) = c("Gene.Set", "Class", "Genes")
    system.time({
        GeneSetsDf <- GeneSetsDf %>% dplyr::group_by(
            .data$Gene.Set, .data$Class) %>%
            dplyr::summarise("Genes" = paste(.data$Genes,
                                        collapse = ","))
    })
    corres = data.frame(
        long_name = c("c1_positional", "c2_curated", "c3_motif", 
                    "c4_computational", "c5_GO", "c6_oncogenic",
                    "c7_immunologic", "hallmark"), short_name = c(
                        paste0("C", seq_len(7)), "H"))
    GeneSetsDf$Class = corres$long_name[
        match(GeneSetsDf$Class, corres$short_name)]
    GeneSetsDf = GeneSetsDf[which(GeneSetsDf$Class %in% GeneSetClasses),]
    GeneSets = lapply(GeneSetsDf$Gene.Set, function(x) {
        unlist(strsplit(
            GeneSetsDf$Genes[which(GeneSetsDf$Gene.Set == x)], split = ","))})
    names(GeneSets) = GeneSetsDf$Gene.Set
    message(paste0(
        "ChromSCape::gene_set_enrichment_analysis_scExp - Selecting ", 
        ref, " genes from Gencode."))
    eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
    GenePool = eval(parse(text = paste0("", ref, ".GeneTSS")))
    GenePool = unique(as.character(GenePool[, "Gene"]))
    
    database_MSIG <- list("GeneSets" = GeneSets,
                        "GeneSetsDf" = GeneSetsDf,
                        "GenePool" = GenePool)
    return(database_MSIG)
}


#' Resutls of hypergeometric gene set enrichment test
#'
#' Run hypergeometric enrichment test and combine significant pathways into
#' a data.frame
#'
#' @param enrichment_qval Adusted p-value threshold above which a pathway is
#' considered significative
#' @param GeneSets List of pathways
#' @param GeneSetsDf Data.frame of pathways
#' @param GenePool Pool of possible genes for testing
#' @param differentialGenes Genes significantly over / under expressed
#'
#' @return A data.frame with pathways passing q.value threshold
#'
results_enrichmentTest <- function(
    differentialGenes, enrichment_qval, GeneSets, GeneSetsDf, GenePool){
    enrich.test = enrichmentTest(gene.sets = GeneSets,
                                mylist = differentialGenes, 
                                possibleIds = GenePool)
    
    enrich.test = data.frame(Gene_set_name = rownames(enrich.test),
                            enrich.test, 
                            check.names = FALSE)
    
    enrich.test = merge(GeneSetsDf[,c("Gene.Set","Class")],
                        enrich.test, by.x = "Gene.Set", by.y = "Gene_set_name",
                        all.y = TRUE, sort = FALSE)  ## Get class of gene set
    
    enrich.test = enrich.test[order(enrich.test$`p-value`), ]
    ind = which(enrich.test$`q-value` <= enrichment_qval)
    if (!length(ind))
    {
        ind = seq_len(10)
    }
    return(enrich.test[ind, ]) 
}


#' Run enrichment tests and combine into list
#' 
#' @param diff Differential list
#' @param enrichment_qval Adusted p-value threshold above which a pathway is
#' considered significative list
#' @param qval.th Differential analysis adjusted p.value threshold
#' @param logFC.th Differential analysis log-fold change threshold
#' @param min.percent Minimum fraction of cells having the feature active to
#' consider it as significantly differential. (0.01)
#' @param annotFeat_long Long annotation
#' @param peak_distance Maximum gene to peak distance
#' @param refined_annotation Refined annotation data.frame if peak calling is 
#' done
#' @param GeneSets List of pathways
#' @param GeneSetsDf Data.frame of pathways
#' @param GenePool Pool of possible genes for testing
#' @param progress A shiny Progress instance to display progress bar. 
#' 
#' @return A list of list of pathway enrichment data.frames for 
#' Both / Over / Under and for each cluster
#'
combine_enrichmentTests <- function(
    diff, enrichment_qval, qval.th, logFC.th, min.percent,
    annotFeat_long, peak_distance, refined_annotation, GeneSets, GeneSetsDf, 
    GenePool, progress = NULL)
    {
    stopifnot(is.data.frame(diff),
              is.numeric(enrichment_qval), is.numeric(qval.th),
                is.numeric(logFC.th), is.data.frame(annotFeat_long),
                is.numeric(peak_distance), is.list(GeneSets),
                is.data.frame(GeneSetsDf), is.character(GenePool))
    groups <- gsub(".*\\.","", grep("qval",colnames(diff), value = TRUE))
    res <- diff
    enr = list(Both = list(), Overexpressed = list(), Underexpressed = list())
    for (i in seq_along(groups)) {
        gp = groups[i]
        if (!is.null(progress)) progress$set(
            detail = paste0("Enrichment for ",gp, "..."),
            value = 0.3 + (0.5 / length(groups)))
        qval.col <- paste("qval", gp, sep = ".")
        logFC.col <- paste("logFC", gp, sep = ".")
        group_activation.col <- paste("group_activation", gp, sep = ".")
        reference_activation.col <- paste("reference_activation", gp, sep = ".")
        
        over = res$ID[which(res[, qval.col] <= qval.th & 
                                res[, logFC.col] > logFC.th &
                                res[, group_activation.col] > min.percent)]
        overG = unique(
            annotFeat_long$Gene[annotFeat_long$distanceToTSS < peak_distance & 
                                            annotFeat_long$ID %in% over])
        under = res$ID[which(res[, qval.col] <= qval.th &
                                res[,logFC.col] < -logFC.th  &
                                 res[, reference_activation.col] > min.percent)]
        underG = unique(
            annotFeat_long$Gene[annotFeat_long$distanceToTSS < peak_distance & 
                                    annotFeat_long$ID %in% under])
        signific = c(over, under)
        significG = unique(c(overG, underG))
        
        if (!is.null(refined_annotation)){
            out <- filter_genes_with_refined_peak_annotation(
                refined_annotation, peak_distance, signific, over, under)
            significG = out$significG
            overG = out$overG
            underG = out$underG
        }
        if (length(significG)) enr$Both[[i]] = results_enrichmentTest(
            differentialGenes = significG, enrichment_qval, GeneSets,
            GeneSetsDf, GenePool)
        if (length(overG)) enr$Overexpressed[[i]] = results_enrichmentTest(
            differentialGenes = overG, enrichment_qval, GeneSets,
            GeneSetsDf, GenePool)
        if (length(underG)) enr$Underexpressed[[i]] = results_enrichmentTest(
            differentialGenes = underG, enrichment_qval, GeneSets,
            GeneSetsDf, GenePool)
    }
return(enr)
}

#' Filter genes based on peak calling refined annotation
#'
#' @param refined_annotation A data.frame containing each gene distance to real 
#' peak
#' @param peak_distance Minimum distance to an existing peak to accept a given 
#' gene 
#' @param signific Indexes of all significantly differential genes 
#' @param over Indexes of all significantly overexpressed genes
#' @param under Indexes of all significantly underexpressed genes
#'
#' @return List of significantly differential, overexpressed 
#' and underexpressed genes close enough to existing peaks
#' 
filter_genes_with_refined_peak_annotation <- function(
    refined_annotation, peak_distance, signific, over, under){
    w_ids <- refined_annotation$window_ID
    p_ids <- refined_annotation$peak_ID
    genes <- refined_annotation$Gene
    distTSS <- refined_annotation$distance
    signific_associated_peak = p_ids[w_ids %in% signific]
    over_associated_peak = p_ids[w_ids %in% over]
    under_associated_peak = p_ids[w_ids %in% under]
    
    signific_associated_gene = genes[p_ids %in% signific_associated_peak &
                                        distTSS < peak_distance]
    over_associated_gene = genes[p_ids %in% over_associated_peak &
                                    distTSS < peak_distance]
    under_associated_gene = genes[p_ids %in% under_associated_peak &
                                    distTSS < peak_distance]
    significG = unique(signific_associated_gene)
    overG = unique(over_associated_gene)
    underG = unique(under_associated_gene)
    
    return(list("significG" = significG, "overG" = overG, "underG" = underG))
}

#' Creates table of enriched genes sets
#'
#' @param scExp  A SingleCellExperiment object containing list of enriched gene
#'   sets.
#' @param set A character vector, either 'Both', 'Overexpressed' or
#'   'Underexpressed'. ('Both')
#' @param enr_class_sel Which classes of gene sets to show. (c('c1_positional',
#'   'c2_curated', ...))
#' @param group The "group" name from differential analysis. Can be the cluster
#' name or the custom name in case of a custom differential analysis.
#'
#' @return A DT::data.table of enriched gene sets.
#' @export
#'
#' @importFrom DT datatable
#' @importFrom tidyr unite
#' 
#' @examples
#' data("scExp")
#' \dontrun{table_enriched_genes_scExp(scExp)}
#' 
table_enriched_genes_scExp <- function(
    scExp, set = "Both", group = "C1", 
    enr_class_sel = c(
        "c1_positional", "c2_curated", "c3_motif", "c4_computational", "c5_GO",
        "c6_oncogenic", "c7_immunologic", "hallmark"))
{
    stopifnot(is(scExp, "SingleCellExperiment"),
            is.character(set), is.character(group), 
            is.character(enr_class_sel))
    if (is.null(scExp@metadata$enr)) 
        stop("ChromSCape::table_enriched_genes_scExp - No GSEA, please ",
        "run gene_set_enrichment_analysis_scExp first.")
    
    if (!set %in% c("Both", "Overexpressed", "Underexpressed")) 
        stop("ChromSCape::table_enriched_genes_scExp - set variable",
        "must be 'Both', 'Overexpressed' or 'Underexpressed'.")
    
    return_df = setNames(
        data.frame(matrix(ncol = 6, nrow = 0)),
        c("Gene_set", "Class", "Num_deregulated_genes", "p.value",
          "adj.p.value", "Deregulated_genes"))
  groups = gsub(".*\\.","",grep("qval",
                                colnames(SingleCellExperiment::rowData(scExp)),
                                value = TRUE))
    if (!group %in% groups) 
        stop("ChromSCape::table_enriched_genes_scExp - Group is not in ",
        "differential analysis")
    
    n = match(group, groups)
    if(n > length(scExp@metadata$enr[[set]]) || 
       is.null(scExp@metadata$enr[[set]][[n]])) return(return_df)
    
    table <- scExp@metadata$enr[[set]][[n]]
    table <- table[which(table[, "Class"] %in% enr_class_sel), ]
    if (is.null(table))
    {
        return(return_df)
    }
    table <- tidyr::unite(table, "dereg_genes", c(
        "Nb_of_deregulated_genes", "Nb_of_genes"), sep = "/")
    colnames(table) <- c("Gene_set", "Class", "Num_deregulated_genes", 
                        "p.value", "adj.p.value", "Deregulated_genes")
    table[, 4] <- round(table[, 4], 9)
    table[, 5] <- round(table[, 5], 9)
    table <- table[order(table$adj.p.value, table$p.value), ]
    return(table)
}




#' Find the TF that are enriched in the differential genes using ChEA3 database
#'
#' @param scExp A SingleCellExperiment object containing list of differential
#'   features.
#' @param ref A reference annotation, either 'hg38' or 'mm10'. ('hg38')
#' @param qval.th Adjusted p-value threshold to define differential features.
#'   (0.01)
#' @param logFC.th Fold change threshold to define differential features. (1)
#' @param min.percent Minimum fraction of cells having the feature active to
#' consider it as significantly differential. (0.01)
#' @param peak_distance Maximum distanceToTSS of feature to gene TSS to consider
#'   associated, in bp. (1000)
#' @param use_peaks Use peak calling method (must be calculated beforehand).
#'   (FALSE)
#' @param progress A shiny Progress instance to display progress bar. 
#' @param verbose A logical to print message or not. (TRUE)
#' 
#' @return Returns a SingleCellExperiment object containing list of enriched
#'   Gene Sets for each cluster, either in depleted features, enriched features
#'   or simply differential features (both).
#'
#' @export
#' @importFrom SingleCellExperiment colData rowData
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom tidyr separate_rows 
#'
#' @examples
#' data("scExp")
#' 
#' scExp = enrich_TF_ChEA3_scExp(
#'  scExp,
#'  ref = "hg38",
#'  qval.th = 0.01,
#'  logFC.th = 1,
#'  min.percent = 0.01)
#' 
#' 
enrich_TF_ChEA3_scExp = function(
    scExp, ref = "hg38", qval.th = 0.01, logFC.th = 1,
    min.percent = 0.01, peak_distance = 1000, use_peaks = FALSE,
    progress=NULL, verbose = TRUE){
  stopifnot(is(scExp, "SingleCellExperiment"), is.character(ref),
            is.numeric(peak_distance), is.numeric(qval.th),
            is.numeric(logFC.th), is.numeric(min.percent))
  if (!is.null(progress)) progress$inc(
    detail = "Loading reference gene sets...", amount = 0.1)
  
  if (!any(grepl("qval", colnames(SingleCellExperiment::rowData(scExp))))) 
    stop("ChromSCape::enrich_for_TF_ChEA3 - No DA, ",
         "please run run_differential_analysis_scExp first.")
  if (is.null(SingleCellExperiment::rowData(scExp)$Gene)) 
    stop("ChromSCape::enrich_for_TF_ChEA3 - No Gene ",
         "annotation, please annotate features with genes using ",
         "feature_annotation_scExp first.")
  if (is.null(use_peaks)) use_peaks = FALSE
  
  if (use_peaks & (!"refined_annotation" %in% names(scExp@metadata))) 
    stop("ChromSCape::enrich_for_TF_ChEA3 - ",
         "When use_peaks is TRUE, metadata must contain refined_annotation ",
         "object.")
  

  refined_annotation = NULL
  if(use_peaks) refined_annotation = scExp@metadata$refined_annotation
  
  if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
    stop("ChromSCape::enrich_for_TF_ChEA3 - scExp ",
         "object must have selected number of clusters.")

  annotFeat_long = as.data.frame(tidyr::separate_rows(
    as.data.frame(SummarizedExperiment::rowRanges(scExp)), 
    .data$Gene, sep = ", "))
  df = data.frame("Rank" = 0,
                  "TF" = 0,
                  "Score" = 0,
                  "Library" = NA,
                  "Overlapping_Genes" = 0,
                  "nTargets" = 0,
                  "totalTargetsTF" = 0,
                  'totalGenesInSet' = 0)
  res <- annotFeat_long
  groups <- gsub(".*\\.","", grep("qval",colnames(res), value = TRUE))
  TF_enrichment = list()
  for (i in seq_along(groups)) {
    enr = list(Differential = data.frame(), Enriched = data.frame(), Depleted = data.frame())
    gp = groups[i]
    if (!is.null(progress)) progress$inc(
      detail = paste0("TF enrichment for ",gp, "..."),
      amount =(0.7 / length(groups)))
    
    if(verbose) cat(paste0("TF enrichment for ",gp, "...\n"))
    qval.col <- paste("qval", gp, sep = ".")
    logFC.col <- paste("logFC", gp, sep = ".")
    group_activation.col <- paste("group_activation", gp, sep = ".")
    reference_activation.col <- paste("reference_activation", gp, sep = ".")
    
    over = res$ID[which(res[, qval.col] <= qval.th & 
                          res[, logFC.col] > logFC.th &
                          res[, group_activation.col] > min.percent)]
    overG = unique(
      annotFeat_long$Gene[annotFeat_long$distanceToTSS < peak_distance & 
                            annotFeat_long$ID %in% over])
    under = res$ID[which(res[, qval.col] <= qval.th &
                           res[,logFC.col] < -logFC.th  &
                           res[, reference_activation.col] > min.percent)]
    underG = unique(
      annotFeat_long$Gene[annotFeat_long$distanceToTSS < peak_distance & 
                            annotFeat_long$ID %in% under])
    signific = c(over, under)
    significG = unique(c(overG, underG))
    
    if (!is.null(refined_annotation)){
      out <- filter_genes_with_refined_peak_annotation(
        refined_annotation, peak_distance, signific, over, under)
      significG = out$significG
      overG = out$overG
      underG = out$underG
    }
    if (length(significG)) enr$Differential = enrich_TF_ChEA3_genes(significG) else
        enr$Differential = df
    if (length(overG)) enr$Enriched = enrich_TF_ChEA3_genes(overG) else
        enr$Differential = df
    if (length(underG)) enr$Depleted = enrich_TF_ChEA3_genes(underG) else
        enr$Differential = df
    TF_enrichment[[gp]] = enr
    
  }
  scExp@metadata$TF_enrichment= TF_enrichment
  return(scExp)
}


#' Find the TF that are enriched in the differential genes using ChEA3 API
#'
#' @param genes A character vector with the name of genes to enrich for TF.
#' 
#' @return Returns a SingleCellExperiment object containing list of enriched
#'   Gene Sets for each cluster, either in depleted features, enriched features
#'   or simply differential features (both).
#'   
#' @references Keenan AB, Torre D, Lachmann A, Leong AK, Wojciechowicz M, Utti V, 
#' Jagodnik K, Kropiwnicki E, Wang Z, Ma'ayan A (2019)
#'  ChEA3: transcription factor enrichment analysis by orthogonal omics integration. 
#'  Nucleic Acids Research. doi: 10.1093/nar/gkz446
#'  +
#' @export
#' @importFrom utils data
#' @examples
#' data(scExp)
#' enrich_TF_ChEA3_genes(head(unlist(strsplit(SummarizedExperiment::rowData(scExp)$Gene, split = ",", fixed = TRUE)), 15))
#' 
enrich_TF_ChEA3_genes = function(genes){
  response = NULL
  if(!requireNamespace("httr", quietly=TRUE)){
    warning("ChromSCape::enrich_TF_ChEA3_genes - In order to access ChEA3 database, you must install httr Package first.",
                            "Run install.packages('httr') in console. Exiting.")
    return()
  }
  if(!requireNamespace("jsonlite", quietly=TRUE)){
    warning("ChromSCape::enrich_TF_ChEA3_genes - In order to access ChEA3 database, you must install jsonlite Package first.",
            "Run install.packages('jsonlite') in console. Exiting.")
    return()
  }
  
    return_df = data.frame("Rank" = 0,
                          "TF" = 0,
                          "Score" = 0,
                          "Library" = NA,
                          "Overlapping_Genes" = 0,
                          "nTargets" = 0,
                          "totalTargetsTF" = 0,
                          'totalGenesInSet' = 0)
    
  if(length(genes) >= 10){
  utils::data("CheA3_TF_nTargets")
  
  #POST to ChEA3 server
  httr::set_config(httr::config(ssl_verifypeer = 0L))
  url = "https://maayanlab.cloud/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = "myQuery", gene_set = genes)
  
  #POST to ChEA3 server
  tryCatch({response = httr::POST(url = url, body = payload, encode = encode)},
           error = function(e) e)
  
  if(is.null(response)) return(return_df)
  json =  httr::content(response, "text")
  
  #results as list of R dataframes
  results =  tryCatch({  jsonlite::fromJSON(json)},
                      error = function(e) return_df) 
  if(length(results) > 1){
      results = results$`Integrated--meanRank`
      results$nTargets = sapply(results$Overlapping_Genes,
                                function(x) length(unlist(strsplit(x, split = ","))))
      results$Score = as.numeric(results$Score)
      results = results[,-1]
      results$totalTargetsTF = CheA3_TF_nTargets$nTargets_TF[
          match(results$TF,CheA3_TF_nTargets$TF)]
      results$totalGenesInSet = length(genes)
  }

  } else {
    results = return_df
    
  }
  return(results)
}
