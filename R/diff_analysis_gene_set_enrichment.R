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
#' p-value threshold (qval.th) and the fold-change threshold (cdiff.th). It
#' outputs a SingleCellExperiment object containing a differential list.
#'
#' @param scExp A SingleCellExperiment object containing consclust with selected
#'   number of cluster.
#' @param de_type Type of comparisons. Either 'one_vs_rest', to compare each
#'   cluster against all others, or 'pairwise' to make 1 to 1 comparisons.
#'   ('one_vs_rest')
#' @param qval.th Adjusted p-value threshold. (0.01)
#' @param cdiff.th Fold change threshold. (1)
#' @param method Wilcoxon or edgerGLM
#' @param block Use batches as blocking factors ?
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
#'  data("scExp")
#' scExp_cf = differential_analysis_scExp(scExp_cf)
#' 
differential_analysis_scExp = function(
    scExp, de_type = "one_vs_rest", method = "wilcox", qval.th = 0.01, 
    cdiff.th = 1, block = NULL)
{
    warning_DA(scExp, de_type, method, qval.th, cdiff.th, block)
    
    if (isFALSE(block)) block = NULL
    nclust = length(unique(SingleCellExperiment::colData(scExp)$cell_cluster))
    if(method == "wilcox"){counts = SingleCellExperiment::normcounts(scExp)
    } else{counts = SingleCellExperiment::counts(scExp)}
    feature <- SingleCellExperiment::rowData(scExp)
    feature = data.frame(ID = feature[, "ID"], chr = feature[, "chr"],
        start = feature[, "start"], end = feature[, "end"])
    affectation = as.data.frame(SingleCellExperiment::colData(scExp))
    diff = list(res = NULL, summary = NULL, groups = NULL, refs = NULL)
    if (de_type == "one_vs_rest"){
        out <- DA_one_vs_rest_fun(affectation, nclust, counts,
                                method, feature, block)
    } else {
        out <- DA_pairwise(affectation,nclust, counts, method, feature, block)
    }
    res <- out$res
    groups <- out$groups
    refs <- out$refs
    diff$summary = matrix(nrow = 3, ncol = length(groups), dimnames = list(
        c("differential", "over", "under"), groups))
    for (k in seq_along(groups)){
        gpsamp = groups[k]
        qval.col <- paste("qval", gpsamp, sep = ".")
        cdiff.col <- paste("cdiff", gpsamp, sep = ".")
        diff$summary["differential", gpsamp] = sum(res[, qval.col] <= qval.th & 
                                                    abs(
                res[, cdiff.col]) > cdiff.th, na.rm = TRUE)
        diff$summary["over", gpsamp] = sum(res[,qval.col] <= qval.th & 
                res[, cdiff.col] > cdiff.th, na.rm = TRUE)
        diff$summary["under", gpsamp] = sum(res[, qval.col] <= qval.th &
                res[, cdiff.col] < -cdiff.th, na.rm = TRUE)
    }
    diff$groups = groups
    diff$refs = refs
    diff$res = res
    scExp@metadata$diff = diff
    return(scExp)
}

#' Warning for differential_analysis_scExp
#'
#' @param scExp A SingleCellExperiment object containing consclust with selected
#'   number of cluster.
#' @param de_type Type of comparisons. Either 'one_vs_rest', to compare each
#'   cluster against all others, or 'pairwise' to make 1 to 1 comparisons.
#'   ('one_vs_rest')
#' @param qval.th Adjusted p-value threshold. (0.01)
#' @param cdiff.th Fold change threshold. (1)
#' @param method Wilcoxon or edgerGLM
#' @param block Use batches as blocking factors ?
#'
#' @return Warnings or Errors if the input are not correct
warning_DA <- function(scExp, de_type, method, qval.th, cdiff.th, block){
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(de_type),
            is.numeric(qval.th), is.numeric(cdiff.th))
    if (!de_type %in% c("one_vs_rest", "pairwise")) 
        stop("ChromSCape::run_differential_analysis_scExp - de_type",
                    "must be 'one_vs_rest' or 'pairwise'.")
    if (!method %in% c("wilcox", "neg.binomial")) 
        stop("ChromSCape::run_differential_analysis_scExp - method must",
                    "be 'wilcox' or 'neg.binomial'.")
    if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::run_differential_analysis_scExp - scExp ",
                    "object must have selected number of clusters.")
    if (FALSE %in% (c("chr", "start", "end") %in% 
                    colnames(SingleCellExperiment::rowData(scExp)))) 
        stop("ChromSCape::run_differential_analysis_scExp - Please run ",
                    "feature_annotation_scExp first.")
}

#' Differential Analysis in 'One vs Rest' mode
#'
#' @param affectation An annotation data.frame with cell_id and cell_cluster
#'   columns
#' @param nclust Number of clusters
#' @param counts Count matrix
#' @param method DA method : Wilcoxon or EdgeR
#' @param feature Feature tables
#' @param block Blocking feature
#'
#' @return A list of results, groups compared and references
#'   
DA_one_vs_rest_fun <- function(affectation,nclust, counts, method, feature,
                            block){
    stopifnot(is.data.frame(affectation),is.integer(nclust),
                is(counts,"dgCMatrix")|is.matrix(counts), is.character(method),
                is.data.frame(feature))
    # compare each cluster to all the rest
    mygps = lapply(seq_len(nclust), function(i)
    {
        affectation[which(affectation$cell_cluster == paste0("C", i)),
                    "cell_id"]
    })
    names(mygps) = paste0("C", seq_len(nclust))
    groups = names(mygps)
    myrefs = lapply(seq_len(nclust), function(i)
    {
        affectation[which(affectation$cell_cluster != paste0("C", i)),
                    "cell_id"]
    })
    names(myrefs) = paste0("notC", seq_len(nclust))
    refs = names(myrefs)
    if(method == "wilcox"){ res = CompareWilcox(
        dataMat = counts, annot = affectation, ref_group = myrefs, 
        groups = mygps, featureTab = feature, block = block)
    } else {
        res = CompareedgeRGLM( dataMat = counts, 
                            annot = affectation, ref_group = myrefs,
                            groups = mygps, featureTab = feature)
        colnames(res)[
            grep("logCPM",colnames(res))] = gsub(
                "logCPM","Count",colnames(res)[grep("logCPM",colnames(res))])
        colnames(res)[grep("log2FC",colnames(res))] = gsub(
            "log2FC","cdiff",colnames(res)[grep("log2FC",colnames(res))])
        
    }
    out <- list("res" = res, "refs" = refs, "groups" = groups)
    return(out)
}

#' Run differential analysis in Pairwise mode
#'
#' @param affectation An annotation data.frame with cell_cluster and cell_id 
#' columns
#' @param nclust Number of clusters 
#' @param counts Count matrix
#' @param feature Feature data.frame
#' @param method DA method, Wilcoxon or edgeR
#' @param block Blocking feature
#'
#' @return A list of results, groups compared and references
#'
DA_pairwise <- function(affectation,nclust, counts,
                        method, feature, block){
    stopifnot(is.data.frame(affectation),is.integer(nclust),
                is(counts,"dgCMatrix")|is.matrix(counts), is.character(method),
                is.data.frame(feature))
    res = feature
    out <- run_pairwise_tests(affectation, nclust, counts,feature,method)
    count_save <- out$count_save
    pairs <- out$pairs
    single_results <- out$single_results
    
    tmp.gp = list(affectation[which(
        affectation$cell_cluster == paste0("C", nclust)), "cell_id"])
    count_save[, paste0("C", nclust)] = apply(
        counts, 1, function(x) mean(x[as.character(tmp.gp[[1]])]))
    combinedTests = scran::combineMarkers(
        de.lists = single_results, pairs = pairs, pval.field = "p.val",
        effect.field = "cdiff", pval.type = "any", log.p.in = FALSE, 
        log.p.out = FALSE, output.field = "stats", full.stats = TRUE)
    for (i in seq_len(as.integer(nclust)))
    {
        cdiffs = vapply(seq_len((as.integer(nclust) - 1)), function(k){
            combinedTests[[paste0("C", i)]][feature$ID, k + 4]$cdiff
        }, FUN.VALUE = rep(0,nrow(feature)))
        res[, paste0("Rank.C", i)] = combinedTests[[paste0("C", i)]][feature$ID,
                                                                    "Top"]
        res[, paste0("Count.C", i)] = as.numeric(count_save[, paste0("C", i)])
        res[, paste0("cdiff.C", i)] = Matrix::rowMeans(cdiffs)
        res[, paste0("pval.C", i)] = combinedTests[[paste0("C", i)]][feature$ID,
                                                                    "p.value"]
        res[, paste0("qval.C", i)] = combinedTests[[paste0("C", i)]][feature$ID,
                                                                    "FDR"]
    }
    groups = paste0("C", seq_len(nclust)) 
    refs = paste0("pairedTest", seq_len(nclust))
    out <- list( "res" = res, "refs" = refs, "groups" = groups)
    return(out)
}


#' Run pairwise tests
#'
#' @param affectation An annotation data.frame with cell_cluster and cell_id 
#' columins
#' @param nclust Number of clusters 
#' @param counts Count matrix
#' @param feature Feature data.frame
#' @param method DA method, Wilcoxon or edgeR
#'
#' @return A list containing objects for DA function
#'
run_pairwise_tests <- function(affectation, nclust, counts, 
                            feature, method){
    stopifnot(is.data.frame(affectation),is.integer(nclust),
                is(counts,"dgCMatrix")|is.matrix(counts), is.character(method))
    count_save = data.frame(ID = feature$ID)
    single_results = list()
    pairs = setNames(data.frame(matrix(ncol = 2, nrow = 0)),
                    c("group", "ref"))
    for (i in seq_len((as.integer(nclust) - 1))){
        mygps = list(affectation[which(
            affectation$cell_cluster == paste0("C", i)), "cell_id"])
        names(mygps) = paste0("C", i)
        groups = names(mygps)
        for (j in (i + 1):as.integer(nclust)){
            myrefs = list(affectation[which(
                affectation$cell_cluster == paste0("C",j)), "cell_id"])
            names(myrefs) = paste0("C", j)
            refs = names(myrefs)
            if(method == "wilcox"){ tmp_result = CompareWilcox(
                dataMat = counts, annot = affectation, ref_group = myrefs, 
                groups = mygps, featureTab = feature)
            } else {
                tmp_result = CompareedgeRGLM(
                    dataMat = counts, annot = affectation, ref_group = myrefs,
                    groups = mygps, featureTab = feature)
                logCPM.cols <- grep("logCPM",colnames(tmp_result))
                log2FC.cols <- grep("log2FC",colnames(tmp_result))
                colnames(tmp_result)[logCPM.cols] = gsub(
                    "logCPM","Count", colnames(tmp_result)[logCPM.cols])
                colnames(tmp_result)[log2FC.cols] = gsub(
                    "log2FC","cdiff",colnames(tmp_result)[log2FC.cols])}
            tmp_result = tmp_result[match(count_save$ID, tmp_result$ID), ]
            rownames(tmp_result) = tmp_result$ID
            if(method == "wilcox"){
                tmp_result[5] = NULL  
                colnames(tmp_result)[5:8] <- c("count", "cdiff",
                                            "p.val", "adj.p.val")} else {
                colnames(tmp_result)[5:8] = c("cdiff", "count",
                                            "p.val", "adj.p.val")}
            count_save[, paste0("C", i)] = tmp_result$count
            single_results = rlist::list.append(single_results, tmp_result)
            pairs[nrow(pairs) + 1, ] = list(groups[1], refs[1])
            tmp_mirror = tmp_result
            tmp_mirror$cdiff = tmp_mirror$cdiff * (-1)
            tmp_mirror$count = 0  
            single_results = rlist::list.append(single_results, tmp_mirror)
            pairs[nrow(pairs) + 1, ] = list(refs[1], groups[1])}}
    out <- list("single_results" = single_results, "pairs" = pairs,
                "count_save" = count_save)
    return(out)}

#' Runs Gene Set Enrichment Analysis on genes associated with differential
#' features
#'
#' This function takes previously calculated differential features and runs
#' hypergeometric test to look for enriched gene sets in the genes associated
#' with differential features, for each cell cluster. This functions takes as
#' input a SingleCellExperiment object with consclust, the type of comparison,
#' either 'one_vs_rest' or 'pairwise', the adjusted p-value threshold (qval.th)
#' and the fold-change threshold (cdiff.th). It outputs a SingleCellExperiment
#' object containing a differential list.
#'
#' @param scExp A SingleCellExperiment object containing list of differential
#'   features.
#' @param ref A reference annotation. ('hg38')
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
#' @param cdiff.th Fold change threshold to define differential features. (1)
#' @param peak_distance Maximum distanceToTSS of feature to gene TSS to consider
#'   associated, in bp. (1000)
#' @param use_peaks Use peak calling method (must be calculated beforehand).
#'   (FALSE)
#' @param GeneSetClasses Which classes of MSIGdb to look for.  
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
#' #Usually recommanding qval.th = 0.01 & cdiff.th = 1 or 2
#' \dontrun{scExp_cf = gene_set_enrichment_analysis_scExp(scExp_cf,
#'  qval.th = 0.4, cdiff.th = 0.3)}
#' 
gene_set_enrichment_analysis_scExp = function(
    scExp, enrichment_qval = 0.1, ref = "hg38", GeneSets = NULL,
    GeneSetsDf = NULL, GenePool = NULL, qval.th = 0.01, cdiff.th = 1, 
    peak_distance = 1000, use_peaks = FALSE, GeneSetClasses = c(
        "c1_positional","c2_curated","c3_motif","c4_computational",
        "c5_GO","c6_oncogenic","c7_immunologic","hallmark"))
{
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(ref),
            is.numeric(enrichment_qval), is.numeric(peak_distance),
            is.numeric(qval.th), is.numeric(cdiff.th), 
            is.character(GeneSetClasses))
    if (is.null(scExp@metadata$diff)) 
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
    
    if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - scExp ",
        "object must have selected number of clusters.")
    
    database_MSIG <- load_MSIGdb(ref,GeneSetClasses)
    GeneSets = database_MSIG$GeneSets
    GeneSetsDf = database_MSIG$GeneSetsDf
    GenePool = database_MSIG$GenePool
    
    nclust = length(unique(SingleCellExperiment::colData(scExp)$cell_cluster))
    
    annotFeat_long = as.data.frame(tidyr::separate_rows(
        as.data.frame(SingleCellExperiment::rowData(scExp)), 
        .data$Gene, sep = ", "))
    enr <- combine_enrichmentTests(
        scExp@metadata$diff, enrichment_qval, qval.th, cdiff.th,
        annotFeat_long,peak_distance, refined_annotation, GeneSets,
        GeneSetsDf, GenePool)
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
    GenePool = unique(as.character(GenePool[, "gene"]))
    
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
    
    enrich.test = merge(subset(GeneSetsDf, select = -Genes), enrich.test, 
                        by.x = "Gene.Set", by.y = "Gene_set_name",
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
#' @param cdiff.th Differential analysis log-fold change threshold
#' @param annotFeat_long Long annotation
#' @param peak_distance Maximum gene to peak distance
#' @param refined_annotation Refined annotation data.frame if peak calling is 
#' done
#' @param GeneSets List of pathways
#' @param GeneSetsDf Data.frame of pathways
#' @param GenePool Pool of possible genes for testing
#'
#' @return A list of list of pathway enrichment data.frames for 
#' Both / Over / Under and for each cluster
#'
combine_enrichmentTests <- function(
    diff, enrichment_qval, qval.th, cdiff.th, annotFeat_long, peak_distance,
    refined_annotation, GeneSets, GeneSetsDf, GenePool)
    {
    stopifnot(is.list(diff), is.numeric(enrichment_qval), is.numeric(qval.th),
              is.numeric(cdiff.th), is.data.frame(annotFeat_long),
              is.numeric(peak_distance), is(refined_annotation,"GRanges"),
              is.list(GeneSets),is.data.frame(GeneSetsDf),
              is.character(GenePool))
    groups <- diff$groups
    res <- diff$res
    enr = list(Both = list(), Overexpressed = list(), Underexpressed = list())
    for (i in seq_along(groups)) {
        gp = groups[i]
        qval.col <- paste("qval", gp, sep = ".")
        cdiff.col <- paste("cdiff", gp, sep = ".")
        signific = res$ID[which(res[, qval.col] <= qval.th & 
                                    abs(res[, cdiff.col]) > cdiff.th)]
        significG = unique(
            annotFeat_long$Gene[annotFeat_long$distanceToTSS < peak_distance & 
                                    annotFeat_long$ID %in% signific])
        over = res$ID[which(res[, qval.col] <= qval.th & 
                                res[, cdiff.col] > cdiff.th)]
        overG = unique(
            annotFeat_long$Gene[annotFeat_long$distanceToTSS < peak_distance & 
                                            annotFeat_long$ID %in% over])
        under = res$ID[which(res[, qval.col] <= qval.th &
                                res[,cdiff.col] < -cdiff.th)]
        underG = unique(
            annotFeat_long$Gene[annotFeat_long$distanceToTSS < peak_distance & 
                                    annotFeat_long$ID %in% under])
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
#' @param cell_cluster Cell cluster. ('C1')
#' @param enr_class_sel Which classes of gene sets to show. (c('c1_positional',
#'   'c2_curated', ...))
#'
#' @return A DT::data.table of enriched gene sets.
#' @export
#'
#' @importFrom DT datatable
#' @importFrom tidyr unite
#' 
#' @examples
#' data("scExp")
#' table_enriched_genes_scExp(scExp_cf)
table_enriched_genes_scExp <- function(
    scExp, set = "Both", cell_cluster = "C1", 
    enr_class_sel = c(
        "c1_positional", "c2_curated", "c3_motif", "c4_computational", "c5_GO",
        "c6_oncogenic", "c7_immunologic", "hallmark"))
{
    stopifnot(is(scExp, "SingleCellExperiment"),
            is.character(set), is.character(cell_cluster), 
            is.character(enr_class_sel))
    if (is.null(scExp@metadata$enr)) 
        stop("ChromSCape::table_enriched_genes_scExp - No GSEA, please ",
        "run gene_set_enrichment_analysis_scExp first.")
    
    if (!set %in% c("Both", "Overexpressed", "Underexpressed")) 
        stop("ChromSCape::table_enriched_genes_scExp - set variable",
        "must be 'Both', 'Overexpressed' or 'Underexpressed'.")
    
    if (!cell_cluster %in% scExp@metadata$diff$groups) 
        stop("ChromSCape::table_enriched_genes_scExp - No GSEA, please",
        "run gene_set_enrichment_analysis_scExp first.")
    
    table <- scExp@metadata$enr[[set]][[match(cell_cluster,
                                            scExp@metadata$diff$groups)]]
    table <- table[which(table[, "Class"] %in% enr_class_sel), ]
    if (is.null(table))
    {
        return(
            setNames(
                data.frame(matrix(ncol = 6, nrow = 0)),
                c("Gene_set", "Class", "Num_deregulated_genes", "p.value",
                "q.value", "Deregulated_genes")))
    }
    table <- tidyr::unite(table, "dereg_genes", c(
        "Nb_of_deregulated_genes", "Nb_of_genes"), sep = "/")
    colnames(table) <- c("Gene_set", "Class", "Num_deregulated_genes", 
                        "p.value", "adj.p.value", "Deregulated_genes")
    table[, 4] <- round(table[, 4], 9)
    table[, 5] <- round(table[, 5], 9)
    table <- table[order(table$adj.p.value, table$p.value), ]
    DT::datatable(table, options = list(dom = "tpi"))
}
