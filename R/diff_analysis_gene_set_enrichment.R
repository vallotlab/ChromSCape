## Authors : Pacôme Prompsy Title : Wrappers & functions to run differential
## analysis on clusters of single-cell based on epigenomic profiling and Gene Set
## Enrichment Analysis on the genes associated with differentially enriched
## features

#' Runs differential analysis between cell clusters
#'
#' Based on clusters of cell defined previously, runs non-parametric Wilcoxon Rank Sum test to
#' find significantly depleted or enriched features, in 'one_vs_rest' mode or 'pairwise' mode.
#' In pairwise mode, each cluster is compared to all other cluster individually, and then pairwise
#' comparisons between clusters are combined to find overall differential features using combineMarkers
#'  function from scran.
#'
#' This functions takes as input a SingleCellExperiment object with consclust, the type of comparison, either 
#' 'one_vs_rest' or 'pairwise', the adjusted p-value threshold (qval.th) and the fold-change threshold (cdiff.th). 
#' It outputs a SingleCellExperiment object containing a differential list.
#' 
#' @param scExp A SingleCellExperiment object containing consclust with selected number of cluster.
#' @param de_type Type of comparisons. Either 'one_vs_rest', to compare each cluster against all others, or 
#' 'pairwise' to make 1 to 1 comparisons. ['one_vs_rest']
#' @param qval.th Adjusted p-value threshold. [0.01]
#' @param cdiff.th Fold change threshold. [1]
#' @param block 
#'
#' @return Returns a SingleCellExperiment object containing a differential list.
#' @export
#'
#' @examples
#' @importFrom Matrix rowMeans
#' @importFrom scran combineMarkers 
#' @importFrom SingleCellExperiment colData normcounts rowData
#' @importFrom rlist list.append
differential_analysis_scExp = function(scExp, de_type = "one_vs_rest", 
                                       method = "wilcox",
                                       qval.th = 0.01, 
                                       cdiff.th = 1, block = NULL)
    {
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(de_type), is.numeric(qval.th), 
        is.numeric(cdiff.th))
    
    if (!de_type %in% c("one_vs_rest", "pairwise")) 
        stop("ChromSCape::run_differential_analysis_scExp - de_type must be 'one_vs_rest' or 'pairwise'.")
    if (!method %in% c("wilcox", "neg.binomial")) 
        stop("ChromSCape::run_differential_analysis_scExp - method must be 'wilcox' or 'neg.binomial'.")
    
    if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::run_differential_analysis_scExp - scExp object must have selected number of clusters.")
    
    if (FALSE %in% (c("chr", "start", "end") %in% colnames(SingleCellExperiment::rowData(scExp)))) 
        stop("ChromSCape::run_differential_analysis_scExp - Please run feature_annotation_scExp first.")
    
    if (isFALSE(block)) 
        block = NULL
    nclust = length(unique(SingleCellExperiment::colData(scExp)$cell_cluster))
    
    if(method == "wilcox"){
        counts = SingleCellExperiment::normcounts(scExp)
    } else {
        counts = SingleCellExperiment::counts(scExp)
    }
    
    feature = data.frame(ID = SingleCellExperiment::rowData(scExp)[, "ID"], chr = SingleCellExperiment::rowData(scExp)[, 
        "chr"], start = SingleCellExperiment::rowData(scExp)[, "start"], end = SingleCellExperiment::rowData(scExp)[, 
        "end"])
    affectation = as.data.frame(SingleCellExperiment::colData(scExp))
    
    diff = list(res = NULL, summary = NULL, groups = NULL, refs = NULL)
    
    if (de_type == "one_vs_rest")
    {
        # compare each cluster to all the rest
        mygps = lapply(1:nclust, function(i)
        {
            affectation[which(affectation$cell_cluster == paste0("C", i)), "cell_id"]
        })
        names(mygps) = paste0("C", 1:nclust)
        groups = names(mygps)
        myrefs = lapply(1:nclust, function(i)
        {
            affectation[which(affectation$cell_cluster != paste0("C", i)), "cell_id"]
        })
        names(myrefs) = paste0("notC", 1:nclust)
        refs = names(myrefs)
        if(method == "wilcox"){ res = CompareWilcox(
            dataMat = counts, annot = affectation, ref = myrefs, 
            groups = mygps, featureTab = feature, block = block)
        } else {
            res = CompareedgeRGLM( dataMat = counts, 
                                        annot = affectation, ref = myrefs,
                                        groups = mygps, featureTab = feature)
            colnames(res)[grep("logCPM",colnames(res))] = gsub("logCPM","Count",
                                                               colnames(res)[grep("logCPM",colnames(res))])
            colnames(res)[grep("log2FC",colnames(res))] = gsub("log2FC","cdiff",
                                                               colnames(res)[grep("log2FC",colnames(res))])
        }
    } else
    {
        # pairwise one-vs-one testing for each cluster
        
        res = feature
        count_save = data.frame(ID = feature$ID)
        single_results = list()
        pairs = setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("group", "ref"))
        for (i in 1:(as.integer(nclust) - 1))
        {
            mygps = list(affectation[which(affectation$cell_cluster == paste0("C", 
                i)), "cell_id"])
            names(mygps) = paste0("C", i)
            groups = names(mygps)
            for (j in (i + 1):as.integer(nclust))
            {
                myrefs = list(affectation[which(affectation$cell_cluster == paste0("C", 
                  j)), "cell_id"])
                names(myrefs) = paste0("C", j)
                refs = names(myrefs)
                if(method == "wilcox") tmp_result = CompareWilcox(
                    dataMat = counts, annot = affectation, ref = myrefs, 
                    groups = mygps, featureTab = feature)
                else {
                    tmp_result = CompareedgeRGLM(
                    dataMat = counts, annot = affectation, ref = myrefs,
                    groups = mygps, featureTab = feature)
                colnames(tmp_result)[grep("logCPM",colnames(tmp_result))] = gsub("logCPM","Count",
                                                                   colnames(tmp_result)[grep("logCPM",colnames(tmp_result))])
                colnames(tmp_result)[grep("log2FC",colnames(tmp_result))] = gsub("log2FC","cdiff",
                                                                   colnames(tmp_result)[grep("log2FC",colnames(tmp_result))])
                }
                
                tmp_result = tmp_result[match(count_save$ID, tmp_result$ID), ]
                rownames(tmp_result) = tmp_result$ID
                if(method == "wilcox"){
                    tmp_result[5] = NULL  #remove rank because it will be added later
                    colnames(tmp_result)[5:8] = c("count", "cdiff", "p.val", "adj.p.val")
                } else {
                    colnames(tmp_result)[5:8] = c("cdiff", "count",  "p.val", "adj.p.val")
                }
                count_save[, paste0("C", i)] = tmp_result$count
                single_results = list.append(single_results, tmp_result)
                pairs[nrow(pairs) + 1, ] = list(groups[1], refs[1])
                tmp_mirror = tmp_result
                tmp_mirror$cdiff = tmp_mirror$cdiff * (-1)
                tmp_mirror$count = 0  #not correct, but doesn't matter because it won't be used
                single_results = list.append(single_results, tmp_mirror)
                pairs[nrow(pairs) + 1, ] = list(refs[1], groups[1])
            }
        }
        # get count for last group as the loop doesn't cover it
        tmp.gp = list(affectation[which(affectation$cell_cluster == paste0("C", 
            nclust)), "cell_id"])
        count_save[, paste0("C", nclust)] = apply(counts, 1, function(x) mean(x[as.character(tmp.gp[[1]])]))
        combinedTests = scran::combineMarkers(de.lists = single_results, pairs = pairs, 
            pval.field = "p.val", effect.field = "cdiff", pval.type = "any", log.p.in = FALSE, 
            log.p.out = FALSE, output.field = "stats", full.stats = TRUE)
        for (i in 1:as.integer(nclust))
        {
            cdiffs = sapply(1:(as.integer(nclust) - 1), function(k)
            {
                combinedTests[[paste0("C", i)]][feature$ID, k + 3]$cdiff
            })
            res[, paste0("Rank.C", i)] = combinedTests[[paste0("C", i)]][feature$ID, 
                "Top"]
            res[, paste0("Count.C", i)] = as.numeric(count_save[, paste0("C", i)])
            res[, paste0("cdiff.C", i)] = Matrix::rowMeans(cdiffs)
            res[, paste0("pval.C", i)] = combinedTests[[paste0("C", i)]][feature$ID, 
                "p.value"]
            res[, paste0("qval.C", i)] = combinedTests[[paste0("C", i)]][feature$ID, 
                "FDR"]
        }
        groups = paste0("C", 1:nclust)  #needed for following code
        refs = paste0("pairedTest", 1:nclust)
    }
    
    diff$summary = matrix(nrow = 3, ncol = length(groups), dimnames = list(c("differential", 
        "over", "under"), groups))
    for (k in 1:length(groups))
    {
        gpsamp = groups[k]
        
        # For log2(x1/x2) > 1 || log2(x1/x2) > -1
        diff$summary["differential", gpsamp] = sum(res[, paste("qval", gpsamp, sep = ".")] <= 
            qval.th & abs(res[, paste("cdiff", gpsamp, sep = ".")]) > cdiff.th, na.rm = T)
        diff$summary["over", gpsamp] = sum(res[, paste("qval", gpsamp, sep = ".")] <= 
            qval.th & res[, paste("cdiff", gpsamp, sep = ".")] > cdiff.th, na.rm = T)
        diff$summary["under", gpsamp] = sum(res[, paste("qval", gpsamp, sep = ".")] <= 
            qval.th & res[, paste("cdiff", gpsamp, sep = ".")] < -cdiff.th, na.rm = T)
        
    }
    
    diff$groups = groups
    diff$refs = refs
    diff$res = res
    
    scExp@metadata$diff = diff
    return(scExp)
}


#' Runs Gene Set Enrichment Analysis on genes associated with differential features
#'
#' This function takes previously calculated differential features and runs hypergeometric test
#' to look for enriched gene sets in the genes associated with differential features, for each cell 
#' cluster.
#' This functions takes as input a SingleCellExperiment object with consclust, the type of comparison, either 
#' 'one_vs_rest' or 'pairwise', the adjusted p-value threshold (qval.th) and the fold-change threshold (cdiff.th). 
#' It outputs a SingleCellExperiment object containing a differential list.
#' @param scExp A SingleCellExperiment object containing list of differential features.
#' @param ref A reference annotation. ['hg38']
#' @param enrichment_qval Adjusted p-value threshold for gene set enrichment. [0.1]
#' @param GeneSets A named list of gene sets. If NULL will automatically load MSigDB list
#' of gene sets for specified reference genome. [NULL]
#' @param GeneSetsDf A dataframe containing gene sets & class of gene sets. If NULL will automatically 
#' load MSigDB dataframe of gene sets for specified reference genome. [NULL]
#' @param GenePool The pool of genes to run enrichment in. If NULL will automatically load 
#' Gencode list of genes fro specified reference genome. [NULL] 
#' @param qval.th Adjusted p-value threshold to define differential features. [0.01]
#' @param cdiff.th Fold change threshold to define differential features. [1]
#' @param peak_distance Maximum distance of feature to gene TSS to consider associated, in bp. [1000]
#' @param use_peaks Use peak calling method (must be calculated beforehand). [F]
#'
#' @return Returns a SingleCellExperiment object containing list of enriched Gene Sets for each cluster, either
#' in depleted features, enriched features or simply differential features (both). 
#' 
#' @export
#'
#' @examples
#' 
#' @importFrom SingleCellExperiment colData normcounts rowData
gene_set_enrichment_analysis_scExp = function(scExp, enrichment_qval = 0.1, ref = "hg38", 
    GeneSets = NULL, GeneSetsDf = NULL, GenePool = NULL, qval.th = 0.01, cdiff.th = 1, 
    peak_distance = 1000, use_peaks = F)
    {
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(ref), is.numeric(enrichment_qval), 
        is.numeric(peak_distance), is.numeric(qval.th), is.numeric(cdiff.th))
    
    if (is.null(scExp@metadata$diff)) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - No DA, please run run_differential_analysis_scExp first.")
    
    if (is.null(SingleCellExperiment::rowData(scExp)$Gene)) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - No Gene Annotation, please annotate features with genes 
         using feature_annotation_scExp first.")
    
    if (is.null(use_peaks)) 
        use_peaks = F
    
    if (use_peaks & (!"refined_annotation" %in% names(scExp@metadata))) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - When use_peaks is TRUE, metadata must contain refined_annotation object.")
    
    if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - scExp object must have selected number of clusters.")
    
    if ((!ref %in% c("hg38", "mm10")) & (is.null(GeneSets) | is.null(GeneSetsDf) | 
        is.null(GenePool))) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - Reference genome (ref) must be 'hg38' or 'mm10' if gene sets 
         not specified.")
    
    if (is.null(GeneSets) | is.null(GeneSetsDf))
    {
        message(paste0("ChromSCape::gene_set_enrichment_analysis_scExp - Selecting ", 
            ref, " MSigDB gene sets."))
        eval(parse(text = paste0("data(", ref, ".MSIG.ls)")))
        eval(parse(text = paste0("data(", ref, ".MSIG.gs)")))
        eval(parse(text = paste0("GeneSets = ", ref, ".MSIG.ls")))
        eval(parse(text = paste0("GeneSetsDf = ", ref, ".MSIG.gs")))
        
    }
    
    if (is.null(GenePool))
    {
        message(paste0("ChromSCape::gene_set_enrichment_analysis_scExp - Selecting ", 
            ref, " genes from Gencode."))
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        GenePool = eval(parse(text = paste0("", ref, ".GeneTSS")))
        GenePool = unique(as.character(GenePool[, "gene"]))
    }
    
    nclust = length(unique(SingleCellExperiment::colData(scExp)$cell_cluster))
    
    enr = list(Both = NULL, Overexpressed = NULL, Underexpressed = NULL)
    
    Both = list()
    Overexpressed = list()
    Underexpressed = list()
    
    groups = scExp@metadata$diff$groups
    res = scExp@metadata$diff$res
    diff = scExp@metadata$diff
    
    annotFeat_long = as.data.frame(tidyr::separate_rows(as.data.frame(SingleCellExperiment::rowData(scExp)), 
        Gene, sep = ", "))
    
    for (i in 1:length(groups))
    {
        gp = groups[i]
        ref = diff$refs[i]
        signific = res$ID[which(res[, paste("qval", gp, sep = ".")] <= qval.th & 
            abs(res[, paste("cdiff", gp, sep = ".")]) > cdiff.th)]
        significG = unique(annotFeat_long$Gene[annotFeat_long$distance < peak_distance & 
            annotFeat_long$ID %in% signific])
        over = res$ID[which(res[, paste("qval", gp, sep = ".")] <= qval.th & res[, 
            paste("cdiff", gp, sep = ".")] > cdiff.th)]
        overG = unique(annotFeat_long$Gene[annotFeat_long$distance < peak_distance & 
            annotFeat_long$ID %in% over])
        under = res$ID[which(res[, paste("qval", gp, sep = ".")] <= qval.th & res[, 
            paste("cdiff", gp, sep = ".")] < -cdiff.th)]
        underG = unique(annotFeat_long$Gene[annotFeat_long$distance < peak_distance & 
            annotFeat_long$ID %in% under])
        if (!is.null(use_peaks))
        {
            if (use_peaks == TRUE)
            {
                refined_annotation = scExp@metadata$refined_annotation
                signific_associated_peak = refined_annotation$peak_ID[refined_annotation$window_ID %in% 
                  signific]
                over_associated_peak = refined_annotation$peak_ID[refined_annotation$window_ID %in% 
                  over]
                under_associated_peak = refined_annotation$peak_ID[refined_annotation$window_ID %in% 
                  under]
                signific_associated_gene = refined_annotation$Gene[refined_annotation$peak_ID %in% 
                  signific_associated_peak & refined_annotation$distance < peak_distance]
                over_associated_gene = refined_annotation$Gene[refined_annotation$peak_ID %in% 
                  over_associated_peak & refined_annotation$distance < peak_distance]
                under_associated_gene = refined_annotation$Gene[refined_annotation$peak_ID %in% 
                  under_associated_peak & refined_annotation$distance < peak_distance]
                significG = unique(signific_associated_gene)
                overG = unique(over_associated_gene)
                underG = unique(under_associated_gene)
            }
        }
        
        if (length(significG))
        {
            enrich.test = enrichmentTest(gene.sets = GeneSets, mylist = significG, 
                possibleIds = GenePool)
            enrich.test = data.frame(Gene_set_name = rownames(enrich.test), enrich.test, 
                check.names = FALSE)
            enrich.test = merge(subset(GeneSetsDf, select = -Genes), enrich.test, 
                by.x = "Gene.Set", by.y = "Gene_set_name", all.y = TRUE, sort = FALSE)  ## Get class of gene set
            enrich.test = enrich.test[order(enrich.test$`p-value`), ]
            ind = which(enrich.test$`q-value` <= enrichment_qval)
            if (!length(ind))
            {
                ind = 1:10
            }
            Both[[i]] = enrich.test[ind, ]
        }
        if (length(overG))
        {
            enrich.test = enrichmentTest(gene.sets = GeneSets, mylist = overG, 
                possibleIds = GenePool)
            enrich.test = data.frame(Gene_set_name = rownames(enrich.test), enrich.test, 
                check.names = FALSE)
            enrich.test = merge(subset(GeneSetsDf, select = -Genes), enrich.test, 
                by.x = "Gene.Set", by.y = "Gene_set_name", all.y = TRUE, sort = FALSE)  ## Get class of gene set
            enrich.test = enrich.test[order(enrich.test$`p-value`), ]
            ind = which(enrich.test$`q-value` <= enrichment_qval)
            if (!length(ind))
            {
                ind = 1:10
            }
            Overexpressed[[i]] = enrich.test[ind, ]
        }
        if (length(underG))
        {
            enrich.test = enrichmentTest(gene.sets = GeneSets, mylist = underG, 
                possibleIds = GenePool)
            enrich.test = data.frame(Gene_set_name = rownames(enrich.test), enrich.test, 
                check.names = FALSE)
            enrich.test = merge(subset(GeneSetsDf, select = -Genes), enrich.test, 
                by.x = "Gene.Set", by.y = "Gene_set_name", all.y = TRUE, sort = FALSE)  ## Get class of gene set
            enrich.test = enrich.test[order(enrich.test$`p-value`), ]
            ind = which(enrich.test$`q-value` <= enrichment_qval)
            if (!length(ind))
            {
                ind = 1:10
            }
            Underexpressed[[i]] = enrich.test[ind, ]
        }
    }
    enr$Both = Both
    enr$Overexpressed = Overexpressed
    enr$Underexpressed = Underexpressed
    scExp@metadata$enr = enr
    return(scExp)
}

#' Creates table of enriched genes sets
#'
#' @param scExp  A SingleCellExperiment object containing list of enriched gene sets.
#' @param set A character vector, either 'Both', 'Overexpressed' or 'Underexpressed'. ['Both']
#' @param cell_cluster Cell cluster. ['C1']
#' @param enr_class_sel Which classes of gene sets to show. [c('c1_positional', 'c2_curated', ...)]
#'
#' @return A DT::data.table of enriched gene sets.
#' @export
#'
#' @importFrom DT datatable
#' @importFrom tidyr unite
table_enriched_genes_scExp <- function(scExp, set = "Both", cell_cluster = "C1", 
    enr_class_sel = c("c1_positional", "c2_curated", "c3_motif", "c4_computational", 
        "c5_GO", "c6_oncogenic", "c7_immunologic", "hallmark"))
        {
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(set), is.character(cell_cluster), 
        is.character(enr_class_sel))
    
    if (is.null(scExp@metadata$enr)) 
        stop("ChromSCape::table_enriched_genes_scExp - No GSEA, please run gene_set_enrichment_analysis_scExp first.")
    
    if (!set %in% c("Both", "Overexpressed", "Underexpressed")) 
        stop("ChromSCape::table_enriched_genes_scExp - set variable must be 'Both', 'Overexpressed' or 'Underexpressed'.")
    
    if (!cell_cluster %in% scExp@metadata$diff$groups) 
        stop("ChromSCape::table_enriched_genes_scExp - No GSEA, please run gene_set_enrichment_analysis_scExp first.")
    
    
    table <- scExp@metadata$enr[[set]][[match(cell_cluster, scExp@metadata$diff$groups)]]
    table <- table[which(table[, "Class"] %in% enr_class_sel), ]
    if (is.null(table))
    {
        return(setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene_set", "Class", 
            "Num_deregulated_genes", "p.value", "q.value", "Deregulated_genes")))
    }
    table <- tidyr::unite(table, "dereg_genes", c("Nb_of_deregulated_genes", "Nb_of_genes"), 
        sep = "/")
    colnames(table) <- c("Gene_set", "Class", "Num_deregulated_genes", "p.value", 
        "adj.p.value", "Deregulated_genes")
    table[, 4] <- round(table[, 4], 9)
    table[, 5] <- round(table[, 5], 9)
    table <- table[order(table$adj.p.value, table$p.value), ]
    DT::datatable(table, options = list(dom = "tpi"))
}
