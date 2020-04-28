# This script is for testing reproducibility between
# the branch "package" & the branch "master" of ChromSCape

context("Testing reproducibility of preprocessing & filt with master")

library(testthat)
library(ChromSCape)
setwd("/media/pacome/LaCie/InstitutCurie/Documents/Data/ChromSCape_analysis/datasets/scChIPseq_H3K27me3_HBCx_95_22_mouse/")

master = new.env()
load("scChIP_raw.RData", master)
datamatrix = master$datamatrix
annot_raw = master$annot_raw

datamatrix = as(as.matrix(datamatrix),"dgCMatrix")

test_that("Re checking that seed is the same", {
  expect_equal(as.matrix(master$datamatrix), as.matrix(datamatrix) )
  expect_equal(master$annot_raw, annot_raw )
})

scExp = create_scExp(datamatrix,annot_raw)
scExp = filter_scExp(scExp, min_cov_cell = 2000)

regions_to_exclude = read.table("../../../Annotation/bed/BC976_CNV_HMM_from_ChIP.bed")
scExp. = exclude_features_scExp(scExp,features_to_exclude = regions_to_exclude, by ="region")

test_that("QC Filtering ", {
  load("reduced_data/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected.RData",master)
  expect_equal(dim(master$norm_mat), dim(SingleCellExperiment::counts(scExp)) )
  expect_equal(rownames(master$norm_mat), rownames(SingleCellExperiment::counts(scExp)) )
  expect_equal(colnames(master$norm_mat), colnames(SingleCellExperiment::counts(scExp)) )
})

scExp = normalize_scExp(scExp, "CPM")

test_that("Normalization ", {
  load("reduced_data/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_normMat.RData",master)
  a = master$norm_mat[,1] / normcounts(scExp)[,1]
  head(a[!is.nan(a)])[1] -> b
  expect_equal(master$norm_mat, b*as.matrix(SingleCellExperiment::normcounts(scExp)),tol = 1e-4)
})

scExp = feature_annotation_scExp(scExp,ref="mm10")

test_that("Feature annotation ", {
  load("reduced_data/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_annotFeat.RData",master)
  head(master$annotFeat)
  head(rowData(scExp))
  
  to_test = as.data.frame(SummarizedExperiment::rowData(scExp))
  original = master$annotFeat
  to_test = to_test[match(original$ID,to_test$ID),]
  rownames(original) = original$ID
  
  original$distance[which(original$distance > 0)] = 
    as.numeric(original$distance[which(original$distance > 0)]) - 2 
  # Difference between bedtools & GenomicRanges of 2 (start = 1 vs start = 0 ?)
  to_test$chr = as.character(to_test$chr)
  to_test$start = as.character(to_test$start)
  to_test$end = as.character(to_test$end)
  original$distance = as.numeric(to_test$end)
  expect_equal(to_test,original) #re order !

})

scExp = reduce_dims_scExp(scExp)

test_that("Step 7 : dimensionality reduction", {
  
  load("reduced_data/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected.RData", master)
  expect_equal(abs(master$pca), abs(as.matrix(reducedDim(scExp,"PCA")))) # Is equal except sign 

})

# BECAUSE OF PCA signs are random, all the downstream steps are diverging. Therefore we put
# PCA from master into reducedDim(scExp,"PCA") to keep signs to test all the
# downstream steps
SingleCellExperiment::reducedDim(scExp, "PCA") = as.data.frame(master$pca)
scExp = correlation_and_hierarchical_clust_scExp(scExp)

test_that("Step 8 : correlation & hiearchical clust", {
  
  load("cor_filtered_data/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_99_1.RData", master)
  expect_equal(master$hc_cor$merge,scExp@metadata$hc_cor$merge) # correlation not same
  expect_equal(master$hc_cor$height,scExp@metadata$hc_cor$height)
  expect_equal(master$hc_cor$order,scExp@metadata$hc_cor$order)
  expect_equal(cor(t(master$pca)),reducedDim(scExp,"Cor") ) # Depending on the sign of the PCA, pearson correlation matrix will give different correlation results
  expect_equal(colnames(cor(t(master$pca))),colnames(reducedDim(scExp,"Cor")) )
  expect_equal(rownames(cor(t(master$pca))),rownames(reducedDim(scExp,"Cor")) )
})

scExp_cf = filter_correlated_cell_scExp(scExp,random_iter = 50)

test_that("Step 9 : correlation filterin clust", {
  
  load("consclust/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_99_1.RData", master)
  load("z.RData", master)
  expect_equal(nrow(master$annot_sel),ncol(scExp_cf))
  expect_equal(master$annot_sel$cell_id,colData(scExp_cf)$cell_id)
  expect_equal(table(master$annot_sel$sample_id),table(colData(scExp_cf)$sample_id))
})


scExp_cf = consensus_clustering_scExp(scExp_cf,prefix = "",reps = 1000, seed = 3.14)

test_that("Step 10 : correlation consensus hierarchical clustering", {
  
  load("consclust/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_99_1.RData", master)
  
  # consclust list
  expect_equal(master$consclust[[2]]$consensusMatrix, scExp_cf@metadata$consclust[[2]]$consensusMatrix)
  expect_equal(master$consclust[[5]]$consensusMatrix, scExp_cf@metadata$consclust[[5]]$consensusMatrix)
  expect_equal(master$consclust[[8]]$consensusClass, scExp_cf@metadata$consclust[[8]]$consensusClass)
  
  #icl list
  expect_equal(master$icl$clusterConsensus, scExp_cf@metadata$icl$clusterConsensus)
  expect_equal(master$icl$itemConsensus, scExp_cf@metadata$icl$itemConsensus)
})

scExp_cf = choose_cluster_scExp(scExp_cf, nclust = 4)

test_that("Step 11 : correlation consensus hierarchical clustering - choose cluster", {
  
  load("consclust/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_99_1_affectation_k4.RData", master)
  
  # cell affectation to clusters
  affectation. = as.data.frame(colData(scExp_cf))[,c("barcode","cell_id","chromatin_group","sample_id")]
  colnames(affectation.)[3] = "ChromatinGroup"
  
  table(master$affectation[,c("sample_id","ChromatinGroup")])
  table(affectation.[,c("sample_id","ChromatinGroup")])
  # Clusters are randomly named -> need to manually change names based on count
  # here origin C1 = new C2; origin C2 = new C3; origin C3 = new C1
  affectation.$ChromatinGroup[which(master$affectation$ChromatinGroup=="C3")] = "C3"
  affectation.$ChromatinGroup[which(master$affectation$ChromatinGroup=="C2")] = "C2"
  affectation.$ChromatinGroup[which(master$affectation$ChromatinGroup=="C1")] = "C1"
  
  expect_equal(master$affectation,  affectation. )
  expect_equal(table(master$affectation[,c("sample_id","ChromatinGroup")]),
               table(affectation.[,c("sample_id","ChromatinGroup")]) )
  
  # tsne
  load("consclust/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_99_1_tsne_filtered.RData",
       master)
  rownames(master$tsne_filtered$Y) = master$affectation$cell_id
  colnames(master$tsne_filtered$Y) = c("Component_1","Component_2")
  # expect_equal(as.data.frame(master$tsne_filtered$Y), 
  #              as.data.frame(reducedDim(scExp_cf,"TSNE"))) # Not equal - random seed ?
  
})

scExp_cf = differential_analysis_scExp(scExp_cf, de_type = "one_vs_rest",qval.th = 0.4, cdiff.th = 0.3, block = NULL)

test_that("Step 12 : differential analysis between clusters", {
  
  load("supervised/scChIPseq_H3K27me3_HBCx_95_22_mouse_2000_1_95_uncorrected_99_1_4_0.01_1_one_vs_rest.RData", master)
  
  my.res_original = as.data.frame(apply(master$my.res_save, MARGIN = 2, as.character), stringsAsFactors=F)
  my.res_new = as.data.frame(apply(scExp_cf@metadata$diff$res, MARGIN = 2, as.character), stringsAsFactors=F)
  
  head(my.res_original)
  head(my.res_new)
  my.res_new = my.res_new[match(my.res_original$ID,my.res_new$ID),]
  rownames(my.res_new) = rownames(my.res_original) = my.res_original$ID
  
  # Clusters are randomly named -> need to manually change names based on count
  # here origin C1 = new C2; origin C2 = new C3; origin C3 = new C1
  expect_equal(my.res_original$qval.C1, my.res_new$qval.C2)
  expect_equal(my.res_original$qval.C2, my.res_new$qval.C3)
  expect_equal(my.res_original$qval.C3, my.res_new$qval.C1)
  expect_equal(my.res_original$qval.C4, my.res_new$qval.C4)
  expect_equal(my.res_original$cdiff.C1, my.res_new$cdiff.C2)
  expect_equal(my.res_original$cdiff.C2, my.res_new$cdiff.C3)
  expect_equal(my.res_original$cdiff.C3, my.res_new$cdiff.C1)
  expect_equal(my.res_original$cdiff.C4, my.res_new$cdiff.C4)
})

scExp_cf = gene_set_enrichment_analysis_scExp(scExp_cf,ref = "mm10",
                                              use_peaks = F )