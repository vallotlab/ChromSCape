# This script is for testing reproducibility between
# the branch "package" & the branch "master" of ChromSCape

context("Testing reproducibility of preprocessing & filt with master")

library(testthat)
library(ChromSCape)
setwd("/media/pacome/LaCie/InstitutCurie/Documents/Data/ChromSCape_Data/Reproducibility_new/Rdata_for_tests/Synthetic_dataset/")
out = create_scDataset_raw(featureType = "window", sparse = T)
datamatrix = out$mat
annot_raw = out$annot

master = new.env()
load("scChIP_raw.RData", master)
test_that("Re checking that seed is the same", {
  
  expect_equal(master$datamatrix, as.matrix(datamatrix) )
  expect_equal(master$annot_raw, annot_raw )
})

scExp = create_scExp(datamatrix,annot_raw)

test_that("Step 1 : creating scExp", {
  load("filter_red_01.RData", master)
  expect_equal(SummarizedExperiment::assay(master$umi), as.matrix(SummarizedExperiment::assay(scExp)) )
  expect_equal(colData(master$umi), colData(scExp))
})

scExp = filter_scExp(scExp)

test_that("Step 2 : filtering ", {
  load("filter_red_02.RData", master)
  expect_equal(master$SelMatCov, as.matrix(SummarizedExperiment::assay(scExp)))
})

regions_to_exclude = read.table("../../Data/Annotation/bed/MM468_5FU3_all_5FU5_initial_CNV_from_ChIP_input.bed")
scExp = exclude_features_scExp(scExp,features_to_exclude = regions_to_exclude, by ="region")

test_that("Step 3 : remove specific features ", {
  load("filter_red_03.RData", master)
  expect_equal(master$SelMatCov, as.matrix(SummarizedExperiment::assay(scExp)))
})

scExp = normalize_scExp(scExp, "CPM")
test_that("Step 4 : Normalize ", {
  load("filter_red_04.RData", master)
  expect_equivalent(master$norm_mat, as.matrix(SingleCellExperiment::normcounts(scExp)) )
})

scExp = feature_annotation_scExp(scExp,ref="hg38")

test_that("Step 5 : feature annotation ", {
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_annotFeat.RData",master)
  to_test = as.numeric(as.data.frame(SummarizedExperiment::rowData(scExp))$distance)
  original = as.numeric(GenomicRanges::makeGRangesFromDataFrame(master$annotFeat,
                                                 keep.extra.columns = T)$distance)
  original[original > 0] = original[original > 0] -2 # Difference between bedtools & GenomicRanges of 2 (start = 1 vs start = 0 ?)
  expect_equal(to_test,original[match(rownames(scExp),master$annotFeat$ID)]) #re order !
  
  expect_equal(master$annotFeat$Gene[match(rownames(scExp),master$annotFeat$ID)],rowData(scExp)$Gene)
  
})

# test_that("Step 6 : batch correction ", {
#   scExp = ChromSCape::filter_scExp(scExp)
#   load("filter_red_02.RData", master)
#   expect_equal(master$SelMatCov, assay(scExp))
# })

scExp = reduce_dims_scExp(scExp)

test_that("Step 7 : dimensionality reduction", {
  
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected.RData", master)
  # expect_equal(as.data.frame(master$pca), reducedDim(scExp,"PCA")) # at a sign
  expect_equal(cor(master$pca[,1], reducedDim(scExp,"PCA")[,1]), -1)
  expect_equal(cor(master$pca[,2], reducedDim(scExp,"PCA")[,2]), 1)
  expect_equal(cor(master$pca[,3], reducedDim(scExp,"PCA")[,3]), 1)

})

# BECAUSE OF PCA signs are random, all the downstream steps are diverging. Therefore we put
# PCA from master into reducedDim(scExp,"PCA") to keep signs to test all the
# downstream steps
SingleCellExperiment::reducedDim(scExp, "PCA") = as.data.frame(master$pca)
scExp = correlation_and_hierarchical_clust_scExp(scExp)

test_that("Step 8 : correlation & hiearchical clust", {
  
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_cor_filtered.RData", master)
  expect_equal(master$hc_cor$merge,scExp@metadata$hc_cor$merge) # correlation not same
  expect_equal(master$hc_cor$height,scExp@metadata$hc_cor$height)
  expect_equal(master$hc_cor$order,scExp@metadata$hc_cor$order)
  expect_equal(cor(t(master$pca)),reducedDim(scExp,"Cor") ) # Depending on the sign of the PCA, pearson correlation matrix will give different correlation results
  expect_equal(colnames(cor(t(master$pca))),colnames(reducedDim(scExp,"Cor")) )
  expect_equal(rownames(cor(t(master$pca))),rownames(reducedDim(scExp,"Cor")) )

})

scExp_cf = filter_correlated_cell_scExp(scExp,random_iter = 50)

test_that("Step 9 : correlation filterin clust", {
  
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_cor_filtered.RData", master)
  load("z.RData", master)
  expect_equal(nrow(master$annot_sel),ncol(scExp_cf))
  expect_equal(master$annot_sel$cell_id,colData(scExp_cf)$cell_id)
  expect_equal(table(master$annot_sel$sample_id),table(colData(scExp_cf)$sample_id))
})


scExp_cf = consensus_clustering_scExp(scExp_cf,prefix = "",reps = 1000, seed = 3.14)

test_that("Step 10 : correlation consensus hierarchical clustering", {
  
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_consclust.RData", master)
  
  # consclust list
  expect_equal(master$consclust[[2]]$consensusMatrix, scExp_cf@metadata$consclust[[2]]$consensusMatrix)
  expect_equal(master$consclust[[5]]$consensusMatrix, scExp_cf@metadata$consclust[[5]]$consensusMatrix)
  expect_equal(master$consclust[[8]]$consensusClass, scExp_cf@metadata$consclust[[8]]$consensusClass)
  
  #icl list
  expect_equal(master$icl$clusterConsensus, scExp_cf@metadata$icl$clusterConsensus)
  expect_equal(master$icl$itemConsensus, scExp_cf@metadata$icl$itemConsensus)
})

scExp_cf = choose_cluster_scExp(scExp_cf, nclust = 2)

test_that("Step 11 : correlation consensus hierarchical clustering - choose cluster", {
  
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_affectation_k2.RData", master)
  
  # cell affectation to clusters
  affectation. = as.data.frame(colData(scExp_cf))[,c("cell_id","sample_id","chromatin_group")]
  colnames(affectation.)[3] = "ChromatinGroup"
  expect_equal(master$affectation[,c(1,3,2)],  affectation. )
  expect_equal(table(master$affectation[,c("sample_id","ChromatinGroup")]),  table(affectation.[,c("sample_id","ChromatinGroup")]) )
  
  # tsne
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_tsne_filtered.RData", master)
  rownames(master$tsne_filtered$Y) = master$affectation$cell_id
  colnames(master$tsne_filtered$Y) = c("Component_1","Component_2")
  expect_equal(as.data.frame(master$tsne_filtered$Y), as.data.frame(reducedDim(scExp_cf,"TSNE")))
  
})

scExp_cf = differential_analysis_scExp(scExp_cf, de_type = "one_vs_rest",qval.th = 0.4, cdiff.th = 0.3, block = NULL)

test_that("Step 12 : differential analysis between clusters", {
  
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_2_0.4_0.3_one_vs_rest.RData", master)
  
  my.res_original = as.data.frame(apply(master$my.res_save, MARGIN = 2, as.character), stringsAsFactors=F)
  my.res_new = as.data.frame(apply(scExp_cf@metadata$diff$res, MARGIN = 2, as.character), stringsAsFactors=F)
  
  head(my.res_original)
  head(my.res_new)
  my.res_new = my.res_new[match(my.res_original$ID,my.res_new$ID),]
  rownames(my.res_new) = rownames(my.res_original) = my.res_original$ID
  
  master$norm_mat = master$norm_mat[,as.character(master$affectation$cell_id)]
  expect_equal(master$norm_mat,as.matrix(SingleCellExperiment::normcounts(scExp_cf)))
  
  expect_equal(my.res_original, my.res_new)
  expect_equal(master$summary_save, scExp_cf@metadata$diff$summary)
  expect_equal(my.res_original, my.res_new)
  
})

scExp_cf = gene_set_enrichment_analysis_scExp(scExp_cf,ref = "hg38", qval.th = 0.4, cdiff.th = 0.3,
                                              use_peaks = F )

test_that("Step 12 : GSEA of genes associated to differential loci:", {
  
  load("Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_2_0.4_0.3_one_vs_rest_GSEA.RData", master)

  # GSEA
  expect_equal(master$Both,scExp_cf@metadata$enr$Both)
  expect_equal(master$Overexpressed[[1]],scExp_cf@metadata$enr$Overexpressed[[1]])
  expect_equal(master$Overexpressed[[2]],scExp_cf@metadata$enr$Overexpressed[[2]])
  expect_equal(master$Underexpressed[[1]],scExp_cf@metadata$enr$Underexpressed[[1]])
  expect_equal(master$Underexpressed[[2]],scExp_cf@metadata$enr$Underexpressed[[2]])

})
