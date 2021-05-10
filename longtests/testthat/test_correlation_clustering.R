context("Testing correlation & clustering results.")

library(testthat)

out = create_scDataset_raw(featureType = "window")
mat = out$mat
annot = out$annot

scExp = create_scExp(mat,annot)

scExp = normalize_scExp(scExp)

scExp = reduce_dims_scExp(scExp)

scExp = correlation_and_hierarchical_clust_scExp(scExp)

test_that("Correlation & hierarchical clustering - Wrong inputs.", {
  expect_error(correlation_and_hierarchical_clust_scExp(NULL))
  expect_error(correlation_and_hierarchical_clust_scExp(scExp, hc_linkage = "ABC"))
  
  scExp. = scExp
  SingleCellExperiment::reducedDim(scExp., "PCA") = NULL  
  expect_error(correlation_and_hierarchical_clust_scExp(scExp.))
  
})

test_that("Correlation & hierarchical clustering - Right inputs.", {
  expect_s4_class(correlation_and_hierarchical_clust_scExp(scExp),
                  "SingleCellExperiment")

  expect_is(scExp@metadata$hc_cor,
                  "hclust")
  expect_type(as.matrix(SingleCellExperiment::reducedDim(scExp,"Cor")), "double")
})


scExp_cf. = filter_correlated_cell_scExp(scExp, random_iter = 50, percent_correlation = 0.5,
                                         corr_threshold = 99)

test_that("Filtering lowly correlated cells - Wrong inputs.", {
  expect_error(filter_correlated_cell_scExp(NULL))
  expect_error(filter_correlated_cell_scExp(scExp, percent_correlation = NULL))
  expect_error(filter_correlated_cell_scExp(scExp, corr_threshold = NULL))
  
  
  scExp. = scExp
  SingleCellExperiment::reducedDim(scExp., "Cor") = NULL  
  expect_error(filter_correlated_cell_scExp(scExp.))
  
})

test_that("Filtering lowly correlated cells - Right inputs.", {
  expect_s4_class(filter_correlated_cell_scExp(scExp),
                  "SingleCellExperiment")
  
  # No filtering :
  expect_equal(dim(filter_correlated_cell_scExp(scExp, corr_threshold = 0,verbose = FALSE )),
               dim(scExp))
  expect_equal(dim(filter_correlated_cell_scExp(scExp, percent_correlation = 0, verbose = FALSE)),
               dim(scExp))
  expect_lt(ncol(filter_correlated_cell_scExp(scExp,percent_correlation = 2, verbose = FALSE)),
                   ncol(scExp))
  
})

scExp_cf = scExp_cf.

test_that("Consensus Clustering - Wrong inputs.", {
  expect_error(consensus_clustering_scExp(NULL))

  expect_error(consensus_clustering_scExp(scExp_cf,maxK = NULL))
})

scExp_cf = consensus_clustering_scExp(scExp_cf)

test_that("Consensus Clustering - Right inputs.", {
  expect_s4_class(scExp_cf,
                  "SingleCellExperiment")

  # No filtering :
  expect_is(scExp_cf@metadata$consclust, "list")
  expect_is(scExp_cf@metadata$consclust[[1]], "matrix")
  expect_is(scExp_cf@metadata$consclust[[2]], "list")
  expect_is(scExp_cf@metadata$icl, "list")
  expect_is(scExp_cf@metadata$icl$clusterConsensus, "matrix")
  expect_is(scExp_cf@metadata$icl$itemConsensus, "data.frame")
  
})


test_that("Choosing clusters - Wrong inputs.", {
  expect_error(choose_cluster_scExp(NULL))
  
  expect_error(choose_cluster_scExp(scExp_cf, nclust = NULL))
  
  
  scExp_cf. = scExp_cf
  scExp_cf.@metadata$consclust = NULL
  expect_error(choose_cluster_scExp(scExp_cf.))
  
})

scExp_cf = choose_cluster_scExp(scExp_cf, nclust = 3)

test_that("Choosing clusters - Right inputs.", {
  expect_s4_class(choose_cluster_scExp(scExp_cf),
                  "SingleCellExperiment")
  
  # No filtering :
  expect_is(SingleCellExperiment::colData(scExp_cf)$cell_cluster, "character")
  expect_is(SingleCellExperiment::colData(scExp_cf)$cell_cluster_color, "character")
  expect_is(scExp_cf@metadata$hc_cor, "hclust")
})


