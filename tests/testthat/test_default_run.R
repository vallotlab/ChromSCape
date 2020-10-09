context("Testing a A to Z run with default parameters")

# Functions for testing purposes

set.seed(47)
out = create_scDataset_raw(featureType = "window",sparse = TRUE,
                        batch_id = factor(c(1,1,2,2)))
mat = out$mat
annot = out$annot
batches = out$batches

#test sparse matrix
test_that("Sparse matrices", {
scExp = create_scExp(mat,annot)
expect_is(SingleCellExperiment::counts(scExp),"dgCMatrix")
scExp = filter_scExp(scExp)
expect_is(SingleCellExperiment::counts(scExp),"dgCMatrix")
scExp = normalize_scExp(scExp,type = "CPM")
expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
scExp=feature_annotation_scExp(scExp)
expect_is(SummarizedExperiment::rowRanges(scExp),"GRanges")
scExp = reduce_dims_scExp(scExp,n = 50,batch_correction = FALSE)
expect_is(SingleCellExperiment::reducedDim(scExp,"PCA"),"data.frame")

scExp = colors_scExp(scExp,annotCol = c("sample_id","batch_id","total_counts"))
plot_reduced_dim_scExp(scExp,reduced_dim = "PCA",color_by = "sample_id")

scExp = correlation_and_hierarchical_clust_scExp(scExp)
expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
scExp = filter_correlated_cell_scExp(scExp)
expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
scExp = consensus_clustering_scExp(scExp)
expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
expect_is(scExp@metadata$consclust,"list")
expect_is(scExp@metadata$consclust[[2]]$consensusClass,"integer")
scExp = choose_cluster_scExp(scExp,nclust = 2)
expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
})



