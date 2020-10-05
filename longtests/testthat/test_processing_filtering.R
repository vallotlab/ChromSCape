context("Testing preprocessing, filtering & reduction functions")

# Functions for testing purposes

set.seed(47)
out = create_scDataset_raw(featureType = "window",sparse = TRUE,
                           batch_id = factor(c(1,1,2,2)))
mat = out$mat
annot = out$annot
batches = out$batches

#test sparse matrix + batch correction
test_that("Sparse matrices + Batch Correction", {
    scExp = create_scExp(mat, annot)
    expect_is(SingleCellExperiment::counts(scExp),"dgCMatrix")
    scExp = filter_scExp(scExp)
    expect_is(SingleCellExperiment::counts(scExp),"dgCMatrix")
    scExp = normalize_scExp(scExp)
    expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
    scExp = reduce_dims_scExp(scExp, n = 50, batch_correction = TRUE,
                              batch_list = list(
                                  "batch_1"=c("sample_1","sample_2"),
                                  "batch_2"=c("sample_3","sample_4")))
    expect_is(SingleCellExperiment::reducedDim(scExp,"PCA"),"data.frame")
    
    scExp = colors_scExp(scExp,annotCol = c(
        "sample_id","batch_id","total_counts"))
    plot_reduced_dim_scExp(scExp,reduced_dim = "PCA",color_by = "sample_id")
    
})

#### create_scExp Wrapper to create the single cell experiment from sparse
#datamatrix and annot, remove Zero count Features then Cells and calculate QC
#Metrics (scater::calculateQCMetrics)

test_that("Wrong input - basic", {
    expect_error(create_scExp(c(),list()))
    expect_error(create_scExp(-1, ))
    expect_error(create_scExp(NA, NA))
})

test_that("Wrong input - advanced", {
    expect_error(create_scExp(annot, mat))
    expect_error(create_scExp(mat, annot, remove_zero_cells = 3))
    expect_error(create_scExp(mat, annot, remove_zero_features = 3))
    expect_error(create_scExp(mat, annot[1:5,]))
})

test_that("Some cells are empty", {
    mat. = mat
    mat.[,sample(1:ncol(mat.),3)] = 0
    annot. = annot
    annot.$total_counts = Matrix::colSums(mat.)
    expect_output(create_scExp(mat.,annot.),
                  "cells with 0 signals were removed." )
})

test_that("Some features are empty", {
    mat. = mat
    mat.[sample(1:nrow(mat.),3),] = 0
    annot. = annot
    annot.$total_counts = Matrix::colSums(mat.)
    expect_output(create_scExp(mat.,annot.),
                  "features with 0 signals were removed." )
})

test_that("Removing chrM - non canonical", {
    mat. = mat
    rownames(mat.)[sample(1:nrow(mat.),3)] = paste0("chrM_1_",1:3)
    expect_output(create_scExp(mat.,annot),
                  "chromosome M regions were removed." )
    no_removal = create_scExp(mat,annot)
    removal = create_scExp(mat.,annot)
    expect_equal(nrow(no_removal), nrow(removal)+3)
})

test_that("Removing chrM - non canonical", {
    mat. = mat
    rownames(mat.)[sample(1:nrow(mat.),3)] = paste0("chrRandom_Unk_1_",1:3)
    expect_output(create_scExp(mat.,annot),
                  "non canonical regions were removed." )
    no_removal = create_scExp(mat,annot)
    removal = create_scExp(mat.,annot)
    expect_equal(nrow(no_removal), nrow(removal)+3)
})

scExp = create_scExp(mat,annot)

#### filter_scExp Function to filter out cells & features from sparse matrix
#based on total count per cell, number of cells 'ON' (count >= 2) in features
#and top covered cells that might be doublets

test_that("Wrong input - basic", {
    expect_error(filter_scExp(c(),list()))
    expect_error(filter_scExp(-1, ))
    expect_error(filter_scExp(NA, NA))
})

test_that("No cell filter doesn't change number cells", {
    
    expect_equal(ncol(filter_scExp(scExp,
                                   percentMin = 0,
                                   quant_removal = 100,
                                   min_cov_cell = 0)), ncol(scExp) )
})

test_that("No feature filter doesn't change number features", {
    
    expect_equal(nrow(filter_scExp(scExp,
                                   percentMin = 0, quant_removal = 100, 
                                   min_cov_cell = 0)), nrow(scExp))
})

test_that("Max cell filters remove all cells", {
    
    expect_equal(
        ncol(filter_scExp(scExp, min_cov_cell = max(
            Matrix::colSums(SingleCellExperiment::counts(scExp))) )), 0)
})

test_that("Max feature filters remove all features", {
    expect_equal(nrow(filter_scExp(scExp,percentMin = 101)),0 )
})


test_that("Verbose is on /off", {
    expect_output(filter_scExp(scExp,verbose=TRUE))
    expect_output(filter_scExp(scExp,verbose=TRUE))
    expect_invisible({scExp. = filter_scExp(scExp,verbose=FALSE)} )
})

test_that("Some cells are empty", {
    mat. = mat
    mat.[sample(1:ncol(mat.),3),] = 0
    annot. = annot
    annot.$total_counts = Matrix::colSums(mat.)
    scExp. = create_scExp(mat., annot.)
    expect_type(filter_scExp(scExp.),typeof(scExp))
})

test_that("Some features are empty", {
    mat. = mat
    mat.[sample(1:nrow(mat.),3),] = 0
    annot. = annot
    annot.$total_counts = Matrix::colSums(mat.)
    scExp. = create_scExp(mat., annot.)
    expect_type(filter_scExp(scExp.),typeof(scExp))
})

#### has_genomic_coordinates
# Function to return TRUE if can find chromosome, FALSE if not
test_that("Is not genomic coordinates", {
    scExp. = scExp 
    rownames(scExp.) = paste0("Gene", 1:nrow(scExp.))
    expect_equal(has_genomic_coordinates(scExp.),FALSE)
    rownames(scExp.) = sample(letters, nrow(scExp.), replace=TRUE)
    expect_equal(has_genomic_coordinates(scExp.),FALSE)
    rownames(scExp.) = NULL
    expect_error(has_genomic_coordinates(scExp.))
})

test_that("Is genomic coordinates", {
    expect_equal(has_genomic_coordinates(scExp),TRUE)
})


#### exclude_features
# Function to exclude genomic coordinates
# or features from a scExp object
test_that("Exclude features ", {
    scExp. = scExp 
    rownames(scExp.) = paste0("Gene", 1:nrow(scExp.))
    
    expect_equal(
        exclude_features_scExp(scExp.,features_to_exclude = data.frame(
            gene=paste0("Gene", 1:10)),by = "feature_name"),
        scExp.[-c(1:10),])
    rownames(scExp.) = rownames(scExp)
    expect_warning(
        exclude_features_scExp(scExp.,features_to_exclude = data.frame(
            gene=paste0("Gene", 1:10)),by = "feature_name"))
})

#### normalize_scExp
# Function to normalize scExp
# by library size, feature size or both

test_that("Normalize features ", {
    expect_s4_class(normalize_scExp(scExp,"TPM"),"SingleCellExperiment")
    expect_s4_class(normalize_scExp(scExp,"CPM"),"SingleCellExperiment")
    expect_s4_class(normalize_scExp(scExp,"RPKM"),"SingleCellExperiment")
    expect_s4_class(normalize_scExp(scExp,"feature_size_only"),
                    "SingleCellExperiment")
    rownames(scExp) = paste0("Gene", 1:nrow(scExp))
    expect_warning(normalize_scExp(scExp))
})

#### feature_annotation
# Function to normalize scExp
# by library size, feature size or both
test_that("Feature annotation wrong input", {
    expect_error(feature_annotation_scExp(scExp,NULL))
    expect_error(feature_annotation_scExp(scExp,data.frame()))
})


test_that("Dimensionality reduction wrong input", {
    expect_error(reduce_dims_scExp(scExp,"PCB"))
    expect_warning(reduce_dims_scExp(scExp,NULL))
    expect_error(reduce_dims_scExp(data.frame()))
})


test_that("Dimensionality reduction right input", {
    set.seed(47)
    scExp. = normalize_scExp(scExp)
    pca_1 = SingleCellExperiment::reducedDim(reduce_dims_scExp(scExp.),"PCA")
    pca_2 = as.data.frame(
        prcomp(Matrix::t(SingleCellExperiment::normcounts(scExp.)),
               retx = TRUE, center = TRUE, scale. = FALSE)$x[,1:50])
    expect_equal(abs(pca_1[,1]), abs(pca_2[,1]))
})
