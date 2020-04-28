context("Testing correlation & clustering results.")

library(testthat)

# Functions for testing purposes
create_scDataset_raw <- function(cells=300,features=600,
                                 featureType = c("window","peak","gene"),
                                sparse=T, nsamp=4, ref = "hg38",batch_id = rep(1,nsamp)) {
  
  stopifnot(featureType %in% c("window","peak","gene"), ref %in% c("mm10","hg38"),
            nsamp >= 1, cells >= nsamp, features >=1, length(batch_id) == nsamp)
  
  stopifnot()
  
  set.seed(47)
  
  # Create cell names
  cell_counts = sapply(split( 1:cells , sample(nsamp, cells , repl = TRUE) ), length)
  cell_names = sample = batches = list()
  for(i in 1:length(cell_counts)) {
    cell_names[[i]] = paste0("sample_",i,"_c", 1:cell_counts[i])
    sample[[i]] = rep(paste0("sample_",i),cell_counts[i])
    batches[[i]] = rep(batch_id[i],cell_counts[i])
  }
  cell_names = as.character(unlist(cell_names))
  sample = as.character(unlist(sample))
  batches = as.numeric(unlist(batches))
  
  # Create feature names
  eval(parse(text = paste0("chr <- ChromSCape::",ref,".chromosomes"))) #load species chromsizes
  chr = GRanges(chr)
  
  if(featureType[1] == "window") {
    chr_ranges = unlist(tileGenome(setNames(width(chr),seqnames(chr)),
                                      ntile = features))[1:features] # ~constant window size
    features_names = paste(as.data.frame(chr_ranges)$seqnames,
                           as.data.frame(chr_ranges)$start,
                           as.data.frame(chr_ranges)$end, sep="_")
  }
  if(featureType[1] == "peak") {
    size_peaks = c(1000,2500,7999,10000,150000,10^6) #Different size of peaks
    peaks = sapply(split( 1:features , sample(length(size_peaks), features , repl = TRUE) ), length)
    chr_ranges_list = GRangesList()
    for(i in 1:length(peaks)){
      chr_ranges = unlist(tileGenome(setNames(width(chr),seqnames(chr)),
                                     tilewidth = size_peaks[i], cut.last.tile.in.chrom = F))
      chr_ranges_list[[i]] = chr_ranges[sample(1:length(chr_ranges),size = peaks[i]),]
    }
    chr_ranges = GenomicRanges::sort.GenomicRanges(unlist(chr_ranges_list))[1:features]
    
    features_names = paste(as.data.frame(chr_ranges)$seqnames,
                           as.data.frame(chr_ranges)$start,
                           as.data.frame(chr_ranges)$end, sep="_")
  }
  if(featureType[1] == "gene") {
    eval(parse(text = paste0("chr <- ChromSCape::",ref,".GeneTSS"))) #load species chromsizes
    features_names = as.character(sample(chr$gene,features,replace = F))
  }
  vec = rpois(cells*features,0.5) #Add count to values > 0, iteratively
  for(i in 1:10) vec[vec >= i] = vec[vec >= i]  +  i^2*rpois(length(vec[vec >= i]),0.5)
  mat = matrix(vec, nrow = features, ncol = cells, 
               dimnames = list( features_names,cell_names))
  annot = data.frame(cell_id = cell_names,
                     sample_id = sample,
                     batch_id = batches,
                     total_counts = Matrix::colSums(mat))
  if(sparse) return(list("mat" =  as(mat,"dgCMatrix"), "annot" = annot)) else return(list("mat" =  mat, "annot" = annot))
}

out = create_scDataset_raw(featureType = "window")
mat = out$mat
annot = out$annot

scExp = create_scExp(mat,annot)

scExp = normalize_scExp(scExp)

scExp = reduce_dims_scExp(scExp)

scExp = correlation_and_hierarchical_clust_scExp(scExp)

test_that("Correlation & hierarchical clustering - Wrong inputs.", {
  expect_error(correlation_and_hierarchical_clust_scExp(NULL))
  expect_error(correlation_and_hierarchical_clust_scExp(scExp, correlation = "ABC"))
  expect_error(correlation_and_hierarchical_clust_scExp(scExp, hc_linkage = "ABC"))
  
  scExp. = scExp
  reducedDim(scExp., "PCA") = NULL  
  expect_error(correlation_and_hierarchical_clust_scExp(scExp.))
  
})

test_that("Correlation & hierarchical clustering - Right inputs.", {
  expect_s4_class(correlation_and_hierarchical_clust_scExp(scExp),
                  "SingleCellExperiment")
  expect_s4_class(correlation_and_hierarchical_clust_scExp(scExp, correlation = "spearman", hc_linkage = "mcquitty"),
                  "SingleCellExperiment")
  expect_s4_class(correlation_and_hierarchical_clust_scExp(scExp, correlation = "kendall", hc_linkage = "median"),
                  "SingleCellExperiment")

  expect_is(scExp@metadata$hc_cor,
                  "hclust")
  expect_type(reducedDim(scExp,"Cor"),
                  "double")
  expect_is(reducedDim(scExp,"Cor"),
                  "matrix")
})


scExp_cf. = filter_correlated_cell_scExp(scExp, random_iter = 50, percent_correlation = 2,
                                         corr_threshold = 99)

test_that("Filtering lowly correlated cells - Wrong inputs.", {
  expect_error(filter_correlated_cell_scExp(NULL))
  expect_error(filter_correlated_cell_scExp(scExp, percent_correlation = NULL))
  expect_error(filter_correlated_cell_scExp(scExp, corr_threshold = NULL))
  
  
  scExp. = scExp
  reducedDim(scExp., "Cor") = NULL  
  expect_error(filter_correlated_cell_scExp(scExp.))
  
})

test_that("Filtering lowly correlated cells - Right inputs.", {
  expect_s4_class(filter_correlated_cell_scExp(scExp),
                  "SingleCellExperiment")
  
  # No filtering :
  expect_equal(dim(filter_correlated_cell_scExp(scExp, corr_threshold = 0,verbose = F )),
               dim(scExp))
  expect_equal(dim(filter_correlated_cell_scExp(scExp, percent_correlation = 0, verbose = F)),
               dim(scExp))
  expect_lt(ncol(filter_correlated_cell_scExp(scExp,percent_correlation = 5, verbose = F)),
                   ncol(scExp))
  
})

scExp_cf = scExp_cf.

test_that("Consensus Clustering - Wrong inputs.", {
  expect_error(consensus_clustering_scExp(NULL))

  expect_error(consensus_clustering_scExp(scExp_cf,maxK = NULL))

  scExp_cf. = scExp_cf
  reducedDim(scExp_cf., "Cor") = NULL  
  expect_error(consensus_clustering_scExp(scExp_cf.))
  
})

scExp_cf = consensus_clustering_scExp(scExp_cf, prefix = "tests/test_scChIP/consensus_test",plot_consclust = "pdf")

test_that("Consensus Clustering - Right inputs.", {
  expect_s4_class(scExp_cf,
                  "SingleCellExperiment")

  expect_equal(file.exists("tests/test_scChIP/consensus_test/consensus.pdf"),TRUE)
  expect_equal(file.exists("tests/test_scChIP/consensus_test/icl001.png"),TRUE)
  if(dir.exists("tests/test_scChIP/consensus_test")) unlink("tests/test_scChIP/consensus_test",recursive = T)
  
  # No filtering :
  expect_is(scExp_cf@metadata$consclust, "list")
  expect_is(scExp_cf@metadata$consclust[[1]], "matrix")
  expect_is(scExp_cf@metadata$consclust[[2]], "list")
  expect_is(scExp_cf@metadata$consclust[[2]]$consensusMatrix, "matrix")
  expect_is(scExp_cf@metadata$consclust[[2]]$ml, "matrix")
  expect_is(scExp_cf@metadata$icl, "list")
  expect_is(scExp_cf@metadata$icl$clusterConsensus, "matrix")
  expect_is(scExp_cf@metadata$icl$itemConsensus, "matrix")
  
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
  expect_is(colData(scExp_cf)$chromatin_group, "character")
  expect_is(colData(scExp_cf)$chromatin_group_color, "character")
  expect_is(reducedDim(scExp_cf,"ConsensusAssociation"), "matrix")
  expect_is(scExp_cf@metadata$hc_consensus_association, "hclust")
})


