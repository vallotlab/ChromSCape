#' Data.frame of chromosome length - hg38
#' 
#' This data frame provides the length of each "canonical" chromosomes of
#'  Homo Sapiens genome build hg38. 
#'
#' @usage data("hg38.chromosomes")
#'
#' @format hg38.chromosomes - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome - character}
#'   \item{start}{Start of the chromosome (bp) - integer}
#'   \item{end}{End of the chromosome (bp) - integer}
#' }
"hg38.chromosomes"

#' Data.frame of cytoBandlocation - hg38
#' 
#' This data frame provides the location of each cytoBands of
#'  Homo Sapiens genome build hg38. 
#'
#' @usage data("hg38.cytoBand")
#'
#' @format hg38.cytoBand - a data frame with 862 rows and 4 variables:
#' \describe{
#'   \item{chr}{Chromosome - character}
#'   \item{start}{Start of the chromosome (bp) - integer}
#'   \item{end}{End of the chromosome (bp) - integer}
#'   \item{cytoBand}{Name of the cytoBand - character}
#' }
"hg38.cytoBand"

#' Data.frame of cytoBandlocation - mm10
#' 
#' This data frame provides the location of each cytoBands of
#'  Homo Sapiens genome build mm10. 
#'
#' @usage data("mm10.cytoBand")
#'
#' @format mm10.cytoBand - a data frame with 862 rows and 4 variables:
#' \describe{
#'   \item{chr}{Chromosome - character}
#'   \item{start}{Start of the chromosome (bp) - integer}
#'   \item{end}{End of the chromosome (bp) - integer}
#'   \item{cytoBand}{Name of the cytoBand - character}
#' }
"mm10.cytoBand"

#' Data.frame of chromosome length - mm10
#'
#'This data frame provides the length of each "canonical" chromosomes of
#'  Mus Musculus (Mouse) genome build mm10. 
#'
#' @usage data("mm10.chromosomes")
#'
#' @format mm10.chromosomes - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome - character}
#'   \item{start}{Start of the chromosome (bp) - integer}
#'   \item{end}{End of the chromosome (bp) - integer}
#' }
"mm10.chromosomes"

#' Data.frame of gene TSS - hg38
#'
#' This dataframe was extracted from Gencode v25 and report the Transcription
#' Start Site of each gene in the Homo Sapiens genome build hg38.
#' 
#' @usage data("hg38.GeneTSS")
#' 
#' @format hg38.GeneTSS - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome - character}
#'   \item{start}{Start of the gene (TSS) - integer}
#'   \item{end}{End of the gene - integer}
#'   \item{gene}{Gene symbol - character}
#' }
"hg38.GeneTSS"

#' Data.frame of gene TSS - mm10
#'
#' This dataframe was extracted from Gencode v25 and report the Transcription
#' Start Site of each gene in the Mus Musculus genome build mm10 (Mouse).
#' 
#' @usage data("mm10.GeneTSS")
#' 
#' @format mm10.GeneTSS - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome name - character}
#'   \item{start}{Start of the gene (TSS) - integer}
#'   \item{end}{End of the gene - integer}
#'   \item{gene}{Gene symbol - character}
#' }
"mm10.GeneTSS"

#' A SingleCellExperiment outputed by ChromSCape
#' 
#' Data from a single-cell ChIP-seq experiment against H3K4me3 active mark from 
#' two cell lines, Jurkat B cells and Ramos T cells from Grosselin et al., 2019.
#' The count matrices, on 5kbp bins, were given to ChromSCape and the filtering
#'  parameter was set to 3% of cells active in regions and subsampled down to 
#'  150 cells per sample. After correlation filtering, the experiment is 
#'  composed of respectively 51 and 55 cells from Jurkat & Ramos and 5499 
#'  5kbp-genomic bins where signal is located.
#' 
#' The scExp is composed of :
#'  * counts and normcounts assays, PCA, UMAP, and Correlation matrix in 
#' reducedDims(scExp)
#'  * Assignation of genes to genomic bins  in rowRanges(scExp)
#'  * Cluster information in colData(scExp) correlation 
#'  * Hierarchical clustering dengogram in metadata$hc_cor
#'  * Consensus clustering raw data in metadata$consclust
#'  * Consensus clustering cluster-consensus and item consensus dataframes
#'  in metadata$icl
#'  * Differential analysis in metadata$diff
#'  * Gene Set Analysis in metadata$enr
#'  
#' @usage data("scExp")
#'  
#' @format scExp - a SingleCellExperiment with 106 cells and 5499 features
#'   (genomic bins) in hg38: \describe{ \item{chr}{A SingleCellExperiment} }
#'
#' @examples
#' data("scExp")
#' plot_reduced_dim_scExp(scExp)
#' plot_reduced_dim_scExp(scExp,color_by = "cell_cluster")
#' plot_heatmap_scExp(scExp)
#' plot_differential_volcano_scExp(scExp,cell_cluster = "C1")
#' plot_differential_summary_scExp(scExp)
"scExp"

#' A data.frame with the number of targets of each TF in ChEA3 
#' 
#' This data.frame was obtained by downloading datasets from ChEA3 database
#' (https://maayanlab.cloud/chea3/) and merging targets for :
#' * ARCHS4_Coexpression
#' * ENCODE_ChIP-seq
#' * Enrichr_Queries
#' * GTEx_Coexpression
#' * Literature_ChIP-seq
#' * ReMap_ChIP-seq
#' 
#' @references Keenan AB, Torre D, Lachmann A, Leong AK, Wojciechowicz M, Utti V, 
#' Jagodnik K, Kropiwnicki E, Wang Z, Ma'ayan A (2019)
#'  ChEA3: transcription factor enrichment analysis by orthogonal omics integration. 
#'  Nucleic Acids Research. doi: 10.1093/nar/gkz446
#' 
#' The data.frame is composed of two columns:
#'  * TF column containing the TF gene names (human)
#'  * nTargets_TF containing the number of targets for this TF in the combined
#'  database.
#'  
#' @usage data("CheA3_TF_nTargets")
#'  
#' @format CheA3_TF_nTargets - a data.frame with 1632 rows (unique TFs) and 
#' 2 columns
#'
#' @examples
#' data("CheA3_TF_nTargets")
#' head(CheA3_TF_nTargets)
#' 
"CheA3_TF_nTargets"
