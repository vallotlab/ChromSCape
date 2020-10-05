#' Data.frame of chromosome length - hg38
#'
#'
#' @format hg38.chromosomes - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome name}
#'   \item{star}{Start}
#'   \item{end}{End}
#' }
"hg38.chromosomes"

#' Data.frame of chromosome length - mm10
#'
#'
#' @format mm10.chromosomes - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome name}
#'   \item{star}{Start}
#'   \item{end}{End}
#' }
"mm10.chromosomes"

#' Data.frame of gene TSS - hg38
#'
#'
#' @format hg38.GeneTSS - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome name}
#'   \item{start}{TSS start}
#'   \item{end}{TSS end}
#'   \item{gene}{Gene symbol}
#' }
"hg38.GeneTSS"

#' Data.frame of gene TSS - mm10
#'
#'
#' @format mm10.GeneTSS - a data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{chr}{Chromosome name}
#'   \item{start}{TSS start}
#'   \item{end}{TSS end}
#'   \item{gene}{Gene symbol}
#' }
"mm10.GeneTSS"

#' A randomly generated SingleCellExperiment normalized and with reduced
#' dimensions (hg38)
#'
#' @format scExp - a SingleCellExperiment with 300 cells and 600 features
#'   (genomic bins) in hg38: \describe{ \item{chr}{A SingleCellExperiment} }
"scExp"

#' A randomly generated SingleCellExperiment normalized and with reduced
#' dimensions, processed until GSA (hg38)
#'
#' @format scExp_cf - a SingleCellExperiment with 300 cells and 600 features
#'   (genomic bins) in hg38 and processed until GSA
#'   : \describe{ \item{chr}{A SingleCellExperiment} }
"scExp_cf"
