% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copy_number.R
\name{plot_gain_or_loss_barplots}
\alias{plot_gain_or_loss_barplots}
\title{Plot Gain or Loss of cytobands of the most variables cytobands}
\usage{
plot_gain_or_loss_barplots(scExp, cells = NULL, top = 20)

plot_gain_or_loss_barplots(scExp, cells = NULL, top = 20)
}
\arguments{
\item{scExp}{A SingleCellExperiment with "logRatio_cytoBand" reducedDim slot
filled. See  \code{\link{calculate_logRatio_CNA}}}

\item{cells}{Cell IDs of the tumor samples to}

\item{top}{Number of most variables cytobands to plot}
}
\value{
Plot the gains/lost in the selected cells of interest as multiple
barplots

Plot the gains/lost in the selected cells of interest as multiple
barplots
}
\description{
Plot Gain or Loss of cytobands of the most variables cytobands

Plot Gain or Loss of cytobands of the most variables cytobands
}
\examples{

data("scExp")
scExp = calculate_CNA(scExp,  control_samples = unique(scExp$sample_id)[1],
ref_genome="hg38", quantiles_to_define_gol = c(0.05,0.95))
plot_gain_or_loss_barplots(scExp, cells = scExp$cell_id[which(
scExp$sample_id \%in\% unique(scExp$sample_id)[2])])


data("scExp")
scExp = calculate_CNA(scExp,  control_samples = unique(scExp$sample_id)[1],
ref_genome="hg38", quantiles_to_define_gol = c(0.05,0.95))
plot_gain_or_loss_barplots(scExp, cells = scExp$cell_id[which(
scExp$sample_id \%in\% unique(scExp$sample_id)[2])])

}
