% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing_filtering_reduction.R
\name{num_cell_scExp}
\alias{num_cell_scExp}
\title{Table of cells}
\usage{
num_cell_scExp(annot, datamatrix)
}
\arguments{
\item{annot}{An annotation of cells. Can be obtain through 'colData(scExp)'.}

\item{datamatrix}{A matrix of cells per regions  before filtering.}
}
\value{
A formatted kable in HTML.
}
\description{
Table of cells
}
\examples{
scExp = create_scExp(create_scDataset_raw()$mat,create_scDataset_raw()$annot)
\dontrun{num_cell_scExp(SingleCellExperiment::colData(scExp))}
}
