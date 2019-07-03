\name{geco.pangenomCoord}
\alias{geco.pangenomCoord}
\title{geco.pangenomCoord}
\description{return pangenomic absolute position, convenient for plots}
\usage{geco.pangenomCoord( x,chrom="chrom",pos="meanPos",startPos="start",   endPos="end",                         cyto=NULL,absPos="absPos",chromNum="chrNum")}
\arguments{
\item{x}{data.frame which columns include genomic position information }
\item{chrom}{the column in \code{x} giving the chromosome position (in 1:22,X,Y)}
\item{pos}{the column in \code{x} giving the base pair position   (optional)}
\item{startPos}{the column in \code{x} giving the start base pair position}
\item{endPos}{the column in \code{x} giving the end base pair position}
\item{cyto}{cytoband object (optional)}
\item{absPos}{name of added column giving absolute pangenomic order}
\item{chromNum}{the column in \code{x} giving the chromosome position with X = 23 & Y =24}
}
\value{\code{x} in genormic order, with absolute position (base pair) added
}
\author{Eric Letouze}
\keyword{graphs}
