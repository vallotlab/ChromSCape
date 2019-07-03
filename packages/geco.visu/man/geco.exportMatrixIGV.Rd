\name{geco.exportMatrixIGV}
\alias{geco.exportMatrixIGV}
\title{geco.exportMatrixIGV}
\description{Function to export a matrix of numerical values in a file readable withIGV}
\usage{geco.exportMatrixIGV(mat=NA,							    feature=NA,							    chrom.col=NA,							    start.col=NA,							    end.col=NA,							    probe.col=NA,							    graphType="points",							    color="0,150,50",							    altColor="0,0,0",							    igvfile=NA							    ){		}
\arguments{
\item{mat}{Matrix of numerical values to export}
\item{feature}{Feature annotation table}
\item{chrom.col}{Name of chromosome column in feature table}
\item{start.col}{Name of start position column in feature table}
\item{end.col}{Name of end position column in feature table}
\item{probe.col}{Name of probe ID column in feature table}
\item{graphType}{Type of graph (bar, points or heatmap)}
\item{color}{RGB color for positive values}
\item{altColor}{RGB color for negative values}
\item{igvfile}{File to which the output should be written}
}
\value{No value returned
}
\author{Eric Letouze}
\keyword{methods}
