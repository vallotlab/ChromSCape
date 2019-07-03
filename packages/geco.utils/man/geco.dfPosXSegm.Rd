\name{geco.dfPosXSegm}
\alias{geco.dfPosXSegm}
\title{geco.dfPosXSegm}
\description{Function to annotate positions in a data frame (dfPos) using segments in another data frame (dfSegm)}
\usage{geco.dfPosXSegm(dfPos=NULL,							dfPos.chrom.col="chrom",							dfPos.pos.col="pos",							dfSegm=NULL,							dfSegm.chrom.col="chrom",							dfSegm.start.col="start",							dfSegm.end.col="end",							colsToAdd=NULL,							namesColsToAdd=NULL						   )}
\arguments{
\item{dfPos}{Data frame with positions to annotate}
\item{dfPos.chrom.col}{Chromosome column in dfPos}
\item{dfPos.pos.col}{Position column in dfPos}
\item{dfSegm}{Data frame with segments to use for annotating dfPos}
\item{dfSegm.chrom.col}{Chromosome column in dfSegm}
\item{dfSegm.start.col}{Start position column in dfSegm}
\item{dfSegm.end.col}{End position column in dfSegm}
\item{colsToAdd}{Names of columns in dfSegm that should be used to annotate dfPos}
\item{namesColsToAdd}{Column names to give in the columns added to dfPos}
}
\value{annotated dfPos
}
\author{Eric Letouze}
\keyword{methods}
