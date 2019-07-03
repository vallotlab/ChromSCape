\name{geco.chromString2num}
\alias{geco.chromString2num}
\title{geco.chromString2num}
\description{convert chromosomes X, Y and MT into numerics 23, 24 and 25 }
\usage{geco.chromString2num(d, chrom="chrom",		  chromNum="chrNum")}
\arguments{
\item{d}{data.frame which columns include genomic position information }
\item{chrom}{the column in \code{d} giving the chromosome position}
}
\value{\code{d} with added \code{chromNum} column
}
\author{Eric Letouze}
\keyword{methods}
