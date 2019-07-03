\name{geco.genomOrder}
\alias{geco.genomOrder}
\title{geco.genomOrder}
\description{required}
\usage{geco.genomOrder(d,chrom="chrom",pos="pos",startPos="start",endPos="end" ,chromNum="chrNum") }
\arguments{
\item{d}{a \code{data.frame}}
\item{chrom}{chromosome column name}
\item{pos}{sequence location column name}
\item{startPos}{the column in \code{d} giving the start base pair position}
\item{endPos}{the column in \code{d} giving the end base pair position}
\item{chromNum}{the name of added column whith chromosome from 1 to 24}
}
\value{\code{d} ordered genomically
}
\author{Eric Letouze}
\keyword{methods}
