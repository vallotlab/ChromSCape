\name{geco.discretize}
\alias{geco.discretize}
\title{geco.discretize}
\description{Discretize a continous variable by specified \code{lim} cut-off(s)}
\usage{geco.discretize(x, lim, quant=FALSE, addlevels=FALSE)}
\arguments{
\item{x}{ a vector}
\item{lim}{ cut-off(s)}
\item{quant}{ if TRUE (default FALSE) the \code{x} is discretize by quantile and \code{lim} is considered as cut-off(s) for quantile, ie 0<lim<1}
\item{addlevels}{ add character levels indicating the cut-offs (i.e. for un cut-off iqq levels=c("<iqq",">=iqq") )}
}
\value{a vector of integers
}
\author{Eric Letouze}
\keyword{utilities}
