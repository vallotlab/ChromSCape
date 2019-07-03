\name{geco.peaks}
\alias{geco.peaks}
\title{geco.peaks}
\description{wrapping of the function 'geco.density' to control the maximum number of peaks}
\usage{geco.peaks(x,percentHighestPeak=.2,maxNbPeaks=NULL,minDeltaBetweenPeaks=.03,deltaApproach=1,doplot=FALSE,...)}
\arguments{
\item{x}{a numeric vector}
\item{doplot}{boolean}
\item{percentHighestPeak}{...}
\item{maxNbPeaks}{...}
\item{minDeltaBetweenPeaks}{...}
\item{deltaApproach}{...}
\item{...}{additionnal paramters to be passed to plot}
}
\value{a list with 3 elements : top (peaks) and down (inter-peaks) points  , number of peaks
}
\author{Eric Letouze}
\keyword{methods}
