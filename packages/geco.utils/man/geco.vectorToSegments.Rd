\name{geco.vectorToSegments}
\alias{geco.vectorToSegments}
\title{geco.vectorToSegments}
\description{Function to convert a numeric vector in a series of segments}
\usage{geco.vectorToSegments(v,replaceNA=999,d=NULL)}
\arguments{
\item{v}{Numeric vector to convert into a matrix of segments}
\item{replaceNA}{Numeric value that temporarily replaces NAs in the function. Should be different from any value in v.}
\item{d}{Data frame with other characteristics to add to the result. Each column should have a unique value for each continuous segment in v.}
}
\value{Data frame giving the coordinates of segment boundaries and values within each segment from v, plus additional characteristics if d was supplied.
}
\author{Eric Letouze}
\keyword{methods}
