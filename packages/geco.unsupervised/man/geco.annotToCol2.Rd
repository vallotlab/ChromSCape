\name{geco.annotToCol2}
\alias{geco.annotToCol2}
\title{geco.annotToCol2}
\description{Function to convert a table of annotations to a table of colors}
\usage{geco.annotToCol2(annotS=NULL, annotT=NULL, missing=c("", NA), anotype=NULL, maxnumcateg=2, categCol=NULL, quantitCol=NULL, plotLegend=T, plotLegendFile=NULL) }
\arguments{
\item{annotT}{Annotation table}
\item{annotS}{selection of interest in Annotation Table}
\item{missing}{Fields that should be considered as missing data}
\item{anotype}{Vector indicating the type of each column ("categ","binary","quantit"). If not provided, categories will be determined automatically}
\item{maxnumcateg}{Maximum number of numerical values to consider a numeric column as categorical}
\item{categcol}{Colors for categorical variables (optional)}
\item{quantitcol}{Colors for quantitive variables (optional)}
\item{plotLegend}{Should the legends be plotted in a file?}
\item{plotLegendFile}{Name of the file where the legends should be plotted}
}
\value{Annotation matrix with color codes
}
\author{Eric Letouze & Celine Vallot}
\keyword{methods}
