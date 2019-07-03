\name{geco.pangenomPlot}
\alias{geco.pangenomPlot}
\title{geco.pangenomPlot}
\description{pangenomic plot}
\usage{geco.pangenomPlot(d=NULL,ycol=NULL, chrom="chrom",pos="meanPos",startPos="start",endPos="end",colorscol=NULL,cyto=NULL,absPos="absPos",rectangleOption=FALSE,plotCentro=TRUE,plotCentro.color="lightgrey",plotChrom=TRUE,plotnew=TRUE,plotyaxis=TRUE,plotxaxis=TRUE,yaxmark=NULL,yaxlab=NULL,chromToPlot=c(1:22,"X","Y"),xlim=NULL,CexAxis=0.8,...)}
\arguments{
\item{d}{data.frame which columns include genomic position information and the wanted y axis}
\item{ycol}{the column in \code{d} to be used as y axis}
\item{chrom}{the column in  \code{d} giving the chromosome position}
\item{pos}{the column in  \code{d} giving the base pair position   (optional)}
\item{startPos}{the column in  \code{d} giving the start base pair position}
\item{endPos}{the column in  \code{d} giving the end base pair position}
\item{colorscol}{(optional) column in \code{d} giving the color for each point from ycol or a unique color name or number (default black)}
\item{cyto}{cytoband object (optional)}
\item{absPos}{the column in  \code{d} giving the absolute base pair position (will be calculated if not present in \code{d})}
\item{rectangleOption}{boolean : if \code{TRUE} rectangles are plotted instead of points (default) - NB : for points, several 'looks' can be obtained (ex using parameter type , pch,...)}
\item{plotCentro}{boolean : should centromeric delimitation be plotted ?}
\item{plotCentro.color}{color of the centromeric delimitation}
\item{plotChrom}{ boolean : should chromosome delimitation be plotted ?}
\item{plotnew}{ boolean : should the a new graph be plotted (default=\code{TRUE}), or should something be added to the current graph (-> \code{plotnew=FALSE})}
\item{plotyaxis}{ boolean : should an Y axis be plotted ? (considered only if \code{plotnew = TRUE})}
\item{yaxmark}{ y axis 'at' parameter (considered only if \code{plotyaxis = TRUE})}
\item{yaxlab}{ y axis 'labels' parameter (considered only if \code{plotyaxis = TRUE})}
\item{chromToPlot}{a vector of the chromosomes to be plotted (default c(1:22,"X","Y"))}
\item{...}{required}
}
\value{\code{d} in genormic order, with absolute position (base pair) added
}
\author{Eric Letouze}
\keyword{graphs}
