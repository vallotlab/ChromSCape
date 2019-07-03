\name{geco.groupsVSannot}
\alias{geco.groupsVSannot}
\title{geco.groupsVSannot}
\description{This function tests the association between all colunms in an annotation table vs 2 or more groups}
\usage{geco.groupsVSannot(groups=NULL,							   annotS=NULL,							   annotT=NULL,							   groupOrder=NA,#							   maxnumcateg=10,							   orientated=F							   )}
\arguments{
\item{groups}{Vector giving the group of each sample (ordered like annot rows)}
\item{annot}{Annotation data frame}
\item{groupOrder}{Ordering of the groups for the outputs (optional)}
\item{orientated}{Should tests for quantitative variables be orientated?}
}
\value{List containing a pval table (each annot vs all groups), a pvalgp matrix (each annot vs each group) and a pvalgpcol matrix (equivalent to pvalgp with a color code)
}
\author{Eric Letouze}
\keyword{methods}
