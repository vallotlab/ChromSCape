\name{geco.getMakeRdInternal}
\alias{geco.getMakeRdInternal}
\title{geco.getMakeRdInternal}
\description{generate a template in the \emph{makeRd} format to ease the production of \file{.Rd} files for internal functions}
\usage{geco.getMakeRdInternal(fun=NULL, author=NULL)}
\arguments{
\item{fun}{an existing \code{R} function: fill the items in the template argument field }
\item{author}{a string: fill the template author field}
}
\value{return a makeRd template to copy/paste as commentary in \code{R} files. 
}
\author{Eric Letouze}
\keyword{methods}
