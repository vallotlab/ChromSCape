\name{geco.getMakeRd}
\alias{geco.getMakeRd}
\title{geco.getMakeRd}
\description{generate a template in the \emph{makeRd} format to ease the production of \file{.Rd} files}
\usage{geco.getMakeRd(fun=NULL, author=NULL, keyword="methods", description="required", value="required")}
\arguments{
\item{fun}{an existing \code{R} function: fill the items in the template argument field }
\item{author}{a string: fill the template author field}
\item{keyword}{a string: fill the template value field}
\item{description}{a string: fill the template description field }
\item{value}{a string: fill the template value field }
}
\value{return a makeRd template to copy/paste as commentary in \code{R} files. 
}
\author{Eric Letouze}
\keyword{methods}
