\name{geco.makeRdPackScan}
\alias{geco.makeRdPackScan}
\title{geco.makeRdPackScan}
\description{generate a \file{.Rd} for all the CIT packages}
\usage{geco.makeRdPackScan(pack, directory=NULL, verbose=FALSE)}
\arguments{
\item{ pack}{the set of package to parse. \code{NULL} by default corresponding to all the CIT packages}
\item{directory}{the top directory of the CIT packages}
\item{verbose}{if \code{TRUE}, it lists the files created. \code{FALSE} by default}
}
\value{generate \file{.Rd} files in the \file{man} directory of each package
}
\author{Eric Letouze}
\keyword{documentation}
