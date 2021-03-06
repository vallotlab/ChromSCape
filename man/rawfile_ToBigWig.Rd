% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage.R
\name{rawfile_ToBigWig}
\alias{rawfile_ToBigWig}
\title{rawfile_ToBigWig : reads in BAM file and write out BigWig coverage file,
normalized and smoothed}
\usage{
rawfile_ToBigWig(
  filename,
  BigWig_filename,
  format = "BAM",
  bin_width = 150,
  n_smoothBin = 5,
  ref = "hg38",
  read_size = 101
)
}
\arguments{
\item{filename}{Path to the BAM file (with index) or BED file}

\item{BigWig_filename}{Path to write the output BigWig file}

\item{format}{File format, either "BAM" or "BED"}

\item{bin_width}{Bin size for coverage}

\item{n_smoothBin}{Number of bins for smoothing values}

\item{ref}{Reference genome.}

\item{read_size}{Length of the reads.}
}
\value{
Writes a BigWig file as output
}
\description{
rawfile_ToBigWig : reads in BAM file and write out BigWig coverage file,
normalized and smoothed
}
