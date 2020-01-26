\name{BiocPyRanges}

\alias{BiocPyRanges}

\alias{pandas}
\alias{pyranges}

\alias{makePyRangesFromGRanges}
\alias{makeGRangesFromPyRanges}

\title{Conveniently access PyRanges functionalities from Bioconductor}

\description{
  PyRanges is a Python module for "efficient comparison of genomic
  intervals in Python".

  \pkg{BiocPyRanges} is an \R package that leverages the power of popular
  \R/Python interface \pkg{reticulate} to make it convenient to access
  PyRanges functionalities from Bioconductor.
  In particular:
  \itemize{
    \item It provides two functions to switch back and forth between
          \link[GenomicRanges]{GRanges} and PyRanges objects.
    \item For maximum convenience it also defines the \code{pandas} and
          \code{pyranges} global symbols which are references to the
	  corresponding Python modules.
  }

  See below for some examples.
}

\usage{
makePyRangesFromGRanges(gr)

makeGRangesFromPyRanges(pyr)
}

\arguments{
  \item{gr}{
    A \link[GenomicRanges]{GRanges} object, or, more generally, it can be
    any \link[GenomicRanges]{GenomicRanges} derivative.
  }
  \item{pyr}{
    A PyRanges object.
  }
}

\value{
  \code{makePyRangesFromGRanges} returns a PyRanges object.

  \code{makeGRangesFromPyRanges} returns a \link[GenomicRanges]{GRanges}
  object.
}

\references{
  PyRanges: \url{https://github.com/biocore-NTNU/pyranges}

  reticulate: \url{https://cran.r-project.org/web/packages/reticulate/}
}

\examples{
gr <- GRanges("chr1", IRanges(1:5, 20))
gr

## Turn into PyRanges object:
pyr <- makePyRangesFromGRanges(gr)
pyr  # PyRanges objects follow the 0-based start convention!

## Turn back into GRanges object:
gr2 <- makeGRangesFromPyRanges(pyr)

## Sanity check:
stopifnot(identical(gr, gr2))
}

