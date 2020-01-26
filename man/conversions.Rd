\name{Bioc-PyRanges-conversions}

\alias{Bioc-PyRanges-conversions}
\alias{conversions}
\alias{makePyRangesFromGRanges}
\alias{makeGRangesFromPyRanges}

\title{Convert between GRanges and PyRanges objects}

\description{
  Two utility functions for converting back and forth between
  \link[GenomicRanges]{GRanges} and PyRanges objects.
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
}

\seealso{
  \itemize{
    \item \link[GenomicRanges]{GRanges} objects defined in
          the \pkg{GenomicRanges} package.

    \item What is the \pkg{\link{BiocPyRanges}} package about.
  }
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

