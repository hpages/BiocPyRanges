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

  PyRanges' documentation: \url{https://biocore-NTNU.github.io/pyranges/}
}

\seealso{
  \itemize{
    \item \link[GenomicRanges]{GRanges} objects defined in
          the \pkg{GenomicRanges} package.

    \item What is the \pkg{\link{BiocPyRanges}} package about.
  }
}

\examples{
## ---------------------------------------------------------------------
## A SIMPLE EXAMPLE
## ---------------------------------------------------------------------

gr <- GRanges("chr1", IRanges(1:3, 20))
pyr <- makePyRangesFromGRanges(gr)
pyr  # PyRanges objects follow the 0-based start convention!!

makeGRangesFromPyRanges(pyr)

## ---------------------------------------------------------------------
## TYPICAL BiocPyRanges WORKFLOW
## ---------------------------------------------------------------------

gr1 <- GRanges("chr1", IRanges(c(1:3, 15, 35, 51), c(20:22, 30, 50, 60)),
               strand=c("+", "+", "-", "+", "+", "+"))
gr1

## The typical BiocPyRanges workflow has 3 steps:

## STEP 1 - Go to PyRanges space (i.e. turn your GRange object into a
##          PyRanges object):
pyr1 <- makePyRangesFromGRanges(gr1)
pyr1

## STEP 2 - Operate in PyRanges space (i.e. perform some operation
##          provided by the PyRanges module on your PyRanges object):
pyr2 <- pyr1$merge()  # equivalent of GenomicRanges::reduce()

## STEP 3 - Come back to R space (i.e. turn result back into GRanges
##          object):
gr2 <- makeGRangesFromPyRanges(pyr2)

## Sanity check:
stopifnot(identical(reduce(gr1, ignore.strand=TRUE), gr2))

## ---------------------------------------------------------------------
## PERFORMANCE
## ---------------------------------------------------------------------

## PyRanges is very fast!
random_GRanges <- function(n, width) {
    gr_starts <- sample(50e6, n, replace=TRUE)
    chroms <- c("chr1", "chr2", "chr3")
    gr_chrom <- sample(factor(chroms, levels=chroms), n, replace=TRUE)
    gr_strand <- sample(levels(strand()), n, replace=TRUE)
    GRanges(gr_chrom, IRanges(gr_starts, width=width, names=seq_len(n)),
            strand=gr_strand)
}

set.seed(123)
gr1 <- random_GRanges(3e6, width=40)  # 3 million random ranges

system.time(gr2a <- reduce(gr1, ignore.strand=TRUE))

pyr1 <- makePyRangesFromGRanges(gr1)
system.time(pyr2 <- pyr1$merge())
gr2b <- makeGRangesFromPyRanges(pyr2)

## Sanity check:
stopifnot(identical(gr2a, gr2b))
}

