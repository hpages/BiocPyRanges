\name{BiocPyRanges-package}

\alias{BiocPyRanges-package}
\alias{BiocPyRanges}
\docType{package}

\alias{install_pandas}
\alias{install_pyranges}
\alias{pandas}
\alias{pyranges}

\title{Conveniently access PyRanges functionalities from Bioconductor}

\description{
  PyRanges is a Python module for "efficient comparison of genomic
  intervals in Python".

  \pkg{BiocPyRanges} is an \R package that leverages the power of popular
  \R/Python interface \pkg{reticulate} to make it convenient to access
  PyRanges functionalities from Bioconductor.
  In particular:
  \itemize{
    \item It provides two utility functions for converting back and forth
          between \link[GenomicRanges]{GRanges} and PyRanges objects.
          See \code{?\link{conversions}}.

    \item For maximum convenience it also defines the \code{pandas} and
          \code{pyranges} global symbols (which are references to the
          corresponding Python modules) as well as utility functions
          to facilitate installation of these modules. See below.
  }
}

\usage{
## Wrappers for py_install() to facilitate installation of the Python
## modules required by BiocPyRanges:

install_pandas(method="auto", conda="auto")
install_pyranges(method="auto", conda="auto")
}

\arguments{
  \item{method, conda}{
    See \code{?\link[reticulate]{py_install}} in the \pkg{reticulate} package.
  }
}

\references{
  PyRanges: \url{https://github.com/biocore-NTNU/pyranges}

  PyRanges' documentation: \url{https://biocore-NTNU.github.io/pyranges/}

  reticulate: \url{https://cran.r-project.org/web/packages/reticulate/}
}

\seealso{
  \itemize{
    \item \code{\link[reticulate]{py_install}} in the \pkg{reticulate}
          package for installing Python modules into a virtual environment
          or Conda environment.

    \item \link[GenomicRanges]{GRanges} objects defined in
          the \pkg{GenomicRanges} package.

    \item \code{?\link{conversions}} for converting back and forth between
          \link[GenomicRanges]{GRanges} and PyRanges objects.
  }
}

