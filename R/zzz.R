### Per recommendation of the reticulate folks, we provide convenient
### wrappers for py_install() to make it easy for the user to install the
### required Python modules. See:
### https://cran.r-project.org/web/packages/reticulate/vignettes/package.html

install_pandas <- function(method="auto", conda="auto")
    reticulate::py_install("pandas", method=method, conda=conda)

install_pyranges <- function(method="auto", conda="auto")
    reticulate::py_install("pyranges", method=method, conda=conda)

### Global references to pandas and pyranges (will be initialized in .onLoad).
pandas <- NULL
pyranges <- NULL

.onLoad <- function(libname, pkgname)
{
    ## Delay loading of Python modules (will only be loaded when
    ## accessed via $). Also per recommendation of the reticulate folks.
    pandas <<- reticulate::import("pandas", delay_load=TRUE)
    pyranges <<- reticulate::import("pyranges", delay_load=TRUE)
}

### Define methods for S3 generics r_to_py() and py_to_r() defined in
### reticulate. Also per recommendation of the reticulate folks.
#r_to_py.GenomicRanges <- function(x, convert=FALSE) makePyRangesFromGRanges(x)
#py_to_r.pyranges.pyranges.PyRanges <- function(x) makeGRangesFromPyRanges(x)

