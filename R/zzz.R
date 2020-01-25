# Global reference to pandas (will be initialized in .onLoad).
pandas <- NULL
pyranges <- NULL

.onLoad <- function(libname, pkgname) {
  pandas <<- reticulate::import("pandas", delay_load=TRUE)
  pyranges <<- reticulate::import("pyranges", delay_load=TRUE)
}

.make_GRanges_from_PyRanges_data_frame <- function(df)
{
    stopifnot(is.data.frame(df))

    df_colnames <- colnames(df)

    ## Core fields: Chromosome, Start, End
    core_fields <- c("Chromosome", "Start", "End")
    core_idx <- match(core_fields, df_colnames)
    missing_idx <- which(is.na(core_idx))
    if (length(missing_idx) != 0L) {
        missing_fields <- core_fields[missing_idx]
        if (length(missing_fields) == 1L) {
            msg <- c("column ", missing_fields, " is missing")
        } else {
            missing_fields <- paste(missing_fields, collapse=", ")
            msg <- c("some columns are missing (", missing_fields, ")")
        }
        stop(wmsg("invalid PyRanges object: ", msg))
    }
    used_idx <- core_idx

    ans_seqnames <- df[[core_idx[[1L]]]]
    ## PyRanges objects follow the 0-based start convention!
    ans_start <- df[[core_idx[[2L]]]] + 1L
    ans_end <- df[[core_idx[[3L]]]]

    ## "Strand" field.
    strand_idx <- match("Strand", df_colnames)
    if (is.na(strand_idx)) {
        ans_strand <- Rle(strand("*"), nrow(df))
    } else {
        ans_strand <- df[[strand_idx]]
        used_idx <- c(used_idx, strand_idx)
    }

    ## "Name" field.
    names_idx <- match("Name", df_colnames)
    if (is.na(names_idx)) {
        ans_names <- NULL
    } else {
        ans_names <- df[[names_idx]]
        used_idx <- c(used_idx, names_idx)
    }

    ## Metadata columns.
    ans_mcols <- df[ , -used_idx, drop=FALSE]

    GRanges(ans_seqnames,
            IRanges(ans_start, ans_end, names=ans_names),
            ans_strand,
            ans_mcols)
}

makeGRangesFromPyRanges <- function(x)
{
    stopifnot(inherits(x, "pyranges.pyranges.PyRanges"))

    df <- x$as_df()

    ## Not safe to use at the moment e.g. will error (with "cannnot determine
    ## seqnames column unambiguously") if 'df' has columns "Chromosome"
    ## and "chromosome".
    #ans <- makeGRangesFromDataFrame(df,
    #                                seqnames.field="Chromosome",
    #                                start.field="Start",
    #                                end.field="End",
    #                                strand.field="Strand",
    #                                keep.extra.columns=TRUE,
    #                                starts.in.df.are.0based=TRUE)
    #ans_mcols <- mcols(ans)
    #if (!is.null(ans_mcols)) {
    #    idx <- match("Name", colnames(ans_mcols))
    #    if (!is.na(idx)) {
    #        names(ans) <- ans_mcols[[idx]]
    #        mcols(ans) <- ans_mcols[-idx]
    #    }
    #}
    #ans

    .make_GRanges_from_PyRanges_data_frame(df)
}

makePyRangesFromGRanges <- function(x)
{
    stopifnot(is(x, "GRanges"))
}

