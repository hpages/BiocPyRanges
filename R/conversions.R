### =========================================================================
### Convert back and forth between GRanges and PyRanges objects
### -------------------------------------------------------------------------


.make_PyRanges_data_frame_from_GRanges <- function(gr)
{
    stopifnot(is(gr, "GenomicRanges"))

    ## Core fields: Chromosome, Start, End
    ## PyRanges objects follow the 0-based start convention!
    df <- data.frame(Chromosome=as.factor(seqnames(gr)),
                     Start=start(gr) - 1L,
                     End=end(gr),
                     stringsAsFactors=FALSE)

    ## "Strand" field.
    gr_strand <- strand(gr)
    if (length(runValue(gr_strand)) != 1L ||
        as.character(runValue(gr_strand)) != "*" )
        df <- cbind(df, Strand=as.factor(gr_strand), stringsAsFactors=FALSE)

    ## "Name" field.
    gr_names <- names(gr)
    if (!is.null(gr_names))
        df <- cbind(df, Name=gr_names, stringsAsFactors=FALSE)

    ## Metadata columns.
    gr_mcols <- mcols(gr, use.names=FALSE)
    if (length(gr_mcols) != 0L)
        df <- cbind(df, as.data.frame(gr_mcols), stringsAsFactors=FALSE)

    df
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

    gr_seqnames <- df[[core_idx[[1L]]]]
    ## PyRanges objects follow the 0-based start convention!
    gr_start <- df[[core_idx[[2L]]]] + 1L
    gr_end <- df[[core_idx[[3L]]]]

    ## "Strand" field.
    strand_idx <- match("Strand", df_colnames)
    if (is.na(strand_idx)) {
        gr_strand <- Rle(strand("*"), nrow(df))
    } else {
        gr_strand <- df[[strand_idx]]
        used_idx <- c(used_idx, strand_idx)
    }

    ## "Name" field.
    names_idx <- match("Name", df_colnames)
    if (is.na(names_idx)) {
        gr_names <- NULL
    } else {
        gr_names <- df[[names_idx]]
        used_idx <- c(used_idx, names_idx)
    }

    ## Metadata columns.
    gr_mcols <- df[ , -used_idx, drop=FALSE]

    GRanges(gr_seqnames,
            IRanges(gr_start, gr_end, names=gr_names),
            gr_strand,
            gr_mcols)
}

makePyRangesFromGRanges <- function(gr)
{
    df <- .make_PyRanges_data_frame_from_GRanges(gr)
    pyranges$PyRanges(df)
}

makeGRangesFromPyRanges <- function(pyr)
{
    stopifnot(inherits(pyr, "pyranges.pyranges.PyRanges"))

    df <- pyr$as_df()

    ## Not safe to use makeGRangesFromDataFrame() at the moment e.g. it
    ## will error (with "cannnot determine seqnames column unambiguously")
    ## if 'df' has columns "Chromosome" and "chromosome".
    #gr <- makeGRangesFromDataFrame(df,
    #                               seqnames.field="Chromosome",
    #                               start.field="Start",
    #                               end.field="End",
    #                               strand.field="Strand",
    #                               keep.extra.columns=TRUE,
    #                               starts.in.df.are.0based=TRUE)
    #gr_mcols <- mcols(gr, use.names=FALSE)
    #if (length(gr_mcols) != 0L) {
    #    idx <- match("Name", colnames(gr_mcols))
    #    if (!is.na(idx)) {
    #        names(gr) <- gr_mcols[[idx]]
    #        mcols(gr) <- gr_mcols[-idx]
    #    }
    #}
    #gr

    .make_GRanges_from_PyRanges_data_frame(df)
}

