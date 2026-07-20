#' Count sgRNAs Using MAGeCK count
#'
#' @description
#' Builds a shell command to count CRISPR sgRNA reads across samples with
#' `mageck count`, producing a guide x sample count matrix. All samples are
#' counted in a single job; each sample may span several lane FASTQs, which are
#' treated as technical replicates (comma-joined). Submit with [fn_submit()].
#'
#' @param fastqs A list of character vectors, one element per sample; each
#'   vector holds that sample's FASTQ paths (e.g. the 4 lanes). Lanes are
#'   comma-joined for MAGeCK (technical replicates of one sample).
#' @param sample.names Character vector of sample labels, same length and order
#'   as `fastqs`. Used for `--sample-label` and the count-matrix columns.
#' @param library.file Path to the MAGeCK library CSV (`sgRNA,sequence,gene`).
#' @param base.name Output prefix name (MAGeCK `-n`).
#' @param control.sgrna Optional path to a control sgRNA ID list for the null
#'   distribution / normalisation (`--control-sgrna`). Default NULL.
#' @param trim5 5' trimming length(s) passed to `--trim-5`; comma-separated
#'   string or "AUTO". Default "23,24".
#' @param output.dir Output directory. Default "db/counts/".
#' @return data.table with columns: counts, counts_normalized, summary, cmd, path.
#' @import data.table
#' @export
fn_mageck_count <- function(fastqs,
                            sample.names,
                            library.file,
                            base.name,
                            control.sgrna = NULL,
                            trim5 = "23,24",
                            output.dir = "db/counts/") {

  # ---- Input validation ----
  if (!is.list(fastqs)) stop("fastqs must be a list of character vectors (one per sample).")
  if (length(fastqs) != length(sample.names))
    stop("fastqs and sample.names must have the same length.")
  if (length(library.file) != 1) stop("Provide exactly one library.file.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Output paths ----
  prefix            <- file.path(output.dir, base.name)
  counts_out        <- paste0(prefix, ".count.txt")
  counts_norm_out   <- paste0(prefix, ".count_normalized.txt")
  summary_out       <- paste0(prefix, ".countsummary.txt")

  # ---- Build FASTQ / label arguments ----
  # Each sample -> its lanes comma-joined; samples space-separated.
  fastq_arg <- paste(vapply(fastqs,
                            function(v) paste(shQuote(v), collapse = ","),
                            character(1)),
                     collapse = " ")
  label_arg <- paste(sample.names, collapse = ",")

  # ---- Build command ----
  cmd <- paste(
    "mageck count",
    "--list-seq", shQuote(library.file),
    "--fastq", fastq_arg,
    "--sample-label", shQuote(label_arg),
    "-n", shQuote(prefix),
    "--trim-5", shQuote(trim5)
  )
  if (!is.null(control.sgrna))
    cmd <- paste(cmd, "--control-sgrna", shQuote(control.sgrna))

  # ---- Return data.table ----
  data.table::data.table(
    counts = counts_out,
    counts_normalized = counts_norm_out,
    summary = summary_out,
    cmd = cmd,
    path = counts_out
  )
}
