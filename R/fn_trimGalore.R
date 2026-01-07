#' Trim FASTQ Reads with Trim Galore
#'
#' @description
#' General-purpose wrapper for Trim Galore. Handles single-end and paired-end reads.
#' Supports adapter trimming (auto-detect or custom), hard clipping, and quality trimming.
#' Keeps read pairs synchronized when filtering.
#'
#' @param fq1 A non-empty character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default NULL.
#' @param output.dir Directory for output FASTQ files. Default "db/fq/".
#' @param clip.5p.r1 Hard trim N bases from 5' end of R1. Default NULL.
#' @param clip.5p.r2 Hard trim N bases from 5' end of R2. Default NULL.
#' @param clip.3p.r1 Hard trim N bases from 3' end of R1. Default NULL.
#' @param clip.3p.r2 Hard trim N bases from 3' end of R2. Default NULL.
#' @param adapter.r1 Adapter sequence for R1. Default NULL (auto-detect).
#' @param adapter.r2 Adapter sequence for R2. Default NULL (auto-detect).
#' @param quality Quality cutoff for trimming. Default 20.
#' @param min.length Minimum read length to keep after trimming. Default 20.
#' @param cores Number of threads. Default 4.
#' @param suffix Suffix for output file.type column names. Default "trim".
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels (fq1.{suffix}, fq2.{suffix}).
#' - `path`: Paths to trimmed FASTQ files.
#' - `cmd`: Shell command to run Trim Galore.
#' - `job.name`: Default name for the job = "trimGalore".
#'
#' @examples
#' # Single-end with default adapter detection
#' cmd <- fn_trimGalore(fq1 = "sample_R1.fq.gz")
#'
#' # Paired-end with default adapter detection
#' cmd <- fn_trimGalore(fq1 = "sample_R1.fq.gz", fq2 = "sample_R2.fq.gz")
#'
#' # Paired-end with 20bp hard trim from 5' of R2 (e.g., RT primer removal)
#' cmd <- fn_trimGalore(fq1 = "sample_R1.fq.gz", fq2 = "sample_R2.fq.gz", clip.5p.r2 = 20, suffix = "clipped")
#'
#' @export
fn_trimGalore <- function(fq1, 
                          fq2 = NULL,
                          output.dir = "db/fq/", 
                          clip.5p.r1 = NULL,
                          clip.5p.r2 = NULL,
                          clip.3p.r1 = NULL,
                          clip.3p.r2 = NULL,
                          adapter.r1 = NULL,
                          adapter.r2 = NULL,
                          quality = 20,
                          min.length = 20,
                          cores = 4,
                          suffix = "trim") {

  # ---- Input validation ----
  if (!is.character(fq1) || length(fq1) == 0)
    stop("`fq1` must be a non-empty character vector.")
  if (!is.null(fq2) && (!is.character(fq2) || length(fq2) == 0))
    stop("`fq2` must be a non-empty character vector if provided.")
  if (!is.null(fq2) && length(fq1) != length(fq2))
    stop("`fq1` and `fq2` must have the same length.")
  if (any(!grepl("\\.(fq|fastq)(\\.gz)?$", c(fq1, fq2))))
    stop("File paths must end with `.fq`, `.fastq`, `.fq.gz`, or `.fastq.gz`")

  # ---- Path handling ----
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Output file paths ----
  if (is.null(fq2)) {
    fq1_out <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", "_trimmed.fq.gz", basename(fq1)))
  } else {
    fq1_out <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", "_val_1.fq.gz", basename(fq1)))
    fq2_out <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", "_val_2.fq.gz", basename(fq2)))
  }

  # ---- Build command ----
  cmd <- paste("trim_galore", "--gzip", "-j", cores, "-q", quality, "--length", min.length, "-o", shQuote(output.dir))
  
  if (!is.null(clip.5p.r1)) cmd <- paste(cmd, "--clip_R1", clip.5p.r1)
  if (!is.null(clip.5p.r2)) cmd <- paste(cmd, "--clip_R2", clip.5p.r2)
  if (!is.null(clip.3p.r1)) cmd <- paste(cmd, "--three_prime_clip_R1", clip.3p.r1)
  if (!is.null(clip.3p.r2)) cmd <- paste(cmd, "--three_prime_clip_R2", clip.3p.r2)
  if (!is.null(adapter.r1)) cmd <- paste(cmd, "-a", shQuote(adapter.r1))
  if (!is.null(adapter.r2)) cmd <- paste(cmd, "-a2", shQuote(adapter.r2))
  
  if (is.null(fq2)) {
    cmd <- paste(cmd, shQuote(fq1))
  } else {
    cmd <- paste(cmd, "--paired", shQuote(fq1), shQuote(fq2))
  }

  # ---- Return data.table ----
  if (is.null(fq2)) {
    data.table::data.table(
      file.type = paste0("fq1.", suffix),
      path = fq1_out,
      cmd = cmd,
      job.name = "trimGalore"
    )
  } else {
    data.table::data.table(
      file.type = c(paste0("fq1.", suffix), paste0("fq2.", suffix)),
      path = c(fq1_out, fq2_out),
      cmd = cmd,
      job.name = "trimGalore"
    )
  }
}