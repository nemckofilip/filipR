#' Trim FASTQ Reads with Trim Galore
#'
#' @description
#' General-purpose wrapper for Trim Galore. Handles single-end and paired-end reads.
#' Supports adapter trimming (auto-detect or custom), hard clipping, and quality trimming.
#' Keeps read pairs synchronized when filtering.
#'
#' @param fq1 Character string. Path to R1 FASTQ file.
#' @param fq2 Character string. Path to R2 FASTQ file. Default NULL (single-end).
#' @param base.name Base name for output files.
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
#'
#' @return A `data.table` with columns: fq1_in, fq2_in, fq1_out, fq2_out, cmd.
#'
#' @examples
#' # Single-end with default adapter detection
#' cmd <- fn_trimGalore(fq1 = "sample_R1.fq.gz", base.name = "sample1")
#'
#' # Paired-end with default adapter detection
#' cmd <- fn_trimGalore(fq1 = "sample_R1.fq.gz", fq2 = "sample_R2.fq.gz", base.name = "sample1")
#'
#' # Paired-end with 20bp hard trim from 5' of R2 (e.g., RT primer removal)
#' cmd <- fn_trimGalore(fq1 = "sample_R1.fq.gz", fq2 = "sample_R2.fq.gz", 
#'                      base.name = "sample1", clip.5p.r2 = 20)
#'
#' @export
fn_trimGalore <- function(fq1, 
                          fq2 = NULL,
                          base.name,
                          output.dir = "db/fq/", 
                          clip.5p.r1 = NULL,
                          clip.5p.r2 = NULL,
                          clip.3p.r1 = NULL,
                          clip.3p.r2 = NULL,
                          adapter.r1 = NULL,
                          adapter.r2 = NULL,
                          quality = 20,
                          min.length = 20,
                          cores = 4) {

  # ---- Input validation ----
  if (length(fq1) != 1) stop("fn_trimGalore processes one sample at a time.")
  if (!is.null(fq2) && length(fq2) != 1) stop("fn_trimGalore processes one sample at a time.")
  if (!grepl("\\.(fq|fastq)(\\.gz)?$", fq1)) stop("fq1 must end with .fq, .fastq, .fq.gz, or .fastq.gz")
  if (!is.null(fq2) && !grepl("\\.(fq|fastq)(\\.gz)?$", fq2)) stop("fq2 must end with .fq, .fastq, .fq.gz, or .fastq.gz")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  is_paired <- !is.null(fq2)

  # ---- Output file paths ----
  if (is_paired) {
    fq1_out <- file.path(output.dir, paste0(base.name, "_trim_R1.fq.gz"))
    fq2_out <- file.path(output.dir, paste0(base.name, "_trim_R2.fq.gz"))
  } else {
    fq1_out <- file.path(output.dir, paste0(base.name, "_trim.fq.gz"))
  }

  # ---- Build command ----
  # Trim Galore names outputs itself, so we use its convention then rename
  # Actually: Trim Galore output naming is fixed based on input filenames.
  # We run trim_galore into output.dir, then rename to our base.name convention.
  
  tg_cmd <- paste("trim_galore", "--gzip", "-j", cores, "-q", quality, 
                  "--length", min.length, "-o", shQuote(output.dir))
  
  if (!is.null(clip.5p.r1)) tg_cmd <- paste(tg_cmd, "--clip_R1", clip.5p.r1)
  if (!is.null(clip.5p.r2)) tg_cmd <- paste(tg_cmd, "--clip_R2", clip.5p.r2)
  if (!is.null(clip.3p.r1)) tg_cmd <- paste(tg_cmd, "--three_prime_clip_R1", clip.3p.r1)
  if (!is.null(clip.3p.r2)) tg_cmd <- paste(tg_cmd, "--three_prime_clip_R2", clip.3p.r2)
  if (!is.null(adapter.r1)) tg_cmd <- paste(tg_cmd, "-a", shQuote(adapter.r1))
  if (!is.null(adapter.r2)) tg_cmd <- paste(tg_cmd, "-a2", shQuote(adapter.r2))
  
  if (is_paired) {
    tg_cmd <- paste(tg_cmd, "--paired", shQuote(fq1), shQuote(fq2))
  } else {
    tg_cmd <- paste(tg_cmd, shQuote(fq1))
  }

  # Trim Galore auto-names outputs based on input filenames, so we rename after
  if (is_paired) {
    tg_fq1 <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", "_val_1.fq.gz", basename(fq1)))
    tg_fq2 <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", "_val_2.fq.gz", basename(fq2)))
    cmd <- paste(tg_cmd, 
                 "&&", "mv", shQuote(tg_fq1), shQuote(fq1_out),
                 "&&", "mv", shQuote(tg_fq2), shQuote(fq2_out))
  } else {
    tg_fq1 <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", "_trimmed.fq.gz", basename(fq1)))
    cmd <- paste(tg_cmd, "&&", "mv", shQuote(tg_fq1), shQuote(fq1_out))
  }

  # ---- Return data.table ----
  data.table::data.table(
    fq1_in = fq1,
    fq2_in = if (is_paired) fq2 else NA_character_,
    fq1_out = fq1_out,
    fq2_out = if (is_paired) fq2_out else NA_character_,
    cmd = cmd,
    path = fq1_out
  )
}