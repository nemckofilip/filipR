#' Filter Paired-End Reads by R2 Motif (cutadapt, mismatch-tolerant)
#'
#' @description
#' Keeps only read pairs where R2 starts with the specified motif.
#' Uses cutadapt for mismatch-tolerant matching. Optionally trims the motif from R2.
#'
#' @param fq1 Character vector of R1 FASTQ files.
#' @param fq2 Character vector of R2 FASTQ files.
#' @param motif Character string of the motif to match at position 1 of R2. Default: "CCATTGTCACGCTCTCTACCGGAACCAG".
#' @param mismatches Integer. Number of mismatches allowed. Default: 2.
#' @param trim Logical. If TRUE, remove the motif from R2. If FALSE, keep motif but still filter. Default: TRUE.
#' @param output.dir Directory for filtered FASTQs. Default: "db/fq/motif_filtered/".
#' @param cores Number of threads. Default: 4.
#' @param suffix Suffix to append to output file names. Default: "motifFiltered".
#'
#' @return A `data.table` with:
#' - `file.type`: fq1/fq2 labels
#' - `path`: paths to filtered FASTQ files
#' - `cmd`: shell command to perform filtering
#' - `job.name`: default job name "filter_motif_R2"
#'
#' @export
fn_filter_motif_R2 <- function(fq1, fq2,
                               motif = "CCATTGTCACGCTCTCTACCGGAACCAG",
                               mismatches = 2,
                               trim = TRUE,
                               output.dir = "db/fq_filtered/",
                               cores = 4,
                               suffix = "filtered") {

  # ---- Input validation ----
  if (length(fq1) != length(fq2)) stop("fq1 and fq2 must have the same length")
  if (!is.character(motif) || nchar(motif) == 0) stop("motif must be a non-empty string")
  if (!is.numeric(mismatches) || mismatches < 0) stop("mismatches must be a non-negative integer")
  
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Calculate error rate from mismatches ----
  error_rate <- mismatches / nchar(motif)

  # ---- Generate output FASTQ paths ----
  fq1_out <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", paste0("_", suffix, "_R1.fq.gz"), basename(fq1)))
  fq2_out <- file.path(output.dir, sub("\\.(fq|fastq)(\\.gz)?$", paste0("_", suffix, "_R2.fq.gz"), basename(fq2)))

  # ---- Build cutadapt command ----
  # -G = 5' adapter on R2 (uppercase G for read2)
  # ^ = anchor to start of read
  # --discard-untrimmed = only keep pairs where motif found
  # --pair-filter=any = discard pair if R2 doesn't match
  
  if (trim) {
    adapter_arg <- paste0("-G ^", motif)
  } else {
    # Use X prefix to mark as non-trimming (linked adapter trick)
    adapter_arg <- paste0("-G ^", motif, "X")
  }

  cmd <- paste(
    "cutadapt",
    "-j", cores,
    adapter_arg,
    "-e", round(error_rate, 3),
    "--discard-untrimmed",
    "--pair-filter=any",
    "-o", shQuote(fq1_out),
    "-p", shQuote(fq2_out),
    shQuote(fq1),
    shQuote(fq2)
  )

  # ---- Return data.table ----
  data.table::data.table(
    file.type = c(paste0("fq1.", suffix), paste0("fq2.", suffix)),
    path = c(fq1_out, fq2_out),
    cmd = cmd,
    job.name = "filter_motif_R2"
  )
}