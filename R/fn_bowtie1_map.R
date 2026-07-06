#' Map Reads with Bowtie1 (PRO-seq)
#'
#' @description
#' Builds a shell command to align reads to a reference index using Bowtie1
#' (`-v` mismatch mode, unique mapping only), piping the output to a sorted BAM.
#' Preserves UMIs in read headers and keeps unaligned reads in the BAM (used by
#' [fn_umi_collapse_proseq()] for total-read statistics). Submit with [fn_submit()].
#'
#' @param fq1 Character. Path(s) to R1 FASTQ file(s). Several files are merged
#'   (comma-separated) during alignment.
#' @param index Path to the Bowtie1 index prefix (e.g. references .../bowtie1/hg38).
#' @param base.name Base name for output files.
#' @param fq2 Character. Path(s) to R2 FASTQ for paired-end data. Default NULL (single-end).
#' @param output.dir Directory for the output BAM. Default 'db/bam/'.
#' @param alignment.stats.output.dir Directory for Bowtie1 alignment stats. Default 'db/alignment_stats/'.
#' @param suffix Suffix appended to output file names (e.g. '_hg38'). Default ''.
#' @param max.mismatch Maximum number of mismatches (`-v`). Default 2.
#' @param max.ins Maximum insert size for paired-end reads (`--maxins`). Default 500.
#' @param cores Total number of threads. Default 8.
#' @return data.table with columns: fq1_in, fq2_in, bam, stats, cmd, path.
#' @import data.table
#' @export
fn_bowtie1_map <- function(fq1,
                           index,
                           base.name,
                           fq2 = NULL,
                           output.dir = "db/bam/",
                           alignment.stats.output.dir = "db/alignment_stats/",
                           suffix = "",
                           max.mismatch = 2,
                           max.ins = 500,
                           cores = 8) {

  # ---- Input validation ----
  if (length(base.name) != 1) stop("fn_bowtie1_map processes one sample at a time.")
  if (!is.null(fq2) && length(fq1) != length(fq2))
    stop("When provided, fq2 files should match fq1 files.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if (!dir.exists(alignment.stats.output.dir)) dir.create(alignment.stats.output.dir, recursive = TRUE)

  is_paired <- !is.null(fq2)

  # ---- Thread balancing ----
  sort_cores <- 2
  bt_cores <- max(1, cores - sort_cores)

  # ---- Output paths ----
  bam_out <- file.path(output.dir, paste0(base.name, suffix, ".bam"))
  stats_out <- file.path(alignment.stats.output.dir, paste0(base.name, suffix, "_stats.txt"))

  # ---- Input files string (comma-separated, merged during alignment) ----
  if (is_paired) {
    files <- paste("-1", paste(shQuote(fq1), collapse = ","),
                   "-2", paste(shQuote(fq2), collapse = ","))
  } else {
    files <- paste(shQuote(fq1), collapse = ",")
  }

  # ---- Build Bowtie1 command ----
  cmd_align <- paste(
    "bowtie -p", bt_cores,
    "-q",                # FASTQ input
    "-v", max.mismatch,  # Max mismatches
    "-m 1",              # Unique mappers only
    "--best --strata",   # Best alignment
    "--sam"              # SAM output
  )
  if (is_paired) cmd_align <- paste(cmd_align, "--maxins", max.ins)
  cmd_align <- paste(cmd_align, shQuote(index), files)

  # ---- Pipeline: Align -> Sort ----
  cmd <- paste(cmd_align, "2>", shQuote(stats_out),
               "| samtools sort -@", sort_cores, "-o", shQuote(bam_out), "-")

  # ---- Return data.table ----
  data.table::data.table(
    fq1_in = paste(fq1, collapse = ","),
    fq2_in = if (is_paired) paste(fq2, collapse = ",") else NA_character_,
    bam = bam_out,
    stats = stats_out,
    cmd = cmd,
    path = bam_out
  )
}
