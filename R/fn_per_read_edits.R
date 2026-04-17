#' Compute per-read editing distribution from a BAM file
#'
#' Runs `samtools view` piped into a Python parser that uses MD tags and CIGAR
#' to count all 12 substitution types per read. Assumes an RF-stranded library
#' (R1 antisense, R2 sense); mismatches are converted to the mRNA/sense-strand
#' coordinate frame using the flip condition: is_r1 != is_reverse.
#'
#' @param bam Character. Path to sorted BAM.
#' @param base.name Base name for output files.
#' @param output.dir Output directory.
#'
#' @return data.table with columns: bam, output, cmd, path.
#' @export
fn_per_read_edits <- function(bam,
                              base.name,
                              output.dir = "db/alignment_stats/per_read_edits/") {

  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("python", "compute_per_read_edits.py", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_per_read_edits.tsv"))

  cmd <- paste(
    "samtools view -F 2308",   # exclude unmapped, secondary, supplementary
    shQuote(bam),
    "| python3", shQuote(script),
    "--sample", shQuote(base.name),
    "--output", shQuote(out)
  )

  data.table::data.table(
    bam    = bam,
    output = out,
    cmd    = cmd,
    path   = out
  )
}
