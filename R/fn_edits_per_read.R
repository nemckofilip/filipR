#' Count A>G and C>T substitutions per read from a BAM file
#'
#' @description
#' Parses the MD tag of each mapped primary read and counts A>G and C>T
#' substitutions alongside the total reference A and C bases in the aligned
#' region. Results are aggregated into a frequency table. Generates a command
#' to run as a SLURM job via [fn_submit()].
#'
#' Output columns: `sample`, `region`, `n_AG`, `n_CT`, `n_ref_A`, `n_ref_C`,
#' `n_reads`. Per-read editing fractions can be computed as `n_AG / n_ref_A`
#' and `n_CT / n_ref_C`.
#'
#' `region` is `"ercc"` for reads on ERCC spike-in chromosomes
#' (sequencing-error baseline) and `"human"` for all other chromosomes.
#'
#' @param bam Character. Path to sorted BAM file.
#' @param base.name Base name for the output file.
#' @param output.dir Output directory.
#'   Default `"db/alignment_stats/edits_per_read/"`.
#'
#' @return A one-row `data.table` with columns: `bam`, `output`, `cmd`, `path`.
#' @export
fn_edits_per_read <- function(bam,
                               base.name,
                               output.dir = "db/alignment_stats/edits_per_read/") {

  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("Rscript", "count_edits_per_read.R", package = "filipR")
  output <- file.path(output.dir, paste0(base.name, "_edits_per_read.tsv"))

  cmd <- paste(
    "Rscript", shQuote(script),
    shQuote(bam),
    shQuote(base.name),
    shQuote(output)
  )

  data.table::data.table(bam = bam, output = output, cmd = cmd, path = output)
}
