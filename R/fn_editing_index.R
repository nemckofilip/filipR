#' Compute per-base editing index from a BAM file
#'
#' Runs `samtools mpileup` piped into a strand-correct Python parser.
#' All 12 substitution types are computed. Human and ERCC chromosomes
#' are tallied separately. The reference FASTA must cover all chromosomes
#' in the BAM (use a merged human+ERCC FASTA when the BAM was aligned to
#' a combined index).
#'
#' @param bam Character. Path to sorted, indexed BAM.
#' @param base.name Base name for output files.
#' @param output.dir Output directory.
#' @param fasta Path to reference FASTA (must have a .fai index).
#' @param min.bq Minimum base quality passed to mpileup. Default 20.
#' @param min.mq Minimum mapping quality passed to mpileup. Default 20.
#' @param max.depth Maximum per-position depth cap. Default 1000000.
#'
#' @return data.table with columns: bam, output, cmd, path.
#' @export
fn_editing_index <- function(bam,
                             base.name,
                             output.dir = "db/alignment_stats/editing_index/",
                             fasta,
                             min.bq     = 20,
                             min.mq     = 20,
                             max.depth  = 1000000) {

  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("python", "compute_editing_index.py", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_editing_index.tsv"))

  cmd <- paste(
    "samtools mpileup",
    "-Q", min.bq,
    "-q", min.mq,
    "--no-BAQ",
    "-d", max.depth,
    "-f", shQuote(fasta),
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
