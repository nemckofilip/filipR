#' Compute per-base editing index from a BAM file
#'
#' Runs `samtools mpileup` piped into a strand-correct Python parser.
#' All 12 substitution types are computed. Human and ERCC chromosomes
#' are tallied separately. The reference FASTA must cover all chromosomes
#' in the BAM (use a merged human+ERCC FASTA when the BAM was aligned to
#' a combined index).
#'
#' For forward-stranded paired-end libraries (e.g. CORALL), set
#' \code{stranded = TRUE}. R1 reads (sense) are processed with normal
#' strand logic; R2 reads (antisense) are processed with swapped strand
#' logic so that both contribute correctly to the transcript-level edit.
#' The two partial TSVs are combined by summing numerators and denominators.
#'
#' @param bam Character. Path to sorted, indexed BAM.
#' @param base.name Base name for output files.
#' @param output.dir Output directory.
#' @param fasta Path to reference FASTA (must have a .fai index).
#' @param stranded Logical. If TRUE, process R1 and R2 separately and
#'   combine (required for forward-stranded PE libraries such as CORALL).
#'   Default TRUE.
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
                             stranded   = TRUE,
                             min.bq     = 20,
                             min.mq     = 20,
                             max.depth  = 1000000) {

  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script_ei   <- system.file("python", "compute_editing_index.py",
                             package = "filipR")
  script_comb <- system.file("python", "combine_editing_index.py",
                             package = "filipR")
  out  <- file.path(output.dir, paste0(base.name, "_editing_index.tsv"))

  mpileup_base <- paste(
    "samtools mpileup",
    "-Q", min.bq,
    "-q", min.mq,
    "--no-BAQ",
    "-d", max.depth,
    "-f", shQuote(fasta)
  )

  if (!stranded) {
    cmd <- paste(
      mpileup_base, shQuote(bam),
      "| python3", shQuote(script_ei),
      "--sample", shQuote(base.name),
      "--output", shQuote(out)
    )
  } else {
    tmp1 <- file.path(output.dir, paste0(base.name, "_r1.tmp.tsv"))
    tmp2 <- file.path(output.dir, paste0(base.name, "_r2.tmp.tsv"))
    cmd <- paste(
      mpileup_base, "--rf 64",  shQuote(bam),
      "| python3", shQuote(script_ei),
      "--sample", shQuote(base.name), "--output", shQuote(tmp1),
      "&&",
      mpileup_base, "--rf 128", shQuote(bam),
      "| python3", shQuote(script_ei), "--swap-strands",
      "--sample", shQuote(base.name), "--output", shQuote(tmp2),
      "&&",
      "python3", shQuote(script_comb),
      shQuote(tmp1), shQuote(tmp2), "--output", shQuote(out),
      "&& rm -f", shQuote(tmp1), shQuote(tmp2)
    )
  }

  data.table::data.table(
    bam    = bam,
    output = out,
    cmd    = cmd,
    path   = out
  )
}
