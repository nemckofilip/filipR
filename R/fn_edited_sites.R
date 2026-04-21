#' Extract per-position edited sites from a BAM file
#'
#' Runs \code{samtools mpileup} on R1 and R2 reads in parallel, parses
#' mismatches for the requested edit types, and combines the two passes into
#' a single BED with columns: chrom, start (0-based), end, name (ref>alt),
#' alt_count, strand, ref_count.
#'
#' Only positions where \code{alt_count > 0} are written. Strand and ref/alt
#' are reported at the transcript level (complement applied to reverse-strand
#' positions).
#'
#' @param bam Character. Path to sorted, indexed BAM.
#' @param base.name Base name for output files.
#' @param output.dir Output directory. Default \code{"db/edited_sites/"}.
#' @param fasta Path to reference FASTA (must have a .fai index).
#' @param edit.types Character vector of edit types to extract,
#'   e.g. \code{c("A>G", "C>T")}. Default \code{c("A>G", "C>T")}.
#' @param stranded Logical. Forward-stranded PE (e.g. CORALL). Default TRUE.
#' @param min.bq Minimum base quality. Default 20.
#' @param min.mq Minimum mapping quality. Default 20.
#' @param max.depth Maximum per-position depth cap. Default 1000000.
#'
#' @return data.table with columns: bam, output, cmd, path.
#' @export
fn_edited_sites <- function(bam,
                             base.name,
                             output.dir = "db/edited_sites/",
                             fasta,
                             edit.types = c("A>G", "C>T"),
                             stranded   = TRUE,
                             min.bq     = 20,
                             min.mq     = 20,
                             max.depth  = 1000000) {

  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script_es   <- system.file("python", "extract_edited_sites.py",
                              package = "filipR")
  script_comb <- system.file("python", "combine_edited_sites.py",
                              package = "filipR")

  out <- file.path(output.dir, paste0(base.name, "_edited_sites.bed"))
  et  <- paste(edit.types, collapse = ",")

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
      "| python3", shQuote(script_es),
      "--edit-types", shQuote(et),
      "--output",     shQuote(out)
    )
  } else {
    tmp1 <- file.path(output.dir, paste0(base.name, "_r1.tmp.bed"))
    tmp2 <- file.path(output.dir, paste0(base.name, "_r2.tmp.bed"))
    cmd <- paste(
      # R1 and R2 mpileup run in parallel
      mpileup_base, "--rf 64",  shQuote(bam),
      "| python3", shQuote(script_es),
      "--edit-types", shQuote(et),
      "--output",     shQuote(tmp1), "&",

      mpileup_base, "--rf 128", shQuote(bam),
      "| python3", shQuote(script_es), "--swap-strands",
      "--edit-types", shQuote(et),
      "--output",     shQuote(tmp2), "&",

      "wait &&",

      "python3", shQuote(script_comb),
      shQuote(tmp1), shQuote(tmp2),
      "--output", shQuote(out),
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
