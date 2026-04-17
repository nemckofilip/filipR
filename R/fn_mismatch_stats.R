#' Count base substitutions from a BAM file, split by human vs ERCC
#'
#' @description
#' Generates a command that streams all mapped reads from a BAM, parses
#' MD tags via a bundled Python script, and writes per-substitution-type
#' counts to two text files (human chromosomes and ERCC spike-ins).
#' Also writes positional count files and total-base count files so that
#' per-position mismatch rates can be computed downstream (rate =
#' mismatch_pos / base_pos).
#'
#' With `rf.stranded = TRUE` (e.g. Lexogen CORALL), substitutions and
#' base counts are reported from the mRNA perspective: reads from
#' minus-strand genes are complement-flipped (flip when R1 XOR
#' reverse-mapped), collapsing e.g. T>C into A>G so that each editing
#' type appears as a single substitution class regardless of gene strand.
#' Positions use R1/R2 identity for the 5-prime-of-mRNA flip.
#'
#' @param bam Path to input BAM file (must be indexed).
#' @param base.name Base name for output files.
#' @param output.dir Output directory.
#'   Default: `"db/alignment_stats/mismatch_stats/"`.
#' @param rf.stranded Logical; use mRNA-perspective counting for RF
#'   stranded PE libraries (e.g. Lexogen CORALL). Default `FALSE`.
#'
#' @return A `data.table` with columns `bam`, `path` (human mismatch
#'   output, used for existence check by [fn_submit()]), `path_human`,
#'   `path_ercc`, `path_human_pos`, `path_ercc_pos`, `path_human_base`,
#'   `path_ercc_base`, `path_human_base_pos`, `path_ercc_base_pos`,
#'   `cmd`.
#'
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' cmd <- fn_mismatch_stats("sample.bam", base.name = "sample1",
#'                          rf.stranded = TRUE)
#' fn_submit(cmd, job.name = "mismatch_sample1", mem = "8G", cores = 2)
#' }
fn_mismatch_stats <- function(bam,
                              base.name,
                              output.dir =
                                "db/alignment_stats/mismatch_stats/",
                              rf.stranded = FALSE) {

  script            <- system.file("python", "parse_md_substitutions.py",
                                   package = "filipR")
  out_human         <- file.path(output.dir,
                                 paste0(base.name, "_mismatch_human.txt"))
  out_ercc          <- file.path(output.dir,
                                 paste0(base.name, "_mismatch_ERCC.txt"))
  out_human_pos     <- file.path(output.dir,
                                 paste0(base.name, "_mismatch_human_pos.txt"))
  out_ercc_pos      <- file.path(output.dir,
                                 paste0(base.name, "_mismatch_ERCC_pos.txt"))
  out_human_base    <- file.path(output.dir,
                                 paste0(base.name, "_base_human.txt"))
  out_ercc_base     <- file.path(output.dir,
                                 paste0(base.name, "_base_ERCC.txt"))
  out_human_base_pos <- file.path(output.dir,
                                  paste0(base.name, "_base_human_pos.txt"))
  out_ercc_base_pos  <- file.path(output.dir,
                                  paste0(base.name, "_base_ERCC_pos.txt"))

  rf_flag <- if (isTRUE(rf.stranded)) "--rf-stranded" else ""

  cmd <- paste(
    "samtools view -F 2308", shQuote(bam),
    "| python3", shQuote(script),
    "--out-human",          shQuote(out_human),
    "--out-ercc",           shQuote(out_ercc),
    "--out-human-pos",      shQuote(out_human_pos),
    "--out-ercc-pos",       shQuote(out_ercc_pos),
    "--out-human-base",     shQuote(out_human_base),
    "--out-ercc-base",      shQuote(out_ercc_base),
    "--out-human-base-pos", shQuote(out_human_base_pos),
    "--out-ercc-base-pos",  shQuote(out_ercc_base_pos),
    rf_flag
  )

  data.table(
    bam                = bam,
    path               = out_human,
    path_human         = out_human,
    path_ercc          = out_ercc,
    path_human_pos     = out_human_pos,
    path_ercc_pos      = out_ercc_pos,
    path_human_base    = out_human_base,
    path_ercc_base     = out_ercc_base,
    path_human_base_pos = out_human_base_pos,
    path_ercc_base_pos  = out_ercc_base_pos,
    cmd                = cmd
  )
}
