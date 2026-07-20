#' Tabulate reference sequence context around edited vs non-edited bases
#'
#' @description
#' Builds a shell command that tabulates the reference sequence context
#' (+/- `flank` nt, oriented to the transcript strand) around each target base
#' in an [fn_edit_calls()] table, split into an edited set and a non-edited set.
#' The non-edited target bases are the matched null - same base, same reads,
#' same regions - so per-position enrichment (edited / background) isolates the
#' deaminase's sequence preference rather than the genome's base composition.
#'
#' A per-unique-site view is produced: one context per `chrom+pos+strand`, so a
#' site counts once no matter how many reads cover it (immune to expression
#' skew). The center position is always the target base by construction; the
#' signal lives entirely in the flanks. Reference bases come from `fasta` and the
#' `-` strand contexts are reverse-complemented so every context reads 5'->3' on
#' the transcript. The tabulation runs in `inst/Rscript/edit_motif.R`. Submit
#' with [fn_submit()].
#'
#' @param edit.calls Path to a `*_edit_calls.tsv.gz` from [fn_edit_calls()].
#' @param base.name Base name for the output file.
#' @param fasta Path to reference FASTA (must have a .fai index).
#' @param output.dir Output directory. Default `"db/alignment_stats/edit_motif/"`.
#' @param flank Number of flanking nt on each side of the target base. Default 4.
#' @param edit.min Coverage/support threshold. A site is called edited if at
#'   least this many reads support the edit; the background is coverage-matched -
#'   sites with at least this many reads that are NEVER edited (0). Under-covered
#'   sites and the ambiguous `1..edit.min-1` edited zone are excluded. Raise it
#'   (e.g. 5) to drop single-read false-positive calls. Default 1.
#' @param bg.cap Max number of unique non-edited (background) sites to retain;
#'   larger sets are randomly subsampled (per-position frequencies stay
#'   unbiased). All edited sites are always kept. Default 2e6.
#'
#' @return A one-row `data.table` with columns: `edit_calls`, `output`, `cmd`, `path`.
#' @export
fn_edit_motif <- function(edit.calls,
                          base.name,
                          fasta,
                          output.dir = "db/alignment_stats/edit_motif/",
                          flank      = 4,
                          edit.min   = 1,
                          bg.cap     = 2e6) {

  if (!file.exists(fasta)) stop("fasta does not exist: ", fasta)
  if (flank < 1) stop("flank must be >= 1.")
  if (edit.min < 1) stop("edit.min must be >= 1.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("Rscript", "edit_motif.R", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_edit_motif.rds"))

  cmd <- paste(
    "Rscript", shQuote(script),
    shQuote(edit.calls),
    shQuote(fasta),
    flank,
    edit.min,
    format(bg.cap, scientific = FALSE),
    shQuote(out)
  )

  data.table::data.table(edit_calls = edit.calls, output = out,
                         cmd = cmd, path = out)
}
