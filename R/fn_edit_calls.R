#' Extract per-observation edit calls from a BAM file
#'
#' @description
#' Walks every read once (via `inst/python/extract_edit_calls.py`, pysam) and
#' emits one row per reference target base it covers, keeping read identity so
#' per-read edit counts and per-site editing both derive from the same table:
#'
#' `read_id`, `chrom`, `pos` (0-based), `strand`, `is_edited` (0/1).
#'
#' `read_id` is an integer per fragment (mates share it), so you can filter
#' reads by edit count (e.g. keep reads with >= 5 edits) and still recover
#' exactly where those edits fall — the two views are linked. The target base
#' is set by `edit.type` at the transcript level (A>G -> genomic A/G on +
#' strand, T/C on -). Reference bases come from `fasta` (BAM sequences are
#' original, not 3N-converted). Submit with [fn_submit()].
#'
#' @param bam Character. Path to a sorted BAM.
#' @param base.name Base name for the output file.
#' @param output.dir Output directory. Default `"db/alignment_stats/edit_calls/"`.
#' @param fasta Path to reference FASTA (must have a .fai index).
#' @param edit.type Single transcript-level edit, e.g. `"A>G"` or `"C>T"`.
#' @param region.bed Optional BED restricting to reads overlapping those
#'   regions (e.g. functional gene bodies). Default `NULL` = whole BAM.
#' @param min.bq Minimum base quality. Default 20.
#' @param min.mq Minimum mapping quality. Default 20.
#'
#' @return A one-row `data.table` with columns: `bam`, `output`, `cmd`, `path`.
#' @export
fn_edit_calls <- function(bam,
                          base.name,
                          output.dir = "db/alignment_stats/edit_calls/",
                          fasta,
                          edit.type,
                          region.bed = NULL,
                          min.bq     = 20,
                          min.mq     = 20) {

  if (length(edit.type) != 1 || !grepl("^[ACGT]>[ACGT]$", edit.type))
    stop("edit.type must be a single edit like 'A>G' or 'C>T'.")
  if (!is.null(region.bed) && !file.exists(region.bed))
    stop("region.bed does not exist: ", region.bed)
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("python", "extract_edit_calls.py", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_edit_calls.tsv.gz"))

  worker <- paste(
    "python3", shQuote(script),
    "--edit-type", shQuote(edit.type),
    "--fasta",     shQuote(fasta),
    "--min-bq",    min.bq,
    "--min-mq",    min.mq,
    "--output",    shQuote(out)
  )

  if (is.null(region.bed)) {
    cmd <- paste(worker, "--bam", shQuote(bam))
  } else {
    # Pre-filter to reads overlapping the regions; needs a BAM index.
    tmp <- file.path(output.dir, paste0(base.name, "_edit_calls.tmp.bam"))
    cmd <- paste(
      "if [ ! -f", paste0(shQuote(bam), ".bai"), "]; then samtools index",
      shQuote(bam), "; fi &&",
      "samtools view -b -q", min.mq, "-L", shQuote(region.bed),
      shQuote(bam), "-o", shQuote(tmp), "&&",
      worker, "--bam", shQuote(tmp), "&&",
      "rm -f", shQuote(tmp)
    )
  }

  data.table::data.table(bam = bam, output = out, cmd = cmd, path = out)
}
