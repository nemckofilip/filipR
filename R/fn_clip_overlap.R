#' Clip overlapping mate-pair regions with bamUtil clipOverlap
#'
#' @description
#' Soft-clips the overlapping portion of paired-end mates so that every
#' genomic position is covered by at most one mate. The lower-quality mate
#' is clipped (bamUtil default). The original CIGAR is preserved in the
#' `OC` tag.
#'
#' The BAM is first name-sorted (so mates of a pair are adjacent), run
#' through `bam clipOverlap --readName` (robust against large RNA-seq
#' insert sizes), then re-sorted by coordinate and indexed. Intermediate
#' files live in a per-job temp directory inside `output.dir` and are
#' removed at the end. Stats from clipOverlap (overlap counts, clipped
#' bases) are captured to a `.stats` file.
#'
#' @param bam Character. Path to a coordinate-sorted BAM (typically the
#'   UMI-deduplicated output).
#' @param base.name Base name for output files.
#' @param output.dir Output directory. Default
#'   `"db/bam/mapped_STAR_dedup_clipped/"`.
#' @param cores Threads for `samtools sort` / `samtools index`. Default 4.
#'
#' @return A one-row `data.table` with columns:
#'   `bam_in`, `bam_clipped`, `stats`, `cmd`, `path`.
#' @export
fn_clip_overlap <- function(bam,
                            base.name,
                            output.dir = "db/bam/mapped_STAR_dedup_clipped/",
                            cores      = 4) {

  if (length(bam) != 1) stop("fn_clip_overlap processes one file at a time.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  bam_clipped <- file.path(output.dir, paste0(base.name, "_clipped.bam"))
  stats       <- file.path(output.dir, paste0(base.name, "_clipped.stats"))

  cmd <- paste(
    "TD=$(mktemp -d -p", shQuote(output.dir), ")", "&&",
    "samtools sort -n -@", cores, "-T \"$TD/ns\"",
      "-o \"$TD/ns.bam\"", shQuote(bam), "&&",
    "bam clipOverlap",
      "--in",  "\"$TD/ns.bam\"",
      "--out", "\"$TD/ns_clipped.bam\"",
      "--readName",
      "--storeOrig OC",
      "--stats",
      "2>", shQuote(stats), "&&",
    "samtools sort -@", cores, "-T \"$TD/cs\"",
      "-o", shQuote(bam_clipped), "\"$TD/ns_clipped.bam\"", "&&",
    "samtools index -@", cores, shQuote(bam_clipped), "&&",
    "rm -rf \"$TD\""
  )

  data.table::data.table(
    bam_in      = bam,
    bam_clipped = bam_clipped,
    stats       = stats,
    cmd         = cmd,
    path        = bam_clipped
  )
}
