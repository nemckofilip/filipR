#' Convert BAM to BED
#'
#' @description
#' Converts a sorted BAM file to BED format using bedtools.
#'
#' @param bam Character string. Path to the input BAM file.
#' @param base.name Base name for output files.
#' @param output.dir Directory for BED output. Default: 'db/bed/'.
#' @param sorted Logical. If TRUE, sort the BED output. Default TRUE.
#'
#' @return A data.table with columns: bam_in, bed, cmd, path.
#' @export
fn_bam_to_bed <- function(bam,
                          base.name,
                          output.dir = "db/bed/",
                          sorted = TRUE) {

  # ---- Input validation ----
  if (length(bam) != 1) stop("fn_bam_to_bed processes one BAM file at a time.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Output path ----
  bed_out <- file.path(output.dir, paste0(base.name, ".bed"))

  # ---- Build command ----
  cmd <- paste("bedtools bamtobed -i", shQuote(bam))
  if (sorted) {
    cmd <- paste(cmd, "| sort -k1,1 -k2,2n")
  }
  cmd <- paste(cmd, ">", shQuote(bed_out))

  # ---- Return data.table ----
  data.table::data.table(
    bam_in = bam,
    bed = bed_out,
    cmd = cmd,
    path = bed_out
  )
}