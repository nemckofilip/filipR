#' Map Reads with minimap2
#'
#' @description
#' Builds a shell command to align reads with minimap2 and returns it without
#' executing. Suited to long noisy reads (Nanopore, PacBio) but also handles
#' short queries against a small reference via the `preset` and `kmer`/`window`
#' arguments. Emits PAF by default, which is what most assignment work needs;
#' set `output = "bam"` for a sorted, indexed BAM.
#'
#' Submit with [fn_submit()], or for a small job run it directly with
#' `system(cmd$cmd)` — like every other wrapper here it only builds the string.
#'
#' @param fq Character. Path to the query FASTQ/FASTA (gzip supported).
#' @param reference Character. Path to the reference FASTA.
#' @param base.name Character. Sample base name used for output file naming.
#' @param output.dir Character. Directory for output. Default: `"db/paf/"`.
#' @param preset Character. minimap2 `-x` preset, e.g. `"map-ont"`,
#'   `"map-pb"`, `"sr"`, `"asm5"`. Default: `"map-ont"`.
#' @param output Character. `"paf"` (default) or `"bam"`. `"bam"` adds `-a`,
#'   pipes through `samtools sort` and indexes the result.
#' @param cores Integer. Threads passed to `-t`. Default: `8`.
#' @param kmer Integer or `NULL`. Passed to `-k`. Lower values (e.g. `11`) help
#'   when queries are short or error-rich. Default: `NULL` (preset default).
#' @param window Integer or `NULL`. Passed to `-w`. Default: `NULL`.
#' @param secondary Logical. Keep secondary alignments. `TRUE` is required if
#'   you want a best-vs-second-best ambiguity margin. Default: `FALSE`.
#' @param max.secondary Integer. Passed to `-N` when `secondary = TRUE`.
#'   Default: `5`.
#' @param min.secondary.ratio Numeric or `NULL`. Passed to `-p`: report
#'   secondaries scoring at least this fraction of the best. Default: `NULL`.
#' @param cigar Logical. Add `-c` to emit base-level CIGAR in PAF. Ignored for
#'   BAM output. Default: `FALSE`.
#' @param extra Character or `NULL`. Extra flags passed verbatim.
#'   Default: `NULL`.
#'
#' @return A `data.table` with columns: `fq_in`, `reference`, `out`, `cmd`,
#'   `path`.
#'
#' @examples
#' \dontrun{
#' # Nanopore reads to a plasmid reference, PAF out, run locally
#' cmd <- fn_minimap2_map(fq = "reads.fastq.gz", reference = "refs.fa",
#'                        base.name = "ppLS1", output.dir = "db/paf/")
#' system(cmd$cmd)
#'
#' # short error-rich inserts, keep secondaries for an ambiguity margin
#' cmd <- fn_minimap2_map(fq = "inserts.fasta", reference = "catalog.fa",
#'                        base.name = "ppLS1", kmer = 11, window = 5,
#'                        secondary = TRUE, max.secondary = 20,
#'                        min.secondary.ratio = 0.4)
#' fn_submit(cmd, job.name = "mm2", logs = "db/logs/minimap2/")
#' }
#' @import data.table
#' @export
fn_minimap2_map <- function(fq,
                            reference,
                            base.name,
                            output.dir          = "db/paf/",
                            preset              = "map-ont",
                            output              = "paf",
                            cores               = 8,
                            kmer                = NULL,
                            window              = NULL,
                            secondary           = FALSE,
                            max.secondary       = 5,
                            min.secondary.ratio = NULL,
                            cigar               = FALSE,
                            extra               = NULL) {

  # ---- Input validation ----
  if (length(fq) != 1) stop("fn_minimap2_map processes one sample at a time.")
  if (!output %in% c("paf", "bam")) stop("output must be 'paf' or 'bam'.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Output path ----
  out <- file.path(output.dir,
                   paste0(base.name, "_", gsub("[^A-Za-z0-9]", "", preset),
                          if (output == "bam") ".bam" else ".paf"))

  # ---- Optional flags ----
  kmer_flag      <- if (!is.null(kmer))   paste("-k", kmer)   else ""
  window_flag    <- if (!is.null(window)) paste("-w", window) else ""
  cigar_flag     <- if (isTRUE(cigar) && output == "paf")     "-c" else ""
  # --secondary=no is the quiet default; secondaries are only useful when a
  # best-vs-second-best comparison is wanted downstream
  secondary_flag <- if (isTRUE(secondary)) {
    paste("--secondary=yes -N", max.secondary,
          if (!is.null(min.secondary.ratio))
            paste("-p", min.secondary.ratio) else "")
  } else "--secondary=no"

  # ---- Build command ----
  cmd <- paste("minimap2",
               "-x", preset,
               "-t", cores,
               kmer_flag, window_flag, cigar_flag, secondary_flag,
               if (output == "bam") "-a" else "",
               if (!is.null(extra)) extra else "",
               shQuote(reference), shQuote(fq))

  if (output == "bam") {
    sort_tmp <- paste0(out, "_tmp")
    cmd <- paste(paste0("rm -f ", shQuote(sort_tmp), "*.bam"), "&&", cmd,
                 "| samtools sort -@ 2 -T", shQuote(sort_tmp),
                 "-o", shQuote(out), "-",
                 "&&", "samtools index", shQuote(out))
  } else {
    cmd <- paste(cmd, ">", shQuote(out))
  }

  # ---- Return ----
  data.table::data.table(
    fq_in     = fq,
    reference = reference,
    out       = out,
    cmd       = cmd,
    path      = out
  )
}
