#' Read the Sequence from a SnapGene .dna File
#'
#' @description
#' Extracts the nucleotide sequence from a SnapGene `.dna` plasmid map without
#' any external tool. The format is a chain of segments, each
#' `[1 byte type][4 byte big-endian length][payload]`; segment type `0` holds
#' the sequence, whose first payload byte is a topology flag rather than a base.
#' Useful for reading cloning maps directly instead of exporting FASTA by hand.
#'
#' @param path Character. Path to a single `.dna` file.
#' @param as.character Logical. If `TRUE` (default) returns an uppercase
#'   character string; if `FALSE` returns a `Biostrings::DNAString`, which
#'   requires the Biostrings package.
#'
#' @return A character string of the full sequence, or a `DNAString` when
#'   `as.character = FALSE`.
#'
#' @examples
#' \dontrun{
#' s <- fn_read_snapgene("db/pools_identity/pHL349.dna")
#' nchar(s)
#' # locate a restriction site
#' gregexpr("GGATCC", s)[[1]]
#' }
#' @export
fn_read_snapgene <- function(path, as.character = TRUE) {

  # ---- Input validation ----
  if (length(path) != 1) stop("fn_read_snapgene reads one file at a time.")
  if (!file.exists(path)) stop("file not found: ", path)

  # ---- Walk the segment chain ----
  raw <- readBin(path, "raw", n = file.info(path)$size)
  i   <- 1L
  seq <- NULL
  while (i < length(raw) - 5L) {
    type <- as.integer(raw[i])
    n    <- sum(as.integer(raw[(i + 1L):(i + 4L)]) * 256^(3:0))
    if (type == 0L) {
      # skip the topology flag byte at i + 5
      seq <- toupper(rawToChar(raw[(i + 6L):(i + 4L + n)]))
      break
    }
    i <- i + 5L + n
  }
  if (is.null(seq)) stop("no sequence segment (type 0) in ", path)

  # ---- Return ----
  if (as.character) return(seq)
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("as.character = FALSE requires the Biostrings package.")
  Biostrings::DNAString(seq)
}
