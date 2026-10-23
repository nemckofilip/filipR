#' Reverse-translate proteins into codon-optimised coding sequences
#'
#' @description Builds a shell command that reverse-translates every protein of a FASTA into a
#' coding sequence with DnaChisel: exact translation, codon usage matched to a species table
#' (default human), GC content window, no recognition site of the listed restriction enzymes,
#' no homopolymers and no AATAAA. Identical proteins get identical DNA, and the result is
#' deterministic (the seed is a hash of the protein). The coding sequence starts with the codon of
#' the first residue and ends with the stop codon `stop.codon`. A record that cannot be designed (letters other than the
#' 20 amino acids, empty) does not stop the run: it gets an empty cds and the reason in `error`.
#' Submit with [fn_submit()].
#' @param fasta Character. Path to the protein FASTA (one record per construct).
#' @param base.name Base name for the output file.
#' @param output.dir Output directory. Default 'db/library/'.
#' @param avoid Character vector of restriction enzyme names (as in Biopython, e.g. "XhoI") whose
#'   sites must not occur. Default none.
#' @param species Codon usage table of python_codon_tables. Default "h_sapiens".
#' @param gc.min,gc.max GC content limits inside the window. Defaults 0.35 and 0.65.
#' @param gc.window Window (nt) of the GC limits. Default 50.
#' @param homopolymer Runs of this many identical bases or more are avoided. Default 6.
#' @param stop.codon Stop codon appended to every coding sequence. Default "TAA".
#' @param prefix,suffix Constant sequences placed before the coding sequence and after the stop
#'   codon. They are not changed, but the restriction sites, homopolymers and AATAAA are also avoided
#'   across the junctions; the returned cds does not contain them. Default none.
#' @param cores Number of worker processes. Default 8.
#' @return data.table with columns: fasta, output, cmd, path. The output TSV has id, aa_len (number
#'   of residues), cds (with the stop codon), gc, max_homopolymer, translation_ok, sites (restriction sites, also counted
#'   across the junctions), constraints_ok, error.
#' @import data.table
#' @export
#' @examples
#' \dontrun{
#' cmd <- fn_reverse_translate("db/library/coat_proteins.faa", "coat", avoid = c("XhoI", "NheI"))
#' fn_submit(cmd, job.name = "reverse_translate", cores = 8, mem = "16G", time = "02:00:00")
#' }
fn_reverse_translate <- function(fasta,
                                 base.name,
                                 output.dir  = "db/library/",
                                 avoid       = character(0),
                                 species     = "h_sapiens",
                                 gc.min      = 0.35,
                                 gc.max      = 0.65,
                                 gc.window   = 50,
                                 homopolymer = 6,
                                 stop.codon  = "TAA",
                                 prefix      = "",
                                 suffix      = "",
                                 cores       = 8) {

  if (!file.exists(fasta)) stop("fasta does not exist: ", fasta)
  if (gc.min >= gc.max) stop("gc.min must be smaller than gc.max.")
  if (!toupper(stop.codon) %in% c("TAA", "TAG", "TGA")) stop("stop.codon must be TAA, TAG or TGA.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  script <- system.file("python", "reverse_translate.py", package = "filipR")
  out    <- file.path(output.dir, paste0(base.name, "_cds.tsv"))

  cmd <- paste(
    "python3", shQuote(script), shQuote(fasta),
    "--output", shQuote(out),
    "--avoid", shQuote(paste(avoid, collapse = ",")),
    "--species", shQuote(species),
    "--gc-min", gc.min, "--gc-max", gc.max, "--gc-window", gc.window,
    "--homopolymer", homopolymer,
    "--stop", shQuote(stop.codon),
    "--prefix", shQuote(prefix), "--suffix", shQuote(suffix),
    "--threads", cores
  )

  data.table::data.table(fasta = fasta, output = out, cmd = cmd, path = out)
}
