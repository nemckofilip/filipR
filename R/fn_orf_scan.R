#' Scan a genome FASTA for stop-codon-free ORFs in all 6 reading frames
#'
#' Generates a shell command to run the ORF scan as a SLURM job via
#' [fn_submit()]. Scans every primary chromosome (`chrN`, `chrNA/B`, `chrX`,
#' `chrY`, `chrM`) in all 6 reading frames for maximal stretches without stop
#' codons. No start codon is required. Scaffolds and alts are excluded
#' automatically. One chromosome is held in memory at a time.
#'
#' @param fasta_path Path to the genome FASTA file. Must be indexed (`.fai`).
#' @param output_dir Directory where the output TSV will be written.
#' @param base_name Base name for the output file, e.g. `"panTro6"` produces
#'   `panTro6_orf_scan_results.tsv`.
#' @param min_codons Minimum ORF length in codons. Default `100L`.
#'
#' @return A one-row `data.table` with columns `fasta_in`, `output`, `cmd`,
#'   and `path` (used by [fn_submit()] for skip-if-exists logic).
#'
#' @examples
#' \dontrun{
#' cmd <- fn_orf_scan(
#'   fasta_path = "/references/primate/panTro6/fasta/panTro6.fa",
#'   output_dir = "db/orf_scan/species_scans",
#'   base_name  = "panTro6"
#' )
#' fn_submit(cmd, job.name = "orf_scan_panTro6", mem = "32G", time = "08:00:00")
#' }
#'
#' @export
fn_orf_scan <- function(fasta_path, output_dir, base_name, min_codons = 100L) {
  script <- system.file("Rscript", "orf_scan.R", package = "filipR")
  output <- file.path(output_dir, paste0(base_name, "_orf_scan_results.tsv"))
  cmd    <- paste("Rscript", shQuote(script),
                  shQuote(fasta_path), shQuote(output), as.character(min_codons))
  data.table(fasta_in = fasta_path, output = output, cmd = cmd, path = output)
}
