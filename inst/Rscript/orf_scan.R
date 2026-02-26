#!/usr/bin/env Rscript
# Called internally by fn_orf_scan(). Do not run directly.
# Args: <fasta_path> <output_file> <min_codons>

library(data.table)
library(Biostrings)
library(Rsamtools)

args        <- commandArgs(trailingOnly = TRUE)
fasta_path  <- args[1L]
output_file <- args[2L]
min_codons  <- as.integer(args[3L])

#==============================================================================
# .fn_find_orfs
#
# Scans all 6 reading frames of one chromosome DNAString. Returns a data.table
# of stop-codon-free stretches >= min_codons with 0-based genomic coordinates.
#==============================================================================

.fn_find_orfs <- function(seq_dna, chrom_name, min_codons) {
  seq_len <- length(seq_dna)
  results <- vector("list", 6L)

  scan_frame <- function(seq, fr, strand) {
    n_cod <- (length(seq) - fr) %/% 3L
    if(n_cod < 1L) return(NULL)

    aa_vec <- strsplit(
      as.character(translate(subseq(seq, fr + 1L, fr + n_cod * 3L),
                             no.init.codon = TRUE, if.fuzzy.codon = "X")),
      "")[[1L]]

    r        <- rle(aa_vec != "*" & aa_vec != "X")
    keep     <- r$values
    if(!any(keep)) return(NULL)

    n_codons <- r$lengths[keep]
    ends     <- cumsum(r$lengths)[keep]
    starts   <- ends - n_codons + 1L
    s0       <- fr + (starts - 1L) * 3L
    e0       <- s0 + n_codons * 3L

    if(strand == "-") { tmp <- seq_len - e0; e0 <- seq_len - s0; s0 <- tmp }

    data.table(chrom = chrom_name, strand = strand, frame = fr,
               start_0based = s0, end_0based = e0, n_codons = n_codons)
  }

  for(fr in 0:2) results[[fr + 1L]] <- scan_frame(seq_dna, fr, "+")
  rc_seq <- reverseComplement(seq_dna)
  for(fr in 0:2) results[[fr + 4L]] <- scan_frame(rc_seq, fr, "-")

  out <- rbindlist(results, use.names = TRUE)
  if(nrow(out) == 0L) return(out)
  out[n_codons >= min_codons]
}

#==============================================================================
# LOAD FASTA — PRIMARY CHROMOSOMES ONLY
# Keeps chrN, chrNA/B (e.g. chr2A in chimp), chrX, chrY, chrM.
#==============================================================================

fa          <- FaFile(fasta_path)
idx         <- scanFaIndex(fa)
chrom_names <- as.character(seqnames(idx))
chrom_lens  <- setNames(width(idx), chrom_names)
chrom_names <- chrom_names[grepl("^chr([0-9]+[AB]?|X|Y|M)$", chrom_names)]
chrom_lens  <- chrom_lens[chrom_names]

#==============================================================================
# SCAN CHROMOSOMES
#==============================================================================

all_results <- vector("list", length(chrom_names))

for(i in seq_along(chrom_names)) {
  chr              <- chrom_names[i]
  seq_dna          <- scanFa(fa, GRanges(chr, IRanges(1, chrom_lens[i])))[[1L]]
  all_results[[i]] <- .fn_find_orfs(seq_dna, chr, min_codons)
  rm(seq_dna)
  if(chrom_lens[i] > 100e6) gc(verbose = FALSE)
}

#==============================================================================
# MERGE, SORT, SAVE
#==============================================================================

results <- rbindlist(all_results, use.names = TRUE)
rm(all_results); gc(verbose = FALSE)

setorder(results, -n_codons)
results[, `:=`(rank = .I, length_nt = n_codons * 3L)]
setcolorder(results, c("rank", "chrom", "strand", "frame",
                        "start_0based", "end_0based", "n_codons", "length_nt"))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(results, output_file, sep = "\t")
