#!/usr/bin/env Rscript
# Called by fn_edit_motif(). Do not run directly.
# From a per-observation edit_calls table, tabulate the reference sequence
# context (+/- flank nt, oriented to the transcript strand) around target bases,
# split into edited vs non-edited (matched background) sets. Non-edited target
# bases are the natural null: same base, same reads, same regions. The center
# position is always the target base by construction; the signal is in the
# flanks. Contexts are fetched from the reference FASTA (BAM reads are original,
# not 3N-converted) and '-' strand contexts are reverse-complemented so every
# context reads 5'->3' on the transcript.
#
# Args: <edit_calls.tsv.gz> <fasta> <flank> <edit_min> <bg_cap> <output.rds>
#
# Per-site view only. Both classes require >= edit_min reads of coverage so they
# are coverage-matched: a site is edited if >= edit_min reads support the edit;
# background is sites with >= edit_min reads that are NEVER edited (0). Sites with
# < edit_min coverage, or the ambiguous 1..edit_min-1 edited zone, are dropped.
#
# Output RDS: list(
#   counts = list(site_edited, site_bg),   # each (2*flank+1) x 4 (A,C,G,T)
#   n      = c(site_edited, site_bg),       # set sizes
#   flank, positions)                       # positions: -flank..+flank

suppressMessages({
  library(data.table)
  library(Rsamtools)
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
})

args     <- commandArgs(trailingOnly = TRUE)
calls    <- args[1]
fasta    <- args[2]
flank    <- as.integer(args[3])
edit_min <- as.integer(args[4])
bg_cap   <- as.numeric(args[5])
out_path <- args[6]

set.seed(1)
L     <- 2L * flank + 1L
BASES <- c("A", "C", "G", "T")

# ---- Stream per-site aggregation (constant memory; the table is ~10^8 rows) ----
# awk collapses observations to one row per (chrom,pos,strand) with the edited
# and total counts, so R never holds the full per-observation table.
awk <- paste0("awk 'NR>1{k=$2\"\\t\"$3\"\\t\"$4; e[k]+=$5; t[k]++} ",
              "END{for(k in e) print k\"\\t\"e[k]\"\\t\"t[k]}'")
sites <- fread(cmd = paste("gzip -dc", shQuote(calls), "|", awk),
               header = FALSE,
               col.names  = c("chrom", "pos", "strand", "n_edited", "n_total"),
               colClasses = c("character", "integer", "character",
                              "integer", "integer"))

# Keep only the two site classes we tabulate, both requiring >= edit_min reads of
# coverage so the sets are coverage-matched: edited = >= edit_min edited reads;
# background = >= edit_min reads but NEVER edited (0). Sites with < edit_min
# coverage, or the ambiguous 1..edit_min-1 edited zone, are dropped.
sites <- sites[n_edited >= edit_min | (n_edited == 0L & n_total >= edit_min)]

# Cap background sites (keep all edited sites); per-position freqs stay unbiased.
bg_idx <- which(sites$n_edited == 0L)
if (length(bg_idx) > bg_cap) {
  drop  <- sample(bg_idx, length(bg_idx) - bg_cap)
  sites <- sites[-drop]
}

# ---- Reference context for every retained unique site (fetch strand-agnostic) ----
faidx  <- scanFaIndex(fasta)
chrlen <- setNames(width(faidx), as.character(seqnames(faidx)))

sites[, `:=`(gstart = pos + 1L - flank, gend = pos + 1L + flank)]  # pos is 0-based
ok <- sites$chrom %in% names(chrlen) &
      sites$gstart >= 1L &
      sites$gend   <= chrlen[sites$chrom]
sites <- sites[ok]

gr   <- GRanges(sites$chrom, IRanges(sites$gstart, sites$gend))
seqs <- scanFa(fasta, gr)                                # always + strand
neg  <- sites$strand == "-"
if (any(neg)) seqs[neg] <- reverseComplement(seqs[neg])  # orient to transcript
kmers <- toupper(as.character(seqs))

good  <- !grepl("[^ACGT]", kmers)                        # drop contexts with N
sites <- sites[good]; kmers <- kmers[good]

# ---- Position x base count matrix over a (weighted) set of contexts ----
base_counts <- function(km, w) {
  m <- matrix(0, nrow = L, ncol = 4L, dimnames = list(NULL, BASES))
  for (j in seq_len(L)) {
    b <- substr(km, j, j)
    a <- tapply(w, factor(b, levels = BASES), sum)
    a[is.na(a)] <- 0
    m[j, ] <- a
  }
  m
}

ed_site <- sites$n_edited >= edit_min                    # confident edited sites
bg_site <- sites$n_edited == 0L                          # never-edited sites

counts <- list(
  site_edited = base_counts(kmers[ed_site], rep(1, sum(ed_site))),
  site_bg     = base_counts(kmers[bg_site], rep(1, sum(bg_site)))
)
n <- c(site_edited = sum(ed_site), site_bg = sum(bg_site))

res <- list(counts = counts, n = n, flank = flank,
            positions = seq(-flank, flank))
saveRDS(res, out_path)
cat("motif contexts:", n["site_edited"], "edited sites,",
    n["site_bg"], "bg sites ->", out_path, "\n")
