#!/usr/bin/env Rscript
# Called by fn_metagene_editing(). Do not run directly.
# Bins per-observation edit calls by position along the mature transcript and
# tallies edited / total target bases per bin, split by gene-length band, giving
# two positional profiles:
#   - dist3p : editing rate vs absolute distance from the 3' end (nt bins)
#   - frac   : editing rate vs fractional 5'->3' position (metagene)
# Gene-length banding lets the plot keep a comparable gene set: each band's
# distance profile is fair only up to the band's shortest gene (capped downstream).
#
# Args: <edit_calls_annot.rds> <out.tsv> <min_gene_len> <bin_size> <max_dist> <n_frac_bins>
# Output TSV columns: len_band, scale, x, n_edited, n_total  (rate = n_edited / n_total).

suppressMessages(library(data.table))

a           <- commandArgs(trailingOnly = TRUE)
in_path     <- a[1]
out_path    <- a[2]
min_gene_len <- as.integer(a[3])
bin_size    <- as.integer(a[4])
max_dist    <- as.integer(a[5])
n_frac_bins <- as.integer(a[6])

d <- readRDS(in_path)
d <- d[gene_len >= min_gene_len]
d[, len_band := cut(gene_len, breaks = c(1000, 3000, 8000, 20000, Inf),
                    labels = c("1-3kb", "3-8kb", "8-20kb", ">20kb"),
                    right = FALSE)]

# ---- 3'-distance profile (0 .. max_dist), per length band ----
da <- d[dist_from_3p >= 0 & dist_from_3p < max_dist]
da[, x := (dist_from_3p %/% bin_size) * bin_size]
prof_dist <- da[, .(n_edited = sum(is_edited), n_total = .N), by = .(len_band, x)]
prof_dist[, scale := "dist3p"]

# ---- fractional 5'->3' metagene, per length band ----
d[, x := pmin(as.integer(frac_5to3 * n_frac_bins), n_frac_bins - 1L)]
prof_frac <- d[, .(n_edited = sum(is_edited), n_total = .N), by = .(len_band, x)]
prof_frac[, `:=`(x = x / n_frac_bins, scale = "frac")]

out <- rbind(prof_dist[, .(len_band, scale, x = as.numeric(x), n_edited, n_total)],
             prof_frac[, .(len_band, scale, x = as.numeric(x), n_edited, n_total)])
fwrite(out, out_path, sep = "\t")
cat("metagene:", nrow(d), "genic obs ->", out_path, "\n")
