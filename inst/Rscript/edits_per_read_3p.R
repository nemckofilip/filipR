#!/usr/bin/env Rscript
# Called by fn_edits_per_read_3p(). Do not run directly.
# Per-read on-target edit count from annotated edit calls, split into all genic
# reads vs reads overlapping the 3'-proximal window (a read is "proximal" if its
# closest target base to the 3' end is within dist_thr nt).
#
# Args: <edit_calls_annot.rds> <out.tsv> <dist_thr>
# Output TSV columns: region ("all" | "prox"), n_edited, n_reads.

suppressMessages(library(data.table))

a        <- commandArgs(trailingOnly = TRUE)
in_path  <- a[1]
out_path <- a[2]
dist_thr <- as.integer(a[3])

d <- readRDS(in_path)

# per read: on-target edit count and closest approach to the 3' end
reads <- d[, .(n_edited = sum(is_edited), min_d = min(dist_from_3p)),
           by = read_id]

all_d  <- reads[,                .(n_reads = .N), by = n_edited][, region := "all"]
prox_d <- reads[min_d < dist_thr, .(n_reads = .N), by = n_edited][, region := "prox"]

out <- rbind(all_d, prox_d)[order(region, n_edited)]
fwrite(out, out_path, sep = "\t")
cat("edits/read 3p:", nrow(reads), "reads,",
    reads[min_d < dist_thr, .N], "proximal ->", out_path, "\n")
