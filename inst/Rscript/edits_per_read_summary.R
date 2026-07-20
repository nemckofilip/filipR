#!/usr/bin/env Rscript
# Called by fn_edits_per_read(). Do not run directly.
# From a per-observation edit_calls table, compute the per-read on-target
# edit-count distribution (how many reads carry 0,1,2,... edits).
# Args: <edit_calls.tsv.gz> <output.rds>
#
# Output RDS: data.table(n_edited, n_reads); total read count in attr "total_reads".

suppressMessages(library(data.table))

a  <- commandArgs(trailingOnly = TRUE)
d  <- fread(a[1], select = c("read_id", "is_edited"))
pr <- d[, .(n_edited = sum(is_edited)), by = read_id]
dist <- pr[, .(n_reads = .N), by = n_edited][order(n_edited)]
attr(dist, "total_reads") <- nrow(pr)
saveRDS(dist, a[2])
cat("reads:", nrow(pr), " max edits/read:", max(pr$n_edited), "->", a[2], "\n")
