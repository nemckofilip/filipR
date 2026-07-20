#!/usr/bin/env Rscript
# Called by fn_annotate_edit_calls(). Do not run directly.
# Maps per-observation edit calls to their gene's mature-transcript coordinate.
# Args: <edit_calls.tsv.gz> <gene_model.rds> <output.rds>
#
# gene_model.rds: exon-union segments (built once by Annotate_edit_calls.R) with
# columns chrom, start, end, strand, gene_id, gene_type, gene_len, cum_start
# and keyed on (chrom, start, end) for foverlaps.
#
# Output RDS: the edit calls joined to gene coordinates (genic observations only):
#   read_id, chrom, pos, strand, is_edited, gene_id, gene_type, gene_len,
#   mature_pos (0-based, 5'->3'), frac_5to3, dist_from_3p

suppressMessages(library(data.table))

args       <- commandArgs(trailingOnly = TRUE)
calls_path <- args[1]
model_path <- args[2]
out_path   <- args[3]

model <- readRDS(model_path)
setkey(model, chrom, start, end)

dt <- fread(calls_path)   # read_id chrom pos strand is_edited

# Annotate the UNIQUE positions (far fewer than observations), then join back.
u <- unique(dt[, .(chrom, pos, strand)])
u[, `:=`(start = pos, end = pos + 1L)]

ov <- foverlaps(u, model, by.x = c("chrom", "start", "end"),
                type = "any", nomatch = 0L)
ov <- ov[i.strand == strand]                       # sense-strand match to gene
ov[, mature_pos := cum_start + (i.start - start)]  # 5'->3' offset (genomic order)
ov[strand == "-", mature_pos := gene_len - 1L - mature_pos]
# One gene per position (drop rare ambiguous multi-gene overlaps).
ov <- unique(ov, by = c("chrom", "pos", "strand"))

ann <- ov[, .(chrom, pos, strand, gene_id, gene_type, gene_len, mature_pos,
              frac_5to3    = mature_pos / gene_len,
              dist_from_3p = gene_len - 1L - mature_pos)]

res <- merge(dt, ann, by = c("chrom", "pos", "strand"))   # genic observations
saveRDS(res, out_path)
cat("annotated", nrow(res), "genic observations of", nrow(dt), "->", out_path, "\n")
