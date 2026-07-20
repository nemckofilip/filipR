#!/usr/bin/env python3
"""Per-observation edit calls from a BAM (called by fn_edit_calls()).

Walks every read once and emits ONE row per reference target base it covers,
keeping the read identity so per-read and per-site views stay linked:

    read_id  chrom  pos  strand  is_edited

  - read_id : integer per fragment (both mates of a pair share it)
  - pos     : 0-based genomic position of the reference target base
  - strand  : transcript strand ('+' if R1 XOR reverse, forward-stranded PE)
  - is_edited: 1 if the read base is the (transcript) edited base, else 0

The target base is defined by --edit-type at the transcript level (e.g. A>G):
on + strand transcripts the genomic ref/alt are A/G, on - strand they are the
complements T/C. Reference bases come from the genome FASTA (BAM read
sequences are original, not 3N-converted). Only the two on-target alleles are
recorded; other mismatches are ignored. Output is gzip TSV.
"""
import argparse
import gzip
import pysam

COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def target_bases(edit_type, strand):
    x, y = edit_type.split(">")
    return (x, y) if strand == "+" else (COMP[x], COMP[y])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--edit-type", required=True)   # e.g. "A>G"
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-mq", type=int, default=20)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fa = pysam.FastaFile(args.fasta)

    rid = {}
    next_id = 0
    min_bq = args.min_bq
    min_mq = args.min_mq

    with gzip.open(args.output, "wt") as out:
        out.write("read_id\tchrom\tpos\tstrand\tis_edited\n")
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mq:
                continue
            rend = read.reference_end
            qseq = read.query_sequence
            if rend is None or qseq is None:
                continue

            strand = "+" if (read.is_read1 != read.is_reverse) else "-"
            refb, altb = target_bases(args.edit_type, strand)

            chrom = read.reference_name
            rstart = read.reference_start
            refseq = fa.fetch(chrom, rstart, rend).upper()
            quals = read.query_qualities

            qn = read.query_name
            myid = rid.get(qn)
            if myid is None:
                myid = next_id
                rid[qn] = next_id
                next_id += 1

            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if refseq[rpos - rstart] != refb:
                    continue
                if quals is not None and quals[qpos] < min_bq:
                    continue
                qb = qseq[qpos].upper()
                if qb == refb:
                    ed = 0
                elif qb == altb:
                    ed = 1
                else:
                    continue
                out.write(f"{myid}\t{chrom}\t{rpos}\t{strand}\t{ed}\n")


if __name__ == "__main__":
    main()
