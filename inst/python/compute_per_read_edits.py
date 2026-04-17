#!/usr/bin/env python3
"""
Per-read editing distribution from samtools view (stdin).

RF-stranded library (dUTP): R1 antisense, R2 sense.
Flip condition: is_r1 != is_reverse.
  - R2 forward / R1 reverse  (+ strand gene reads) : no flip -> ref->read as-is
  - R2 reverse / R1 forward  (- strand gene reads) : flip    -> comp(ref)->comp(read)
This puts all mismatches in the mRNA/sense-strand coordinate frame.

ERCC chromosomes (ERCC-*) are tallied separately.

Output: one row per sample x region x edit_type with:
  n_reads        total reads (all, not just edited)
  mean_edits     mean edits of this type per read (over all reads)
  frac_ge1       fraction of reads with >= 1 edit
  frac_ge3       fraction of reads with >= 3 edits
  frac_hyperedit   fraction of reads with >= 5 edits AND edits/aligned_length >= 0.05

Usage:
  samtools view -F 2308 sample.bam | \
    python3 compute_per_read_edits.py --sample ID --output out.tsv
"""

import sys, re, argparse
from collections import defaultdict

BASES = 'ACGT'
_CT   = str.maketrans('ACGTacgt', 'TGCAtgca')

def _c(b): return b.upper().translate(_CT)

EDIT_TYPES = [f'{x}>{y}' for x in BASES for y in BASES if x != y]


def aligned_bases(seq, cigar):
    """Read bases at each reference-consuming (M/=/X) position, in order."""
    out, pos = [], 0
    for m in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        n, op = int(m.group(1)), m.group(2)
        if op in 'M=X':
            out.extend(seq[pos:pos + n])
            pos += n
        elif op in 'IS':
            pos += n
    return out


def parse_md(md, aln):
    """Yield (ref_base, read_base) for each mismatch."""
    ref_pos = i = 0
    while i < len(md):
        c = md[i]
        if c.isdigit():
            j = i
            while j < len(md) and md[j].isdigit(): j += 1
            ref_pos += int(md[i:j]); i = j
        elif c == '^':
            i += 1
            while i < len(md) and md[i].isalpha(): i += 1
        elif c.isalpha():
            if ref_pos < len(aln):
                yield c.upper(), aln[ref_pos].upper()
            ref_pos += 1; i += 1
        else:
            i += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample',  required=True)
    ap.add_argument('--output',  required=True)
    args = ap.parse_args()

    # counts[region][edit_type] = [total_edits, n_reads, ge1, ge3, hyperedit]
    counts = {r: {et: [0, 0, 0, 0, 0] for et in EDIT_TYPES}
              for r in ('human', 'ercc')}

    for line in sys.stdin:
        if line.startswith('@'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 10: continue

        flag  = int(f[1])
        chrom = f[2]
        cigar = f[5]
        seq   = f[9]
        if cigar == '*' or not seq: continue

        md = next((t[5:] for t in f[11:] if t.startswith('MD:Z:')), None)
        if md is None: continue

        region     = 'ercc' if chrom.startswith('ERCC-') else 'human'
        is_reverse = bool(flag & 0x10)
        is_r1      = bool(flag & 0x40)
        flip       = (is_r1 != is_reverse)   # RF library: flip = - strand gene read

        aln       = aligned_bases(seq, cigar)
        edit_cnts = defaultdict(int)

        for ref_b, read_b in parse_md(md, aln):
            if ref_b not in BASES or read_b not in BASES or ref_b == read_b:
                continue
            if flip:
                ref_b, read_b = _c(ref_b), _c(read_b)
            edit_cnts[f'{ref_b}>{read_b}'] += 1

        for et in EDIT_TYPES:
            n  = edit_cnts[et]
            c  = counts[region][et]
            c[0] += n            # total edits
            c[1] += 1            # n_reads
            if n >= 1:  c[2] += 1
            if n >= 3:  c[3] += 1
            if n >= 5 and len(aln) > 0 and n / len(aln) >= 0.05: c[4] += 1

    with open(args.output, 'w') as f:
        f.write('sample\tregion\tedit_type\tn_reads\tmean_edits\tfrac_ge1\tfrac_ge3\tfrac_hyperedit\n')
        for region in ('human', 'ercc'):
            for et in EDIT_TYPES:
                tot, n, ge1, ge3, ge10 = counts[region][et]
                if n == 0:
                    f.write(f'{args.sample}\t{region}\t{et}\t0\tNA\tNA\tNA\tNA\n')
                else:
                    f.write(f'{args.sample}\t{region}\t{et}\t{n}\t'
                            f'{tot/n:.6g}\t{ge1/n:.6g}\t{ge3/n:.6g}\t{ge10/n:.6g}\n')


if __name__ == '__main__':
    main()
