#!/usr/bin/env python3
"""
Per-base editing index from samtools mpileup (stdin).

Strand-correct for all 12 substitution types. At each pileup position (ref = R):
  - Forward reads  -> edit types  R -> Y           (+ strand transcripts)
  - Reverse reads  -> edit types  comp(R) -> comp(Y)  (- strand transcripts)
    A reverse read showing pileup base B means the - strand RNA has comp(B)
    instead of comp(R), so the RNA-level edit is comp(R) -> comp(B).

ERCC chromosomes (ERCC-*) are tallied separately from human.

Usage:
  samtools mpileup -Q 0 -q 0 --no-BAQ -d 1000000 -f ref.fa sample.bam | \
    python3 compute_editing_index.py --sample ID --output out.tsv
"""

import sys, re, argparse
from collections import defaultdict

BASES = 'ACGT'
_CT   = str.maketrans('ACGTacgt', 'TGCAtgca')

def _c(b): return b.translate(_CT)

EDIT_TYPES = [f'{x}>{y}' for x in BASES for y in BASES if x != y]


def parse_pileup(s):
    """Return (fwd_mismatches, fwd_n, rev_mismatches, rev_n).
    Mismatches are the observed base (uppercase). Rev mismatches are + strand equivalent."""
    fwd_mm, rev_mm = [], []
    fwd_n = rev_n = 0
    i = 0
    while i < len(s):
        c = s[i]
        if   c == '^':                          i += 2          # ^<mapq>
        elif c == '$':                          i += 1
        elif c in '+-':                                         # indel
            i += 1
            j  = i
            while j < len(s) and s[j].isdigit(): j += 1
            i  = j + int(s[i:j])
        elif c == '*':                          i += 1          # deletion
        elif c == '.':   fwd_n += 1;            i += 1          # fwd match
        elif c == ',':   rev_n += 1;            i += 1          # rev match
        elif c in BASES: fwd_n += 1; fwd_mm.append(c);   i += 1
        elif c.upper() in BASES:
            rev_n += 1; rev_mm.append(c.upper()); i += 1
        else:                                   i += 1
    return fwd_mm, fwd_n, rev_mm, rev_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample',  required=True)
    ap.add_argument('--output',  required=True)
    args = ap.parse_args()

    fwd_cov = {'human': defaultdict(int), 'ercc': defaultdict(int)}
    rev_cov = {'human': defaultdict(int), 'ercc': defaultdict(int)}
    fwd_mm  = {'human': defaultdict(lambda: defaultdict(int)),
               'ercc':  defaultdict(lambda: defaultdict(int))}
    rev_mm  = {'human': defaultdict(lambda: defaultdict(int)),
               'ercc':  defaultdict(lambda: defaultdict(int))}

    for line in sys.stdin:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 5:
            continue
        chrom, ref, pileup = parts[0], parts[2].upper(), parts[4]
        if ref not in BASES:
            continue

        region = 'ercc' if chrom.startswith('ERCC-') else 'human'
        cr     = _c(ref)

        fwd_list, fwd_n, rev_list, rev_n = parse_pileup(pileup)

        fwd_cov[region][ref] += fwd_n
        for b in fwd_list:
            fwd_mm[region][ref][b] += 1

        rev_cov[region][cr] += rev_n
        for b in rev_list:
            rev_mm[region][cr][_c(b)] += 1

    with open(args.output, 'w') as f:
        f.write('sample\tregion\tedit_type\tnumerator\tdenominator\tediting_index\n')
        for region in ('human', 'ercc'):
            for et in EDIT_TYPES:
                x, y  = et[0], et[2]
                num   = fwd_mm[region][x][y] + rev_mm[region][x][y]
                denom = fwd_cov[region][x]   + rev_cov[region][x]
                ei    = num / denom if denom > 0 else float('nan')
                f.write(f'{args.sample}\t{region}\t{et}\t{num}\t{denom}\t{ei:.8g}\n')


if __name__ == '__main__':
    main()
