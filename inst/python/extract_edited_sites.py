#!/usr/bin/env python3
"""
Per-position edited sites from samtools mpileup (stdin).

Outputs one row per (chrom, pos, strand, ref, alt) where alt_count > 0
for the requested edit types. Strand and ref/alt are at the transcript level.

Usage (R1 pass):
  samtools mpileup --rf 64 ... bam | python3 extract_edited_sites.py \
      --edit-types A>G,C>T --output r1.bed

Usage (R2 pass, forward-stranded PE):
  samtools mpileup --rf 128 ... bam | python3 extract_edited_sites.py \
      --edit-types A>G,C>T --swap-strands --output r2.bed
"""

import sys, argparse
from collections import Counter

BASES = 'ACGT'
_CT   = str.maketrans('ACGTacgt', 'TGCAtgca')

def _c(b): return b.translate(_CT)

def parse_pileup(s):
    fwd_mm, rev_mm = [], []
    fwd_n = rev_n = 0
    i = 0
    while i < len(s):
        c = s[i]
        if   c == '^':              i += 2
        elif c == '$':              i += 1
        elif c in '+-':
            i += 1; j = i
            while j < len(s) and s[j].isdigit(): j += 1
            i = j + int(s[i:j])
        elif c == '*':              i += 1
        elif c == '.':   fwd_n += 1; i += 1
        elif c == ',':   rev_n += 1; i += 1
        elif c in BASES: fwd_n += 1; fwd_mm.append(c);        i += 1
        elif c.upper() in BASES:
                         rev_n += 1; rev_mm.append(c.upper()); i += 1
        else:            i += 1
    return fwd_mm, fwd_n, rev_mm, rev_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output',       required=True)
    ap.add_argument('--edit-types',   required=True,
                    help='Comma-separated, e.g. A>G,C>T')
    ap.add_argument('--swap-strands', action='store_true',
                    help='Use for R2 in forward-stranded PE libraries')
    args = ap.parse_args()

    et_pairs = [(t[0], t[2]) for t in args.edit_types.split(',')]

    with open(args.output, 'w') as out:
        out.write('chrom\tstart\tend\tname\talt_count\tstrand\tref_count\n')

        for line in sys.stdin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            chrom = parts[0]
            pos   = int(parts[1]) - 1  # 0-based
            ref   = parts[2].upper()
            if ref not in BASES:
                continue
            cr = _c(ref)

            fwd_list, fwd_n, rev_list, rev_n = parse_pileup(parts[4])

            if not args.swap_strands:
                # fwd reads → + strand (transcript ref = ref)
                # rev reads → - strand (transcript ref = cr, mismatch complemented)
                fwd_cnt = Counter(fwd_list)
                rev_cnt = Counter(_c(b) for b in rev_list)

                for r, a in et_pairs:
                    if ref == r and fwd_n > 0:
                        alt_c = fwd_cnt.get(a, 0)
                        if alt_c > 0:
                            ref_c = fwd_n - sum(fwd_cnt.values())
                            out.write(
                                f'{chrom}\t{pos}\t{pos+1}\t{r}>{a}'
                                f'\t{alt_c}\t+\t{ref_c}\n'
                            )
                    if cr == r and rev_n > 0:
                        alt_c = rev_cnt.get(a, 0)
                        if alt_c > 0:
                            ref_c = rev_n - sum(rev_cnt.values())
                            out.write(
                                f'{chrom}\t{pos}\t{pos+1}\t{r}>{a}'
                                f'\t{alt_c}\t-\t{ref_c}\n'
                            )
            else:
                # R2 swap: R2 fwd → effective - strand; R2 rev → effective + strand
                eff_fwd_cnt = Counter(rev_list)                 # R2 rev → + strand
                eff_rev_cnt = Counter(_c(b) for b in fwd_list) # R2 fwd → - strand

                for r, a in et_pairs:
                    if ref == r and rev_n > 0:
                        alt_c = eff_fwd_cnt.get(a, 0)
                        if alt_c > 0:
                            ref_c = rev_n - sum(eff_fwd_cnt.values())
                            out.write(
                                f'{chrom}\t{pos}\t{pos+1}\t{r}>{a}'
                                f'\t{alt_c}\t+\t{ref_c}\n'
                            )
                    if cr == r and fwd_n > 0:
                        alt_c = eff_rev_cnt.get(a, 0)
                        if alt_c > 0:
                            ref_c = fwd_n - sum(eff_rev_cnt.values())
                            out.write(
                                f'{chrom}\t{pos}\t{pos+1}\t{r}>{a}'
                                f'\t{alt_c}\t-\t{ref_c}\n'
                            )


if __name__ == '__main__':
    main()
