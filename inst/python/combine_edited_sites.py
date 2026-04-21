#!/usr/bin/env python3
"""
Combine per-position edited site BEDs from R1 and R2 passes.
Sums ref_count and alt_count for matching (chrom, pos, strand, name).

Usage:
  python3 combine_edited_sites.py r1.bed r2.bed --output combined.bed
"""

import argparse
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inputs', nargs='+')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    counts = defaultdict(lambda: [0, 0])

    for path in args.inputs:
        with open(path) as fh:
            next(fh)
            for line in fh:
                p = line.rstrip('\n').split('\t')
                # BED columns: chrom, start, end, name, alt_count, strand, ref_count
                key = (p[0], int(p[1]), p[3], p[5])  # chrom, start, name, strand
                counts[key][0] += int(p[6])  # ref_count
                counts[key][1] += int(p[4])  # alt_count

    with open(args.output, 'w') as fh:
        fh.write('chrom\tstart\tend\tname\talt_count\tstrand\tref_count\n')
        for (chrom, start, name, strand), (ref_c, alt_c) in sorted(
                counts.items(), key=lambda x: (x[0][0], x[0][1])):
            fh.write(
                f'{chrom}\t{start}\t{start+1}\t{name}'
                f'\t{alt_c}\t{strand}\t{ref_c}\n'
            )


if __name__ == '__main__':
    main()
