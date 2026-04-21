#!/usr/bin/env python3
"""
Sum numerators and denominators across two or more editing_index TSV files.
Recomputes editing_index = numerator / denominator for each region x edit_type.

Usage:
  python3 combine_editing_index.py r1.tsv r2.tsv --output combined.tsv
"""

import argparse
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inputs', nargs='+')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    counts = defaultdict(lambda: [0, 0])
    sample = None

    for path in args.inputs:
        with open(path) as fh:
            next(fh)
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if sample is None:
                    sample = p[0]
                key = (p[1], p[2])
                counts[key][0] += int(p[3])
                counts[key][1] += int(p[4])

    with open(args.output, 'w') as fh:
        fh.write('sample\tregion\tedit_type\tnumerator\tdenominator\tediting_index\n')
        for (region, et), (num, denom) in sorted(counts.items()):
            ei = num / denom if denom > 0 else float('nan')
            fh.write(f'{sample}\t{region}\t{et}\t{num}\t{denom}\t{ei:.8g}\n')


if __name__ == '__main__':
    main()
