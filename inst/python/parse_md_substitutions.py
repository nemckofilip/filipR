#!/usr/bin/env python3
"""
Reads samtools view output from stdin, counts base substitutions via MD tag + CIGAR.
Splits counts into human chromosomes (chr*) vs ERCC spike-ins (ERCC-*).

Output files (all tab-separated):
  --out-human / --out-ercc                    REF>ALT\tcount
  --out-human-pos / --out-ercc-pos            REF>ALT\tread_pos\tcount
  --out-human-base / --out-ercc-base          ref_b\tcount
  --out-human-base-pos / --out-ercc-base-pos  ref_b\tread_pos\tcount

With --rf-stranded (e.g. Lexogen CORALL, PE only):
  Substitutions and base counts are reported from the mRNA perspective.
  CORALL is FR stranded: R1 is sense, R2 is antisense to the mRNA.
  Reads from minus-strand genes are complement-flipped (flip condition:
  is_r1 == is_reverse) so T>C becomes A>G, G>A becomes C>T, etc.
  Positions use R2 identity for the 5-prime-of-mRNA end (R2 = antisense,
  so R2 position 0 is closest to the 5' end of the mRNA fragment).
"""
import sys
import re
import argparse
from collections import defaultdict

COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def aligned_bases(seq, cigar):
    """Return read bases at reference-consuming positions (M/=/X), in order."""
    result = []
    pos = 0
    for m in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        n, op = int(m.group(1)), m.group(2)
        if op in ('M', '=', 'X'):
            result.extend(seq[pos:pos + n])
            pos += n
        elif op in ('I', 'S'):
            pos += n
    return result


def write_counts(counts, path):
    with open(path, 'w') as fh:
        for key, count in sorted(counts.items()):
            fh.write(f"{key}\t{count}\n")


def write_pos_counts(counts, path):
    with open(path, 'w') as fh:
        for (key, pos), count in sorted(counts.items()):
            fh.write(f"{key}\t{pos}\t{count}\n")


parser = argparse.ArgumentParser()
parser.add_argument('--out-human',          required=True)
parser.add_argument('--out-ercc',           required=True)
parser.add_argument('--out-human-pos',      default=None)
parser.add_argument('--out-ercc-pos',       default=None)
parser.add_argument('--out-human-base',     default=None)
parser.add_argument('--out-ercc-base',      default=None)
parser.add_argument('--out-human-base-pos', default=None)
parser.add_argument('--out-ercc-base-pos',  default=None)
parser.add_argument('--rf-stranded',        action='store_true', default=False)
args = parser.parse_args()

counts_human          = defaultdict(int)
counts_ercc           = defaultdict(int)
counts_human_pos      = defaultdict(int)
counts_ercc_pos       = defaultdict(int)
counts_human_base     = defaultdict(int)
counts_ercc_base      = defaultdict(int)
counts_human_base_pos = defaultdict(int)
counts_ercc_base_pos  = defaultdict(int)

for line in sys.stdin:
    if line.startswith('@'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 10:
        continue
    flag = int(f[1])
    if flag & (4 | 256 | 2048):   # unmapped, secondary, supplementary
        continue

    rname      = f[2]
    is_ercc    = rname.startswith('ERCC-')
    is_reverse = bool(flag & 0x10)
    is_r1      = bool(flag & 0x40)

    if args.rf_stranded:
        # flip when gene is on minus strand: for FR (CORALL), R1 maps sense,
        # so minus-strand reads have is_r1 == is_reverse
        flip_sub = (is_r1 == is_reverse)
        # position 0 = 5' end of mRNA fragment: R2 is antisense in FR
        flip_pos = not is_r1
    else:
        flip_sub = False
        flip_pos = is_reverse

    cigar, seq = f[5], f[9]
    md = next((t[5:] for t in f[11:] if t.startswith('MD:Z:')), None)
    if not md:
        continue

    aln     = aligned_bases(seq, cigar)
    aln_len = len(aln)

    ref_pos, i = 0, 0
    while i < len(md):
        c = md[i]
        if c.isdigit():
            j = i
            while j < len(md) and md[j].isdigit():
                j += 1
            n = int(md[i:j])
            # match region: read base == ref base at each position
            for k in range(ref_pos, min(ref_pos + n, aln_len)):
                rb   = aln[k].upper()
                if flip_sub:
                    rb = COMP.get(rb, rb)
                spos = (aln_len - 1 - k) if flip_pos else k
                if rb in 'ACGT':
                    if is_ercc:
                        counts_ercc_base[rb] += 1
                        counts_ercc_base_pos[(rb, spos)] += 1
                    else:
                        counts_human_base[rb] += 1
                        counts_human_base_pos[(rb, spos)] += 1
            ref_pos += n
            i = j
        elif c == '^':
            i += 1
            while i < len(md) and md[i].isalpha():
                i += 1
        elif c.isalpha():
            ref_b  = c.upper()
            spos   = (aln_len - 1 - ref_pos) if flip_pos else ref_pos
            if ref_pos < aln_len:
                read_b = aln[ref_pos].upper()
                if flip_sub:
                    ref_b_out  = COMP.get(ref_b,  ref_b)
                    read_b_out = COMP.get(read_b, read_b)
                else:
                    ref_b_out  = ref_b
                    read_b_out = read_b
                if ref_b_out in 'ACGT':
                    if is_ercc:
                        counts_ercc_base[ref_b_out] += 1
                        counts_ercc_base_pos[(ref_b_out, spos)] += 1
                    else:
                        counts_human_base[ref_b_out] += 1
                        counts_human_base_pos[(ref_b_out, spos)] += 1
                if (ref_b_out in 'ACGT' and read_b_out in 'ACGT'
                        and ref_b_out != read_b_out):
                    sub = f"{ref_b_out}>{read_b_out}"
                    if is_ercc:
                        counts_ercc[sub] += 1
                        counts_ercc_pos[(sub, spos)] += 1
                    else:
                        counts_human[sub] += 1
                        counts_human_pos[(sub, spos)] += 1
            ref_pos += 1
            i += 1
        else:
            i += 1

write_counts(counts_human, args.out_human)
write_counts(counts_ercc,  args.out_ercc)

if args.out_human_pos:
    write_pos_counts(counts_human_pos, args.out_human_pos)
if args.out_ercc_pos:
    write_pos_counts(counts_ercc_pos, args.out_ercc_pos)
if args.out_human_base:
    write_counts(counts_human_base, args.out_human_base)
if args.out_ercc_base:
    write_counts(counts_ercc_base, args.out_ercc_base)
if args.out_human_base_pos:
    write_pos_counts(counts_human_base_pos, args.out_human_base_pos)
if args.out_ercc_base_pos:
    write_pos_counts(counts_ercc_base_pos, args.out_ercc_base_pos)
