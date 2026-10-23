#!/usr/bin/env python3
"""
Reverse-translate proteins into codon-optimised coding sequences with DnaChisel.

Every distinct protein is processed once (identical proteins get identical DNA). The result is
deterministic: the random seed is a hash of the protein sequence. Constraints: exact translation,
GC content window, no recognition site of the listed restriction enzymes, no homopolymer of
--homopolymer or more bases, no AATAAA. Objective: codon usage matching the species table.
The CDS starts with the codon of the first residue (M = ATG) and ends with the stop codon --stop (default TAA).
--prefix / --suffix: constant sequences placed before the CDS / after the stop codon; they are not changed, but the
constraints (restriction sites, homopolymers, AATAAA) are checked on the whole sequence, so sites
that span a junction are avoided too. The output CDS does not contain them.
A record that cannot be designed (letters other than the 20 amino acids, empty) does not stop the
run: it gets an empty cds, FALSE flags and the reason in the column error.

Usage:
  python3 reverse_translate.py proteins.faa --output cds.tsv --avoid XhoI,XbaI --threads 16
Output TSV: id, aa_len, cds (with the stop codon), gc, max_homopolymer, translation_ok, sites, constraints_ok, error
"""

import argparse
import os
import random
import re
import sys
import zlib
from multiprocessing import Pool

import numpy as np

AA = set('ACDEFGHIKLMNPQRSTVWY')
ARGS = None
ENZYMES = []


def init(args, enzymes):
    global ARGS, ENZYMES
    ARGS, ENZYMES = args, enzymes


def read_fasta(path):
    ids, seqs = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                ids.append(line[1:].split()[0] if line[1:].split() else '')
                seqs.append([])
            elif line and seqs:
                seqs[-1].append(line)
    return ids, [''.join(s).upper().rstrip('*') for s in seqs]


def site_of(enzyme):
    from Bio.Restriction import Restriction_Dictionary as rd
    return rd.rest_dict[enzyme]['site']


def failed(msg):
    return {'cds': '', 'gc': float('nan'), 'max_homopolymer': 0, 'translation_ok': False, 'sites': 0,
            'constraints_ok': False, 'error': msg}


def design(protein):
    if not protein:
        return protein, failed('empty')
    bad = sorted(set(protein) - AA)
    if bad:
        return protein, failed('invalid residues ' + ''.join(bad))
    from Bio.Restriction import RestrictionBatch
    from Bio.Seq import Seq
    from dnachisel import (AvoidPattern, CodonOptimize, DnaOptimizationProblem, AvoidChanges,
                           EnforceGCContent, EnforceTranslation, Location, reverse_translate)
    seed = zlib.crc32(protein.encode())
    random.seed(seed)
    np.random.seed(seed)
    pre, stop, suf = ARGS.prefix.upper(), ARGS.stop.upper(), ARGS.suffix.upper()
    suf = stop + suf                                                   # the stop codon is constant: it joins the suffix
    a, b = len(pre), len(pre) + 3 * len(protein)
    cds_loc = Location(a, b)
    constraints = [EnforceTranslation(location=cds_loc),
                   EnforceGCContent(mini=ARGS.gc_min, maxi=ARGS.gc_max, window=ARGS.gc_window, location=cds_loc)]
    if pre:
        constraints.append(AvoidChanges(location=Location(0, a)))
    if suf:
        constraints.append(AvoidChanges(location=Location(b, b + len(suf))))
    # each pattern is only avoided where it overlaps the CDS (the region widened by the pattern length - 1): sites that lie entirely in the constant flanks (e.g. cloning sites) are allowed
    def avoid(pattern, length):
        m = length - 1
        return AvoidPattern(pattern, location=Location(max(0, a - m), min(b + len(suf), b + m)))
    site_len = {e: len(site_of(e)) for e in ENZYMES}
    constraints += [avoid(e + '_site', site_len[e]) for e in ENZYMES]
    constraints += [avoid('%dx%s' % (ARGS.homopolymer, x), ARGS.homopolymer) for x in 'ACGT']
    constraints += [avoid('AATAAA', 6)]
    errors = []
    try:
        problem = DnaOptimizationProblem(sequence=pre + reverse_translate(protein) + suf, constraints=constraints,
                                         objectives=[CodonOptimize(species=ARGS.species, method='match_codon_usage', location=cds_loc)],
                                         logger=None)
    except Exception as e:
        return protein, failed('problem not built: %s: %s' % (type(e).__name__, e))
    try:
        problem.resolve_constraints()
    except Exception as e:
        errors.append('resolve_constraints: ' + type(e).__name__)
    try:                                                                # optimise even if some constraint could not be resolved
        problem.optimize()
    except Exception as e:
        errors.append('optimize: ' + type(e).__name__)
    full = problem.sequence
    dna = full[a:b]                                                     # the designed coding sequence without the stop codon
    n_sites = lambda x: sum(len(e.search(Seq(x))) for e in RestrictionBatch(ENZYMES)) if x else 0
    sites = n_sites(full) - n_sites(pre) - n_sites(suf)                  # sites that are not already in the constant flanks
    run = max((len(m.group(0)) for m in re.finditer(r'(A+|C+|G+|T+)', full)), default=0)
    cds = dna + stop
    return protein, {'cds': cds, 'gc': (cds.count('G') + cds.count('C')) / len(cds), 'max_homopolymer': run,
                     'translation_ok': str(Seq(dna).translate()) == protein, 'sites': sites,
                     'constraints_ok': bool(problem.all_constraints_pass()), 'error': ';'.join(errors)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fasta')
    ap.add_argument('--output', required=True)
    ap.add_argument('--avoid', default='')
    ap.add_argument('--species', default='h_sapiens')
    ap.add_argument('--gc-min', type=float, default=0.35)
    ap.add_argument('--gc-max', type=float, default=0.65)
    ap.add_argument('--gc-window', type=int, default=50)
    ap.add_argument('--homopolymer', type=int, default=6)
    ap.add_argument('--prefix', default='')
    ap.add_argument('--stop', default='TAA')
    ap.add_argument('--suffix', default='')
    ap.add_argument('--threads', type=int, default=1)
    args = ap.parse_args()
    enzymes = [e for e in args.avoid.split(',') if e]
    ids, prots = read_fasta(args.fasta)
    uniq = sorted(set(prots))
    print('%d sequences, %d distinct proteins' % (len(ids), len(uniq)), file=sys.stderr)
    with Pool(args.threads, initializer=init, initargs=(args, enzymes)) as pool:
        res = dict(pool.imap_unordered(design, uniq, chunksize=20))
    nfail = sum(1 for p in prots if res[p]['error'] or not res[p]['constraints_ok'])
    if nfail:
        print('WARNING: %d records with an error or an unmet constraint (see the column error)' % nfail, file=sys.stderr)
    with open(args.output + '.part', 'w') as out:
        out.write('id\taa_len\tcds\tgc\tmax_homopolymer\ttranslation_ok\tsites\tconstraints_ok\terror\n')
        for i, p in zip(ids, prots):
            r = res[p]
            out.write('%s\t%d\t%s\t%.4f\t%d\t%s\t%d\t%s\t%s\n' % (i, len(p), r['cds'], r['gc'], r['max_homopolymer'], r['translation_ok'], r['sites'], r['constraints_ok'], r['error']))
    os.replace(args.output + '.part', args.output)


if __name__ == '__main__':
    main()
