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
--max-donor / --max-acceptor: splice sites are forbidden: every donor (5' splice site) and acceptor (3' splice site) window inside the CDS and stop codon must score
below these MaxEntScan scores (maxentscan.py, a port of the Perl scripts). They are enforced together with the codon optimisation: DnaChisel first resolves the constraints, then optimises
the codon usage without breaking them.
--splice-max/--splice-step: if the splice limits cannot be met together with the other constraints, both limits are raised by --splice-step until --splice-max.
--splice-fallback: a protein that still fails is designed WITHOUT splice limits (all other constraints kept).
splice_design = limit used (e.g. 7, 8.5), 'unconstrained' (fallback) or 'none' (no limits given).
A record that cannot be designed (letters other than the 20 amino acids, empty) does not stop the
run: it gets an empty cds, FALSE flags and the reason in the column error.

Usage:
  python3 reverse_translate.py proteins.faa --output cds.tsv --avoid XhoI,XbaI --threads 16
Output TSV: id, aa_len, cds (with the stop codon), gc, max_homopolymer, translation_ok, sites, constraints_ok, max_donor, max_acceptor, splice_design, error
"""

import argparse
import os
import random
import re
import sys
import zlib
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import maxentscan

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


def make_splice_spec():
    from dnachisel import Specification, SpecEvaluation, Location

    class AvoidSpliceSites(Specification):
        """No donor / acceptor window inside the location scores at or above max_donor / max_acceptor."""
        best_possible_score = 0

        def __init__(self, max_donor=None, max_acceptor=None, location=None, boost=1.0):
            self.max_donor, self.max_acceptor = max_donor, max_acceptor
            self.location = Location.from_data(location)
            self.boost = boost

        def initialized_on_problem(self, problem, role=None):
            return self._copy_with_full_span_if_no_location(problem)

        def evaluate(self, problem):
            loc = self.location
            don, acc = maxentscan.sites(problem.sequence, loc.start, loc.end)
            bad = []
            if self.max_donor is not None:
                bad += [Location(p - 2, p + 7) for p, sc in don if sc >= self.max_donor]
            if self.max_acceptor is not None:
                bad += [Location(p - 19, p + 4) for p, sc in acc if sc >= self.max_acceptor]
            return SpecEvaluation(self, problem, score=-len(bad), locations=bad or [loc],
                                  message='%d splice sites at or above the limits' % len(bad) if bad else 'no strong splice sites')

        def localized(self, location, problem=None, with_righthand=True):
            # a mutation inside `location` can only change the windows that overlap it (acceptor windows are the longest: 23 nt)
            lo = max(self.location.start, location.start - 22)
            hi = min(self.location.end, location.end + 22)
            if hi <= lo:
                return None
            return self.__class__(self.max_donor, self.max_acceptor, location=Location(lo, hi), boost=self.boost)

        def label_parameters(self):
            return ['max_donor=%s' % self.max_donor, 'max_acceptor=%s' % self.max_acceptor]

    return AvoidSpliceSites


def failed(msg):
    return {'cds': '', 'gc': float('nan'), 'max_homopolymer': 0, 'translation_ok': False, 'sites': 0,
            'constraints_ok': False, 'max_donor': float('nan'), 'max_acceptor': float('nan'), 'splice_design': 'none', 'error': msg}


def design(protein):
    # splice limits start at --max-donor/--max-acceptor and rise by --splice-step up to --splice-max; if none can be met the protein is designed without them
    if ARGS.max_donor is None and ARGS.max_acceptor is None:
        res = attempt(protein, None)
        res['splice_design'] = 'none'
        return protein, res
    d0, a0 = ARGS.max_donor, ARGS.max_acceptor
    top = ARGS.splice_max if ARGS.splice_max is not None else 0
    n_steps = int(round(max(0, top - min(x for x in (d0, a0) if x is not None)) / ARGS.splice_step))
    for k in range(n_steps + 1):
        delta = k * ARGS.splice_step
        res = attempt(protein, (None if d0 is None else d0 + delta, None if a0 is None else a0 + delta))
        res['splice_design'] = '%g' % (min(x for x in (d0, a0) if x is not None) + delta)
        if res['constraints_ok'] or res['cds'] == '':
            return protein, res
    if ARGS.splice_fallback:
        res = attempt(protein, None)
        res['splice_design'] = 'unconstrained'
    return protein, res


def attempt(protein, splice):
    if not protein:
        return failed('empty')
    bad = sorted(set(protein) - AA)
    if bad:
        return failed('invalid residues ' + ''.join(bad))
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
    if splice:
        constraints.append(make_splice_spec()(max_donor=splice[0], max_acceptor=splice[1], location=Location(a, b + len(stop))))
    errors = []
    try:
        problem = DnaOptimizationProblem(sequence=pre + reverse_translate(protein) + suf, constraints=constraints,
                                         objectives=[CodonOptimize(species=ARGS.species, method='match_codon_usage', location=cds_loc)],
                                         logger=None)
    except Exception as e:
        return failed('problem not built: %s: %s' % (type(e).__name__, e))
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
    don, acc = maxentscan.sites(cds)
    max_d = max((sc for _, sc in don), default=float('nan'))
    max_a = max((sc for _, sc in acc), default=float('nan'))
    return {'cds': cds, 'gc': (cds.count('G') + cds.count('C')) / len(cds), 'max_homopolymer': run,
                     'translation_ok': str(Seq(dna).translate()) == protein, 'sites': sites,
                     'constraints_ok': bool(problem.all_constraints_pass()), 'max_donor': max_d, 'max_acceptor': max_a, 'error': ';'.join(errors)}


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
    ap.add_argument('--splice-max', type=float, default=None)
    ap.add_argument('--splice-step', type=float, default=0.5)
    ap.add_argument('--splice-fallback', action='store_true')
    ap.add_argument('--max-donor', type=float, default=None)
    ap.add_argument('--max-acceptor', type=float, default=None)
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
        out.write('id\taa_len\tcds\tgc\tmax_homopolymer\ttranslation_ok\tsites\tconstraints_ok\tmax_donor\tmax_acceptor\tsplice_design\terror\n')
        for i, p in zip(ids, prots):
            r = res[p]
            out.write('%s\t%d\t%s\t%.4f\t%d\t%s\t%d\t%s\t%.2f\t%.2f\t%s\t%s\n' % (i, len(p), r['cds'], r['gc'], r['max_homopolymer'], r['translation_ok'], r['sites'], r['constraints_ok'], r['max_donor'], r['max_acceptor'], r['splice_design'], r['error']))
    os.replace(args.output + '.part', args.output)


if __name__ == '__main__':
    main()
