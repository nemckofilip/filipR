"""
Python port of MaxEntScan (Yeo & Burge 2004; the Perl scripts score5.pl and score3.pl of the bioconda
package 'maxentscan'), using the model tables of that package. Same numbers as the Perl scripts.

donor(s9)    score of a 9-mer: 3 exon + 6 intron nt, the donor GT at positions 4-5 (1-based)
acceptor(s23) score of a 23-mer: 20 intron + 3 exon nt, the acceptor AG at positions 19-20
Only 9-mers with GT / 23-mers with AG at these positions are splice sites; other windows still get a score.
"""

import glob
import math
import os
import sys

_DIR = None
_ME2X5 = _SEQ5 = _ACC = None
_BGD = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
_D1 = {'A': 0.004, 'C': 0.0032, 'G': 0.9896, 'T': 0.0032}
_D2 = {'A': 0.0034, 'C': 0.0039, 'G': 0.0042, 'T': 0.9884}
_A1 = {'A': 0.9903, 'C': 0.0032, 'G': 0.0034, 'T': 0.0030}
_A2 = {'A': 0.0027, 'C': 0.0037, 'G': 0.9905, 'T': 0.0030}
_FOUR = (1, 4, 16, 64, 256, 1024, 4096, 16384)
_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def _load():
    global _DIR, _ME2X5, _SEQ5, _ACC
    if _ME2X5 is not None:
        return
    hits = glob.glob(os.path.join(sys.prefix, 'share', 'maxentscan-*'))
    if not hits:
        raise RuntimeError('maxentscan data not found in %s (conda install -c bioconda maxentscan)' % sys.prefix)
    _DIR = hits[0]
    _ME2X5 = [float(l) for l in open(os.path.join(_DIR, 'me2x5')) if l.strip()]
    _SEQ5 = {l.strip(): n for n, l in enumerate(open(os.path.join(_DIR, 'splicemodels', 'splice5sequences'))) if l.strip()}
    _ACC = [[float(l) for l in open(os.path.join(_DIR, 'splicemodels', 'me2x3acc%d' % k)) if l.strip()] for k in range(1, 10)]


def _hash(seq):
    n = len(seq)
    return sum(_IDX[c] * _FOUR[n - i - 1] for i, c in enumerate(seq))


def donor(s):
    _load()
    s = s.upper()
    cons = _D1[s[3]] * _D2[s[4]] / (_BGD[s[3]] * _BGD[s[4]])
    return math.log2(cons * _ME2X5[_SEQ5[s[0:3] + s[5:9]]])


def acceptor(s):
    _load()
    s = s.upper()
    r = s[:18] + s[20:23]
    sc = [_ACC[0][_hash(r[0:7])], _ACC[1][_hash(r[7:14])], _ACC[2][_hash(r[14:21])], _ACC[3][_hash(r[4:11])], _ACC[4][_hash(r[11:18])],
          _ACC[5][_hash(r[4:7])], _ACC[6][_hash(r[7:11])], _ACC[7][_hash(r[11:14])], _ACC[8][_hash(r[14:18])]]
    cons = _A1[s[18]] * _A2[s[19]] / (_BGD[s[18]] * _BGD[s[19]])
    return math.log2(cons * sc[0] * sc[1] * sc[2] * sc[3] * sc[4] / (sc[5] * sc[6] * sc[7] * sc[8]))


def sites(seq, start=0, end=None):
    """Donor and acceptor scores of every GT / AG window that lies completely inside seq[start:end].
    Returns two lists of (junction position, score): for a donor the 0-based index of the last exon base,
    for an acceptor the 0-based index of the last intron base (the G of AG)."""
    seq = seq.upper()
    end = len(seq) if end is None else end
    don, acc = [], []
    for i in range(start, end - 8):
        if seq[i + 3:i + 5] == 'GT':
            don.append((i + 2, donor(seq[i:i + 9])))
    for i in range(start, end - 22):
        if seq[i + 18:i + 20] == 'AG':
            acc.append((i + 19, acceptor(seq[i:i + 23])))
    return don, acc
