"""Pairwise sequence alignment, shared by the scripts that need it.

Here rather than in either caller because both `find_structures.py` (is this
entry's protein the same protein?) and `derive_start_sites.py` (where is the
conserved site?) need it, and a second copy would drift from the first.

No dependency: this example set is meant to run from a clean checkout.
"""

from __future__ import annotations

import numpy as np

__all__ = ["AA3", "AA_ORDER", "BLOSUM62", "align", "seq_of"]

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def seq_of(atoms, chain):
    """``(one-letter sequence, residue numbers)`` for one polymer chain."""
    res = {}
    for x in atoms.atoms:
        if x.chain == chain and not x.hetatm and x.resname.upper() in AA3:
            res.setdefault(x.resseq, x.resname.upper())
    nums = sorted(res)
    return "".join(AA3[res[n]] for n in nums), nums


#: BLOSUM62, the standard substitution matrix, as rows in AA_ORDER.
#:
#: Identity scoring was enough for the 52-74 % identity mappings and not for
#: SQS: the paper's SQS is the *Priestia megaterium* enzyme and the only
#: cross-referenced structure is HUMAN squalene synthase, where an
#: identity-scored local alignment found a block of just 35 residues -- under
#: the minimum, so the mapping was refused. A substitution matrix is what
#: makes a distant homolog alignable at all, and refusing a site for want of
#: one is refusing it for the wrong reason.
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
BLOSUM62 = {}
for _i, _row in enumerate(_B62.strip().split("\n")):
    for _j, _v in enumerate(_row.split()):
        BLOSUM62[(AA_ORDER[_i], AA_ORDER[_j])] = int(_v)


def _score(x, y):
    return BLOSUM62.get((x, y), -1)


def align(a, b, *, gap=-10, extend=-1):
    """Smith-Waterman LOCAL alignment, returning index pairs and identity.

    Local, not global, and that is the whole point. A global alignment forces
    the ends to meet, so when the homolog covers only part of the model it
    spreads the correspondence across the whole chain and the result is not an
    alignment at all. REPI is 900 residues and its homolog 1ZGD is 308 -- the
    reductase domain of a P450-reductase fusion -- and global alignment of
    those two produced a "site" with a superposition RMSD of **31.6 A**, which
    the first version of this script printed as a result. A number that large
    is not a near miss; it is a mapping onto the wrong part of the protein.

    Scored with BLOSUM62 and affine gaps. Written here rather than taken from
    a dependency so this set runs from a clean checkout.
    """
    n, m = len(a), len(b)
    sc = np.zeros((n + 1, m + 1), dtype=float)
    ptr = np.zeros((n + 1, m + 1), dtype=np.int8)
    best, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = sc[i], sc[i - 1]
        for j in range(1, m + 1):
            d = prev[j - 1] + _score(ai, b[j - 1])
            u = prev[j] + (extend if ptr[i - 1, j] == 1 else gap)
            left = row[j - 1] + (extend if ptr[i, j - 1] == 2 else gap)
            v, p = d, 0
            if u > v:
                v, p = u, 1
            if left > v:
                v, p = left, 2
            if v <= 0:
                row[j], ptr[i, j] = 0.0, 3
                continue
            row[j], ptr[i, j] = v, p
            if v > best:
                best, bi, bj = v, i, j
    pairs = []
    i, j = bi, bj
    while i > 0 and j > 0 and ptr[i, j] != 3 and sc[i, j] > 0:
        p = ptr[i, j]
        if p == 0:
            i, j = i - 1, j - 1
            pairs.append((i, j))
        elif p == 1:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    ident = sum(1 for i, j in pairs if a[i] == b[j])
    return pairs, (ident / max(1, len(pairs)))



