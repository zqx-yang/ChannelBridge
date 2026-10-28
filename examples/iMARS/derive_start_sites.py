#!/usr/bin/env python3
"""Where a channel starts, for the inputs that have no ligand to start from.

Fourteen of the 24 enzymes here carry their channelled species, or a usable
analogue of it, in the file: the start point is where that sits. **Eleven do
not** -- the five AlphaFold3 predictions, the five AlphaFold DB models, and
one experimental entry (3S9V) that is apo -- and for those the start point has
to come from the residues lining the site (`caver config --start-cavity`).

Which residues those are is the decision this script makes, and the one place
where it would be easiest to simply assert an answer. Asserting catalytic
residue numbers from memory is exactly the class of error this example set
already contains: a number that looks authoritative, is never checked, and
turns out to be about a different protein. So every site here is **located in
the model itself and then checked against something independent**, by one of
three routes, chosen by family:

**1. The type III PKS triad, by geometry.** Chalcone-synthase-family enzymes
hold a Cys-His-Asn triad. It is found by its geometry -- the Cys SG closest to
both a His NE2 and an Asn ND2 -- not by residue number, and the number that
comes out is then compared with the literature. On the CHS model this returns
**Cys164**, the catalytic cysteine of *Medicago* CHS2 by number, with His at
3.7 A and Asn at 4.2 A; on CUS it returns Cys174. The first version of this
search ranked candidates by burial instead and returned Cys341 on CHS -- a
buried cysteine that is not the catalytic one -- so the ranking matters and
is part of the method, not a detail.

**2. The class I terpene synthase metal motifs, by pattern and geometry.**
These enzymes bind a trinuclear Mg/Mn cluster between an aspartate-rich
``DDxxD`` motif and an ``NSE/DTE`` motif across the cleft. Both are found by
pattern and then confirmed by distance to the metal in holo members of the
same family: on taxadiene synthase (3P5P) the search returns exactly one of
each, ``DDMAD`` at 613 and ``NDTKTYQAE`` at 757, **2.2 and 3.7 A from the
Mg**; on limonene synthase (2ONG) ``DDIYD`` at 352 and ``DDLGTSVEE`` at 496,
**2.2 and 2.4 A from the Mn**. One hit each in 745 and 543 residues, no false
positives, and the metal was never used to find them.

**A fourth route was tried for the prenyltransferase fold and refuted.** The
squalene-synthase and trans-prenyltransferase folds hold their diphosphate
between two aspartate-rich motifs rather than on a DDxxD/NSE pair, and the
pattern matches too often to use alone -- six times in the SQS model -- so the
pair was chosen by geometry, as the two whose carboxylates come closest. It
validated on two holo members (3WEF 3.1 A from the ligand, 2E8T 2.8 A, the
latter over two decoys) and is still **wrong**: mapped through GPPS's own holo
homolog, the second motif of the real site is ``DDILD`` at 307 in the model --
3KRF chain D's ``DDILDVT`` at 221, 3.4 A from the ligand -- and the geometric
criterion picks ``DIASEGD`` at 245 instead. Two validations agreeing did not
make the criterion right; a third case with an independent answer broke it.
Both prenyltransferases have a holo homolog, so they use route 3, and the
motif-pair code is gone rather than left in as an option.

**3. A holo homolog, by alignment.** For the rest, the site is mapped from a
homologous holo structure: align it to the model and take the model residues
contacting its ligand. The homologs used are the very entries this set
originally mistook for the enzymes themselves -- which is what they are
legitimately evidence for. A homolog is not a stand-in for the structure; it
is perfectly good evidence about where a conserved site is.

    ./derive_start_sites.py              # report, with the evidence
    ./derive_start_sites.py --tsv        # write start_sites.tsv
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, HERE)

from seqalign import AA3, align, seq_of  # noqa: E402,F401

from channelbridge.geometry import superpose  # noqa: E402
from channelbridge.io import read_mmcif, read_pdb  # noqa: E402

#: The aspartate-rich motif and the NSE/DTE motif of a class I terpene
#: synthase. Validated against the metal in two holo members; see the module
#: docstring for the numbers.
DDXXD = re.compile(r"DD..D")
NSE = re.compile(r"[ND]D..[ST]...E")

#: Residues of the model within this distance of a mapped ligand count as
#: lining its site.
CONTACT = 5.0

#: How each apo input's site is located, and why that route. The family is the
#: reason; it is recorded because "which residues line the site" is a claim
#: about the protein, not about the file.
FAMILY = {
    "CHS": ("pks", "type III PKS (chalcone synthase); Cys-His-Asn triad"),
    "CUS": ("pks", "type III PKS (curcuminoid synthase); as CHS"),
    "STS": ("pks", "type III PKS (stilbene synthase); as CHS"),
    "AS": ("terpene", "class I diterpene synthase (abietadiene synthase); "
                      "trinuclear metal cluster"),
    "NES": ("terpene", "class I sesquiterpene synthase (nerolidol synthase)"),
    "GPPS": ("homolog", "trans-prenyltransferase (geranyl diphosphate "
                        "synthase); mapped, because the motif route is "
                        "refuted for this fold -- see the module docstring"),
    "GGPPS": ("homolog", "trans-prenyltransferase (GGPP synthase); as GPPS"),
    "4CL": ("homolog", "ANL-superfamily adenylate-forming enzyme; no simple "
                       "sequence motif marks the site, so it is mapped"),
    "SQS": ("homolog", "squalene synthase; mapped from a holo homolog"),
    "RAS": ("homolog", "acyltransferase/reductase; mapped"),
    "REPI": ("homolog", "NADPH-dependent reductase; mapped"),
}

#: ``enzyme -> (pdb id, chain, ligand resname)`` for the homolog route. The
#: ligand is the one whose contacts define the site; it is checked to be
#: present before use.
HOMOLOG = {
    "4CL": ("5BSM", "A", "ATP"),
    # 3KRF is a heterotetramer and only chain D has the substrate analogue
    # bound, so the chain is part of the decision, not a default.
    "GPPS": ("3KRF", "D", "DST"),
    "GGPPS": ("2E8T", "A", "FPS"),
    "SQS": ("3WEF", "A", "FPS"),
    "RAS": ("5FAL", "A", "COA"),
    "REPI": ("1ZGD", "A", "NAP"),
}

ENTRIES = os.path.join(HERE, "_entries")

#: AlphaFold DB model version, defined once (find_structures.py
#: downloads it, build_inputs.py and this module name the file).
#: Three copies of "6" used to drift independently.
from find_structures import AFDB_VERSION  # noqa: E402


def read_any(path):
    return read_pdb(path) if path.lower().endswith(".pdb") else read_mmcif(path)

def atom_xyz(atoms, chain, resnums, elements=None):
    out = []
    want = set(resnums)
    for x in atoms.atoms:
        if x.chain != chain or x.resseq not in want or x.hetatm:
            continue
        if elements and (x.element or "").upper() not in elements:
            continue
        out.append(x.xyz)
    return np.asarray(out, dtype=float)


def pks_triad(atoms, chain):
    """The Cys-His-Asn triad of a type III PKS, found by geometry.

    Ranked by how close the Cys SG is to BOTH a His NE2 and an Asn ND2, which
    is the catalytic geometry. Ranking by burial instead returns the most
    buried cysteine, which on the CHS model is Cys341 and wrong.
    """
    his_res = [(x.resseq, np.asarray(x.xyz)) for x in atoms.atoms
               if x.chain == chain and not x.hetatm
               and x.resname.upper() == "HIS" and x.name == "NE2"]
    asn_res = [(x.resseq, np.asarray(x.xyz)) for x in atoms.atoms
               if x.chain == chain and not x.hetatm
               and x.resname.upper() == "ASN" and x.name == "ND2"]
    best = None
    for x in atoms.atoms:
        if (x.chain != chain or x.hetatm or x.resname.upper() != "CYS"
                or x.name != "SG"):
            continue
        p = np.asarray(x.xyz)
        if not his_res or not asn_res:
            continue
        hn, hp = min(his_res, key=lambda r: np.linalg.norm(r[1] - p))
        an, ap = min(asn_res, key=lambda r: np.linalg.norm(r[1] - p))
        dh = float(np.linalg.norm(hp - p))
        da = float(np.linalg.norm(ap - p))
        if best is None or dh + da < best[0]:
            best = (dh + da, x.resseq, dh, da, p, hn, an)
    return best


def atom_xyz_named(atoms, chain, resname, name):
    return np.asarray([x.xyz for x in atoms.atoms
                       if x.chain == chain and not x.hetatm
                       and x.resname.upper() == resname and x.name == name],
                      dtype=float)


def terpene_motifs(atoms, chain):
    """``[(label, motif, first residue, residue numbers)]``, by pattern."""
    s, nums = seq_of(atoms, chain)
    out = []
    for label, rx in (("DDxxD", DDXXD), ("NSE/DTE", NSE)):
        for m in rx.finditer(s):
            sel = nums[m.start():m.end()]
            out.append((label, m.group(0), sel[0], sel))
    return out

#: A mapping whose superposition is worse than this is refused rather than
#: reported. A site derived from a bad superposition sends CAVER to the wrong
#: place, and it does it silently: the output is a tunnel, not an error.
MAX_MAP_RMSD = 3.0
#: And it must rest on at least this many residue pairs.
MIN_MAP_PAIRS = 50


def _superpose_trimmed(src, dst, *, rounds=5, keep=2.5):
    """Superpose, dropping the pairs that do not fit, and report what was kept.

    A homolog agrees with the model over its conserved core and diverges in
    loops; including the loops drags the fit and moves the mapped ligand. Each
    round refits on the pairs within *keep* times the median deviation.
    """
    idx = np.arange(len(src))
    R, t, rmsd = superpose(src, dst)
    for _ in range(rounds):
        dev = np.linalg.norm((src[idx] @ R.T + t) - dst[idx], axis=1)
        cut = max(1.0, keep * float(np.median(dev)))
        sel = idx[dev <= cut]
        if len(sel) < MIN_MAP_PAIRS or len(sel) == len(idx):
            break
        idx = sel
        R, t, rmsd = superpose(src[idx], dst[idx])
    return R, t, rmsd, len(idx)


def homolog_site(model, model_chain, hom, hom_chain, ligand):
    """Model residues lining the site that *ligand* occupies in the homolog.

    The homolog is superposed on the model over their aligned residues, so the
    ligand lands where the conserved site is in the model's own frame, and the
    residues within CONTACT of it are returned with the superposition RMSD and
    the sequence identity so the mapping can be judged.
    """
    ms, mnums = seq_of(model, model_chain)
    hs, hnums = seq_of(hom, hom_chain)
    pairs, ident = align(hs, ms)
    hxyz, mxyz, hn, mn = [], [], [], []
    hca = {x.resseq: x.xyz for x in hom.atoms
           if x.chain == hom_chain and x.name == "CA" and not x.hetatm}
    mca = {x.resseq: x.xyz for x in model.atoms
           if x.chain == model_chain and x.name == "CA" and not x.hetatm}
    for i, j in pairs:
        a, b = hnums[i], mnums[j]
        if a in hca and b in mca:
            hxyz.append(hca[a])
            mxyz.append(mca[b])
            hn.append(a)
            mn.append(b)
    if len(hxyz) < MIN_MAP_PAIRS:
        return {"refused": "only %d residue pair(s) aligned, under the %d "
                           "minimum" % (len(hxyz), MIN_MAP_PAIRS)}
    R, t, rmsd, n_kept = _superpose_trimmed(np.asarray(hxyz), np.asarray(mxyz))
    if rmsd > MAX_MAP_RMSD:
        return {"refused": "superposition RMSD %.2f A over %d kept pair(s) "
                           "exceeds the %.1f A limit, so the mapped site would "
                           "not be this protein's site"
                           % (rmsd, n_kept, MAX_MAP_RMSD),
                "identity": ident, "rmsd": rmsd, "n_aligned": len(hxyz)}
    lig = np.asarray([x.xyz for x in hom.atoms
                      if x.hetatm and x.resname.upper() == ligand.upper()],
                     dtype=float)
    if not len(lig):
        return {"refused": "%s is not in the homolog" % ligand}
    lig_moved = lig @ R.T + t
    lining = []
    for resnum in mnums:
        pts = atom_xyz(model, model_chain, [resnum])
        if not len(pts):
            continue
        if np.linalg.norm(pts[:, None, :] - lig_moved[None, :, :], axis=2).min() <= CONTACT:
            lining.append(resnum)
    return {"identity": ident, "rmsd": rmsd, "n_aligned": len(hxyz),
            "n_kept": n_kept, "lining": lining,
            "centre": lig_moved.mean(axis=0)}


def structure_for(row):
    s, uni = row["structure_source"], row["uniprot"]
    if s == "alphafill":
        return os.path.join(ENTRIES, "alphafill_%s.cif" % uni)
    if s == "afdb":
        return os.path.join(ENTRIES, "afdb_%s_v%s.cif" % (uni, AFDB_VERSION))
    if s in ("pdb", "pdb-other"):
        return os.path.join(ENTRIES, "%s.cif" % row["use_entry"])
    if s == "af3":
        return os.path.join(ENTRIES, "af3",
                            "af3_%s.cif" % row["enzyme"].replace("'", "p"))
    return None


def main(argv=None):
    from find_structures import read_decisions

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tsv", action="store_true")
    ap.add_argument("--only", action="append")
    args = ap.parse_args(argv)

    rows = read_decisions(os.path.join(HERE, "structures.tsv"))
    if args.only:
        rows = [r for r in rows if r["enzyme"] in args.only]
    if args.tsv:
        print("\t".join(("enzyme", "family", "route", "chain", "residues",
                         "evidence")))
    bad = 0
    for r in rows:
        name = r["enzyme"]
        if name not in FAMILY:
            continue
        family, why = FAMILY[name]
        path = structure_for(r)
        if not path or not os.path.exists(path):
            print("%-7s structure not found: %s" % (name, path), file=sys.stderr)
            bad += 1
            continue
        st = read_any(path)
        chain = sorted({x.chain for x in st.atoms if not x.hetatm})[0]
        residues, evidence = [], ""
        # Which route actually found the site. It used to be a second copy
        # of `family`, so the column could never say anything `family` did
        # not -- and the useful distinction is WHICH homolog a site came
        # from, which only this records.
        route = {"pks": "pks-triad"}.get(family, family)
        if family == "pks":
            t = pks_triad(st, chain)
            if t is None:
                evidence = "no Cys-His-Asn triad found"
                bad += 1
            else:
                _, cys, dh, da, _, his_n, asn_n = t
                # All THREE triad residues, not just the Cys: a single residue
                # gives --start-cavity one atom cluster to centre on, and a
                # centroid inside an atom is the failure that machinery exists
                # to avoid. The triad brackets the cavity.
                residues = sorted({cys, his_n, asn_n})
                evidence = ("catalytic triad by geometry: Cys%d, His%d at "
                            "%.1f A, Asn%d at %.1f A"
                            % (cys, his_n, dh, asn_n, da))
        elif family == "terpene":
            route = "terpene-motifs"
            ms = terpene_motifs(st, chain)
            if not ms:
                evidence = "neither DDxxD nor NSE/DTE found"
                bad += 1
            else:
                for _label, _motif, _first, sel in ms:
                    residues.extend(sel)
                evidence = "; ".join("%s %s at %d" % (lab, mo, f)
                                     for lab, mo, f, _ in ms)
        else:
            pid, hch, lig = HOMOLOG[name]
            route = "homolog-%s" % pid
            hpath = None
            for cand in ("%s/%s.pdb" % (ENTRIES, pid), "%s/%s.cif" % (ENTRIES, pid)):
                if os.path.exists(cand):
                    hpath = cand
                    break
            if hpath is None:
                evidence = "homolog %s not fetched" % pid
                bad += 1
            else:
                got = homolog_site(st, chain, read_any(hpath), hch, lig)
                if got is None or "refused" in got:
                    evidence = "REFUSED: %s" % (
                        got["refused"] if got else "could not map %s from %s"
                        % (lig, pid))
                    bad += 1
                else:
                    residues = got["lining"]
                    evidence = ("%d residue(s) lining %s mapped from %s chain %s "
                                "(%.0f%% identity, %d pairs aligned, %d kept, "
                                "RMSD %.2f A)"
                                % (len(residues), lig, pid, hch,
                                   100 * got["identity"], got["n_aligned"],
                                   got["n_kept"], got["rmsd"]))
        # No slice. An earlier version kept 40 residues while `evidence` on
        # the same row counted all of them, so a larger site would have
        # described one cavity and seeded CAVER from another.
        rng = ",".join(str(x) for x in residues)
        if args.tsv:
            print("\t".join((name, family, route, chain, rng, evidence)))
        else:
            print("%-7s %-8s chain %-3s %s" % (name, family, chain, evidence))
            print("        %s: %s" % (why, rng or "(none)"))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
