#!/usr/bin/env python3
"""Check the predicted 4CL against the one experimental structure of it.

4CL is the only enzyme in this set taken from a fresh AlphaFold3 prediction
rather than from a database model, and the reason is specific: its channel
depends on the angle between its two domains, and that is the one thing this
project measured AlphaFold2 to be unreliable at. A prediction chosen for that
reason has to be checked, so this does the check in two parts -- and the two
parts are not equally strong, which is the point of reporting them apart.

**Part 1, the alignment.** 3TSY is the only experimental structure containing
this protein. It is a 4CL::stilbene-synthase FUSION at 3.1 A, and the last 102
residues of 4CL were never modelled in it. So it can confirm the fold and the
geometry of the large N-terminal domain, residue for residue, against real
density -- and it says nothing whatever about where the C-terminal domain
sits, because that is the part it is missing. The verification is partial by
construction and the missing part is exactly the part in question. Saying so
is not a hedge; a reader who takes a 1 A RMSD here as "the model is right"
would be taking it as evidence about a domain that was never compared.

**Part 2, the prediction's own PAE.** For the inter-domain angle, which part 1
cannot reach, the available evidence is AlphaFold3's own predicted aligned
error between the two domains. The domain split is taken from the PAE matrix
itself rather than asserted from a sequence boundary, so the number describes
the blocks the prediction actually has. The yardstick is this project's own
measurement on iMARS fusions: inter-domain PAE of 28 A against a 31 A ceiling
was a draw from the prior and gave a 0.36 DO spread, while intra-domain PAE
was 3.0-5.7 A.

    ./verify_4cl.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, HERE)

from channelbridge.geometry import superpose  # noqa: E402
from channelbridge.io import read_mmcif  # noqa: E402

#: The fusion, and the chain 4CL is in it.
FUSION = os.path.join(HERE, "_entries", "3TSY.cif")
FUSION_CHAIN = "A"
#: The prediction, as `make_af3_inputs.py --collect` leaves it.
PREDICTION = os.path.join(HERE, "_entries", "af3", "af3_4CL.cif")
CONFIDENCES = os.path.join(HERE, "_entries", "af3", "af3_4CL_confidences.json")

#: A residue-name agreement below this means the two numberings were not
#: brought into register and no RMSD computed from them means anything.
MIN_AGREEMENT = 0.95


def ca_by_residue(atoms, chain):
    """``{resseq: (resname, xyz)}`` for the CA atoms of one chain."""
    out = {}
    for a in atoms.atoms:
        if a.chain == chain and a.name == "CA" and not a.hetatm:
            out.setdefault(a.resseq, (a.resname.upper(), np.asarray(a.xyz)))
    return out


def best_offset(fusion_ca, pred_ca):
    """The numbering offset between the two files, found rather than assumed.

    RCSB says 3TSY's entity positions 24-584 are this protein's residues
    1-561, which implies an offset of -23. That is a statement about the
    entity numbering and this file's author-numbering need not match it, so
    the offset is chosen by residue-name agreement and the agreement is
    reported. An offset that merely maximises a bad agreement is caught by
    MIN_AGREEMENT rather than used.
    """
    best = (0.0, None, 0)
    for off in range(-200, 201):
        hits = tot = 0
        for seq, (name, _) in fusion_ca.items():
            other = pred_ca.get(seq + off)
            if other is None:
                continue
            tot += 1
            hits += (other[0] == name)
        if tot >= 50 and hits / tot > best[0]:
            best = (hits / tot, off, tot)
    return best


#: A block smaller than this is a terminus, not a domain.
MIN_BLOCK = 100
#: A block whose residues are not confidently placed relative to EACH OTHER
#: is not a domain either, whatever its PAE to the rest.
MAX_INTRA = 8.0


def domain_split(pae, token_res):
    """Split the chain in two where the PAE says two bodies are, or None.

    A domain is a block that is placed confidently **within itself** and not
    confidently **relative to the other**. Scoring only the second half --
    maximising inter-block minus mean intra-block PAE -- looks right and is
    not: a floppy 80-residue terminus has a high PAE to everything including
    itself, so it wins. Measured on the four models already predicted here,
    that version put CUS's boundary 322 residues from the truth with a
    *higher* score than CHS's correct one, and invented a split in
    single-chain GGPPS. Requiring both blocks to be internally confident and
    scoring by ``inter - max(intra_a, intra_b)`` puts both dimers' splits on
    their true chain boundary exactly.

    The score is a weak discriminator and is reported as a number against
    the controls in :func:`calibration`, never as a verdict.
    """
    n = len(pae)
    if n < 2 * MIN_BLOCK:
        return None
    best = None
    for k in range(MIN_BLOCK, n - MIN_BLOCK + 1):
        a, b = slice(0, k), slice(k, n)
        ia, ib = pae[a, a].mean(), pae[b, b].mean()
        if max(ia, ib) > MAX_INTRA:
            continue
        inter = 0.5 * (pae[a, b].mean() + pae[b, a].mean())
        score = inter - max(ia, ib)
        if best is None or score > best[0]:
            best = (score, k, ia, ib, inter)
    if best is None:
        return None
    score, k, ia, ib, inter = best
    return {"at": int(token_res[k]), "index": k, "score": score,
            "intra_a": ia, "intra_b": ib, "inter": inter,
            "n_a": k, "n_b": n - k}


def calibration(af3_dir):
    """The split score on the models whose answer is known, computed now.

    A dimer has a two-body split whose position is known exactly -- the
    chain boundary -- and a single-chain single-domain enzyme has nothing to
    find, so the two together give a scale. Computing it instead of quoting
    it keeps that scale honest when a model is repredicted.

    How many of each there are is NOT asserted, because it changed: CUS was
    predicted as a dimer until its own deposited entry (3ALE) was found, and
    the calibration quietly went from two dimer controls to one. The caller
    prints the counts so a thin calibration is visible rather than implied.
    """
    import glob
    rows = []
    for path in sorted(glob.glob(os.path.join(af3_dir, "af3_*_confidences.json"))):
        name = os.path.basename(path)[len("af3_"):-len("_confidences.json")]
        if name == "4CL":
            continue
        try:
            with open(path, encoding="utf-8") as fh:
                conf = json.load(fh)
        except (OSError, ValueError):
            continue
        pae = np.asarray(conf["pae"], dtype=float)
        chains = conf["token_chain_ids"]
        uniq = sorted(set(chains))
        truth = chains.index(uniq[1]) if len(uniq) > 1 else None
        s = domain_split(pae, conf["token_res_ids"])
        rows.append((name, len(uniq), s, truth))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fusion", default=FUSION)
    ap.add_argument("--prediction", default=PREDICTION)
    ap.add_argument("--confidences", default=CONFIDENCES)
    args = ap.parse_args(argv)

    for path, what in ((args.fusion, "the fusion structure"),
                       (args.prediction, "the prediction")):
        if not os.path.exists(path):
            raise SystemExit(
                "%s is not there: %s\nRun `./make_af3_inputs.py --collect` "
                "once the AF3 job has finished." % (what, path))

    fusion = read_mmcif(args.fusion)
    pred = read_mmcif(args.prediction)
    pred_chain = sorted({a.chain for a in pred.atoms if not a.hetatm})[0]
    f_ca = ca_by_residue(fusion, FUSION_CHAIN)
    p_ca = ca_by_residue(pred, pred_chain)
    print("fusion      %s chain %s: %d residues with CA, numbered %d-%d"
          % (os.path.basename(args.fusion), FUSION_CHAIN, len(f_ca),
             min(f_ca), max(f_ca)))
    print("prediction  %s chain %s: %d residues, numbered %d-%d"
          % (os.path.basename(args.prediction), pred_chain, len(p_ca),
             min(p_ca), max(p_ca)))

    agree, off, tot = best_offset(f_ca, p_ca)
    if off is None or agree < MIN_AGREEMENT:
        raise SystemExit(
            "the two numberings could not be brought into register (best "
            "agreement %.0f%% at offset %s over %d residues). An RMSD from "
            "this would be a number about the wrong residues."
            % (100 * agree, off, tot))
    print("\nregister    offset %+d, residue names agree on %d of %d (%.0f%%)"
          % (off, round(agree * tot), tot, 100 * agree))

    pairs = [(seq, seq + off) for seq in sorted(f_ca)
             if seq + off in p_ca and p_ca[seq + off][0] == f_ca[seq][0]]
    mob = np.asarray([p_ca[j][1] for _, j in pairs])
    tgt = np.asarray([f_ca[i][1] for i, _ in pairs])
    _, _, rmsd = superpose(mob, tgt)
    lo, hi = pairs[0][1], pairs[-1][1]
    print("            compared on %d CA pairs, prediction residues %d-%d"
          % (len(pairs), lo, hi))
    print("\nPART 1  fold and the observed domain")
    print("  CA RMSD %.2f A over %d residues" % (rmsd, len(pairs)))

    # Where the agreement is worst, which is more informative than the mean.
    R, t, _ = superpose(mob, tgt)
    dev = np.linalg.norm((mob @ R.T + t) - tgt, axis=1)
    worst = np.argsort(dev)[::-1][:5]
    print("  per-residue deviation: median %.2f, 90th pct %.2f, max %.2f A"
          % (np.median(dev), np.percentile(dev, 90), dev.max()))
    print("  worst: %s" % ", ".join(
        "%s%d %.1f A" % (p_ca[pairs[i][1]][0], pairs[i][1], dev[i]) for i in worst))

    missing = sorted(set(range(min(p_ca), max(p_ca) + 1)) - {j for _, j in pairs})
    runs = []
    for n in missing:
        if runs and n == runs[-1][-1] + 1:
            runs[-1].append(n)
        else:
            runs.append([n])
    print("  NOT COMPARED: %d residue(s) of the prediction have no observed "
          "counterpart" % len(missing))
    for r in sorted(runs, key=len, reverse=True)[:4]:
        print("     %d-%d (%d residues)" % (r[0], r[-1], len(r)))
    print("  The largest of those runs is 4CL's C-terminal region, which the\n"
          "  fusion never modelled. Part 1 says nothing about where it sits,\n"
          "  and that is the domain the channel depends on.")

    if not os.path.exists(args.confidences):
        print("\nPART 2  skipped: no %s" % os.path.basename(args.confidences))
        return 0
    with open(args.confidences, encoding="utf-8") as fh:
        conf = json.load(fh)
    pae = np.asarray(conf["pae"], dtype=float)
    token_res = conf["token_res_ids"]
    print("\nPART 2  the inter-domain angle, which part 1 cannot reach")
    print("  PAE matrix %d x %d, mean %.1f A" % (pae.shape[0], pae.shape[1], pae.mean()))
    split = domain_split(pae, token_res)
    if split is None:
        print("  no two-block split passes: either the chain is too short or "
              "no pair of blocks is each confidently placed within itself")
        return 0
    print("  best two-block split at residue %d (%d + %d residues), taken "
          "from the PAE itself" % (split["at"], split["n_a"], split["n_b"]))
    print("  intra-block PAE %.1f / %.1f A   INTER-block PAE %.1f A   "
          "score %.1f" % (split["intra_a"], split["intra_b"],
                          split["inter"], split["score"]))

    print("\n  The score is a weak discriminator, so here it is on the models "
          "whose\n  answer is known (dimers have a real two-body split at the "
          "chain\n  boundary; the single-chain enzymes have none):")
    rows = calibration(os.path.dirname(args.confidences))
    n_dimer = sum(1 for _n, nch, _s, _t in rows if nch > 1)
    print("    (%d control(s): %d dimer, %d single-chain. A thin calibration "
          "is a weak scale -- say so rather than reading a verdict off it.)"
          % (len(rows), n_dimer, len(rows) - n_dimer))
    for name, nch, s, truth in rows:
        if s is None:
            print("    %-6s %d chain(s): no split passes" % (name, nch))
            continue
        where = ("found the chain boundary exactly" if truth is not None
                 and s["index"] == truth else
                 "%d residues from the chain boundary" % abs(s["index"] - truth)
                 if truth is not None else "nothing to find")
        print("    %-6s %d chain(s): score %.1f -- %s"
              % (name, nch, s["score"], where))
    print("\n  For the inter-block number itself, the scale from this "
          "project's own\n  measurement: AF2 fusion models had inter-domain "
          "PAE 28 A against a\n  31 A ceiling -- a draw from the prior, 0.36 "
          "DO spread -- while\n  intra-domain PAE was 3.0-5.7 A.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
