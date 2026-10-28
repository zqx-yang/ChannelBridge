#!/usr/bin/env python3
"""Choose the channel-start ligand in each AlphaFill model, with evidence.

AlphaFill is not a cofactor service. It transplants **everything** it can
find from any homologous holo structure, and for these ten enzymes that is
15 to 82 ligand copies across 11 to 30 chemical types per model. The
inventory for one monoterpene synthase model includes CHES buffer (3CX/3FX),
a polyvinylpyrrolidone oligomer (QQ4), the antiviral baloxavir (E4Z) and
three bisphosphonate drugs (ibandronate, pamidronate, alendronate) alongside
the prenyl diphosphates. Picking by clash count or donor identity alone lands
on buffer.

So the ligand is chosen by **what the enzyme's channel actually carries**,
from the paper's own substrate and product columns, and only then checked
against what AlphaFill has. `CHOICE` records each decision and its reason;
this script verifies the resname is really in the fetched model and prints
AlphaFill's evidence for that transplant (donor entry, local identity, RMSD
and clash count) so the choice is auditable rather than asserted.

One choice is a deliberate surrogate. GES, MYS and OCS consume geranyl
diphosphate (GPP, C10), which AlphaFill did not transplant into any of the
three; it did transplant **GST**, geranyl S-thiolodiphosphate, the
non-hydrolysable C10 analogue crystallographers use in its place. That is a
real substitution and is marked as one. The obvious-looking alternative, DST,
is a C5 dimethylallyl analogue -- one isoprene unit short -- and would put
the start point in the wrong subsite.

    ./pick_ligands.py            # report and verify
    ./pick_ligands.py --tsv      # write the decision columns
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from find_structures import CLASH_PER_ATOM_MAX, ENTRIES, alphafill, read_decisions  # noqa: E402

#: enzyme -> (resname, exact|surrogate, why)
CHOICE = {
    # SAM-dependent methyltransferases: the methyl donor marks the catalytic
    # site the alkaloid is channelled to.
    "4'OMT": ("SAM", "exact", "S-adenosylmethionine, the methyl donor of this "
                              "4'-O-methyltransferase"),
    "7OMT": ("SAM", "exact", "as 4'OMT; a SAM-dependent O-methyltransferase"),
    "RNMT": ("SAM", "exact", "as 4'OMT; a SAM-dependent N-methyltransferase"),
    # trans-prenyltransferase: FPP is what leaves ERG20 and what the
    # downstream sesquiterpene synthases receive, so it is the cargo.
    "ERG20": ("FPP", "exact", "farnesyl diphosphate, ERG20's product and the "
                              "intermediate channelled onward (GPP, its "
                              "substrate, is also present)"),
    "ADS": ("FPP", "exact", "farnesyl diphosphate, the C15 substrate of "
                            "amorphadiene synthase"),
    "BBS": ("FPP", "exact", "farnesyl diphosphate, the C15 substrate of "
                            "bisabolene/bisabolol synthase"),
    # Monoterpene synthases consume GPP (C10). AlphaFill has no GPP in any of
    # the three models; GST is the C10 thiolo analogue.
    "GES": ("GST", "surrogate", "geranyl S-thiolodiphosphate, the C10 "
                                "non-hydrolysable analogue of GPP, which "
                                "AlphaFill did not transplant here"),
    "MYS": ("GST", "surrogate", "as GES"),
    "OCS": ("GST", "surrogate", "as GES"),
    # Cytochrome P450: the haem is the catalytic centre.
    "CTS": ("HEM", "exact", "protoporphyrin IX / Fe, the catalytic centre of "
                            "this P450"),
}

#: Codes seen in these models that are NOT cargo, with what they are. Kept so
#: that a future reader can see what picking by clash count would have chosen.
NOT_CARGO_SEEN = {
    "3CX": "CHES buffer", "3FX": "CHES buffer", "QQ4": "PVP oligomer",
    "E4Z": "baloxavir acid (antiviral)", "BFQ": "ibandronate",
    "210": "pamidronate", "AHD": "alendronate", "GDP": "guanosine diphosphate",
    "DPO": "bare diphosphate", "POP": "pyrophosphate (leaving group)",
    "FAR": "farnesyl, no diphosphate", "SNE": "sabinene (a product)",
    "DST": "C5 dimethylallyl analogue, one isoprene short of GPP",
}


#: The clash-per-atom threshold is defined once, in `find_structures`, because
#: it decides two things that must not drift apart: whether AlphaFill wins as
#: a SOURCE there, and whether its transplant can give a start point here.
#: Above it, the point comes from the residues lining the site instead
#: (`caver config --start-cavity`) -- the same machinery built for the
#: catalytic triads whose centre of gravity lands inside an atom.


def ligand_atom_count(path, resname):
    """Heavy atoms in one copy of *resname*, for normalising the clash count."""
    best = 0
    seen = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue
            f = line.split()
            if len(f) < 9 or f[5].strip().upper() != resname:
                continue
            key = tuple(f[6:9])
            seen[key] = seen.get(key, 0) + 1
    for n in seen.values():
        best = max(best, n)
    return best


def resnames_in(path):
    """HETATM residue names present in a fetched model."""
    out = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("HETATM"):
                f = line.split()
                if len(f) > 5:
                    out.add(f[5].strip().upper())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tsv", action="store_true")
    ap.add_argument("--no-evidence", action="store_true",
                    help="skip the AlphaFill query (offline). The clash "
                         "count comes from that query and decides start_mode, "
                         "so without it start_mode is written `unknown` and "
                         "build_inputs.py will refuse the row rather than "
                         "guess")
    args = ap.parse_args(argv)

    rows = [r for r in read_decisions(os.path.join(HERE, "structures.tsv"))
            if r["structure_source"] == "alphafill"]
    out, bad = [], 0
    if args.tsv:
        print("\t".join(("enzyme", "uniprot", "model", "channel_ligand",
                         "exactness", "donor", "donor_identity", "clashes",
                         "ligand_atoms", "clashes_per_atom", "start_mode",
                         "reason")))
    for r in rows:
        name, uni = r["enzyme"], r["uniprot"]
        model = os.path.join(ENTRIES, "alphafill_%s.cif" % uni)
        pick = CHOICE.get(name)
        if pick is None:
            print("%-7s NO CHOICE RECORDED" % name, file=sys.stderr)
            bad += 1
            continue
        lig, exactness, why = pick
        if not os.path.exists(model):
            print("%-7s model not fetched: %s" % (name, model), file=sys.stderr)
            bad += 1
            continue
        present = resnames_in(model)
        if lig not in present:
            print("%-7s %s IS NOT IN %s -- the choice is wrong or the model "
                  "changed" % (name, lig, os.path.basename(model)), file=sys.stderr)
            bad += 1
            continue
        donor, ident, clashes = "", "", ""
        if not args.no_evidence:
            fill = alphafill(uni)
            for h in (fill or {}).get("ligands", []):
                if h["ligand"] == lig:
                    donor = h["donor"]
                    ident = "%.3f" % h["identity"]
                    clashes = str(h["clashes"])
                    break
        natoms = ligand_atom_count(model, lig)
        # The clash count decides the start mode, so no clash count means no
        # decision -- NOT a clash count of zero. With `--no-evidence` the
        # AlphaFill query is skipped and `clashes` stays empty; reading that
        # as 0.0 silently turned the three worst transplants in this set
        # (GES, MYS, OCS at 2.5-5.2 clashes per atom) from `cavity` into
        # `ligand`, which would have started CAVER from the centroid of a
        # transplant that clashes into the protein at 0.56 A.
        if clashes == "":
            per_atom_s, start_mode = "", "unknown"
        else:
            per_atom = (float(clashes) / natoms) if natoms else 0.0
            per_atom_s = "%.1f" % per_atom
            start_mode = ("cavity" if per_atom > CLASH_PER_ATOM_MAX
                          else "ligand")
        if args.tsv:
            print("\t".join((name, uni, os.path.basename(model), lig, exactness,
                             donor, ident, clashes, str(natoms),
                             per_atom_s, start_mode, why)))
        else:
            print("%-7s %-4s %-10s donor %-5s identity %-6s clashes %-4s "
                  "/%2d atoms = %4.1f -> %s"
                  % (name, lig, exactness, donor or "-", ident or "-",
                     clashes or "-", natoms, per_atom, start_mode.upper()))
            print("        %s" % why)
            if start_mode == "cavity":
                print("        the transplant clashes too hard to mark a real "
                      "cavity, so the start point comes from the residues "
                      "lining it (--start-cavity), not its centroid")
            junk = sorted(present & set(NOT_CARGO_SEEN))
            if junk:
                print("        also in this model, and not cargo: %s"
                      % ", ".join("%s (%s)" % (j, NOT_CARGO_SEEN[j]) for j in junk[:5]))
        out.append(name)
    if not args.tsv:
        print("\n%d of %d assigned" % (len(out), len(rows)))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
