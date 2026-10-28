#!/usr/bin/env python3
"""The one table the pipeline runs from: what to prepare, and how.

`structures.tsv` says which coordinates each enzyme is represented by.
`af_ligands.tsv` says which ligand marks the channel in the AlphaFill models.
`derive_start_sites.py` says which residues line the site where there is no
ligand at all. This merges the three into `inputs.tsv`, one row per enzyme,
with everything `chbridge prep` and `chbridge caver config` need -- and it
checks each decision against the file rather than trusting the table it came
from.

Two decisions live here because they are not in any of those three:

**Which chains.** A monomer is the cheap answer and is wrong wherever the
active site is completed across the interface. The type III PKS cavity is
built partly from the partner subunit -- verified directly on the CHS model,
where the catalytic Cys164's 12 A neighbourhood contains the *partner's*
Met137 -- so those are prepared as dimers. A class I terpene synthase or a
*trans*-prenyltransferase carries a complete site in one subunit, so those are
monomers even where the crystal has four chains.

**What to trim.** A disordered terminus is a wall to a sphere-rolling
algorithm exactly as a water is. GGPPS's model has an 18-residue N-terminal
tail that swings 9 A across AF3's samples, around a core of 326 of 358, and it
sits on the outside of the protein where a channel mouth would be.

    ./build_inputs.py            # report
    ./build_inputs.py --tsv      # write inputs.tsv
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, HERE)

from find_structures import AFDB_VERSION, ENTRIES, read_decisions  # noqa: E402

from channelbridge.io import read_mmcif, read_pdb  # noqa: E402

#: Fewer residues than this and a chain is not a protomer.
MIN_PROTOMER = 30

#: Chains to keep, and why. Absent means "the first polymer chain", which is
#: only right for a single-chain file and is asserted, not assumed: a file
#: with more chains than this says is an error, not a default.
CHAINS = {
    # Type III PKS: the cavity is completed by the partner subunit.
    "CHS": (["A", "B"], "type III PKS; partner Met137 lines the Cys164 site"),
    "CUS": (["A", "B"], "type III PKS, as CHS"),
    "STS": (["A"], "type III PKS, but the source is a single-chain model; the "
                   "dimer partner is absent and the site is therefore "
                   "incomplete -- recorded as a known limitation"),
    "BAS": (["A", "B"], "type III PKS; 3A5R is the dimer"),
    # One complete site per subunit: a monomer is enough.
    "AS": (["A"], "class I diterpene synthase; 3S9V has 4 chains and each "
                  "carries a complete site"),
    "TS": (["A"], "class I diterpene synthase; single chain"),
    "LIC": (["A"], "class I terpene synthase; 2ONG is a dimer but each "
                   "subunit's site is complete"),
    "GPPS": (["A"], "trans-prenyltransferase"),
    "GGPPS": (["A"], "trans-prenyltransferase"),
    "SQS": (["A"], "squalene synthase; site within one subunit"),
    "4CL": (["A"], "ANL-superfamily enzyme; monomeric"),
}

#: Residues to drop before measuring, and why. Trimming changes which atoms
#: exist, so each one is named.
TRIM = {
    "GGPPS": ([(1, 18)], "an 18-residue N-terminal tail swinging 9 A across "
                         "AF3's samples, outside a core of 326 of 358, in the "
                         "region a channel mouth would occupy"),
}

#: Cargo for the sources `af_ligands.tsv` does not cover: the deposited
#: entries. Each is the species whose site gives the start point, and each is
#: checked to be in the file.
EXTRA_CARGO = {
    "BAS": ("HC4", "4-coumaric acid, the substrate, bound in both subunits"),
    "LIC": ("F3P", "2-fluorolinalyl diphosphate: the deposited component in "
                   "BOTH 2ONG and 2ONH despite 2ONG's title saying "
                   "fluorogeranyl -- the analogue isomerised in the crystal"),
    "TS": ("A3C", "the substrate analogue in 3P5P"),
}

#: Ions and metals to keep because they are part of the site, not furniture.
KEEP_ALSO = {
    "LIC": (["MN"], "the catalytic trinuclear manganese cluster"),
    "TS": (["MG"], "the catalytic trinuclear magnesium cluster"),
}


#: Residues this close to the cargo count as lining its site.
CONTACT = 5.0


def lining_residues(atoms, chains, resname):
    """Residues of *chains* within CONTACT of any copy of *resname*."""
    import numpy as np

    lig = np.asarray([x.xyz for x in atoms.atoms
                      if x.hetatm and x.resname.upper() == resname.upper()],
                     dtype=float)
    if not len(lig):
        return []
    keep = set(chains)
    by_res, out = {}, []
    for x in atoms.atoms:
        if x.hetatm or x.chain not in keep:
            continue
        by_res.setdefault((x.chain, x.resseq), []).append(x.xyz)
    for (_ch, num), pts in sorted(by_res.items(), key=lambda kv: kv[0][1]):
        p = np.asarray(pts, dtype=float)
        if np.linalg.norm(p[:, None, :] - lig[None, :, :], axis=2).min() <= CONTACT:
            out.append(num)
    return sorted(set(out))


def read_tsv(path):
    rows, header = [], None
    if not os.path.exists(path):
        return rows
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = f
            continue
        f += [""] * (len(header) - len(f))
        rows.append(dict(zip(header, f)))
    return rows


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


def read_any(path):
    return read_pdb(path) if path.lower().endswith(".pdb") else read_mmcif(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tsv", action="store_true")
    args = ap.parse_args(argv)

    ligands = {r["enzyme"]: r for r in read_tsv(os.path.join(HERE, "af_ligands.tsv"))}
    sites = {r["enzyme"]: r for r in read_tsv(os.path.join(HERE, "start_sites.tsv"))}
    if args.tsv:
        print("\t".join(("enzyme", "source", "file", "chains", "cargo",
                         "keep_also", "trim", "start_mode", "start_residues",
                         "note")))
    problems = []
    for r in read_decisions(os.path.join(HERE, "structures.tsv")):
        name = r["enzyme"]
        path = structure_for(r)
        if not path or not os.path.exists(path):
            problems.append("%s: no structure at %s" % (name, path))
            continue
        st = read_any(path)
        # Only chains long enough to BE a protomer count as one. AlphaFill can
        # write a transplanted cofactor as ATOM records in a chain of its own:
        # 7OMT's model has a single SAH residue as "chain CF", which makes a
        # one-protomer file look like a dimer and turns a chain decision into
        # a question about a ligand.
        per_chain = {}
        for x in st.atoms:
            if not x.hetatm:
                per_chain.setdefault(x.chain, set()).add(x.resseq)
        present = sorted(c for c, rs in per_chain.items()
                         if len(rs) >= MIN_PROTOMER)
        stubs = sorted(c for c, rs in per_chain.items()
                       if len(rs) < MIN_PROTOMER)
        het = {x.resname.upper() for x in st.atoms if x.hetatm}

        chains, chain_why = CHAINS.get(name, (None, ""))
        if chains is None:
            if len(present) != 1:
                problems.append(
                    "%s: %d polymer chains (%s) and no decision recorded -- "
                    "which protomer is a choice, not a default"
                    % (name, len(present), "+".join(present)))
                continue
            chains, chain_why = [present[0]], "single-chain file"
            if stubs:
                chain_why += ("; ignored chain(s) %s, under %d residues and so "
                              "not a protomer" % ("+".join(stubs), MIN_PROTOMER))
        missing = [c for c in chains if c not in present]
        if missing:
            problems.append("%s: chain(s) %s are not in %s (it has %s)"
                            % (name, "+".join(missing), os.path.basename(path),
                               "+".join(present)))
            continue

        cargo, cargo_why = "", ""
        if name in EXTRA_CARGO:
            cargo, cargo_why = EXTRA_CARGO[name]
        elif name in ligands:
            cargo = ligands[name].get("channel_ligand", "")
            cargo_why = "%s (%s)" % (ligands[name].get("reason", "")[:60],
                                     ligands[name].get("exactness", ""))
        if cargo and cargo not in het:
            problems.append("%s: cargo %s is not in %s"
                            % (name, cargo, os.path.basename(path)))
            continue

        also, also_why = KEEP_ALSO.get(name, ([], ""))
        also = [x for x in also if x in het]
        trim, trim_why = TRIM.get(name, ([], ""))

        start_mode, start_res = "", ""
        if cargo:
            mode = (ligands.get(name, {}) or {}).get("start_mode", "")
            if mode == "unknown":
                # `pick_ligands.py --no-evidence` could not fetch the clash
                # count, so nothing decided whether this transplant can BE a
                # start point. Guessing either way is how a channel comes to
                # be measured from inside a side chain.
                problems.append(
                    "%s: start_mode is `unknown` -- af_ligands.tsv was "
                    "written without the AlphaFill clash evidence. Re-run "
                    "./pick_ligands.py --tsv with network access." % name)
                continue
            start_mode = mode or "ligand"
            if start_mode == "cavity":
                # A transplant too clashing to BE the start point still marks
                # the right residues: its centroid may sit inside a side chain
                # while the pocket it was placed into is the real one. So the
                # lining comes from its own contacts, in this same file, and
                # the start point is derived from those instead.
                start_res = ",".join(
                    str(r) for r in lining_residues(st, chains, cargo))
                if not start_res:
                    problems.append(
                        "%s: start_mode is cavity but %s has no residues "
                        "within %.1f A of it" % (name, cargo, CONTACT))
                    continue
        else:
            start_mode = "cavity"
            s = sites.get(name)
            if s and s.get("residues"):
                start_res = s["residues"]
            else:
                problems.append(
                    "%s: apo and no start site derived -- run "
                    "./derive_start_sites.py --tsv first" % name)
                continue

        note = "; ".join(x for x in (chain_why, cargo_why, also_why, trim_why) if x)
        if args.tsv:
            print("\t".join((
                name, r["structure_source"], os.path.relpath(path, HERE),
                "+".join(chains), cargo, "+".join(also),
                ",".join("%d-%d" % t for t in trim), start_mode, start_res,
                " ".join(note.split()))))
        else:
            print("%-7s %-9s %-24s chains %-5s %s"
                  % (name, r["structure_source"], os.path.basename(path)[:24],
                     "+".join(chains),
                     ("cargo %s" % cargo) if cargo
                     else "start-cavity from %d residue(s)"
                          % len(start_res.split(",")) if start_res else "?"))
            if also or trim:
                print("        %s%s" % ("keep also %s  " % "+".join(also) if also else "",
                                        "trim %s" % ",".join("%d-%d" % t for t in trim)
                                        if trim else ""))
    if problems:
        print("\n%d PROBLEM(S):" % len(problems), file=sys.stderr)
        for p in problems:
            print("  %s" % p, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
