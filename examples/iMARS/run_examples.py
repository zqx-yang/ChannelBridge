#!/usr/bin/env python3
"""Build and run the iMARS example set: stages 0 to 7, both tunnel branches.

The iMARS paper (Wang et al., Cell 188:1349, 2025) reports 36 enzyme fusion
constructs. Collapsed on the enzyme PAIR -- a linker variant is not a
different pair, and ChannelBridge does not fuse anything -- they are the 19
pairs in ``pairs.tsv``, over the 24 enzymes in ``enzymes.tsv``.

What this is for: exercising the pipeline on a set nobody curated for it,
and showing what input preparation decides. Every enzyme is now represented
by coordinates of its own protein, which leaves **2 of the 24** on an entry
ChannelsDB has annotated (AS and TS); the other 22 take the CAVER branch.
Two pairs are therefore mixed, with one branch on each side.

    ./run_examples.py fetch            # the site-transfer donors + annotated entries
    ./run_examples.py prep             # one chain per enzyme
    ./run_examples.py pairs            # write the 19 pair.yaml files
    ./run_examples.py tunnels          # stage 0 or 2a, then 2b
    ./run_examples.py report
    ./run_examples.py all --only iMARS_01_4CL-STS

Nothing here is bundled: the structures are deposited entries, fetched on
demand, and the directories they land in are not in version control. What IS
committed is the two tables, the pair files this writes, and this script --
which is to say, every decision, and none of the data.

CAVER is needed for all 38 enzyme sides -- the two annotated ones end up on
a seeded CAVER run as well -- so point $CAVER at it. Network access is needed
for `find_structures.py --fetch`, `fetch` and the ChannelsDB reports.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # the repository
SRC = os.path.join(ROOT, "src")
ENTRIES = os.path.join(HERE, "_entries")
RCSB = "https://files.rcsb.org/download/%s.pdb"

#: Which chain id each side gets in the PREPARED file. Most of these entries
#: have only chain A, so without renaming both enzymes of a pair would be
#: chain A and `chbridge validate` rejects the pair -- the pose stage needs
#: them distinguishable. Renaming is not moving: the ChannelsDB branch is
#: unaffected, because a database sphere chain carries coordinates and no
#: chain id.
CHAIN_OF = {"upstream": "A", "downstream": "B"}

#: How close a fetched chain must come to the prepared protein to be a tunnel
#: OF it. A report annotates the whole deposited assembly and does not say
#: which protomer each tunnel belongs to, so a chain can be perfectly valid
#: and 32 A from the chain you kept.
ATTACHED_MAX = 5.0

if os.path.isdir(os.path.join(SRC, "channelbridge")) and SRC not in sys.path:
    sys.path.insert(0, SRC)


# --------------------------------------------------------------- the tables

def _rows(name):
    """Read one of the TSV tables, dropping the commentary."""
    path = os.path.join(HERE, name)
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


#: Chain ids handed to each side of a pair, in order. Two enzymes of one
#: complex must not share a chain id -- `write_complex` refuses it -- and a
#: dimeric side needs more than one, so each side gets its own alphabet.
CHAIN_IDS = {"upstream": "ACEG", "downstream": "BDFH"}


def inputs():
    """`inputs.tsv`, the merged table input prep produces."""
    return {r["enzyme"]: r for r in _rows("inputs.tsv")}


def enzymes():
    """The enzyme table, with `inputs.tsv` layered over it.

    `enzymes.tsv` still carries what the FIRST version of this set used, which
    is kept deliberately as the record of what changed. Everything about which
    coordinates to measure now comes from `inputs.tsv`, so the two are merged
    here rather than in each caller: a caller that read `enzymes.tsv` directly
    would silently measure the superseded structure.
    """
    base = {r["enzyme"]: dict(r) for r in _rows("enzymes.tsv")}
    merged = inputs()
    for name, row in base.items():
        got = merged.get(name)
        if got is None:
            continue
        row["source_file"] = os.path.join(HERE, got["file"])
        # What is MEASURED, which for 19 of the 24 is not `entry`. `entry`
        # is the homolog the structure search first settled on and is kept
        # as the record of that search; a later pass replaced most of them
        # with models. Anything that reports provenance -- report.tsv, the
        # pair files -- has to read these two fields and not that one, or it
        # credits the measurement to coordinates nobody opened. It did: the
        # report named 5BSM for a tunnel computed on af3_4CL.cif.
        row["structure"] = os.path.basename(got["file"])
        row["structure_source"] = got["source"]
        row["chains"] = [c for c in got["chains"].split("+") if c]
        row["channel_ligand"] = got["cargo"]
        row["keep_also"] = [c for c in got["keep_also"].split("+") if c]
        row["trim"] = got["trim"]
        row["start_mode"] = got["start_mode"]
        row["start_residues"] = got["start_residues"]
        # The spec cavity_start_point takes, in the PREPARED file's numbering.
        row["start_spec"] = got["start_residues"]
        # ChannelsDB coverage is per deposited ENTRY. The annotation was
        # fetched for `entry`; it describes this enzyme's structure only if
        # that is still the entry being measured. For 19 of the 21 enzymes
        # that used the branch it is not -- they are now models, or a
        # different entry -- so the branch is withdrawn rather than left
        # pointing at tunnels of a protein the set no longer uses.
        same_entry = (got["source"] in ("pdb", "pdb-other")
                      and str(row.get("entry", "")).upper()
                      == os.path.splitext(row["structure"])[0].upper())
        #: Set when `entry` is no longer the structure being measured, so
        #: every reader can say which it is looking at instead of guessing.
        row["superseded_entry"] = "" if same_entry else row.get("entry", "")
        row["branch_changed"] = ""
        if row.get("branch") == "channelsdb" and not same_entry:
            row["branch"] = "caver"
            row["branch_changed"] = (
                "ChannelsDB annotated %s, which is no longer this enzyme's "
                "structure" % row.get("entry"))
    return base


def _provenance(e):
    """One clause naming the coordinates a result was computed on."""
    out = "%s (%s)" % (e.get("structure") or "?", e.get("structure_source") or "?")
    if e.get("superseded_entry"):
        out += ", superseding the %s of the first structure search" % (
            e["superseded_entry"])
    return out


def _chains_for(e, side):
    """``{old chain: new chain}`` for this side of a pair."""
    ids = CHAIN_IDS[side]
    kept = e.get("chains")
    if not kept:
        raise SystemExit(
            "%s has no `chains` in inputs.tsv. There is a `chain` column in "
            "enzymes.tsv and it is the FIRST version's answer, so it is not "
            "used as a fallback: run ./build_inputs.py --tsv > inputs.tsv."
            % e.get("enzyme", "?"))
    if len(kept) > len(ids):
        raise SystemExit("%s keeps %d chains; only %d ids are reserved for the "
                         "%s side" % (e["enzyme"], len(kept), len(ids), side))
    return {old: ids[i] for i, old in enumerate(kept)}


def pairs():
    return _rows("pairs.tsv")


def _selected(args):
    got = pairs()
    if args.only:
        want = set(args.only)
        got = [p for p in got if p["prefix"] in want]
        missing = want - {p["prefix"] for p in got}
        if missing:
            raise SystemExit("no such example(s): %s" % ", ".join(sorted(missing)))
    return got


def _enzyme_dir(prefix, what):
    d = os.path.join(HERE, prefix, what)
    os.makedirs(d, exist_ok=True)
    return d


def _prepped(prefix, enz):
    return os.path.join(HERE, prefix, "input", "%s.pdb" % enz["name"])


def _run(cmd, **kw):
    """A chbridge call, as a subprocess, so the example uses the real CLI."""
    env = dict(os.environ)
    env["PYTHONPATH"] = SRC + os.pathsep + env.get("PYTHONPATH", "")
    full = [sys.executable, "-m", "channelbridge.cli"] + [str(c) for c in cmd]
    return subprocess.run(full, env=env, capture_output=True, text=True, **kw)


def _say(proc, indent="      "):
    for line in (proc.stdout or "").splitlines() + (proc.stderr or "").splitlines():
        if line.strip():
            print(indent + line)


# ------------------------------------------------------------------- stages

def cmd_fetch(args):
    """The deposited entries this set still reads, by name.

    Two kinds, and neither is "every entry in `enzymes.tsv`" -- which is
    what this used to fetch. That column names the homologs the FIRST
    structure search settled on, so it downloaded twenty files of which
    twelve were read by nothing, and the six that mattered arrived by
    accident: edit the superseded table and `derive_start_sites.py` would
    silently lose its donors. A fourth reader of data nothing should be
    reading.

    * the holo homologs `derive_start_sites.py` transfers a site from, which
      that module names in `HOMOLOG`;
    * the entries whose ChannelsDB annotation is still fetched, which is the
      two enzymes left on that branch.

    Everything else each enzyme is measured ON comes from
    `find_structures.py --fetch` (deposited entries and models) or from an
    AlphaFold3 run, per `inputs.tsv`.
    """
    import urllib.request

    sys.path.insert(0, HERE)
    from derive_start_sites import HOMOLOG

    os.makedirs(ENTRIES, exist_ok=True)
    donors = {pid for pid, _ch, _lig in HOMOLOG.values()}
    annotated = {e["entry"] for e in enzymes().values()
                 if e.get("branch") == "channelsdb" and e.get("entry")}
    want = sorted(donors | annotated)
    print("  %d site-transfer donor(s): %s" % (len(donors), " ".join(sorted(donors))))
    print("  %d still-annotated entr(y/ies): %s"
          % (len(annotated), " ".join(sorted(annotated)) or "-"))
    for pid in want:
        dst = os.path.join(ENTRIES, "%s.pdb" % pid)
        if os.path.exists(dst) and os.path.getsize(dst) > 0 and not args.force:
            print("  have %s" % pid)
            continue
        print("  fetch %s" % pid, flush=True)
        with urllib.request.urlopen(RCSB % pid, timeout=60) as fh:
            data = fh.read()
        if len(data) < 1000:
            raise SystemExit("%s came back as %d bytes" % (pid, len(data)))
        with open(dst, "wb") as out:
            out.write(data)
    print("%d entries in %s" % (len(want), ENTRIES))
    return 0


def cmd_prep(args):
    """Stage 1, per enzyme per example: keep one chain and the ligands.

    The same enzyme is prepped separately under each example that uses it.
    That is deliberate: an example should be deletable on its own, and the
    duplication is a few hundred kilobytes of a file that is not committed.
    """
    enz = enzymes()
    bad = 0
    for pair in _selected(args):
        print("%s" % pair["prefix"])
        for side in ("upstream", "downstream"):
            e = enz[pair[side]]
            src = e.get("source_file")
            if not src or not os.path.exists(src):
                raise SystemExit(
                    "%s: no %s. Run ./find_structures.py --fetch and "
                    "./build_inputs.py --tsv first" % (e["enzyme"], src))
            out = _prepped(pair["prefix"], e)
            _enzyme_dir(pair["prefix"], "input")
            mapping = _chains_for(e, side)
            cmd = ["prep", src,
                   "--chains", ",".join(mapping),
                   "--out", out,
                   "--relabel", ",".join("%s:%s" % kv for kv in mapping.items()),
                   "--ligand-copies", "best"]
            keep = [v for v in (e["channel_ligand"], *(e.get("keep_also") or [])) if v]
            if keep:
                cmd += ["--keep", ",".join(keep)]
            if e.get("trim"):
                cmd += ["--drop-resi", e["trim"]]
            proc = _run(cmd)
            if proc.returncode != 0:
                bad += 1
                print("   !! %-6s %s FAILED" % (e["enzyme"],
                                                os.path.basename(src)))
                _say(proc)
                continue
            print("   %-6s %-26s chain %s -> %s (chain %s)%s" % (
                e["enzyme"], os.path.basename(src), "+".join(mapping),
                os.path.basename(out), "+".join(mapping.values()),
                "  trimmed %s" % e["trim"] if e.get("trim") else ""))
    return 1 if bad else 0


def _site_selection(e):
    """What `site`/`vent` say in the pair file, and why they say it.

    With a ligand the site is the ligand. Without one it is the catalytic
    residues the start point was derived from -- which is the honest
    statement: nobody has told us where the reactive atom is, and a guess
    written as a selection would look like knowledge.
    """
    if e["channel_ligand"]:
        return "resn %s" % e["channel_ligand"]
    if e["start_spec"]:
        nums = "+".join(s.split(":")[-1] for s in e["start_spec"].split(","))
        return "resi %s" % nums
    return "polymer"


def cmd_pairs(args):
    """Write one pair file per example."""
    enz = enzymes()
    for pair in _selected(args):
        prefix = pair["prefix"]
        os.makedirs(os.path.join(HERE, prefix), exist_ok=True)
        up, dn = enz[pair["upstream"]], enz[pair["downstream"]]
        text = _pair_yaml(pair, up, dn)
        path = os.path.join(HERE, prefix, "pair.yaml")
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        print("  wrote %s/pair.yaml  (%s / %s, %s)"
              % (prefix, up["enzyme"], dn["enzyme"], _branches(up, dn)))
    return 0


def _branches(up, dn):
    """Which branch each side of a pair actually takes.

    Computed from the merged table, not read from `pairs.tsv`'s `branches`
    column. That column was written when the set first ran and stopped being
    true when `inputs.tsv` replaced 19 of the 24 structures: it still called
    GGPPS-TS a `channelsdb` pair although GGPPS had become an AF3 model with
    no annotation to fetch. A label that can go stale is not kept beside the
    thing that determines it.
    """
    kinds = {up["branch"], dn["branch"]}
    if len(kinds) > 1:
        return "mixed"
    return kinds.pop()


def _pair_yaml(pair, up, dn):
    out = [
        "# ChannelBridge pair file -- %s (upstream) to %s (downstream)."
        % (up["enzyme"], dn["enzyme"]),
        "#",
        "# GENERATED by run_examples.py from inputs.tsv, enzymes.tsv and",
        "# pairs.tsv; edit those, not this. One of the %s pairs of the iMARS"
        % len(pairs()),
        "# fusion set.",
        "#",
        "# Neither structure is the iMARS construct. The paper's enzymes are",
        "# plant and yeast proteins and almost none has a structure of its",
        "# own, so each side here is a homologous deposited entry or a model:",
        "#   upstream   %s" % _provenance(up),
        "#   downstream %s" % _provenance(dn),
        "# enzymes.tsv records the search that found the homolog and its",
        "# score, which is how much the homology is worth.",
        "#",
        "# Tunnel source: %s." % _branch_sentence(up, dn),
        "",
        "name: %s" % pair["prefix"],
        "",
        "# What travels between them: the upstream enzyme's product.",
        "channelled: %s" % (up["product"] or "unknown"),
        "",
        "# O_IC is measured at the UPSTREAM exit only, so this is not cosmetic.",
        "upstream: 1",
        "",
        "targets:",
    ]
    for tid, e, side in ((1, up, "upstream"), (2, dn, "downstream")):
        out += _target_yaml(tid, e, CHAIN_OF[side])
    out += [
        "",
        "# Normalisation is what makes the geometry independent of the CAVER",
        "# build or the database that produced the chain. See",
        "# docs/caver-builds.md for the measurement behind each number.",
        "tunnel:",
        "  normalise: true",
        "  spacing: 1.0",
        "  burial_min: 15",
        "  burial_radius: 8.0",
        "  persist: 2.0",
        "  oic_arm: 1.0",
        "",
        "# The published DO coefficients. Nothing here is tuned per pair.",
        "do_coefficients: imars-fig2c",
        "",
        "# Left at the package defaults on purpose: this set is for exercising",
        "# the tunnel stages, and a per-pair search grid would be a claim",
        "# about these enzymes that nobody here has earned.",
        "",
    ]
    return "\n".join(out) + "\n"


def _branch_sentence(up, dn):
    if up["branch"] == dn["branch"]:
        return ("ChannelsDB for both" if up["branch"] == "channelsdb"
                else "CAVER for both")
    fetched = up if up["branch"] == "channelsdb" else dn
    computed = dn if fetched is up else up
    return ("MIXED -- %s is annotated in ChannelsDB and its tunnels are "
            "fetched, %s is not and CAVER is run on it"
            % (fetched["enzyme"], computed["enzyme"]))


def _target_yaml(tid, e, chain):
    out = [
        "  - id: %d" % tid,
        "    name: %s                 # %s" % (e["name"], e["enzyme"]),
        "    file: input/%s.pdb" % e["name"],
        "    chain: %s                       # %s chain %s, renamed at prep"
        % (chain, e["structure"], "+".join(e["chains"])),
    ]
    if e.get("superseded_entry"):
        out.append("    # Not %s: that is the homolog the first structure"
                   % e["superseded_entry"])
        out.append("    # search settled on, and this file replaced it.")
    if e["channel_ligand"]:
        out.append("    ligands: [%s]" % e["channel_ligand"])
        out.append("    channel_ligand: %s" % e["channel_ligand"])
    sel = _site_selection(e)
    out.append('    site: "%s"' % sel)
    out.append('    vent: "%s"' % sel)
    if not e["channel_ligand"]:
        out.append("    # No ligand marks this site, so site/vent name the")
        out.append("    # catalytic residues instead. Narrow them to the")
        out.append("    # reactive atom if you know which it is.")
    if e["branch"] == "caver":
        lig = (e["channel_ligand"] or "").strip().upper()
        out.append("    caver:")
        # The cargo belongs in `exclude` whenever there is one, whether or
        # not it is also the start point: left in, it sits in the channel
        # whose bottleneck is being measured. This file listed only the
        # waters, which was right for the enzymes that start from a ligand
        # (the runner excluded it there) and wrong for the ones that keep a
        # cargo and start from lining residues.
        out.append("      exclude: [%s]" % ", ".join(_caver_exclude(lig)))
        out.append("      probe_radius: 1.4")
        out.append("      shell_radius: 10")
        out.append("      shell_depth: 6")
        out.append("      clustering_threshold: 3.5")
        out.append("      seed: 1")
        if e.get("start_mode") == "ligand" and lig:
            out.append("      # Start point: the centroid of ONE %s copy, in"
                       % lig)
            out.append("      # the chain whose channel this is. Not the mean")
            out.append("      # over copies -- that lands between two sites and")
            out.append("      # in neither. See run_examples._cargo_seed.")
        else:
            out.append("      # Derived, not written by hand: the centroid of the")
            out.append("      # residues below sits inside an atom, and CAVER then")
            out.append("      # finds zero tunnels without saying why. See")
            out.append("      # channelbridge.caver.cavity_start_point.")
            out.append("      #   chbridge caver config --structure input/%s.pdb \\"
                       % e["name"])
            out.append("      #       --start-cavity %s" % e["start_spec"])
    else:
        # The guard above is `branch`, and `branch` is withdrawn whenever
        # the structure changed -- so `entry` IS what is being measured
        # here. Asserted rather than inherited: this is the one place a
        # superseded `entry` would still be printed as provenance.
        assert not e.get("superseded_entry"), (
            "%s is on the ChannelsDB branch but its entry %s is superseded "
            "by %s" % (e["enzyme"], e.get("entry"), e.get("structure")))
        out.append("    # ChannelsDB branch: the tunnels come from entry %s"
                   % e["entry"])
        out.append("    # (%s tunnels, %s). Nothing is computed, and nothing"
                   % (e["cdb_tunnels"], e["cdb_category"]))
        out.append("    # may move the coordinates -- their chain is in this")
        out.append("    # entry's frame. Where the fetched chain stops before")
        out.append("    # the surface, run_examples re-runs CAVER SEEDED from")
        out.append("    # its buried end; report.tsv's `resolved` column says")
        out.append("    # which of the two produced each row.")
    if e["note"]:
        for chunk in _wrap(e["note"], 66):
            out.append("    # %s" % chunk)
    return out


def _wrap(text, width):
    words, line, out = text.split(), "", []
    for w in words:
        if line and len(line) + 1 + len(w) > width:
            out.append(line)
            line = w
        else:
            line = (line + " " + w).strip()
    if line:
        out.append(line)
    return out


def cmd_tunnels(args):
    """Stage 0 or 2a per enzyme, then 2b (normalise) for both."""
    enz = enzymes()
    results = []
    for pair in _selected(args):
        print("%s  (%s)" % (pair["prefix"],
                            _branches(enz[pair["upstream"]],
                                      enz[pair["downstream"]])))
        for side in ("upstream", "downstream"):
            e = enz[pair[side]]
            try:
                row = (_channelsdb_side if e["branch"] == "channelsdb"
                       else _caver_side)(pair["prefix"], e, side, args)
            except Exception as exc:
                row = dict(ok=False, detail="%s: %s" % (type(exc).__name__, exc))
            row.update(prefix=pair["prefix"], side=side, enzyme=e["enzyme"],
                       branch=e["branch"],
                       branch_changed=e.get("branch_changed", ""),
                       structure=e.get("structure", ""),
                       source=e.get("structure_source", ""),
                       superseded_entry=e.get("superseded_entry", ""))
            results.append(row)
            print("   %-9s %-6s %-24s %s" % (
                e["branch"], e["enzyme"], e.get("structure", "?"),
                row.get("detail", "ok")))
    _write_report(results)
    return 0 if all(r.get("ok") for r in results) else 1


def _channelsdb_side(prefix, e, side, args):
    """Stage 0: fetch the entry's tunnels and keep the best of its category.

    Four filters, in this order, and each is here because the obvious rule
    without it picked something wrong on this set:

    1. **At least 5 A long.** Vavra et al. flag a first tunnel shorter than
       that as a sign the pocket is simply on the surface; it is already
       ``channelbridge.caver.SHORT_TUNNEL`` on the CAVER branch. ChannelsDB
       has two-sphere, 1.0 A "tunnels" -- 5UV1 has nothing else in its
       Procognate category -- and taking one is not a conservative choice,
       it is a wrong one.
    2. **It must be attached to the chain we prepared.** A report annotates
       the tunnels of the whole deposited assembly, and which protomer they
       belong to is not stated. Prepping chain A of 3KRF and taking its
       top Procognate tunnel gave a chain sitting **32 A from the protein at
       both ends** -- correct coordinates, in a frame whose protein had been
       deleted. It passed every other filter: a detached chain trivially
       "leaves the protein". Four entries here annotate a protomer other
       than A (`enzymes.tsv` records which), and this filter is what makes
       a wrong chain a failure instead of a plausible number.
    3. **The chain must leave the protein.** Everything downstream measures
       an exit and a mouth direction, so a chain that stops while still
       buried has no exit to find and its last sphere is not one. MEASURED,
       with the same burial test the normalise stage uses, rather than
       inferred from which program produced the chain.
    4. **Then the widest bottleneck**, ties broken by the shorter chain.
       Not the longest: ChannelsDB reports routes up to 84 A through these
       entries, and the longest thing leaving an active site is whatever
       path the algorithm could thread furthest, not the channel a
       metabolite uses.

    Filter 2 is where the two methods separate, which is worth seeing
    rather than assuming, so the per-method counts are reported: the
    CAVER-derived chains are truncated at the database's own surface
    definition and largely stop inside the protein, while the MOLE chains
    run out into solvent.

    And be clear what none of this is. CAVER's own ranking is by *priority*
    (throughput), and nothing in a ChannelsDB report carries it, so the two
    branches do not choose on the same criterion. A mixed pair has one
    tunnel chosen by throughput and one by bottleneck -- a difference
    between the branches, not a property of the enzymes.
    """
    import numpy as np

    from channelbridge.caver import SHORT_TUNNEL
    from channelbridge.channelsdb import fetch, to_tunnel
    from channelbridge.channelsdb import tunnels as db_tunnels
    from channelbridge.io.pdb import read_pdb
    from channelbridge.tunnel import BURIAL_SURFACE_MIN, burial_at

    out = _enzyme_dir(prefix, "output")
    cache = os.path.join(HERE, "_entries")
    report = fetch(e["entry"], cache_dir=cache)
    got = db_tunnels(report)
    if not got:
        return dict(ok=False, resolved="channelsdb",
                    detail="no tunnels in the report")
    # CATEGORIES is already in the order this pipeline prefers, so the first
    # tunnel's category is the best one on offer.
    best_cat = got[0].category
    same = [t for t in got if t.category == best_cat]
    methods = sorted({t.method for t in same})

    long_enough = [t for t in same if t.length >= SHORT_TUNNEL]
    if not long_enough:
        return dict(ok=False, resolved="channelsdb",
                    detail="%d tunnels (%d in %s, from %s), none reaching "
                           "%.0f A; the longest is %.1f A, which is a surface "
                           "pocket rather than a channel"
                           % (len(got), len(same), best_cat, "+".join(methods),
                              SHORT_TUNNEL, max(t.length for t in same)))

    prot = np.array([a.xyz for a in read_pdb(_prepped(prefix, e)).atoms
                     if not a.hetatm], dtype=float)

    attached = [
        t for t in long_enough
        if float(np.linalg.norm(np.asarray(t.xyz, dtype=float)[:, None, :]
                                - prot[None, :, :], axis=2).min()) <= ATTACHED_MAX
    ]
    if not attached:
        return dict(ok=False, resolved="wrong-chain",
                    detail="%d tunnels (%d in %s, %d over %.0f A) and NONE "
                           "comes within %.0f A of chain %s. They belong to "
                           "another protomer of %s: prepare the chain the "
                           "annotation is for (see the chain column of "
                           "enzymes.tsv)"
                           % (len(got), len(same), best_cat, len(long_enough),
                              SHORT_TUNNEL, ATTACHED_MAX, CHAIN_OF[side],
                              e["entry"]))

    emerges, by_method = [], {}
    for t in attached:
        by_method.setdefault(t.method, [0, 0])
        by_method[t.method][1] += 1
        if float(burial_at(t.xyz, prot)[[0, -1]].min()) < BURIAL_SURFACE_MIN:
            by_method[t.method][0] += 1
            emerges.append(t)
    shape = ", ".join("%s %d/%d leave the protein" % (m, n, d)
                      for m, (n, d) in sorted(by_method.items()))
    if not emerges:
        # No chain reaches the surface, so none of them HAS an exit, and
        # D_IC and O_IC are both measured at one. Do not measure it: hand
        # the route over to CAVER, seeded from where the database put the
        # site. See _seeded_caver_side.
        pick = max(attached, key=lambda t: (round(t.bottleneck, 2), -t.length))
        row = _seeded_caver_side(prefix, e, side, pick, args)
        row["detail"] = ("%d tunnels (%d in %s, %d attached, %d over %.0f A); "
                         "none reaches the surface (%s), so CAVER was run "
                         "instead -- %s"
                         % (len(got), len(same), best_cat, len(attached),
                            len(long_enough), SHORT_TUNNEL, shape,
                            row["detail"]))
        row["picked"] = pick.label
        row["n"] = len(got)
        return row

    pick = max(emerges, key=lambda t: (round(t.bottleneck, 2), -t.length))
    tun = to_tunnel(pick)
    dst = os.path.join(out, "%s_tunnel.pdb" % e["name"])
    tun.write_caver_pdb(dst)
    norm = _normalise(prefix, e, dst, args)
    return dict(ok=norm["ok"], n=len(got), picked=pick.label,
                resolved="channelsdb",
                detail="%d tunnels (%d in %s, %d attached); %d over %.0f A, "
                       "of which %s; took %s, %.1f A, bottleneck %.2f A; %s"
                       % (len(got), len(same), best_cat, len(attached),
                          len(long_enough), SHORT_TUNNEL, shape, pick.label,
                          pick.length, pick.bottleneck, norm["detail"]))


def _seeded_caver_side(prefix, e, side, db_tunnel, args):
    """Run CAVER on the route ChannelsDB found, because the fetched chain
    stops before the surface.

    This is what the branch should do with an unusable chain, and it is not
    a workaround. The two sources are good at different things. ChannelsDB
    knows *what a tunnel is for* -- its categories say whether a route leads
    to the cognate ligand, a cofactor or the catalytic site, which CAVER
    cannot tell you and which is the hard part of choosing a tunnel. What it
    does not give is a chain that reaches bulk solvent, because it truncates
    by its own surface definition. CAVER gives exactly the converse: run it
    with a large enough shell radius and the chain leaves the protein, but
    nothing in its output says which of its clusters is the substrate route.

    So: take the start point from the database and the geometry from CAVER.
    The database chain's buried end sits in the site it was annotated for,
    which is the piece that was otherwise guesswork -- on the CAVER branch
    proper it has to come from a ligand or from transferred catalytic
    residues. Then CHECK that CAVER found the same route, with
    ``channelsdb.compare``; no superposition is needed because prep never
    moves coordinates, so both are already in this entry's frame.

    Judge that check by ``Match.recovery()``, not ``Match.verdict()``. Ours
    is longer by design, so the exits cannot coincide and ``verdict``'s
    coverage term counts our extra spheres as misses: read the wrong way
    round it called six recovered routes "different tunnel", one of them an
    exact match over 100 % of the database chain. ``recovery()`` asks the
    question that is actually being asked -- how much of the chain they
    published lies on the chain we computed.

    "different route" is still a real outcome, and on this set it happens.
    When it does, the annotation has not carried over and the tunnel is
    CAVER's own choice with no category attached to it, which is worth
    knowing rather than smoothing over.
    """
    import numpy as np

    from channelbridge.caver import (
        CaverParams,
        read_summary,
        select,
        write_config,
    )
    from channelbridge.channelsdb import compare
    from channelbridge.io.pdb import read_pdb
    from channelbridge.tunnel import Tunnel, burial_at

    struct = _prepped(prefix, e)
    out = _enzyme_dir(prefix, "output")
    atoms = read_pdb(struct)
    prot = np.array([a.xyz for a in atoms.atoms if not a.hetatm], dtype=float)

    # Seed order matters, and the database chain comes FIRST. Its buried end
    # is where the annotated route starts, so a CAVER run from there is "the
    # route they identified, computed far enough to reach solvent" -- which
    # is the whole reason to come here rather than run CAVER plainly. The
    # ligand centroid is a better-behaved start point but a different
    # question: seeding from it instead took the number of enzymes whose
    # CAVER route matched the annotation from 6 down to 2, because CAVER
    # then finds whatever leaves the ligand site, annotation or not.
    #
    # So the ligand is a FALLBACK, for when the chain end yields nothing.
    # It rescues taxadiene synthase (3P5P) and squalene synthase (3WEF),
    # which gave zero clusters from the chain end at probe 1.4, 1.1 and 0.9
    # and three each from their bound ligand.
    # The cargo is excluded on BOTH seeds, not only on the one that starts
    # from it: a cargo left in the structure sits in the channel whose width
    # is being measured whichever point the search starts from.
    lig = (e["channel_ligand"] or "").upper()
    exclude = _caver_exclude(lig)
    b = burial_at(db_tunnel.xyz, prot)
    seeds = [(db_tunnel.xyz[0] if b[0] >= b[-1] else db_tunnel.xyz[-1],
              "the buried end of %s" % db_tunnel.label)]
    # ONE copy of the cargo. Not the mean over all of them -- see
    # `_cargo_seed`, which is where this was fixed on the other branch and
    # where it had to be fixed again here. Taxadiene synthase is the case:
    # 3P5P holds two A3C copies in one chain, 8.20 A apart, and the pooled
    # centroid of the two sat 4.10 A from either site. CAVER ran from it
    # without complaint and returned tunnels of the space between them.
    seed = _cargo_seed(atoms, lig, _first_chain(e, side)) if lig else None
    if seed is not None:
        seeds.append((seed[0], seed[1]))

    # One run directory per seed. Sharing one would make the fallback
    # uncacheable: the config is compared by content, the first seed's
    # config is written first, so a set whose answer came from the second
    # seed would look stale on every later pass and recompute both. Keeping
    # them apart also leaves the failed attempt on disk, which is the
    # evidence for "the chain end yields nothing here".
    base = os.path.join(out, "caver_%s" % e["name"])
    kept, how, run_dir, tried = [], None, None, []
    for i, (seed_xyz, label) in enumerate(seeds):
        tag = "" if i == 0 else "_seed%d" % i
        run_dir = base + tag
        params = CaverParams(probe_radius=1.4, shell_radius=10.0,
                             shell_depth=6.0, clustering_threshold=3.5,
                             seed=1, exclude=tuple(exclude))
        conf = write_config(
            os.path.join(out, "%s_caver%s.conf" % (e["name"], tag)),
            params, [float(v) for v in seed_xyz])
        tried.append(label)
        try:
            # Keyed on the structure AND this config, so a corrected
            # starting point is never answered out of the cache computed
            # from the old one.
            _run_caver(struct, conf, run_dir, force=args.force)
        except RuntimeError:
            # CAVER writes no summary at all when the point it was given is
            # outside the structure or inside an atom. That is a reason to
            # try the next seed, not to fail the enzyme -- unless there is
            # no next seed.
            if i + 1 < len(seeds):
                continue
            raise
        kept = select(read_summary(run_dir), keep=3, floor=0.55).kept
        how = label
        if kept:
            break
    if not kept:
        return dict(ok=False, resolved="seeded-caver-failed",
                    detail="CAVER found no tunnel above the priority floor "
                           "from %s -- the site is enclosed in this protomer. "
                           "Measured: nothing at probe 1.4, 1.1 or 0.9 "
                           "either. Neither source gives a usable exit here; "
                           "if the route runs along a subunit interface, keep "
                           "more chains" % " or ".join(tried))
    kept_atoms = _check_loaded(run_dir, struct, exclude)
    top = kept[0]
    cluster = _find_caver_cluster(run_dir, top.cluster_id)
    match = compare(Tunnel.from_caver_cluster(cluster), db_tunnel)
    norm = _normalise(prefix, e, cluster, args)
    return dict(ok=norm["ok"], resolved="channelsdb-seeded CAVER",
                recovery=match.recovery(), loaded=kept_atoms,
                detail="seeded from %s; %s; %d clusters kept, took cl%d "
                       "(priority %.3f); %.0f %% of the database chain lies "
                       "on ours (exit gap %.1f A, ours runs further by "
                       "design) -> %s; %s"
                       % (how, kept_atoms, len(kept), top.cluster_id,
                          top.priority, 100 * match.db_coverage,
                          match.exit_gap, match.recovery(), norm["detail"]))


def _start_spec_for(e, side):
    """``start_spec`` with the chain the PREPARED file actually uses.

    The table names residues in the deposited entry's chain; prep renames
    that chain so the two enzymes of a pair are distinguishable. Carrying
    the old id over meant CAVER was handed residues in a chain that no
    longer exists -- which surfaced as "none of the residues is in the
    structure", and would have been a silent miss had the ids happened to
    exist in the renamed chain too.
    """
    out = []
    want = (_chains_for(e, side) or {}).get(
        (e.get("chains") or ["A"])[0], CHAIN_OF[side])
    for raw in (e.get("start_spec") or "").split(","):
        raw = raw.strip()
        out.append("%s:%s" % (want, raw.split(":")[-1]) if raw else raw)
    return [r for r in out if r]


def _first_chain(e, side):
    """The PREPARED chain id of the protomer whose channel is measured."""
    return (_chains_for(e, side) or {}).get(
        (e.get("chains") or ["A"])[0], CHAIN_OF[side])


def _cargo_copies(atoms, lig):
    """``{(chain, resseq): [xyz, ...]}`` for every copy of *lig* in a file.

    Grouped by residue, because a pooled centroid over copies that face each
    other is not in either of them. The het flag is not required: AlphaFill
    writes a transplanted cofactor as ATOM records, so a cargo can be a
    perfectly ordinary-looking polymer record.
    """
    copies = {}
    for a in atoms.atoms:
        if (a.resname or "").upper() == lig:
            copies.setdefault((a.chain, a.resseq), []).append(a.xyz)
    return copies


def _cargo_seed(atoms, lig, first_chain):
    """One cargo copy's centroid as a start point, and which copy it is.

    ONE copy, never the mean over all of them, and this is the second place
    the same error was found. A kept dimer has a copy per subunit: BAS's two
    HC4 sit 29.8 A apart and their pooled centroid lands 1.40 A from the
    nearest protein atom -- inside an atom, between the subunits -- where
    CAVER finds nothing, which is how it was noticed the first time. The
    second time it was silent instead: taxadiene synthase's entry has two
    A3C copies 8.20 A apart in ONE chain, the pooled centroid sits 4.10 A
    from each of them, in neither site, and CAVER ran happily from it and
    returned tunnels of a place that is not the active site. A start point
    that is wrong by 4 A does not fail, it answers a different question.

    So the copy is chosen, by the chain whose channel is being measured and
    then by residue number. Where one chain holds two copies that choice is
    arbitrary and says so in its own sentence: there is nothing in the file
    that prefers one, and a tiebreak dressed up as a reason would be worse
    than an admitted coin toss.

    Returns ``(xyz, how)``, or ``None`` when the cargo is not in the file.
    """
    import numpy as np

    copies = _cargo_copies(atoms, lig)
    if not copies:
        return None
    want = sorted(copies, key=lambda k: (k[0] != first_chain, k[0], k[1]))[0]
    xyz = [float(v) for v in np.asarray(copies[want]).mean(axis=0)]
    same_chain = [k for k in copies if k[0] == want[0]]
    if len(same_chain) > 1:
        how = ("%s %d of %d copies in chain %s (%d in this chain; which one "
               "is arbitrary -- nothing in the file prefers either, and the "
               "lower residue number wins)"
               % (lig, want[1], len(copies), want[0], len(same_chain)))
    elif len(copies) > 1:
        how = ("the %s centroid in chain %s (%d copies present, one per "
               "protomer)" % (lig, want[0], len(copies)))
    else:
        how = "the %s centroid in chain %s" % (lig, want[0])
    return xyz, how


def _caver_exclude(lig):
    """Residue names CAVER must not load, for a side carrying cargo *lig*.

    The cargo is excluded whenever there IS one, not only when it also
    supplies the start point. The old rule was gated on ``start_mode ==
    "ligand"``, so the four enzymes that keep their cargo and start from
    lining residues instead (GES, MYS, OCS here, and taxadiene synthase)
    were measured with nothing but waters named -- which is to say, with the
    cargo asked to sit in the channel whose width was being measured.

    On this set it did no damage, and that is the point: the GST transplants
    carry their donor's occupancy of 0.75 and CAVER never loaded them, and
    the A3C was removed as a starting-point residue. Both are accidents of
    this data. An occupancy-1.00 transplant in cavity mode would have
    blocked its own channel and reported a bottleneck for it. The rule now
    matches the intent, and `_check_loaded` measures that it held.
    """
    exclude = ["HOH", "WAT"]
    if lig and lig not in exclude:
        exclude.append(lig)
    return exclude


def _check_loaded(run_dir, struct, exclude):
    """Did CAVER load what we meant to hand it? Returns a sentence.

    One-sided, and the bound is deliberately the only one that holds. CAVER
    cannot measure more than the file minus what we excluded by name, so a
    SURPLUS is an error: something meant to be kept out of the measurement
    is in it -- a misspelled residue name, a keyword CAVER ignored -- which
    is the failure this exists to catch and the one no other output shows.

    A DEFICIT is reported, not judged, because CAVER removes atoms of its
    own accord and does not say so. Measured here, on this set:

    * **Starting-point residues.** A residue holding or touching the start
      point is dropped; that is how taxadiene synthase's 58 A3C atoms left
      the measurement before anyone named them.
    * **A clashing transplant.** The three AlphaFill GST models lose all 19
      of their ligand atoms, and that was first mistaken for an occupancy
      rule, because they carry their donor's occupancy of 0.75. It is not:
      taxadiene synthase's 74 half-occupancy altloc atoms are all loaded.
      The GST models clash into the protein -- 0.56 A closest approach
      against A3C's 2.77 A -- and the removal is order-dependent, so moving
      the ligand to the top of the file brings 9 of its 19 atoms back and
      drops 10 protein atoms instead. Whatever the internal rule is, it is
      not one to predict from the file, so it is not predicted here.

    This is why the cargo is now excluded BY NAME rather than left to be
    absent: on these four enzymes it happened to be absent anyway, and an
    exclusion that depends on a clash being bad enough is not an exclusion.
    """
    from channelbridge.caver import read_log
    from channelbridge.io import read_pdb

    atoms = read_pdb(struct).atoms
    names = {n.upper() for n in exclude}
    ceiling = [a for a in atoms if (a.resname or "").upper() not in names]
    got = read_log(run_dir)
    cut = len(atoms) - len(ceiling)
    if got.atoms > len(ceiling):
        raise RuntimeError(
            "CAVER loaded %d atoms from %s but the file has only %d once the "
            "%d atoms of %s are taken out. Something that was meant to be "
            "kept out of the measurement is in it."
            % (got.atoms, os.path.basename(struct), len(ceiling), cut,
               " ".join(sorted(names))))
    short = len(ceiling) - got.atoms
    return ("loaded %d of %d atoms (%d excluded by name%s)"
            % (got.atoms, len(ceiling), cut,
               ", %d dropped by CAVER itself" % short if short else ""))


def _caver_side(prefix, e, side, args):
    """Stage 2a: derive a start point, run CAVER, select, then normalise."""
    from channelbridge.caver import (
        CaverParams,
        cavity_start_point,
        read_summary,
        select,
        write_config,
    )

    struct = _prepped(prefix, e)
    if not os.path.exists(struct):
        return dict(ok=False, resolved="caver",
                    detail="not prepped; run `prep` first")
    out = _enzyme_dir(prefix, "output")
    # Where the channel starts. With a cargo in the file it is where that
    # sits; without one, the point is derived from the lining residues --
    # their centre of gravity lands inside an atom for a catalytic triad,
    # which is why cavity_start_point exists. Either way the cargo is
    # excluded by name, so it cannot block the channel being measured.
    from channelbridge.io import read_pdb

    lig = (e.get("channel_ligand") or "").strip().upper()
    exclude = _caver_exclude(lig)
    if e.get("start_mode") == "ligand" and lig:
        seed = _cargo_seed(read_pdb(struct), lig, _first_chain(e, side))
        if seed is None:
            return dict(ok=False, resolved="caver",
                        detail="start_mode is ligand but %s is not in the "
                               "prepared file" % lig)
        start, how = seed
        offset = burial = None
    else:
        spec = _start_spec_for(e, side)
        if not spec:
            return dict(ok=False, resolved="caver",
                        detail="no cargo and no start residues; run "
                               "./build_inputs.py --tsv")
        found = cavity_start_point(struct, spec, probe_radius=1.4)
        start, offset, burial = found.xyz, found.offset, found.burial
        how = "%d lining residue(s)" % len(spec)
    params = CaverParams(probe_radius=1.4, shell_radius=10.0, shell_depth=6.0,
                         clustering_threshold=3.5, seed=1,
                         exclude=tuple(exclude))
    conf = write_config(os.path.join(out, "%s_caver.conf" % e["name"]),
                        params, start)
    run_dir = os.path.join(out, "caver_%s" % e["name"])
    _run_caver(struct, conf, run_dir, force=args.force)
    kept_atoms = _check_loaded(run_dir, struct, exclude)
    sel = select(read_summary(run_dir), keep=3, floor=0.55)
    if not sel.kept:
        return dict(ok=False, resolved="caver",
                    detail="CAVER found no tunnel above the floor")
    top = sel.kept[0]
    cluster = _find_caver_cluster(run_dir, top.cluster_id)
    norm = _normalise(prefix, e, cluster, args)
    where = ("start from %s, %.1f A off the residue centroid (burial %d)"
             % (how, offset, burial) if offset is not None
             else "start from %s" % how)
    return dict(ok=norm["ok"], n=len(sel.kept), resolved="caver",
                loaded=kept_atoms,
                detail="%s; %s; %d clusters kept, took cl%d (priority %.3f); %s"
                       % (where, kept_atoms, len(sel.kept), top.cluster_id,
                          top.priority, norm["detail"]))


def _is_stale(output, *inputs):
    """True if *output* must be rebuilt: missing, or older than any input.

    Existence alone is the wrong test, and it cost three stages in one run.
    Re-preparing the set on corrected structures left yesterday's CAVER
    summaries, patch tables and pose searches in place, so every stage
    reported results computed on the HOMOLOG coordinates as though they were
    the new ones -- the exact substitution the re-prep existed to undo,
    reintroduced by a cache. The patch stage was worse: it judged SUCCESS by
    the same test, so a command that had just failed looked like one whose
    work was already done.

    A cache keyed on "is there output" cannot see that its input changed.
    Keyed on timestamps it can, and ONE helper serves every stage so that the
    next stage added cannot quietly get this wrong again.
    """
    if not os.path.exists(output):
        return True
    try:
        out_at = os.path.getmtime(output)
        return any(os.path.getmtime(i) > out_at
                   for i in inputs if os.path.exists(i))
    except OSError:
        return True


#: Where a finished CAVER run keeps the config it was actually run with.
CAVER_CONF_USED = "chbridge_conf_used.txt"


def _caver_is_stale(run_dir, struct, conf):
    """True if CAVER has to run again for this structure AND this config.

    Timestamps are not enough here, and this is the second half of the same
    lesson as `_is_stale`. That one keyed the cache on "is there output",
    could not see a changed input, and reported results computed on
    superseded coordinates. Keying on the structure's mtime fixed that and
    left the other input uncovered: the SETTINGS. Correct the starting point
    of a run -- which is the whole content of the taxadiene synthase fix --
    and the structure has not changed, so the old clusters are reused, the
    tunnel is still the one from the wrong point, and the report prints the
    corrected start point beside it. The cache would have hidden its own fix.

    So the config is compared by CONTENT, not by mtime: it is rewritten on
    every pass, so its timestamp always looks newer and would force a rerun
    of all 24 every time. A finished run stores the config it used; a run
    whose stored copy differs from the one about to be used is stale.
    """
    if _is_stale(os.path.join(run_dir, "summary.txt"), struct):
        return True
    used = os.path.join(run_dir, CAVER_CONF_USED)
    if not os.path.exists(used):
        return True          # a run from before this check: redo it once
    try:
        with open(used, encoding="utf-8") as a, open(conf, encoding="utf-8") as b:
            return a.read() != b.read()
    except OSError:
        return True


def _run_caver(struct, conf, run_dir, force=False, report=lambda *a: None):
    """Run CAVER unless a finished run of this structure AND config exists.

    Stamps the config into the run directory AFTER CAVER succeeds, so an
    interrupted run is stale rather than falsely complete.
    """
    import shutil

    from channelbridge.caver import run_caver

    if not (force or _caver_is_stale(run_dir, struct, conf)):
        return False
    run_caver(struct, conf, run_dir, report=report)
    shutil.copyfile(conf, os.path.join(run_dir, CAVER_CONF_USED))
    return True


def _find_caver_cluster(run_dir, cluster_id):
    """A single structure makes CAVER write clusters_timeless, not clusters."""
    import glob

    for sub in ("clusters", "clusters_timeless"):
        hits = sorted(glob.glob(os.path.join(
            run_dir, "data", sub, "tun_cl_%03d*.pdb" % cluster_id)))
        if hits:
            return hits[0]
    raise FileNotFoundError("no cluster %d under %s" % (cluster_id, run_dir))


def _normalise(prefix, e, cluster_pdb, args):
    """Stage 2b, through the real CLI, so the example tests what ships.

    Normalisation is the default; ``--raw-only`` is what opts out. Returns
    the numbers rather than a sentence, so a failure here fails the row:
    a tunnel that was fetched but could not be normalised is not a result,
    and reporting it as one is how a broken stage stays hidden.
    """
    out = _enzyme_dir(prefix, "output")
    dst = os.path.join(out, "%s_normalised.pdb" % e["name"])
    cmd = ["tunnel", cluster_pdb, "--out-pdb", dst,
           "--structure", _prepped(prefix, e)]
    proc = _run(cmd)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-1:]
        return dict(ok=False,
                    detail="normalise FAILED (%s)" % ("; ".join(tail) or "no output"))
    row = verdict = ""
    for line in (proc.stdout or "").splitlines():
        parts = line.split()
        if len(parts) >= 8 and parts[1] == "normalised":
            row, verdict = parts, parts[-1]
    if not row:
        return dict(ok=False, detail="normalise printed no 'normalised' row")
    return dict(ok=True, span=row[4], bottleneck=row[5], burial=row[6],
                verdict=verdict,
                detail="normalised: %s spheres, span %s A, bottleneck %s A, "
                       "burial %s, %s" % (row[2], row[4], row[5], row[6], verdict))


def _write_report(results):
    path = os.path.join(HERE, "report.tsv")
    # `branch` is where the enzyme STARTED (what enzymes.tsv says);
    # `resolved` is what actually produced the tunnel, which differs
    # whenever a fetched chain turned out not to reach the surface.
    #
    # `structure` and `source` are the coordinates the row was computed on.
    # They replace a bare `entry`, which named the homolog of the FIRST
    # structure search for all 24 rows and so misattributed every result
    # that a model has since superseded -- four of them to a PDB entry that
    # no stage in this run opened. `superseded_entry` keeps that history
    # where it cannot be mistaken for provenance, and `branch_changed` says
    # why a row that `enzymes.tsv` calls `channelsdb` ran CAVER: it was
    # computed and then dropped on the floor, never written to any column.
    cols = ("prefix", "side", "enzyme", "branch", "branch_changed",
            "resolved", "recovery", "structure", "source", "superseded_entry",
            "loaded", "ok", "n", "detail")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in results:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    ok = sum(1 for r in results if r.get("ok"))
    print("\n%d of %d enzyme sides got a normalised tunnel; wrote %s"
          % (ok, len(results), path))


def cmd_patches(args):
    """Stage 3: classify every residue and rank candidate epitopes.

    Needs PyMOL, which is run as a subprocess. Judge it by the printed
    lines, not the exit code: PyMOL raises ``SystemExit: 0`` through
    ``execfile`` during teardown, so a clean run still prints a traceback.
    """
    enz = enzymes()
    bad = 0
    for pair in _selected(args):
        print("%s" % pair["prefix"])
        for side in ("upstream", "downstream"):
            e = enz[pair[side]]
            out = _enzyme_dir(pair["prefix"], "output")
            tunnel = os.path.join(out, "%s_normalised.pdb" % e["name"])
            cmd = ["patches", _prepped(pair["prefix"], e),
                   "--chain", CHAIN_OF[side], "--name", e["name"],
                   "--out-dir", out, "--pymol", args.pymol]
            # The stage runs without a tunnel -- there is then no mouth, so
            # no CHANNEL or SHADOW class. Worth having anyway for the two
            # enzymes whose site turned out to be enclosed: the
            # classification is still true, it just cannot say which surface
            # sits near a channel. Stage 4 is skipped for those pairs.
            if os.path.exists(tunnel):
                cmd += ["--tunnel", tunnel]
            else:
                print("   .. %-6s no tunnel (see report.tsv): classifying "
                      "without a mouth, so no CHANNEL/SHADOW class"
                      % e["enzyme"])
            if e["channel_ligand"]:
                cmd += ["--ligands", e["channel_ligand"]]
            table = os.path.join(out, "%s_hotspot_patches.tsv" % e["name"])
            proc = _run(cmd)
            lines = [ln.strip() for ln in (proc.stdout or "").splitlines()
                     if "[patches]" in ln]
            # The exit code first, then a table that is actually from THIS
            # run. Existence alone let a stale table from the previous prep
            # stand in for a command that had just failed.
            kept = (proc.returncode == 0
                    and not _is_stale(table, _prepped(pair["prefix"], e), tunnel))
            if not kept:
                bad += 1
                print("   !! %-6s FAILED" % e["enzyme"])
                _say(proc)
                continue
            head = next((ln for ln in lines if "patch(es) kept" in ln), "")
            print("   %-6s %s" % (e["enzyme"],
                                  head.replace("[patches]", "").strip() or "ok"))
    return 1 if bad else 0


def cmd_map(args):
    """Stage 4: the hotspot map, one PyMOL session per pair.

    With no pose search (stage 5 is not built) this takes the pre-sampling
    route: the two enzymes are placed at the closest, best-oriented corner
    of the default grid. The map is for deciding which epitopes to aim at,
    and that does not need the final pose -- but the D_IC and O_IC written
    into the session are the grid corner, NOT a search result, and should
    not be quoted as one.
    """
    enz = enzymes()
    bad = 0
    for pair in _selected(args):
        prefix = pair["prefix"]
        up, dn = enz[pair["upstream"]], enz[pair["downstream"]]
        out = _enzyme_dir(prefix, "output")
        need = []
        cmd = ["hotspots", "--prefix", prefix,
               "--up-pdb", _prepped(prefix, up),
               "--dn-pdb", _prepped(prefix, dn),
               "--up-chains", CHAIN_OF["upstream"],
               "--dn-chains", CHAIN_OF["downstream"],
               "--up-name", up["name"], "--dn-name", dn["name"],
               "--out-dir", out, "--pymol", args.pymol]
        need_all = []
        for flag, e in (("up", up), ("dn", dn)):
            tun = os.path.join(out, "%s_normalised.pdb" % e["name"])
            cls = os.path.join(out, "%s_residue_classes.tsv" % e["name"])
            pat = os.path.join(out, "%s_hotspot_patches.tsv" % e["name"])
            need_all += [tun, cls, pat]
            need += [p for p in (tun, cls, pat) if not os.path.exists(p)]
            cmd += ["--tunnel-%s" % flag, tun,
                    "--classes-%s" % flag, cls,
                    "--patches-%s" % flag, pat]
        if need:
            missing = os.path.basename(need[0])
            if missing.endswith("_normalised.pdb"):
                # Not a missing step: `tunnels` ran and could not produce a
                # tunnel for this enzyme. A map without one side's channel
                # would show a classification with no CHANNEL class and read
                # as if the enzyme had none, so it is not built.
                print("%-22s -- no tunnel for %s, so no map; the enclosed-site "
                      "failure is in report.tsv"
                      % (prefix, missing[:-len("_normalised.pdb")]))
                continue
            bad += 1
            print("%-22s !! missing %s; run `tunnels` then `patches` first"
                  % (prefix, missing))
            continue
        pse = os.path.join(out, "%s_hotspot-map.pse" % prefix)
        # Timestamps, not existence. `_is_stale` was written to abolish the
        # existence cache and this call site was missed: a map built from
        # yesterday's tunnel survived a corrected `tunnels` run and was
        # reported as current. The inputs are the ones named in `cmd`.
        if not (_is_stale(pse, *need_all) or args.force):
            print("%-22s have the map already" % prefix)
            continue
        proc = _run(cmd)
        if not os.path.exists(pse):
            bad += 1
            print("%-22s FAILED" % prefix)
            _say(proc)
            continue
        grid = next((ln.strip() for ln in (proc.stdout or "").splitlines()
                     if "D_IC" in ln and "DO" in ln), "")
        print("%-22s %s" % (prefix, grid or "wrote the map"))
    return 1 if bad else 0


#: How many feasible poses each pair's search keeps. The grid is ordered
#: best-DO-first, so this truncates without biasing the answer -- and the
#: whole 108,624-point grid on nineteen pairs is hours for an example set
#: whose job is to show the stage works on uncurated input.
POSES_WANTED = 5

#: Complexes written per pair. The tables carry every feasible pose; the
#: PDBs are for looking at, and two is enough to see what the search chose.
COMPLEXES_WANTED = 2


def _superseded_poses(out_dir, prefix):
    """Files a previous pose search left in *out_dir*, newest run excluded.

    Everything the pose stage writes per pose: the complexes, the occupancy
    copies and the .pml scripts. The tables are rewritten in place and are
    not listed.
    """
    import glob

    out = []
    for pat in ("%s_geom*_complex.pdb", "%s_geom*_occupancy.pdb",
                "%s_geom*.pml"):
        out += glob.glob(os.path.join(out_dir, pat % prefix))
    return sorted(out)


def cmd_poses(args):
    """Stage 5: search each pair for a pose a bispecific binder could hold.

    Needs no PyMOL -- this stage reads what stages 1-3 wrote and is pure
    numpy. A pair whose enzyme has no tunnel is skipped for the same reason
    stage 4 skips it: without a channel there is no D_IC to search on.
    """
    enz = enzymes()
    bad = 0
    rows = []
    for pair in _selected(args):
        prefix = pair["prefix"]
        up, dn = enz[pair["upstream"]], enz[pair["downstream"]]
        out = _enzyme_dir(prefix, "output")
        need, srcs = [], []
        cmd = ["pose", "--prefix", prefix, "--out-dir", out,
                         "--up-pdb", _prepped(prefix, up),
                         "--dn-pdb", _prepped(prefix, dn),
                         "--up-chain", CHAIN_OF["upstream"],
                         "--dn-chain", CHAIN_OF["downstream"],
                         "--up-name", up["name"], "--dn-name", dn["name"],
                         "--max-ok", str(POSES_WANTED),
                         "--top-n", str(COMPLEXES_WANTED)]
        for flag, e in (("up", up), ("dn", dn)):
            tun = os.path.join(out, "%s_normalised.pdb" % e["name"])
            pat = os.path.join(out, "%s_hotspot_patches.tsv" % e["name"])
            need += [p for p in (tun, pat) if not os.path.exists(p)]
            srcs += [tun, pat, _prepped(prefix, e)]
            cmd += ["--tunnel-%s" % flag, tun, "--patches-%s" % flag, pat]
        ligs = [e["channel_ligand"] for e in (up, dn) if e["channel_ligand"]]
        if ligs:
            cmd += ["--ligands", ",".join(ligs)]
        if need:
            missing = os.path.basename(need[0])
            if missing.endswith("_normalised.pdb"):
                print("%-22s -- no tunnel for %s, so nothing to search on; "
                      "the enclosed-site failure is in report.tsv"
                      % (prefix, missing[:-len("_normalised.pdb")]))
                continue
            bad += 1
            print("%-22s !! missing %s; run `patches` first" % (prefix, missing))
            continue
        table = os.path.join(out, "%s_poses.tsv" % prefix)
        if not _is_stale(table, *srcs) and not args.force:
            print("%-22s have the search already" % prefix)
            rows.append((prefix, _pose_summary(table)))
            continue
        # The previous search's complexes first. The pose stage encodes the
        # geometry in the file NAME, so a re-run writes
        # `geom001_DO0.856_...` beside the old `geom001_DO0.874_...` and
        # nothing overwrites anything: one pose id, two poses, two
        # structures. Stages 6 and 7 then pick by DO, so the SUPERSEDED pose
        # wins wherever it scored higher -- which is how a figure and a
        # design spec came to be built on a structure this set no longer
        # uses. `figures.scan_output` refuses the collision now; this stops
        # making it.
        stale = _superseded_poses(out, prefix)
        for path in stale:
            os.remove(path)
        if stale:
            print("%-22s cleared %d file(s) of the previous search"
                  % (prefix, len(stale)))
        proc = _run(cmd)
        if not os.path.exists(table):
            bad += 1
            print("%-22s FAILED" % prefix)
            _say(proc)
            continue
        best = next((ln.strip() for ln in (proc.stdout or "").splitlines()
                     if ln.startswith("best:")), "no feasible pose")
        print("%-22s %s" % (prefix, best.replace("best: ", "")))
        rows.append((prefix, _pose_summary(table)))
    if rows:
        print("\nfeasible poses per pair (max %d kept)" % POSES_WANTED)
        for prefix, summary in rows:
            print("  %-22s %s" % (prefix, summary))
    return 1 if bad else 0


def _pose_summary(table):
    """One line about a finished search, read back from its own table."""
    with open(table, encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    i_status, i_do, i_d = (head.index("status"), head.index("DO_clamped"),
                           head.index("D_IC"))
    ok = [r for r in rows if r[i_status] == "ok"]
    if not ok:
        counts = {}
        for r in rows:
            counts[r[i_status]] = counts.get(r[i_status], 0) + 1
        return "none feasible of %d evaluated (%s)" % (
            len(rows), ", ".join("%s %d" % kv for kv in sorted(counts.items())))
    best = max(ok, key=lambda r: float(r[i_do]))
    return ("%d feasible of %d evaluated; best DO %s at D_IC %s A"
            % (len(ok), len(rows), best[i_do], best[i_d]))


def shared_entries():
    """Enzymes this set represents by the SAME deposited coordinates.

    Each enzyme here stands in for a plant or yeast protein that has no
    structure of its own, so two enzymes can land on one entry -- and when
    they are close paralogs they usually do. Then the two are not two
    examples: every tunnel, patch and pose is identical, and a pair built
    from one is the same computation as a pair built from the other.

    The set already makes this decision once, deliberately: 7OMT's only
    annotated candidate was 5ICE, which 4'OMT uses, so 7OMT was sent to the
    CAVER branch instead. That judgement was never applied to the other 23,
    and three collisions survived it. This function is how they are found
    rather than noticed, which is what stage 5 did: three pairs of pairs
    returned identical evaluated counts and identical best poses.

    Keyed on `inputs.tsv`, which is what every stage reads. Keyed on
    `enzymes.tsv`'s `(entry, chain)` -- as it was -- it answered for the
    FIRST version of this set and went on answering after a later pass gave
    the collided enzymes models of their own. It reported four shared
    structures and "16 distinct computations" when all 24 inputs were by
    then distinct and the answer was 19; `duplicate_pairs` passed that on,
    and the comparison against the paper dropped three real pairs as
    duplicates of each other.

    The key is the whole computation, not just the file: two enzymes that
    read one file but keep different chains, trim differently or start from
    different residues are two different measurements of it.
    """
    groups = {}
    for name, r in inputs().items():
        key = (r["file"], r["chains"], r["cargo"], r["keep_also"], r["trim"],
               r["start_mode"], r["start_residues"])
        groups.setdefault(key, []).append(name)
    return {k: sorted(v) for k, v in groups.items() if len(v) > 1}


def duplicate_pairs():
    """Pairs that are the same computation, because both sides collide."""
    same = {}
    for names in shared_entries().values():
        for n in names:
            same[n] = names[0]
    seen, dups = {}, []
    for p in pairs():
        key = (same.get(p["upstream"], p["upstream"]),
               same.get(p["downstream"], p["downstream"]))
        if key in seen:
            dups.append((seen[key], p["prefix"]))
        else:
            seen[key] = p["prefix"]
    return dups


def cmd_check(args):
    """What in this set is not as independent as its names suggest."""
    shared = shared_entries()
    print("enzymes measured by the SAME computation (inputs.tsv)")
    if not shared:
        print("  none")
    for key, names in sorted(shared.items(), key=lambda kv: kv[1]):
        print("  %-38s chains %-4s %s"
              % (os.path.basename(key[0]), key[1], " = ".join(names)))
    dups = duplicate_pairs()
    print("\npairs that are therefore the same computation")
    if not dups:
        print("  none")
    for first, second in dups:
        print("  %-22s == %s" % (first, second))
    distinct = len(pairs()) - len(dups)
    print("\n%d pair file(s), %d distinct computation(s)" % (len(pairs()), distinct))
    return 0


def cmd_figures(args):
    """Stage 6: one session and one editable SVG per pair, top pose.

    Two commands per pair, in that order, because choosing the camera is a
    human act: `chbridge session` builds the .pse, `chbridge render` draws
    from whatever session it is pointed at -- including one you opened,
    rotated and saved.

    The arguments below are not boilerplate. Each tunnel is in the frame of
    the PREPPED enzyme, not of the pose, so `tunnel_*_ref` is what carries it
    onto the placement; without it the session draws both channels where the
    unposed enzymes had them and labels the figure with a geometry belonging
    to no placement at all (measured: D_IC 38.0 A against the pose's 6.0 A,
    with nothing failing). `chbridge session` now checks its own geometry
    against the pose's and refuses the mismatch, which is why this is a note
    rather than a trap.
    """
    enz = enzymes()
    bad = 0
    for pair in _selected(args):
        prefix = pair["prefix"]
        up, dn = enz[pair["upstream"]], enz[pair["downstream"]]
        out = _enzyme_dir(prefix, "output")
        fig = os.path.join(out, "figures")
        tun = {k: os.path.join(out, "%s_normalised.pdb" % e["name"])
               for k, e in (("up", up), ("dn", dn))}
        missing = [p for p in tun.values() if not os.path.exists(p)]
        if missing:
            print("%-22s -- no tunnel yet; run `tunnels` first" % prefix)
            continue
        ligands = []
        for e in (up, dn):
            for r in [e["channel_ligand"], *(e.get("keep_also") or [])]:
                if r and r not in ligands:
                    ligands.append(r)
        pse_args = [
            "up_chains=%s" % CHAIN_OF["upstream"],
            "dn_chains=%s" % CHAIN_OF["downstream"],
            "tunnel_up=%s" % tun["up"], "tunnel_dn=%s" % tun["dn"],
            "tunnel_up_ref=%s" % _prepped(prefix, up),
            "tunnel_dn_ref=%s" % _prepped(prefix, dn),
            "label=%s" % prefix,
        ]
        # The hotspot patch tables, so the figure shows the epitope pair THIS
        # pose is built on, in red. Which patch that is comes from the pose's
        # own file name and `chbridge session` fills it in.
        for key, e in (("up", up), ("dn", dn)):
            table = os.path.join(out, "%s_hotspot_patches.tsv" % e["name"])
            if os.path.exists(table):
                pse_args.append("hotspots_%s=%s" % (key, table))
        for key, e in (("up", up), ("dn", dn)):
            if e["channel_ligand"]:
                pse_args.append("channel_ligand_%s=%s" % (key, e["channel_ligand"]))
        if ligands:
            pse_args.append("ligands=%s" % ",".join(ligands))

        cmd = ["session", out, "--prefix", prefix, "--pse", "top1",
               "--render", "pair", "--out-dir", fig, "--pymol", args.pymol]
        for a in pse_args:
            cmd += ["--pse-arg", a]
        pse = os.path.join(fig, "%s_geom001_render.pse" % prefix)
        if _is_stale(pse, tun["up"], tun["dn"], _prepped(prefix, up),
                     _prepped(prefix, dn)) or args.force:
            proc = _run(cmd)
            if proc.returncode or not os.path.exists(pse):
                bad += 1
                print("%-22s SESSION FAILED" % prefix)
                _say(proc)
                continue
            geom = next((ln.strip() for ln in (proc.stdout or "").splitlines()
                         if "geometry " in ln and "D_IC" in ln), "")
            print("%-22s pse  %s" % (prefix, geom.split("geometry")[-1].strip()))
        else:
            print("%-22s have the session already" % prefix)

        svg = os.path.join(fig, "%s_geom001_render_pair.svg" % prefix)
        if _is_stale(svg, pse) or args.force:
            proc = _run(["render", pse, "--prefix", prefix, "--pose", "geom001",
                         "--variant", "pair", "--out-dir", fig,
                         "--width", str(args.width),
                         "--height", str(args.height), "--pymol", args.pymol])
            if proc.returncode or not os.path.exists(svg):
                bad += 1
                print("%-22s RENDER FAILED" % prefix)
                _say(proc)
                continue
            print("%-22s svg  %s (%.2f MB)"
                  % (prefix, os.path.basename(svg),
                     os.path.getsize(svg) / 1e6))
        else:
            print("%-22s have the svg already" % prefix)
        # Copied out of `output/`, which is gitignored, into a committed
        # folder: the README displays these, and a README that points at
        # files nobody who clones the repo will have is not a README. Outside
        # the staleness branch on purpose, so re-running after the copy was
        # added publishes the figures that already exist.
        pub = os.path.join(HERE, "figures")
        os.makedirs(pub, exist_ok=True)
        dst_svg = os.path.join(pub, "%s.svg" % prefix)
        if _is_stale(dst_svg, svg):
            shutil.copyfile(svg, dst_svg)
            print("%-22s published %s" % (prefix, os.path.relpath(dst_svg, HERE)))
        # ...and the flat preview beside it. The SVG is the deliverable --
        # every line, arc and label in it is editable -- but it carries the
        # proteins as an embedded raster, and the sanitiser GitHub puts in
        # front of an SVG in a README drops embedded images. A README that
        # displays the SVG directly would therefore show the overlay floating
        # on nothing. So the page shows the PNG and links the SVG.
        png = os.path.join(fig, "intermediate",
                           "%s_geom001_svgcheck_pair.png" % prefix)
        dst_png = os.path.join(pub, "%s.png" % prefix)
        if os.path.exists(png) and _is_stale(dst_png, png):
            shutil.copyfile(png, dst_png)
    return 1 if bad else 0


def cmd_generate(args):
    """Stage 7: binder-design input from each pair's best pose.

    The pose is what is being handed on, so the two enzymes are conditioned
    where the search put them and must stay there. See
    `channelbridge.generators`, and the warning every generated file carries.
    """
    bad = 0
    for pair in _selected(args):
        prefix = pair["prefix"]
        out = _enzyme_dir(prefix, "output")
        table = os.path.join(out, "%s_poses.tsv" % prefix)
        if not os.path.exists(table):
            print("%-22s -- no pose search yet; run `poses` first" % prefix)
            continue
        gen = os.path.join(out, "generate")
        proc = _run(["generate", out, "--prefix", prefix, "--pose", "top1",
                     "--tool", ",".join(args.tool), "--out-dir", gen])
        if proc.returncode:
            bad += 1
            print("%-22s FAILED" % prefix)
            _say(proc)
            continue
        wrote = [ln.strip() for ln in (proc.stdout or "").splitlines()
                 if ln.strip().startswith("wrote ")]
        print("%-22s %d file(s) -> %s"
              % (prefix, len(wrote), os.path.relpath(gen, HERE)))
        for ln in (proc.stdout or "").splitlines():
            if "note:" in ln or "hotspots" in ln:
                print("      %s" % ln.strip())
    return 1 if bad else 0


def cmd_report(args):
    """Summarise report.tsv: what produced each tunnel, and what failed.

    Grouped by `resolved` rather than `branch`, because `branch` is only
    where a side started: most annotated sides end up on a seeded CAVER run.
    """
    path = os.path.join(HERE, "report.tsv")
    if not os.path.exists(path):
        raise SystemExit("no report.tsv yet; run `tunnels` first")
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    print("how each of the %d enzyme sides got its tunnel" % len(rows))
    by = {}
    for r in rows:
        by.setdefault(r.get("resolved") or "?", []).append(r)
    for how, got in sorted(by.items(), key=lambda kv: -len(kv[1])):
        ok = sum(1 for r in got if r["ok"] == "True")
        print("  %-26s %2d%s" % (how, len(got),
                                 "" if ok == len(got) else "  (%d failed)"
                                 % (len(got) - ok)))

    # the normalise verdict is the last clause of the detail line
    verdicts = {}
    for r in rows:
        if r["ok"] != "True":
            continue
        verdicts.setdefault(r["detail"].rsplit(", ", 1)[-1], []).append(r)
    if verdicts:
        print("\nnormalised-tunnel verdict (only `mouth` gives a usable exit)")
        for v, got in sorted(verdicts.items(), key=lambda kv: -len(kv[1])):
            print("  %-26s %2d" % (v, len(got)))

    seeded = [r for r in rows if r.get("recovery")]
    if seeded:
        print("\ndid a seeded CAVER run recover the annotated route? "
              "(%d sides, %d distinct enzymes)"
              % (len(seeded), len({r["enzyme"] for r in seeded})))
        got = {}
        for r in seeded:
            got.setdefault(r["recovery"], set()).add(r["enzyme"])
        for v, enz in sorted(got.items(), key=lambda kv: -len(kv[1])):
            print("  %-26s %2d  %s" % (v, len(enz), " ".join(sorted(enz))))

    bad = [r for r in rows if r["ok"] != "True"]
    if bad:
        print("\n%d side(s) produced no tunnel:" % len(bad))
        for r in bad:
            print("  %-22s %-6s %s" % (r["prefix"], r["enzyme"], r["detail"]))
    return 1 if bad else 0


def cmd_all(args):
    worst = 0
    for step in (cmd_fetch, cmd_prep, cmd_pairs, cmd_tunnels, cmd_patches,
                 cmd_map, cmd_poses):
        print("\n== %s" % step.__name__.replace("cmd_", ""))
        rc = step(args)
        worst = max(worst, rc)
        if rc and step in (cmd_fetch, cmd_prep, cmd_pairs):
            return rc          # nothing later can work without these
    return worst


ACTIONS = {"fetch": cmd_fetch, "prep": cmd_prep, "pairs": cmd_pairs,
           "tunnels": cmd_tunnels, "patches": cmd_patches, "map": cmd_map,
           "poses": cmd_poses, "figures": cmd_figures,
           "generate": cmd_generate, "report": cmd_report,
           "check": cmd_check, "all": cmd_all}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("action", choices=tuple(ACTIONS))
    ap.add_argument("--only", action="append",
                    help="one example prefix; repeatable")
    ap.add_argument("--force", action="store_true",
                    help="redo work that already has an output")
    ap.add_argument("--tool", action="append", default=[],
                    help="generator to write input for; repeatable, "
                         "comma-separated. Default RFD4-Proteina and "
                         "RFdiffusion3; `chbridge generate --list-tools` "
                         "shows them all")
    ap.add_argument("--width", type=int, default=800,
                    help="render width in px (default 800; 2000 for print, "
                         "which costs about six times as long)")
    ap.add_argument("--height", type=int, default=600)
    ap.add_argument("--pymol", default=os.environ.get("PYMOL", "pymol"),
                    help="the PyMOL to run for stages 3, 4 and 6 "
                         "(default: $PYMOL, else `pymol`)")
    args = ap.parse_args(argv)
    args.tool = [t for v in args.tool for t in str(v).split(",") if t.strip()]
    if not args.tool:
        args.tool = ["RFD4-Proteina", "RFdiffusion3"]
    return ACTIONS[args.action](args)


if __name__ == "__main__":
    raise SystemExit(main())
