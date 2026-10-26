#!/usr/bin/env python3
"""Build and run the iMARS example set: stages 0 to 4, both tunnel branches.

The iMARS paper (Wang et al., Cell 188:1349, 2025) reports 36 enzyme fusion
constructs. Collapsed on the enzyme PAIR -- a linker variant is not a
different pair, and ChannelBridge does not fuse anything -- they are the 19
pairs in ``pairs.tsv``, over the 24 enzymes in ``enzymes.tsv``.

What this is for: exercising both ways of getting a tunnel on a set nobody
curated for the purpose. Twenty-one of the 24 enzymes have a ChannelsDB-
annotated homolog and three do not, so three pairs run one branch on one
enzyme and the other branch on the other, inside a single pair file.

    ./run_examples.py fetch            # the 20 deposited entries (~15 MB)
    ./run_examples.py prep             # one chain per enzyme
    ./run_examples.py pairs            # write the 19 pair.yaml files
    ./run_examples.py tunnels          # stage 0 or 2a, then 2b
    ./run_examples.py report
    ./run_examples.py all --only iMARS_01_4CL-STS

Nothing here is bundled: the structures are deposited entries, fetched on
demand, and the directories they land in are not in version control. What IS
committed is the two tables, the pair files this writes, and this script --
which is to say, every decision, and none of the data.

CAVER is needed for the three ``caver``-branch enzymes only; point $CAVER at
it. The other 21 need network access instead.
"""

from __future__ import annotations

import argparse
import csv
import os
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


def enzymes():
    return {r["enzyme"]: r for r in _rows("enzymes.tsv")}


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
    """The deposited entries. One file per entry, shared between examples."""
    import urllib.request

    os.makedirs(ENTRIES, exist_ok=True)
    want = sorted({e["entry"] for e in enzymes().values()})
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
            src = os.path.join(ENTRIES, "%s.pdb" % e["entry"])
            if not os.path.exists(src):
                raise SystemExit("no %s; run `fetch` first" % src)
            out = _prepped(pair["prefix"], e)
            _enzyme_dir(pair["prefix"], "input")
            cmd = ["prep", src, "--chains", e["chain"], "--out", out,
                   "--relabel", CHAIN_OF[side]]
            keep = [v for v in (e["channel_ligand"],) if v]
            if keep:
                cmd += ["--keep", ",".join(keep)]
            proc = _run(cmd)
            if proc.returncode != 0:
                bad += 1
                print("   !! %-6s %s FAILED" % (e["enzyme"], e["entry"]))
                _say(proc)
                continue
            print("   %-6s %s chain %s -> %s (chain %s)" % (
                e["enzyme"], e["entry"], e["chain"], os.path.basename(out),
                CHAIN_OF[side]))
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
              % (prefix, up["enzyme"], dn["enzyme"], pair["branches"]))
    return 0


def _pair_yaml(pair, up, dn):
    out = [
        "# ChannelBridge pair file -- %s (upstream) to %s (downstream)."
        % (up["enzyme"], dn["enzyme"]),
        "#",
        "# GENERATED by run_examples.py from enzymes.tsv and pairs.tsv; edit",
        "# those, not this. One of the %s pairs of the iMARS fusion set."
        % len(pairs()),
        "#",
        "# Neither structure is the iMARS construct: both are homologous",
        "# deposited entries (see enzymes.tsv for the search that found them,",
        "# and the score, which is how much the homology is worth).",
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
        % (chain, e["entry"], e["chain"]),
    ]
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
        out.append("    caver:")
        out.append("      exclude: [HOH, WAT]")
        out.append("      probe_radius: 1.4")
        out.append("      shell_radius: 10")
        out.append("      shell_depth: 6")
        out.append("      clustering_threshold: 3.5")
        out.append("      seed: 1")
        out.append("      # Derived, not written by hand: the centroid of the")
        out.append("      # residues below sits inside an atom, and CAVER then")
        out.append("      # finds zero tunnels without saying why. See")
        out.append("      # channelbridge.caver.cavity_start_point.")
        out.append("      #   chbridge caver config --structure input/%s.pdb \\"
                   % e["name"])
        out.append("      #       --start-cavity %s" % e["start_spec"])
    else:
        out.append("    # ChannelsDB branch: the tunnels come from entry %s"
                   % e["entry"])
        out.append("    # (%s tunnels, %s). Nothing is computed, and nothing"
                   % (e["cdb_tunnels"], e["cdb_category"]))
        out.append("    # may move the coordinates -- their chain is in this")
        out.append("    # entry's frame.")
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
        print("%s  (%s)" % (pair["prefix"], pair["branches"]))
        for side in ("upstream", "downstream"):
            e = enz[pair[side]]
            try:
                row = (_channelsdb_side if e["branch"] == "channelsdb"
                       else _caver_side)(pair["prefix"], e, side, args)
            except Exception as exc:
                row = dict(ok=False, detail="%s: %s" % (type(exc).__name__, exc))
            row.update(prefix=pair["prefix"], side=side, enzyme=e["enzyme"],
                       branch=e["branch"], entry=e["entry"])
            results.append(row)
            print("   %-9s %-6s %-5s %s" % (
                e["branch"], e["enzyme"], e["entry"],
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
        run_caver,
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
    b = burial_at(db_tunnel.xyz, prot)
    seeds = [(db_tunnel.xyz[0] if b[0] >= b[-1] else db_tunnel.xyz[-1],
              "the buried end of %s" % db_tunnel.label, ["HOH", "WAT"])]
    lig = (e["channel_ligand"] or "").upper()
    lig_xyz = np.array([a.xyz for a in atoms.atoms
                        if (a.resname or "").upper() == lig], dtype=float) \
        if lig else np.empty((0, 3))
    if len(lig_xyz):
        # excluded, or the ligand blocks its own channel
        seeds.append((lig_xyz.mean(axis=0), "the %s centroid" % lig,
                      ["HOH", "WAT", lig]))

    run_dir = os.path.join(out, "caver_%s" % e["name"])
    kept, how, tried = [], None, []
    for seed, how, exclude in seeds:
        params = CaverParams(probe_radius=1.4, shell_radius=10.0,
                             shell_depth=6.0, clustering_threshold=3.5,
                             seed=1, exclude=tuple(exclude))
        conf = write_config(os.path.join(out, "%s_caver.conf" % e["name"]),
                            params, [float(v) for v in seed])
        done = os.path.exists(os.path.join(run_dir, "summary.txt"))
        if not (done and not args.force and not tried):
            run_caver(struct, conf, run_dir, report=lambda *a: None)
        kept = select(read_summary(run_dir), keep=3, floor=0.55).kept
        tried.append(how)
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
    top = kept[0]
    cluster = _find_caver_cluster(run_dir, top.cluster_id)
    match = compare(Tunnel.from_caver_cluster(cluster), db_tunnel)
    norm = _normalise(prefix, e, cluster, args)
    return dict(ok=norm["ok"], resolved="channelsdb-seeded CAVER",
                recovery=match.recovery(),
                detail="seeded from %s; %d clusters kept, took cl%d "
                       "(priority %.3f); %.0f %% of the database chain lies "
                       "on ours (exit gap %.1f A, ours runs further by "
                       "design) -> %s; %s"
                       % (how, len(kept), top.cluster_id, top.priority,
                          100 * match.db_coverage, match.exit_gap,
                          match.recovery(), norm["detail"]))


def _start_spec_for(e, side):
    """``start_spec`` with the chain the PREPARED file actually uses.

    The table names residues in the deposited entry's chain; prep renames
    that chain so the two enzymes of a pair are distinguishable. Carrying
    the old id over meant CAVER was handed residues in a chain that no
    longer exists -- which surfaced as "none of the residues is in the
    structure", and would have been a silent miss had the ids happened to
    exist in the renamed chain too.
    """
    want = CHAIN_OF[side]
    out = []
    for raw in e["start_spec"].split(","):
        raw = raw.strip()
        out.append("%s:%s" % (want, raw.split(":")[-1]) if raw else raw)
    return [r for r in out if r]


def _caver_side(prefix, e, side, args):
    """Stage 2a: derive a start point, run CAVER, select, then normalise."""
    from channelbridge.caver import (
        CaverParams,
        cavity_start_point,
        read_summary,
        run_caver,
        select,
        write_config,
    )

    struct = _prepped(prefix, e)
    if not os.path.exists(struct):
        return dict(ok=False, resolved="caver",
                    detail="not prepped; run `prep` first")
    out = _enzyme_dir(prefix, "output")
    params = CaverParams(probe_radius=1.4, shell_radius=10.0, shell_depth=6.0,
                         clustering_threshold=3.5, seed=1)
    found = cavity_start_point(struct, _start_spec_for(e, side),
                               probe_radius=params.probe_radius)
    conf = write_config(os.path.join(out, "%s_caver.conf" % e["name"]),
                        params, found.xyz)
    run_dir = os.path.join(out, "caver_%s" % e["name"])
    if not (os.path.exists(os.path.join(run_dir, "summary.txt")) and not args.force):
        run_caver(struct, conf, run_dir, report=lambda *a: None)
    sel = select(read_summary(run_dir), keep=3, floor=0.55)
    if not sel.kept:
        return dict(ok=False, resolved="caver",
                    detail="CAVER found no tunnel above the floor")
    top = sel.kept[0]
    cluster = _find_caver_cluster(run_dir, top.cluster_id)
    norm = _normalise(prefix, e, cluster, args)
    return dict(ok=norm["ok"], n=len(sel.kept), resolved="caver",
                detail="start %.1f A off the residue centroid (burial %d); "
                       "%d clusters kept, took cl%d (priority %.3f); %s"
                       % (found.offset, found.burial, len(sel.kept),
                          top.cluster_id, top.priority, norm["detail"]))


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
    cols = ("prefix", "side", "enzyme", "branch", "resolved", "recovery",
            "entry", "ok", "n", "detail")
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
            proc = _run(cmd)
            lines = [ln.strip() for ln in (proc.stdout or "").splitlines()
                     if "[patches]" in ln]
            kept = os.path.exists(os.path.join(
                out, "%s_hotspot_patches.tsv" % e["name"]))
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
        for flag, e in (("up", up), ("dn", dn)):
            tun = os.path.join(out, "%s_normalised.pdb" % e["name"])
            cls = os.path.join(out, "%s_residue_classes.tsv" % e["name"])
            pat = os.path.join(out, "%s_hotspot_patches.tsv" % e["name"])
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
        if os.path.exists(pse) and not args.force:
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
                 cmd_map):
        print("\n== %s" % step.__name__.replace("cmd_", ""))
        rc = step(args)
        worst = max(worst, rc)
        if rc and step in (cmd_fetch, cmd_prep, cmd_pairs):
            return rc          # nothing later can work without these
    return worst


ACTIONS = {"fetch": cmd_fetch, "prep": cmd_prep, "pairs": cmd_pairs,
           "tunnels": cmd_tunnels, "patches": cmd_patches, "map": cmd_map,
           "report": cmd_report, "all": cmd_all}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("action", choices=tuple(ACTIONS))
    ap.add_argument("--only", action="append",
                    help="one example prefix; repeatable")
    ap.add_argument("--force", action="store_true",
                    help="redo work that already has an output")
    ap.add_argument("--pymol", default=os.environ.get("PYMOL", "pymol"),
                    help="the PyMOL to run for stages 3 and 4 "
                         "(default: $PYMOL, else `pymol`)")
    args = ap.parse_args(argv)
    return ACTIONS[args.action](args)


if __name__ == "__main__":
    raise SystemExit(main())
