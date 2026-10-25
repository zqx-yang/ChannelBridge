"""Stage 3 inside PyMOL: measure the residues, then classify and cluster.

    pymol -cq run_patches.py -- structure=enzyme.pdb chain=B name=GDH \
        tunnel=cl1_norm.pdb ligands=NAP,BGC out_dir=. [avoid="resi 145+158"]
    pymol -cq run_patches.py -- --selftest

Only the MEASUREMENTS happen here -- solvent accessibility, secondary
structure, CB positions, distances to the channel mouth. Every decision made
from them lives in :mod:`channelbridge.patches`, which imports no PyMOL and
is tested on its own. This file should stay boring.

Two measurements are easy to get subtly wrong, so they are spelled out:

* SASA is computed TWICE -- once on the whole structure and once on the
  isolated chain. The difference is what the rest of the assembly buries,
  which is what separates an IFACE residue from a merely buried one. Taking
  only the assembly value silently classifies every subunit contact as
  BURIED and hands the search a protein with no interface.
* `dot_solvent` must be on. Without it `get_area` returns the solvent
  EXCLUDED surface, a different quantity on a different scale, and every
  exposure threshold in PatchParams would be meaningless against it.
"""

from __future__ import annotations

import os
import sys

import numpy as np
from pymol import cmd

# A file PyMOL execfile()s gets __file__ = pymol/__init__.py, NOT this
# script; sys.argv[0] is the real path. Same workaround as run_pse.py, and
# the same reason: it cost a run to find once already.
_here = os.path.dirname(os.path.abspath(sys.argv[0]))
_src = os.path.dirname(os.path.dirname(_here))          # .../src
if os.path.isdir(os.path.join(_src, "channelbridge")) and _src not in sys.path:
    sys.path.insert(0, _src)

from channelbridge.patches import (  # noqa: E402
    PatchParams,
    Residue,
    classify,
    find_patches,
    grow_avoid,
    write_patch_table,
    write_residue_table,
)
from channelbridge.tunnel import fitted_axis  # noqa: E402


def _args(argv):
    out = {}
    for a in argv:
        if "=" in a:
            k, v = a.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _sasa(obj, chain):
    """Per-residue solvent accessible area of one chain of *obj*.

    ONE `get_area` call with ``load_b=1``, which writes each atom's area into
    its B-factor, then summed per residue. Calling `get_area` once per
    residue instead is the obvious way and recomputes the dot surface of the
    whole object every time: on a 260-residue chain that was minutes per
    enzyme, against under a second this way, for identical numbers.

    `dot_solvent` must be 1. At 0 `get_area` returns the solvent EXCLUDED
    surface -- a different quantity on a different scale -- and every
    exposure threshold in PatchParams would be measured against the wrong
    thing while still looking plausible.
    """
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", 3)
    sel = "%s and chain %s and polymer" % (obj, chain)
    cmd.get_area(sel, load_b=1)
    out = {}
    cmd.iterate(sel, "out[resi] = out.get(resi, 0.0) + b", space={"out": out})
    return out


def _resi_set(sel):
    acc = set()
    cmd.iterate("byres (%s)" % sel, "acc.add(resi)", space={"acc": acc})
    return acc


def measure(structure, chain, *, tunnel=None, ligands=(), params=None,
            axis_spheres=5):
    """Load a structure and return ``(residues, core, channel)``."""
    params = params or PatchParams()
    cmd.delete("all")
    cmd.load(structure, "prot")
    cmd.remove("solvent")
    cmd.dss("prot")

    meta = {}
    cmd.iterate("prot and chain %s and polymer and name CA" % chain,
                "meta[resi] = (resn, ss)", space={"meta": meta})
    if not meta:
        raise SystemExit("no CA atoms in chain %r of %s" % (chain, structure))

    sasa_all = _sasa("prot", chain)
    cmd.create("mono", "prot and chain %s" % chain)
    sasa_mono = _sasa("mono", chain)
    cmd.delete("mono")

    lig_sel = " or ".join("resn %s" % r for r in ligands) if ligands else None
    core = _resi_set("prot and polymer and chain %s within %g of "
                     "(prot and chain %s and (%s))"
                     % (chain, params.ligand_r, chain, lig_sel)) if lig_sel else set()

    # Several tunnels may be given. The FIRST is the primary: its geometry
    # is what the table reports and what CHANNEL is measured from. All of
    # them contribute a shadow cylinder, because a binder capping a secondary
    # mouth blocks the channel just as surely as one capping the primary.
    tunnels = [t for t in str(tunnel or "").replace(",", " ").split() if t]
    openings = []
    channel = set()
    for i, path in enumerate(tunnels):
        cmd.load(path, "tun")
        xyz = cmd.get_coords("tun")
        if xyz is None or len(xyz) < 2:
            raise SystemExit("tunnel %s has fewer than two spheres" % path)
        mouth = np.asarray(xyz[-1], float)
        # The FITTED axis, over the outer spheres -- not the last two. A
        # single pair of CAVER spheres sits ~1 A apart and its direction is
        # dominated by placement jitter; the two differ by 9-12 degrees on
        # these tunnels, which is enough to move residues in and out of the
        # SHADOW class. (O_IC itself still uses the short arm: see
        # Tunnel.axis vs Tunnel.fitted_axis.)
        axis = fitted_axis(xyz, axis_spheres, source=path)
        openings.append((mouth, axis))
        if i == 0:
            channel = _resi_set("prot and polymer and chain %s within %g of tun"
                                % (chain, params.channel_r))
        cmd.delete("tun")

    residues = {}
    for resi, (resn, ss) in meta.items():
        cb = cmd.get_coords("prot and chain %s and resi %s and name CB"
                            % (chain, resi))
        if cb is None or not len(cb):
            cb = cmd.get_coords("prot and chain %s and resi %s and name CA"
                                % (chain, resi))
        cb = np.asarray(cb[0], float) if cb is not None and len(cb) else np.zeros(3)
        d = th = ax = perp = 0.0
        per_opening = []
        for j, (mouth, axis) in enumerate(openings):
            v = cb - mouth
            o_ax = float(v @ axis)
            o_perp = float(np.sqrt(max(0.0, float(v @ v) - o_ax * o_ax)))
            per_opening.append((o_ax, o_perp))
            if j:
                continue
            ax, perp = o_ax, o_perp
            d = float(np.linalg.norm(v))
            if d > 1e-9:
                th = float(np.degrees(np.arccos(
                    np.clip(float(v @ axis) / d, -1.0, 1.0))))
        a_all = float(sasa_all.get(resi, 0.0))
        residues[resi] = Residue(
            resi=resi, resn=resn, ss=ss or "L", sasa=a_all,
            dsasa=float(sasa_mono.get(resi, a_all)) - a_all,
            cb=tuple(float(x) for x in cb),
            dist_mouth=d, angle_mouth=th, axial=ax, perp=perp,
            openings=per_opening)
    return residues, core, channel


def run(a):
    params = PatchParams()
    for key in vars(params):
        if key in a:
            setattr(params, key, type(getattr(params, key))(a[key]))
    chain = a.get("chain", "A")
    name = a.get("name") or os.path.splitext(os.path.basename(a["structure"]))[0]
    ligands = [x for x in a.get("ligands", "").replace(",", " ").split() if x]
    residues, core, channel = measure(
        a["structure"], chain, tunnel=a.get("tunnel"), ligands=ligands,
        params=params, axis_spheres=int(a.get("axis_spheres", 5)))

    named = set()
    if a.get("avoid"):
        for spec in a["avoid"].split(";"):
            spec = spec.strip()
            if spec:
                named |= _resi_set("prot and chain %s and (%s)" % (chain, spec))
        if not named:
            print("[patches] WARNING avoid= matched nothing: %s" % a["avoid"])
    avoid = grow_avoid(named, residues, params)

    classes = {r: classify(v, in_core=r in core, in_channel=r in channel,
                           in_avoid=r in avoid, params=params)
               for r, v in residues.items()}
    counts = {}
    for c in classes.values():
        counts[c] = counts.get(c, 0) + 1
    print("[patches] %s chain %s: %d residues, %s"
          % (name, chain, len(residues),
             ", ".join("%s %d" % (k, counts[k]) for k in sorted(counts))))
    if named:
        print("[patches] avoid: %d named, %d excluded after growing %.0f A "
              "from the exposed ones" % (len(named), len(avoid), params.avoid_r))

    patches = find_patches(residues, classes, params)
    print("[patches] %d patch(es) kept of %d ALLOWED residues"
          % (len(patches), counts.get("ALLOWED", 0)))
    for p in patches:
        print("    %-9s score %5.2f  geom %.2f  %2d res, %d anchor(s), "
              "%3.0f A from the mouth" % (p.label, p.score, p.geom_score,
                                          p.n_res, p.anchors, p.dist_mouth))

    out_dir = a.get("out_dir", ".")
    os.makedirs(out_dir, exist_ok=True)
    rt = write_residue_table(
        os.path.join(out_dir, "%s_residue_classes.tsv" % name), chain,
        residues, classes)
    pt = write_patch_table(
        os.path.join(out_dir, "%s_hotspot_patches.tsv" % name), patches)
    print("[patches] wrote %s" % rt)
    print("[patches] wrote %s" % pt)
    return 0


def _selftest():
    import tempfile

    ok = [True]

    def check(cond, what):
        print("  %s %s" % ("ok  " if cond else "FAIL", what))
        if not cond:
            ok[0] = False

    tmp = tempfile.mkdtemp()
    pdb = os.path.join(tmp, "t.pdb")
    # two layers of residues: an outer shell that is exposed and an inner one
    # that is not, so the classification has something to separate
    lines = []
    n = 0
    for i, (x, resn) in enumerate([(0.0, "PHE"), (3.8, "LEU"), (7.6, "ILE"),
                                   (11.4, "TRP"), (15.2, "MET"), (19.0, "ALA")]):
        for nm, dz in (("N", 0.0), ("CA", 1.2), ("C", 2.4), ("O", 3.0), ("CB", 1.8)):
            n += 1
            lines.append(
                "ATOM  %5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C"
                % (n, nm, resn, i + 1, x, 0.0, dz))
    open(pdb, "w").write("\n".join(lines) + "\nEND\n")

    res, _core, _chan = measure(pdb, "A")
    check(len(res) == 6, "six residues measured (%d)" % len(res))
    check(all(r.sasa > 0 for r in res.values()),
          "every residue has a positive SASA -- dot_solvent is on")
    check(all(r.cb != (0.0, 0.0, 0.0) for r in res.values()),
          "every residue found a CB or CA")

    params = PatchParams()
    classes = {r: classify(v, in_core=False, in_channel=False, in_avoid=False,
                           params=params) for r, v in res.items()}
    check(set(classes.values()) <= set(
        ("LIGAND_CORE", "CHANNEL", "SHADOW", "AVOID", "IFACE", "BURIED",
         "ALLOWED")), "every label is one of the seven")

    # the avoid growth, which is the part most likely to be got wrong
    exposed = max(res, key=lambda r: res[r].sasa)
    grown = grow_avoid([exposed], res, params)
    check(len(grown) > 1, "an exposed named residue spreads its exclusion")
    buried = dict(res)
    buried[exposed] = Residue(resi=exposed, resn=res[exposed].resn,
                              sasa=1.0, cb=res[exposed].cb)
    check(grow_avoid([exposed], buried, params) == {exposed},
          "a buried one protects only itself")

    out = write_patch_table(os.path.join(tmp, "p.tsv"),
                            find_patches(res, classes, params))
    head = open(out).readline().rstrip("\n").split("\t")
    check(head[0] == "rank" and head[-1] == "sequence",
          "the patch table keeps the original column names")

    print("[selftest] %s" % ("PASSED" if ok[0] else "FAILED"))
    return 0 if ok[0] else 1


def main(argv):
    if "--selftest" in argv:
        return _selftest()
    a = _args(argv)
    if "structure" not in a:
        print(__doc__)
        return 2
    return run(a)


# Called unconditionally: a file PyMOL execfile()s does NOT get
# __name__ == "__main__", so an `if __name__` guard would silently do
# nothing and exit 0.
raise SystemExit(main(sys.argv[1:]))
