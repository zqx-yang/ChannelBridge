"""Validate one of our tunnels against ChannelsDB's, inside PyMOL.

    pymol -cq run_validate.py -- ours=tun_cl_001.pdb structure=mine.pdb \
        template=5BSR.pdb report=channelsdb_5bsr.json \
        [method=align|super|cealign] [select="chain A and resi 10-320"] \
        [template_select=...] [tol=2.0]
    pymol -cq run_validate.py -- --selftest

PyMOL is here only for the SUPERPOSITION; the comparison itself is pure and
lives in :mod:`channelbridge.channelsdb`.

**The alignment method is a real choice and the default is not always right.**

``align``   sequence alignment, then iterative outlier rejection. The default.
            Needs sequence similarity, and degrades where identity is low --
            which is exactly where a homologous template is used.
``super``   sequence-independent structural superposition. Better below
            ~30 % identity.
``cealign`` combinatorial extension, purely structural. Best for remote
            homologs and the slowest.

**And what you align on matters more than which method.** A whole-chain fit
answers "is this the same fold". A tunnel comparison needs "is this the same
SITE", and on a homolog those differ by a few angstroms -- the tunnel
comparison then inherits that error and reports a disagreement that is
really a misalignment. Pass ``select``/``template_select`` to fit on the
site's own residues when you know them.

A fusion protein needs this too, for a different reason: aligning a
single-domain template against a two-domain model lets PyMOL choose which
domain to land on, and it will pick the better-scoring one. Restricting the
template to its own chain took five of our enzymes from 22-30 A rms to under
2 A.
"""

from __future__ import annotations

import os
import sys

import numpy as np
from pymol import cmd

_here = os.path.dirname(os.path.abspath(sys.argv[0]))
_src = os.path.dirname(os.path.dirname(_here))
if os.path.isdir(os.path.join(_src, "channelbridge")) and _src not in sys.path:
    sys.path.insert(0, _src)

from channelbridge.channelsdb import best_match, read_report, tunnels  # noqa: E402
from channelbridge.tunnel import Tunnel  # noqa: E402

METHODS = ("align", "super", "cealign")


def _args(argv):
    out = {}
    for a in argv:
        if "=" in a:
            k, v = a.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def superpose_structures(template, structure, method="align",
                         template_select=None, select=None):
    """Move *template* onto *structure*; return its 4x4 matrix and the rms.

    The template moves, not ours, because the database tunnels are in the
    template's frame and have to come with it.
    """
    if method not in METHODS:
        raise SystemExit("method must be one of %s, got %r" % (", ".join(METHODS), method))
    cmd.delete("all")
    cmd.load(template, "tpl")
    cmd.load(structure, "mine")
    t_sel = "tpl and polymer and name CA"
    m_sel = "mine and polymer and name CA"
    if template_select:
        t_sel += " and (%s)" % template_select
    if select:
        m_sel += " and (%s)" % select
    if cmd.count_atoms(t_sel) < 3 or cmd.count_atoms(m_sel) < 3:
        raise SystemExit(
            "fewer than three CA atoms to align (template %d, structure %d); "
            "check the selections" % (cmd.count_atoms(t_sel), cmd.count_atoms(m_sel)))
    if method == "cealign":
        res = cmd.cealign(m_sel, t_sel)
        rms, n = float(res["RMSD"]), int(res["alignment_length"])
    else:
        fn = cmd.align if method == "align" else cmd.super
        out = fn(t_sel, m_sel, cycles=5, transform=1)
        rms, n = float(out[0]), int(out[1])
    return np.array(cmd.get_object_matrix("tpl")).reshape(4, 4), rms, n


def run(a):
    matrix, rms, n = superpose_structures(
        a["template"], a["structure"], a.get("method", "align"),
        a.get("template_select"), a.get("select"))
    R, t = matrix[:3, :3], matrix[:3, 3]
    print("[validate] %s superposition: rms %.2f A over %d atoms"
          % (a.get("method", "align"), rms, n))
    if rms > 5.0:
        print("[validate] WARNING rms %.2f A is too large to compare tunnels "
              "against. Try method=cealign, or restrict the selections to one "
              "domain or to the site." % rms)
    db = tunnels(read_report(a["report"]))
    if not db:
        print("[validate] the report has no tunnels")
        return 1
    ours = Tunnel.from_caver_cluster(a["ours"])
    prot = cmd.get_coords("mine and polymer")
    if prot is not None and len(prot) and a.get("truncate", "1") not in ("0", "no"):
        try:
            ours = ours.oriented(prot).truncate_at_surface(prot)
            print("[validate] our tunnel truncated at the surface (%d spheres)"
                  % len(ours.xyz))
        except Exception as exc:
            print("[validate] could not truncate (%s); comparing the raw chain" % exc)
    best = best_match(ours, db, R, t, tol=float(a.get("tol", 2.0)))
    print("[validate] best match: %s" % best.db_label)
    print("    exit gap   %6.2f A   <- what distinguishes one cluster from another"
          % best.exit_gap)
    print("    hausdorff  %6.2f A" % best.hausdorff)
    print("    coverage   %6.0f %%" % (100 * best.coverage))
    print("    verdict:   %s" % best.verdict())
    if best.verdict() != "same tunnel":
        print("[validate] NOT a confirmation. A database tunnel computed on a "
              "different structure need not correspond to yours: different "
              "start point, different probe, different surface. Treat this as "
              "'no independent support', not as 'your tunnel is wrong'.")
    return 0


def _selftest():
    ok = [True]

    def check(c, what):
        print("  %s %s" % ("ok  " if c else "FAIL", what))
        if not c:
            ok[0] = False

    from channelbridge.channelsdb import DbTunnel, compare, superpose

    # a superposition recovered from a known transform
    rng = np.random.default_rng(0)
    P = rng.normal(size=(30, 3)) * 10
    th = np.radians(37.0)
    R0 = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    t0 = np.array([5.0, -2.0, 1.0])
    R, _t, rms = superpose(P, P @ R0.T + t0)
    check(rms < 1e-6, "an exact transform is recovered exactly (rms %.2e)" % rms)
    check(np.allclose(R, R0, atol=1e-6), "...and it is the right rotation")

    check(np.linalg.det(R) > 0, "the fit is a rotation, not a reflection")
    try:
        superpose(P[:2], P[:2])
        check(False, "two atoms should be refused")
    except ValueError:
        check(True, "fewer than three paired atoms is refused")
    try:
        superpose(P, P[:10])
        check(False, "mismatched arrays should be refused")
    except ValueError:
        check(True, "mismatched arrays are refused")

    # the comparison, in a frame that needs the transform
    pts = np.array([[float(i), 0, 0] for i in range(10)])
    ours = Tunnel(xyz=pts, radii=np.full(10, 2.0), source="t")
    moved = pts @ R0.T + t0
    db = DbTunnel("CofactorTunnels", "Caver", "x", moved, np.full(10, 2.0),
                  np.full(10, 2.1), np.arange(10.0), {})
    check(compare(ours, db).exit_gap > 1.0, "unsuperposed, it does not match")
    check(compare(ours, db, R0.T, -t0 @ R0).exit_gap < 1e-6 or
          compare(ours, db, np.linalg.inv(R0), -np.linalg.inv(R0) @ t0).exit_gap < 1e-6,
          "with the inverse transform it matches exactly")
    print("[selftest] %s" % ("PASSED" if ok[0] else "FAILED"))
    return 0 if ok[0] else 1


def main(argv):
    if "--selftest" in argv:
        return _selftest()
    a = _args(argv)
    missing = [k for k in ("ours", "structure", "template", "report") if k not in a]
    if missing:
        print(__doc__)
        print("missing: %s" % ", ".join(missing))
        return 2
    return run(a)


raise SystemExit(main(sys.argv[1:]))
