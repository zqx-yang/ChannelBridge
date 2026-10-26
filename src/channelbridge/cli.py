"""The ``chbridge`` command line.

Subcommands land here as the stages are ported out of the P5 prototype
scripts. What is wired up now is the part that is independent of PyMOL:

``chbridge validate``
    Read a pair file, report what it says and what is wrong with it.
``chbridge tunnel``
    Normalise a CAVER tunnel cluster and print its geometry. This is the
    command to reach for when two CAVER runs disagree.
``chbridge do``
    The DO score for a given D_IC and O_IC.
``chbridge figures``
    Sessions, editable SVGs and structures for chosen poses. Needs PyMOL,
    which it runs as a subprocess rather than importing.

Every subcommand exits non-zero on failure and prints to stderr, so it can
be used in a shell pipeline.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from . import __version__

EPILOG = """\
examples:
  chbridge validate examples/gdh_fnr/pair.yaml
  chbridge tunnel caver/output/out_GDH_primary_sr14 \\
      --structure caver/in/GDH/gdh.pdb --exclude HOH,WAT,BGC,NAP
  chbridge tunnel a/out_X b/out_X --structure s.pdb --exclude HOH,WAT --compare
  chbridge do --dic 16 --oic 150
  chbridge figures out-complex_GDH_FNR --svg top3 --render both --prefix run1
  chbridge session out-complex_GDH_FNR --pse rank1 --prefix run1
  chbridge render figures/run1_geom001_render.pse --prefix run1 --variant both
"""


# --------------------------------------------------------------------- helpers


def _fail(msg: str) -> int:
    print("error: %s" % msg, file=sys.stderr)
    return 2


def _split_csv(value: Optional[str]) -> List[str]:
    if not value:
        return []
    return [v.strip().upper() for v in value.split(",") if v.strip()]


def _unique_labels(paths: Sequence[str]) -> List[str]:
    """Shortest tail of each path that tells them apart.

    Two CAVER runs of the same config by different builds have the same
    ``out_<tag>`` name, so labelling by basename alone would print two
    identical rows.
    """
    parts = [Path(p).resolve().parts for p in paths]
    depth = 1
    while depth <= max(len(p) for p in parts):
        labels = ["/".join(p[-depth:]) for p in parts]
        if len(set(labels)) == len(labels):
            return labels
        depth += 1
    return [str(Path(p).resolve()) for p in paths]  # pragma: no cover


def _find_cluster(caver_out: Path, cluster: Optional[int]) -> Path:
    """Locate a tunnel cluster PDB under a CAVER output folder."""
    if caver_out.is_file():
        return caver_out
    # CAVER names this folder by whether the run had a time dimension.
    # "clusters" is for a trajectory; a SINGLE structure -- the only kind
    # this pipeline runs -- gets "clusters_timeless", and its files carry a
    # trailing snapshot index (tun_cl_001_1.pdb). Looking only in
    # "clusters" made every static run look like a run that found nothing.
    candidates = [caver_out / "data" / "clusters",
                  caver_out / "data" / "clusters_timeless"]
    d = next((c for c in candidates if c.is_dir()), None)
    if d is None:
        raise FileNotFoundError(
            "%s has neither data/clusters/ nor data/clusters_timeless/. Either "
            "CAVER found no pathway (check the starting point is not inside an "
            "atom, and read warnings.txt -- CAVER exits 0 either way), or this "
            "is not a CAVER output folder." % caver_out
        )
    found = sorted(d.glob("tun_cl_*.pdb"))
    if not found:
        raise FileNotFoundError("%s holds no tun_cl_*.pdb" % d)
    if cluster is None:
        return found[0]
    want = sorted(d.glob("tun_cl_%03d*.pdb" % cluster))
    if not want:
        raise FileNotFoundError(
            "no tun_cl_%03d*.pdb in %s; this run has clusters %s"
            % (cluster, d, ", ".join(sorted({p.name.split("_")[2] for p in found})))
        )
    return want[0]


# ------------------------------------------------------------------- validate


def cmd_validate(args) -> int:
    from .config import ConfigError, load_pair_file

    try:
        pair = load_pair_file(args.pair_file)
    except ConfigError as exc:
        return _fail(str(exc))
    print(pair.summary())
    problems = pair.validate_against_disk()
    if problems:
        print()
        for p in problems:
            print("PROBLEM  %s" % p, file=sys.stderr)
        return 1
    print("\nok: pair file is valid and every target file is present")
    return 0


# --------------------------------------------------------------------- tunnel


def cmd_tunnel(args) -> int:
    from .io.pdb import read_pdb
    from .tunnel import Tunnel

    exclude = _split_csv(args.exclude)
    protein = None
    if args.structure:
        structure = Path(args.structure)
        if not structure.is_file():
            return _fail("structure not found: %s" % structure)
        arr = read_pdb(structure)
        if exclude:
            arr = arr.exclude_resnames(exclude)
        protein = arr.xyz
        if len(protein) == 0:
            return _fail(
                "every atom of %s was excluded by --exclude %s"
                % (structure.name, ",".join(exclude))
            )

    rows = []
    labels = _unique_labels(args.caver_out)
    for raw_path, label in zip(args.caver_out, labels):
        try:
            path = _find_cluster(Path(raw_path), args.cluster)
            tun = Tunnel.from_caver_cluster(path)
        except (FileNotFoundError, ValueError) as exc:
            return _fail(str(exc))
        entries = [("raw", tun)]
        if protein is not None and not args.raw_only:
            entries.append(
                (
                    "normalised",
                    tun.normalised(
                        protein,
                        burial_min=args.burial_min,
                        spacing=args.spacing,
                        radius=args.burial_radius,
                    ),
                )
            )
        for kind, t in entries:
            d = t.describe(protein)
            d["run"] = label
            d["kind"] = kind
            try:
                d["arm_chord_A"] = round(t.arm_chord(args.oic_arm), 3)
            except ValueError:
                d["arm_chord_A"] = float("nan")
            rows.append(d)

    width = max([len("run")] + [len(r["run"]) for r in rows])
    cols = [
        ("run", "%-" + str(width) + "s", "run"),
        ("kind", "%-11s", "kind"),
        ("spheres", "%7s", "spheres"),
        ("mean_spacing_A", "%9s", "spacing"),
        ("span_A", "%8s", "span"),
        ("bottleneck_A", "%11s", "bottleneck"),
        ("burial_exit", "%7s", "burial"),
        ("oic_arm_spread_deg", "%7s", "spread"),
        ("verdict", "%-18s", "verdict"),
    ]
    print("  ".join(fmt % head for _, fmt, head in cols))
    for r in rows:
        print("  ".join(fmt % r.get(key, "-") for key, fmt, _ in cols))

    print()
    for r in rows:
        print(
            "%s / %s: exit %8.3f %8.3f %8.3f   O_IC arm chord %.3f A"
            % (r["run"], r["kind"], r["exit_x"], r["exit_y"], r["exit_z"], r["arm_chord_A"])
        )
        print("    %s" % r["history"])

    if args.compare and len(rows) >= 2:
        import numpy as np

        print()
        groups = {}
        for r in rows:
            groups.setdefault(r["kind"], []).append(r)
        for kind, rs in groups.items():
            if len(rs) < 2:
                continue
            pts = [np.array([r["exit_x"], r["exit_y"], r["exit_z"]]) for r in rs]
            worst = max(
                float(np.linalg.norm(pts[i] - pts[j]))
                for i in range(len(pts))
                for j in range(i + 1, len(pts))
            )
            print("%-11s largest exit separation between runs: %6.3f A" % (kind, worst))
        print(
            "\nA large raw separation that shrinks when normalised is the sphere-"
            "sampling\nartefact described in channelbridge.tunnel, not a different "
            "tunnel."
        )

    if args.out_pdb:
        if protein is None:
            return _fail("--out-pdb needs --structure: without it there is nothing to normalise")
        if len(args.caver_out) != 1:
            return _fail(
                "--out-pdb writes one file, so give exactly one CAVER output "
                "(got %d)" % len(args.caver_out)
            )
        tun = Tunnel.from_caver_cluster(_find_cluster(Path(args.caver_out[0]), args.cluster))
        norm = tun.normalised(
            protein,
            burial_min=args.burial_min,
            spacing=args.spacing,
            radius=args.burial_radius,
        )
        n = norm.write_caver_pdb(args.out_pdb)
        print("\n[out] %s  (%d spheres, CAVER tunnel format)" % (args.out_pdb, n))
        print("      %s" % "; ".join(norm.history))

    if protein is None:
        print(
            "\nnote: without --structure the chain cannot be oriented or truncated "
            "by burial,\nso the exit and span are whatever the CAVER build wrote.",
            file=sys.stderr,
        )
    return 0


# ------------------------------------------------------------------------- do


def cmd_do(args) -> int:
    from .config import DO_SETS
    from .geometry import do_score

    coeffs = DO_SETS[args.coefficients]
    try:
        value = do_score(args.dic, args.oic, coeffs)
    except ValueError as exc:
        return _fail(str(exc))
    print("D_IC %.3f A   O_IC %.2f deg   DO %.4f   [%s]" % (args.dic, args.oic, value, coeffs.name))
    if args.all_sets:
        for name, c in sorted(DO_SETS.items()):
            if name != coeffs.name:
                print(
                    "  %-14s DO %.4f  (difference %+.4f)"
                    % (name, c(args.dic, args.oic), c(args.dic, args.oic) - value)
                )
    print(
        "  1 A of D_IC is worth %.1f deg of O_IC at these coefficients"
        % coeffs.degrees_per_angstrom
    )
    return 0


# ----------------------------------------------------------------------- main


def cmd_figures(args) -> int:
    from .figures import Runner, build, collect_poses, scan_output
    from .render.select import resolve
    from .render.style import RenderStyle
    from .render.svg import SvgStyle

    try:
        geoms = scan_output(args.output_dir)
    except FileNotFoundError as exc:
        return _fail(str(exc))
    if not geoms:
        return _fail("no geom*_complex.pdb in %s" % args.output_dir)
    print("%d pose(s) in %s" % (len(geoms), args.output_dir))

    style = RenderStyle()
    if args.pair_file:
        from .config import ConfigError, load_pair_file

        try:
            pair = load_pair_file(args.pair_file)
        except (ConfigError, OSError) as exc:
            return _fail(str(exc))
        style = getattr(pair, "render", None) or style
    try:
        svg_style = SvgStyle.from_dict(getattr(style, "svg", None))
    except ValueError as exc:
        return _fail(str(exc))

    existing = [g.pose_id for g in geoms if g.pse]
    try:
        plan = resolve(
            collect_poses(geoms),
            pse=args.pse, svg=args.svg, pdb=args.pdb,
            render=args.render,
            existing_pse=existing, rebuild=args.rebuild,
        )
    except ValueError as exc:
        return _fail(str(exc))
    print(plan.summary())
    if plan.nothing_to_do:
        print("nothing selected; use --pse/--svg/--pdb (e.g. --svg top3)")
        return 0

    out_dir = args.out_dir or os.path.join(args.output_dir, "figures")
    runner = Runner(pymol=args.pymol, dry_run=args.dry_run)
    made = build(plan, geoms, out_dir, runner=runner, width=args.width, height=args.height,
                 pse_args=_split_csv(args.pse_arg), style=style, svg_style=svg_style,
                 prefix=args.prefix)
    for kind in ("pse", "svg", "pdb"):
        for path in made[kind]:
            print("  %-4s %s" % (kind, path))
    for bad in made["failed"]:
        print("  FAILED %s" % bad, file=sys.stderr)
    return 1 if made["failed"] else 0


def cmd_session(args) -> int:
    """Build the .pse only. Rendering is a separate step on purpose.

    Separating them means you can open the session, rotate it to a frame you
    like, save it, and render THAT -- or hand `chbridge render` a session you
    built yourself. A command that renders as soon as it builds forecloses both.
    """
    from .figures import Runner, collect_poses, scan_output, session_name
    from .render.select import resolve

    try:
        geoms = scan_output(args.output_dir)
    except FileNotFoundError as exc:
        return _fail(str(exc))
    if not geoms:
        return _fail("no geom*_complex.pdb in %s" % args.output_dir)
    try:
        plan = resolve(collect_poses(geoms), pse=args.pse or "top1",
                       render=args.render, rebuild=True)
    except ValueError as exc:
        return _fail(str(exc))
    out_dir = args.out_dir or os.path.join(args.output_dir, "figures")
    runner = Runner(pymol=args.pymol, dry_run=args.dry_run)
    by_id = {g.pose_id: g for g in geoms}
    made = []
    for pid in plan.pse:
        g = by_id.get(pid)
        if g is None:
            print("  skip %s: not in the output folder" % pid, file=sys.stderr)
            continue
        dst = os.path.join(out_dir, session_name(args.prefix, pid))
        print("  [%s] -> %s" % (pid, dst))
        rc = runner.run(runner.script("run_pse.py"),
                        ["complex=%s" % g.complex_pdb, "out=%s" % dst,
                         "render=%s" % args.render, *_split_csv(args.pse_arg)])
        if rc:
            print("  FAILED %s: run_pse exited %d" % (pid, rc), file=sys.stderr)
        else:
            made.append(dst)
    for m in made:
        print("  pse  %s" % m)
    print("\nnow render it, or rotate it first and save, then:"
          "\n  chbridge render <session.pse> --prefix %s" % args.prefix)
    return 0 if made or args.dry_run else 1


def cmd_render(args) -> int:
    """Render scenes and assemble SVGs from an EXISTING .pse."""
    from .figures import Runner, helper_name, svg_name
    from .render.style import RenderStyle
    from .render.svg import SvgStyle, build_svg

    pse = args.session
    if not os.path.isfile(pse):
        return _fail("no such session: %s" % pse)
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(pse))
    inter = os.path.join(out_dir, "intermediate")
    if not args.dry_run:
        os.makedirs(inter, exist_ok=True)
    runner = Runner(pymol=args.pymol, dry_run=args.dry_run)
    pose = args.pose or "geom001"
    style = RenderStyle()
    svg_style = SvgStyle.from_dict(getattr(style, "svg", None))

    variants = ["pair", "assembly"] if args.variant == "both" else [args.variant]
    made = []
    for v in variants:
        sfx = "_assembly" if v == "assembly" else ""
        helpers = {}
        ok = True
        for h in ("mask", "mouth", "footprint"):
            path = os.path.join(inter, helper_name(args.prefix, pose, h + sfx))
            helpers[h] = path
            if os.path.exists(path) and not args.dry_run:
                continue
            if runner.run(runner.script("run_png.py"),
                          ["pse=%s" % pse, "scene=%s" % (h + sfx), "out=%s" % path,
                           "w=%d" % args.width, "h=%d" % args.height]):
                ok = False
        base = os.path.join(inter, helper_name(args.prefix, pose,
                                               "base_%s_transparent" % v))
        if not (os.path.exists(base) and not args.dry_run):
            if runner.run(runner.script("run_png.py"),
                          ["pse=%s" % pse, "scene=base_%s" % v, "out=%s" % base,
                           "w=%d" % args.width, "h=%d" % args.height,
                           "transparent=1"]):
                ok = False
        if not ok:
            print("  FAILED %s: a render step failed" % v, file=sys.stderr)
            continue
        dst = os.path.join(out_dir, svg_name(args.prefix, pose, v))
        if args.dry_run:
            print("  would assemble %s" % dst)
            made.append(dst)
            continue
        try:
            axis_px = None
            if args.axis_px:
                try:
                    axis_px = [float(x) for x in args.axis_px.split(",")]
                except ValueError:
                    return _fail("--axis-px wants 'x,y', got %r" % args.axis_px)
                if len(axis_px) != 2:
                    return _fail("--axis-px wants two numbers, got %r" % args.axis_px)
            build_svg(base, helpers["mask"], helpers["mouth"], helpers["footprint"],
                      dst, {"D_IC": args.dic, "O_IC": args.oic, "DO": args.do},
                      axis_px=axis_px,
                      style=style, svg_style=svg_style,
                      title="%s %s %s" % (args.prefix, pose, v),
                      preview_png=os.path.join(
                          inter, helper_name(args.prefix, pose, "svgcheck_" + v)),
                      report=lambda m: print("    " + m))
            made.append(dst)
        except Exception as exc:
            print("  FAILED %s SVG: %s" % (v, exc), file=sys.stderr)
    for m in made:
        print("  svg  %s" % m)
    return 0 if made else 1


def cmd_prep(args) -> int:
    """Stage 1: select the atoms the rest of the pipeline measures."""
    from .prep import PrepParams, prepare, write_pdb

    params = PrepParams(
        chains=_split_csv(args.chains), keep_ligands=_split_csv(args.keep),
        drop_ligands=_split_csv(args.drop), keep_waters=args.keep_waters,
        keep_altlocs=args.keep_altlocs, relabel=args.relabel)
    try:
        atoms, report = prepare(args.structure, params)
    except (OSError, ValueError) as exc:
        return _fail(str(exc))
    print(report.summary())
    out = args.out or os.path.splitext(os.path.basename(args.structure))[0] + "_prep.pdb"
    write_pdb(out, atoms, header=[
        "prepared by ChannelBridge from %s" % os.path.basename(args.structure),
        "chains %s%s" % ("+".join(report.chains_kept),
                         " renamed to %s" % report.relabelled_to
                         if report.relabelled_to else ""),
        "ligands kept %s" % (", ".join(sorted(report.ligands_kept)) or "none"),
        "coordinates are unchanged",
    ])
    print("wrote %s" % out)
    return 0


def cmd_channelsdb(args) -> int:
    """Tunnels already computed for a deposited PDB entry."""
    from .channelsdb import fetch, read_report, to_tunnel, tunnels

    try:
        report = (read_report(args.pdb_id) if os.path.isfile(args.pdb_id)
                  else fetch(args.pdb_id, cache_dir=args.cache))
    except LookupError as exc:
        return _fail(str(exc))
    except OSError as exc:
        return _fail("could not reach ChannelsDB: %s" % exc)
    found = tunnels(report, method=args.method) if args.method else tunnels(report)
    if not found:
        return _fail("%s is in ChannelsDB but has no tunnels in the requested "
                     "categories" % args.pdb_id)
    print("%d tunnel(s)" % len(found))
    for t in found:
        print("  %-28s length %6.1f A  bottleneck %.2f A (free %.2f)"
              % (t.label, t.length, t.bottleneck, float(t.free_radii.min())))
    print()
    print("Categories say what a tunnel leads to: Procognate = the cognate")
    print("ligand site (usually the substrate channel), Cofactor = a bound")
    print("cofactor -- which is the relevant one when the shuttled species IS")
    print("the cofactor, as for a dehydrogenase/reductase pair.")
    print("No energetics here: ChannelsDB stores geometry, not CaverDock profiles.")
    if args.validate:
        missing = [k for k in ("structure", "template") if not getattr(args, k)]
        if missing:
            return _fail("--validate also needs --%s" % " and --".join(missing))
        from .figures import Runner

        runner = Runner(pymol=args.pymol, dry_run=False)
        passthrough = ["ours=%s" % args.validate, "structure=%s" % args.structure,
                       "template=%s" % args.template,
                       "report=%s" % (args.pdb_id if os.path.isfile(args.pdb_id)
                                      else _cached_report(args, report)),
                       "method=%s" % args.align_method]
        for key, value in (("select", args.align_select),
                           ("template_select", args.template_select)):
            if value:
                passthrough.append("%s=%s" % (key, value))
        return runner.run(runner.script("run_validate.py"), passthrough)

    if args.write:
        os.makedirs(args.write, exist_ok=True)
        for t in found:
            dst = os.path.join(args.write, "%s_%s.pdb" % (args.pdb_id.lower(), t.label))
            to_tunnel(t).write_caver_pdb(dst)
            print("  wrote %s" % dst)
    return 0


def _cached_report(args, report):
    """Write the fetched report somewhere the PyMOL shim can read it."""
    import json
    import tempfile

    out = os.path.join(args.cache or tempfile.mkdtemp(),
                       "channelsdb_%s.json" % args.pdb_id.lower())
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(out):
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(report, fh)
    return out


def cmd_caver(args) -> int:
    """Stage 2: write a CAVER config, or choose tunnels from its output."""
    from .caver import CaverParams, caverdock_note, read_summary, select, write_config

    if args.action == "config":
        params = CaverParams()
        for key in ("probe_radius", "shell_radius", "shell_depth",
                    "clustering_threshold", "seed"):
            value = getattr(args, key, None)
            if value is not None:
                setattr(params, key, value)
        if args.exclude:
            params.exclude = tuple(_split_csv(args.exclude))
        if args.start_residues:
            params.start_residues = tuple(_split_csv(args.start_residues))
        if args.start_atoms:
            try:
                params.start_atoms = tuple(int(v) for v in _split_csv(args.start_atoms))
            except ValueError:
                return _fail("--start-atoms wants PDB atom serial numbers, got %r"
                             % args.start_atoms)
        point = None
        if args.start:
            try:
                point = [float(v) for v in _split_csv(args.start)]
            except ValueError:
                return _fail("--start wants three numbers, got %r" % args.start)
        if args.start_cavity:
            if not args.structure:
                return _fail("--start-cavity needs --structure to measure against")
            if point is not None:
                return _fail("give --start or --start-cavity, not both: one is a "
                             "point you chose, the other derives one")
            from .caver import cavity_start_point

            try:
                found = cavity_start_point(
                    args.structure, _split_csv(args.start_cavity),
                    probe_radius=params.probe_radius)
            except (OSError, ValueError) as exc:
                return _fail(str(exc))
            point = list(found.xyz)
            print(found.summary())
            if found.burial < 80:
                print("  WARNING burial %d is low for a catalytic site; this "
                      "point may be in bulk solvent rather than a pocket"
                      % found.burial)
        out = args.out or "caver_config.txt"
        try:
            write_config(out, params, point)
        except ValueError as exc:
            return _fail(str(exc))
        print("wrote %s" % out)
        print("  probe %.2f A (CAVER's own default is 0.90, which finds "
              "channels no water can pass)" % params.probe_radius)
        return 0

    if args.action == "run":
        from .caver import run_caver

        for need in ("structure", "config", "output_dir"):
            if not getattr(args, need, None):
                return _fail("'caver run' needs --structure, --config and an "
                             "output directory")
        try:
            out = run_caver(args.structure, args.config, args.output_dir,
                            caver=args.caver, java=args.java, heap=args.heap)
        except (FileNotFoundError, RuntimeError) as exc:
            return _fail(str(exc))
        print("CAVER finished; output in %s" % out)
        args.output_dir = out

    try:
        clusters = read_summary(args.output_dir)
    except (OSError, ValueError) as exc:
        return _fail(str(exc))
    sel = select(clusters, keep=args.keep, floor=args.floor)
    print("%d cluster(s) found; carrying %d forward"
          % (len(clusters), len(sel.kept)))
    print(sel.summary())
    if sel.kept and not args.quiet:
        print()
        print(caverdock_note(sel))
    return 0 if sel.kept else 1


def cmd_patches(args) -> int:
    """Stage 3: classify one enzyme's residues and rank its epitopes."""
    from .figures import Runner

    runner = Runner(pymol=args.pymol, dry_run=args.dry_run)
    out_dir = args.out_dir or "."
    pass_through = ["structure=%s" % args.structure, "chain=%s" % args.chain,
                    "out_dir=%s" % out_dir]
    for key, value in (("name", args.name), ("tunnel", args.tunnel),
                       ("avoid", args.avoid)):
        if value:
            pass_through.append("%s=%s" % (key, value))
    if args.ligands:
        pass_through.append("ligands=%s" % ",".join(_split_csv(args.ligands)))
    # every PatchParams field is overridable on the command line; only the
    # ones actually given are passed, so the defaults stay in one place
    for key in ("ligand_r", "channel_r", "shadow_r", "shadow_back", "exposed",
                "iface_cut", "seed_sasa", "patch_r", "patch_min", "patch_sep",
                "n_patches", "axis_lo", "axis_hi", "reach_lo", "reach_hi",
                "avoid_r", "avoid_min_sasa"):
        value = getattr(args, key, None)
        if value is not None:
            pass_through.append("%s=%s" % (key, value))
    rc = runner.run(runner.script("run_patches.py"), pass_through)
    if rc:
        return _fail("run_patches exited %d" % rc)
    return 0


def cmd_hotspots(args) -> int:
    """The overall hotspot map: every candidate epitope of both enzymes."""
    from .figures import Runner, scan_output

    style = None
    if args.pair_file:
        from .config import ConfigError, load_pair_file

        try:
            style = load_pair_file(args.pair_file).render
        except (ConfigError, OSError) as exc:
            return _fail(str(exc))

    complex_pdb = args.complex
    computed = bool(args.up_pdb and args.dn_pdb)
    if not complex_pdb and not computed:
        try:
            geoms = scan_output(args.output_dir)
        except FileNotFoundError as exc:
            return _fail(str(exc))
        wanted = args.pose or (geoms[0].pose_id if geoms else None)
        match = [g for g in geoms if g.pose_id == wanted]
        if not match:
            return _fail("pose %r is not in %s (have: %s)"
                         % (wanted, args.output_dir,
                            ", ".join(g.pose_id for g in geoms[:8]) or "none"))
        complex_pdb = match[0].complex_pdb
        print("positions from %s" % os.path.basename(complex_pdb))

    # One map per enzyme PAIR, not per pose: <prefix>_hotspot-map.pse.
    out = args.out
    if not out:
        from .figures import hotspot_map_name
        base = args.out_dir or (os.path.join(args.output_dir, "figures")
                                if args.output_dir else ".")
        out = os.path.join(base, hotspot_map_name(args.prefix))
        if not args.dry_run:
            os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        print("writing %s" % out)
    runner = Runner(pymol=args.pymol, dry_run=args.dry_run)
    pass_through = [
        "out=%s" % out,
        "up_chains=%s" % ",".join(_split_csv(args.up_chains)),
        "dn_chains=%s" % ",".join(_split_csv(args.dn_chains)),
        "up_name=%s" % args.up_name, "dn_name=%s" % args.dn_name,
        "binder_radius=%s" % args.binder_radius,
    ]
    for key, value in (("complex", complex_pdb),
                       ("up_pdb", args.up_pdb), ("dn_pdb", args.dn_pdb),
                       ("tunnel_up", args.tunnel_up), ("tunnel_dn", args.tunnel_dn),
                       ("classes_up", args.classes_up), ("classes_dn", args.classes_dn),
                       ("patches_up", args.patches_up), ("patches_dn", args.patches_dn)):
        if value:
            pass_through.append("%s=%s" % (key, value))
    if computed:
        pass_through += [
            "exclude=%s" % ",".join(_split_csv(args.exclude)),
            "binder_length=%s" % args.binder_length,
            "binder_step=%s" % args.binder_step,
            "d_min=%s" % args.d_min, "d_max=%s" % args.d_max, "d_step=%s" % args.d_step,
            "clash_cut=%s" % args.clash_cut,
        ]
    if args.ligands:
        pass_through.append("ligands=%s" % ",".join(_split_csv(args.ligands)))
    if style is not None and getattr(style, "hotspots", None):
        # the palette is read inside PyMOL, where the pair file is not
        # available, so the one setting that changes the drawing order of the
        # shades is passed explicitly and the rest come from the defaults
        if style.hotspots.get("best_is_darkest"):
            pass_through.append("best_is_darkest=1")
    rc = runner.run(runner.script("run_hotspots.py"), pass_through)
    if rc:
        return _fail("the hotspot map exited %d" % rc)
    print("wrote %s" % out)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="chbridge",
        description="Channel-aware pose search for bispecific binder design.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version="chbridge %s" % __version__)
    sub = p.add_subparsers(dest="command", metavar="<command>")

    v = sub.add_parser("validate", help="check a pair file and print what it says")
    v.add_argument("pair_file", help="a .yaml/.yml/.json pair file")
    v.set_defaults(func=cmd_validate)

    t = sub.add_parser(
        "tunnel",
        help="normalise CAVER tunnel cluster(s) and print the geometry",
        description=(
            "Reads one or more CAVER output folders (or tun_cl_*.pdb files) and "
            "reports the tunnel geometry raw and normalised. With two or more "
            "runs and --compare it reports how far apart their exits are, which "
            "is how to tell a sphere-sampling difference from a different tunnel."
        ),
    )
    t.add_argument("caver_out", nargs="+", help="out_<tag> folder(s), or a tun_cl_*.pdb")
    t.add_argument(
        "--structure",
        help="the PDB CAVER searched. Required to orient and truncate by burial.",
    )
    t.add_argument(
        "--exclude",
        default="HOH,WAT",
        help="comma-separated residue names CAVER excluded (default: HOH,WAT). "
        "Must match the config's exclude_residue_names or burial is counted "
        "against atoms CAVER never searched.",
    )
    t.add_argument("--cluster", type=int, help="cluster number (default: the first found)")
    t.add_argument("--spacing", type=float, default=1.0, help="resample spacing, A (default 1.0)")
    t.add_argument("--burial-min", type=int, default=15, help="exit burial threshold (default 15)")
    t.add_argument("--burial-radius", type=float, default=8.0, help="burial radius, A (default 8.0)")
    t.add_argument("--oic-arm", type=float, default=1.0, help="O_IC arm, A (default 1.0)")
    t.add_argument("--raw-only", action="store_true", help="skip normalisation")
    t.add_argument(
        "--out-pdb",
        help="write the normalised chain here in CAVER tunnel format, to hand "
        "to a tool that reads CAVER output directly. Needs --structure and a "
        "single input.",
    )
    t.add_argument("--compare", action="store_true", help="report exit separation between runs")
    t.set_defaults(func=cmd_tunnel)

    f = sub.add_parser(
        "figures",
        help="sessions, SVGs and structures for chosen poses",
        description=(
            "Picks poses out of a pose-search output folder and builds figures "
            "for them. Each of --pse/--svg/--pdb takes the same selection "
            "language: a count (5), a top-N (top5), a rank (rank3, #3), a "
            "range (rank2-4), an id (geom007), 'all', 'none', or a "
            "comma-separated mixture. An SVG implies its session, and the "
            "ray-traces are run one at a time."
        ),
    )
    f.add_argument("output_dir", help="a pose-search output folder (holds geom*_complex.pdb)")
    f.add_argument("--prefix", default="fig",
                   help="goes at the front of every file this step writes (e.g. run1)")
    f.add_argument("--pse", help="which poses get a PyMOL session")
    f.add_argument("--svg", help="which poses get an editable SVG")
    f.add_argument("--pdb", help="which poses get their complex copied out")
    f.add_argument(
        "--render", default="pair",
        help="which render variant the SVG shows: none, pair, assembly, both "
        "(default: pair)",
    )
    f.add_argument("--out-dir", help="where the figures go (default: <output_dir>/figures)")
    f.add_argument("--pair-file", help="a pair file, for its render: style section")
    f.add_argument("--width", type=int, default=1200, help="render width, px (default 1200)")
    f.add_argument("--height", type=int, default=900, help="render height, px (default 900)")
    f.add_argument("--pymol", default=None, help="the pymol binary to use")
    f.add_argument(
        "--pse-arg", action="append", default=[],
        help="extra key=value passed through to the session builder; repeatable",
    )
    f.add_argument("--rebuild", action="store_true", help="rebuild sessions that exist")
    f.add_argument("--dry-run", action="store_true", help="print the steps, run nothing")
    f.set_defaults(func=cmd_figures)

    hs = sub.add_parser(
        "hotspots",
        help="one session showing every candidate epitope of both enzymes",
        description=(
            "Builds the overall hotspot map: the residue classification "
            "(surface, avoid, shadow, channel) painted the SAME way on both "
            "enzymes, every patch as a shade of one hue and labelled, and a "
            "binder envelope at each patch. The enzymes stand in a real pose's "
            "relative position (the best pose by default)."
        ),
    )
    hs.add_argument("output_dir", nargs="?", help="a pose-search output folder")
    hs.add_argument("--prefix", default="fig",
                   help="goes at the front of every file this step writes (e.g. run1)")
    hs.add_argument("--up-pdb", help="prepared upstream enzyme (the pre-sampling route)")
    hs.add_argument("--dn-pdb", help="prepared downstream enzyme")
    hs.add_argument("--tunnel-up", help="normalised upstream tunnel chain")
    hs.add_argument("--tunnel-dn", help="normalised downstream tunnel chain")
    hs.add_argument("--exclude", default="HOH,WAT",
                    help="residue names CAVER excluded (default HOH,WAT)")
    hs.add_argument("--binder-length", type=float, default=6.0,
                    help="shortest binder worth drawing, A (default 6, the pose "
                         "stage's binder gap)")
    hs.add_argument("--binder-step", type=float, default=2.0,
                    help="spacing of the envelope blobs, A (default 2)")
    hs.add_argument("--d-min", type=float, default=6.0, help="smallest D_IC to try, A")
    hs.add_argument("--d-max", type=float, default=36.0, help="largest D_IC to try, A")
    hs.add_argument("--d-step", type=float, default=1.0, help="D_IC step, A")
    hs.add_argument("--clash-cut", type=float, default=3.0,
                    help="closest allowed interprotein approach, A (default 3)")
    hs.add_argument("--complex", help="use this complex PDB instead of computing a placement")
    hs.add_argument("--pose", help="which pose fixes the positions (default: the best)")
    hs.add_argument("--out",
                    help="where the session goes "
                         "(default: <out-dir>/<prefix>_hotspot-map.pse)")
    hs.add_argument("--out-dir",
                    help="folder for the default name (default: the current "
                         "directory, or <output_dir>/figures if one is given)")
    hs.add_argument("--up-chains", default="", help="upstream chains, comma-separated")
    hs.add_argument("--dn-chains", default="", help="downstream chains, comma-separated")
    hs.add_argument("--up-name", default="up", help="name for the upstream enzyme's objects")
    hs.add_argument("--dn-name", default="dn", help="name for the downstream enzyme's objects")
    hs.add_argument("--classes-up", help="<NAME>_residue_classes.tsv for the upstream enzyme")
    hs.add_argument("--classes-dn", help="<NAME>_residue_classes.tsv for the downstream enzyme")
    hs.add_argument("--patches-up", help="<NAME>_hotspot_patches.tsv for the upstream enzyme")
    hs.add_argument("--patches-dn", help="<NAME>_hotspot_patches.tsv for the downstream enzyme")
    hs.add_argument("--ligands", help="ligand residue names to show as sticks")
    hs.add_argument("--binder-radius", type=float, default=11.0,
                    help="binder envelope radius, A (default 11, the pose stage's)")
    hs.add_argument("--pair-file", help="a pair file, for its render.hotspots: section")
    hs.add_argument("--pymol", default=None, help="the pymol binary to use")
    hs.add_argument("--dry-run", action="store_true", help="print the step, run nothing")
    hs.set_defaults(func=cmd_hotspots)

    se = sub.add_parser("session", help="build the .pse for chosen poses (no rendering)")
    se.add_argument("output_dir", help="a pose-search output folder")
    se.add_argument("--prefix", default="fig", help="front of every file written (e.g. run1)")
    se.add_argument("--pse", help="which poses (default top1); same selection language")
    se.add_argument("--render", default="both", help="none, pair, assembly or both")
    se.add_argument("--out-dir", help="where sessions go (default <output_dir>/figures)")
    se.add_argument("--pse-arg", action="append", default=[],
                    help="extra key=value for the session builder; repeatable")
    se.add_argument("--pymol", default=None)
    se.add_argument("--dry-run", action="store_true")
    se.set_defaults(func=cmd_session)

    rn = sub.add_parser(
        "render", help="render scenes and assemble SVGs from an existing .pse",
        description="Takes any session -- ours, or one you rotated and saved "
                    "yourself -- and writes <prefix>_<pose>_render_<variant>.svg.")
    rn.add_argument("session", help="the .pse to render")
    rn.add_argument("--prefix", default="fig", help="front of every file written")
    rn.add_argument("--pose", default="geom001", help="pose id used in the names")
    rn.add_argument("--variant", default="both", help="pair, assembly or both")
    rn.add_argument("--dic", type=float, default=0.0, help="D_IC to label, A")
    rn.add_argument("--oic", type=float, default=0.0, help="O_IC to label, degrees")
    rn.add_argument("--do", type=float, default=0.0, help="DO to label")
    rn.add_argument(
        "--axis-px", default=None,
        help="the projected P2->P1 direction as 'x,y' (image y downward), "
             "which the session reports. REQUIRED for a rolled frame: "
             "without it the angle marker falls back to assuming the "
             "P1-P2-P3 plane is face-on, and draws the arm at an angle the "
             "picture does not contain")
    rn.add_argument("--out-dir", help="default: beside the session")
    rn.add_argument("--width", type=int, default=1200)
    rn.add_argument("--height", type=int, default=900)
    rn.add_argument("--pymol", default=None)
    rn.add_argument("--dry-run", action="store_true")
    rn.set_defaults(func=cmd_render)

    pr = sub.add_parser(
        "prep", help="stage 1: select chains, ligands and altlocs",
        description="Keeps the atoms the rest of the pipeline measures and "
                    "drops the rest. Coordinates are never moved -- not "
                    "re-centred, not symmetry-expanded. If you took a tunnel "
                    "from ChannelsDB it is in the deposited entry's frame, "
                    "and moving the protein silently invalidates it.")
    pr.add_argument("structure", help="a deposited or modelled structure, PDB")
    pr.add_argument("--chains", help="chains to keep, comma separated "
                                     "(default: the first polymer chain)")
    pr.add_argument("--keep", help="ligand resnames to KEEP: cofactors, "
                                   "prosthetic groups, and the channelled "
                                   "species (whose site gives the start point)")
    pr.add_argument("--drop", help="extra resnames to drop")
    pr.add_argument("--keep-waters", dest="keep_waters", action="store_true",
                    help="rarely right: a water in a channel is a wall")
    pr.add_argument("--keep-altlocs", dest="keep_altlocs", action="store_true",
                    help="keep every alternate location rather than the "
                         "highest-occupancy one. Two copies of a side chain "
                         "narrow a channel that is not narrow")
    pr.add_argument("--relabel", help="rename the kept chain to this id. A "
                                     "pair built from two single-chain "
                                     "entries has both enzymes on chain A, "
                                     "and the pose stage cannot tell them "
                                     "apart; renaming a chain does not move "
                                     "it, so this is safe on the ChannelsDB "
                                     "branch too")
    pr.add_argument("--out", help="output PDB (default: <name>_prep.pdb)")
    pr.set_defaults(func=cmd_prep)

    cd = sub.add_parser(
        "channelsdb", help="tunnels already computed for a deposited PDB entry",
        description="Fetches precomputed tunnels from ChannelsDB 2.0, by both "
                    "MOLE and CAVER, annotated by what they lead to. Faster "
                    "than running CAVER and needs no CAVER install -- but "
                    "coverage is per PDB ENTRY, not per protein, and the "
                    "coordinates are that entry's.")
    cd.add_argument("pdb_id", help="a 4-character PDB id, or a report.json on disk")
    cd.add_argument("--method", choices=("Caver", "MOLE"),
                    help="keep only one method's tunnels (default: both, so "
                         "you can see where they disagree)")
    cd.add_argument("--cache", help="directory to cache downloaded reports in")
    cd.add_argument("--write", help="write each tunnel as a CAVER-format PDB here")
    cd.add_argument("--validate", metavar="TUNNEL.PDB",
                    help="one of YOUR CAVER clusters: superpose the database's "
                         "structure onto yours and report whether any of its "
                         "tunnels is the same tunnel")
    cd.add_argument("--structure", help="your structure, for --validate")
    cd.add_argument("--template", help="the database entry's structure, for --validate")
    cd.add_argument("--align-method", dest="align_method", default="align",
                    choices=("align", "super", "cealign"),
                    help="align = sequence-based (default, needs similarity); "
                         "super = sequence-independent; cealign = structural, "
                         "best for remote homologs")
    cd.add_argument("--align-select", dest="align_select",
                    help="PyMOL selection in YOUR structure to fit on. Fitting "
                         "on the site beats fitting on the whole fold: a "
                         "whole-chain fit answers 'same fold', a tunnel "
                         "comparison needs 'same site'")
    cd.add_argument("--template-select", dest="template_select",
                    help="the matching selection in the template -- restrict it "
                         "to one chain when your structure is a fusion")
    cd.add_argument("--pymol", default=None)
    cd.set_defaults(func=cmd_channelsdb)

    cv = sub.add_parser(
        "caver", help="stage 2: configure CAVER, or pick tunnels from its output",
        description="`caver config` writes a CAVER configuration with "
                    "ChannelBridge's defaults (a 1.4 A water probe, not "
                    "CAVER's 0.9) and your overrides. `caver select` reads a "
                    "finished CAVER output folder and applies the published "
                    "rule: keep the top three clusters above priority 0.55. "
                    "Then proceed with the highest-priority one, or validate "
                    "the survivors with CaverDock.")
    cv.add_argument("action", choices=("config", "run", "select"))
    cv.add_argument("output_dir", nargs="?",
                    help="for 'select': a CAVER output folder (or its summary.txt)")
    cv.add_argument("--out", help="for 'config': where the file goes")
    cv.add_argument("--start", help="for 'config': starting point as 'x,y,z'")
    cv.add_argument("--start-residues", dest="start_residues",
                    help="for 'config': CAVER residue ids giving the starting "
                         "point, comma separated ('B:128,B:303'). Use this "
                         "when no ligand marks the site. Name the chain: "
                         "without one CAVER silently assumes the first")
    cv.add_argument("--start-atoms", dest="start_atoms",
                    help="for 'config': PDB atom serial numbers, comma "
                         "separated. Combines with the other two")
    cv.add_argument("--start-cavity", dest="start_cavity",
                    help="for 'config': DERIVE the starting point from these "
                         "residues (needs --structure) -- the nearest point to "
                         "their centroid that can hold the probe. Use this "
                         "rather than --start-residues for a catalytic triad: "
                         "the centroid of residues facing each other lands "
                         "inside an atom, and CAVER then finds zero tunnels "
                         "without saying why")
    cv.add_argument("--probe-radius", dest="probe_radius", type=float, default=None,
                    help="default 1.4 (water). CAVER's own default is 0.9, "
                         "which roughly doubles the cluster count and leaves "
                         "the top-ranked tunnel unchanged")
    cv.add_argument("--shell-radius", dest="shell_radius", type=float, default=None)
    cv.add_argument("--shell-depth", dest="shell_depth", type=float, default=None)
    cv.add_argument("--clustering-threshold", dest="clustering_threshold",
                    type=float, default=None)
    cv.add_argument("--seed", type=int, default=None,
                    help="fixed by default: CAVER's search is stochastic and "
                         "the cluster NUMBERING can change without it")
    cv.add_argument("--exclude", help="residue names to remove, comma separated")
    cv.add_argument("--keep", type=int, default=3,
                    help="how many clusters to carry forward (default 3)")
    cv.add_argument("--floor", type=float, default=0.55,
                    help="minimum CAVER priority (default 0.55, published)")
    cv.add_argument("--structure", help="for 'run': the prepared structure")
    cv.add_argument("--config", help="for 'run': the CAVER config to use")
    cv.add_argument("--caver", help="CAVER directory (default: $CAVER)")
    cv.add_argument("--java", default="java")
    cv.add_argument("--heap", default="6g", help="JVM heap (default 6g)")
    cv.add_argument("--quiet", action="store_true",
                    help="omit the note on what CaverDock would and would not settle")
    cv.set_defaults(func=cmd_caver)

    pa = sub.add_parser(
        "patches", help="stage 3: classify residues and rank the epitopes",
        description="Writes <name>_residue_classes.tsv and "
                    "<name>_hotspot_patches.tsv, which the hotspot map and "
                    "the pose search both read. Edit the patch table by hand "
                    "or supply your own: nothing downstream cares who wrote "
                    "it, which is the point of the stage being optional.")
    pa.add_argument("structure", help="one prepared enzyme, as PDB")
    pa.add_argument("--chain", default="A", help="which chain to classify")
    pa.add_argument("--name", help="prefix for the two tables (default: the file stem)")
    pa.add_argument("--tunnel", help="its normalised CAVER tunnel; without one "
                                     "there is no mouth, so no CHANNEL or SHADOW class")
    pa.add_argument("--ligands", help="cofactor resnames, comma separated")
    pa.add_argument("--avoid", help="PyMOL selection(s) to keep binders away "
                                    "from, ';' separated (e.g. \"resi 145+158\")")
    pa.add_argument("--out-dir", help="where the two tables go (default: here)")
    for key, kind, helptext in (
            ("ligand_r", float, "within this of a cofactor is LIGAND_CORE"),
            ("channel_r", float, "within this of the tunnel is CHANNEL"),
            ("shadow_r", float, "radius of the no-block cylinder at the mouth"),
            ("shadow_back", float, "how far inside the mouth it still reaches"),
            ("exposed", float, "SASA at which a residue counts as exposed"),
            ("iface_cut", float, "dSASA above which a residue is IFACE"),
            ("seed_sasa", float, "a seed anchor needs this much exposed area"),
            ("patch_r", float, "patch members lie within this of the seed"),
            ("patch_min", int, "fewer members than this is not a patch"),
            ("patch_sep", float, "minimum centroid separation of kept patches"),
            ("n_patches", int, "how many patches to keep"),
            ("axis_lo", float, "best band from the channel axis, low end"),
            ("axis_hi", float, "best band from the channel axis, high end"),
            ("reach_lo", float, "best band from the mouth, low end"),
            ("reach_hi", float, "best band from the mouth, high end"),
            ("avoid_r", float, "how far an exposed avoided residue spreads"),
            ("avoid_min_sasa", float, "below this it protects only itself")):
        pa.add_argument("--" + key.replace("_", "-"), dest=key, type=kind,
                        default=None, help=helptext)
    pa.add_argument("--pymol", default=None)
    pa.add_argument("--dry-run", action="store_true")
    pa.set_defaults(func=cmd_patches)

    d = sub.add_parser("do", help="the DO score for one D_IC and O_IC")
    d.add_argument("--dic", type=float, required=True, help="inter-channel distance, A")
    d.add_argument("--oic", type=float, required=True, help="inter-channel orientation, degrees")
    d.add_argument(
        "--coefficients",
        default="imars-fig2c",
        choices=["imars-fig2c", "imars-exact", "imars-st2a"],
        help="which constant set: imars-fig2c (the paper's stated values, default), imars-exact (the same derivation unrounded; reproduces their Table S2F), imars-st2a (reproduces their Table S2A, which uses a set that does not derive from their own regression)",
    )
    d.add_argument("--all-sets", action="store_true", help="also show the other set")
    d.set_defaults(func=cmd_do)

    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    try:
        return int(args.func(args) or 0)
    except KeyboardInterrupt:  # pragma: no cover
        return 130
    except BrokenPipeError:  # pragma: no cover - `| head`
        return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
