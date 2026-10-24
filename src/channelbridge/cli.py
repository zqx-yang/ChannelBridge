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
  chbridge figures out-complex_GDH_FNR --svg top3 --render both --cross-section monomer
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
    d = caver_out / "data" / "clusters"
    if not d.is_dir():
        raise FileNotFoundError(
            "%s has no data/clusters/. Either CAVER found no pathway (check the "
            "starting point is not inside an atom, and read warnings.txt -- CAVER "
            "exits 0 either way), or this is not a CAVER output folder." % caver_out
        )
    found = sorted(d.glob("tun_cl_*.pdb"))
    if not found:
        raise FileNotFoundError("%s holds no tun_cl_*.pdb" % d)
    if cluster is None:
        return found[0]
    want = d / ("tun_cl_%03d.pdb" % cluster)
    if not want.is_file():
        raise FileNotFoundError(
            "%s not found; this run has clusters %s"
            % (want.name, ", ".join(p.stem.split("_")[-1] for p in found))
        )
    return want


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
            render=args.render, cross_section=args.cross_section,
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
                 pse_args=_split_csv(args.pse_arg), style=style, svg_style=svg_style)
    for kind in ("pse", "svg", "pdb"):
        for path in made[kind]:
            print("  %-4s %s" % (kind, path))
    for bad in made["failed"]:
        print("  FAILED %s" % bad, file=sys.stderr)
    return 1 if made["failed"] else 0


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

    if not args.out:
        return _fail("--out is required: say where the session goes")
    out = args.out
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
    f.add_argument("--pse", help="which poses get a PyMOL session")
    f.add_argument("--svg", help="which poses get an editable SVG")
    f.add_argument("--pdb", help="which poses get their complex copied out")
    f.add_argument(
        "--render", default="monomer",
        help="which render variant the SVG shows: none, monomer, multimer, both "
        "(default: monomer)",
    )
    f.add_argument(
        "--cross-section", default="none",
        help="also assemble the cross-section SVG: none, monomer, multimer, both "
        "(default: none)",
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
    hs.add_argument("--out", help="where the session goes")
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
