"""Turn a pose-search output folder into sessions, SVGs and structures.

The pose search writes one ``geom<NNN>_...`` set per surviving pose. Which of
them become figures is a separate question from the search, asked long after
it: hence a separate command rather than an option of the search.

What this module adds over calling the PyMOL scripts by hand:

* it finds the poses and their files by scanning the output folder, so a pose
  is named by its id (``geom007``) or its rank (``top3``) rather than by a
  90-character file name;
* it knows that an SVG needs a session, and that a session needs renders, so
  ``--svg top3`` is enough;
* it runs the ray-traces **one at a time**. They are minutes to tens of
  minutes each on a shared machine, and two at once is slower than two in a
  row as well as unfair to whoever else is on the node.

The PyMOL steps are subprocesses, not imports: PyMOL owns its interpreter,
and a figure run should not be able to take the rest of the tool down with it.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

__all__ = [
    "DEFAULT_PYMOL",
    "SCENES_FOR_SVG",
    "GeomFiles",
    "collect_poses",
    "default_pymol",
    "helper_renders",
    "scan_output",
    "svg_scene_name",
]

#: Where to find PyMOL, in order: the ``CHBRIDGE_PYMOL`` environment variable,
#: whatever ``pymol`` is on PATH, then the path it lives at on the cluster this
#: was developed on. Hard-coding that last path alone made the tool unusable
#: anywhere else -- including inside its own dev container, where PyMOL comes
#: from the distribution and lands in /usr/bin.
_CLUSTER_PYMOL = "/opt/conda/envs-share/pymol-script/bin/pymol"


def default_pymol() -> str:
    env = os.environ.get("CHBRIDGE_PYMOL")
    if env:
        return env
    found = shutil.which("pymol")
    if found:
        return found
    return _CLUSTER_PYMOL


#: Kept as a module constant for callers that want it eagerly; prefer
#: :func:`default_pymol`, which re-reads the environment each call.
DEFAULT_PYMOL = default_pymol()

#: Flat helper renders every SVG of a frame shares. They are cheap (minutes)
#: and identical for every render variant of the same pose, so they are made
#: once per pose rather than once per SVG.
SCENES_FOR_SVG = ("mask", "mouth", "footprint")

_GEOM = re.compile(r"^(geom\d+)")


@dataclass
class GeomFiles:
    """One pose's files, as found on disk."""

    pose_id: str
    stem: str                                  #: the long name without a suffix
    directory: str
    do: Optional[float] = None
    complex_pdb: Optional[str] = None
    occupancy_pdb: Optional[str] = None
    pse: Optional[str] = None
    d_ic: Optional[float] = None
    o_ic: Optional[float] = None
    patch_up: str = ""
    patch_dn: str = ""

    @property
    def values(self) -> Dict[str, float]:
        out = {}
        for key, v in (("D_IC", self.d_ic), ("O_IC", self.o_ic), ("DO", self.do)):
            if v is not None:
                out[key] = float(v)
        return out


def _from_stem(stem: str) -> Dict[str, object]:
    """Pull what the pose stage encoded in the file name.

    The name is the only place some of this is written per pose (the tables
    hold it too, but a figure run should work on a folder someone copied
    without them), so it is parsed rather than required.
    """
    out: Dict[str, object] = {}
    m = re.search(r"_DO([0-9.]+)", stem)
    if m:
        out["do"] = float(m.group(1))
    m = re.search(r"_d(\d+(?:\.\d+)?)", stem)
    if m:
        out["d_ic"] = float(m.group(1))
    m = re.search(r"_o(\d+(?:\.\d+)?)", stem)
    if m:
        out["o_ic"] = float(m.group(1))
    m = re.search(r"_([A-Za-z0-9]+)p(\d+)-([A-Za-z0-9]+)p(\d+)_", stem)
    if m:
        out["patch_up"] = m.group(2)
        out["patch_dn"] = m.group(4)
    return out


def scan_output(directory) -> List[GeomFiles]:
    """Every pose in a sampling output folder, in file order.

    A pose counts as present when it has a ``_complex.pdb``: that is what
    every downstream step reads, and a folder can hold a ``.pml`` for a pose
    whose structure was never written.
    """
    directory = str(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError("no such output folder: %s" % directory)
    found: Dict[str, GeomFiles] = {}
    for name in sorted(os.listdir(directory)):
        m = _GEOM.match(name)
        if not m:
            continue
        pid = m.group(1)
        path = os.path.join(directory, name)
        if name.endswith("_complex.pdb"):
            stem = name[: -len("_complex.pdb")]
            g = found.setdefault(pid, GeomFiles(pid, stem, directory))
            g.stem = stem
            g.complex_pdb = path
            for k, v in _from_stem(stem).items():
                setattr(g, k, v)
        elif name.endswith("_occupancy.pdb"):
            g = found.setdefault(pid, GeomFiles(pid, name[: -len("_occupancy.pdb")], directory))
            g.occupancy_pdb = path
        elif name.endswith(".pse"):
            g = found.setdefault(pid, GeomFiles(pid, name[: -len(".pse")], directory))
            g.pse = path
    return [g for _, g in sorted(found.items()) if g.complex_pdb]


def collect_poses(geoms: Sequence[GeomFiles]):
    """``(pose_id, DO)`` pairs for the selection language."""
    return [(g.pose_id, g.do if g.do is not None else float("-inf")) for g in geoms]


def svg_scene_name(variant: str, cross_section: bool = False) -> str:
    """The scene an SVG of this variant is rendered from."""
    base = "render_cs" if cross_section else "render"
    return "%s_%s" % (base, variant) if not cross_section else base


def helper_renders(out_dir: str, pose_id: str, width: int, height: int) -> Dict[str, str]:
    return {
        scene: os.path.join(out_dir, "intermediate",
                            "%s_%s_bgFFFFFF.png" % (pose_id, scene))
        for scene in SCENES_FOR_SVG
    }


# ------------------------------------------------------------------- running


@dataclass
class Runner:
    """Runs the PyMOL steps, one at a time, and says what it is doing."""

    #: ``None`` means the default, so a caller can pass an unset flag straight
    #: through without having to know what the default is.
    pymol: Optional[str] = None
    dry_run: bool = False
    report: object = print
    package_root: str = field(default_factory=lambda: os.path.dirname(os.path.abspath(__file__)))

    @property
    def binary(self) -> str:
        return self.pymol or default_pymol()

    def script(self, name: str) -> str:
        return os.path.join(self.package_root, "render", name)

    def run(self, script: str, args: Sequence[str]) -> int:
        cmd = [self.binary, "-cq", script, "--", *args]
        self.report("    $ %s" % " ".join(cmd))
        if self.dry_run:
            return 0
        if not os.path.exists(self.binary):
            self.report("    PyMOL not found at %s; skipping" % self.binary)
            return 127
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True)
        for line in (proc.stdout or "").splitlines():
            self.report("      " + line)
        return proc.returncode

    def copy(self, src: str, dst: str) -> None:
        self.report("    copy %s -> %s" % (os.path.basename(src), dst))
        if not self.dry_run:
            os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
            shutil.copy2(src, dst)


def _ensure_dirs(out_dir: str, dry_run: bool) -> None:
    if dry_run:
        return
    for sub in ("", "intermediate"):
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)


def build(plan, geoms: Sequence[GeomFiles], out_dir: str, *, runner: Optional[Runner] = None,
          width: int = 1200, height: int = 900, pse_args: Sequence[str] = (),
          style=None, svg_style=None, report=print) -> Dict[str, List[str]]:
    """Carry out a :class:`~channelbridge.render.select.FigurePlan`.

    Returns what was written, by kind. Nothing here fails the whole run: a
    pose whose render fails is reported and the rest continue, because a
    figure run of twenty poses should not be lost to one of them.
    """
    runner = runner or Runner(report=report)
    by_id = {g.pose_id: g for g in geoms}
    made: Dict[str, List[str]] = {"pse": [], "svg": [], "pdb": [], "failed": []}
    _ensure_dirs(out_dir, runner.dry_run)

    for pid in plan.pdb:
        g = by_id.get(pid)
        if g is None or not g.complex_pdb:
            made["failed"].append("%s: no complex PDB" % pid)
            continue
        dst = os.path.join(out_dir, "%s_complex.pdb" % g.stem)
        runner.copy(g.complex_pdb, dst)
        made["pdb"].append(dst)

    for pid in plan.pse:
        g = by_id.get(pid)
        if g is None:
            made["failed"].append("%s: not in the output folder" % pid)
            continue
        dst = g.pse or os.path.join(out_dir, "%s.pse" % g.stem)
        report("  [%s] session -> %s" % (pid, dst))
        rc = runner.run(runner.script("run_pse.py"),
                        ["complex=%s" % g.complex_pdb, "out=%s" % dst, *pse_args])
        if rc:
            made["failed"].append("%s: run_pse exited %d" % (pid, rc))
        else:
            made["pse"].append(dst)
            g.pse = dst

    for pid in plan.svg:
        g = by_id.get(pid)
        if g is None or not g.pse:
            made["failed"].append("%s: no session to render from" % pid)
            continue
        report("  [%s] SVG" % pid)
        helpers = helper_renders(out_dir, pid, width, height)
        ok = True
        for scene, path in helpers.items():
            if os.path.exists(path) and not runner.dry_run:
                report("    [have] %s" % os.path.basename(path))
                continue
            rc = runner.run(runner.script("run_png.py"),
                            ["pse=%s" % g.pse, "scene=%s" % scene, "out=%s" % path,
                             "w=%d" % width, "h=%d" % height])
            if rc:
                made["failed"].append("%s: %s render exited %d" % (pid, scene, rc))
                ok = False
                break
        if not ok:
            continue
        for scene in plan_scene_names(plan):
            # the SVG is assembled on the base raster of that variant, NOT on
            # the render scene: the render scene already carries the 3D
            # channel and the 3D D_IC/O_IC/DO annotation, and the SVG adds
            # both as vectors, so building on it prints every label twice.
            raster_scene = _svg_raster_scene(scene)
            base = os.path.join(out_dir, "intermediate",
                                "%s_%s_transparent.png" % (pid, raster_scene))
            if not (os.path.exists(base) and not runner.dry_run):
                rc = runner.run(runner.script("run_png.py"),
                                ["pse=%s" % g.pse, "scene=%s" % raster_scene,
                                 "out=%s" % base, "w=%d" % width, "h=%d" % height,
                                 "transparent=1"])
                if rc:
                    made["failed"].append("%s: %s render exited %d" % (pid, scene, rc))
                    continue
            dst = os.path.join(out_dir, "%s_%s_%s.svg" % (g.stem, scene, "rotated0"))
            if runner.dry_run:
                report("    would assemble %s" % dst)
                made["svg"].append(dst)
                continue
            from .render.svg import build_svg

            try:
                build_svg(base, helpers["mask"], helpers["mouth"], helpers["footprint"],
                          dst, g.values, style=style, svg_style=svg_style,
                          title="%s %s" % (g.pose_id, scene), report=lambda m: report("    " + m))
                made["svg"].append(dst)
            except Exception as exc:
                made["failed"].append("%s: %s SVG failed: %s" % (pid, scene, exc))
    return made


def _svg_raster_scene(render_scene: str) -> str:
    """The annotation-free scene an SVG of *render_scene* is built on."""
    if render_scene == "render_cs":
        return "base_cs"
    return render_scene.replace("render_", "base_")


def plan_scene_names(plan) -> List[str]:
    """The render scenes an SVG is wanted for, from the plan's variants."""
    names: List[str] = []
    for variant in ("monomer", "multimer"):
        if plan.render in (variant, "both"):
            names.append("render_" + variant)
    if plan.cross_section != "none":
        names.append("render_cs")
    return names or ["render_monomer"]


def main(argv=None) -> int:  # pragma: no cover - thin
    from .cli import main as cli_main

    return cli_main(["figures", *(argv if argv is not None else sys.argv[1:])])
