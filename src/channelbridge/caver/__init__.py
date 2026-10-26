"""Stage 2: run CAVER, then choose which tunnels to carry forward.

Two jobs, and the second is the one with a literature answer.

**Running it.** CAVER needs a configuration file and a starting point. The
defaults here are ours, chosen for a substrate channel rather than for
CAVER's own general-purpose defaults, and every one is overridable. The
difference matters: CAVER defaults to a 0.9 A probe, which threads gaps no
water could pass and roughly doubles the number of candidate clusters, while
changing nothing about the top-ranked one. Measured on our own data -- the
same six domains at probe 0.9 and 1.4 gave cluster 1 an identical span,
bottleneck and exit at every shell radius, and 51 against 33 clusters.

**Choosing.** "Take the top-ranked tunnel" is what iMARS does and it is not
safe: across our pilot domains CAVER's priority for cluster 2 was 91.6-99.9 %
of cluster 1's, a near-tie that a different seed or a 0.1 A side-chain shift
could swap -- and the swap moves D_IC by tens of angstroms.

Measured on the worked example, to show how little it takes. The same GDH
structure, the same config and the same seed, run twice with the start point
derived two ways -- the centroid of the bound glucose, and the centroid of
the six protein residues lining it, which sit 4.05 A apart. Clusters 1 and 2
swapped. They were tied at priority 0.773, so nothing in the ranking said
which was which, and the "top tunnel" exit moved 10.02 A, worth 0.08 in DO.
The route itself was not lost: cluster 2 of the second run exits 0.00 A from
cluster 1 of the first. Only its rank changed.

Two things follow. A start point taken from lining residues is NOT
interchangeable with one taken from a ligand, so record which you used. And
a tie in priority is not a tie in consequence -- :func:`select` reports the
margin for this reason.

The rule implemented here is published and large. Vavra et al. (2024) ran
CaverDock on ~13,158 experimental enzyme structures with cognate-like ligands
and ~29,752 energy profiles, and found the top-ranked tunnel is the
energetically most favourable one in **75 %** of cases, the second in a
further **19 %** -- over 93 % in the top two. They recommend inspecting
tunnels with priority >= 0.55 and, absent prior information, considering the
first three.

So: keep the first three above the priority floor, and say plainly that the
first is a 75 % bet rather than a certainty. A user who needs better than
that runs CaverDock on the survivors (stage 2b), which is the only thing that
distinguishes them physically -- and should read
:func:`channelbridge.caver.caverdock_note` first, because the method's own
run-to-run spread can exceed the differences it is being asked to resolve.

One thing their paper rules OUT, which we had proposed and withdrew: gating
on bottleneck radius. *"the energy maximum was not in the geometrical
bottleneck in 50% of analysed tunnels."* Geometry finds the true barrier
half the time, so a bottleneck gate is a coin toss dressed as a criterion.

Vavra O, Tyzack J, Haddadi F, Stourac J, Damborsky J, Mazurenko S,
Thornton JM, Bednar D. "Large-scale annotation of biochemically relevant
pockets and tunnels in cognate enzyme-ligand complexes." J. Cheminform.
16:114 (2024). doi:10.1186/s13321-024-00907-z
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    "CaverParams",
    "CavityPoint",
    "Cluster",
    "Selection",
    "caverdock_note",
    "cavity_start_point",
    "find_caver",
    "parse_residue_id",
    "read_summary",
    "run_caver",
    "select",
    "write_config",
]


@dataclass
class CaverParams:
    """CAVER settings, with ChannelBridge's defaults and why they differ.

    CAVER's own defaults are tuned for finding every pathway; ours are tuned
    for finding the one a substrate uses. Both are reasonable, they answer
    different questions, and the one you want depends on yours -- so every
    value here is a flag.
    """

    #: 1.4 A, a water probe. CAVER defaults to 0.9, which finds channels no
    #: water could pass: on our six pilot domains it gave 51 clusters against
    #: 33 at 1.4, with cluster 1 identical in span, bottleneck and exit. Use
    #: 0.9 to reproduce a run made with CAVER's defaults; use 1.4 to ask
    #: which channels a solvated substrate could actually traverse.
    probe_radius: float = 1.4
    #: How far out from the starting point CAVER looks. Large enough that
    #: the chain certainly leaves the protein: the exit is then recovered by
    #: truncating back to the surface (`Tunnel.truncate_at_surface`), which
    #: is what makes the geometry independent of this number. Measured: raw
    #: D_IC drifts up to 19 A across shell radii 3-16, the normalised value
    #: stops moving once the chain clears the surface.
    shell_radius: float = 10.0
    shell_depth: float = 6.0
    #: Average-link clustering cut. CAVER's default; 3.5 -> 4.3 merges
    #: pathways, 3.5 -> lower splits one into branches.
    clustering_threshold: float = 3.5
    #: Fixed, so a rerun reproduces the clusters. CAVER's tunnel search is
    #: stochastic; without this the cluster numbering can change under you.
    seed: int = 1
    #: Residues excluded from the protein before the search -- waters and
    #: the ligand whose site you are starting from, or the ligand blocks its
    #: own channel.
    exclude: Tuple[str, ...] = ("HOH", "WAT")
    #: Residues whose centre of gravity gives the starting point, as CAVER
    #: residue identifiers: ``"128"``, ``"B:146"``, ``"20A"`` with an
    #: insertion code. Without a chain the FIRST chain is assumed, which is
    #: silent and wrong on a multimer -- so name the chain.
    #:
    #: This is the route to use when the site is not marked by a ligand. Of
    #: the 24 enzymes in ``examples/iMARS``, ten have no ligand at all in
    #: their template entry: for those there is no centroid to start from
    #: and the catalytic residues are the only handle.
    start_residues: Tuple[str, ...] = ()
    #: Atoms whose centre of gravity contributes to the starting point, by
    #: PDB serial number. Precise, and fragile in exactly the way serial
    #: numbers are: prep renumbers atoms, so take these from the file you
    #: will actually run CAVER on.
    start_atoms: Tuple[int, ...] = ()

    def as_config(self) -> Dict[str, str]:
        out = {
            "probe_radius": "%.2f" % self.probe_radius,
            "shell_radius": "%.1f" % self.shell_radius,
            "shell_depth": "%.1f" % self.shell_depth,
            "clustering_threshold": "%.1f" % self.clustering_threshold,
            "seed": str(int(self.seed)),
        }
        if self.exclude:
            out["exclude_residue_names"] = " ".join(self.exclude)
        if self.start_residues:
            out["starting_point_residue"] = " ".join(self.start_residues)
        if self.start_atoms:
            out["starting_point_atom"] = " ".join(str(int(a)) for a in self.start_atoms)
        return out


def write_config(path: str, params: Optional[CaverParams] = None,
                 start_point: Optional[Sequence[float]] = None,
                 extra: Optional[Dict[str, str]] = None) -> str:
    """Write a CAVER configuration file.

    There are three ways to say where to start, they COMBINE rather than
    override, and CAVER weights each entity equally: a point
    (``start_point``), residues (``params.start_residues``) and atoms
    (``params.start_atoms``). At least one is required -- CAVER refuses to
    run without a starting point -- so this refuses to write a config that
    cannot run, rather than leaving the failure to CAVER.

    The starting point is written into the file and echoed in a comment,
    because CAVER moves a starting point that does not fit its probe and
    does not say by how much: the initial point is only the seed for a
    search over Voronoi vertices within ``max_distance`` (3 A by default).
    CAVER's own documentation asks for the start point to be reported with
    the output for exactly this reason.
    """
    params = params or CaverParams()
    if start_point is None and not params.start_residues and not params.start_atoms:
        raise ValueError(
            "no starting point: give start_point, or CaverParams.start_residues "
            "(e.g. ('B:128',)), or start_atoms. CAVER has no default and will "
            "not run without one, and a config written now would fail later "
            "with an error that does not mention the starting point.")
    lines = [
        "# CAVER configuration written by ChannelBridge.",
        "# Defaults are ChannelBridge's, not CAVER's: probe 1.4 A (water)",
        "# rather than CAVER's 0.9, which finds channels no water can pass.",
    ]
    if start_point is not None:
        xyz = [float(v) for v in start_point]
        if len(xyz) != 3:
            raise ValueError("start_point needs three numbers, got %r" % (start_point,))
        lines.append("# starting point as given: %.3f %.3f %.3f" % tuple(xyz))
        lines.append("starting_point_coordinates %.3f %.3f %.3f" % tuple(xyz))
    for key, value in sorted(params.as_config().items()):
        lines.append("%s %s" % (key, value))
    for key, value in sorted((extra or {}).items()):
        lines.append("%s %s" % (key, value))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


@dataclass
class Cluster:
    """One CAVER tunnel cluster, as its summary reports it."""

    cluster_id: int
    n_tunnels: int = 0
    avg_bottleneck: float = 0.0
    max_bottleneck: float = 0.0
    avg_length: float = 0.0
    priority: float = 0.0

    @property
    def label(self) -> str:
        return "cl%d" % self.cluster_id


#: Published floor: "inspect the tunnels with priority scores above 0.55",
#: which is the average priority of the third tunnel in their dataset.
PRIORITY_FLOOR = 0.55
#: "focusing on the first three tunnels" when there is no prior information.
KEEP_TOP = 3
#: "If the first tunnel is shorter than 5 A, the pocket may be on the surface
#: and tunnel analysis may not be relevant."
SHORT_TUNNEL = 5.0


def _num(parts, want, name, default=0.0):
    try:
        return float(parts[want[name]])
    except (KeyError, ValueError, IndexError):
        return default


@dataclass
class CavityPoint:
    """A derived CAVER starting point, with the evidence that it is usable."""

    xyz: Tuple[float, float, float]
    #: Distance to the nearest protein atom CENTRE. The free sphere here has
    #: radius ``free - atom_radius``.
    free: float
    #: How far this sits from the centroid of the named residues.
    offset: float
    #: Protein atoms within 10 A. The check that this is a pocket and not
    #: bulk solvent: a buried site runs 150-250, open solvent under ~50.
    burial: int

    def summary(self) -> str:
        return ("start point %.3f %.3f %.3f  (free sphere %.2f A, %.2f A from "
                "the residue centroid, burial %d atoms within 10 A)"
                % (self.xyz[0], self.xyz[1], self.xyz[2],
                   self.free, self.offset, self.burial))


def parse_residue_id(text: str) -> Tuple[Optional[str], int, str]:
    """A CAVER residue identifier -> ``(chain, number, insertion code)``.

    ``"A:157"``, ``"157"``, ``"20A"``, ``"B:20A"``. A missing chain comes
    back as ``None``, which is CAVER's "the first chain" -- silent and wrong
    on a multimer, so callers should say so.
    """
    s = str(text).strip()
    if not s:
        raise ValueError("empty residue identifier")
    chain = None
    if ":" in s:
        chain, s = s.split(":", 1)
        chain = chain.strip() or None
        s = s.strip()
    icode = ""
    if s and s[-1].isalpha():
        icode, s = s[-1], s[:-1]
    try:
        return chain, int(s), icode
    except ValueError:
        raise ValueError(
            "cannot read %r as a CAVER residue identifier; expected forms are "
            "157, A:157, 20A, B:20A" % (text,)) from None


def cavity_start_point(structure, residues: Sequence[str], *,
                       probe_radius: float = 1.4, reach: float = 8.0,
                       step: float = 0.5, atom_radius: float = 1.6) -> CavityPoint:
    """A startable point in the cavity beside the residues you name.

    Why this exists, measured. CAVER will take the residues themselves
    (``starting_point_residue``), but it then starts from their centre of
    gravity over ALL their atoms -- and for a catalytic triad, whose
    residues point at each other from three sides, that centre lands inside
    an atom. On stilbene synthase (3A5S, Cys157/His296/Asn329) it sat 0.61 A
    from the nearest atom centre and on curcuminoid synthase (3ALE) 0.49 A.
    CAVER then looks for a Voronoi vertex within ``max_distance`` (3 A) and
    gives up: both runs found **zero** tunnels, at probe 0.9 and 1.4 and at
    shell radius 3, 5 and 10. Nothing in the output said why. Moving the
    start into the cavity -- 3.7 A and 4.0 A respectively -- gave 4 clusters
    each, the best at priority 0.77 and 0.80.

    So: take the centroid of the named residues, then return the NEAREST
    point to it that could hold the probe. Nearest, not roomiest: the
    roomiest point within reach is on the way out of the pocket, and the
    question here is where the catalytic site is, not where the widest
    cavity is.

    ``free`` is measured to atom CENTRES, so a point qualifies when it is at
    least ``probe_radius + atom_radius`` away. ``atom_radius`` is one generic
    heavy-atom radius rather than per-element: this picks a starting point
    for a search that then does its own geometry properly, and a tenth of an
    angstrom here changes nothing.

    Raises if no point qualifies -- which means the named residues do not
    line an opening big enough for the probe, and is worth hearing about
    before CAVER reports it as zero tunnels.
    """
    import numpy as np

    from ..io.pdb import read_pdb

    atoms = read_pdb(structure) if isinstance(structure, str) else structure
    want = [parse_residue_id(r) for r in residues]
    if not want:
        raise ValueError("no residues given to derive a start point from")

    prot = np.array([a.xyz for a in atoms.atoms if not a.hetatm], dtype=float)
    if not len(prot):
        raise ValueError("the structure has no polymer atoms")

    keys = {(c, n, i) for c, n, i in want}
    sel = np.array(
        [a.xyz for a in atoms.atoms
         if not a.hetatm
         and ((a.chain, a.resseq, (a.icode or "").strip()) in keys
              or (None, a.resseq, (a.icode or "").strip()) in keys)],
        dtype=float)
    if not len(sel):
        raise ValueError(
            "none of the residues %s is in the structure"
            % ", ".join(str(r) for r in residues))
    centre = sel.mean(axis=0)

    need = float(probe_radius) + float(atom_radius)
    g = np.arange(-reach, reach + 1e-9, step)
    grid = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3) + centre
    offs = np.linalg.norm(grid - centre, axis=1)
    grid, offs = grid[offs <= reach], offs[offs <= reach]

    best = None
    for i in range(0, len(grid), 20000):
        chunk = grid[i:i + 20000]
        d = np.linalg.norm(chunk[:, None, :] - prot[None, :, :], axis=2)
        free = d.min(axis=1)
        ok = free >= need
        if not ok.any():
            continue
        idx = np.flatnonzero(ok)
        j = idx[int(np.argmin(offs[i:i + 20000][idx]))]
        cand = (float(offs[i:i + 20000][j]), chunk[j], float(free[j]),
                int((d[j] <= 10.0).sum()))
        if best is None or cand[0] < best[0]:
            best = cand
    if best is None:
        raise ValueError(
            "no point within %.1f A of the centroid of %s can hold a %.2f A "
            "probe (needs %.2f A clearance to an atom centre). The residues "
            "may not line an opening, or the probe may be too large for this "
            "site -- CAVER would report this as zero tunnels without saying "
            "why." % (reach, ", ".join(str(r) for r in residues),
                      probe_radius, need))
    offset, xyz, free, burial = best
    return CavityPoint(xyz=(float(xyz[0]), float(xyz[1]), float(xyz[2])),
                       free=free, offset=offset, burial=burial)


def find_caver(hint: Optional[str] = None) -> str:
    """Where CAVER's jar is. ``$CAVER`` or an explicit path.

    Not bundled, and not installable by us: CAVER 3.0 is GPL-3 software from
    Masaryk University distributed under a registration we are not entitled
    to pass on. Download it (free for academic use) and point ``$CAVER`` at
    the directory holding ``caver.jar`` and ``lib/``.
    """
    root = hint or os.environ.get("CAVER")
    if not root:
        raise FileNotFoundError(
            "CAVER was not found. Set $CAVER to the directory holding "
            "caver.jar, lib/ and bin/, or pass --caver. It is not bundled: "
            "CAVER is GPL-3 software from Masaryk University, free for "
            "academic use from caver.cz, and not ours to redistribute.")
    jar = root if root.endswith(".jar") else os.path.join(root, "caver.jar")
    if not os.path.isfile(jar):
        raise FileNotFoundError("no caver.jar at %s" % jar)
    home = os.path.dirname(os.path.abspath(jar))
    if not os.path.isfile(os.path.join(home, "bin", "atom_radii.csv")):
        raise FileNotFoundError(
            "found %s but no bin/atom_radii.csv beside it. CAVER is passed "
            "its own directory with -home and reads its atom radii from "
            "there; point $CAVER at the unpacked caver/ directory itself, "
            "not at a copy of the jar." % jar)
    return jar


def run_caver(structure: str, config: str, out_dir: str, *,
              caver: Optional[str] = None, java: str = "java",
              heap: str = "6g", timeout: Optional[float] = None,
              report=print) -> str:
    """Run CAVER on one structure. Returns the output directory.

    CAVER wants a DIRECTORY of structures, not a file, and treats each one
    as a snapshot of the same system -- which is how its priority averaging
    is defined. For a single static structure that directory holds one file,
    so priority equals throughput and the absence penalty never applies.

    The structure is COPIED into a scratch input directory rather than
    linked, because CAVER writes beside its input and a pipeline should not
    scatter output through the folder holding your prepared structures.
    """
    import shutil
    import subprocess
    import tempfile

    jar = find_caver(caver)
    if not os.path.isfile(structure):
        raise FileNotFoundError("no structure at %s" % structure)
    if not os.path.isfile(config):
        raise FileNotFoundError("no CAVER config at %s" % config)
    os.makedirs(out_dir, exist_ok=True)
    work = tempfile.mkdtemp(prefix="chbridge_caver_")
    in_dir = os.path.join(work, "in")
    os.makedirs(in_dir)
    shutil.copy2(structure, os.path.join(in_dir, os.path.basename(structure)))
    # -jar, not -cp plus a main class: the jar's manifest names the entry
    # point (caver.ui.Launcher) and resolves lib/ through its own Class-Path,
    # so naming a class here is both wrong and unnecessary. -home is NOT
    # optional -- CAVER reads bin/atom_radii.csv relative to it, and without
    # it every atom gets a default radius.
    home = os.path.dirname(os.path.abspath(jar))
    cmd = [java, "-Xmx%s" % heap, "-jar", jar, "-home", home,
           "-pdb", in_dir, "-conf", config, "-out", out_dir]
    report("  $ %s" % " ".join(cmd))
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=timeout, stdin=subprocess.DEVNULL)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "could not run %r. CAVER is pure Java; any of 8, 11 or 21 works." % java
        ) from exc
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-5:]
        raise RuntimeError("CAVER exited %d:\n  %s"
                           % (proc.returncode, "\n  ".join(tail)))
    summary = os.path.join(out_dir, "summary.txt")
    if not os.path.isfile(summary):
        raise RuntimeError(
            "CAVER exited 0 but wrote no summary.txt in %s. It does that when "
            "the starting point is outside the structure or inside an atom -- "
            "check the warnings in %s." % (out_dir, out_dir))
    return out_dir


def read_summary(path: str) -> List[Cluster]:
    """Parse CAVER's ``summary.txt`` into clusters, best priority first.

    The table is fixed-width under an ``ID No No_snaps ...`` header; parsing
    from that header rather than from a line number means a CAVER version
    that adds a column does not silently shift the priority into another
    field.
    """
    if os.path.isdir(path):
        path = os.path.join(path, "summary.txt")
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    header = None
    for i, line in enumerate(lines):
        if line.split()[:3] == ["ID", "No", "No_snaps"]:
            header = line.split()
            start = i + 1
            break
    if header is None:
        raise ValueError(
            "%s has no cluster table: expected a header beginning "
            "'ID No No_snaps'. Is this a CAVER summary.txt?" % path)
    want = {name: header.index(name) for name in
            ("Avg_BR", "Max_BR", "Avg_L", "Priority") if name in header}
    missing = {"Avg_BR", "Avg_L", "Priority"} - set(want)
    if missing:
        raise ValueError("%s: summary table has no %s column(s)"
                         % (path, ", ".join(sorted(missing))))
    out = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < len(header) or not parts[0].isdigit():
            break
        out.append(Cluster(cluster_id=int(parts[0]), n_tunnels=int(parts[1]),
                           avg_bottleneck=_num(parts, want, "Avg_BR"),
                           max_bottleneck=_num(parts, want, "Max_BR"),
                           avg_length=_num(parts, want, "Avg_L"),
                           priority=_num(parts, want, "Priority")))
    out.sort(key=lambda c: -c.priority)
    return out


@dataclass
class Selection:
    """Which tunnels to carry forward, and what the choice rests on."""

    kept: List[Cluster] = field(default_factory=list)
    dropped: List[Cluster] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def best(self) -> Optional[Cluster]:
        return self.kept[0] if self.kept else None

    def summary(self) -> str:
        lines = []
        if not self.kept:
            lines.append("no tunnel passed the filter")
        for i, c in enumerate(self.kept, 1):
            lines.append("  %-5s priority %.3f  bottleneck %.2f A  length %.1f A%s"
                         % (c.label, c.priority, c.avg_bottleneck, c.avg_length,
                            "   <- highest priority" if i == 1 else ""))
        if self.dropped:
            lines.append("  dropped %d below priority %.2f or outside the top %d"
                         % (len(self.dropped), PRIORITY_FLOOR, KEEP_TOP))
        lines.extend("  note: " + n for n in self.notes)
        return "\n".join(lines)


def select(clusters: Sequence[Cluster], *, keep: int = KEEP_TOP,
           floor: float = PRIORITY_FLOOR) -> Selection:
    """Apply the published rule: the top *keep* above the priority *floor*.

    Deliberately NOT "take cluster 1". On our pilot domains cluster 2's
    priority was 91.6-99.9 % of cluster 1's, and the published benchmark puts
    cluster 1 energetically first in 75 % of cases -- a good bet, not a fact.
    Carrying three keeps the 93 % that the first two cover, plus margin,
    and leaves the choice between them to CaverDock or to the user.
    """
    ranked = sorted(clusters, key=lambda c: -c.priority)
    above = [c for c in ranked if c.priority >= floor]
    sel = Selection(kept=above[:keep],
                    dropped=[c for c in ranked if c not in above[:keep]])
    if not ranked:
        sel.notes.append("CAVER found no tunnel clusters at all")
        return sel
    if not above:
        sel.notes.append(
            "every cluster is below priority %.2f; the highest is %.3f. The "
            "starting point may not be in a channel, or the structure may "
            "have no route out at this probe radius" % (floor, ranked[0].priority))
    if len(above) > keep:
        sel.notes.append(
            "%d clusters cleared the floor; keeping %d. Published: the top "
            "tunnel is energetically best in 75 %% of cases and the second in "
            "a further 19 %%" % (len(above), keep))
    if sel.kept and sel.kept[0].avg_length < SHORT_TUNNEL:
        sel.notes.append(
            "the highest-priority tunnel is only %.1f A long; the pocket may "
            "be on the surface, and tunnel analysis may not be meaningful "
            "here (Vavra 2024)" % sel.kept[0].avg_length)
    if len(sel.kept) >= 2:
        a, b = sel.kept[0].priority, sel.kept[1].priority
        if b >= 0.9 * a:
            sel.notes.append(
                "the top two are within %.1f %% of each other (%.3f vs %.3f): "
                "the ranking is a near-tie that a different seed could swap, "
                "so do not treat the first as determined"
                % (100 * (1 - b / a) if a else 0.0, a, b))
    return sel


def caverdock_note(selection: Selection) -> str:
    """What running CaverDock on a selection would and would not settle.

    Printed rather than buried, because the honest answer is not "run it and
    find out". Our own measurement: five runs of one identical input, varying
    only the random seed, gave an activation energy of 2.00, 3.60, 3.30, 3.30
    and 4.80 kcal/mol -- sd 1.00 at the default exhaustiveness. The
    cluster-to-cluster differences we were trying to resolve were 0.30-0.70,
    i.e. inside half a standard deviation of the method's own noise.
    """
    n = len(selection.kept)
    return "\n".join([
        "CaverDock would rank these %d by transport energetics rather than" % n,
        "geometry, which is the physically meaningful criterion. Before",
        "spending on it, two measured cautions:",
        "  * at the default exhaustiveness its run-to-run spread (sd ~1.0",
        "    kcal/mol on repeated identical runs) can exceed the differences",
        "    between clusters. Replicate, or raise exhaustiveness, or the",
        "    answer is noise.",
        "  * it needs a tunnel as input and cannot find one: it ranks what",
        "    CAVER gives it.",
        "The published alternative is to accept the priority ranking: top",
        "tunnel energetically best in 75 % of cases, top two in 93 %.",
    ])
