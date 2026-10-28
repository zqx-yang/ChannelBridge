"""Stage 3: residue classification and hotspot patch selection.

Where a binder may sit on each enzyme, and which patches of that surface are
worth carrying into the pose search.

The split in this module is deliberate. Measuring a residue -- its solvent
accessibility, its secondary structure, where its CB sits -- needs a
structure package, and in this project that means shelling out to PyMOL.
Deciding what those measurements MEAN is arithmetic, and arithmetic belongs
somewhere it can be tested without a renderer. So:

* :class:`Residue` holds the measurements, whoever made them;
* :func:`classify` turns them into one of the seven classes;
* :func:`find_patches` turns the ALLOWED ones into ranked patches;
* ``render/run_patches.py`` is the thin PyMOL shim that fills in the
  measurements and calls these.

Everything here is pure. The tests drive it with hand-built residues.

The classes are a PARTITION: every residue gets exactly one, and the order
of the tests in :func:`classify` is the priority, not an accident --

``LIGAND_CORE``
    holds the cofactor; touching it breaks the enzyme
``CHANNEL``
    lines the tunnel, and is exposed
``SHADOW``
    sits in the cylinder in front of the mouth, where a binder would cap the
    channel it is supposed to leave open
``AVOID``
    named by the user, or near enough to a named exposed residue that a
    binder centred there would cover it
``IFACE``
    buried by another chain of the assembly
``BURIED``
    not exposed enough for anything to bind
``ALLOWED``
    everything else: the usable surface

A hotspot patch is a cluster of ALLOWED residues around an anchor seed, and
the score rewards anchors, secondary structure, hydrophobicity, reach from
the mouth and exposure, while penalising side-chain entropy.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

__all__ = [
    "ANCHOR",
    "CLASS_ORDER",
    "HYDROPHOBIC",
    "N_CHI",
    "Patch",
    "PatchParams",
    "Residue",
    "classify",
    "find_patches",
    "geometry_score",
    "grow_avoid",
    "patch_score",
    "pick_hotspots",
    "write_patch_table",
    "write_residue_table",
]

#: Rotatable side-chain dihedrals. Used as an entropy penalty: a patch of
#: long floppy side chains is a worse epitope than one of short rigid ones.
#: MSE (selenomethionine) counts as MET.
N_CHI: Dict[str, int] = {
    "GLY": 0, "ALA": 0, "PRO": 1, "SER": 1, "CYS": 1, "THR": 1, "VAL": 1,
    "ILE": 2, "LEU": 2, "ASP": 2, "ASN": 2, "PHE": 2, "TYR": 2,
    "TRP": 2, "HIS": 2, "MET": 3, "GLU": 3, "GLN": 3, "LYS": 4,
    "ARG": 4, "MSE": 3,
}
HYDROPHOBIC = frozenset("ALA VAL LEU ILE MET PHE TRP TYR CYS PRO MSE".split())
#: Residues a designed binder can get a real grip on. A patch without one of
#: these is not a hotspot however well it scores on geometry.
ANCHOR = frozenset("PHE TYR TRP ILE LEU MET MSE".split())

#: Priority order. `classify` tests them in exactly this order.
CLASS_ORDER = ("LIGAND_CORE", "CHANNEL", "SHADOW", "AVOID",
               "IFACE", "BURIED", "ALLOWED")


@dataclass
class Residue:
    """One residue, as measured. Distances in angstrom, areas in A^2."""

    resi: str
    resn: str
    #: DSSP-style code; anything outside "HS" counts as loop.
    ss: str = "L"
    #: Solvent accessible area in the ASSEMBLY -- what a binder actually sees.
    sasa: float = 0.0
    #: How much area the rest of the assembly buries: monomer minus assembly.
    dsasa: float = 0.0
    #: CB, or CA for glycine.
    cb: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    #: Distance from the channel mouth to this residue.
    dist_mouth: float = 0.0
    #: Angle at the mouth between its outward axis and this residue.
    angle_mouth: float = 0.0
    #: Signed distance along the channel axis; positive is outside the mouth.
    axial: float = 0.0
    #: Distance from the channel axis.
    perp: float = 0.0
    #: ``(axial, perp)`` at EVERY opening, primary first. An enzyme may have
    #: more than one mouth, and a binder capping a secondary one blocks the
    #: channel just as effectively -- so SHADOW is tested against all of
    #: them, while the table reports the primary. Residues 100 and 103 of
    #: the worked example's GDH are SHADOW on exactly this basis: 22 A off the primary axis,
    #: inside a secondary mouth's cylinder.
    openings: List[Tuple[float, float]] = field(default_factory=list)

    @property
    def n_chi(self) -> int:
        return N_CHI.get(self.resn, 3)

    @property
    def is_anchor(self) -> bool:
        return self.resn in ANCHOR

    @property
    def is_hydrophobic(self) -> bool:
        return self.resn in HYDROPHOBIC

    @property
    def structured(self) -> bool:
        return self.ss in ("H", "S")


@dataclass
class PatchParams:
    """Every threshold the stage uses, with the values the worked example was run at.

    These are the defaults of the original script, kept exactly: changing one
    silently would make old and new patch tables incomparable.
    """

    #: Within this of the cofactor is LIGAND_CORE.
    ligand_r: float = 5.0
    #: Within this of the channel is CHANNEL, when no CAVER residues are given.
    channel_r: float = 12.0
    #: Radius of the cylinder in front of the mouth a binder must not block.
    shadow_r: float = 13.0
    #: How far back INSIDE the mouth that cylinder still reaches.
    shadow_back: float = 10.0
    #: SASA at or above which a residue counts as exposed.
    exposed: float = 15.0
    #: dSASA above which a residue is buried by another chain.
    iface_cut: float = 15.0
    #: A seed anchor needs at least this much exposed area.
    seed_sasa: float = 25.0
    #: Patch members are ALLOWED residues within this of the seed's CB.
    patch_r: float = 10.0
    #: Fewer members than this is not a patch.
    patch_min: int = 3
    #: Two kept patches must be at least this far apart, centroid to centroid.
    patch_sep: float = 9.0
    #: How many patches to keep per enzyme.
    n_patches: int = 6
    #: Band of distances from the channel axis that scores best.
    axis_lo: float = 16.0
    axis_hi: float = 32.0
    #: Band of distances from the mouth that scores best.
    reach_lo: float = 15.0
    reach_hi: float = 35.0
    #: A named AVOID residue spreads its exclusion this far...
    avoid_r: float = 12.0
    #: ...but only if it is at least this exposed. A buried catalytic residue
    #: protects only itself: nothing on the surface can reach it, and growing
    #: a 12 A ball from it would remove a third of a usable face for nothing.
    avoid_min_sasa: float = 15.0


def _dist(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(sum((float(x) - float(y)) ** 2 for x, y in zip(a, b)))


def grow_avoid(named: Iterable[str], residues: Dict[str, Residue],
               params: PatchParams,
               radius: Optional[Dict[str, float]] = None) -> set:
    """Residues excluded because the user named them, or they sit beside one.

    A hotspot is the CENTRE of a binder footprint, not its edge, so keeping a
    binder off a named residue means keeping hotspots further away than the
    residue itself. The exclusion is therefore grown.

    It is grown only from EXPOSED named residues. That asymmetry matters: a
    catalytic residue with 5 A^2 of exposed area sits under the surface where
    no epitope can reach it, and growing a ball from it would delete usable
    surface to protect something already unreachable.
    """
    radius = radius or {}
    out: set = set()
    for r in named:
        res = residues.get(r)
        if res is None:
            continue
        out.add(r)
        if res.sasa < params.avoid_min_sasa:
            continue                      # buried: protects only itself
        grow = float(radius.get(r, params.avoid_r))
        for other, o in residues.items():
            if _dist(o.cb, res.cb) <= grow:
                out.add(other)
    return out


def classify(res: Residue, *, in_core: bool, in_channel: bool,
             in_avoid: bool, params: PatchParams) -> str:
    """Which of the seven classes this residue is. Exactly one."""
    exposed = res.sasa >= params.exposed
    cylinders = res.openings or [(res.axial, res.perp)]
    shadow = any(ax > -params.shadow_back and pp < params.shadow_r
                 for ax, pp in cylinders)
    if in_core:
        return "LIGAND_CORE"
    if in_channel and exposed:
        return "CHANNEL"
    if shadow and exposed:
        return "SHADOW"
    if in_avoid:
        return "AVOID"
    if res.dsasa > params.iface_cut:
        return "IFACE"
    if not exposed:
        return "BURIED"
    return "ALLOWED"


def geometry_score(dist_mouth: float, perp: float, params: PatchParams) -> float:
    """How well placed a patch is, 0 to 1.

    Zero inside the shadow cylinder -- a binder there caps the channel, and
    no amount of good chemistry redeems that. Otherwise the mean of two
    plateau functions: one on the distance from the channel axis, one on the
    distance from the mouth. Both are flat across a band and fall off
    linearly outside it, because "far enough to clear the mouth, near enough
    to bridge" has a range rather than an optimum.
    """
    if perp < params.shadow_r:
        return 0.0
    lo, hi = params.axis_lo, params.axis_hi
    s_p = 1.0 if lo <= perp <= hi else max(
        0.0, 1.0 - abs(perp - (lo if perp < lo else hi)) / 12.0)
    dlo, dhi = params.reach_lo, params.reach_hi
    s_d = 1.0 if dlo <= dist_mouth <= dhi else max(
        0.0, 1.0 - abs(dist_mouth - (dlo if dist_mouth < dlo else dhi)) / 15.0)
    return 0.5 * s_p + 0.5 * s_d


@dataclass
class Patch:
    """One candidate epitope."""

    seed: str
    members: List[str]
    score: float
    geom_score: float
    centroid: Tuple[float, float, float]
    n_res: int
    anchors: int
    ss_frac: float
    phob_frac: float
    mean_n_chi: float
    mean_sasa: float
    dist_mouth: float
    angle_mouth: float
    perp: float
    rank: int = 0
    sequence: List[str] = field(default_factory=list)

    @property
    def label(self) -> str:
        return "patch_%02d" % self.rank


def patch_score(*, anchors: int, ss_frac: float, phob_frac: float,
                geom: float, mean_sasa: float, mean_n_chi: float) -> float:
    """The weighted sum, exactly as the original script computed it.

    Anchors saturate at three: a fourth buys nothing, because a binder
    cannot grip more than a few. Side-chain entropy is the only negative
    term.
    """
    return (1.2 * min(anchors, 3) / 3.0
            + 1.5 * ss_frac
            + 1.0 * phob_frac
            + 1.5 * geom
            + 0.6 * min(mean_sasa / 80.0, 1.0)
            - 0.30 * mean_n_chi)


def find_patches(residues: Dict[str, Residue], classes: Dict[str, str],
                 params: Optional[PatchParams] = None) -> List[Patch]:
    """Rank the hotspot patches of one enzyme, best first.

    Seeds are exposed anchor residues in a helix or strand; a patch is the
    ALLOWED residues within ``patch_r`` of a seed. Patches are then taken in
    score order, skipping any whose centroid is within ``patch_sep`` of one
    already kept, so the result is spread over the surface instead of being
    six views of the same spot.
    """
    params = params or PatchParams()
    allowed = {r: v for r, v in residues.items() if classes.get(r) == "ALLOWED"}
    seeds = [r for r, v in allowed.items()
             if v.is_anchor and v.structured and v.sasa >= params.seed_sasa]

    cand: List[Patch] = []
    for s in seeds:
        members = [r for r, v in allowed.items()
                   if _dist(v.cb, allowed[s].cb) <= params.patch_r]
        if len(members) < params.patch_min:
            continue
        if not any(allowed[m].is_anchor for m in members):
            continue
        n = float(len(members))
        ss_frac = sum(1 for m in members if allowed[m].structured) / n
        phob = sum(1 for m in members if allowed[m].is_hydrophobic) / n
        chi = sum(allowed[m].n_chi for m in members) / n
        anchors = sum(1 for m in members if allowed[m].is_anchor)
        cen = tuple(sum(allowed[m].cb[i] for m in members) / n for i in range(3))
        d = sum(allowed[m].dist_mouth for m in members) / n
        th = sum(allowed[m].angle_mouth for m in members) / n
        perp = sum(allowed[m].perp for m in members) / n
        sasa = sum(allowed[m].sasa for m in members) / n
        g = geometry_score(d, perp, params)
        cand.append(Patch(
            seed=s, members=sorted(members, key=_as_int), score=patch_score(
                anchors=anchors, ss_frac=ss_frac, phob_frac=phob, geom=g,
                mean_sasa=sasa, mean_n_chi=chi),
            geom_score=g, centroid=cen, n_res=int(n), anchors=anchors,
            ss_frac=ss_frac, phob_frac=phob, mean_n_chi=chi, mean_sasa=sasa,
            dist_mouth=d, angle_mouth=th, perp=perp,
            sequence=["%s%s" % (allowed[m].resn, m)
                      for m in sorted(members, key=_as_int)]))

    cand.sort(key=lambda p: -p.score)
    kept: List[Patch] = []
    for p in cand:
        if all(_dist(p.centroid, k.centroid) > params.patch_sep for k in kept):
            p.rank = len(kept) + 1
            kept.append(p)
        if len(kept) >= params.n_patches:
            break
    return kept


def _as_int(resi: str) -> int:
    try:
        return int("".join(c for c in str(resi) if c.isdigit() or c == "-"))
    except ValueError:
        return 0


_PATCH_COLUMNS = ("rank", "score", "geom_score", "seed", "n_res", "anchors",
                  "ss_pct", "hydrophobic_pct", "mean_n_chi", "mean_sasa",
                  "dist_mouth", "angle_mouth", "perp_offset",
                  "hotspot_residues", "sequence")


#: How many residues of a patch become hotspots. A patch is 7-14 residues;
#: a generator wants the handful the binder should actually lock onto, not
#: the whole epitope. Five is the worked example's value and BindCraft's
#: advice ("several surface residues within a radial patch").
N_HOTSPOTS = 5

#: Below this solvent-accessible area a residue cannot be contacted by a
#: binder, so it is not a hotspot however good its type.
HOTSPOT_MIN_SASA = 15.0


def pick_hotspots(members: Sequence[str], info: Mapping[str, Tuple[str, float]],
                  *, n: int = N_HOTSPOTS,
                  min_sasa: float = HOTSPOT_MIN_SASA,
                  anchors_only: bool = True) -> List[Tuple[str, str]]:
    """``n`` residues of one patch to hand a generator, in sequence order.

    *members* are the patch's residue ids; *info* maps a residue id to
    ``(resname, sasa)``. Returns ``[(resi, resname)]``, which may be SHORTER
    than *n* -- see below.

    Two rules, and they pull against each other.

    **Anchors only.** A hotspot should be one of the big, few-rotamer
    hydrophobics a binder actually locks onto (:data:`ANCHOR`: F, Y, W, I, L,
    M), which is also BindCraft's advice for choosing a site. A patch of 7-14
    residues contains plenty that are not: the whole patch handed over as
    hotspots includes glycines, alanines and prolines, and those are the
    three residue types RFdiffusion3's ``TIP`` keyword is not defined for --
    so an unfiltered patch forces a per-residue atom table that exists only
    to describe residues that should not have been hotspots.

    **Spread.** Hotspots clustered on one turn of a helix tell a generator
    far less about where the binder should sit than ones that span the patch,
    so the patch is cut into ``n`` runs along the sequence and the best
    candidate of each run is taken. Picking the top ``n`` by score alone
    collapses onto whichever element happens to be most exposed.

    A run with no anchor contributes nothing rather than its best
    non-anchor, so a patch with three usable anchors yields three hotspots.
    That is the honest answer: the alternative is to pad the list with
    residues the binder cannot grip, which reads as five anchors to anyone
    downstream. ``anchors_only=False`` restores the older behaviour of
    taking the best of every run whatever its type.

    Anything under *min_sasa* is too buried to be contacted and is never
    returned, anchor or not.
    """
    have = [m for m in members if m in info
            and float(info[m][1]) >= float(min_sasa)
            and (not anchors_only or (info[m][0] or "").upper() in ANCHOR)]
    if not have:
        return []
    have.sort(key=_as_int)

    def score(m: str) -> float:
        resn, sasa = info[m]
        s = min(float(sasa), 120.0) / 120.0
        if (resn or "").upper() in ANCHOR:
            s += 0.6
        return s

    n = max(1, int(n))
    out: List[str] = []
    for i in range(n):
        run = have[i * len(have) // n:(i + 1) * len(have) // n]
        if run:
            best = max(run, key=score)
            if best not in out:
                out.append(best)
    out.sort(key=_as_int)
    return [(m, info[m][0]) for m in out]


def write_patch_table(path: str, patches: Sequence[Patch]) -> str:
    """The ``<NAME>_hotspot_patches.tsv`` the later stages read.

    Column names and order match the original script's output exactly, so a
    table written here is interchangeable with one written there.
    """
    lines = ["\t".join(_PATCH_COLUMNS)]
    for p in patches:
        lines.append("\t".join(str(v) for v in (
            p.rank, "%.2f" % p.score, "%.2f" % p.geom_score, p.seed, p.n_res,
            p.anchors, "%.0f" % (100 * p.ss_frac), "%.0f" % (100 * p.phob_frac),
            "%.2f" % p.mean_n_chi, "%.0f" % p.mean_sasa,
            "%.1f" % p.dist_mouth, "%.1f" % p.angle_mouth, "%.1f" % p.perp,
            "+".join(p.members), " ".join(p.sequence))))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


_RESIDUE_COLUMNS = ("chain", "resi", "resn", "ss", "class", "sasa",
                    "dsasa_assembly", "dist_mouth", "angle_mouth", "axial",
                    "perp_offset", "n_chi")


def write_residue_table(path: str, chain: str, residues: Dict[str, Residue],
                        classes: Dict[str, str]) -> str:
    """The ``<NAME>_residue_classes.tsv`` the map session paints from."""
    lines = ["\t".join(_RESIDUE_COLUMNS)]
    for resi in sorted(residues, key=_as_int):
        r = residues[resi]
        lines.append("\t".join(str(v) for v in (
            chain, resi, r.resn, r.ss or "L", classes.get(resi, "ALLOWED"),
            "%.1f" % r.sasa, "%.1f" % r.dsasa, "%.1f" % r.dist_mouth,
            "%.1f" % r.angle_mouth, "%.1f" % r.axial, "%.1f" % r.perp,
            r.n_chi)))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# --------------------------------------------------- a partner's footprint
#
# A binding partner that is NOT part of the design -- ferredoxin on FNR -- is
# still a reason not to put an epitope somewhere. It enters as its own
# complex, superposed onto the target, rather than as an extra chain in the
# target file: adding the partner to the structure would also make it a
# steric obstacle in every pose, which is a different and much stronger
# claim than "do not put a binder here".

@dataclass(frozen=True)
class Footprint:
    """Where a partner protein would sit on the target."""

    #: Residue ids of the target the partner covers.
    residues: Tuple[str, ...]
    #: How well the partner complex's own copy of the target superposed.
    rmsd: float
    #: How many residues the superposition was solved on.
    n_matched: int
    contact: float

    @property
    def selection(self) -> str:
        """As a PyMOL selection, for the classification stage to consume."""
        return "resi " + "+".join(self.residues) if self.residues else ""

    def warnings(self) -> List[str]:
        out = []
        if self.rmsd > 3.0:
            out.append(
                "the partner complex superposed at %.2f A RMSD over %d residues, "
                "which is high for a superposition whose job is to place a "
                "partner: check the align chain is the right one before "
                "trusting this footprint" % (self.rmsd, self.n_matched))
        if not self.residues:
            out.append(
                "no target residue lies within %.1f A of the partner, so the "
                "footprint is empty -- either the chains are wrong or the two "
                "do not touch" % self.contact)
        return out


def _ca_by_residue(atoms, chain: Optional[str] = None) -> Dict[str, Tuple[float, float, float]]:
    out = {}
    for a in atoms:
        if a.hetatm or a.name.strip() != "CA":
            continue
        if chain is not None and a.chain != chain:
            continue
        out["%d%s" % (a.resseq, a.icode.strip())] = a.xyz
    return out


def partner_footprint(target_atoms, target_chain: Optional[str], partner_atoms,
                      align_chain: str, partner_chain: str, *,
                      contact: float = 5.0) -> Footprint:
    """Target residues a partner protein would cover.

    *partner_atoms* holds a complex of the target (``align_chain``) with its
    partner (``partner_chain``). That copy of the target is superposed onto
    ours, the same transform is applied to the partner, and every target
    residue within *contact* angstroms of it is returned.

    The correspondence is by residue number, which is what makes this cheap
    and is the normal case for a homolog of the same protein; residues the
    two do not share are simply not used, and ``n_matched`` says how many
    were. A high ``rmsd`` means the chains were probably mispaired, so it is
    reported rather than swallowed -- see :meth:`Footprint.warnings`.
    """
    import numpy as np

    from .geometry import superpose

    ours = _ca_by_residue(target_atoms, target_chain)
    theirs = _ca_by_residue(partner_atoms, align_chain)
    shared = sorted(set(ours) & set(theirs), key=_as_int)
    if len(shared) < 3:
        raise ValueError(
            "the partner complex's chain %s shares %d residue number(s) with "
            "the target%s, too few to superpose on. The two are probably not "
            "the same protein, or are numbered differently."
            % (align_chain, len(shared),
               " chain %s" % target_chain if target_chain else ""))
    rot, trans, rmsd = superpose(
        np.array([theirs[r] for r in shared]),
        np.array([ours[r] for r in shared]))

    partner_xyz = np.array([a.xyz for a in partner_atoms
                            if a.chain == partner_chain and not a.hetatm], float)
    if not len(partner_xyz):
        raise ValueError(
            "the partner complex has no polymer atoms in chain %s. Chains "
            "present: %s" % (partner_chain, ", ".join(sorted(
                {a.chain for a in partner_atoms if a.chain.strip()}))))
    moved = partner_xyz @ rot.T + trans

    hits = set()
    cut2 = float(contact) ** 2
    for a in target_atoms:
        if a.hetatm or (target_chain is not None and a.chain != target_chain):
            continue
        d2 = ((moved - np.asarray(a.xyz, float)) ** 2).sum(1)
        if d2.min() <= cut2:
            hits.add("%d%s" % (a.resseq, a.icode.strip()))
    return Footprint(residues=tuple(sorted(hits, key=_as_int)), rmsd=rmsd,
                     n_matched=len(shared), contact=float(contact))
