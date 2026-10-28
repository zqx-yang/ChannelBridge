"""Stage 5: the pose search.

The method, and the only stage that cannot be handed an answer from outside.
Everything before it describes two enzymes; this searches for a relative
placement of them that a bispecific binder could actually hold, and scores
each candidate in iMARS terms.

What a pose is
--------------

A pose is a rigid-body transform of the **downstream** enzyme; the upstream
one never moves. It is parameterised the way the metric is, not the way a
docking program would:

``d_ic``
    separation of the two channel exits, held exactly.
``o_ic``
    the iMARS angle P1-P2-P3 at the **upstream** exit. 180 degrees puts the
    downstream exit straight out along the upstream channel; smaller values
    swing it off that axis by ``180 - o_ic``.
``azimuth``
    where on that cone the downstream enzyme sits. Degenerate at
    ``o_ic = 180``, where the cone is a point, and dropped there.
``tilt``
    how far the downstream channel turns away from pointing back at the
    upstream exit. ``tilt = 0`` means aimed straight back, which is the one
    orientation ``O_IC_rev`` calls 180 degrees.
``tilt_azimuth``
    which way it tilts. See `The orientation the default grid does not
    sample`_ -- this axis exists because leaving it out is a real gap, and
    defaults to unswept because sweeping it costs twelve times the grid.
``spin``
    roll about the downstream channel. Changes nothing about D_IC or O_IC and
    everything about which face of the downstream enzyme is presented, so it
    is what decides its epitope.

Why the angles are built on the O_IC arm
----------------------------------------

The direction a pose is built along is ``Tunnel.axis``: exit minus the point
one O_IC arm back, the same two points :func:`channelbridge.geometry.o_ic`
measures the result with. :meth:`~channelbridge.tunnel.Tunnel.fitted_axis` is
steadier and is 9-12 degrees away from it on real tunnels, which would leave
the grid asking for one O_IC and recording another. Building on the scored
direction instead makes "requested" and "realised" the same number, and
:func:`search` asserts it.

The same axis is used for the downstream enzyme, which is why ``tilt = 0``
gives ``O_IC_rev = 180`` exactly rather than approximately.

Three filters, cheapest first
-----------------------------

A pose has to survive all of them to be feasible:

1. **clash** -- the two enzymes may not interpenetrate.
2. **binder gap** -- they must be at least ``binder_gap`` apart, because
   something has to run between them.
3. **a bridgeable epitope pair** -- some patch on one enzyme and some patch
   on the other must be the right distance apart *and* have a clear straight
   line between them. A geometry with no such pair scores well and cannot be
   built, which is the whole reason this stage exists.
4. **occupancy**, where the upstream enzyme is a homomultimer: the pose has
   to work on every protomer at once. Run last, on poses that already passed
   the rest, because it is by far the most expensive.

The grid is ordered best-DO-first, so the first feasible pose found is the
best one, and ``max_ok`` can stop the scan without biasing it.

The orientation the default grid does not sample
------------------------------------------------

O_IC constrains where the downstream exit sits and says nothing whatever
about the downstream enzyme's own orientation, which is three free
parameters. ``tilt`` and ``spin`` are two of them; the third is which way
the channel tilts, ``tilt_azimuth``. The default is a single value, so the
default grid reaches the downstream channel directions on one arc rather
than over the sphere.

That default is for grid size, not because the axis does not matter, and
``examples/iMARS`` measures what it costs. Sweep it when a pair comes back
with no feasible pose: it is the cheapest dimension to add, because nothing
about D_IC, O_IC or DO changes along it -- only the clash and epitope tests.

Pure numpy, no PyMOL
--------------------

The reference implementation ran inside PyMOL because the stages around it
need SASA and a session. This one reads what those stages wrote -- prepared
PDBs, normalised tunnels, the patch tables -- so it is importable, testable
and fast. Two consequences worth knowing:

* protomer copies are superposed here (Kabsch over matching CA atoms)
  instead of by ``cmd.super``, and from the structure's own protomers rather
  than an assumed symmetry operator, so an imperfect assembly still gets the
  copies it actually has.
* distances from a point cloud to a line segment are solved analytically
  rather than by sampling the segment every angstrom. That is both faster
  and strictly more correct: a 1 A sample can step over a thin obstruction.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import DO_PUBLISHED, DOCoefficients, do_score, o_ic, superpose
from .reference_pose import min_interprotein_distance, rotation_about

__all__ = [
    "STATUSES",
    "PatchPair",
    "PatchSite",
    "Pose",
    "PoseParams",
    "PoseSide",
    "SearchResult",
    "equivalent_chains",
    "path_is_clear",
    "pose_grid",
    "pose_transform",
    "read_patch_sites",
    "search",
    "segment_cloud_distance",
    "segment_segment_distance",
    "symmetry_transforms",
    "write_complex",
    "write_patch_pair_best",
    "write_patch_pairs",
    "write_poses",
]

#: Every verdict a pose can get. ``ok`` is the only feasible one.
STATUSES = ("ok", "clash", "no_binder_gap", "no_bridgeable_pair", "occupancy_clash")

#: Below this D_IC the two exits coincide and the O_IC vertex is degenerate,
#: so the angle the pose was built to is not a measurable property of it.
MIN_DIC = 0.5


# --------------------------------------------------------------- geometry

def _unit(v, what: str = "vector") -> np.ndarray:
    v = np.asarray(v, float).reshape(3)
    n = float(np.linalg.norm(v))
    if n < 1e-9:
        raise ValueError("%s has zero length" % what)
    return v / n


def _frame(u, helper) -> np.ndarray:
    """Right-handed orthonormal frame with ``e1 = u``, as rows."""
    e1 = _unit(u, "frame axis")
    h = _unit(helper, "frame helper")
    if abs(float(e1 @ h)) > 0.95:
        h = np.array([1.0, 0.0, 0.0]) if abs(e1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e2 = h - e1 * float(h @ e1)
    e2 = e2 / float(np.linalg.norm(e2))
    return np.stack([e1, e2, np.cross(e1, e2)])


def _perp(v, axis) -> np.ndarray:
    u = _unit(axis, "axis")
    v = np.asarray(v, float).reshape(3)
    return v - u * float(v @ u)


def _dihedral(v1, axis, v2) -> float:
    """Signed angle between the parts of *v1* and *v2* across *axis*."""
    a = _unit(axis, "dihedral axis")
    n1, n2 = _perp(v1, a), _perp(v2, a)
    if float(np.linalg.norm(n1)) < 1e-6 or float(np.linalg.norm(n2)) < 1e-6:
        return float("nan")
    return math.degrees(math.atan2(float(np.cross(n1, n2) @ a), float(n1 @ n2)))


def segment_cloud_distance(p, q, xyz, *, t_lo: float = 0.0,
                           t_hi: float = 1.0) -> float:
    """Closest approach between the segment ``p->q`` and a point cloud.

    Solved per point by projecting onto the segment and clamping, so it is
    exact -- no sampling step to step over. ``t_lo``/``t_hi`` restrict the
    parameter range, which is how the end margins of :func:`path_is_clear`
    are applied.
    """
    p = np.asarray(p, float).reshape(3)
    q = np.asarray(q, float).reshape(3)
    xyz = np.atleast_2d(np.asarray(xyz, float))
    if xyz.size == 0:
        return float("inf")
    v = q - p
    L2 = float(v @ v)
    if L2 < 1e-12:
        return float(np.linalg.norm(xyz - p, axis=1).min())
    t = ((xyz - p) @ v) / L2
    t = np.clip(t, float(t_lo), float(t_hi))
    closest = p[None, :] + t[:, None] * v[None, :]
    return float(np.linalg.norm(xyz - closest, axis=1).min())


def path_is_clear(p, q, xyz, clearance: float = 3.5, margin: float = 7.0) -> bool:
    """Is the straight line between two surface points free of protein?

    Both endpoints sit *on* the surface, so ``margin`` angstroms at each end
    are not tested -- without that every patch reports itself blocked by its
    own residues. A segment shorter than the two margins has nothing left to
    test and counts as clear.
    """
    p = np.asarray(p, float).reshape(3)
    q = np.asarray(q, float).reshape(3)
    L = float(np.linalg.norm(q - p))
    if L < 1e-6:
        return False
    if L <= 2.0 * float(margin):
        return True
    f = float(margin) / L
    return segment_cloud_distance(p, q, xyz, t_lo=f, t_hi=1.0 - f) >= float(clearance)


def segment_segment_distance(p1, q1, p2, q2) -> float:
    """Closest approach between two line SEGMENTS (not infinite lines)."""
    p1 = np.asarray(p1, float).reshape(3)
    q1 = np.asarray(q1, float).reshape(3)
    p2 = np.asarray(p2, float).reshape(3)
    q2 = np.asarray(q2, float).reshape(3)
    d1, d2 = q1 - p1, q2 - p2
    r = p1 - p2
    a, e, f = float(d1 @ d1), float(d2 @ d2), float(d2 @ r)
    if a < 1e-9 and e < 1e-9:
        return float(np.linalg.norm(r))
    if a < 1e-9:
        sc, tc = 0.0, min(1.0, max(0.0, f / e))
    else:
        c = float(d1 @ r)
        if e < 1e-9:
            tc, sc = 0.0, min(1.0, max(0.0, -c / a))
        else:
            b = float(d1 @ d2)
            den = a * e - b * b
            sc = min(1.0, max(0.0, (b * f - c * e) / den)) if den > 1e-9 else 0.0
            tc = (b * sc + f) / e
            if tc < 0.0:
                tc, sc = 0.0, min(1.0, max(0.0, -c / a))
            elif tc > 1.0:
                tc, sc = 1.0, min(1.0, max(0.0, (b - c) / a))
    return float(np.linalg.norm((p1 + sc * d1) - (p2 + tc * d2)))


# ------------------------------------------------- the proximity question
#
# The scan asks one question tens of thousands of times: is anything in the
# downstream enzyme closer than `cut` to anything in the upstream one? Asked
# exactly that is 5000 x 3000 distances a pose, which is minutes of grid at
# tens of milliseconds each; a voxel lookup is one array index per atom.
#
# A voxel lookup is also approximate, and the usual way of using one -- a
# conservative prefilter, exact check on the survivors -- is the wrong way
# round here. Nearly every pose in a best-DO-first grid interpenetrates, so
# "conservative" would send nearly every pose to the exact check and save
# nothing. What is needed is a cheap DEFINITE reject.
#
# So each threshold gets two grids, eroded and dilated by the half-diagonal
# of a voxel, which bounds the only approximation left (the query point is
# snapped to the nearest voxel centre; the marking below is exact). A hit in
# the eroded grid means definitely closer; a miss in the dilated grid means
# definitely not; between the two the exact distance is measured. The verdict
# is therefore exact everywhere, and the exact measurement runs on the few
# poses that are genuinely near the boundary.

def _clearance_grid(xyz, clearance: float, spacing: float):
    """Voxels whose CENTRE lies within *clearance* of an atom.

    Marked exactly -- each atom's own sub-voxel offset is kept and the
    distance to each candidate voxel centre is tested -- rather than by
    growing a ball of whole voxels around a snapped atom position, which
    would add a second half-voxel of slop to every answer.
    """
    xyz = np.atleast_2d(np.asarray(xyz, float))
    if xyz.size == 0 or xyz.shape[-1] != 3:
        raise ValueError("expected (N, 3) coordinates, got %r" % (xyz.shape,))
    spacing = float(spacing)
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    clearance = max(0.0, float(clearance))
    pad = clearance + 2.0 * spacing
    origin = xyz.min(axis=0) - pad
    shape = tuple(math.ceil(e / spacing) + 1
                  for e in (xyz.max(axis=0) + pad - origin))
    grid = np.zeros(shape, dtype=bool)
    base = np.floor((xyz - origin) / spacing).astype(int)
    r = math.ceil(clearance / spacing) + 1
    offs = np.arange(-r, r + 1)
    for i in offs:
        for j in offs:
            for k in offs:
                if (i * i + j * j + k * k) * spacing * spacing > (clearance + spacing * 1.74) ** 2:
                    continue
                cand = base + np.array([i, j, k])
                centre = origin + cand * spacing
                keep = ((centre - xyz) ** 2).sum(1) <= clearance * clearance
                if not keep.any():
                    continue
                sel = cand[keep]
                ok = np.ones(len(sel), bool)
                for ax in range(3):
                    ok &= (sel[:, ax] >= 0) & (sel[:, ax] < shape[ax])
                sel = sel[ok]
                if len(sel):
                    grid[sel[:, 0], sel[:, 1], sel[:, 2]] = True
    return grid, origin, spacing


def _hits_mask(grid_tuple, xyz) -> np.ndarray:
    """Which of *xyz* snap onto a marked voxel."""
    grid, origin, spacing = grid_tuple
    idx = np.rint((np.asarray(xyz, float) - origin) / spacing).astype(int)
    out = np.ones(len(idx), bool)
    for ax in range(3):
        out &= (idx[:, ax] >= 0) & (idx[:, ax] < grid.shape[ax])
    if not out.any():
        return out
    sel = idx[out]
    out[out] = grid[sel[:, 0], sel[:, 1], sel[:, 2]]
    return out


def _hits_any(grid_tuple, xyz) -> bool:
    """Does any of *xyz* snap onto a marked voxel?"""
    return bool(_hits_mask(grid_tuple, xyz).any())


#: How wide a neighbourhood :func:`_min_distance` trusts a voxel grid to
#: bound. Comfortably above any clash or binder-gap cut-off, so the narrowed
#: computation is exact for every distance the verdict turns on.
NEAR_RADIUS = 12.0


def _min_distance(up_xyz, near_grid, moved) -> float:
    """The exact closest approach, computed on the atoms that can matter.

    The honest version of this is every upstream atom against every
    downstream one: 4,000 x 3,000 distances, about 0.3 s, and the scan needs
    it thousands of times. Almost none of those pairs can be the closest,
    though, so the downstream atoms are first narrowed to those within
    :data:`NEAR_RADIUS` of the upstream enzyme (one voxel lookup each) and
    the upstream atoms to a box around what is left.

    Still exact wherever it matters: the closest pair survives both filters
    whenever the true distance is under the radius. When nothing at all is
    within the radius the narrowing has nothing to offer, and the full
    computation runs rather than a bound being reported as a measurement.
    """
    up_xyz = np.asarray(up_xyz, float)
    moved = np.asarray(moved, float)
    mask = _hits_mask(near_grid, moved)
    if not mask.any():
        return min_interprotein_distance(up_xyz, moved)
    sub = moved[mask]
    lo = sub.min(axis=0) - NEAR_RADIUS
    hi = sub.max(axis=0) + NEAR_RADIUS
    keep = np.all((up_xyz >= lo) & (up_xyz <= hi), axis=1)
    if not keep.any():
        return min_interprotein_distance(up_xyz, moved)
    return min_interprotein_distance(up_xyz[keep], sub)


class _Proximity:
    """Is anything closer than *cut*? Yes, no, or 'measure it exactly'."""

    def __init__(self, xyz, cut: float, spacing: float):
        half = 0.5 * spacing * math.sqrt(3.0)
        self.cut = float(cut)
        # eroded: a hit here is definitely closer than cut
        self._lo = (_clearance_grid(xyz, self.cut - half, spacing)
                    if self.cut > half else None)
        # dilated: a miss here is definitely not closer than cut
        self._hi = _clearance_grid(xyz, self.cut + half, spacing)

    def test(self, xyz) -> Optional[bool]:
        if self._lo is not None and _hits_any(self._lo, xyz):
            return True
        if not _hits_any(self._hi, xyz):
            return False
        return None


# -------------------------------------------------------- what goes in

@dataclass(frozen=True)
class PatchSite:
    """One candidate epitope, as stage 5 needs it.

    ``centroid`` is the mean CB (CA for glycine) of the member residues --
    recomputed from the structure rather than read from a column, so a patch
    table written by any version of stage 3 can be used.
    """

    rank: int
    centroid: np.ndarray
    residues: Tuple[str, ...] = ()
    score: float = 0.0

    @property
    def label(self) -> str:
        return "p%d" % self.rank


@dataclass
class PoseSide:
    """One enzyme, ready to be posed."""

    name: str
    #: Protein heavy atoms, the clash and path tests' cloud.
    xyz: np.ndarray
    #: P2, the channel exit.
    exit: np.ndarray
    #: P1, one O_IC arm inside the channel.
    axis_point: np.ndarray
    patches: Sequence[PatchSite] = ()
    #: Whole-enzyme centroid, for the ``phi`` dihedral only.
    centroid: Optional[np.ndarray] = None
    #: Carried through so :func:`write_complex` can write a real PDB.
    atoms: object = None
    #: ``(chain, R, t)`` per protomer of a homomultimer, identity first.
    symmetry: Sequence[Tuple[str, np.ndarray, np.ndarray]] = ()

    def __post_init__(self):
        self.xyz = np.atleast_2d(np.asarray(self.xyz, float))
        self.exit = np.asarray(self.exit, float).reshape(3)
        self.axis_point = np.asarray(self.axis_point, float).reshape(3)
        if self.centroid is None:
            self.centroid = self.xyz.mean(axis=0)
        else:
            self.centroid = np.asarray(self.centroid, float).reshape(3)

    @property
    def axis(self) -> np.ndarray:
        """Outward channel direction over the O_IC arm."""
        return _unit(self.exit - self.axis_point, "%s channel axis" % self.name)

    @property
    def tag_name(self) -> str:
        """The name as it appears inside a machine-read tag.

        Letters and digits only. The epitope-pair tag is
        ``<up>p<N>-<dn>p<M>``, which the figure stage parses out of a file
        name, and a separator inside an enzyme name splits the field: the
        iMARS set calls its enzymes ``enz_4CL``, and with the underscore
        kept both patch ids were silently lost from every file name.
        ``best_pair`` keeps the real names for people to read.
        """
        return re.sub(r"[^A-Za-z0-9]", "", self.name) or "x"

    @classmethod
    def from_tunnel(cls, name, xyz, tunnel, *, arm: float = 1.0, **kw) -> "PoseSide":
        """Build from a :class:`~channelbridge.tunnel.Tunnel`."""
        return cls(name=name, xyz=xyz, exit=tunnel.exit,
                   axis_point=tunnel.axis_point(arm), **kw)


@dataclass(frozen=True)
class PoseParams:
    """The grid, the filters and the budget."""

    # --- the grid
    dic_min: float = 6.0
    dic_max: float = 36.0
    dic_step: float = 1.0
    oic: Tuple[float, ...] = (180, 165, 150, 135, 120, 105, 90)
    azimuth: Tuple[float, ...] = tuple(range(0, 360, 30))
    tilt: Tuple[float, ...] = (0, 30, 60, 90)
    #: Which way the downstream channel tilts. One value by default; see the
    #: module docstring on what that leaves unsampled.
    tilt_azimuth: Tuple[float, ...] = (0.0,)
    spin: Tuple[float, ...] = tuple(range(0, 360, 30))

    # --- the filters
    clash_cut: float = 3.0
    binder_gap: float = 6.0
    bridge_lo: float = 15.0
    bridge_hi: float = 50.0
    path_clearance: float = 3.5
    path_margin: float = 7.0

    # --- occupancy
    occupancy: str = "auto"
    binder_radius: float = 11.0
    occupancy_gap: float = 4.0

    # --- budget and bookkeeping
    max_ok: int = 50
    grid_spacing: float = 1.0
    oic_arm: float = 1.0
    do_coefficients: DOCoefficients = DO_PUBLISHED

    def __post_init__(self):
        if self.dic_step <= 0:
            raise ValueError("dic_step must be positive, got %r" % (self.dic_step,))
        if self.dic_min > self.dic_max:
            raise ValueError("dic_min %g is above dic_max %g" % (self.dic_min, self.dic_max))
        for angle in self.oic:
            if not 0.0 <= float(angle) <= 180.0:
                raise ValueError("oic %g is not in [0, 180] degrees" % float(angle))
        if self.occupancy not in ("auto", "on", "off"):
            raise ValueError("occupancy must be auto, on or off, got %r" % (self.occupancy,))
        if self.bridge_lo > self.bridge_hi:
            raise ValueError("bridge_lo %g is above bridge_hi %g"
                             % (self.bridge_lo, self.bridge_hi))

    @classmethod
    def from_spec(cls, spec, **overrides) -> "PoseParams":
        """From a pair file's ``pose:`` section (:class:`config.PoseSpec`)."""
        kw = dict(
            dic_min=spec.dic_min, dic_max=spec.dic_max, dic_step=spec.dic_step,
            oic=tuple(spec.oic), azimuth=tuple(spec.azimuth),
            tilt=tuple(spec.tilt), spin=tuple(spec.spin),
            binder_gap=spec.binder_gap, max_ok=spec.max_ok,
        )
        kw.update(overrides)
        return cls(**kw)

    @property
    def dic_values(self) -> List[float]:
        n = math.floor((self.dic_max - self.dic_min) / self.dic_step + 1e-9) + 1
        return [self.dic_min + i * self.dic_step for i in range(max(0, n))]


@dataclass(frozen=True)
class GridPoint:
    d_ic: float
    o_ic: float
    azimuth: float
    tilt: float
    tilt_azimuth: float
    spin: float


def pose_grid(params: PoseParams) -> List[GridPoint]:
    """Every grid point, best DO first.

    Best-DO-first means the scan starts at the tightest, most face-to-face
    geometry the score can reward and walks outward, so the first feasible
    pose is the best one and ``max_ok`` truncates without biasing.

    Degenerate axes are collapsed rather than scanned: ``azimuth`` where the
    O_IC cone is a point (``o_ic = 180``), ``tilt_azimuth`` where there is no
    tilt to steer (``tilt = 0``).
    """
    dics = [d for d in params.dic_values if d >= MIN_DIC]
    if not dics:
        raise ValueError(
            "no usable D_IC: every value is below %.1f A, where the two exits "
            "coincide and the O_IC vertex is degenerate" % MIN_DIC)
    grid = []
    for d in dics:
        for o in params.oic:
            azis = params.azimuth if abs(float(o) - 180.0) > 1e-6 else params.azimuth[:1]
            for az in azis:
                for ti in params.tilt:
                    taz = params.tilt_azimuth if abs(float(ti)) > 1e-6 \
                        else params.tilt_azimuth[:1]
                    for ta in taz:
                        for sp in params.spin:
                            grid.append(GridPoint(float(d), float(o), float(az),
                                                  float(ti), float(ta), float(sp)))
    grid.sort(key=lambda g: (-do_score(g.d_ic, g.o_ic, params.do_coefficients),
                             g.d_ic, -g.o_ic, g.azimuth, g.tilt, g.tilt_azimuth, g.spin))
    return grid


def pose_transform(up: PoseSide, dn: PoseSide, point: GridPoint
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """The rigid-body transform placing *dn* at one grid point.

    Returns ``(R, t)``, applied as ``xyz @ R.T + t``.
    """
    m_a = up.exit
    u_a = up.axis
    frame = _frame(u_a, [0.0, 0.0, 1.0])

    # where the downstream exit goes: swing u_a off itself by (180 - o_ic),
    # about a perpendicular rolled around u_a by the azimuth. The angle the
    # metric measures is at m_a between -u_a and this direction, so swinging
    # by (180 - o_ic) is what makes the realised O_IC the requested one.
    perp = rotation_about(u_a, point.azimuth) @ frame[1]
    w = _unit(rotation_about(perp, 180.0 - point.o_ic) @ u_a, "O_IC direction")
    m_b_target = m_a + w * point.d_ic

    # which way the downstream channel points: back along w at tilt 0, then
    # tilted, with tilt_azimuth choosing the plane it tilts in.
    back = -w
    frame_b = _frame(back, frame[2])
    tilt_axis = rotation_about(back, point.tilt_azimuth) @ frame_b[1]
    u_b_target = _unit(rotation_about(tilt_axis, point.tilt) @ back, "tilt direction")

    # rotate the downstream body frame onto the target one, then spin. Its
    # reference perpendicular is "toward my own centroid", which is a
    # property of the enzyme, so spin 0 means the same thing every run.
    ref_perp = _perp(dn.centroid - dn.exit, dn.axis)
    if float(np.linalg.norm(ref_perp)) < 1e-6:
        # the centroid sits on the channel axis: any perpendicular will do,
        # and spin sweeps the whole circle anyway
        ref_perp = _frame(dn.axis, [0.0, 0.0, 1.0])[1]
    src = _frame(dn.axis, ref_perp)
    dst = _frame(u_b_target, frame[2])
    rot = dst.T @ src
    rot = rotation_about(u_b_target, point.spin) @ rot
    trans = m_b_target - rot @ dn.exit
    return rot, trans


def score_pose(up: PoseSide, dn: PoseSide, rot, trans,
               coeffs: DOCoefficients = DO_PUBLISHED) -> Dict[str, float]:
    """The full geometric report for one placement, in iMARS terms plus extras.

    ``O_IC`` is the published one-sided angle. ``O_IC_rev`` is the same angle
    measured at the *other* exit: the metric does not see it, and a pose can
    aim one channel perfectly while the other points away. See
    ``docs/pipeline.md`` on what the score does not see.
    """
    rot = np.asarray(rot, float).reshape(3, 3)
    trans = np.asarray(trans, float).reshape(3)
    m_a, p1_a = up.exit, up.axis_point
    m_b = rot @ dn.exit + trans
    p1_b = rot @ dn.axis_point + trans

    d = float(np.linalg.norm(m_b - m_a))
    angle = o_ic(p1_a, m_a, m_b)
    angle_rev = o_ic(p1_b, m_b, m_a)
    raw = do_score(d, angle, coeffs)
    conn = m_b - m_a
    phi = _dihedral(_perp(up.centroid - m_a, conn), conn,
                    _perp((rot @ dn.centroid + trans) - m_b, conn))
    return {"D_IC": d, "O_IC": angle, "O_IC_rev": angle_rev,
            "DO": raw, "DO_clamped": max(0.0, min(1.0, raw)),
            "phi": phi, "theta_up": 180.0 - angle, "theta_dn": 180.0 - angle_rev}


# ----------------------------------------------------------- occupancy
#
# An obligate homomultimer presents a binder once per protomer, so a pose
# that works on one protomer has to be realisable on all of them at once:
# that many copies of the downstream enzyme and that many binders have to
# coexist around one particle. Nothing in D_IC, O_IC or the single-pose
# clash test sees this. The single-pose test already covers one downstream
# enzyme against the WHOLE upstream assembly, but never one downstream copy
# against another.

def equivalent_chains(atoms, ref_chain: str, min_identity: float = 0.9) -> List[str]:
    """Chains of *atoms* whose sequence matches *ref_chain*: the protomers.

    Compared residue number by residue number, which is what makes it cheap
    and is safe here: the protomers of a deposited homomultimer are numbered
    alike. A chain that matches on fewer than *min_identity* of the shared
    positions is a different protein that happens to share numbering.
    """
    def seq(ch):
        sel = atoms.select(chain=ch)
        return {"%d%s" % (a.resseq, a.icode.strip()): a.resname
                for a in sel if a.name.strip() == "CA" and not a.hetatm}

    ref = seq(ref_chain)
    if not ref:
        return [ref_chain]
    chains = sorted({a.chain for a in atoms if a.chain.strip()})
    keep = []
    for ch in chains:
        q = seq(ch)
        shared = set(ref) & set(q)
        if not shared:
            continue
        same = sum(1 for r in shared if ref[r] == q[r])
        if len(shared) >= min_identity * len(ref) and same >= min_identity * len(shared):
            keep.append(ch)
    return keep or [ref_chain]


def symmetry_transforms(atoms, ref_chain: str, chains: Sequence[str],
                        *, warn_rmsd: float = 2.0, report=None
                        ) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    """``(chain, R, t)`` mapping the reference protomer's frame onto each copy.

    Derived from the structure's own protomers by superposition, not from an
    assumed symmetry operator: on a clean D2 tetramer this reproduces the
    crystallographic symmetry, and on an imperfect assembly it is still the
    transform that actually relates the copies.
    """
    def ca(ch):
        sel = atoms.select(chain=ch)
        return {"%d%s" % (a.resseq, a.icode.strip()): np.asarray(a.xyz, float)
                for a in sel if a.name.strip() == "CA" and not a.hetatm}

    ref = ca(ref_chain)
    out = [(ref_chain, np.eye(3), np.zeros(3))]
    for ch in chains:
        if ch == ref_chain:
            continue
        q = ca(ch)
        shared = sorted(set(ref) & set(q))
        if len(shared) < 3:
            if report:
                report("occupancy: chain %s shares %d CA with %s, too few to "
                       "superpose; skipped" % (ch, len(shared), ref_chain))
            continue
        p = np.array([ref[r] for r in shared])
        t = np.array([q[r] for r in shared])
        rot, tr, rmsd = superpose(p, t)
        if rmsd > warn_rmsd and report:
            report("occupancy: protomer %s superposes on %s at %.2f A RMSD, so "
                   "the copies will be approximate" % (ref_chain, ch, rmsd))
        out.append((ch, rot, tr))
    return out


def occupancy_check(xyz_posed, binder_a, binder_b, syms, up_xyz, params: PoseParams
                    ) -> Tuple[bool, str, Dict[str, float]]:
    """Can this pose be realised on every protomer at once?

    The binder is modelled as a capsule of ``binder_radius`` about the line
    joining the two epitope centroids. We do not have its structure, and a
    capsule of a minibinder's own width is the honest envelope for "something
    has to run between these two patches".

    A binder is not tested against the two proteins it binds: that overlap is
    the design's job to resolve, and by symmetry it is identical for every
    copy, so it carries no information about occupancy. What is tested is
    every copy against every *other* copy.
    """
    if len(syms) < 2:
        return True, "", {}
    r_b, gap = params.binder_radius, params.occupancy_gap
    copies, axes = [], []
    for _ch, rot, tr in syms:
        copies.append(np.asarray(xyz_posed, float) @ rot.T + tr)
        axes.append((rot @ np.asarray(binder_a, float) + tr,
                     rot @ np.asarray(binder_b, float) + tr))

    m = {"min_copy_copy": float("inf"), "min_binder_binder": float("inf"),
         "min_binder_copy": float("inf"), "min_copy_up": float("inf")}
    for c in copies:
        m["min_copy_up"] = min(m["min_copy_up"], min_interprotein_distance(c, up_xyz))
    for i in range(len(copies)):
        for j in range(i + 1, len(copies)):
            m["min_copy_copy"] = min(
                m["min_copy_copy"], min_interprotein_distance(copies[i], copies[j]))
            m["min_binder_binder"] = min(
                m["min_binder_binder"],
                segment_segment_distance(axes[i][0], axes[i][1], axes[j][0], axes[j][1]))
        for j in range(len(copies)):
            if i != j:
                m["min_binder_copy"] = min(
                    m["min_binder_copy"],
                    segment_cloud_distance(axes[i][0], axes[i][1], copies[j]))

    if m["min_copy_copy"] < gap:
        return False, "downstream copies overlap (%.1f A)" % m["min_copy_copy"], m
    if m["min_binder_binder"] < 2 * r_b:
        return False, ("binders overlap (axes %.1f A apart, need %.1f)"
                       % (m["min_binder_binder"], 2 * r_b)), m
    if m["min_binder_copy"] < r_b + gap:
        return False, ("a binder runs through another copy (%.1f A from its axis, "
                       "need %.1f)" % (m["min_binder_copy"], r_b + gap)), m
    return True, "", m


# -------------------------------------------------------------- the scan

@dataclass
class Pose:
    """One evaluated grid point."""

    pose: int
    status: str
    d_ic: float
    o_ic: float
    o_ic_rev: float
    do: float
    do_clamped: float
    phi: float
    theta_up: float
    theta_dn: float
    requested: GridPoint
    rotation: np.ndarray
    translation: np.ndarray
    min_distance: float = float("nan")
    n_bridgeable: int = 0
    pairs_clear: Tuple[str, ...] = ()
    pairs_blocked: Tuple[str, ...] = ()
    best_pair: str = "-"
    best_pair_tag: str = "nopair"
    best_pair_dist: float = float("nan")
    best_up_patch: Optional[int] = None
    best_dn_patch: Optional[int] = None
    occupancy: Dict[str, float] = field(default_factory=dict)
    occupancy_note: str = ""
    #: ``geomNNN`` for a feasible pose, empty for a rejected one.
    geom: str = ""

    @property
    def feasible(self) -> bool:
        return self.status == "ok"

    def apply(self, xyz) -> np.ndarray:
        return np.asarray(xyz, float) @ self.rotation.T + self.translation

    @property
    def tag(self) -> str:
        """The file-name stem: what this pose is, readable without the table."""
        stem = self.geom or "pose%05d" % self.pose
        return ("%s_DO%.3f_d%05.1f_o%03.0f_%s_az%03.0f_ti%03.0f_sp%03.0f_%s"
                % (stem, self.do_clamped, self.d_ic, self.o_ic, self.best_pair_tag,
                   self.requested.azimuth, self.requested.tilt,
                   self.requested.spin, self.status))

    def describe(self) -> str:
        return ("%s  D_IC %5.1f A  O_IC %5.1f deg  DO %.4f  gap %4.1f A  %s"
                % (self.geom or "pose%05d" % self.pose, self.d_ic, self.o_ic,
                   self.do_clamped,
                   self.min_distance if self.min_distance == self.min_distance else -1.0,
                   self.status))


@dataclass(frozen=True)
class PatchPair:
    """One epitope pair considered at one pose."""

    pose: int
    up_patch: int
    dn_patch: int
    distance: float
    clear: bool
    up_residues: str = ""
    dn_residues: str = ""


@dataclass
class SearchResult:
    """Everything the scan found, and what it did not reach."""

    poses: List[Pose]
    pairs: List[PatchPair]
    up_name: str
    dn_name: str
    params: PoseParams
    n_grid: int = 0
    truncated: bool = False
    occupancy_chains: Tuple[str, ...] = ()

    @property
    def feasible(self) -> List[Pose]:
        return [p for p in self.poses if p.feasible]

    def best(self) -> Optional[Pose]:
        ok = self.feasible
        return min(ok, key=lambda p: (-p.do_clamped, p.d_ic)) if ok else None

    def counts(self) -> Dict[str, int]:
        return {s: sum(1 for p in self.poses if p.status == s) for s in STATUSES
                if any(p.status == s for p in self.poses)}

    def summary(self) -> str:
        ok = self.feasible
        lines = ["%d of %d grid point(s) evaluated, %d feasible%s"
                 % (len(self.poses), self.n_grid, len(ok),
                    " (stopped at max_ok)" if self.truncated else "")]
        for status, n in sorted(self.counts().items(), key=lambda kv: -kv[1]):
            lines.append("  %-20s %d" % (status, n))
        best = self.best()
        if best is None:
            lines.append("no feasible pose: loosen binder_gap, widen bridge_hi, "
                         "raise dic_max, or sweep tilt_azimuth")
        else:
            lines.append("best: %s" % best.describe())
            lines.append("  best epitope pair %s at %.1f A"
                         % (best.best_pair, best.best_pair_dist))
            floor = min(p.d_ic for p in ok)
            lines.append("  tightest feasible D_IC %.1f A (DO %.4f there at O_IC 180)"
                         % (floor, do_score(floor, 180.0, self.params.do_coefficients)))
            if abs(floor - self.params.dic_min) < 1e-6 and len(self.params.dic_values) > 1:
                lines.append("  NOTE that is dic_min, so the real floor is lower; "
                             "re-run lower to find where the two stop fitting")
            lines.extend("  " + line for line in self.tilt_azimuth_note() or [])
        return "\n".join(lines)

    def tilt_azimuth_note(self) -> List[str]:
        """Whether the unswept tilt direction is limiting this answer.

        ``tilt_azimuth`` is the third degree of freedom of the downstream
        enzyme's orientation and defaults to a single value, so saying "sweep
        it if you get stuck" is not advice anyone can act on. The question is
        answerable from the result: the axis is degenerate at ``tilt = 0``,
        so if every feasible pose has no tilt, sweeping it cannot change the
        answer, and if some do, it was searching a fixed tilt direction.

        Measured on the iMARS set: with the default grid every feasible pose
        on those pairs had tilt 0 and sweeping changed nothing at three times
        the scan; on a thinned grid (one azimuth, spin every 60 degrees,
        D_IC every 2 A) the tilted poses mattered and sweeping was worth
        0.063 DO and 8 A of D_IC.
        """
        ok = self.feasible
        if not ok or len(self.params.tilt_azimuth) > 1:
            return []
        tilted = [p for p in ok if abs(p.requested.tilt) > 1e-9]
        if not tilted:
            return ["every feasible pose has tilt 0, where tilt_azimuth is "
                    "degenerate: sweeping it cannot change this answer"]
        return ["NOTE %d of %d feasible pose(s) tilt the downstream channel "
                "and tilt_azimuth was not swept, so only one tilt direction "
                "was tried. It is the cheapest axis to add -- no scored "
                "quantity changes along it." % (len(tilted), len(ok))]


def _number_feasible(poses: Sequence[Pose]) -> List[Pose]:
    """Rank by DO and issue ``geomNNN`` to the feasible poses only.

    Numbering rejections alongside successes made ``geom050145`` the best
    design in a 100k grid, which says nothing useful. A rejected pose keeps
    its raw grid number and no geom id, so it is visibly not a candidate.
    """
    ranked = sorted(poses, key=lambda p: (-p.do_clamped, p.d_ic))
    ok = [p for p in ranked if p.feasible]
    width = max(3, len(str(len(ok))))
    for p in ranked:
        p.geom = ""
    for i, p in enumerate(ok, 1):
        p.geom = "geom%0*d" % (width, i)
    return ranked


def search(up: PoseSide, dn: PoseSide, params: Optional[PoseParams] = None,
           *, report=None, progress: int = 0) -> SearchResult:
    """Scan the grid and return every evaluated pose.

    The upstream enzyme does not move. Filters run cheapest-first and the
    grid is ordered best-DO-first, so ``params.max_ok`` stops the scan at the
    point where nothing better can lie below.
    """
    params = params or PoseParams()
    say = report or (lambda _m: None)

    syms: List[Tuple[str, np.ndarray, np.ndarray]] = []
    if params.occupancy != "off":
        syms = list(up.symmetry)
        if len(syms) > 1:
            say("occupancy: %s has %d protomer(s) (%s), so every feasible pose "
                "is also tested with %s on all of them at once, plus a binder "
                "of radius %.1f A between each epitope pair"
                % (up.name, len(syms), ",".join(s[0] for s in syms), dn.name,
                   params.binder_radius))
        elif params.occupancy == "on":
            raise ValueError(
                "occupancy is 'on' but %s has no equivalent protomers to tile. "
                "There is nothing to test; use 'auto' or 'off'." % up.name)
        else:
            syms = []

    # Built once, because the upstream enzyme never moves.
    clash = _Proximity(up.xyz, params.clash_cut, params.grid_spacing)
    gap = _Proximity(up.xyz, params.binder_gap, params.grid_spacing)
    near = _clearance_grid(up.xyz, NEAR_RADIUS, params.grid_spacing)

    grid = pose_grid(params)
    say("scanning %d grid point(s), best DO first: D_IC %.0f-%.0f A step %.0f, "
        "O_IC %s deg" % (len(grid), params.dic_min, params.dic_max, params.dic_step,
                         ",".join("%g" % o for o in params.oic)))
    say("best DO this grid could reach: %.4f at D_IC %.0f, O_IC %.0f"
        % (do_score(grid[0].d_ic, grid[0].o_ic, params.do_coefficients),
           grid[0].d_ic, grid[0].o_ic))

    poses: List[Pose] = []
    pairs: List[PatchPair] = []
    n_ok = 0
    truncated = False

    for pid, point in enumerate(grid, 1):
        if params.max_ok and n_ok >= params.max_ok:
            truncated = True
            say("stopping after %d feasible pose(s) (max_ok); %d of %d grid "
                "points evaluated. The grid is best-DO-first, so nothing "
                "better lies below." % (n_ok, pid - 1, len(grid)))
            break
        if progress and pid % progress == 0:
            say("  ... %d/%d evaluated, %d feasible" % (pid, len(grid), n_ok))

        rot, trans = pose_transform(up, dn, point)
        moved = dn.xyz @ rot.T + trans
        row = score_pose(up, dn, rot, trans, params.do_coefficients)
        pose = Pose(pose=pid, status="ok", requested=point,
                    rotation=rot, translation=trans,
                    d_ic=row["D_IC"], o_ic=row["O_IC"], o_ic_rev=row["O_IC_rev"],
                    do=row["DO"], do_clamped=row["DO_clamped"], phi=row["phi"],
                    theta_up=row["theta_up"], theta_dn=row["theta_dn"])

        # 1 and 2: clash, then the binder gap. `min_distance` is left NaN on
        # a pose rejected without measurement -- the verdict is exact, the
        # number was never needed, and writing a voxel bound into a column
        # headed "distance" would be a worse answer than leaving it out.
        verdict = clash.test(moved)
        if verdict is True:
            pose.status = "clash"
            poses.append(pose)
            continue
        if verdict is None:
            pose.min_distance = _min_distance(up.xyz, near, moved)
            if pose.min_distance < params.clash_cut:
                pose.status = "clash"
                poses.append(pose)
                continue
        if pose.min_distance != pose.min_distance:          # not measured yet
            if gap.test(moved) is True:
                pose.status = "no_binder_gap"
                poses.append(pose)
                continue
            # either near the gap boundary or clear of it; a pose that gets
            # this far is reported with a real distance either way
            pose.min_distance = _min_distance(up.xyz, near, moved)
        if pose.min_distance < params.binder_gap:
            pose.status = "no_binder_gap"
            poses.append(pose)
            continue

        # 3: a bridgeable epitope pair. Every pair that clears is recorded,
        # not just the best: one geometry usually bridges several, and
        # reporting only the best made it look as though each patch needed
        # its own, lower-DO pose.
        both = np.vstack([up.xyz, moved])
        clear_tags, blocked_tags = [], []
        best_d = float("nan")
        for p1 in up.patches:
            for p2 in dn.patches:
                c2 = rot @ p2.centroid + trans
                d = float(np.linalg.norm(p1.centroid - c2))
                if not (params.bridge_lo <= d <= params.bridge_hi):
                    continue
                clear = path_is_clear(p1.centroid, c2, both,
                                      params.path_clearance, params.path_margin)
                tag = "%s%s-%s%s" % (up.tag_name, p1.label,
                                     dn.tag_name, p2.label)
                if clear:
                    clear_tags.append(tag)
                    if not best_d == best_d or d < best_d:
                        best_d = d
                        pose.best_pair = "%s_%s + %s_%s" % (up.name, p1.label,
                                                            dn.name, p2.label)
                        pose.best_pair_tag = tag
                        pose.best_up_patch, pose.best_dn_patch = p1.rank, p2.rank
                else:
                    blocked_tags.append(tag)
                pairs.append(PatchPair(pose=pid, up_patch=p1.rank, dn_patch=p2.rank,
                                       distance=d, clear=clear,
                                       up_residues="+".join(p1.residues),
                                       dn_residues="+".join(p2.residues)))
        pose.n_bridgeable = len(clear_tags)
        pose.pairs_clear = tuple(clear_tags)
        pose.pairs_blocked = tuple(blocked_tags)
        pose.best_pair_dist = best_d
        if not clear_tags:
            pose.status = "no_bridgeable_pair"
            poses.append(pose)
            continue

        # 4: occupancy, last because it is the dearest and only poses that
        # already passed everything else can reach it.
        if syms:
            up_patch = _patch_by_rank(up.patches, pose.best_up_patch)
            dn_patch = _patch_by_rank(dn.patches, pose.best_dn_patch)
            fits, why, metrics = occupancy_check(
                moved, up_patch.centroid, rot @ dn_patch.centroid + trans,
                syms, up.xyz, params)
            pose.occupancy = metrics
            if not fits:
                pose.status = "occupancy_clash"
                pose.occupancy_note = why
                poses.append(pose)
                continue

        poses.append(pose)
        n_ok += 1

    ranked = _number_feasible(poses)
    result = SearchResult(poses=ranked, pairs=pairs, up_name=up.name, dn_name=dn.name,
                          params=params, n_grid=len(grid), truncated=truncated,
                          occupancy_chains=tuple(s[0] for s in syms))
    say(result.summary())
    return result


def _patch_by_rank(patches: Sequence[PatchSite], rank: Optional[int]) -> PatchSite:
    for p in patches:
        if p.rank == rank:
            return p
    raise KeyError("no patch of rank %r" % (rank,))


# ------------------------------------------------------------ reading in

def read_patch_sites(table, atoms, *, chain: Optional[str] = None,
                     max_patches: Optional[int] = None) -> List[PatchSite]:
    """Patches from a stage-3 ``*_hotspot_patches.tsv``, centroids remeasured.

    The table lists each patch's residues but not its centroid, so the
    centroid is recomputed here as the mean CB (CA for glycine) of those
    residues in *atoms*. That is the same definition stage 3 used, and
    recomputing rather than reading a column means a table from any version
    of stage 3 -- including the original script's -- can be used.
    """
    rows = []
    with open(table, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_rank = header.index("rank")
            i_res = header.index("hotspot_residues")
        except ValueError:
            raise ValueError(
                "%s is not a hotspot patch table: it has no 'rank' and "
                "'hotspot_residues' columns (found %s)"
                % (table, ", ".join(header) or "nothing")) from None
        i_score = header.index("score") if "score" in header else None
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(header) - len(f))
            rows.append((int(f[i_rank]),
                         tuple(r for r in f[i_res].split("+") if r),
                         float(f[i_score]) if i_score is not None and f[i_score] else 0.0))

    sel = atoms.select(chain=chain) if chain else atoms
    cb: Dict[str, np.ndarray] = {}
    ca: Dict[str, np.ndarray] = {}
    for a in sel:
        if a.hetatm:
            continue
        key = "%d%s" % (a.resseq, a.icode.strip())
        name = a.name.strip()
        if name == "CB":
            cb[key] = np.asarray(a.xyz, float)
        elif name == "CA":
            ca[key] = np.asarray(a.xyz, float)

    out = []
    for rank, residues, score in sorted(rows):
        pts = [cb.get(r, ca.get(r)) for r in residues]
        have = [p for p in pts if p is not None]
        if len(have) != len(residues):
            missing = [r for r, p in zip(residues, pts) if p is None]
            raise ValueError(
                "patch %d of %s names residue(s) %s, which are not in the "
                "structure%s. The patch table and the structure are not from "
                "the same prep." % (rank, table, ",".join(missing),
                                    " chain %s" % chain if chain else ""))
        out.append(PatchSite(rank=rank, centroid=np.mean(have, axis=0),
                             residues=residues, score=score))
        if max_patches and len(out) >= max_patches:
            break
    return out


# ------------------------------------------------------------- writing out

_POSE_COLUMNS = (
    "geom_id", "pose", "status", "D_IC", "O_IC", "DO", "DO_clamped", "O_IC_rev",
    "dic_requested", "oic_requested", "azimuth", "tilt", "tilt_azimuth", "spin",
    "theta_upstream", "theta_downstream", "phi", "min_interprotein_dist",
    "binder_gap_ok", "n_bridgeable", "best_pair_dist", "best_pair",
    "pairs_clear", "pairs_blocked", "copy_copy", "binder_binder", "binder_copy",
    "occupancy_note",
)


def _num(value, fmt="%.1f") -> str:
    if value is None:
        return "-"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if f != f or f == float("inf"):
        return "-"
    return fmt % f


def write_poses(path, result: SearchResult) -> str:
    """``poses.tsv``: every evaluated pose, best DO first."""
    lines = ["\t".join(_POSE_COLUMNS)]
    for p in result.poses:
        r = p.requested
        lines.append("\t".join((
            p.geom, str(p.pose), p.status,
            _num(p.d_ic, "%.2f"), _num(p.o_ic), _num(p.do, "%.4f"),
            _num(p.do_clamped, "%.4f"), _num(p.o_ic_rev),
            _num(r.d_ic), _num(r.o_ic), _num(r.azimuth, "%.0f"),
            _num(r.tilt, "%.0f"), _num(r.tilt_azimuth, "%.0f"), _num(r.spin, "%.0f"),
            _num(p.theta_up), _num(p.theta_dn), _num(p.phi),
            _num(p.min_distance), str(p.min_distance >= result.params.binder_gap),
            str(p.n_bridgeable), _num(p.best_pair_dist), p.best_pair,
            "+".join(p.pairs_clear) or "-", "+".join(p.pairs_blocked) or "-",
            _num(p.occupancy.get("min_copy_copy")),
            _num(p.occupancy.get("min_binder_binder")),
            _num(p.occupancy.get("min_binder_copy")),
            p.occupancy_note)))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return str(path)


def write_patch_pairs(path, result: SearchResult) -> str:
    """``patch_pairs.tsv``: every epitope pair each pose could reach."""
    geom_of = {p.pose: p.geom for p in result.poses}
    lines = ["\t".join(("geom_id", "pose", "%s_patch" % result.up_name,
                        "%s_patch" % result.dn_name, "centroid_dist", "path_clear"))]
    for r in sorted(result.pairs, key=lambda q: (geom_of.get(q.pose) or "zzz", q.distance)):
        lines.append("\t".join((geom_of.get(r.pose, ""), str(r.pose),
                                str(r.up_patch), str(r.dn_patch),
                                "%.1f" % r.distance, str(r.clear))))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return str(path)


def write_patch_pair_best(path, result: SearchResult) -> str:
    """``patch_pair_best.tsv``: the best geometry for each epitope pair.

    ``patch_pairs.tsv`` lists every combination a pose could reach, but a
    pose's FILE NAME carries only its own best pair, so reading names
    suggests a patch's best geometry is the first geom named after it --
    when an earlier, higher-DO geom may bridge that pair too. Most listed
    pairs are also ``path_clear False``: candidates, not routes. This table
    answers the question directly, clear pairs only.
    """
    geom_of = {p.pose: p.geom for p in result.poses}
    do_of = {p.pose: p.do_clamped for p in result.poses}
    best: Dict[Tuple[int, int], Tuple[str, float, int]] = {}
    seen: Dict[Tuple[int, int], int] = {}
    for r in result.pairs:
        if not r.clear or not geom_of.get(r.pose):
            continue
        key = (r.up_patch, r.dn_patch)
        seen[key] = seen.get(key, 0) + 1
        row = (geom_of[r.pose], r.distance, r.pose)
        if key not in best or row[0] < best[key][0]:
            best[key] = row
    lines = ["\t".join(("%s_patch" % result.up_name, "%s_patch" % result.dn_name,
                        "best_geom", "pose", "centroid_dist", "n_geoms_clear", "DO"))]
    for key in sorted(best):
        geom, dist, pid = best[key]
        lines.append("\t".join((str(key[0]), str(key[1]), geom, str(pid),
                                "%.1f" % dist, str(seen[key]),
                                "%.4f" % do_of.get(pid, 0.0))))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return str(path)


def write_complex(path, up: PoseSide, dn: PoseSide, pose: Pose) -> str:
    """One pose as a PDB: the upstream enzyme, plus the downstream one moved.

    Needs both sides' ``atoms``. Chain ids are left exactly as prepped, and a
    collision is refused rather than silently merged -- two enzymes in one
    chain is not a complex anything downstream can read.
    """
    import dataclasses

    from .io.pdb import write_pdb

    if up.atoms is None or dn.atoms is None:
        raise ValueError("write_complex needs both sides' atoms; build the "
                         "PoseSide with atoms=read_pdb(...)")
    up_chains = {a.chain for a in up.atoms}
    dn_chains = {a.chain for a in dn.atoms}
    clash = up_chains & dn_chains
    if clash:
        raise ValueError(
            "%s and %s both use chain %s, so the two would merge into one "
            "chain. Re-prep the downstream enzyme with a different chain id "
            "(`chbridge prep --relabel`)."
            % (up.name, dn.name, "+".join(sorted(clash))))
    moved = []
    for a in dn.atoms:
        xyz = pose.rotation @ np.asarray(a.xyz, float) + pose.translation
        moved.append(dataclasses.replace(a, xyz=tuple(float(v) for v in xyz)))
    write_pdb(path, list(up.atoms) + moved)
    return str(path)


def write_sessions(out_dir, up: PoseSide, dn: PoseSide, result: SearchResult,
                   *, which: str = "top", top_n: int = 8, prefix: str = "",
                   report=None) -> List[str]:
    """Write a complex PDB per chosen pose. Returns the paths written.

    ``which`` is ``top`` (the best *top_n* feasible poses), ``success``
    (every feasible one) or ``all`` (including rejections, for inspection).
    """
    say = report or (lambda _m: None)
    if which == "all":
        chosen: Iterable[Pose] = result.poses
    elif which == "success":
        chosen = result.feasible
    elif which == "top":
        chosen = result.feasible[:top_n]
    else:
        raise ValueError("which must be top, success or all, got %r" % (which,))
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for pose in chosen:
        name = "%s%s_complex.pdb" % (prefix, pose.tag)
        path = os.path.join(out_dir, name)
        write_complex(path, up, dn, pose)
        written.append(path)
        say("wrote %s" % name)
    return written
