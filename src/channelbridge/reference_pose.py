"""A reference placement of the two enzymes, computed without a pose search.

The hotspot map is shown **before** the DO sampling, so that the epitopes can
be judged before committing to a search. That rules out taking the relative
position from the search's own output: at this point there is none.

So the placement is constructed from the two channels alone:

* O_IC is one-sided -- the angle P1-P2-P3 at the **upstream** exit -- so
  ``O_IC = 180`` says only that the downstream exit lies on the ray leaving
  the upstream channel. That fixes P3 to ``P2 + d * u`` and nothing else.
* The downstream enzyme's own orientation is therefore **free**: O_IC says
  nothing about it. It is scanned, over the same tilt/azimuth/spin grid the
  pose search uses, with "channel aimed straight back at the upstream exit"
  as tilt 0. Constraining it to that one orientation is tempting and wrong --
  on the p6a pair it costs 10 A of D_IC (18 A instead of 8 A) for a
  requirement the metric does not make.
* ``d`` is raised from the minimum until a spin exists that does not make the
  two enzymes interpenetrate, and the first such ``(d, spin)`` is the
  placement. That is "the best clash-free geometry at O_IC = 180": DO falls
  monotonically with d, so the smallest clash-free d is the highest-DO one.

This is deliberately *not* the pose search. There is no bridgeable-pair test,
no binder-gap test and no occupancy test -- those are what the search is for,
and a figure that refused to be drawn because no binder fits would hide
exactly the fact worth seeing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np

#: How big the dilation ball may get, in voxels. See :func:`occupancy_grid`.
MAX_SHELL_VOXELS = 50_000

__all__ = [
    "Placement",
    "clash_free_distance",
    "clear_view_about_axis",
    "min_interprotein_distance",
    "occupancy_grid",
    "place_at_oic180",
    "rotation_about",
    "rotation_taking",
]


def rotation_taking(a, b) -> np.ndarray:
    """The shortest rotation taking unit vector *a* onto unit vector *b*."""
    a = np.asarray(a, float).reshape(3)
    b = np.asarray(b, float).reshape(3)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-9 or nb < 1e-9:
        raise ValueError("cannot rotate a zero-length vector")
    a, b = a / na, b / nb
    v = np.cross(a, b)
    c = float(a @ b)
    s = float(np.linalg.norm(v))
    if s < 1e-12:
        if c > 0:
            return np.eye(3)
        # antiparallel: a half turn about any axis perpendicular to a
        tmp = np.array([1.0, 0.0, 0.0])
        if abs(float(a @ tmp)) > 0.9:
            tmp = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, tmp)
        return rotation_about(axis / np.linalg.norm(axis), 180.0)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], float)
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def rotation_about(axis, degrees: float) -> np.ndarray:
    """Right-handed rotation about *axis* by *degrees* (Rodrigues)."""
    axis = np.asarray(axis, float).reshape(3)
    n = float(np.linalg.norm(axis))
    if n < 1e-9:
        raise ValueError("rotation axis has zero length")
    k = axis / n
    th = math.radians(float(degrees))
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]], float)
    return np.eye(3) + math.sin(th) * kx + (1 - math.cos(th)) * (kx @ kx)


def occupancy_grid(xyz, clearance: float, spacing: float = 1.0):
    """A boolean voxel grid of everything within *clearance* of an atom.

    Used instead of an all-pairs distance for the placement scan: hundreds of
    trial placements against thousands of atoms is tens of millions of
    distances, while a grid lookup is one array index per atom. It is
    conservative to within one voxel, so the placement it returns is checked
    exactly afterwards by :func:`min_interprotein_distance`.

    Returns ``(grid, origin, spacing)``.
    """
    xyz = np.atleast_2d(np.asarray(xyz, float))
    if xyz.size == 0 or xyz.shape[-1] != 3:
        raise ValueError("expected (N, 3) coordinates, got %r" % (xyz.shape,))
    spacing = float(spacing)
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    pad = float(clearance) + spacing
    origin = xyz.min(axis=0) - pad
    extent = xyz.max(axis=0) + pad - origin
    shape = tuple(math.ceil(e / spacing) + 1 for e in extent)
    grid = np.zeros(shape, dtype=bool)
    idx = np.rint((xyz - origin) / spacing).astype(int)
    grid[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    # grow by the clearance, as a ball in voxels
    r = math.ceil(float(clearance) / spacing)
    if r > MAX_SHELL_VOXELS ** (1 / 3.0) * 0.62:
        # the ball is built as one array roll per offset, so its cost is
        # cubic in the radius. A clash cut-off is a few angstroms; anything
        # that would need a ball this big is a misuse, not a slow run.
        raise ValueError(
            "clearance %.1f A at spacing %.1f A needs a shell of radius %d "
            "voxels, which is too large to build; clearance is a clash "
            "cut-off (a few A), not a search radius"
            % (float(clearance), spacing, r))
    if r > 0:
        offs = np.arange(-r, r + 1)
        ball = [(i, j, k) for i in offs for j in offs for k in offs
                if i * i + j * j + k * k <= r * r]
        grown = np.zeros_like(grid)
        for i, j, k in ball:
            grown |= np.roll(np.roll(np.roll(grid, i, 0), j, 1), k, 2)
        grid = grown
    return grid, origin, spacing


def _hits(grid, origin, spacing, xyz) -> int:
    """How many of *xyz* fall in an occupied voxel."""
    idx = np.rint((np.asarray(xyz, float) - origin) / spacing).astype(int)
    inside = np.ones(len(idx), bool)
    for ax in range(3):
        inside &= (idx[:, ax] >= 0) & (idx[:, ax] < grid.shape[ax])
    if not inside.any():
        return 0
    sel = idx[inside]
    return int(grid[sel[:, 0], sel[:, 1], sel[:, 2]].sum())


def min_interprotein_distance(a, b, chunk: int = 2048) -> float:
    """Closest heavy-atom approach between two coordinate sets, exactly."""
    a = np.atleast_2d(np.asarray(a, float))
    b = np.atleast_2d(np.asarray(b, float))
    if a.size == 0 or b.size == 0:
        return float("inf")
    best = float("inf")
    for i in range(0, len(a), chunk):
        block = a[i:i + chunk]
        d = np.sqrt(((block[:, None, :] - b[None, :, :]) ** 2).sum(-1))
        best = min(best, float(d.min()))
    return best


@dataclass(frozen=True)
class Placement:
    """Where the downstream enzyme goes, and what that placement is."""

    rotation: np.ndarray
    translation: np.ndarray
    d_ic: float
    o_ic: float
    spin: float
    min_distance: float
    do: float
    tilt: float = 0.0
    azimuth: float = 0.0

    def apply(self, xyz) -> np.ndarray:
        return np.asarray(xyz, float) @ self.rotation.T + self.translation

    def describe(self) -> str:
        return ("D_IC %.1f A, O_IC %.0f deg, DO %.4f; tilt %.0f / azimuth %.0f "
                "/ spin %.0f deg, closest approach %.2f A"
                % (self.d_ic, self.o_ic, self.do, self.tilt, self.azimuth,
                   self.spin, self.min_distance))


def place_at_oic180(up_xyz, up_exit, up_axis_point, dn_xyz, dn_exit, dn_axis_point,
                    *, d_min: float = 6.0, d_max: float = 36.0, d_step: float = 1.0,
                    spins: Optional[Sequence[float]] = None,
                    tilts: Optional[Sequence[float]] = None,
                    azimuths: Optional[Sequence[float]] = None,
                    clash_cut: float = 3.0,
                    grid_spacing: float = 1.0, do_coefficients=None,
                    report=None) -> Optional[Placement]:
    """The smallest clash-free D_IC with the two channels collinear.

    ``up_axis_point`` and ``dn_axis_point`` are the points one O_IC arm inside
    each channel (``Tunnel.axis_point``); the exits are ``Tunnel.exit``. The
    upstream enzyme does not move. ``tilts``/``azimuths``/``spins`` default to
    the pose search's own grid, so the answer is comparable with its output.

    Returns ``None`` when no spin at any ``d`` up to ``d_max`` is clash-free,
    which is a real answer about the pair rather than an error.
    """
    from .geometry import do_score

    up_xyz = np.atleast_2d(np.asarray(up_xyz, float))
    dn_xyz = np.atleast_2d(np.asarray(dn_xyz, float))
    up_exit = np.asarray(up_exit, float).reshape(3)
    dn_exit = np.asarray(dn_exit, float).reshape(3)
    u = up_exit - np.asarray(up_axis_point, float).reshape(3)
    nu = float(np.linalg.norm(u))
    if nu < 1e-9:
        raise ValueError("the upstream channel has no direction at its exit")
    u = u / nu
    v = dn_exit - np.asarray(dn_axis_point, float).reshape(3)
    nv = float(np.linalg.norm(v))
    if nv < 1e-9:
        raise ValueError("the downstream channel has no direction at its exit")
    v = v / nv

    # "aimed straight back at the upstream exit" is the tilt-0 orientation;
    # every other orientation is a tilt away from it, turned by an azimuth
    # about the shared axis, and then spun about the downstream channel.
    back = rotation_taking(v, -u)
    if spins is None:
        spins = list(range(0, 360, 30))
    if tilts is None:
        tilts = [0.0, 30.0, 60.0, 90.0]
    if azimuths is None:
        azimuths = list(range(0, 360, 30))
    # two axes perpendicular to u, to tilt within
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(float(u @ tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(u, tmp)
    e1 /= np.linalg.norm(e1)
    grid, origin, spacing = occupancy_grid(up_xyz, clash_cut, grid_spacing)

    orientations = []
    for tilt in tilts:
        for az in (azimuths if float(tilt) != 0.0 else [0.0]):
            tilt_axis = rotation_about(u, float(az)) @ e1
            aim = rotation_about(tilt_axis, float(tilt)) @ back
            for spin in spins:
                orientations.append((float(tilt), float(az), float(spin),
                                     rotation_about(aim @ v, float(spin)) @ aim))
    if report:
        report("scanning %d orientation(s) per D_IC step" % len(orientations))

    d = float(d_min)
    while d <= float(d_max) + 1e-9:
        p3 = up_exit + d * u
        best = None
        for tilt, az, spin, rot in orientations:
            trans = p3 - rot @ dn_exit
            moved = dn_xyz @ rot.T + trans
            n_hit = _hits(grid, origin, spacing, moved)
            if best is None or n_hit < best[0]:
                best = (n_hit, tilt, az, spin, rot, trans, moved)
            if n_hit == 0:
                break
        n_hit, tilt, az, spin, rot, trans, moved = best
        if report:
            report("  D_IC %5.1f A: best tilt %3.0f az %3.0f spin %3.0f, "
                   "%d atom(s) inside the clash shell" % (d, tilt, az, spin, n_hit))
        if n_hit == 0:
            exact = min_interprotein_distance(up_xyz, moved)
            if exact >= clash_cut:
                return Placement(
                    rotation=rot, translation=trans, d_ic=d, o_ic=180.0,
                    tilt=tilt, azimuth=az, spin=spin, min_distance=exact,
                    do=float(do_score(d, 180.0, do_coefficients)
                             if do_coefficients is not None else do_score(d, 180.0)),
                )
        d += float(d_step)
    return None


def clash_free_distance(*args, **kwargs) -> Optional[float]:
    """Just the D_IC of :func:`place_at_oic180`, for a quick report."""
    p = place_at_oic180(*args, **kwargs)
    return None if p is None else p.d_ic


def patch_outward_normal(protein_xyz, centroid, radius: float = 12.0) -> np.ndarray:
    """Which way a binder on this patch would stick out.

    Away from the local protein mass, not from the whole enzyme's centroid: on
    an elongated or concave enzyme the two differ by tens of degrees, and the
    local one is what a binder actually sees. Falls back to the whole-protein
    direction when the patch has no neighbourhood to speak of.
    """
    xyz = np.atleast_2d(np.asarray(protein_xyz, float))
    c = np.asarray(centroid, float).reshape(3)
    d = xyz - c
    near = (d ** 2).sum(1) <= float(radius) ** 2
    ref = xyz[near] if near.sum() >= 10 else xyz
    out = c - ref.mean(axis=0)
    n = float(np.linalg.norm(out))
    if n < 1e-9:
        out = c - xyz.mean(axis=0)
        n = float(np.linalg.norm(out))
        if n < 1e-9:
            return np.array([0.0, 0.0, 1.0])
    return out / n


def binder_blobs(centroid, direction, length: float = 6.0, step: float = 2.0,
                 radius: float = 11.0) -> Tuple[np.ndarray, np.ndarray]:
    """The spheres a minimum-length binder on this patch would occupy.

    The pose stage models a binder as a capsule of ``radius`` between two
    epitope centroids and needs ``binder_gap`` (6 A by default) of room for
    it. Drawing that same minimum span outward from each patch -- as the same
    chain of spheres every ``step`` A that the pose session uses -- shows how
    far a binder would reach from each epitope, on the same terms, before any
    pose has been chosen.
    """
    c = np.asarray(centroid, float).reshape(3)
    direction = np.asarray(direction, float).reshape(3)
    n = float(np.linalg.norm(direction))
    if n < 1e-9:
        raise ValueError("binder direction has zero length")
    direction = direction / n
    if float(length) < 0:
        raise ValueError("binder length must not be negative")
    count = max(2, int(float(length) / float(step)) + 1)
    pts = np.array([c + direction * (float(length) * f)
                    for f in np.linspace(0.0, 1.0, count)])
    return pts, np.full(len(pts), float(radius))


def clear_view_about_axis(xyz_a, xyz_b, axis, *, n: int = 36,
                          grid: float = 3.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Camera axes that put *axis* across the screen and the two sets apart.

    The channel axis is laid along screen x, so the upstream enzyme is on the
    left and the downstream on the right -- the house order. That leaves one
    free angle, the roll about the axis, and the default (whatever PyMOL
    happens to have) routinely puts one enzyme in front of the other. So the
    roll is scanned and the one whose projections overlap least is kept.

    Returns ``(x, y, z)`` camera axes in model coordinates, with z pointing
    at the viewer.
    """
    a = np.atleast_2d(np.asarray(xyz_a, float))
    b = np.atleast_2d(np.asarray(xyz_b, float))
    x = np.asarray(axis, float).reshape(3)
    nx = float(np.linalg.norm(x))
    if nx < 1e-9:
        raise ValueError("the view axis has zero length")
    x = x / nx
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(float(x @ tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(x, tmp)
    e1 /= np.linalg.norm(e1)
    best = None
    for k in range(max(1, int(n))):
        roll = 360.0 * k / max(1, int(n))
        z = rotation_about(x, roll) @ e1            # toward the viewer
        y = np.cross(z, x)
        pa = np.stack([a @ x, a @ y], 1)
        pb = np.stack([b @ x, b @ y], 1)
        cells_a = {tuple(v) for v in np.floor(pa / grid).astype(int)}
        cells_b = {tuple(v) for v in np.floor(pb / grid).astype(int)}
        overlap = len(cells_a & cells_b) / max(1, min(len(cells_a), len(cells_b)))
        if best is None or overlap < best[0]:
            best = (overlap, x.copy(), y.copy(), z.copy())
    return best[1], best[2], best[3]
