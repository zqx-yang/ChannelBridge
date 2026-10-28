"""Geometry for the cross-section and the binder envelope.

Kept apart from the PyMOL code so it can be tested without a graphics
session: these are the parts where an error is silent in a picture.

**The cross-section cuts by deleting atoms, not with a clipping plane.** A
clipping plane also cuts the tunnel spheres: a thick slab leaves protein
covering the channel, a thin one turns the spheres into hollow rings. So a
*copy* of each enzyme has the near half removed and is drawn solid, while the
original stays in place at high transparency. That is what keeps the enzyme
from looking broken into pieces -- the removed half is still faintly there.
"""

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np

__all__ = [
    "capsule_points",
    "cut_mask",
    "enzyme_section_point",
    "parse_hotspot_residues",
    "patch_centroid",
    "section_plane",
]


def section_plane(p2, p3, view_direction):
    """The shared cutting orientation for a channel figure.

    Returns ``(point, normal)``: a plane through the midpoint of the two
    channel exits, facing the camera.

    **Use this for the normal, not usually for the point.** One plane at the
    midpoint cuts both enzymes only if they sit at the same depth along the
    view axis, and they generally do not: on the worked GDH-FNR pose the shared
    plane left GDH with 1296 of 1971 atoms but FNR with 234 of 2352 -- 90 %
    of the downstream enzyme deleted, leaving a sliver that reads as no
    section at all. :func:`enzyme_section_point` places each enzyme's plane at
    its own channel instead.
    """
    p2 = np.asarray(p2, float).reshape(3)
    p3 = np.asarray(p3, float).reshape(3)
    n = np.asarray(view_direction, float).reshape(3)
    norm = float(np.linalg.norm(n))
    if norm < 1e-9:
        raise ValueError("view_direction has zero length; no cutting plane is defined")
    return (p2 + p3) / 2.0, n / norm


def enzyme_section_point(exit_xyz, coords, normal, keep_min: float = 0.45,
                         keep_max: float = 0.75):
    """Where to cut ONE enzyme so the result reads as a section.

    The plane starts at that enzyme's own channel exit and is then clamped so
    that the fraction of atoms left standing lies in
    ``[keep_min, keep_max]``. Both bounds are needed, and each fixes a real
    failure seen on this data:

    * **too little left.** A plane at the midpoint between the two exits cut
      GDH to 1296/1971 atoms but FNR to 234/2352 -- 90 % of the downstream
      enzyme gone, leaving a sliver that reads as no section at all.
    * **too much left.** Cutting at an exit that lies beyond the enzyme's own
      atoms removes nothing, and a section that removes nothing is not a
      section -- it is the uncut enzyme with extra objects in the file.

    Within those bounds the plane sits as close to the exit as it can, so the
    cut passes through the channel the figure is about.
    """
    xyz = np.atleast_2d(np.asarray(coords, float))
    if xyz.size == 0:
        raise ValueError("no coordinates to section")
    n = np.asarray(normal, float).reshape(3)
    nn = float(np.linalg.norm(n))
    if nn < 1e-9:
        raise ValueError("normal has zero length")
    n = n / nn
    keep_min, keep_max = float(keep_min), float(keep_max)
    for name, v in (("keep_min", keep_min), ("keep_max", keep_max)):
        if not 0.0 < v < 1.0:
            raise ValueError("%s must be a fraction in (0, 1), got %r" % (name, v))
    if keep_min > keep_max:
        raise ValueError("keep_min %r is above keep_max %r" % (keep_min, keep_max))
    exit_xyz = np.asarray(exit_xyz, float).reshape(3)
    depth = xyz @ n
    # Deeper plane -> more removed. quantile(keep) is the depth with exactly
    # that fraction of atoms behind it, so keep_min gives the deepest plane
    # allowed and keep_max the shallowest.
    deepest = float(np.quantile(depth, keep_min))
    shallowest = float(np.quantile(depth, keep_max))
    at_exit = float(exit_xyz @ n)
    chosen = min(max(at_exit, deepest), shallowest)
    return exit_xyz + n * (chosen - at_exit)


def cut_mask(coords, point, normal, offset: float = 0.0) -> np.ndarray:
    """Which atoms sit on the camera side of the plane, and so are removed.

    ``offset`` moves the plane along its normal: positive cuts *deeper* into
    the structure (removes more), negative leaves more material in front. The
    sign is fixed this way so that "cut deeper" is a positive number in a
    config file.
    """
    xyz = np.atleast_2d(np.asarray(coords, float))
    if xyz.size == 0:
        return np.zeros(0, dtype=bool)
    if xyz.shape[-1] != 3:
        raise ValueError("coords must be (N, 3), got %r" % (xyz.shape,))
    point = np.asarray(point, float).reshape(3)
    normal = np.asarray(normal, float).reshape(3)
    nn = float(np.linalg.norm(normal))
    if nn < 1e-9:
        raise ValueError("normal has zero length")
    normal = normal / nn
    return ((xyz - point) @ normal) > -float(offset)


def capsule_points(a, b, radius: float, n_along: int = 24, n_ring: int = 16):
    """A capsule between two points, as a cloud that renders as an envelope.

    This stands in for the binder that does not exist yet: the pose stage
    models one as a capsule of ``--binder-radius`` between the two epitope
    centroids and rejects poses whose capsules collide. Drawing the same
    shape shows what the occupancy check actually tested.

    Returns ``(points, radii)`` suitable for one sphere per point.
    """
    a = np.asarray(a, float).reshape(3)
    b = np.asarray(b, float).reshape(3)
    if radius <= 0:
        raise ValueError("capsule radius must be positive, got %r" % (radius,))
    axis = b - a
    length = float(np.linalg.norm(axis))
    if length < 1e-6:
        # degenerate: a sphere at a
        return a[None, :], np.array([radius])
    u = axis / length
    # any two vectors perpendicular to u
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(float(u @ tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    v = np.cross(u, tmp)
    v /= np.linalg.norm(v)
    w = np.cross(u, v)

    n_along = max(2, int(n_along))
    n_ring = max(3, int(n_ring))
    pts = []
    # the shaft, as rings of small spheres on the surface
    small = max(0.4, radius * 0.18)
    for t in np.linspace(0.0, 1.0, n_along):
        c = a + u * (t * length)
        for k in range(n_ring):
            th = 2 * np.pi * k / n_ring
            pts.append(c + radius * (np.cos(th) * v + np.sin(th) * w))
    # the two caps
    for end, sign in ((a, -1.0), (b, 1.0)):
        for phi in np.linspace(0.15, np.pi / 2, 4):
            r = radius * np.cos(phi)
            off = sign * radius * np.sin(phi)
            for k in range(n_ring):
                th = 2 * np.pi * k / n_ring
                pts.append(end + u * off + r * (np.cos(th) * v + np.sin(th) * w))
    pts = np.asarray(pts, float)
    return pts, np.full(len(pts), small)


def parse_hotspot_residues(text) -> List[int]:
    """Residue numbers from a hotspot-patch cell such as ``33+49+58+59``.

    The pose stage writes patches this way. Tolerates spaces, commas and
    stray text, and keeps the order, because a patch is reported as a set
    rather than a range.
    """
    if text is None:
        return []
    s = str(text).strip()
    if not s:
        return []
    out, seen = [], set()
    tok = ""
    for ch in s + " ":
        if ch.isdigit():
            tok += ch
        else:
            if tok:
                v = int(tok)
                if v not in seen:
                    seen.add(v)
                    out.append(v)
                tok = ""
    return out


def patch_centroid(coords_by_resi, residues: Sequence[int]) -> Optional[np.ndarray]:
    """Centroid of the residues of a patch that are present in the structure.

    Returns ``None`` when none of them are, rather than a meaningless origin:
    an epitope that is not in the rendered chain set should disappear from the
    figure, not collapse to (0, 0, 0).
    """
    pts = [coords_by_resi[r] for r in residues if r in coords_by_resi]
    if not pts:
        return None
    return np.asarray(pts, float).mean(axis=0)
