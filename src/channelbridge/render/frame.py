"""Cameras for channel figures, and which poses get a .pse or an SVG.

The default frame, ``rotated0``, is the **angle frame** of the p5c figures: it
puts the iMARS construction in the plane of the page so the figure can be
measured off the image.

    P2 at the centre of the view
    the ray P2->P1 along the horizontal
    the plane P1-P2-P3 face-on

so the segment P2-P3 is drawn at the *true* O_IC and its length is the *true*
D_IC, then the whole thing is turned 180 degrees about the vertical so the
**upstream** enzyme sits on the LEFT (p5c had green on the right once, because
"x axis = O_IC 0" was taken literally, and that was not what was wanted).
After the flip the segment P2-P3 rises at ``180 - O_IC`` from the horizontal,
which is the arithmetic the SVG's angle arc relies on.

Only the camera moves. The coordinates of a posed complex are never
transformed: a figure must stay comparable with the PDB the pose stage wrote.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

__all__ = [
    "Frame",
    "angle_frame",
    "parse_selection",
    "pymol_view",
    "rotated",
    "select_poses",
    "spread_in_view",
]


@dataclass(frozen=True)
class Frame:
    """A camera: three orthonormal axes in model coordinates, plus a centre.

    ``axes`` rows are the camera x (right), y (up) and z (toward the viewer)
    in model coordinates. PyMOL stores them column-wise, which
    :func:`pymol_view` handles -- getting that transpose wrong is the classic
    way to end up with a mirrored figure.
    """

    axes: np.ndarray
    centre: np.ndarray
    name: str = "rotated0"

    def __post_init__(self):
        a = np.asarray(self.axes, dtype=float).reshape(3, 3)
        c = np.asarray(self.centre, dtype=float).reshape(3)
        if not np.allclose(a @ a.T, np.eye(3), atol=1e-6):
            raise ValueError("camera axes are not orthonormal")
        if not np.isclose(np.linalg.det(a), 1.0, atol=1e-6):
            raise ValueError(
                "camera axes are left-handed (det %.3f); the figure would be "
                "mirrored" % np.linalg.det(a)
            )
        object.__setattr__(self, "axes", a)
        object.__setattr__(self, "centre", c)

    @property
    def right(self) -> np.ndarray:
        return self.axes[0]

    @property
    def up(self) -> np.ndarray:
        return self.axes[1]

    @property
    def toward_viewer(self) -> np.ndarray:
        return self.axes[2]

    def project(self, points) -> np.ndarray:
        """Model points to ``(right, up)`` image coordinates, centre removed."""
        p = np.atleast_2d(np.asarray(points, dtype=float))
        rel = p - self.centre
        return np.column_stack([rel @ self.right, rel @ self.up])


def _unit(v, what):
    v = np.asarray(v, dtype=float).reshape(3)
    n = float(np.linalg.norm(v))
    if n < 1e-9:
        raise ValueError("%s has zero length; the frame is undefined" % what)
    return v / n


def angle_frame(p1, p2, p3, *, centre=None, upstream_left=True) -> Frame:
    """The ``rotated0`` angle frame from the iMARS construction.

    Parameters
    ----------
    p1, p2, p3:
        The iMARS points: P1 one O_IC arm inside the upstream exit, P2 the
        upstream exit, P3 the downstream exit.
    centre:
        What the camera looks at. Defaults to the midpoint of P2 and P3, so
        both mouths sit in frame; pass the complex centroid to show more of
        the proteins.
    upstream_left:
        Apply the 180-degree turn that puts the upstream enzyme on the left.
        Switching it off reproduces the raw "P2->P1 along +x" frame.
    """
    p1 = np.asarray(p1, dtype=float).reshape(3)
    p2 = np.asarray(p2, dtype=float).reshape(3)
    p3 = np.asarray(p3, dtype=float).reshape(3)

    x = _unit(p1 - p2, "P2->P1 (the O_IC arm)")
    y_in_plane = p3 - p2
    normal = np.cross(x, y_in_plane)
    if float(np.linalg.norm(normal)) < 1e-6:
        raise ValueError(
            "P1, P2 and P3 are collinear (O_IC is 0 or 180 degrees), so the "
            "plane they define is undefined. Pass an explicit frame, or use "
            "rotated() from a frame built on another pose."
        )
    z = _unit(normal, "plane normal")
    y = np.cross(z, x)

    if upstream_left:
        # 180 degrees about the vertical: x -> -x, z -> -z, y unchanged.
        x, z = -x, -z

    c = np.asarray(p2 + p3, dtype=float) / 2.0 if centre is None else np.asarray(centre, float)
    return Frame(np.array([x, y, z]), c, name="rotated0")


def rotated(frame: Frame, degrees: float) -> Frame:
    """*frame* turned *degrees* about its own vertical axis.

    The naming convention is ``rotated<N>`` with N relative to ``rotated0``,
    always an integer.
    """
    th = np.radians(float(degrees))
    up = frame.up
    k = np.array(
        [[0, -up[2], up[1]], [up[2], 0, -up[0]], [-up[1], up[0], 0]], dtype=float
    )
    r = np.eye(3) + np.sin(th) * k + (1 - np.cos(th)) * (k @ k)
    name = "rotated0" if round(degrees) % 360 == 0 else "rotated%d" % round(degrees)
    return Frame((r @ frame.axes.T).T, frame.centre, name=name)


def pymol_view(frame: Frame, *, slab: float = 3000.0, distance: float = -150.0) -> List[float]:
    """The 18 numbers ``cmd.set_view`` wants.

    PyMOL stores the rotation **column-wise**, so the axes matrix is
    transposed before flattening. The caller normally follows this with
    ``cmd.zoom(sel, buffer=0, complete=0)``, which keeps the rotation and
    fits the selection, then ``cmd.clip('slab', 3000)``.
    """
    v = list(frame.axes.T.flatten())
    v += [0.0, 0.0, float(distance)]          # camera position
    v += [float(c) for c in frame.centre]      # origin of rotation
    v += [1.0, float(slab), -20.0]             # near, far, orthoscopic flag slot
    return v


# --------------------------------------------------------------- what to draw


def parse_selection(spec: Optional[str], available: Sequence[str]) -> List[str]:
    """Turn a user's ``--pse`` / ``--svg`` value into a list of pose ids.

    Accepts, in order of precedence:

    ====================  ==========================================
    ``None`` or ``""``    nothing
    ``"all"``             every pose, best DO first
    an integer ``"8"``    the best 8 by DO
    a list ``"g1,g3"``    exactly those ids, in the order given
    ====================  ==========================================

    *available* must already be ordered best-DO-first, which is how the pose
    stage writes its geom ids. The counts for .pse and SVG are independent on
    purpose: a reviewer wants many sessions to click through and one or two
    publication SVGs.
    """
    if spec is None:
        return []
    s = str(spec).strip()
    if not s:
        return []
    if s.lower() == "all":
        return list(available)
    if s.lower() == "none":
        return []
    try:
        n = int(s)
    except ValueError:
        # not an integer, so treat it as a comma-separated list of ids
        wanted = [t.strip() for t in s.split(",") if t.strip()]
        unknown = [w for w in wanted if w not in available]
        if unknown:
            raise ValueError(
                "unknown pose id(s) %s; available: %s"
                % (", ".join(unknown), ", ".join(available) or "(none)")
            ) from None
        return wanted
    if n < 0:
        raise ValueError("a count cannot be negative, got %r" % (spec,))
    return list(available[:n])


def select_poses(poses, n_pse=None, n_svg=None):
    """Which poses get a session and which get an SVG.

    *poses* is a sequence of ``(pose_id, do_score)``. Returns
    ``(pse_ids, svg_ids)``. Ranking is by DO, descending, with the pose id as
    the tie-break so the order is stable between runs.
    """
    ranked = [p for p, _ in sorted(poses, key=lambda t: (-float(t[1]), str(t[0])))]
    return parse_selection(n_pse, ranked), parse_selection(n_svg, ranked)


def spread_in_view(anchors, x_axis, y_axis, *, min_sep=16.0,
                   iterations: int = 60, step: float = 0.5):
    """Push label anchors apart **in the camera plane** until they stop overlapping.

    Labels are placed in 3D but read in 2D, so two anchors that are far apart
    in space can still print on top of each other. Offsetting along a world
    axis does not help and offsetting along each patch's own surface normal
    only helps when the patches face different ways -- patches clustered on
    one face of an enzyme share a normal and stack anyway.

    So the separation is done where the collision happens: project onto the
    camera's two screen axes, relax the points apart, and add the resulting
    displacement back in the camera plane. Depth is never changed, so a label
    stays at its patch's distance from the camera and in front of the same
    surface.

    ``min_sep`` may be a single number or ``(along_x, along_y)``. Two numbers
    are usually what you want: a label is a line of text, far wider than it is
    tall, so two anchors at the same height need much more horizontal room
    than two stacked vertically. A circular criterion either leaves wide
    labels overlapping or flings short ones needlessly far apart.

    Returns an ``(N, 3)`` array of displacements to add to ``anchors``.
    """
    pts = np.atleast_2d(np.asarray(anchors, float))
    if pts.size == 0:
        return np.zeros((0, 3))
    if pts.shape[-1] != 3:
        raise ValueError("anchors must be (N, 3), got %r" % (pts.shape,))
    x = _unit(x_axis, "x_axis")
    y = _unit(y_axis, "y_axis")
    uv = np.stack([pts @ x, pts @ y], axis=1)
    n = len(uv)
    if n < 2:
        return np.zeros((n, 3))
    sep = np.asarray(min_sep, float).reshape(-1)
    if sep.size == 1:
        sep = np.array([float(sep[0]), float(sep[0])])
    elif sep.size != 2:
        raise ValueError("min_sep takes one number or two, got %r" % (min_sep,))
    if (sep <= 0).any():
        raise ValueError("min_sep must be positive, got %r" % (min_sep,))
    # work in units of the required separation, so the criterion is the unit
    # circle and the relaxation stays isotropic in that space
    uv = uv / sep
    moved = uv.copy()
    min_sep = 1.0
    for _ in range(int(iterations)):
        worst = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                d = moved[j] - moved[i]
                dist = float(np.linalg.norm(d))
                if dist >= min_sep:
                    continue
                if dist < 1e-6:
                    # exactly coincident: separate along an arbitrary but
                    # deterministic direction rather than dividing by zero.
                    # The distance stays 0 -- setting it to 1 here made the
                    # push (min_sep - dist) evaluate to nothing once min_sep
                    # was normalised to 1, so coincident labels never moved.
                    direction = np.array([1.0, 0.0]) if i % 2 == 0 else np.array([0.0, 1.0])
                    dist = 0.0
                else:
                    direction = d / dist
                push = direction * (min_sep - dist) * float(step)
                moved[i] -= push
                moved[j] += push
                worst = max(worst, min_sep - dist)
        if worst < 1e-3:
            break
    delta = (moved - uv) * sep
    return delta[:, 0, None] * x[None, :] + delta[:, 1, None] * y[None, :]
