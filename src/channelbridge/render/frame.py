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
from typing import List, Optional, Sequence, Tuple

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


def _put_upstream_left(x, z, p1, p2, upstream_left: bool):
    """Turn the frame 180 degrees about y unless the upstream BODY is left.

    House rule: the upstream (green) enzyme is always on the left. Judged on
    P1, not on the exits: P1 lies one arm inside the upstream channel, so it
    points at that enzyme's body, while the two exits sit together at the
    interface and say little about which protein is where.

    Checked rather than assumed, because each frame builds x from a different
    vector -- the O_IC arm points INTO the upstream protein, the D_IC line
    points AWAY from it -- so a fixed "apply the flip" means opposite things
    in the two, which is how ic_line_frame shipped mirrored.
    """
    if not upstream_left:
        return x, z
    p1 = np.asarray(p1, float).reshape(3)
    p2 = np.asarray(p2, float).reshape(3)
    if float((p1 - p2) @ x) > 0.0:          # the upstream body is to the right
        return -x, -z
    return x, z


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

    # 180 degrees about the vertical if needed: x -> -x, z -> -z, y unchanged.
    x, z = _put_upstream_left(x, z, p1, p2, upstream_left)

    c = np.asarray(p2 + p3, dtype=float) / 2.0 if centre is None else np.asarray(centre, float)
    return Frame(np.array([x, y, z]), c, name="rotated0")


def ic_line_frame(p1, p2, p3, *, roll: float = 0.0, centre=None,
                  upstream_left: bool = True) -> Frame:
    """A frame with the **D_IC line along x**, rolled about it by *roll*.

    ``angle_frame`` lays the O_IC ARM along x and holds the P1-P2-P3 plane
    face-on, so O_IC is drawn at its true value. That is the right frame for
    reading the angle, and it leaves no freedom at all: every other choice is
    fixed. When a channel's drawn volume happens to lie across the D_IC line
    in that view, there is nothing to be done about it.

    This frame gives up one thing to buy that freedom. The line joining the
    exits is the x axis, and the camera may roll about it by any angle; the
    line stays horizontal at every roll, while the tunnels swing out of the
    way. What is lost is the face-on plane, so O_IC is **foreshortened** --
    use :func:`drawn_oic` to find out by how much, and pass the projected
    arm to the SVG assembler as ``axis_px`` rather than letting it assume
    the true angle.

    How much foreshortening a given roll costs depends on the pose: near
    O_IC 180 the arm is almost antiparallel to the line and barely moves,
    near 90 it moves a lot. So the roll is chosen per pose against a stated
    error budget, never fixed as a house constant.
    """
    p1 = np.asarray(p1, float).reshape(3)
    p2 = np.asarray(p2, float).reshape(3)
    p3 = np.asarray(p3, float).reshape(3)
    x = _unit(p3 - p2, "P2->P3 (the D_IC line)")
    seed = p1 - p2
    seed = seed - float(seed @ x) * x
    if float(np.linalg.norm(seed)) < 1e-6:
        raise ValueError(
            "P1, P2 and P3 are collinear (O_IC is 0 or 180 degrees), so there "
            "is no in-plane direction to roll from. Use angle_frame()."
        )
    seed = _unit(seed, "in-plane perpendicular")
    th = np.radians(float(roll))
    y = _unit(np.cos(th) * seed + np.sin(th) * np.cross(x, seed), "camera y")
    z = np.cross(x, y)
    x, z = _put_upstream_left(x, z, p1, p2, upstream_left)
    c = (p2 + p3) / 2.0 if centre is None else np.asarray(centre, float)
    return Frame(np.array([x, y, z]), c, name="ic%+03d" % round(float(roll)))


def arm_frame(p1, p2, p3, *, roll: float = 0.0, centre=None,
              upstream_left: bool = True) -> Frame:
    """``angle_frame`` with the camera free to roll about the **O_IC arm**.

    x is P2->P1, the upstream channel's last two spheres, exactly as in
    :func:`angle_frame`; at ``roll=0`` this IS that frame. Rolling turns the
    camera about that arm, so the arm stays horizontal while the rest of the
    picture -- the D_IC line, and the drawn channels -- swings.

    The cost is larger than rolling about the D_IC line, and in a second
    place: the line leaves the view plane, so **both** the drawn O_IC and the
    drawn LENGTH of the D_IC line shrink. Report both (:func:`drawn_oic`,
    :func:`drawn_dic_fraction`) and budget both.
    """
    p1 = np.asarray(p1, float).reshape(3)
    p2 = np.asarray(p2, float).reshape(3)
    p3 = np.asarray(p3, float).reshape(3)
    x = _unit(p1 - p2, "P2->P1 (the O_IC arm)")
    seed = p3 - p2
    seed = seed - float(seed @ x) * x
    if float(np.linalg.norm(seed)) < 1e-6:
        raise ValueError(
            "P1, P2 and P3 are collinear (O_IC is 0 or 180 degrees), so there "
            "is no in-plane direction to roll from. Use angle_frame()."
        )
    seed = _unit(seed, "in-plane perpendicular")
    th = np.radians(float(roll))
    y = _unit(np.cos(th) * seed + np.sin(th) * np.cross(x, seed), "camera y")
    z = np.cross(x, y)
    x, z = _put_upstream_left(x, z, p1, p2, upstream_left)
    c = (p2 + p3) / 2.0 if centre is None else np.asarray(centre, float)
    return Frame(np.array([x, y, z]), c, name="arm%+03d" % round(float(roll)))


def drawn_dic_fraction(frame: Frame, p2, p3) -> float:
    """Drawn length of the D_IC line over its true length.

    1.0 when the line lies in the view plane. Anything less means the figure
    shows a shorter line than the number beside it: say so, or do not roll.
    """
    p2 = np.asarray(p2, float).reshape(3)
    p3 = np.asarray(p3, float).reshape(3)
    v = p3 - p2
    n = float(np.linalg.norm(v))
    if n < 1e-9:
        return float("nan")
    u = float(v @ frame.axes[0])
    w = float(v @ frame.axes[1])
    return float(np.hypot(u, w) / n)


def drawn_oic(frame: Frame, p1, p2, p3) -> float:
    """The O_IC a viewer would measure off this frame, in degrees.

    Equal to the true O_IC only when the P1-P2-P3 plane is face-on. Report it
    beside the true value whenever they differ: a figure whose angle cannot
    be measured off it, without saying so, is worse than one that does not
    draw the angle at all.
    """
    p1 = np.asarray(p1, float).reshape(3)
    p2 = np.asarray(p2, float).reshape(3)
    p3 = np.asarray(p3, float).reshape(3)
    x, y = frame.axes[0], frame.axes[1]
    a = np.array([float((p1 - p2) @ x), float((p1 - p2) @ y)])
    b = np.array([float((p3 - p2) @ x), float((p3 - p2) @ y)])
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9:
        return float("nan")
    return float(np.degrees(np.arccos(np.clip(a @ b / (na * nb), -1.0, 1.0))))


def clearest_roll(p1, p2, p3, points, *, about: str = "arm",
                  max_angle_error: float = 2.0, min_dic_fraction: float = 0.95,
                  step: float = 15.0, band: Tuple[float, float] = (0.12, 0.88),
                  centre=None, upstream_left: bool = True):
    """Pick the roll that keeps *points* furthest off the drawn D_IC line.

    ``about`` is the axis the camera rolls around:

    ``"arm"``
        the O_IC arm stays on x, as in :func:`angle_frame`, of which roll 0
        is exactly a copy. The D_IC line leaves the view plane, so both the
        drawn angle and the drawn length of the line shrink.
    ``"line"``
        the D_IC line stays on x. Its length is then always true and only
        the angle foreshortens.

    Clearance is the perpendicular distance in the IMAGE from the projected
    D_IC segment, not from a camera axis -- the two are the same thing only
    when the line is on x, which is why measuring against the axis gave a
    frame-dependent answer and flattered ``"arm"`` by a factor of five.

    Only the stretch of the segment strictly between the exits is judged
    (``band``): a channel is entitled to sit at its own mouth.

    ``max_angle_error`` and ``min_dic_fraction`` cap what the roll may cost
    the figure's own measurements. Rolls outside either budget are not
    considered, so a tidier picture can never be bought with a number the
    reader would misread.

    Returns ``(roll, frame, clearance, drawn, true, dic_fraction)``.
    """
    builders = {"arm": arm_frame, "line": ic_line_frame}
    if about not in builders:
        raise ValueError("about must be 'arm' or 'line', got %r" % (about,))
    build = builders[about]
    pts = np.atleast_2d(np.asarray(points, float))
    p2a = np.asarray(p2, float).reshape(3)
    p3a = np.asarray(p3, float).reshape(3)
    if float(np.linalg.norm(p3a - p2a)) < 1e-6:
        raise ValueError("the two exits coincide; there is no line to clear")
    true = drawn_oic(angle_frame(p1, p2, p3), p1, p2, p3)
    lo, hi = float(band[0]), float(band[1])
    best = None
    for i in range(round(360.0 / float(step))):
        roll = i * float(step)
        f = build(p1, p2, p3, roll=roll, centre=centre,
                  upstream_left=upstream_left)
        shown = drawn_oic(f, p1, p2, p3)
        frac = drawn_dic_fraction(f, p2a, p3a)
        if abs(shown - true) > float(max_angle_error):
            continue
        if frac < float(min_dic_fraction):
            continue
        a = np.array([(p2a - f.centre) @ f.axes[0], (p2a - f.centre) @ f.axes[1]])
        b = np.array([(p3a - f.centre) @ f.axes[0], (p3a - f.centre) @ f.axes[1]])
        d = b - a
        denom = float(d @ d)
        q = np.stack([(pts - f.centre) @ f.axes[0], (pts - f.centre) @ f.axes[1]], axis=1)
        t = ((q - a) @ d) / denom
        inside = (t > lo) & (t < hi)
        if inside.any():
            foot = a[None, :] + t[inside, None] * d[None, :]
            clear = float(np.linalg.norm(q[inside] - foot, axis=1).min())
        else:
            clear = float("inf")
        if best is None or clear > best[2]:
            best = (roll, f, clear, shown, true, frac)
    if best is None:
        raise ValueError(
            "no roll meets the budgets (O_IC within %.1f deg of %.1f, D_IC at "
            "least %.0f%% of true); widen them or use angle_frame()"
            % (max_angle_error, true, 100 * min_dic_fraction))
    return best


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


def pymol_view(frame: Frame, *, slab: float = 3000.0, distance: float = -150.0,
               orthoscopic: bool = False) -> List[float]:
    """The 18 numbers ``cmd.set_view`` wants.

    PyMOL stores the rotation **column-wise**, so the axes matrix is
    transposed before flattening. The caller normally follows this with
    ``cmd.zoom(sel, buffer=0, complete=0)``, which keeps the rotation and
    fits the selection, then ``cmd.clip('slab', 3000)``.

    The last number is the field of view, and its SIGN is the orthoscopic
    flag: negative means perspective. Perspective is right for a figure of
    one complex, where it reads as depth. It is wrong for a map whose job is
    to compare patches: it makes a near patch larger than a far one of the
    same size, and -- the reason this option exists -- it means moving a
    label along the view axis drifts it across the screen, so lifting a label
    clear of the surface also slides it away from the patch it names. Under
    an orthographic camera that lift costs nothing on screen.
    """
    v = list(frame.axes.T.flatten())
    v += [0.0, 0.0, float(distance)]          # camera position
    v += [float(c) for c in frame.centre]      # origin of rotation
    v += [1.0, float(slab), 20.0 if orthoscopic else -20.0]
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


def lift_to_front(anchors, toward_viewer, coords, margin: float = 2.0):
    """Pull label anchors in front of the structure, along the view axis.

    A label sits at its own anchor's depth. One naming a patch on the FAR side
    of an enzyme therefore renders behind the near surface and washes out to
    unreadable -- a collision in depth, which :func:`spread_in_view` cannot
    fix because it deliberately never changes depth.

    PyMOL's ``float_labels`` is the obvious answer and does not work here: it
    applies to the interactive viewport, not to ray-traced output, which is
    what a figure is. So the anchors are moved instead.

    Only anchors that are actually behind get moved, and only as far as they
    must be, because the camera is perspective: shifting a point along the
    view axis drifts its screen position outward, so an unnecessary lift
    pulls a label away from the patch it names.

    Returns an ``(N, 3)`` array of displacements to add to ``anchors``.
    """
    pts = np.atleast_2d(np.asarray(anchors, float))
    if pts.size == 0:
        return np.zeros((0, 3))
    if pts.shape[-1] != 3:
        raise ValueError("anchors must be (N, 3), got %r" % (pts.shape,))
    t = _unit(toward_viewer, "toward_viewer")
    body = np.atleast_2d(np.asarray(coords, float))
    if body.size == 0:
        return np.zeros((len(pts), 3))
    front = float((body @ t).max()) + float(margin)
    depth = pts @ t
    need = np.clip(front - depth, 0.0, None)
    return need[:, None] * t[None, :]


def spread_in_view(anchors, x_axis, y_axis, *, min_sep=16.0,
                   iterations: int = 60, step: float = 0.5, bounds=None):
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

    ``bounds`` is an optional ``(u_lo, u_hi, v_lo, v_hi)`` box in the same
    screen axes. Anchors are kept inside it, so a label pushed apart from its
    neighbours cannot be pushed off the edge of the frame instead -- which is
    what happened to the right-hand labels of the p6a hotspot map, where
    separation worked but one name still ran off the picture.

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
        if bounds is not None:
            lo = np.array([bounds[0], bounds[2]], float) / sep
            hi = np.array([bounds[1], bounds[3]], float) / sep
            moved = np.clip(moved, lo, hi)
        if worst < 1e-3:
            break
    delta = (moved - uv) * sep
    return delta[:, 0, None] * x[None, :] + delta[:, 1, None] * y[None, :]
