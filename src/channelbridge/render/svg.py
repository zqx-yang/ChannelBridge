"""The editable SVG: a vector overlay on a raster of the proteins.

Why an SVG at all
-----------------
The proteins have to be a raster -- there is no honest vector form of a
ray-traced molecular surface -- but everything the figure *claims* (where the
channel runs, how far apart the exits are, what angle they make, what the
numbers are) is a measurement, and a measurement in a figure should stay
editable geometry and editable text rather than becoming pixels. So the
proteins stay an embedded PNG and every line, arc, outline and label on top of
them is an SVG element with an id, which Illustrator and Inkscape can select,
recolour and move.

Why the overlay is traced in 2D rather than rendered in 3D
----------------------------------------------------------
PyMOL *can* draw the channel in 3D, and the ``render`` scenes do. Two things
go wrong when that 3D object is the final figure:

* ``ray_trace_mode 1`` (the outline mode the rest of the style depends on)
  outlines the tunnel's bumpy surface as a field of dots and creases, not as
  one clean edge;
* the part of a channel that reaches into solvent ray-traces as a red smear
  over nothing, because the tunnel is a cloud of spheres and solvent is empty.

Tracing the silhouette from a flat black helper render instead, and clipping
it to the enzymes' own footprint, gives exactly one closed curve per channel
and *guarantees* that no red is drawn outside a protein. The outline is then
cut open at the mouth, so the channel reads as something that opens onto the
surface rather than as a closed blob.

What this module needs
----------------------
The helper renders of the ``mask``, ``mouth`` and ``footprint`` scenes, plus
one transparent render of whichever ``render*`` scene is being assembled.
They all share a camera, so one set of helpers serves every render scene of
the same frame. SciPy, Matplotlib and Pillow are imported lazily: they are
``channelbridge[figures]`` extras, and the pose search must not need them.
"""

from __future__ import annotations

import base64
import io
import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .style import DEFAULT_STYLE, RenderStyle, parse_color, rgb_to_hex

__all__ = [
    "AngleMarker",
    "Mouth",
    "Placement",
    "SvgStyle",
    "angle_marker",
    "build_svg",
    "content_box",
    "find_mouths",
    "open_runs",
    "path_closed",
    "path_open",
    "place_block",
    "svg_document",
    "trace_contours",
    "tunnel_silhouette",
    "write_preview",
]

#: Every pixel constant below is quoted at this render width and scaled from
#: it, so the same numbers give the same-looking figure at preview size and at
#: print size.
REFERENCE_WIDTH = 2000.0


# ------------------------------------------------------------------ optional deps


def _imaging():
    """SciPy / Matplotlib / Pillow, or one error that says what to install."""
    try:
        import matplotlib
        from PIL import Image
        from scipy import ndimage as ndi

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
    except ImportError as exc:  # pragma: no cover - exercised by hand
        raise ImportError(
            "the SVG assembler needs scipy, matplotlib and pillow; "
            "install them with  pip install 'channelbridge[figures]'  "
            "(original error: %s)" % (exc,)
        ) from exc
    return ndi, plt, Image


# ----------------------------------------------------------------------- style


def _default_label_order() -> List[str]:
    return ["D_IC", "O_IC", "DO"]


@dataclass(frozen=True)
class SvgStyle:
    """Sizes and opacities of the vector layer. All of it overridable.

    Colours come from the :class:`~channelbridge.render.style.RenderStyle`
    (roles ``tunnel`` and ``annotation``) so that the SVG cannot drift away
    from the 3D scenes it is assembled from; everything that has no meaning in
    3D -- a stroke width in pixels, a font size, the radius of the angle arc --
    lives here.
    """

    #: Opacity of the channel fill. Not ``1 - transparency['tunnel']``: the 3D
    #: tunnel is seen through two translucent surfaces and needs to be faint,
    #: the flat 2D fill sits on top of everything and would disappear.
    fill_opacity: float = 0.30
    stroke_width: float = 5.0
    font_size: float = 46.0
    #: Radius of the O_IC arc, and the length of the drawn upstream-axis arm.
    arc_radius: float = 110.0
    arm_length: float = 150.0
    exit_dot_radius: float = 8.0
    #: How far a label block must stay clear of the channel line.
    label_clearance: float = 30.0
    #: Fraction of the content's larger side added as a margin around it.
    margin: float = 0.035
    #: How far past a mouth's radius the outline is erased, and how far behind
    #: the mouth the erasure still reaches (so the opening is not a slit).
    mouth_open_factor: float = 2.3
    mouth_open_back: float = 0.3
    #: Smoothing of the traced silhouette, in reference pixels.
    smooth: float = 1.5
    #: Smallest traced area kept, in reference pixels squared. Below this a
    #: contour is a speck of a neighbouring chain's surface, not a channel.
    min_area: float = 200.0
    font_family: str = "Helvetica, Arial, sans-serif"
    label_color: str = "black"
    #: Which measurements to print, in order. An empty list draws the geometry
    #: with no numbers, for a figure whose caption carries them.
    labels: List[str] = field(default_factory=_default_label_order)

    _ALLOWED = (
        "fill_opacity", "stroke_width", "font_size", "arc_radius", "arm_length",
        "exit_dot_radius", "label_clearance", "margin", "mouth_open_factor",
        "mouth_open_back", "smooth", "min_area", "font_family", "label_color",
        "label_colour", "labels",
    )
    _KNOWN_LABELS = ("D_IC", "O_IC", "DO")

    def __post_init__(self):
        for name in ("fill_opacity",):
            v = float(getattr(self, name))
            if not 0.0 <= v <= 1.0:
                raise ValueError("%s must be in [0, 1], got %r" % (name, v))
        for name in ("stroke_width", "font_size", "arc_radius", "arm_length",
                     "exit_dot_radius", "smooth", "mouth_open_factor"):
            if float(getattr(self, name)) <= 0:
                raise ValueError("%s must be positive, got %r" % (name, getattr(self, name)))
        for name in ("label_clearance", "margin", "min_area", "mouth_open_back"):
            if float(getattr(self, name)) < 0:
                raise ValueError("%s must not be negative, got %r" % (name, getattr(self, name)))
        parse_color(self.label_color)
        unknown = [str(k) for k in self.labels if str(k) not in self._KNOWN_LABELS]
        if unknown:
            raise ValueError(
                "unknown label(s) %s; known: %s"
                % (", ".join(repr(u) for u in unknown), ", ".join(self._KNOWN_LABELS))
            )

    @classmethod
    def from_dict(cls, raw: Optional[Mapping], where: str = "render.svg") -> "SvgStyle":
        raw = dict(raw or {})
        unknown = sorted(set(raw) - set(cls._ALLOWED))
        if unknown:
            raise ValueError(
                "%s: unknown key(s) %s; allowed: %s"
                % (where, ", ".join(repr(u) for u in unknown), ", ".join(sorted(cls._ALLOWED)))
            )
        if "label_color" in raw and "label_colour" in raw:
            raise ValueError("%s: give only one of label_color, label_colour" % (where,))
        kw = {k: raw[k] for k in raw if k in cls.__dataclass_fields__}
        if "label_colour" in raw:
            kw["label_color"] = raw["label_colour"]
        if "labels" in raw:
            labels = raw["labels"]
            if isinstance(labels, str):
                labels = [p.strip() for p in labels.replace(",", " ").split() if p.strip()]
            kw["labels"] = [str(x) for x in labels]
        for k in ("fill_opacity", "stroke_width", "font_size", "arc_radius", "arm_length",
                  "exit_dot_radius", "label_clearance", "margin", "mouth_open_factor",
                  "mouth_open_back", "smooth", "min_area"):
            if k in kw:
                kw[k] = float(kw[k])
        try:
            return cls(**kw)
        except TypeError as exc:
            raise ValueError("%s: %s" % (where, exc)) from None


# ------------------------------------------------------------------- raster bits


def _gray(path) -> np.ndarray:
    """A helper render as a 0-1 greyscale array, 0 black."""
    _, _, Image = _imaging()
    if not os.path.exists(path):
        raise FileNotFoundError(
            "missing helper render %s; render the mask/mouth/footprint scenes "
            "of this frame first" % (path,)
        )
    return np.asarray(Image.open(path).convert("L"), dtype=float) / 255.0


def tunnel_silhouette(mask, footprint, scale: float = 1.0, occlusion=None,
                      smooth: float = 2.0, dilate: int = 2) -> np.ndarray:
    """Where the channel may be painted: the mask, clipped to the enzymes.

    ``mask`` and ``footprint`` are flat renders (black on white) as 0-1
    arrays. The clip is the whole point: an unclipped silhouette of a 3D
    tunnel leaves red over solvent, because the sphere cloud reaches past the
    surface. The footprint is dilated by a couple of pixels first so that the
    channel is not shaved along the enzyme's own anti-aliased edge.

    ``occlusion`` is a flat render in which any *opaque* chain is black and
    the channel objects grey. Where such a chain is in front, the channel must
    not be drawn, or the overlay floats over something solid.
    """
    ndi, _, _ = _imaging()
    mask = np.asarray(mask, float)
    footprint = np.asarray(footprint, float)
    if mask.shape != footprint.shape:
        raise ValueError(
            "mask %r and footprint %r come from different renders; they must "
            "share a camera and a size" % (mask.shape, footprint.shape)
        )
    inside = ndi.binary_dilation(footprint < 0.5, iterations=max(1, round(float(dilate) * scale)))
    soft = ndi.gaussian_filter((mask < 0.5).astype(float), smooth * scale) > 0.5
    sil = soft & inside
    if occlusion is not None:
        occ = np.asarray(occlusion, float)
        if occ.shape != mask.shape:
            raise ValueError("occlusion render has size %r, mask %r" % (occ.shape, mask.shape))
        sil &= ~ndi.binary_dilation(occ < 0.3, iterations=1)
    return sil


@dataclass(frozen=True)
class Mouth:
    """One channel exit as it appears in the picture."""

    key: str                      #: ``"up"`` or ``"dn"``
    centre: np.ndarray            #: (x, y) in pixels
    radius: float                 #: of the disc, in pixels
    forward: np.ndarray           #: unit vector out of the protein, in pixels

    def __str__(self):
        return "%s at (%.0f, %.0f) r=%.0f" % (
            self.key, self.centre[0], self.centre[1], self.radius)


def find_mouths(mouth, silhouette, min_pixels: int = 20,
                dark_below: float = 0.25) -> Dict[str, Mouth]:
    """The two exit discs, told apart by their grey level.

    The ``mouth`` helper scene draws the upstream exit black and the
    downstream exit grey precisely so that a flat render identifies which is
    which without a second file. The outward direction is estimated from the
    channel itself: the silhouette near the disc lies mostly *behind* the
    exit, so the opposite of that offset points out of the protein.
    """
    ndi, _, _ = _imaging()
    mouth = np.asarray(mouth, float)
    sil = np.asarray(silhouette, bool)
    if mouth.shape != sil.shape:
        raise ValueError("mouth render %r and silhouette %r differ" % (mouth.shape, sil.shape))
    h, w = mouth.shape
    yy, xx = np.mgrid[0:h, 0:w]
    labels, n = ndi.label(mouth < 0.85)
    found: Dict[str, List[Mouth]] = {}
    for i in range(1, n + 1):
        blob = labels == i
        if blob.sum() < min_pixels:
            continue
        centre = np.array([xx[blob].mean(), yy[blob].mean()])
        radius = float(np.sqrt(blob.sum() / np.pi))
        key = "up" if mouth[blob].mean() < dark_below else "dn"
        near = ((xx - centre[0]) ** 2 + (yy - centre[1]) ** 2 <= (3 * radius) ** 2) & sil
        if near.sum() < 3:
            # no channel beside this disc to take a direction from: fall back
            # to pointing away from the picture's centre, which is where the
            # other enzyme is, rather than inventing a direction of zero
            back = np.array([w / 2.0, h / 2.0]) - centre
        else:
            back = np.array([xx[near].mean(), yy[near].mean()]) - centre
        norm = float(np.linalg.norm(back))
        forward = -back / norm if norm > 1e-9 else np.array([1.0, 0.0])
        found.setdefault(key, []).append(Mouth(key, centre, radius, forward))
    out = {}
    for key, blobs in found.items():
        blobs.sort(key=lambda m: -m.radius)
        out[key] = blobs[0]
    missing = [k for k in ("up", "dn") if k not in out]
    if missing:
        raise ValueError(
            "the mouth render shows no %s exit (found %d disc(s)); it should "
            "carry one black disc (upstream) and one grey disc (downstream)"
            % (" and no ".join(missing), n)
        )
    return out


def trace_contours(silhouette, scale: float = 1.0, smooth: float = 1.5,
                   min_area: float = 200.0, min_points: int = 20) -> List[np.ndarray]:
    """One closed polygon per channel, in image coordinates (y down)."""
    ndi, plt, _ = _imaging()
    sil = np.asarray(silhouette, bool)
    h, w = sil.shape
    field_ = ndi.gaussian_filter(sil.astype(float), smooth * scale)
    cs = plt.contour(np.arange(w), np.arange(h), field_, [0.5])
    try:
        segs = cs.allsegs[0]
    finally:
        plt.close("all")
    polys = [np.asarray(p, float) for p in segs if len(p) > min_points]
    keep = [p for p in polys if _area(p) > min_area * scale * scale]
    return keep


def _area(poly) -> float:
    p = np.asarray(poly, float)
    return 0.5 * abs(
        float(np.dot(p[:, 0], np.roll(p[:, 1], -1)) - np.dot(p[:, 1], np.roll(p[:, 0], -1)))
    )


def open_runs(poly, mouths: Sequence[Mouth], factor: float = 2.3,
              back: float = 0.3, min_run: int = 5) -> List[np.ndarray]:
    """The contour cut open where it crosses a mouth.

    A channel that opens onto the surface should be drawn as a curve with two
    free ends at the opening, not as a closed blob: the blob reads as a cavity.
    Points within ``factor`` radii of a mouth centre *and* on its outward side
    are dropped, and what remains is returned as polylines.

    The array is rotated so that the first kept point follows a dropped one,
    otherwise a run that straddles the start of the contour would come back
    as two pieces with a visible gap at an arbitrary place.
    """
    p = np.asarray(poly, float)
    bad = np.zeros(len(p), bool)
    for m in mouths:
        d = p - m.centre
        bad |= ((d ** 2).sum(1) <= (factor * m.radius) ** 2) & (d @ m.forward >= -back * m.radius)
    if not bad.any():
        return [np.vstack([p, p[:1]])]
    if bad.all():
        return []
    i0 = int(np.argmax(bad))
    p = np.roll(p, -i0, axis=0)
    bad = np.roll(bad, -i0)
    runs, cur = [], []
    for q, b in zip(p, bad):
        if b:
            if len(cur) > min_run:
                runs.append(np.asarray(cur))
            cur = []
        else:
            cur.append(q)
    if len(cur) > min_run:
        runs.append(np.asarray(cur))
    return runs


def path_closed(poly, step: int = 2) -> str:
    p = np.asarray(poly, float)[::step]
    return "M" + " L".join("%.1f,%.1f" % (q[0], q[1]) for q in p) + " Z"


def path_open(poly, step: int = 2) -> str:
    p = np.asarray(poly, float)[::step]
    return "M" + " L".join("%.1f,%.1f" % (q[0], q[1]) for q in p)


# -------------------------------------------------------------- label placement


@dataclass(frozen=True)
class Placement:
    """Where a block of text goes, and how much protein it covers."""

    x: float
    y: float
    overlap: float


def place_block(alpha, p_up, p_dn, width: float, height: float, above: bool,
                clearance: float = 30.0, scale: float = 1.0,
                step: float = 25.0, reach: float = 0.4) -> Placement:
    """Put a text block where it hides the least protein, on one side of the line.

    "Above" and "below" are taken perpendicular to the channel line, not in
    screen y: a near-vertical line has no above, and the p5c script's
    ``slope``-based test degenerates there. Using the line's own normal --
    oriented so that it points up on the screen -- keeps the two blocks on
    opposite sides whatever the angle, which is what the figure needs.

    Every corner of the block must clear the line by ``clearance``, so a block
    never leans over the measurement it annotates.
    """
    a = np.asarray(alpha, bool)
    h, w = a.shape
    p_up = np.asarray(p_up, float)
    p_dn = np.asarray(p_dn, float)
    axis = p_dn - p_up
    n = float(np.linalg.norm(axis))
    if n < 1e-9:
        raise ValueError("the two exits project to the same point; no line to place against")
    axis = axis / n
    normal = np.array([axis[1], -axis[0]])
    if normal[1] > 0:                        # y grows downwards: up on screen is -y
        normal = -normal
    mid = (p_up + p_dn) / 2.0
    want = 1.0 if above else -1.0
    span = reach * max(w, h)
    grid = np.arange(-span, span + 1, max(1.0, step * scale))
    best = None
    for dy in grid:
        for dx in grid:
            x0 = mid[0] + dx - width / 2.0
            y0 = mid[1] + dy - height / 2.0
            if x0 < 20 * scale or y0 < 20 * scale:
                continue
            if x0 + width > w - 20 * scale or y0 + height > h - 20 * scale:
                continue
            corners = np.array([[x0, y0], [x0 + width, y0],
                                [x0, y0 + height], [x0 + width, y0 + height]])
            signed = (corners - p_up) @ normal
            if (want * signed).min() < clearance * scale:
                continue
            patch = a[int(y0):int(y0 + height), int(x0):int(x0 + width)]
            overlap = float(patch.mean()) if patch.size else 1.0
            cost = overlap * 20.0 + math.hypot(dx, dy) / max(span, 1.0) * 2.0
            if best is None or cost < best[0]:
                best = (cost, x0, y0, overlap)
    if best is None:
        # nowhere clear: put it at the midpoint, offset to the right side, and
        # let the caller see a 100 % overlap rather than get no label at all
        off = mid + want * normal * (clearance * scale + height)
        return Placement(float(off[0] - width / 2), float(off[1] - height / 2), 1.0)
    return Placement(float(best[1]), float(best[2]), float(best[3]))


# --------------------------------------------------------------- angle geometry


@dataclass(frozen=True)
class AngleMarker:
    """The O_IC marker: the upstream axis arm, its arrowhead, and the arc."""

    arm_start: np.ndarray
    vertex: np.ndarray
    arrowhead: np.ndarray        #: three points
    arc_start: np.ndarray
    arc_end: np.ndarray
    arc_mid: np.ndarray
    radius: float
    exact: bool                  #: False when the arm came from the analytic rule


def angle_marker(p_up, p_dn, oic_deg: float, radius: float, arm: float,
                 scale: float = 1.0, axis_px=None) -> AngleMarker:
    """The angle marker at the upstream exit, after iMARS Fig. 1D.

    O_IC is the angle P1-P2-P3: P2 is the upstream exit, P3 the downstream
    exit, P1 a point back along the upstream channel. In the ``rotated0``
    frame the plane P1-P2-P3 is face-on, so the ray P2->P1 is the ray P2->P3
    turned by exactly ``-oic_deg`` in image angles (y down) and the arc spans
    the true angle. That is what is drawn when ``axis_px`` is not given, and
    the marker is flagged ``exact=False`` to say the arm was *derived* from
    the number rather than measured off the render.

    Pass ``axis_px`` -- the projected P2->P1 direction in pixels -- for any
    other camera, where that rule does not hold.
    """
    p_up = np.asarray(p_up, float)
    p_dn = np.asarray(p_dn, float)
    radius = float(radius) * scale
    arm = float(arm) * scale
    th3 = math.atan2(p_dn[1] - p_up[1], p_dn[0] - p_up[0])
    exact = axis_px is not None
    if exact:
        v = np.asarray(axis_px, float).reshape(2)
        if float(np.linalg.norm(v)) < 1e-9:
            raise ValueError("axis_px has zero length")
        th1 = math.atan2(v[1], v[0])
    else:
        th1 = th3 - math.radians(float(oic_deg))
    u1 = np.array([math.cos(th1), math.sin(th1)])
    u3 = np.array([math.cos(th3), math.sin(th3)])
    perp = np.array([-u1[1], u1[0]])
    head = np.array([
        p_up,
        p_up + 30 * scale * u1 + 11 * scale * perp,
        p_up + 30 * scale * u1 - 11 * scale * perp,
    ])
    th_mid = th3 + _signed_delta(th3, th1) / 2.0
    return AngleMarker(
        arm_start=p_up + arm * u1,
        vertex=p_up,
        arrowhead=head,
        arc_start=p_up + radius * u3,
        arc_end=p_up + radius * u1,
        arc_mid=p_up + radius * np.array([math.cos(th_mid), math.sin(th_mid)]),
        radius=radius,
        exact=exact,
    )


def _signed_delta(a: float, b: float) -> float:
    """``b - a`` wrapped into (-pi, pi]: the short way round."""
    d = (b - a) % (2 * math.pi)
    return d - 2 * math.pi if d > math.pi else d


def _arc_flags(marker: AngleMarker) -> Tuple[int, int]:
    """SVG ``large-arc-flag`` and ``sweep-flag`` for the marker's arc."""
    a = math.atan2(marker.arc_start[1] - marker.vertex[1],
                   marker.arc_start[0] - marker.vertex[0])
    b = math.atan2(marker.arc_end[1] - marker.vertex[1],
                   marker.arc_end[0] - marker.vertex[0])
    d = _signed_delta(a, b)
    return (1 if abs(d) > math.pi else 0), (1 if d > 0 else 0)


# ------------------------------------------------------------------ framing


def content_box(alpha, silhouette, points: Sequence[Sequence[float]] = (),
                margin: float = 0.035) -> Tuple[float, float, float, float]:
    """The viewBox: everything drawn, plus a margin.

    A render is mostly empty -- the camera was fitted to the complex, not to
    the page -- so cropping to the content is what makes the complex fill the
    figure. Cropping is done in the SVG's viewBox rather than on the raster so
    that the embedded image keeps its own pixel coordinates and the overlay's
    numbers stay the render's numbers.
    """
    a = np.asarray(alpha, bool)
    s = np.asarray(silhouette, bool)
    h, w = a.shape
    ys, xs = np.where(a | s)
    if len(xs) == 0:
        x0, x1, y0, y1 = 0.0, float(w), 0.0, float(h)
    else:
        x0, x1, y0, y1 = float(xs.min()), float(xs.max()), float(ys.min()), float(ys.max())
    for pt in points:
        x0, x1 = min(x0, float(pt[0])), max(x1, float(pt[0]))
        y0, y1 = min(y0, float(pt[1])), max(y1, float(pt[1]))
    m = float(margin) * max(x1 - x0, y1 - y0)
    vx, vy = max(0.0, x0 - m), max(0.0, y0 - m)
    return vx, vy, min(float(w), x1 + m) - vx, min(float(h), y1 + m) - vy


# ------------------------------------------------------------------- the document


def _subscript(head: str, sub: str, tail: str, fs: float) -> str:
    return (
        '<tspan>%s</tspan><tspan dy="%.1f" font-size="%.1f">%s</tspan>'
        '<tspan dy="%.1f">%s</tspan>'
        % (head, 0.22 * fs, 0.68 * fs, sub, -0.22 * fs, _escape(tail))
    )


def _escape(text) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _embed(path) -> str:
    _, _, Image = _imaging()
    img = Image.open(path).convert("RGBA")
    buf = io.BytesIO()
    img.save(buf, "PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii")


@dataclass
class _Layout:
    """Everything positioned, before any text is written."""

    width: int
    height: int
    scale: float
    mouths: Dict[str, Mouth]
    polys: List[np.ndarray]
    silhouette: np.ndarray
    alpha: np.ndarray
    marker: Optional[AngleMarker]
    blocks: Dict[str, Placement]
    box: Tuple[float, float, float, float]


def svg_document(layout: _Layout, raster_b64: str, values: Mapping,
                 style: RenderStyle, svg_style: SvgStyle, title: str = "") -> str:
    """The SVG text. Pure: everything it needs is already measured."""
    k = layout.scale
    tunnel = rgb_to_hex(style.palette.rgb("tunnel"))
    annot = rgb_to_hex(style.palette.rgb("annotation"))
    label_rgb = rgb_to_hex(parse_color(svg_style.label_color))
    fs = svg_style.font_size * k
    sw = svg_style.stroke_width * k
    vx, vy, vw, vh = layout.box
    up, dn = layout.mouths["up"], layout.mouths["dn"]
    out: List[str] = []
    out.append(
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        'width="%.0f" height="%.0f" viewBox="%.1f %.1f %.1f %.1f">' % (vw, vh, vx, vy, vw, vh)
    )
    out.append("<title>%s</title>" % _escape(title or "channel geometry"))
    out.append(
        '<rect id="background" x="%.1f" y="%.1f" width="%.1f" height="%.1f" '
        'fill="#FFFFFF" fill-opacity="0"/>' % (vx, vy, vw, vh)
    )
    out.append(
        '<image id="proteins" x="0" y="0" width="%d" height="%d" '
        'xlink:href="data:image/png;base64,%s"/>' % (layout.width, layout.height, raster_b64)
    )
    out.append('<g id="tunnels">')
    for i, p in enumerate(layout.polys):
        out.append(
            '<path id="tunnel_fill_%d" d="%s" fill="%s" fill-opacity="%.2f" stroke="none"/>'
            % (i, path_closed(p), tunnel, svg_style.fill_opacity)
        )
    for i, p in enumerate(layout.polys):
        for j, run in enumerate(open_runs(p, list(layout.mouths.values()),
                                          svg_style.mouth_open_factor,
                                          svg_style.mouth_open_back)):
            out.append(
                '<path id="tunnel_outline_%d_%d" d="%s" fill="none" stroke="%s" '
                'stroke-width="%.1f" stroke-linecap="round" stroke-linejoin="round"/>'
                % (i, j, path_open(run), tunnel, sw)
            )
    out.append("</g>")
    out.append(
        '<g id="channel">'
        '<line id="channel_line" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
        'stroke-width="%.1f" stroke-linecap="round"/>'
        '<circle id="exit_upstream" cx="%.1f" cy="%.1f" r="%.1f" fill="%s"/>'
        '<circle id="exit_downstream" cx="%.1f" cy="%.1f" r="%.1f" fill="%s"/></g>'
        % (up.centre[0], up.centre[1], dn.centre[0], dn.centre[1], annot, sw,
           up.centre[0], up.centre[1], svg_style.exit_dot_radius * k, annot,
           dn.centre[0], dn.centre[1], svg_style.exit_dot_radius * k, annot)
    )
    m = layout.marker
    if m is not None:
        large, sweep = _arc_flags(m)
        out.append(
            '<g id="angle" fill="none" stroke="%s" stroke-width="%.1f" stroke-linecap="round" '
            'stroke-linejoin="round">'
            '<line id="segment_X" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
            '<polygon id="segment_X_arrowhead" points="%s" fill="%s" stroke="none"/>'
            '<path id="angle_arc" d="M%.1f,%.1f A%.1f,%.1f 0 %d %d %.1f,%.1f"/></g>'
            % (annot, sw, m.arm_start[0], m.arm_start[1], m.vertex[0], m.vertex[1],
               " ".join("%.1f,%.1f" % (p[0], p[1]) for p in m.arrowhead), annot,
               m.arc_start[0], m.arc_start[1], m.radius, m.radius, large, sweep,
               m.arc_end[0], m.arc_end[1])
        )
        # the two segment names, beside the things they name
        normal = _line_normal(up.centre, dn.centre)
        lx = m.vertex + np.array([-26 * k, 0.0]) - np.array([0.0, m.radius + 22 * k])
        ly = (up.centre + dn.centre) / 2.0 + normal * 30 * k
        out.append(
            '<g id="segment_names" font-family="%s" font-size="%.1f" font-style="italic" '
            'fill="%s">'
            '<text id="label_X" x="%.1f" y="%.1f"><tspan>O</tspan>'
            '<tspan dy="%.1f" font-size="%.1f">IC</tspan></text>'
            '<text id="label_Y" x="%.1f" y="%.1f"><tspan>D</tspan>'
            '<tspan dy="%.1f" font-size="%.1f">IC</tspan></text></g>'
            % (svg_style.font_family, 0.85 * fs, annot,
               lx[0], lx[1], 0.2 * fs, 0.58 * fs,
               ly[0], ly[1], 0.2 * fs, 0.58 * fs)
        )
    if layout.blocks:
        out.append(
            '<g id="labels" font-family="%s" font-size="%.1f" fill="%s">'
            % (svg_style.font_family, fs, label_rgb)
        )
        above = layout.blocks.get("above")
        if above is not None:
            line = 0
            for key in svg_style.labels:
                if key == "DO":
                    continue
                line += 1
                y = above.y + fs + (line - 1) * 1.1 * fs
                if key == "D_IC":
                    body = _subscript("D", "IC", " = %.1f Å" % float(values["D_IC"]), fs)
                else:
                    body = _subscript("O", "IC", " = %.1f°" % float(values["O_IC"]), fs)
                out.append('<text id="label_%s" x="%.1f" y="%.1f">%s</text>'
                           % (key.replace("_", ""), above.x, y, body))
        below = layout.blocks.get("below")
        if below is not None and "DO" in svg_style.labels:
            out.append('<text id="label_DO" x="%.1f" y="%.1f">DO = %.3f</text>'
                       % (below.x, below.y + fs, float(values["DO"])))
        out.append("</g>")
    out.append("</svg>")
    return "\n".join(out) + "\n"


def _line_normal(a, b) -> np.ndarray:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = b - a
    n = float(np.linalg.norm(d))
    if n < 1e-9:
        return np.array([0.0, -1.0])
    d = d / n
    normal = np.array([d[1], -d[0]])
    return normal if normal[1] < 0 else -normal


# -------------------------------------------------------------------- the driver


def build_svg(base_png, mask_png, mouth_png, footprint_png, out_svg, values: Mapping,
              style: Optional[RenderStyle] = None, svg_style: Optional[SvgStyle] = None,
              occlusion_png=None, axis_px=None, title: str = "",
              preview_png=None, report=print) -> str:
    """Assemble one editable SVG and write it.

    ``values`` carries ``D_IC``, ``O_IC`` and ``DO``. The helper renders must
    all come from the same frame as ``base_png``; they are reused across every
    ``render*`` scene of that frame, because the camera and the channels do
    not change between them.
    """
    style = style or DEFAULT_STYLE
    svg_style = svg_style or SvgStyle.from_dict(getattr(style, "svg", None))
    _, _, Image = _imaging()
    rgba = Image.open(base_png).convert("RGBA")
    w, h = rgba.size
    scale = w / REFERENCE_WIDTH
    mask = _gray(mask_png)
    if mask.shape != (h, w):
        raise ValueError(
            "the helper renders are %d x %d but %s is %d x %d; render them at "
            "the same size" % (mask.shape[1], mask.shape[0], os.path.basename(base_png), w, h)
        )
    occ = _gray(occlusion_png) if occlusion_png else None
    sil = tunnel_silhouette(mask, _gray(footprint_png), scale, occ)
    mouths = find_mouths(_gray(mouth_png), sil)
    report("exits: " + "; ".join(str(mouths[k]) for k in ("up", "dn")))
    polys = trace_contours(sil, scale, svg_style.smooth, svg_style.min_area)
    report("traced %d channel outline(s) over %d silhouette pixel(s)"
           % (len(polys), int(sil.sum())))
    if not polys:
        raise ValueError(
            "no channel outline survived tracing; the mask render may be empty "
            "or the footprint may not overlap it"
        )
    alpha = np.asarray(rgba)[..., 3] > 20
    marker = None
    want_marker = ("O_IC" in values and "O_IC" in svg_style.labels) or axis_px is not None
    if want_marker:
        marker = angle_marker(mouths["up"].centre, mouths["dn"].centre,
                              float(values.get("O_IC", 0.0)), svg_style.arc_radius,
                              svg_style.arm_length, scale, axis_px)
        if not marker.exact:
            report("O_IC arm drawn from the reported angle (rotated0 rule), not measured "
                   "off the render")
    fs = svg_style.font_size * scale
    blocks: Dict[str, Placement] = {}
    n_above = len([k for k in svg_style.labels if k in ("D_IC", "O_IC")])
    if n_above:
        blocks["above"] = place_block(
            alpha, mouths["up"].centre, mouths["dn"].centre, 380 * scale,
            (1.2 + 1.1 * n_above) * fs, True, svg_style.label_clearance, scale)
        report("label block (%d line(s), above the line) at (%.0f, %.0f), protein overlap %.0f %%"
               % (n_above, blocks["above"].x, blocks["above"].y, 100 * blocks["above"].overlap))
    if "DO" in svg_style.labels:
        blocks["below"] = place_block(
            alpha, mouths["up"].centre, mouths["dn"].centre, 300 * scale, 1.3 * fs,
            False, svg_style.label_clearance, scale)
        report("DO label (below the line) at (%.0f, %.0f), protein overlap %.0f %%"
               % (blocks["below"].x, blocks["below"].y, 100 * blocks["below"].overlap))
    pts: List[Sequence[float]] = [mouths["up"].centre, mouths["dn"].centre]
    if marker is not None:
        pts += [marker.arm_start, marker.arc_start, marker.arc_end, marker.arc_mid]
    for name, p in blocks.items():
        wblock = 380 * scale if name == "above" else 300 * scale
        hblock = (1.2 + 1.1 * n_above) * fs if name == "above" else 1.3 * fs
        pts += [(p.x, p.y), (p.x + wblock, p.y + hblock)]
    box = content_box(alpha, sil, pts, svg_style.margin)
    report("frame (viewBox): x %.0f y %.0f w %.0f h %.0f of %d x %d" % (*box, w, h))
    layout = _Layout(w, h, scale, mouths, polys, sil, alpha, marker, blocks, box)
    text = svg_document(layout, _embed(base_png), values, style, svg_style, title)
    with open(out_svg, "w", encoding="utf-8") as fh:
        fh.write(text)
    report("wrote %s (%.2f MB)" % (out_svg, len(text) / 1e6))
    if preview_png:
        write_preview(preview_png, base_png, layout, values, style, svg_style)
        report("preview %s" % preview_png)
    return text


def write_preview(out_png, base_png, layout: _Layout, values: Mapping,
                  style: RenderStyle, svg_style: SvgStyle,
                  background=(255, 255, 255)) -> str:
    """A raster of the finished SVG, drawn with Pillow.

    Not a general SVG renderer -- it draws the same layers from the same
    layout, which is the point. It answers the question a finished SVG cannot
    be checked for on a machine with no rasteriser: did the outline open at
    the right place, is the arc on the ligand-free side, did a label land on
    a protein. Cropped to the same viewBox, so what it shows is what the SVG
    shows.
    """
    _, _, Image = _imaging()
    from PIL import ImageDraw, ImageFont

    k = layout.scale
    base = Image.open(base_png).convert("RGBA")
    canvas = Image.new("RGBA", base.size, (*background, 255))
    canvas.alpha_composite(base)
    tunnel = tuple(round(c * 255) for c in style.palette.rgb("tunnel"))
    annot = tuple(round(c * 255) for c in style.palette.rgb("annotation"))
    sw = max(1, round(svg_style.stroke_width * k))

    fill_layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    fd = ImageDraw.Draw(fill_layer)
    for poly in layout.polys:
        fd.polygon([tuple(float(x) for x in q) for q in poly],
                   fill=(*tunnel, round(255 * svg_style.fill_opacity)))
    canvas.alpha_composite(fill_layer)

    d = ImageDraw.Draw(canvas)
    mouths = list(layout.mouths.values())
    for poly in layout.polys:
        for run in open_runs(poly, mouths, svg_style.mouth_open_factor,
                             svg_style.mouth_open_back):
            d.line([tuple(float(x) for x in q) for q in run], fill=(*annot, 255), width=sw)
    up, dn = layout.mouths["up"], layout.mouths["dn"]
    d.line([tuple(up.centre), tuple(dn.centre)], fill=(*annot, 255), width=sw)
    for m in (up, dn):
        r = svg_style.exit_dot_radius * k
        d.ellipse([m.centre[0] - r, m.centre[1] - r, m.centre[0] + r, m.centre[1] + r],
                  fill=(*annot, 255))
    mk = layout.marker
    if mk is not None:
        d.line([tuple(mk.arm_start), tuple(mk.vertex)], fill=(*annot, 255), width=sw)
        d.polygon([tuple(p) for p in mk.arrowhead], fill=(*annot, 255))
        a = math.degrees(math.atan2(mk.arc_start[1] - mk.vertex[1],
                                    mk.arc_start[0] - mk.vertex[0]))
        b = math.degrees(math.atan2(mk.arc_end[1] - mk.vertex[1],
                                    mk.arc_end[0] - mk.vertex[0]))
        if (b - a) % 360 > 180:
            a, b = b, a
        d.arc([mk.vertex[0] - mk.radius, mk.vertex[1] - mk.radius,
               mk.vertex[0] + mk.radius, mk.vertex[1] + mk.radius],
              start=a, end=b, fill=(*annot, 255), width=sw)

    fs = max(8, round(svg_style.font_size * k))
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", fs)
    except OSError:
        font = ImageFont.load_default()
    label_rgb = tuple(round(c * 255) for c in parse_color(svg_style.label_color))
    above = layout.blocks.get("above")
    if above is not None:
        line = 0
        for key in svg_style.labels:
            if key == "DO" or key not in values:
                continue
            txt = ("D_IC = %.1f Å" % float(values["D_IC"]) if key == "D_IC"
                   else "O_IC = %.1f°" % float(values["O_IC"]))
            d.text((above.x, above.y + line * 1.1 * fs), txt, fill=(*label_rgb, 255), font=font)
            line += 1
    below = layout.blocks.get("below")
    if below is not None and "DO" in values:
        d.text((below.x, below.y), "DO = %.3f" % float(values["DO"]),
               fill=(*label_rgb, 255), font=font)

    vx, vy, vw, vh = layout.box
    canvas = canvas.crop((int(vx), int(vy), int(vx + vw), int(vy + vh)))
    canvas.convert("RGB").save(out_png)
    return out_png
