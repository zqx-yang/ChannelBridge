"""The iMARS annotation: the D_IC line, the O_IC arm, the arc and the labels.

Its own module because TWO sessions draw it -- the pose session and the
hotspot map -- and a second copy would drift. It did: the map had a
hand-rolled version with different object names, so the same geometry read
differently depending on which file you opened.

Degenerate cases are handled rather than avoided. At O_IC 180 the three
points are collinear and there is no arc to draw; the function omits it
instead of drawing a bend that is not there. The hotspot map's reference
placement is always at O_IC 180, so that is the normal case there, not an
edge case.
"""

from __future__ import annotations

import numpy as np

from .style import RenderStyle

__all__ = ["build_annotation"]


def build_annotation(cmd, geom, st=None, label="", frame=None):
    """Dashed D_IC line, the O_IC arm arrow, the angle arc and the labels.

    ``geom`` needs ``p1``, ``p2``, ``p3``, ``d_ic``, ``o_ic`` and ``do``.
    Returns the object names it made, grouped under ``IC_annotation``.
    """
    from pymol.cgo import CONE, CYLINDER

    st = st or RenderStyle()
    rgb = [float(c) for c in st.palette.rgb("annotation")]
    p1, p2, p3 = geom["p1"], geom["p2"], geom["p3"]
    made = []

    cmd.delete("IC_mouths")
    for p in (p2, p3):
        cmd.pseudoatom("IC_mouths", pos=[float(p[0]), float(p[1]), float(p[2])], vdw=0.7)
    cmd.hide("everything", "IC_mouths")
    cmd.show("spheres", "IC_mouths")
    cmd.set_color("c_ann", rgb)
    cmd.color("c_ann", "IC_mouths")
    made.append("IC_mouths")

    cmd.distance("IC_line", "IC_mouths and index 1", "IC_mouths and index 2")
    cmd.hide("labels", "IC_line")
    for k, v in (("dash_gap", 0.35), ("dash_length", 0.5), ("dash_radius", 0.14)):
        cmd.set(k, v, "IC_line")
    cmd.color("c_ann", "IC_line")
    made.append("IC_line")

    # segment X: the O_IC arm, as an arrow ending at the upstream exit P2.
    # Its DIRECTION is the real one -- the chain's last two spheres, which at
    # the normalised 1 A spacing is exactly "second-to-last to last" -- and
    # only its length is drawn long enough to see.
    d = (p1 - p2) / float(np.linalg.norm(p1 - p2))
    tail, neck = p2 + d * 7.0, p2 + d * 2.4
    cmd.load_cgo([CYLINDER, *tail, *neck, 0.22, *rgb, *rgb,
                  CONE, *neck, *p2, 0.55, 0.0, *rgb, *rgb, 1.0, 1.0], "IC_arm")
    made.append("IC_arm")

    # The downstream channel's own axis is NOT drawn (user, 2026-10-25). It
    # is computed and reported as `o_ic_dn` so it can become a third term in
    # a future score, but a second arm in the figure would imply the metric
    # uses it, and it does not.

    # the O_IC arc at P2, slerped between the two arms
    u1 = (p1 - p2) / np.linalg.norm(p1 - p2)
    u3 = (p3 - p2) / np.linalg.norm(p3 - p2)
    phi = float(np.arccos(max(-1.0, min(1.0, float(u1 @ u3)))))
    # Drawn only where there IS a bend. The slerp divides by sin(phi), so it
    # degenerates at BOTH ends: near 0 there is no angle to show, and near
    # 180 sin(phi) is ~1e-16 and every intermediate point collapses onto the
    # vertex -- which draws a spike through P2 that reads as an arc. Guarding
    # only the low end looked right for years of poses at 135-165 deg and
    # failed the moment the hotspot map asked for its reference placement,
    # which is at exactly 180.
    if phi >= np.pi - 1e-3:
        # Exactly 180: the arm and the line are antiparallel, so the plane
        # the angle turns in is UNDEFINED and the slerp divides by ~1e-16.
        # The angle is still real and worth drawing, so a plane is CHOSEN --
        # the view plane, where it reads -- and the arc is a half turn in it.
        # This is a presentational choice, not a measurement: the figure
        # cannot tell you which way a 180 degree angle opens, because nothing
        # can.
        n = (np.asarray(frame.toward_viewer, float) if frame is not None
             else np.array([0.0, 0.0, 1.0]))
        n = n - float(n @ u1) * u1
        if float(np.linalg.norm(n)) < 1e-6:
            n = np.cross(u1, [1.0, 0.0, 0.0])
        n = n / float(np.linalg.norm(n))
        pts = [p2 + 4.0 * (np.cos(a) * u1 + np.sin(a) * n)
               for a in np.linspace(0.0, np.pi, 40)]
        cgo = []
        for a, b in zip(pts[:-1], pts[1:]):
            cgo += [CYLINDER, *a, *b, 0.16, *rgb, *rgb]
        cmd.load_cgo(cgo, "IC_angle")
        made.append("IC_angle")
    elif 1e-3 < phi < np.pi - 1e-3:
        pts = [p2 + 4.0 * (np.sin((1 - t) * phi) * u1 + np.sin(t * phi) * u3) / np.sin(phi)
               for t in np.linspace(0, 1, 40)]
        cgo = []
        for a, b in zip(pts[:-1], pts[1:]):
            cgo += [CYLINDER, *a, *b, 0.16, *rgb, *rgb]
        cmd.load_cgo(cgo, "IC_angle")
        made.append("IC_angle")

    # labels, offset perpendicular to the line in the view plane. PyMOL labels
    # cannot do subscripts, so they read DIC / OIC / DO; the SVG writes the
    # real subscripts.
    # Offsets must separate the labels ON SCREEN. A world axis projected off
    # the line does not: the first preview stacked all four labels in one
    # illegible pile. Use the camera's own up vector when a frame is known.
    mid = (p2 + p3) / 2.0
    ln = (p3 - p2) / float(np.linalg.norm(p3 - p2))
    if frame is not None:
        up = np.asarray(frame.up, float)
    else:
        up = np.array([0.0, 0.0, 1.0])
    up = up - (up @ ln) * ln
    if float(np.linalg.norm(up)) < 1e-6:
        up = np.cross(ln, [1.0, 0.0, 0.0])
    up = up / float(np.linalg.norm(up))
    # ...and far enough along it to clear the PROTEINS. A fixed 6 A step
    # separated the four labels from each other but left all of them lying
    # over the surfaces, where black text on a shaded translucent protein is
    # barely readable. The offset is measured instead: push past the complex's
    # own projected extent along `up`, so a label cannot land on a protein
    # whatever the pose or the zoom.
    step = 6.0
    above = below = 2.0 * step
    try:
        body = cmd.get_coords("cplx and polymer")
    except Exception:
        body = None
    if body is not None and len(body):
        reach = (body - mid) @ up
        margin = 4.0
        above = float(reach.max()) + margin
        below = -float(reach.min()) + margin
    texts = [("IC_label_D", mid + up * above, "DIC = %.1f A" % geom["d_ic"]),
             ("IC_label_O", mid + up * (above + step), "OIC = %.1f deg" % geom["o_ic"]),
             ("IC_label_DO", mid - up * below, "DO = %.3f" % geom["do"])]
    if label:
        texts.append(("IC_label_tag", mid - up * (below + step), label))
    for name, pos, txt in texts:
        cmd.pseudoatom(name, pos=[float(pos[0]), float(pos[1]), float(pos[2])])
        cmd.hide("everything", name)
        cmd.label(name, '"%s"' % txt)
        made.append(name)
    cmd.set("label_size", 14)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("float_labels", "on")
    cmd.group("IC_annotation", " ".join(made))
    return made

