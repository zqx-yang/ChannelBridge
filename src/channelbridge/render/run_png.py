# pymol -cq run_png.py -- pse=<file> out=<png> [scene=review|base|mask|mouth|footprint] [w=2000] [h=1500] [bg=FFFFFF] [transparent=1]
#
# Ray-traced PNG of one scene of a session built by run_pse.py. The .pse
# carries the style, so this only recalls a scene, switches the flat helper
# scenes to flat lighting, and renders.
#
# Called unconditionally: a file PyMOL execfile()s does not get
# __name__ == "__main__".
import sys

from pymol import cmd

a = dict(x.split("=", 1) for x in sys.argv[1:] if "=" in x)
for k in ("pse", "out"):
    if k not in a:
        raise SystemExit("missing %s=; see the header of this file" % k)

FLAT = ("mask", "mouth", "footprint")
# Scenes rendered WITHOUT the black outline. PyMOL scenes store which objects
# are visible, NOT global settings, so ray_trace_mode cannot differ per scene
# inside one session -- the renderer has to set it per scene, as it already
# does for the flat helpers. The cross-section needs this because a cut
# surface is nothing but creases and the outline mode draws a black edge
# along every one of them, burying the thing the section exists to show.
# p5d renders its own cross_section at ray_trace_mode 0 for the same reason.
NO_OUTLINE = ("classes", "occupancy", "map", "envelopes")

#: Scenes whose labels must be drawn in front of the surfaces. A label sits at
#: its own anchor's depth, so one naming a patch on the far side of an enzyme
#: renders behind the near surface and washes out. The session sets this too,
#: but a setting is global, not per scene, and a render is free to have been
#: given a session built before the setting existed.
FLOAT_LABELS = ("map", "classes", "envelopes")
cmd.load(a["pse"])
scene = a.get("scene")
if scene:
    cmd.scene(scene, "recall")
    if scene in FLAT:
        # ray_trace_mode 1 would outline the bumpy tunnel surface as creases;
        # these renders exist to be thresholded, so light them flat.
        for k_, v_ in (("ray_trace_mode", 0), ("ambient", 1.0), ("direct", 0.0),
                       ("reflect", 0.0), ("specular", 0.0), ("spec_reflect", 0.0),
                       ("ray_shadows", 0), ("ray_trace_fog", 0), ("depth_cue", 0),
                       ("antialias", 2)):
            cmd.set(k_, v_)
    if scene in FLOAT_LABELS:
        # float_labels helps only in the viewport -- a ray trace honours
        # depth regardless -- so the session ALSO lifts the anchors. The
        # orthographic camera is what makes that lift free on screen, and it
        # is a setting, so it is re-applied here: a scene stores object
        # visibility, never globals.
        cmd.set("float_labels", 1)
        cmd.set("orthoscopic", 1)
    if scene in NO_OUTLINE:
        for k_, v_ in (("ray_trace_mode", 0), ("ray_interior_color", "grey50"),
                       ("two_sided_lighting", 1), ("ambient", 0.14),
                       ("direct", 0.45), ("reflect", 0.45), ("specular", 1.0),
                       ("ray_shadows", 1), ("depth_cue", 1), ("antialias", 2)):
            cmd.set(k_, v_)
else:
    # the annotation lives in the .pse only, never in a render
    for o in cmd.get_names("all"):
        if o.startswith("IC_"):
            cmd.disable(o)

if a.get("transparent"):
    cmd.set("ray_opaque_background", 0)
hx = a.get("bg", "FFFFFF").lstrip("#")
cmd.set_color("bg_custom", [int(hx[i:i + 2], 16) / 255.0 for i in (0, 2, 4)])
cmd.bg_color("bg_custom")
cmd.png(a["out"], int(a.get("w", 2000)), int(a.get("h", 1500)), dpi=200, ray=1)
print("[png] %s  scene=%s  %sx%s" % (a["out"], scene or "(as saved)",
                                     a.get("w", 2000), a.get("h", 1500)))
