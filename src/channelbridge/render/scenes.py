"""Build ONE .pse per camera frame, with every variant as an object + scene.

Run under PyMOL, which is where ``cmd`` exists::

    PYTHONPATH=<repo>/src /opt/conda/envs-share/pymol-script/bin/pymol -cq \
        -m channelbridge.render.scenes -- \
        complex=geom001_..._complex.pdb \
        tunnel_up=.../tun_cl_001.pdb tunnel_dn=.../tun_cl_002.pdb \
        up_chains=B,C,D,E dn_chains=N \
        out=.../p6a_geom001_rotated0.pse

    # no inputs, no network: builds a synthetic pair and checks the scenes
    ... -m channelbridge.render.scenes -- --selftest

Why one session with scenes rather than a folder of files: the figure pipeline
needs the same camera for the colour render and for three flat helper renders
(silhouette, exit markers, enzyme footprint). Separate .pse files drift apart
the moment the camera is touched. The scenes are

======================  ====================================================
``review``              proteins + ligands + translucent tunnels + the
                        D_IC / O_IC / DO annotation. What opens.
``base``                proteins + ligands only -- what the 2D overlay is
                        laid over.
``mask``                the tunnel volumes, flat black -> silhouette
``mouth``               one disc per exit, upstream black / downstream grey
                        -> where the outline is left open, and which is which
``footprint``           both enzymes, flat black -> the region outside which
                        no tunnel colour may be drawn
======================  ====================================================

**The tunnels stay as objects** (``tunnel_up`` / ``tunnel_dn``), drawn at
CAVER's own per-sphere radii, so the session is usable as a reference on its
own and not only as a step toward a figure.

Only the camera moves. The coordinates of a posed complex are never
transformed -- a figure has to stay comparable with the PDB the pose stage
wrote.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .cut import (
    capsule_points,
    parse_hotspot_residues,
    patch_centroid,
)
from .frame import Frame, angle_frame, clearest_roll, pymol_view
from .select import parse_variant
from .style import RAY_SETTINGS, RenderStyle, lighten

__all__ = ["SCENE_NAMES", "SceneSpec", "build_session", "main"]

#: Always built. "classes" and "occupancy" are the inspection scenes;
#: base/mask/mouth/footprint are the flat helpers the SVG is traced from.
# "classes" is deliberately NOT here. The per-pose session answers "what does
# this pose look like"; the residue classification belongs to the hotspot map,
# which is built once per pair and shows every patch (user, 2026-10-24). Pass
# classes=1 to put it back in a pose session.
BASE_SCENES = ("occupancy", "base", "mask", "mouth", "footprint")
#: Rasters an SVG is assembled on: the same frame as the matching render
#: scene, but with the 3D channel and the 3D annotation off, because the SVG
#: supplies both as vectors.
SVG_BASE_SCENES = ("base_pair", "base_assembly")
#: Helper scenes exist per variant because the assembly is framed on the whole
#: assembly while the pair is framed on the focus pair. An SVG traces its
#: overlay from these, so a helper rendered with the WRONG camera puts the
#: channel outline and the D_IC/O_IC lines somewhere other than the proteins --
#: which is exactly what happened to the first assembly SVG.
HELPER_SCENES = ("mask", "mouth", "footprint")
#: Built per variant: render_pair, render_assembly.
VARIANT_SCENES = ("render",)
SCENE_NAMES = BASE_SCENES + VARIANT_SCENES + SVG_BASE_SCENES


def _cmd():
    """The PyMOL command object, with a useful message when absent."""
    try:
        from pymol import cmd
    except ImportError:  # pragma: no cover - depends on the interpreter
        raise SystemExit(
            "this module needs PyMOL. Run it as:\n"
            "  PYTHONPATH=<repo>/src /opt/conda/envs-share/pymol-script/bin/pymol "
            "-cq -m channelbridge.render.scenes -- <args>"
        ) from None
    return cmd


class SceneSpec:
    """What goes into one session: the structure, the tunnels, the roles."""

    def __init__(
        self,
        complex_pdb: str,
        *,
        tunnel_up: Optional[str] = None,
        tunnel_dn: Optional[str] = None,
        tunnel_up_ref: Optional[str] = None,
        tunnel_dn_ref: Optional[str] = None,
        up_chains: Sequence[str] = (),
        dn_chains: Sequence[str] = (),
        ligand_resn: Sequence[str] = (),
        style: Optional[RenderStyle] = None,
        oic_arm: float = 1.0,
        oic_error_budget: Optional[float] = None,
        roll_about: str = "arm",
        label: str = "",
        focus_chains: Sequence[str] = (),
        zoom: str = "complex",
        zoom_buffer: float = 10.0,
        occupancy_pdb: Optional[str] = None,
        occ_up_chains: Sequence[str] = (),
        occ_dn_chains: Sequence[str] = (),
        hotspots_up: Optional[str] = None,
        hotspots_dn: Optional[str] = None,
        patch_up: Optional[str] = None,
        patch_dn: Optional[str] = None,
        binder_radius: float = 11.0,
        render: str = "pair",
        classes: bool = False,
        channel_ligand_up: Sequence[str] = (),
        channel_ligand_dn: Sequence[str] = (),
        ligand_wrap: float = 3.0,
    ):
        self.complex_pdb = complex_pdb
        self.tunnel_up = tunnel_up
        self.tunnel_dn = tunnel_dn
        # the structure each tunnel was computed on, if the pose moved that
        # target; the tunnel is carried onto the complex by superposition
        self.tunnel_up_ref = tunnel_up_ref
        self.tunnel_dn_ref = tunnel_dn_ref
        self.up_chains = [c for c in up_chains if c]
        self.dn_chains = [c for c in dn_chains if c]
        self.ligand_resn = [r for r in ligand_resn if r]
        self.style = style or RenderStyle()
        self.oic_arm = float(oic_arm)
        #: None keeps the angle frame, where O_IC is drawn exactly true and
        #: there is no freedom left. A number switches to the D_IC-line frame
        #: and allows that many degrees of O_IC error in exchange for rolling
        #: the channels off the line.
        self.oic_error_budget = (None if oic_error_budget is None
                                 else float(oic_error_budget))
        #: Which axis the camera rolls about: "arm" keeps the O_IC arm on x
        #: (the angle frame's convention), "line" keeps the D_IC line on x.
        self.roll_about = str(roll_about)
        self.label = label
        # Only these chains are drawn, if given. A GDH tetramer plus eight
        # ligand chains is 15 transparent objects and cost 22 min of CPU for
        # an 800x600 preview; the tunnel-bearing protomer plus the partner is
        # the same figure for a fraction of that.
        self.focus_chains = [c for c in focus_chains if c]
        # "complex" frames both enzymes whole, the p5c look and the default;
        # "channel" crops tight to the two mouths. With focus_chains trimming
        # the assembly to the tunnel-bearing protomer plus its partner, the
        # whole-complex view is no longer dominated by symmetry copies.
        self.zoom = zoom
        self.zoom_buffer = float(zoom_buffer)
        # the full assembly with symmetry copies, for the assembly variants
        self.occupancy_pdb = occupancy_pdb
        # every protomer of each enzyme in the assembly. One binder per
        # upstream protomer, so the occupancy scene must show them ALL --
        # a single capsule cannot test whether neighbouring binders collide,
        # which is the only reason the scene exists.
        self.occ_up_chains = [c for c in occ_up_chains if c] or list(self.up_chains)
        self.occ_dn_chains = [c for c in occ_dn_chains if c] or list(self.dn_chains)
        # <NAME>_hotspot_patches.tsv from the pose stage
        self.hotspots_up = hotspots_up
        self.hotspots_dn = hotspots_dn
        # which patch of each target the pose bridges, e.g. "4" from GDHp4
        self.patch_up = patch_up
        self.patch_dn = patch_dn
        self.binder_radius = float(binder_radius)
        self.render = parse_variant(render, "render")
        # The species that actually travels the channel, per enzyme. Only this
        # one is wrapped into the drawn volume: wrapping every cofactor would
        # paint red over prosthetic groups that the channel never reaches.
        # the classification scene, off by default (see BASE_SCENES)
        self.classes = bool(classes)
        self.channel_ligand_up = [r for r in channel_ligand_up if r]
        self.channel_ligand_dn = [r for r in channel_ligand_dn if r]
        self.ligand_wrap = float(ligand_wrap)

    def chain_sel(self, which: str) -> str:
        chains = self.up_chains if which == "up" else self.dn_chains
        if not chains:
            raise ValueError("no chains given for the %sstream enzyme" % which)
        return "chain " + "+".join(chains)

    def channel_ligand_sel(self, which: str) -> str:
        """Selection for the channelled species of one enzyme, or ``""``.

        Falls back to that enzyme's chains intersected with every ligand when
        no channelled species was named, which is right for a target with a
        single cofactor and over-generous for one with several -- hence the
        explicit option.
        """
        resn = self.channel_ligand_up if which == "up" else self.channel_ligand_dn
        if resn:
            return "resn " + "+".join(resn)
        if not self.ligand_resn:
            return ""
        try:
            return "(%s) and (%s)" % (self.chain_sel(which), self.ligand_sel)
        except ValueError:
            return ""

    @property
    def ligand_sel(self) -> str:
        return "resn " + "+".join(self.ligand_resn) if self.ligand_resn else "none"


# ------------------------------------------------------------------ geometry


def read_tunnel_spheres(path) -> Tuple[np.ndarray, np.ndarray]:
    """Sphere centres and radii from a CAVER tunnel-cluster PDB.

    Independent of :mod:`channelbridge.io.pdb` on purpose: this module runs
    inside PyMOL's interpreter, and a figure script should not fail because
    an unrelated import moved.
    """
    pts, rad = [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            b = line[60:66].strip()
            o = line[54:60].strip()
            try:
                rad.append(float(b) if b else (float(o) if o else 1.0))
            except ValueError:
                rad.append(1.0)
    if len(pts) < 2:
        raise ValueError(
            "%s holds %d sphere(s); a tunnel needs at least 2. An empty cluster "
            "file usually means CAVER found no pathway." % (path, len(pts))
        )
    return np.asarray(pts, float), np.asarray(rad, float)


def read_ca(path, chain) -> np.ndarray:
    """CA coordinates of one chain, in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain or line[12:16].strip() != "CA":
                continue
            out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(out, float)


def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Rigid transform taking *mobile* onto *target*. Returns ``(R, t, rmsd)``.

    No reflection: the determinant is corrected, so a near-degenerate set
    cannot mirror the tunnel into the protein.
    """
    mobile = np.asarray(mobile, float)
    target = np.asarray(target, float)
    if mobile.shape != target.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError(
            "kabsch needs two matching (N, 3) arrays, got %r and %r"
            % (mobile.shape, target.shape)
        )
    if len(mobile) < 3:
        raise ValueError("kabsch needs at least 3 points, got %d" % len(mobile))
    cm, ct = mobile.mean(0), target.mean(0)
    h = (mobile - cm).T @ (target - ct)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = ct - r @ cm
    rmsd = float(np.sqrt((((r @ mobile.T).T + t - target) ** 2).sum(axis=1).mean()))
    return r, t, rmsd


def pose_transform(ref_pdb, complex_pdb, chain):
    """The transform the pose stage applied to a target.

    The pose stage leaves the upstream target in place and moves the
    downstream one, so a tunnel computed on the ORIGINAL target has to be
    carried onto the posed complex before it can be drawn. The atom
    correspondence is one-to-one (same file lineage, same numbering), so this
    is exact and the reported RMSD should be ~0; anything else means the
    wrong reference was given.
    """
    a = read_ca(ref_pdb, chain)
    b = read_ca(complex_pdb, chain)
    if len(a) == 0 or len(b) == 0:
        raise ValueError(
            "chain %r has no CA atoms in %s"
            % (chain, ref_pdb if len(a) == 0 else complex_pdb)
        )
    if len(a) != len(b):
        raise ValueError(
            "chain %r has %d CA in the reference and %d in the complex; the "
            "reference does not match. Give the target the tunnel was computed "
            "on." % (chain, len(a), len(b))
        )
    return kabsch(a, b)


def arc_point(xyz: np.ndarray, arm: float) -> np.ndarray:
    """The point *arm* angstroms inward from the outer end, along the chain.

    The same fixed-arc-length rule as
    :meth:`channelbridge.tunnel.Tunnel.axis_point`, repeated here so a figure
    measures O_IC over the same baseline the score used.
    """
    seg = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(s[-1])
    if arm >= total:
        return xyz[0].copy()
    target = total - arm
    return np.array([np.interp(target, s, xyz[:, k]) for k in range(3)])


# ----------------------------------------------------------------- building


def _pseudo_spheres(cmd, name, pts, radii, rgb, transparency):
    """One object of spheres at their true radii."""
    cmd.delete(name)
    for _i, (p, r) in enumerate(zip(np.atleast_2d(pts), np.atleast_1d(radii))):
        cmd.pseudoatom(name, pos=[float(p[0]), float(p[1]), float(p[2])],
                       vdw=float(max(r, 0.1)), state=1)
    cmd.hide("everything", name)
    cmd.show("spheres", name)
    cmd.set_color("c_" + name, [float(c) for c in rgb])
    cmd.color("c_" + name, name)
    cmd.set("sphere_transparency", float(transparency), name)
    return cmd.count_atoms(name)


def _apply(cmd, settings):
    for k, v in settings.items():
        cmd.set(k, v)


def build_session(spec: SceneSpec, out_path: str, frame: Optional[Frame] = None) -> Dict:
    """Build the session and save it. Returns a report dict."""
    cmd = _cmd()
    st = spec.style

    cmd.reinitialize()
    cmd.set("auto_zoom", 0)          # else every added object jumps the view
    cmd.load(spec.complex_pdb, "cplx_full")
    if spec.focus_chains:
        cmd.create("cplx", "cplx_full and chain %s" % "+".join(spec.focus_chains))
        cmd.delete("cplx_full")
    else:
        cmd.set_name("cplx_full", "cplx")

    report = {"complex": os.path.basename(spec.complex_pdb), "objects": [], "scenes": []}

    # ---- colours and representation
    #
    # Two families of the same structure, because the scenes have different
    # jobs. "classes" and "occupancy" are for INSPECTION -- where the surfaces
    # and hotspots actually sit -- and the render scenes are ILLUSTRATION.
    # Separate objects rather than per-scene settings: what a PyMOL scene
    # stores of object settings is version-dependent, an object identity is not.
    #
    # The inspection copy is translucent too, just less so (role "inspect"):
    # a solid surface hides the cofactors, and a hotspot patch only means
    # something next to the ligand the channel serves.
    inspect_t = st.transparency_for("inspect")
    for role, which in (("upstream", "up"), ("downstream", "dn")):
        sel = "cplx and %s and polymer" % spec.chain_sel(which)
        cmd.set_color("c_" + role, [float(c) for c in st.palette.rgb(role)])
        cmd.color("c_" + role, sel)
        cmd.set("transparency", inspect_t, sel)
    cmd.set("transparency", inspect_t, "cplx")     # object level too, so it is checkable
    if spec.ligand_resn:
        cmd.set_color("c_ligand", [float(c) for c in st.palette.rgb("ligand")])
        cmd.color("c_ligand", "cplx and (%s)" % spec.ligand_sel)
    cmd.hide("everything", "all")
    cmd.show("surface", "cplx and polymer")
    if spec.ligand_resn:
        cmd.show("sticks", "cplx and (%s)" % spec.ligand_sel)
        cmd.set("stick_radius", st.stick_radius)
        cmd.set("stick_ball", 0)
    cmd.bg_color(st.background if st.background == "white" else "0x" + st.background)
    _apply(cmd, RAY_SETTINGS)

    # the translucent twin used by the render scenes
    cmd.create("cplx_il", "cplx")
    cmd.set("transparency", st.transparency_for("upstream"), "cplx_il")
    for role, which in (("upstream", "up"), ("downstream", "dn")):
        cmd.set("transparency", st.transparency_for(role),
                "cplx_il and %s and polymer" % spec.chain_sel(which))
    report["objects"].append(
        "cplx (transparency %.2f, for classes/occupancy) + cplx_il "
        "(%.2f, for the render scenes)"
        % (inspect_t, st.transparency_for("upstream")))

    # ---- tunnels as objects, at CAVER's own radii
    #
    # TWO sets are kept, and conflating them would be a quiet disaster:
    #   `tun`   = the CAVER sphere chain alone. Every NUMBER comes from this
    #             -- P1, P2, P3, D_IC, O_IC, the mouth marker.
    #   `drawn` = what is painted: the chain PLUS the channelled ligand
    #             wrapped at `ligand_wrap`. A CAVER chain starts at the
    #             starting point and runs outward, so on its own it leaves the
    #             cofactor sitting outside the red volume -- but the ligand is
    #             IN the channel, and p5c draws the union for exactly this
    #             reason. Visual only: it must never move an exit.
    tun = {}
    drawn = {}
    for key, path, ref in (("up", spec.tunnel_up, spec.tunnel_up_ref),
                           ("dn", spec.tunnel_dn, spec.tunnel_dn_ref)):
        if not path:
            continue
        xyz, rad = read_tunnel_spheres(path)
        if ref:
            chains = spec.up_chains if key == "up" else spec.dn_chains
            rot, trans, rmsd = pose_transform(ref, spec.complex_pdb, chains[0])
            xyz = (rot @ xyz.T).T + trans
            report.setdefault("superposed", []).append(
                "%s tunnel carried onto the pose by chain %s, rmsd %.4f A"
                % (key, chains[0], rmsd))
            if rmsd > 0.05:
                report.setdefault("warnings", []).append(
                    "%s tunnel superposition rmsd %.3f A is not ~0; is %s the "
                    "target the tunnel was computed on?" % (key, rmsd, ref))
        tun[key] = (xyz, rad)
        d_xyz, d_rad = xyz, rad
        wrap_sel = spec.channel_ligand_sel(key)
        if wrap_sel and spec.ligand_wrap > 0:
            try:
                lig = np.asarray(cmd.get_coords("cplx and (%s)" % wrap_sel), float)
            except Exception:
                lig = None
            if lig is not None and len(lig):
                d_xyz = np.vstack([xyz, lig])
                d_rad = np.concatenate([rad, np.full(len(lig), spec.ligand_wrap)])
                report["objects"].append(
                    "tunnel_%s wraps %d ligand atom(s) of (%s) at %.1f A"
                    % (key, len(lig), wrap_sel, spec.ligand_wrap))
            else:
                report.setdefault("warnings", []).append(
                    "tunnel_%s: no atoms matched the channelled ligand (%s), so "
                    "the drawn channel will not reach the cofactor"
                    % (key, wrap_sel))
        drawn[key] = (d_xyz, d_rad)
        n = _pseudo_spheres(cmd, "tunnel_" + key, d_xyz, d_rad,
                            st.palette.rgb("tunnel"), st.transparency_for("tunnel"))
        report["objects"].append("tunnel_%s (%d spheres drawn, %d from CAVER)"
                                 % (key, n, len(xyz)))

    # ---- the iMARS construction, if both tunnels are present
    geom = None
    if "up" in tun and "dn" in tun:
        up_xyz = tun["up"][0]
        dn_xyz = tun["dn"][0]
        p2 = up_xyz[-1]
        p1 = arc_point(up_xyz, spec.oic_arm)
        p3 = dn_xyz[-1]
        # the downstream channel's OWN axis, the same arm length back along
        # its chain. iMARS never uses it -- see o_ic_dn below -- but without
        # it the figure cannot show which way the receiving mouth opens.
        q1 = arc_point(dn_xyz, spec.oic_arm)
        d_ic = float(np.linalg.norm(p3 - p2))
        v1, v3 = p1 - p2, p3 - p2
        cosv = float(np.dot(v1, v3) / (np.linalg.norm(v1) * np.linalg.norm(v3)))
        o_ic = float(np.degrees(np.arccos(max(-1.0, min(1.0, cosv)))))
        # The MIRROR of O_IC, measured at the downstream exit. O_IC is
        # one-sided: it asks only whether the DONOR channel is aimed at the
        # acceptor's exit, and says nothing about which way the acceptor
        # opens. A receiving mouth pointing 90 degrees sideways, or straight
        # back into its own protein, scores exactly like one facing the
        # donor. On the p6a GDH-FNR pose: O_IC 165 deg, this 90 deg.
        # Reported, never scored -- DO's published coefficients cover D_IC
        # and O_IC only, and adding a term here would not be iMARS.
        w1, w3 = q1 - p3, p2 - p3
        cosw = float(np.dot(w1, w3) / (np.linalg.norm(w1) * np.linalg.norm(w3)))
        o_ic_dn = float(np.degrees(np.arccos(max(-1.0, min(1.0, cosw)))))
        from ..geometry import do_score
        geom = {"d_ic": d_ic, "o_ic": o_ic, "do": do_score(d_ic, o_ic),
                "o_ic_dn": o_ic_dn, "p1": p1, "p2": p2, "p3": p3, "q1": q1}
        report["geometry"] = {k: round(v, 4) for k, v in
                              (("d_ic", d_ic), ("o_ic", o_ic), ("do", geom["do"]),
                               ("o_ic_downstream_not_scored", o_ic_dn))}

    # ---- camera: the angle frame, else whatever PyMOL fitted
    if frame is None and geom is not None:
        centre = np.asarray(cmd.get_extent("cplx")).mean(axis=0)
        try:
            if spec.oic_error_budget is not None:
                # Roll the camera until the drawn channels stop lying
                # across the D_IC line. The axis is the O_IC ARM by default,
                # so the arm stays on x exactly as in the angle frame and
                # roll 0 is that frame. Costs some O_IC fidelity, capped by
                # the budget and REPORTED either way -- an angle a reader
                # cannot measure off the figure must not be presented as if
                # they could.
                pts = [cmd.get_coords(n) for n in ("tunnel_up", "tunnel_dn")
                       if n in cmd.get_names("objects")]
                pts = np.vstack(pts) if pts else np.zeros((0, 3))
                roll, frame, clear, shown, true, frac = clearest_roll(
                    geom["p1"], geom["p2"], geom["p3"], pts,
                    about=spec.roll_about,
                    max_angle_error=spec.oic_error_budget, step=5.0,
                    centre=centre)
                geom["oic_drawn"] = shown
                report["frame"] = {
                    "mode": "roll about the %s" % spec.roll_about,
                    "roll_deg": round(roll, 1),
                    "clearance_A": round(clear, 2),
                    "o_ic_true": round(true, 1), "o_ic_drawn": round(shown, 1),
                    "o_ic_error": round(abs(shown - true), 1),
                    "d_ic_drawn_pct": round(100 * frac, 1)}
            else:
                frame = angle_frame(geom["p1"], geom["p2"], geom["p3"],
                                    centre=centre)
                report["frame"] = {"mode": "angle", "o_ic_drawn": "true"}
        except ValueError as exc:
            report["frame_warning"] = str(exc)
            frame = None
    if frame is not None and geom is not None:
        # The projected P2->P1 direction, for the SVG's angle marker. Under
        # the angle frame it can be derived from the reported O_IC; under a
        # rolled frame it cannot, and assuming it would draw the arm at an
        # angle the picture does not contain.
        arm = np.asarray(geom["p1"], float) - np.asarray(geom["p2"], float)
        ax = float(arm @ frame.axes[0])
        ay = -float(arm @ frame.axes[1])        # image y grows downward
        n = float(np.hypot(ax, ay))
        if n > 1e-9:
            report["axis_px"] = [round(ax / n, 6), round(ay / n, 6)]
    zoom_sel = "cplx"
    if spec.zoom == "channel" and tun:
        names = [n for n in ("tunnel_up", "tunnel_dn") if n in cmd.get_names("objects")]
        if names:
            zoom_sel = " or ".join(names)
    if frame is not None:
        cmd.set_view(pymol_view(frame))
        cmd.zoom(zoom_sel, buffer=spec.zoom_buffer, complete=0)
        cmd.clip("slab", 3000)
        report["frame"] = frame.name
    else:
        cmd.reset()
        cmd.zoom(zoom_sel, buffer=spec.zoom_buffer, complete=0)
        report["frame"] = "fitted"
    report["zoom"] = "%s (buffer %.0f A)" % (zoom_sel, spec.zoom_buffer)
    view = cmd.get_view()

    # ---- annotation objects (review only; never rendered to png)
    ann = []
    if geom is not None:
        ann = _annotation(cmd, geom, st, spec.label, frame)
        report["objects"].extend(ann)

    # ---- flat helper objects
    helpers = _helpers(cmd, spec, tun, st, drawn)
    report["objects"].extend(helpers)

    # ---- the four named scenes plus the flat helpers
    tunnels = [n for n in ("tunnel_up", "tunnel_dn") if n in cmd.get_names("objects")]
    classes_objs = _classes(cmd, spec, st, report) if spec.classes else []
    occ_objs = _occupancy(cmd, spec, st, report)

    sets = {
        **({"classes": ["cplx", *classes_objs]} if spec.classes else {}),
        "occupancy": (occ_objs or ["cplx"]),
        "base": ["cplx_il"],
        "mask": [n for n in helpers if n.startswith("mask_")],
        "mouth": [n for n in helpers if n.startswith("mouth_")],
        "footprint": [n for n in helpers if n.startswith("footprint")],
        # the same three again, to be stored with the assembly's camera
        "mask_assembly": [n for n in helpers if n.startswith("mask_")],
        "mouth_assembly": [n for n in helpers if n.startswith("mouth_")],
        "footprint_assembly": [n for n in helpers if n.startswith("footprint")],
    }
    # render, per requested variant
    assembly_view = None
    variants = {"pair": ("cplx", ""), "assembly": ("cplx_assembly", "_assembly")}
    # The cross-section was dropped (user, 2026-10-25): it was not meaningfully
    # different from render_pair, and the cut copy coloured parts of the FNR
    # ligand black. render_pair shows the tunnels and the relative positions,
    # which is what the scene was for. Its builder went with it; the cutting
    # geometry stays in render/cut.py, which is tested on its own.
    for base_name, want in (("render", spec.render),):
        if want == "none":
            continue
        which = ("pair", "assembly") if want == "both" else (want,)
        for v in which:
            src, _sfx = variants[v]
            if src not in cmd.get_names("objects"):
                report.setdefault("warnings", []).append(
                    "%s_%s skipped: no %s object (give occupancy= for the assembly)"
                    % (base_name, v, src))
                continue
            key = base_name if want != "both" else "%s_%s" % (base_name, v)
            if v == "assembly" and src in cmd.get_names("objects"):
                # The assembly is far bigger than the focus pair, so the pair's
                # zoom clips it: the first assembly SVG came out with protomers
                # cut off at the frame edge. Re-fit to the assembly, keeping
                # the rotation, and store that view with the assembly scenes.
                cmd.set_view(view)
                cmd.zoom(src, buffer=6, complete=1)
                cmd.clip("slab", 3000)
                assembly_view = cmd.get_view()
            if base_name == "render":
                sets[key] = [src + "_il", *tunnels, *ann]
                # The same frame WITHOUT the 3D channel and annotation. This
                # is what an SVG is assembled on: the SVG draws the channel as
                # a traced vector outline and the numbers as editable text, so
                # a raster that already carries both comes out with every
                # label twice and a red smear where the 3D tunnel reaches into
                # solvent.
                sets["base_" + v] = [src + "_il"]


    for name in list(sets):
        objs = [o for o in sets.get(name, []) if o in cmd.get_names("objects")]
        if not objs:
            continue
        cmd.disable("all")
        for o in objs:
            cmd.enable(o)
        # every *_assembly scene -- base, render AND the helpers an SVG traces
        # from -- must share the assembly camera, or the vector overlay lands
        # in a different place from the raster it is drawn over
        use_assembly = assembly_view and ("assembly" in name or name.endswith("_assembly"))
        cmd.set_view(assembly_view if use_assembly else view)
        cmd.scene(name, "store", _scene_note(name))
        report["scenes"].append("%s (%s)" % (name, ", ".join(objs)))

    first = next((s for s in ("render", "render_pair", "classes", "base")
                  if s in cmd.get_scene_list()), None)
    if first:
        cmd.scene(first, "recall")
    cmd.set_view(view)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    cmd.save(out_path)
    report["out"] = out_path
    report["needs_occlusion_clipping"] = st.needs_occlusion_clipping
    return report


def _scene_note(name):
    return {
        "review": "proteins + ligands + tunnels + D_IC/O_IC/DO annotation",
        "base": "proteins and ligands only; the 2D overlay is laid over this",
        "mask": "tunnel volumes, flat black -> silhouette",
        "mouth": "exit markers, upstream black / downstream grey",
        "footprint": "both enzymes, flat black -> the no-draw boundary",
        "base_pair": "proteins and ligands only -- the raster an SVG is built on",
        "base_assembly": "the assembly, proteins and ligands only -- SVG raster",
    }.get(name, name)


def read_hotspot_patches(path):
    """Rank -> residue list, from a ``<NAME>_hotspot_patches.tsv``."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            i_rank = hdr.index("rank")
            i_res = hdr.index("hotspot_residues")
        except ValueError:
            return out
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_rank, i_res):
                continue
            out[str(f[i_rank]).strip()] = parse_hotspot_residues(f[i_res])
    return out


def _patch_key(patches, wanted):
    """The key of *wanted* in a patch table, tolerating ``"4"`` vs ``"04"``.

    The pose stage writes the patch into a file name (``GDHp4``) and the table
    writes it as a rank; both are strings whose spelling has drifted before,
    so the lookup compares numbers when it can and text otherwise.
    """
    if wanted is None:
        return None
    w = str(wanted).strip()
    if w in patches:
        return w
    try:
        n = int(w)
    except ValueError:
        return None
    for key in patches:
        try:
            if int(str(key).strip()) == n:
                return key
        except ValueError:
            continue
    return None


def _ca_by_resi(cmd, sel):
    m = cmd.get_model("(%s) and name CA" % sel)
    return {a.resi_number: tuple(a.coord) for a in m.atom}


def _classes(cmd, spec, st, report):
    """Surfaces by role with THIS pose's hotspot patch picked out.

    Only the patch the pose is built on is coloured, not every candidate: the
    question this scene answers is "where would the binder bind in this
    pose", and six red patches on one enzyme answer a different one. The
    overall map of every patch is a separate session
    (:mod:`channelbridge.render.hotspots`).
    """
    made = []
    for which, path, focus in (("up", spec.hotspots_up, spec.patch_up),
                               ("dn", spec.hotspots_dn, spec.patch_dn)):
        patches = read_hotspot_patches(path)
        if not patches:
            continue
        try:
            chains = spec.chain_sel(which)
        except ValueError:
            continue
        chosen = patches
        if focus:
            key = _patch_key(patches, focus)
            if key is None:
                report["objects"].append(
                    "WARNING: %sstream patch %r is not in %s; colouring every patch"
                    % (which, focus, os.path.basename(str(path))))
            else:
                chosen = {key: patches[key]}
        resi = sorted({r for rs in chosen.values() for r in rs})
        if not resi:
            continue
        name = "hotspots_" + which
        cmd.create(name, "cplx and %s and resi %s and polymer"
                   % (chains, "+".join(str(r) for r in resi)))
        if cmd.count_atoms(name) == 0:
            cmd.delete(name)
            continue
        cmd.hide("everything", name)
        cmd.show("surface", name)
        cmd.set_color("c_hot", [float(c) for c in st.palette.rgb("annotation")])
        cmd.color("c_hot", name)
        # the patch itself stays solid against the translucent surface around
        # it: it is the one thing this scene is pointing at
        cmd.set("transparency", 0.0, name)
        made.append(name)
        report["objects"].append(
            "%s (patch %s of %d, %d residues)"
            % (name, ", ".join(sorted(chosen)), len(patches), len(resi)))
    return made


def _segment_distance(p1, q1, p2, q2):
    """Closest distance between two line segments -- the binder clash test."""
    p1, q1, p2, q2 = (np.asarray(x, float) for x in (p1, q1, p2, q2))
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e, f = float(d1 @ d1), float(d2 @ d2), float(d2 @ r)
    if a < 1e-9 and e < 1e-9:
        return float(np.linalg.norm(r))
    if a < 1e-9:
        s, tt = 0.0, max(0.0, min(1.0, f / e))
    else:
        c = float(d1 @ r)
        if e < 1e-9:
            tt, s = 0.0, max(0.0, min(1.0, -c / a))
        else:
            b = float(d1 @ d2)
            den = a * e - b * b
            s = max(0.0, min(1.0, (b * f - c * e) / den)) if den > 1e-9 else 0.0
            tt = (b * s + f) / e
            if tt < 0.0:
                tt, s = 0.0, max(0.0, min(1.0, -c / a))
            elif tt > 1.0:
                tt, s = 1.0, max(0.0, min(1.0, (b - c) / a))
    return float(np.linalg.norm((p1 + d1 * s) - (p2 + d2 * tt)))


def _occupancy(cmd, spec, st, report):
    """The assembly, plus the binder envelope the pose stage tested."""
    made = []
    if spec.occupancy_pdb and os.path.isfile(spec.occupancy_pdb):
        cmd.load(spec.occupancy_pdb, "cplx_assembly")
        cmd.hide("everything", "cplx_assembly")
        cmd.show("surface", "cplx_assembly and polymer")
        # Every protomer gets ITS OWN ENZYME'S colour, in a different shade,
        # rather than a uniform grey for "the copies". Grey copies read as a
        # different species and, worse, a grey copy sitting in front of the
        # coloured protomer hides the one the figure is about. Shades of the
        # same hue keep the pair legible and still tell the copies apart.
        for role, which in (("upstream", "up"), ("downstream", "dn")):
            chains = spec.occ_up_chains if which == "up" else spec.occ_dn_chains
            chains = [c for c in chains if c]
            if not chains:
                continue
            base = st.palette.rgb(role)
            # The first chain is the one the pose was built on: full colour.
            # Every OTHER protomer of that enzyme gets ONE shared lighter
            # shade of the same hue -- not a per-copy gradient. A gradient
            # said something the figure does not mean (the copies are not
            # ordered, and copy 4 is no further from the pose than copy 2),
            # and at the pale end of it the fourth upstream and fourth
            # downstream protomer both washed out to near-white and stopped
            # being tellable apart. One step keeps every copy in its own
            # enzyme's hue and clearly behind the protomer that matters.
            pale = lighten(base, st.copy_lighten)
            cmd.set_color("c_%s_copy" % role, [float(c) for c in pale])
            for i, ch in enumerate(chains):
                sel = "cplx_assembly and chain %s and polymer" % ch
                if cmd.count_atoms(sel) == 0:
                    continue
                if i == 0:
                    key = "c_%s_%s" % (role, ch.lower())
                    cmd.set_color(key, [float(c) for c in base])
                else:
                    key = "c_%s_copy" % role
                cmd.color(key, sel)
            report["objects"].append(
                "cplx_assembly %sstream: %s full colour, %s one lighter shade"
                % (which, chains[0], "+".join(chains[1:]) or "(none)"))
        if spec.ligand_resn:
            # the assembly carries its cofactors too; without this the
            # assembly render loses every ligand the pair render shows
            lig_sel = "cplx_assembly and (%s)" % spec.ligand_sel
            if cmd.count_atoms(lig_sel):
                cmd.set_color("c_ligand", [float(c) for c in st.palette.rgb("ligand")])
                cmd.color("c_ligand", lig_sel)
                cmd.show("sticks", lig_sel)
                report["objects"].append(
                    "cplx_assembly ligands (%d atoms, sticks)" % cmd.count_atoms(lig_sel))
        cmd.set("transparency", st.transparency_for("inspect"),
                "cplx_assembly and polymer")                        # the inspection copy
        cmd.create("cplx_assembly_il", "cplx_assembly")
        cmd.set("transparency", st.transparency_for("copy"), "cplx_assembly_il and polymer")
        made.append("cplx_assembly")
        report["objects"].append("cplx_assembly (%d atoms)" % cmd.count_atoms("cplx_assembly"))

    # One binder per upstream protomer. Each is paired with the NEAREST
    # downstream copy rather than by chain order, because the assembly's
    # chain naming does not have to follow the pairing the pose implies.
    src = "cplx_assembly" if "cplx_assembly" in cmd.get_names("objects") else "cplx"
    patches_up = read_hotspot_patches(spec.hotspots_up)
    patches_dn = read_hotspot_patches(spec.hotspots_dn)
    resi_up = patches_up.get(str(spec.patch_up)) if spec.patch_up else None
    resi_dn = patches_dn.get(str(spec.patch_dn)) if spec.patch_dn else None
    if resi_up and resi_dn:
        up_cents, dn_cents = {}, {}
        for ch in spec.occ_up_chains:
            c = patch_centroid(_ca_by_resi(cmd, "%s and chain %s and polymer" % (src, ch)), resi_up)
            if c is not None:
                up_cents[ch] = c
        for ch in spec.occ_dn_chains:
            c = patch_centroid(_ca_by_resi(cmd, "%s and chain %s and polymer" % (src, ch)), resi_dn)
            if c is not None:
                dn_cents[ch] = c
        pairs = []
        for uch, uc in up_cents.items():
            if not dn_cents:
                break
            dch = min(dn_cents, key=lambda k: float(np.linalg.norm(dn_cents[k] - uc)))
            pairs.append((uch, dch, uc, dn_cents[dch]))
        all_pts, all_rad = [], []
        for _uch, _dch, a, b in pairs:
            pts, rad = capsule_points(a, b, spec.binder_radius)
            all_pts.append(pts)
            all_rad.append(rad)
        if all_pts:
            _pseudo_spheres(cmd, "binder_envelope", np.vstack(all_pts),
                                np.concatenate(all_rad), st.palette.rgb("binder"), 0.6)
            made.append("binder_envelope")
            report["objects"].append(
                "binder_envelope (%d binder(s), r=%.1f A: %s)"
                % (len(pairs), spec.binder_radius,
                   ", ".join("%s-%s" % (u, d) for u, d, _, _ in pairs)))
            # closest approach between neighbouring binder axes: the clash test
            if len(pairs) > 1:
                worst = min(
                    _segment_distance(pairs[i][2], pairs[i][3], pairs[j][2], pairs[j][3])
                    for i in range(len(pairs)) for j in range(i + 1, len(pairs)))
                report["binder_binder_min"] = round(worst, 2)
                report["objects"].append(
                    "  closest binder-binder axis distance %.1f A (capsules touch below %.1f)"
                    % (worst, 2 * spec.binder_radius))
    elif spec.patch_up or spec.patch_dn:
        report.setdefault("warnings", []).append(
            "binder envelope not drawn: needs both hotspot tables and both patch ids")
    return made


def _annotation(cmd, geom, st, label, frame=None):
    """Kept as a thin alias: the drawing lives in render/annotation.py now,
    because the hotspot map draws the same thing and had grown its own copy."""
    from .annotation import build_annotation

    return build_annotation(cmd, geom, st, label, frame)


def _helpers(cmd, spec, tun, st, drawn=None):
    """Flat-black objects the 2D overlay is traced from.

    The MASK traces the drawn volume (chain + wrapped ligand), so the SVG's
    silhouette matches the 3D render. The MOUTH marker stays on the CAVER
    chain's own last sphere: it marks the exit, which the ligand wrap must
    not be allowed to move.
    """
    made = []
    drawn = drawn or tun
    for key, (xyz, rad) in tun.items():
        m_xyz, m_rad = drawn.get(key, (xyz, rad))
        _pseudo_spheres(cmd, "mask_" + key, m_xyz, m_rad, (0.0, 0.0, 0.0), 0.0)
        made.append("mask_" + key)
        # upstream black, downstream grey: tells the two exits apart in the render
        grey = (0.0, 0.0, 0.0) if key == "up" else (0.5, 0.5, 0.5)
        _pseudo_spheres(cmd, "mouth_" + key, xyz[-1][None, :],
                        [float(np.clip(rad[-1], 1.0, 1.8))], grey, 0.0)
        made.append("mouth_" + key)
    sel = []
    for which in ("up", "dn"):
        try:
            sel.append(spec.chain_sel(which))
        except ValueError:
            pass
    if sel:
        cmd.create("footprint", "cplx and (%s) and polymer" % " or ".join(sel))
        cmd.hide("everything", "footprint")
        cmd.show("surface", "footprint")
        cmd.color("black", "footprint")
        cmd.set("transparency", 0.0, "footprint")
        made.append("footprint")
    return made


# --------------------------------------------------------------------- entry


def _csv(value) -> List[str]:
    return [x for x in str(value or "").replace(",", " ").split() if x]


def _parse_argv(argv) -> Dict[str, str]:
    out = {}
    for a in argv:
        if a == "--selftest":
            out["selftest"] = "1"
        elif "=" in a:
            k, v = a.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    args = _parse_argv(argv)
    if args.get("selftest"):
        return selftest()
    missing = [k for k in ("complex", "out") if k not in args]
    if missing:
        print("missing argument(s): %s" % ", ".join(missing), file=sys.stderr)
        print(__doc__, file=sys.stderr)
        return 2
    spec = SceneSpec(
        complex_pdb=args["complex"],
        tunnel_up=args.get("tunnel_up"),
        tunnel_dn=args.get("tunnel_dn"),
        tunnel_up_ref=args.get("tunnel_up_ref"),
        tunnel_dn_ref=args.get("tunnel_dn_ref"),
        up_chains=args.get("up_chains", "").split(","),
        dn_chains=args.get("dn_chains", "").split(","),
        ligand_resn=args.get("ligands", "").split(","),
        oic_arm=float(args.get("oic_arm", 1.0)),
        oic_error_budget=(float(args["oic_error_budget"])
                          if args.get("oic_error_budget") else None),
        roll_about=args.get("roll_about", "arm"),
        label=args.get("label", ""),
        focus_chains=[c for c in args.get("focus_chains", "").split(",") if c],
        zoom=args.get("zoom", "complex"),
        zoom_buffer=float(args.get("zoom_buffer", 10.0)),
        occupancy_pdb=args.get("occupancy"),
        occ_up_chains=[c for c in args.get("occ_up_chains", "").split(",") if c],
        occ_dn_chains=[c for c in args.get("occ_dn_chains", "").split(",") if c],
        hotspots_up=args.get("hotspots_up"),
        hotspots_dn=args.get("hotspots_dn"),
        patch_up=args.get("patch_up"),
        patch_dn=args.get("patch_dn"),
        binder_radius=float(args.get("binder_radius", 11.0)),
        render=args.get("render", "pair"),
        classes=str(args.get("classes", "")).lower() in ("1", "true", "yes"),
        channel_ligand_up=_csv(args.get("channel_ligand_up")),
        channel_ligand_dn=_csv(args.get("channel_ligand_dn")),
        ligand_wrap=float(args.get("ligand_wrap", 3.0)),
    )
    rep = build_session(spec, args["out"])
    for k in ("complex", "frame", "zoom", "out", "needs_occlusion_clipping"):
        if k in rep:
            print("[pse] %-26s %s" % (k, rep[k]))
    if "geometry" in rep:
        g = rep["geometry"]
        print("[pse] %-26s D_IC %.1f A, O_IC %.1f deg, DO %.4f"
              % ("geometry", g["d_ic"], g["o_ic"], g["do"]))
    for s in rep.get("superposed", []):
        print("      align   %s" % s)
    for w in rep.get("warnings", []):
        print("      WARN    %s" % w)
    for o in rep["objects"]:
        print("      object  %s" % o)
    for s in rep["scenes"]:
        print("      scene   %s" % s)
    return 0


def selftest():
    """Build a synthetic pair and check the session, with no input files."""
    cmd = _cmd()
    import tempfile

    ok = True

    def check(c, m):
        nonlocal ok
        if not c:
            print("[selftest] FAIL %s" % m)
            ok = False

    with tempfile.TemporaryDirectory() as d:
        # two little "enzymes" 30 A apart, each with a ligand, and a tunnel
        # running out of each toward the other
        lines, serial = [], 0
        for ch, cx in (("B", -15.0), ("N", 15.0)):
            for i in range(40):
                serial += 1
                x = cx + (i % 4) * 2.0 - 3.0
                y = ((i // 4) % 5) * 2.0 - 4.0
                z = (i // 20) * 8.0 - 4.0     # depth, so a cut plane has something to cut
                lines.append("ATOM  %5d  CA  ALA %1s%4d    %8.3f%8.3f%8.3f  1.00 20.00           C"
                             % (serial, ch, i + 1, x, y, z))
            serial += 1
            lines.append("HETATM%5d  C1  LIG %1s%4d    %8.3f%8.3f%8.3f  1.00 20.00           C"
                         % (serial, ch, 900, cx, 0.0, 0.0))
        pdb = os.path.join(d, "cplx.pdb")
        open(pdb, "w").write("\n".join(lines) + "\nEND\n")

        def tunnel(path, x0, x1, y=0.0):
            rows = []
            for i, t in enumerate(np.linspace(0, 1, 9)):
                rows.append("ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f      %6.2f"
                            % (i + 1, x0 + t * (x1 - x0), y, 0.0, 2.0))
            open(path, "w").write("\n".join(rows) + "\nEND\n")

        t_up, t_dn = os.path.join(d, "up.pdb"), os.path.join(d, "dn.pdb")
        tunnel(t_up, -15.0, -5.0)           # exits at x = -5, pointing +x
        tunnel(t_dn, 15.0, 5.0, y=6.0)      # exits at (5, 6), so P1-P2-P3 is a
                                            # real triangle and the cut plane exists
        out = os.path.join(d, "t_rotated0.pse")

        spec = SceneSpec(pdb, tunnel_up=t_up, tunnel_dn=t_dn,
                         up_chains=["B"], dn_chains=["N"], ligand_resn=["LIG"],
                         label="selftest", zoom="complex")
        rep = build_session(spec, out)

        check(os.path.isfile(out) and os.path.getsize(out) > 0, "a .pse was written")
        check(any(o.startswith("tunnel_up (") and "9 from CAVER" in o
                  for o in rep["objects"]), "upstream tunnel kept as an object")
        check(any(o.startswith("tunnel_dn (") and "9 from CAVER" in o
                  for o in rep["objects"]), "downstream tunnel kept as an object")
        check(any("wraps" in o and "ligand atom" in o for o in rep["objects"]),
              "the drawn channel wraps the channelled ligand (p5c rule)")
        g = rep.get("geometry", {})
        # P2 = (-5, 0, 0), P3 = (5, 6, 0): D_IC = |(10, 6, 0)| = 11.66,
        # O_IC = angle between (-1, 0, 0) and (10, 6, 0) = 149.0 deg
        check(abs(g.get("d_ic", 0) - 11.6619) < 0.01,
              "D_IC is the exit-to-exit distance, got %s" % g.get("d_ic"))
        check(abs(g.get("o_ic", 0) - 149.036) < 0.5,
              "O_IC is the angle at the upstream exit, got %s" % g.get("o_ic"))
        check(rep.get("frame") == "rotated0", "the angle frame was built")

        names = cmd.get_scene_list()
        for s in ("occupancy", "base", "mask", "mouth", "footprint",
                  "render"):
            check(s in names, "scene %r stored (have %s)" % (s, sorted(names)))
        check("classes" not in names,
              "a pose session carries no classification scene by default "
              "(it belongs to the hotspot map)")

        # each scene must show exactly the right objects
        for scene, want, forbid in (
            ("base", {"cplx_il"}, {"mask_up", "footprint", "IC_line"}),
            ("mask", {"mask_up", "mask_dn"}, {"cplx", "footprint"}),
            ("mouth", {"mouth_up", "mouth_dn"}, {"cplx", "mask_up"}),
            ("footprint", {"footprint"}, {"cplx", "mask_up"}),
            ("render", {"cplx_il", "tunnel_up", "tunnel_dn", "IC_line"}, {"mask_up"}),
        ):
            cmd.scene(scene, "recall")
            on = set(cmd.get_names("objects", enabled_only=1))
            check(want <= on, "scene %s shows %s, got %s" % (scene, sorted(want), sorted(on)))
            check(not (forbid & on), "scene %s must not show %s" % (scene, sorted(forbid & on)))

        # the cross-section must carry BOTH a solid cut copy and a faint ghost,
        # otherwise the enzymes look sliced into pieces
        check("render_cs" not in cmd.get_scene_list(),
              "the cross-section scene is no longer built (dropped by the user)")

        # the inspection scenes must be LESS transparent than the render ones,
        # but not solid: a solid surface hides the cofactors
        solid_t = float(cmd.get("transparency", "cplx"))
        illus_t = float(cmd.get("transparency", "cplx_il"))
        check(0.0 < solid_t < illus_t,
              "classes/occupancy use the inspection copy, translucent but less "
              "so than the render copy (%s vs %s)" % (solid_t, illus_t))

        # a tunnel file with a single sphere has no direction
        bad = os.path.join(d, "bad.pdb")
        open(bad, "w").write("ATOM      1  H   FIL T   1       0.000   0.000   0.000      2.00\nEND\n")
        try:
            read_tunnel_spheres(bad)
            check(False, "a one-sphere tunnel should raise")
        except ValueError:
            pass

        # the O_IC arm is a fixed arc length, not one sphere step
        xyz, _ = read_tunnel_spheres(t_up)
        p = arc_point(xyz, 1.0)
        check(abs(float(np.linalg.norm(xyz[-1] - p)) - 1.0) < 1e-6,
              "the O_IC arm is 1.0 A of arc length")

        # THE invariant of the ligand wrap: it changes what is DRAWN and
        # nothing that is MEASURED. If wrapping ever moved an exit, every
        # D_IC and O_IC in the project would silently shift.
        bare = SceneSpec(pdb, tunnel_up=t_up, tunnel_dn=t_dn,
                         up_chains=["B"], dn_chains=["N"], ligand_resn=["LIG"],
                         label="selftest", zoom="complex", ligand_wrap=0.0)
        rep_bare = build_session(bare, os.path.join(d, "t_bare.pse"))
        gb = rep_bare.get("geometry", {})
        for key in ("d_ic", "o_ic"):
            check(abs(float(g.get(key, 0)) - float(gb.get(key, -1))) < 1e-9,
                  "%s is identical with and without the ligand wrap (%s vs %s)"
                  % (key, g.get(key), gb.get(key)))
        check(not any("wraps" in o for o in rep_bare["objects"]),
              "ligand_wrap=0 draws the CAVER chain alone")

    print("[selftest] %s -- session written; tunnels kept as objects at CAVER radii "
          "and the drawn volume wraps the channelled ligand WITHOUT moving D_IC or "
          "O_IC; the angle frame built; scenes stored (occupancy, base, mask, mouth, "
          "footprint, render, base_*) each showing exactly the right "
          "objects, and no classification scene unless asked for; inspection scenes "
          "translucent but less so than the render scenes; no cross-section "
          "scene is built; one-sphere tunnel "
          "refused; fixed-length O_IC arm"
          % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
