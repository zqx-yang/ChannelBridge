"""The overall hotspot map: one session showing every candidate epitope.

This is the *other* half of the figure set. A pose session (
:mod:`channelbridge.render.scenes`) answers "what does this one pose look
like"; this session answers "what did the pose search have to choose from".

Three decisions make it a separate product rather than another scene:

* **Both enzymes are painted the same.** In a pose figure the two enzymes are
  told apart by colour, because which is upstream is the point. Here the
  subject is the residue *classification* -- surface, avoid, shadow, channel
  -- and giving the two enzymes different surface colours would make the same
  class look like two classes. So one class, one colour, on both.
* **Every patch is coloured, in one hue.** The patches of an enzyme are
  ranked, and a rank is an ordering, not a set of unrelated categories, so
  they are shades of a single colour rather than a rainbow. Each is labelled
  where it sits, so the shade does not have to be read as a number.
* **Each patch gets its own binder envelope.** A patch is only a candidate if
  a binder can actually sit on it, so the sphere a binder would occupy is
  drawn at every patch centroid, as its own object. Together with the pose's
  own capsule this shows at a glance which pairs could be bridged.

The geometry is taken from one pose's complex (``geom001`` by default), so
the two enzymes stand in a real relative position rather than side by side.

The pure functions here are importable anywhere; :func:`build` needs PyMOL.
"""

from __future__ import annotations

import itertools
import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .style import RenderStyle, lighten, parse_color

__all__ = [
    "CLASS_ORDER",
    "ClassPalette",
    "HotspotMapSpec",
    "Patch",
    "build",
    "patch_gradient",
    "read_patch_table",
    "read_residue_classes",
]

#: Painted in this order, so a residue in two classes takes the later one.
#: The ordering is deliberate: a channel-wall residue that is also solvent
#: exposed is a channel residue first -- the channel is what the figure is
#: about -- and a hotspot patch overrides everything, because it is drawn as
#: its own object on top.
CLASS_ORDER = (
    "BURIED", "ALLOWED", "IFACE", "LIGAND_CORE", "AVOID", "SHADOW", "CHANNEL",
)


def _default_class_colors() -> Dict[str, str]:
    # The p5d map script's palette, which the user's eye is already trained
    # on, minus the two distinctions this session deliberately does not make:
    # upstream vs downstream, and hydrophobic vs other exposed surface.
    return {
        "BURIED": "wheat",
        "ALLOWED": "cyan",        # solvent-exposed and usable
        "IFACE": "grey30",        # buried by another chain of the assembly
        "LIGAND_CORE": "purple",  # the residues holding the cofactor
        "AVOID": "grey60",
        "SHADOW": "pink",         # inside the cylinder a binder must not block
        "CHANNEL": "magenta",     # the tunnel wall
    }


@dataclass(frozen=True)
class ClassPalette:
    """Colours for the residue classes, the patches and the envelopes."""

    classes: Dict[str, str] = field(default_factory=_default_class_colors)
    #: The hue the patch shades are made from.
    patch_base: str = "red"
    #: Lightest and darkest shade, as a fraction of white mixed in. The first
    #: patch gets ``patch_light``, the last ``patch_dark``.
    patch_light: float = 0.45
    patch_dark: float = 0.0
    #: ``False`` makes rank 1 the palest and the last rank the deepest, which
    #: reads as "increasing shade" down the list. ``True`` puts the strongest
    #: colour on the best-ranked patch instead.
    best_is_darkest: bool = False
    ligand: str = "orange"
    binder: str = "palegreen"
    label: str = "black"
    #: How far the BACKGROUND classes are mixed towards white, so that the
    #: patches are the only saturated thing on the surface. At full
    #: saturation the classification drowns them, which is what the first
    #: version of this session did. Lightening *everything* is worse: it
    #: pushes wheat, grey60 and pink together until BURIED, AVOID and SHADOW
    #: cannot be told apart (0.14 in RGB at 0.45). So only the classes that
    #: mean "nothing special here" are lightened, and the three that are
    #: constraints the user must see -- SHADOW, CHANNEL, LIGAND_CORE -- keep
    #: their full colour.
    class_lighten: float = 0.35
    background_classes: Tuple[str, ...] = ("BURIED", "ALLOWED", "IFACE", "AVOID")
    #: The flat colour the patch-focused scenes paint everything that is not
    #: a patch, so the patches are the only thing competing for attention.
    neutral: str = "grey80"
    #: The channel volume and its mouth, as in the per-pose session.
    channel_volume: str = "magenta"
    mouth: str = "deepteal"
    #: The tunnel sits INSIDE a translucent surface, so it has to be nearly
    #: opaque or it disappears; the binder envelope sits outside and in front
    #: of things worth seeing, so it is faint. Both are p5d's values.
    tunnel_transparency: float = 0.10
    binder_transparency: float = 0.70

    _ALLOWED = ("classes", "patch_base", "patch_light", "patch_dark",
                "best_is_darkest", "ligand", "binder", "label",
                "channel_volume", "mouth", "tunnel_transparency",
                "binder_transparency", "class_lighten", "background_classes",
                "neutral")

    def __post_init__(self):
        for name, value in self.classes.items():
            if name not in CLASS_ORDER:
                raise ValueError(
                    "unknown residue class %r; known: %s"
                    % (name, ", ".join(CLASS_ORDER))
                )
            parse_color(value)
        for name in ("patch_base", "ligand", "binder", "label",
                     "channel_volume", "mouth", "neutral"):
            parse_color(getattr(self, name))
        for name in ("tunnel_transparency", "binder_transparency"):
            v = float(getattr(self, name))
            if not 0.0 <= v <= 1.0:
                raise ValueError("%s must be in [0, 1], got %r" % (name, v))
        bad = [c for c in self.background_classes if c not in CLASS_ORDER]
        if bad:
            raise ValueError(
                "unknown background class(es) %s; known: %s"
                % (", ".join(repr(b) for b in bad), ", ".join(CLASS_ORDER)))
        for name in ("patch_light", "patch_dark", "class_lighten"):
            v = float(getattr(self, name))
            if not 0.0 <= v <= 1.0:
                raise ValueError("%s must be in [0, 1], got %r" % (name, v))

    @classmethod
    def from_dict(cls, raw: Optional[Mapping], where: str = "render.hotspots") -> "ClassPalette":
        raw = dict(raw or {})
        unknown = sorted(set(raw) - set(cls._ALLOWED))
        if unknown:
            raise ValueError(
                "%s: unknown key(s) %s; allowed: %s"
                % (where, ", ".join(repr(u) for u in unknown), ", ".join(sorted(cls._ALLOWED)))
            )
        kw = dict(raw)
        if "classes" in kw:
            merged = _default_class_colors()
            given = dict(kw["classes"] or {})
            bad = sorted(set(given) - set(CLASS_ORDER))
            if bad:
                raise ValueError(
                    "%s: unknown residue class(es) %s; known: %s"
                    % (where, ", ".join(repr(b) for b in bad), ", ".join(CLASS_ORDER))
                )
            merged.update(given)
            kw["classes"] = merged
        # numbers only: `neutral` is a colour and `background_classes` a list,
        # and both were briefly in this loop, where float() rejected them
        for k in ("patch_light", "patch_dark", "tunnel_transparency",
                  "binder_transparency", "class_lighten"):
            if k in kw:
                kw[k] = float(kw[k])
        if "background_classes" in kw:
            v = kw["background_classes"]
            if isinstance(v, str):
                v = [x.strip() for x in v.replace(",", " ").split() if x.strip()]
            kw["background_classes"] = tuple(str(x) for x in v)
        if "best_is_darkest" in kw:
            kw["best_is_darkest"] = bool(kw["best_is_darkest"])
        try:
            return cls(**kw)
        except TypeError as exc:
            raise ValueError("%s: %s" % (where, exc)) from None

    def rgb(self, name: str) -> Tuple[float, float, float]:
        return parse_color(self.classes[name] if name in self.classes else name)

    def final_class_rgb(self, name: str) -> Tuple[float, float, float]:
        """The colour a class is actually drawn in, lightening included."""
        rgb = parse_color(self.classes[name])
        return lighten(rgb, self.class_lighten) if name in self.background_classes else rgb

    def indistinguishable(self, n_patches: int = 6, min_sep: float = 0.30):
        """Pairs of drawn colours too close to tell apart on a shaded surface.

        Reported rather than corrected: the class colours are the ones the
        user's eye is trained on from the p5d map, so a clash is something to
        be told about, not something to silently repaint.
        """
        cols = {c: self.final_class_rgb(c) for c in CLASS_ORDER if c in self.classes}
        for i, shade in enumerate(patch_gradient(n_patches, self.patch_base,
                                                 self.patch_light, self.patch_dark,
                                                 self.best_is_darkest), 1):
            cols["patch_%02d" % i] = shade
        out = []
        for a, b in itertools.combinations(sorted(cols), 2):
            # two shades of the gradient being close is the gradient working,
            # not a clash: they are ordered, adjacent and individually labelled
            if a.startswith("patch_") and b.startswith("patch_"):
                continue
            d = math.dist(cols[a], cols[b])
            if d < float(min_sep):
                out.append((a, b, round(d, 3)))
        return sorted(out, key=lambda r: r[2])


def patch_gradient(n: int, base="red", light: float = 0.62, dark: float = 0.0,
                   best_is_darkest: bool = False) -> List[Tuple[float, float, float]]:
    """``n`` shades of one hue, pale to deep.

    Mixing with white rather than scaling towards black: a dark red and a
    dark grey are hard to tell apart on a shaded surface, while a pale red
    still reads as red. ``n == 1`` gets the deep end -- a lone patch should
    not be drawn in the faintest shade available.
    """
    if int(n) <= 0:
        return []
    rgb = np.asarray(parse_color(base), float)
    light = float(light)
    dark = float(dark)
    if int(n) == 1:
        mix = [dark]
    else:
        mix = list(np.linspace(light, dark, int(n)))
    if best_is_darkest:
        mix = list(reversed(mix))
    out = []
    for m in mix:
        out.append(tuple(float(c + (1.0 - c) * m) for c in rgb))
    return out


# ------------------------------------------------------------------- the tables


def read_residue_classes(path) -> Dict[str, Dict[int, str]]:
    """``<NAME>_residue_classes.tsv`` as ``{class: {resi: resn}}``.

    Keyed by class rather than by residue so that painting is one selection
    per class instead of one per residue: a 300-residue enzyme painted
    residue by residue takes PyMOL minutes and the session megabytes.
    """
    out: Dict[str, Dict[int, str]] = {}
    if not path or not os.path.isfile(path):
        raise FileNotFoundError("no residue-class table at %r" % (path,))
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            i_resi, i_resn, i_cls = hdr.index("resi"), hdr.index("resn"), hdr.index("class")
        except ValueError:
            raise ValueError(
                "%s does not look like a residue-class table (columns: %s)"
                % (path, ", ".join(hdr))
            ) from None
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_resi, i_resn, i_cls):
                continue
            try:
                resi = int(f[i_resi])
            except ValueError:
                continue
            out.setdefault(f[i_cls].strip(), {})[resi] = f[i_resn].strip()
    return out


@dataclass(frozen=True)
class Patch:
    """One candidate epitope."""

    rank: str
    residues: List[int]
    score: Optional[float] = None
    sequence: str = ""

    @property
    def label(self) -> str:
        """``patch_01``: zero-padded so a sorted listing stays in rank order."""
        try:
            return "patch_%02d" % int(self.rank)
        except ValueError:
            return "patch_%s" % self.rank


def read_patch_table(path) -> List[Patch]:
    """``<NAME>_hotspot_patches.tsv`` in rank order."""
    if not path or not os.path.isfile(path):
        raise FileNotFoundError("no hotspot-patch table at %r" % (path,))
    rows: List[Patch] = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            i_rank = hdr.index("rank")
            i_res = hdr.index("hotspot_residues")
        except ValueError:
            raise ValueError(
                "%s does not look like a hotspot-patch table (columns: %s)"
                % (path, ", ".join(hdr))
            ) from None
        i_score = hdr.index("score") if "score" in hdr else None
        i_seq = hdr.index("sequence") if "sequence" in hdr else None
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_rank, i_res):
                continue
            resi = []
            for tok in f[i_res].replace("+", " ").replace(",", " ").split():
                try:
                    resi.append(int(tok))
                except ValueError:
                    continue
            score = None
            if i_score is not None and len(f) > i_score:
                try:
                    score = float(f[i_score])
                except ValueError:
                    score = None
            rows.append(Patch(
                rank=str(f[i_rank]).strip(),
                residues=resi,
                score=score,
                sequence=f[i_seq].strip() if i_seq is not None and len(f) > i_seq else "",
            ))
    rows.sort(key=lambda p: (_rank_number(p.rank), p.rank))
    return rows


def _rank_number(rank) -> float:
    try:
        return float(str(rank).strip())
    except ValueError:
        return float("inf")


# -------------------------------------------------------------------- the spec


@dataclass
class HotspotMapSpec:
    """What the map session is built from."""

    #: Either a pre-posed complex, OR the two prepared enzymes plus their
    #: tunnels, in which case the relative position is computed here (see
    #: :mod:`channelbridge.reference_pose`). The second form is what the
    #: pipeline uses: the map is shown BEFORE the DO sampling, to judge the
    #: epitopes before committing to a search, so there is no pose to borrow.
    complex_pdb: Optional[str] = None
    up_pdb: Optional[str] = None
    dn_pdb: Optional[str] = None
    tunnel_up: Optional[str] = None
    tunnel_dn: Optional[str] = None
    #: Residue names CAVER excluded, so burial is not counted against atoms it
    #: never searched, and the clash test is not run against cofactors.
    exclude_resn: Sequence[str] = ()
    #: The pose search's own binder gap: the shortest binder worth drawing.
    binder_length: float = 6.0
    binder_step: float = 2.0
    d_min: float = 6.0
    d_max: float = 36.0
    d_step: float = 1.0
    clash_cut: float = 3.0
    oic_arm: float = 1.0
    #: How much room a label needs, in A, along the screen's x and y. A label
    #: is a line of text -- far wider than it is tall -- so these differ by a
    #: lot; a single circular distance either leaves wide labels overlapping
    #: or throws short ones much too far from the patch they name.
    label_sep_x: float = 34.0
    label_sep_y: float = 11.0
    up_name: str = "up"
    dn_name: str = "dn"
    up_chains: Sequence[str] = ()
    dn_chains: Sequence[str] = ()
    classes_up: Optional[str] = None
    classes_dn: Optional[str] = None
    patches_up: Optional[str] = None
    patches_dn: Optional[str] = None
    ligand_resn: Sequence[str] = ()
    binder_radius: float = 11.0
    palette: ClassPalette = field(default_factory=ClassPalette)
    #: How see-through the surfaces are. The same "inspect" value the pose
    #: session's classes/occupancy scenes use, so the two sessions of one
    #: figure set look alike and one override changes both.
    transparency: float = field(
        default_factory=lambda: RenderStyle().transparency_for("inspect"))
    #: Chains to paint from each table. The class tables are written for one
    #: protomer; the other protomers of a homo-oligomer are identical, so the
    #: same residue numbers are painted on all of them.
    label: str = ""

    def chains(self, which: str) -> List[str]:
        return [c for c in (self.up_chains if which == "up" else self.dn_chains) if c]

    def name(self, which: str) -> str:
        return self.up_name if which == "up" else self.dn_name


def patch_centroids(coords_by_resi: Mapping[int, Sequence[float]],
                    patches: Sequence[Patch]) -> Dict[str, np.ndarray]:
    """Centroid of each patch's residues that are present in the structure.

    A patch whose residues are all missing is left out rather than collapsed
    to the origin -- a green ball at (0, 0, 0) is worse than no ball.
    """
    out = {}
    for p in patches:
        pts = [coords_by_resi[r] for r in p.residues if r in coords_by_resi]
        if not pts:
            continue
        out[p.label] = np.asarray(pts, float).mean(axis=0)
    return out


# -------------------------------------------------------------------- the build


def _cmd():
    try:
        from pymol import cmd
    except ImportError as exc:  # pragma: no cover - only outside PyMOL
        raise ImportError(
            "the hotspot map is built inside PyMOL; run it with\n"
            "  pymol -cq %s -- <key=value> ..." % (__file__,)
        ) from exc
    return cmd


def _place(spec, report):
    """The reference placement, or ``None`` when a posed complex was given."""
    from ..io.pdb import read_pdb
    from ..reference_pose import place_at_oic180
    from ..tunnel import Tunnel

    if not (spec.up_pdb and spec.dn_pdb and spec.tunnel_up and spec.tunnel_dn):
        raise ValueError(
            "to compute a placement the map needs up_pdb, dn_pdb, tunnel_up and "
            "tunnel_dn; give complex_pdb instead to use a pose that already exists"
        )
    ex = list(spec.exclude_resn)
    up = read_pdb(spec.up_pdb).exclude_resnames(ex).xyz if ex else read_pdb(spec.up_pdb).xyz
    dn = read_pdb(spec.dn_pdb).exclude_resnames(ex).xyz if ex else read_pdb(spec.dn_pdb).xyz
    tu = Tunnel.from_caver_cluster(spec.tunnel_up)
    td = Tunnel.from_caver_cluster(spec.tunnel_dn)
    report("placing at O_IC 180 deg, smallest clash-free D_IC from %.0f A" % spec.d_min)
    placement = place_at_oic180(
        up, tu.exit, tu.axis_point(spec.oic_arm),
        dn, td.exit, td.axis_point(spec.oic_arm),
        d_min=spec.d_min, d_max=spec.d_max, d_step=spec.d_step,
        clash_cut=spec.clash_cut, report=None,
    )
    if placement is None:
        raise ValueError(
            "no clash-free placement at O_IC 180 between %.0f and %.0f A; the two "
            "channels cannot be brought face to face at any sampled distance"
            % (spec.d_min, spec.d_max)
        )
    report("  " + placement.describe())
    return placement, tu, td


def _pymol_matrix(rotation, translation):
    """PyMOL's 16-element TTT matrix for a rotation and a translation."""
    m = np.eye(4)
    m[:3, :3] = np.asarray(rotation, float)
    m[:3, 3] = np.asarray(translation, float)
    return [float(x) for x in m.flatten()]


def _tunnel_object(cmd, name, path, rgb, transparency, transform=None):
    """A CAVER sphere chain as its own object, at CAVER's own radii."""
    from .scenes import read_tunnel_spheres

    xyz, rad = read_tunnel_spheres(path)
    if transform is not None:
        rot, trans = transform
        xyz = xyz @ np.asarray(rot, float).T + np.asarray(trans, float)
    for i, (p, r) in enumerate(zip(xyz, rad), 1):
        cmd.pseudoatom(name, pos=[float(x) for x in p], resn="TUN", resi=str(i),
                       name="T%d" % i, elem="C", b=float(r))
    cmd.alter(name, "vdw = b")
    cmd.rebuild(name)
    cmd.show_as("spheres", name)
    cmd.set_color("c_" + name.lower(), [float(c) for c in rgb])
    cmd.color("c_" + name.lower(), name)
    cmd.set("sphere_transparency", transparency, name)
    return xyz, rad


def build(spec: HotspotMapSpec, out_pse: str, report=print) -> Dict:
    """Build the map session and save it. Returns what was drawn.

    Two routes in: a complex that already exists, or the two enzymes and their
    tunnels, which are placed here. The second is the pipeline's -- the map is
    shown before the DO sampling, so there is no pose yet to borrow.
    """
    from ..reference_pose import (
        binder_blobs,
        clear_view_about_axis,
        patch_outward_normal,
    )
    from .frame import Frame, pymol_view, spread_in_view
    from .scenes import _apply, _ca_by_resi, _pseudo_spheres
    from .style import INSPECT_SETTINGS

    cmd = _cmd()
    cmd.reinitialize()
    cmd.set("auto_zoom", 0)
    pal = spec.palette
    info: Dict[str, object] = {"objects": [], "patches": {}, "envelopes": {},
                               "placement": None}

    # ---- the structures, in a relative position
    placement = None
    frame = None
    up_obj, dn_obj = spec.up_name, spec.dn_name
    if spec.complex_pdb:
        cmd.load(spec.complex_pdb, "cplx")
        info["complex"] = os.path.basename(spec.complex_pdb)
        for which, obj in (("up", up_obj), ("dn", dn_obj)):
            chains = spec.chains(which)
            if not chains:
                raise ValueError("a posed complex needs %s_chains" % which)
            cmd.create(obj, "cplx and chain %s" % "+".join(chains))
        cmd.delete("cplx")
    else:
        placement, _tun_up, _tun_dn = _place(spec, report)
        info["placement"] = placement.describe()
        info["complex"] = "%s + %s, placed here" % (os.path.basename(spec.up_pdb),
                                                    os.path.basename(spec.dn_pdb))
        cmd.load(spec.up_pdb, up_obj)
        cmd.load(spec.dn_pdb, dn_obj)
        cmd.transform_object(dn_obj, _pymol_matrix(placement.rotation,
                                                   placement.translation))

    cmd.hide("everything", "all")
    for obj in (up_obj, dn_obj):
        cmd.show("surface", "%s and polymer" % obj)

    # ---- the classification, painted IN PLACE on each enzyme
    #
    # Not as separate surface objects: a surface built from ten residues alone
    # renders as a free-floating blob that z-fights with the parent surface,
    # which is why the first version of this session was unreadable. Colouring
    # the parent's own atoms keeps one continuous surface.
    #
    # The same colour for the same class on BOTH enzymes: the subject here is
    # the classification, and two surface colours for one class would read as
    # two classes.
    painted: Dict[str, int] = {}
    for path, obj in ((spec.classes_up, up_obj), (spec.classes_dn, dn_obj)):
        if not path:
            continue
        table = read_residue_classes(path)
        for cls in CLASS_ORDER:
            resi = sorted(table.get(cls, {}))
            colour = pal.classes.get(cls)
            if not resi or colour is None:
                continue
            key = "c_cls_%s" % cls.lower()
            rgb = parse_color(colour)
            if cls in pal.background_classes:
                rgb = lighten(rgb, pal.class_lighten)
            cmd.set_color(key, [float(c) for c in rgb])
            cmd.color(key, "%s and resi %s and polymer"
                      % (obj, "+".join(str(r) for r in resi)))
            painted[cls] = painted.get(cls, 0) + len(resi)
    for cls in CLASS_ORDER:
        if cls in painted:
            info["objects"].append("%-12s %-8s %d residue(s) over both enzymes"
                                   % (cls, pal.classes.get(cls), painted[cls]))

    n_up = len(info["patches"].get(spec.up_name, [])) if info.get("patches") else 6
    for a, b, d in pal.indistinguishable(max(6, n_up)):
        report("WARNING: %s and %s are %.2f apart in colour and may not be "
               "tellable apart on a shaded surface" % (a, b, d))

    if spec.ligand_resn:
        sel = "(%s or %s) and resn %s" % (up_obj, dn_obj, "+".join(spec.ligand_resn))
        cmd.set_color("c_lig", [float(c) for c in parse_color(pal.ligand)])
        cmd.color("c_lig", sel)
        cmd.show("sticks", sel)
        cmd.set("stick_radius", 0.45)
        cmd.set("stick_ball", 0)

    for obj in (up_obj, dn_obj):
        cmd.set("transparency", spec.transparency, "%s and polymer" % obj)
        cmd.set("transparency", spec.transparency, obj)
    cmd.bg_color("white")
    # INSPECT_SETTINGS, not RAY_SETTINGS: see the note on both in style.py.
    # The outline mode is what made every previous version of this session
    # come out dark.
    _apply(cmd, INSPECT_SETTINGS)
    info["objects"].append("%s + %s (surfaces at transparency %.2f, ligands as sticks)"
                           % (up_obj, dn_obj, spec.transparency))

    # ---- the tunnels, as objects. Nearly opaque: they are inside a
    # translucent surface and a translucent tunnel inside it disappears.
    tunnel_objs: List[str] = []
    mouth_objs: List[str] = []
    exits: Dict[str, np.ndarray] = {}
    axis_vector = None
    for which, path in (("up", spec.tunnel_up), ("dn", spec.tunnel_dn)):
        if not path:
            continue
        name = "%s_tunnel" % spec.name(which)
        tf = None
        if which == "dn" and placement is not None:
            tf = (placement.rotation, placement.translation)
        xyz, rad = _tunnel_object(cmd, name, path, parse_color(pal.channel_volume),
                                  pal.tunnel_transparency, tf)
        tunnel_objs.append(name)
        mouth = "%s_mouth" % spec.name(which)
        _pseudo_spheres(cmd, mouth, xyz[-1][None, :],
                        [float(np.clip(rad[-1], 1.0, 1.8))],
                        tuple(float(c) for c in parse_color(pal.mouth)), 0.0)
        mouth_objs.append(mouth)
        exits[which] = np.asarray(xyz[-1], float)
        info["objects"].append("%s (%d spheres) + %s" % (name, len(xyz), mouth))
    if "up" in exits and "dn" in exits:
        info["objects"].append(
            "channel axis %.1f A between the two mouths"
            % float(np.linalg.norm(exits["dn"] - exits["up"])))

    # ---- the camera: the channel axis across the screen, upstream left, and
    # the roll chosen so the two enzymes overlap as little as possible. The
    # default view puts one enzyme in front of the other, which is what made
    # the first version of this session unreadable.
    up_coords = cmd.get_coords("%s and polymer" % up_obj)
    dn_coords = cmd.get_coords("%s and polymer" % dn_obj)
    # Across the screen goes the CENTROID offset, not the channel axis. The
    # two mouths both sit at the interface, so laying the channel axis
    # horizontal leaves the bulk of each enzyme hanging off to one side and
    # the two still overlapping -- which is what the first camera did.
    axis_vector = dn_coords.mean(axis=0) - up_coords.mean(axis=0)
    if float(np.linalg.norm(axis_vector)) < 1e-6:
        # no tunnels given and the two enzymes sit on top of each other:
        # there is no axis to lay across the screen, so let PyMOL choose
        cam_y = np.array([0.0, 1.0, 0.0])
        cmd.orient("%s or %s" % (up_obj, dn_obj))
        report("WARNING: no channel axis and no separation between the two "
               "enzymes; falling back to PyMOL's own orientation")
        info["objects"].append("camera: PyMOL orient (no channel axis available)")
    else:
        cam_x, cam_y, cam_z = clear_view_about_axis(up_coords, dn_coords, axis_vector)
        frame = Frame(axes=np.stack([cam_x, cam_y, cam_z]),
                      centre=np.vstack([up_coords, dn_coords]).mean(axis=0))
        cmd.set_view(pymol_view(frame))
        info["objects"].append(
            "camera: the two enzymes across the screen, upstream left, roll "
            "chosen for least projected overlap")

    # ---- a neutral twin of each enzyme.
    #
    # The classification is a seven-colour map; the patches are what the user
    # is here to judge. Shown together, the rainbow wins and the patches are
    # lost in it -- which is what the first two versions of this session did.
    # So the patch-focused scenes get a copy whose surface is one flat grey,
    # with only the patches coloured on it, while "classes" and "map" keep
    # the full classification. Two objects rather than one, because a PyMOL
    # scene does not reliably store per-atom colour.
    plain = {}
    for which, obj in (("up", up_obj), ("dn", dn_obj)):
        twin = obj + "_plain"
        cmd.create(twin, obj)
        cmd.hide("everything", twin)
        cmd.show("surface", "%s and polymer" % twin)
        cmd.set_color("c_plain", [float(c) for c in parse_color(pal.neutral)])
        cmd.color("c_plain", "%s and polymer" % twin)
        if spec.ligand_resn:
            cmd.color("c_lig", "%s and resn %s" % (twin, "+".join(spec.ligand_resn)))
            cmd.show("sticks", "%s and resn %s" % (twin, "+".join(spec.ligand_resn)))
        cmd.set("transparency", spec.transparency, "%s and polymer" % twin)
        cmd.set("transparency", spec.transparency, twin)
        plain[which] = twin
    info["objects"].append("%s + %s (neutral %s twins, for the patch scenes)"
                           % (plain["up"], plain["dn"], pal.neutral))

    # ---- the patches, in place, in one hue, each labelled, each with the
    # envelope a minimum-length binder on it would occupy
    patch_objs: List[str] = []
    env_objs: List[str] = []
    label_objs: List[str] = []
    pending_labels: List[tuple] = []
    for which, path, obj in (("up", spec.patches_up, up_obj),
                             ("dn", spec.patches_dn, dn_obj)):
        if not path:
            continue
        patches = read_patch_table(path)
        shades = patch_gradient(len(patches), pal.patch_base, pal.patch_light,
                                pal.patch_dark, pal.best_is_darkest)
        coords = _ca_by_resi(cmd, obj)
        centroids = patch_centroids(coords, patches)
        surface = cmd.get_coords("%s and polymer" % obj)
        drawn = []
        for patch, shade in zip(patches, shades):
            if not patch.residues:
                continue
            tag = "%s_%s" % (spec.name(which), patch.label)
            sel = "%s and resi %s and polymer" % (obj, "+".join(str(r) for r in patch.residues))
            if cmd.count_atoms(sel) == 0:
                continue
            key = "c_%s" % tag.lower()
            cmd.set_color(key, [float(c) for c in shade])
            # ONLY on the neutral twin. Painting the patches onto the
            # classified object as well destroys both readings at once: the
            # classification gets red holes in it and the patches compete
            # with six other colours. The classified objects stay pure.
            cmd.color(key, sel.replace(obj, plain[which], 1))
            # a named selection, not an object: the colour is already on the
            # parent surface, and this is so the patch can be picked later
            cmd.select(tag, sel)
            patch_objs.append(tag)
            drawn.append(tag)
            centre = centroids.get(patch.label)
            if centre is None:
                continue
            # outward along the LOCAL surface normal, so the envelope shows
            # how a binder on this epitope would extend regardless of where
            # the partner enzyme happens to be
            normal = patch_outward_normal(surface, centre)
            lab = tag + "_label"
            # Short text, lifted off the surface along the patch's OWN normal
            # and clear of its envelope. Anchoring at the centroid piles the
            # labels of neighbouring patches into one illegible stack, and
            # lifting them all along a single axis just moves the stack.
            # Placed after the loop: the anchors of ALL patches are needed
            # at once to push them apart in the camera plane (see below).
            pending_labels.append(
                (lab, np.asarray(centre, float) + normal * (spec.binder_length + 8.0),
                 "%s %s" % (spec.name(which), patch.label.replace("patch_", "p"))))
            pts, radii = binder_blobs(centre, normal, spec.binder_length,
                                      spec.binder_step, spec.binder_radius)
            env = tag + "_binder"
            _pseudo_spheres(cmd, env, pts, radii,
                            tuple(float(c) for c in parse_color(pal.binder)),
                            pal.binder_transparency)
            env_objs.append(env)
            info["envelopes"][tag] = {"centroid": [round(float(x), 2) for x in centre],
                                      "blobs": len(pts)}
        info["patches"][spec.name(which)] = drawn
        report("%s: %d patch(es), shades %s"
               % (spec.name(which), len(drawn),
                  " -> ".join("#%02x%02x%02x" % tuple(round(float(c) * 255) for c in s)
                              for s in shades)))
    # ---- the labels, pushed apart where the collision actually happens
    if pending_labels:
        anchors = np.array([a for _, a, _ in pending_labels])
        if frame is not None:
            anchors = anchors + spread_in_view(
                anchors, cam_x, cam_y, min_sep=(spec.label_sep_x, spec.label_sep_y))
        cmd.set_color("c_lab", [float(c) for c in parse_color(pal.label)])
        for (lab, _, text), pos in zip(pending_labels, anchors):
            cmd.pseudoatom(lab, pos=[float(v) for v in pos], label=text)
            cmd.set("label_size", 18, lab)
            cmd.set("label_color", "c_lab", lab)
            label_objs.append(lab)
        report("%d label(s) placed, separated in the camera plane" % len(label_objs))
    cmd.deselect()

    # ---- the classification scene gets the two enzymes PULLED APART.
    #
    # They are painted alike by design and the reference placement puts them
    # 3.4 A apart, so in the shared frame they fuse into one mass and no
    # class can be traced to an enzyme. Nothing in a classification depends
    # on the relative pose, so this scene -- and only this scene -- separates
    # them. The scene message says so, so it cannot be mistaken for a pose.
    split_objs = []
    if frame is not None:
        gap = 0.0
        for obj in (up_obj, dn_obj):
            proj = cmd.get_coords("%s and polymer" % obj) @ cam_x
            gap = max(gap, float(proj.max() - proj.min()))
        for obj, sign in ((up_obj, -1.0), (dn_obj, 1.0)):
            twin = obj + "_split"
            cmd.create(twin, obj)
            shift = cam_x * sign * gap * 0.32
            m = np.eye(4)
            m[:3, 3] = shift
            cmd.transform_object(twin, [float(v) for v in m.flatten()])
            split_objs.append(twin)
        info["objects"].append(
            "%s (classification, pulled %.0f A apart so the two enzymes can be "
            "told apart; NOT a pose)" % (" + ".join(split_objs), gap * 0.64))

    # ---- scenes
    cmd.zoom("%s or %s" % (up_obj, dn_obj), buffer=8, complete=1)
    cmd.clip("slab", 3000)   # zoom keeps the rotation set above and fits the pair
    proteins = [up_obj, dn_obj]
    neutral = [plain["up"], plain["dn"]]
    note = placement.describe() if placement is not None else ""
    scenes = {
        # the whole story, on the neutral surface so the patches lead
        "map": neutral + tunnel_objs + mouth_objs + label_objs + env_objs,
        # the classification alone -- no patches over it, and the two enzymes
        # separated so each one's classes can be read
        "classes": (split_objs or proteins),
        # the epitopes on a NEUTRAL surface, labelled: the scene for judging
        # which patches to carry into the sampling, with nothing competing
        "patches": neutral + tunnel_objs + mouth_objs + label_objs,
        # can a minimum binder sit on each of them, again on neutral
        "envelopes": neutral + env_objs + label_objs,
    }
    for key, objs in scenes.items():
        cmd.disable("all")
        for o in objs:
            cmd.enable(o)
        cmd.scene(key, "store", message=note if key == "map" else None)
        info["objects"].append("scene %-10s %d object(s)" % (key, len(objs)))
    cmd.disable("all")
    for o in scenes["map"]:
        cmd.enable(o)
    cmd.scene("map", "recall")
    out_dir = os.path.dirname(os.path.abspath(out_pse))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    cmd.save(out_pse)
    info["pse"] = out_pse
    report("saved %s" % out_pse)
    return info


# --------------------------------------------------------------------- entry


def _parse_argv(argv) -> Dict[str, str]:
    out = {}
    for a in argv:
        if a == "--selftest":
            out["selftest"] = "1"
        elif "=" in a:
            k, v = a.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _chains(value) -> List[str]:
    return [c for c in str(value or "").replace(",", " ").split() if c]


def main(argv=None) -> int:
    import sys

    args = _parse_argv(list(sys.argv[1:] if argv is None else argv))
    if args.get("selftest"):
        return selftest()
    missing = [k for k in ("out",) if k not in args]
    if "complex" not in args and not all(k in args for k in ("up_pdb", "dn_pdb")):
        missing.append("complex (or up_pdb + dn_pdb + tunnel_up + tunnel_dn)")
    if missing:
        print("missing argument(s): %s\n"
              "usage: pymol -cq hotspots.py -- complex=<geom001_complex.pdb> "
              "out=<map.pse> up_chains=B,C,D,E dn_chains=N "
              "classes_up=<GDH_residue_classes.tsv> ..." % ", ".join(missing))
        return 2
    spec = HotspotMapSpec(
        complex_pdb=args.get("complex"),
        up_pdb=args.get("up_pdb"), dn_pdb=args.get("dn_pdb"),
        tunnel_up=args.get("tunnel_up"), tunnel_dn=args.get("tunnel_dn"),
        exclude_resn=_chains(args.get("exclude")),
        binder_length=float(args.get("binder_length", 6.0)),
        binder_step=float(args.get("binder_step", 2.0)),
        d_min=float(args.get("d_min", 6.0)), d_max=float(args.get("d_max", 36.0)),
        d_step=float(args.get("d_step", 1.0)),
        clash_cut=float(args.get("clash_cut", 3.0)),
        oic_arm=float(args.get("oic_arm", 1.0)),
        up_name=args.get("up_name", "up"),
        dn_name=args.get("dn_name", "dn"),
        up_chains=_chains(args.get("up_chains")),
        dn_chains=_chains(args.get("dn_chains")),
        classes_up=args.get("classes_up"),
        classes_dn=args.get("classes_dn"),
        patches_up=args.get("patches_up"),
        patches_dn=args.get("patches_dn"),
        ligand_resn=_chains(args.get("ligands")),
        binder_radius=float(args.get("binder_radius", 11.0)),
        palette=ClassPalette(
            best_is_darkest=str(args.get("best_is_darkest", "")).lower()
            in ("1", "true", "yes"),
        ),
        label=args.get("label", ""),
    )
    info = build(spec, args["out"])
    for line in info["objects"]:
        print("  " + str(line))
    return 0


def _synthetic(d, chain, x0, resn="ALA"):
    """A small block of CA atoms, as a stand-in for an enzyme."""
    lines, n = [], 1
    for resi in range(1, 31):
        for k in range(3):
            lines.append(
                "ATOM  %5d  CA  %s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C"
                % (n, resn, chain, resi, x0 + (resi % 10) * 2.0,
                   (resi // 10) * 3.0 + k * 1.5, (resi % 4) * 2.0))
            n += 1
    path = os.path.join(d, "%s.pdb" % chain)
    open(path, "w").write("\n".join(lines) + "\nEND\n")
    return path


def _synthetic_tunnel(d, name, start, direction, n=6):
    """A CAVER-format sphere chain running out of a block."""
    start = np.asarray(start, float)
    direction = np.asarray(direction, float)
    direction = direction / np.linalg.norm(direction)
    lines = []
    for i in range(n):
        p = start + direction * (i * 1.5)
        lines.append("ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f      %4.2f"
                     % (i + 1, i + 1, p[0], p[1], p[2], 2.0))
    path = os.path.join(d, "%s.pdb" % name)
    open(path, "w").write("\n".join(lines) + "\nEND\n")
    return path


def selftest() -> int:
    """Build a map from synthetic structures, both routes, inside PyMOL."""
    import tempfile

    cmd = _cmd()
    ok = [True]

    def check(cond, what):
        print(("  ok   " if cond else "  FAIL ") + what)
        if not cond:
            ok[0] = False

    with tempfile.TemporaryDirectory() as d:
        up_pdb = _synthetic(d, "A", 0.0)
        dn_pdb = _synthetic(d, "B", 0.0)
        # the upstream channel leaves the +x face, the downstream one its own
        t_up = _synthetic_tunnel(d, "tun_up", [14.0, 0.0, 2.0], [1.0, 0.0, 0.0])
        t_dn = _synthetic_tunnel(d, "tun_dn", [4.0, 0.0, 2.0], [-1.0, 0.0, 0.0])
        cls = os.path.join(d, "A_classes.tsv")
        open(cls, "w").write(
            "chain\tresi\tresn\tss\tclass\n"
            + "".join("A\t%d\tALA\tL\t%s\n" % (r, c) for r, c in
                      [(1, "AVOID"), (2, "SHADOW"), (3, "CHANNEL"), (4, "ALLOWED")]))
        pt = os.path.join(d, "A_patches.tsv")
        open(pt, "w").write(
            "rank\tscore\thotspot_residues\tsequence\n"
            "1\t4.5\t5+6+7\tALA ALA ALA\n2\t3.9\t10+11\tALA ALA\n3\t3.1\t15+16\tALA ALA\n")

        print("-- route 1: the two enzymes, placed here")
        spec = HotspotMapSpec(
            up_pdb=up_pdb, dn_pdb=dn_pdb, tunnel_up=t_up, tunnel_dn=t_dn,
            up_name="E1", dn_name="E2", classes_up=cls, patches_up=pt,
            binder_radius=4.0, binder_length=6.0, d_min=4.0, d_max=40.0,
            clash_cut=2.0)
        out = os.path.join(d, "map.pse")
        info = build(spec, out)
        check(os.path.isfile(out), "the session is written")
        check(info["placement"] and "O_IC 180" in str(info["placement"]),
              "the placement is computed, not borrowed: %s" % info["placement"])
        objs = set(cmd.get_names("objects"))
        check({"E1", "E2"} <= objs, "one object per enzyme (%s)" % sorted(objs)[:4])
        check("E1_tunnel" in objs and "E2_tunnel" in objs,
              "both tunnels are objects")
        check("E1_mouth" in objs, "the mouth is its own object")

        # the patches are coloured ON the enzyme, not built as separate
        # surfaces -- the bug that made the first version unreadable
        check("E1_patch_01" not in objs,
              "a patch is a selection, not a floating surface object")
        check("E1_patch_01" in set(cmd.get_names("selections")),
              "...but it is still selectable by name")
        check("E1_patch_01_binder" in objs, "every patch gets an envelope")
        check(len(info["envelopes"]) == 3, "three envelopes (%d)" % len(info["envelopes"]))
        blobs = {k: v["blobs"] for k, v in info["envelopes"].items()}
        check(set(blobs.values()) == {4},
              "a 6 A binder at 2 A steps is 4 blobs everywhere (%s)" % blobs)

        for key in ("map", "classes", "patches", "envelopes"):
            check(key in (cmd.get_scene_list() or []), "scene %s stored" % key)
        cmd.scene("classes", "recall")
        on = set(cmd.get_names("objects", enabled_only=1))
        check(not any(o.endswith("_binder") for o in on),
              "the classes scene drops the envelopes (%s)" % sorted(on))
        check(all(o.endswith("_split") for o in on),
              "the classes scene shows the pulled-apart pair (%s)" % sorted(on))
        # the tunnels were built on the un-split objects, so showing them here
        # would draw each channel offset from the enzyme it belongs to
        check("E1_tunnel" not in on,
              "...and drops the tunnels, which no longer line up with it")
        cmd.scene("envelopes", "recall")
        on = set(cmd.get_names("objects", enabled_only=1))
        check(any(o.endswith("_binder") for o in on), "the envelopes scene shows envelopes")

        t = float(cmd.get("transparency", "E1"))
        check(0.0 < t < 1.0, "the enzymes are translucent, so ligands show (%s)" % t)
        tt = float(cmd.get("sphere_transparency", "E1_tunnel"))
        check(tt < 0.3, "the tunnel is nearly opaque inside it (%s)" % tt)
        bt = float(cmd.get("sphere_transparency", "E1_patch_01_binder"))
        check(bt > 0.5, "the envelope is faint (%s)" % bt)

        print("-- route 2: a complex that already exists")
        both = os.path.join(d, "both.pdb")
        open(both, "w").write(
            open(up_pdb).read().replace("END\n", "")
            + open(dn_pdb).read().replace(" B", " B"))
        spec2 = HotspotMapSpec(complex_pdb=both, up_chains=["A"], dn_chains=["B"],
                               up_name="E1", dn_name="E2", classes_up=cls,
                               patches_up=pt, binder_radius=4.0)
        info2 = build(spec2, os.path.join(d, "map2.pse"))
        check(info2["placement"] is None, "a given complex is not re-placed")
        check(len(info2["patches"].get("E1", [])) == 3, "patches still found")

        print("-- missing inputs are an error, not a half-drawn figure")
        try:
            build(HotspotMapSpec(up_pdb=up_pdb), os.path.join(d, "never.pse"))
            check(False, "should have raised")
        except ValueError as exc:
            check("tunnel_up" in str(exc), "...and names what is missing: %s" % exc)

    print("SELFTEST " + ("PASSED" if ok[0] else "FAILED"))
    return 0 if ok[0] else 1
