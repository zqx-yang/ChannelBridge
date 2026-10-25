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

from .style import lighten, parse_color

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
    #: The designed binder's house colour (user, standing).
    binder: str = "#0065bd"
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
class MapObjects:
    """The named object groups a map session is built from.

    Eight bare ``List[str]`` locals threaded through a 400-line function is
    what this replaces. The groups are not interchangeable -- a scene is
    defined by which of THESE it shows -- so giving them a type makes the
    scene definitions read as what they mean, and makes it impossible to
    append a tunnel to the ligand list without noticing.
    """

    #: One translucent whole-enzyme surface per enzyme; everything else is
    #: solid and sits on top of it.
    base: List[str] = field(default_factory=list)
    #: One per residue class per enzyme. Together they partition the surface.
    classes: List[str] = field(default_factory=list)
    patches: List[str] = field(default_factory=list)
    tunnels: List[str] = field(default_factory=list)
    mouths: List[str] = field(default_factory=list)
    ligands: List[str] = field(default_factory=list)
    #: Minimum-binder envelopes, one per patch.
    envelopes: List[str] = field(default_factory=list)
    #: Patch names. Built, and off unless asked for.
    labels: List[str] = field(default_factory=list)
    #: The iMARS annotation. Built, and off: the map is about the epitopes.
    annotation: List[str] = field(default_factory=list)

    @property
    def structure(self) -> List[str]:
        """Everything that is always shown: the surface and what is in it."""
        return self.base + self.tunnels + self.mouths + self.ligands

    def scene(self, *, classes: bool = True, patches: bool = True,
              envelopes: bool = False, labels: bool = False) -> List[str]:
        """The object list for one scene, by what it shows."""
        out = list(self.structure)
        if classes:
            out += self.classes
        if patches:
            out += self.patches
        if envelopes:
            out += self.envelopes
        if labels:
            out += self.labels
        return out


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
    # A label is a line of text, so it needs far more room across than down.
    # These are the gaps between label ANCHORS, in Angstrom, and they carry a
    # margin on top of the text's own width: the session renders in
    # perspective (the view's field-of-view slot is -20), so a label nearer
    # the camera prints larger and shifted from where an orthographic
    # projection puts it, and separation measured on that projection comes up
    # a little short. At 34 the two most crowded FNR labels still overprinted.
    label_sep_x: float = 46.0
    label_sep_y: float = 13.0
    #: Keep labels within the drawn pair plus this margin, in Angstrom, so a
    #: label pushed clear of its neighbours is not pushed off the picture.
    label_margin: float = 10.0
    #: How far in front of the nearest atom a lifted label sits, in Angstrom.
    label_lift: float = 2.0
    #: Probe radius for the patch surfaces. PyMOL's default is 1.4; a little
    #: larger lifts a patch clear of the coincident class surface beneath it.
    patch_probe: float = 1.7
    #: Patch labels are BUILT but start switched off (user, 2026-10-25).
    #: Ten names on two crowded faces could not be placed without some of
    #: them landing away from the patch they name, and every patch is its own
    #: object with its own name, so the object list is the legend. Switch
    #: them on per scene if you want them in a figure.
    show_labels: bool = False
    #: What to call each enzyme. Every object in the session is named from
    #: these -- `<up_name>_patch_01`, `<up_name>_patch_01_binder`,
    #: `<up_name>_channel`, `<up_name>_tunnel` -- so giving the real names
    #: makes the object list readable without a legend. The defaults say
    #: which side each is rather than pretending to a name we do not have.
    up_name: str = "Eup"
    dn_name: str = "Edn"
    up_chains: Sequence[str] = ()
    dn_chains: Sequence[str] = ()
    classes_up: Optional[str] = None
    classes_dn: Optional[str] = None
    patches_up: Optional[str] = None
    patches_dn: Optional[str] = None
    ligand_resn: Sequence[str] = ()
    binder_radius: float = 11.0
    palette: ClassPalette = field(default_factory=ClassPalette)
    #: How see-through the surfaces are. **Solid by default**, and that is
    #: deliberate: this session exists to be toggled. You switch a patch
    #: object off and on to see exactly which bit of surface it covers, and
    #: through a translucent surface you cannot tell whether what you are
    #: looking at is the patch in front or the class behind it.
    #:
    #: They were translucent for a while, to keep the cofactors visible
    #: inside a solid surface. That reason is gone: every class, patch,
    #: ligand, mouth and tunnel is now its OWN object, so the way to see
    #: inside is to hide what is in front, not to see through it.
    transparency: float = 0.0
    #: The whole-enzyme base, the one translucent thing in the session. Every
    #: element drawn on top of it is solid, so switching one on reads as
    #: "this patch, here, on this surface".
    base_transparency: float = 0.55
    #: Chains to paint from each table. The class tables are written for one
    #: protomer; the other protomers of a homo-oligomer are identical, so the
    #: same residue numbers are painted on all of them.
    label: str = ""

    def chains(self, which: str) -> List[str]:
        return [c for c in (self.up_chains if which == "up" else self.dn_chains) if c]

    def name(self, which: str) -> str:
        return self.up_name if which == "up" else self.dn_name


def _protomer_transforms(cmd, obj, primary, chains, report):
    """Rigid transforms taking the primary protomer onto each of the others.

    The upstream enzyme is an oligomer and the downstream one, as prepared,
    is a single chain -- so the assembly's extra downstream copies do not
    exist in any file and have to be generated. Each is the SAME placement,
    carried onto another protomer by that protomer's own symmetry operator,
    which is what the pose search does when it writes an occupancy model.

    Fitted on matched CA atoms rather than assumed to be a crystal operator,
    because a prepared structure may have been renumbered or trimmed.
    """
    out = []
    def ca(ch):
        m = cmd.get_model("%s and chain %s and polymer and name CA" % (obj, ch))
        return {a.resi: np.asarray(a.coord, float) for a in m.atom}
    ref = ca(primary)
    if not ref:
        return out
    for ch in chains:
        if ch == primary:
            continue
        q = ca(ch)
        shared = sorted(set(ref) & set(q))
        if len(shared) < 20:
            report("  chain %s shares only %d CA with %s; no copy made"
                   % (ch, len(shared), primary))
            continue
        P = np.array([ref[r] for r in shared])
        Q = np.array([q[r] for r in shared])
        pc, qc = P.mean(0), Q.mean(0)
        U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        t = qc - R @ pc
        rms = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
        out.append((ch, R, t))
        report("  protomer %s: %d CA, rms %.2f A" % (ch, len(shared), rms))
    return out


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
    from .frame import Frame, lift_to_front, pymol_view, spread_in_view
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
    symmetry: List[Tuple[str, np.ndarray, np.ndarray]] = []
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
        cmd.load(spec.up_pdb, up_obj + "_all")
        # The PAIR is one protomer of each. Everything -- classes, patches,
        # envelopes -- is built on it, and the assembly is then made by
        # replicating that pair with the upstream oligomer's own symmetry, so
        # the two scenes show the same things and differ only in how many
        # copies. Building the classes across all chains instead would give
        # an assembly with no pair inside it to show on its own.
        up_all_chains = [c for c in cmd.get_chains(up_obj + "_all") if c.strip()]
        primary = (spec.chains("up") or up_all_chains or [""])[0]
        cmd.create(up_obj, "%s_all and chain %s" % (up_obj, primary)
                   if primary else up_obj + "_all")
        cmd.load(spec.dn_pdb, dn_obj)
        cmd.transform_object(dn_obj, _pymol_matrix(placement.rotation,
                                                   placement.translation))
        symmetry[:] = _protomer_transforms(cmd, up_obj + "_all", primary,
                                           up_all_chains, report)
        cmd.delete(up_obj + "_all")

    cmd.hide("everything", "all")

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
    # ONE OBJECT PER CLASS, per enzyme.
    #
    # The p5d original colours the classes in place on one surface and keeps
    # 31 named SELECTIONS beside it (GDH_surface, GDH_avoid, GDH_patch_01 ...).
    # A selection can be picked but not hidden: hiding it does not hide the
    # parent's surface. With seven translucent classes on one object every
    # colour is seen through every colour behind it, which is what made the
    # elements hard to tell apart.
    #
    # Objects fix that, and the classes make it clean: they PARTITION the
    # residues, so the union of the per-class objects is exactly the whole
    # surface. All on = the complete enzyme; switch one off and only that part
    # disappears. No double-drawing, no gaps.
    painted: Dict[str, int] = {}
    objs = MapObjects()
    class_objs = objs.classes
    for which, path, obj in (("up", spec.classes_up, up_obj),
                             ("dn", spec.classes_dn, dn_obj)):
        if not path:
            continue
        table = read_residue_classes(path)
        for cls in CLASS_ORDER:
            resi = sorted(table.get(cls, {}))
            colour = pal.classes.get(cls)
            if not resi or colour is None:
                continue
            name = "%s_%s" % (spec.name(which), cls.lower())
            cmd.create(name, "%s and resi %s and polymer"
                       % (obj, "+".join(str(r) for r in resi)))
            if cmd.count_atoms(name) == 0:
                cmd.delete(name)
                continue
            cmd.hide("everything", name)
            cmd.show("surface", name)
            rgb = parse_color(colour)
            if cls in pal.background_classes:
                rgb = lighten(rgb, pal.class_lighten)
            key = "c_cls_%s" % cls.lower()
            cmd.set_color(key, [float(c) for c in rgb])
            cmd.color(key, name)
            cmd.set("transparency", spec.transparency, name)
            class_objs.append(name)
            painted[cls] = painted.get(cls, 0) + len(resi)
    for cls in CLASS_ORDER:
        if cls in painted:
            info["objects"].append("%-12s %-8s %d residue(s), own object per enzyme"
                                   % (cls, pal.classes.get(cls), painted[cls]))

    # ligands as objects, one per enzyme, so they can be toggled like the rest
    ligand_objs = objs.ligands
    if spec.ligand_resn:
        cmd.set_color("c_lig", [float(c) for c in parse_color(pal.ligand)])
        for which, obj in (("up", up_obj), ("dn", dn_obj)):
            name = "%s_ligands" % spec.name(which)
            cmd.create(name, "%s and resn %s" % (obj, "+".join(spec.ligand_resn)))
            if cmd.count_atoms(name) == 0:
                cmd.delete(name)
                continue
            cmd.hide("everything", name)
            cmd.show("sticks", name)
            cmd.color("c_lig", name)
            cmd.set("stick_radius", 0.45, name)
            cmd.set("stick_ball", 0, name)
            ligand_objs.append(name)
        if ligand_objs:
            info["objects"].append("ligands as objects: " + ", ".join(ligand_objs))

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
    tunnel_objs = objs.tunnels
    mouth_objs = objs.mouths
    exits: Dict[str, np.ndarray] = {}
    arms: Dict[str, np.ndarray] = {}
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
        if which == "up" and len(xyz) >= 2:
            # one sphere back along the chain: the O_IC arm's far point
            arms["p1"] = np.asarray(xyz[-2], float)
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
        cmd.set_view(pymol_view(frame, orthoscopic=True))
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
    base = {}
    for which, obj in (("up", up_obj), ("dn", dn_obj)):
        twin = obj + "_base"
        cmd.create(twin, obj)
        cmd.hide("everything", twin)
        cmd.show("surface", "%s and polymer" % twin)
        cmd.set_color("c_plain", [float(c) for c in parse_color(pal.neutral)])
        cmd.color("c_plain", "%s and polymer" % twin)
        if spec.ligand_resn:
            cmd.color("c_lig", "%s and resn %s" % (twin, "+".join(spec.ligand_resn)))
            cmd.show("sticks", "%s and resn %s" % (twin, "+".join(spec.ligand_resn)))
        # TRANSLUCENT, and the only translucent thing in the session. It is
        # the whole enzyme, drawn so that a solid element switched on over it
        # reads as "here, on this surface" -- which is the one question this
        # session exists to answer. Every element on top of it is solid.
        cmd.set("transparency", spec.base_transparency, "%s and polymer" % twin)
        cmd.set("transparency", spec.base_transparency, twin)
        base[which] = twin
    info["objects"].append("%s + %s (translucent %s base, transparency %.2f)"
                           % (base["up"], base["dn"], pal.neutral,
                              spec.base_transparency))

    # ---- the patches, in place, in one hue, each labelled, each with the
    # envelope a minimum-length binder on it would occupy
    patch_objs = objs.patches
    env_objs = objs.envelopes
    label_objs = objs.labels
    pending_labels: List[tuple] = []
    for which, path, obj in (("up", spec.patches_up, up_obj),
                             ("dn", spec.patches_dn, dn_obj)):
        if not path:
            continue
        patches = read_patch_table(path)
        shades = patch_gradient(len(patches), pal.patch_base, pal.patch_light,
                                pal.patch_dark, pal.best_is_darkest)
        # EVERY protomer, not just one. The patch objects already span all
        # chains (the selection has no chain term), but the envelope was
        # placed at a single centroid taken from a resi->coord map, where
        # residue numbers collide across chains and one protomer silently
        # won. The point of the assembly view is to see whether binders on
        # neighbouring copies collide, which needs all of them.
        chains = [c for c in cmd.get_chains(obj) if c.strip()] or [""]
        per_chain = {}
        for ch in chains:
            where = "%s and chain %s" % (obj, ch) if ch else obj
            if cmd.count_atoms("%s and polymer" % where) == 0:
                continue
            per_chain[ch] = (patch_centroids(_ca_by_resi(cmd, where), patches),
                             cmd.get_coords("%s and polymer" % where))
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
            # An OBJECT, not a named selection. A selection can be picked but
            # not hidden -- hiding it leaves the parent's surface drawn -- so
            # with everything on one object each element is seen through every
            # element behind it. As an object a patch can be switched off.
            cmd.create(tag, sel)
            if cmd.count_atoms(tag) == 0:
                cmd.delete(tag)
                continue
            cmd.hide("everything", tag)
            cmd.show("surface", tag)
            cmd.color(key, tag)
            cmd.set("transparency", 0.0, tag)   # the patches are the subject
            # A patch and the class object under it are built from the SAME
            # atoms, so their surfaces are coincident and which one is drawn
            # is up to depth-buffer luck -- switching a patch off and on in
            # the classes scene could leave it behind the class surface.
            # A slightly larger probe inflates the patch surface just outside
            # the class one, so the patch always wins. It is a rendering
            # radius, not a coordinate: nothing moves.
            cmd.set("solvent_radius", spec.patch_probe, tag)
            # The base underneath is NOT tinted. It was, so that the scene
            # still read with the patch objects off -- which meant switching
            # a patch off changed nothing visible, and toggling, the whole
            # point of the session, silently did not work.
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
            # One envelope object per patch, carrying a blob chain for EVERY
            # protomer that has this patch -- the same object, so it toggles
            # with the patch it belongs to.
            all_pts, all_rad, placed = [], [], 0
            for _ch, (cens, surf) in sorted(per_chain.items()):
                c = cens.get(patch.label)
                if c is None:
                    continue
                n_ch = patch_outward_normal(surf, c)
                pts_ch, rad_ch = binder_blobs(c, n_ch, spec.binder_length,
                                              spec.binder_step, spec.binder_radius)
                all_pts.extend(pts_ch)
                all_rad.extend(rad_ch)
                placed += 1
            if not all_pts:
                fallback = binder_blobs(centre, normal, spec.binder_length,
                                        spec.binder_step, spec.binder_radius)
                all_pts, all_rad, placed = list(fallback[0]), list(fallback[1]), 1
            pts, radii = np.asarray(all_pts, float), list(all_rad)
            env = tag + "_binder"
            _pseudo_spheres(cmd, env, pts, radii,
                            tuple(float(c) for c in parse_color(pal.binder)),
                            pal.binder_transparency)
            env_objs.append(env)
            info["envelopes"][tag] = {"centroid": [round(float(x), 2) for x in centre],
                                      "blobs": len(pts), "protomers": placed}
        info["patches"][spec.name(which)] = drawn
        report("%s: %d patch(es), shades %s"
               % (spec.name(which), len(drawn),
                  " -> ".join("#%02x%02x%02x" % tuple(round(float(c) * 255) for c in s)
                              for s in shades)))
    # ---- the labels, pushed apart where the collision actually happens
    if pending_labels:
        anchors = np.array([a for _, a, _ in pending_labels])
        if frame is not None:
            # The frame is zoomed to the two enzymes, so anything beyond their
            # extent is at or past the edge of the picture. Bounding the
            # relaxation by that box keeps a crowded label beside its own
            # enzyme instead of trading an overlap for a label in the margin.
            both = np.vstack([up_coords, dn_coords])
            u, v = both @ cam_x, both @ cam_y
            m = float(spec.label_margin)
            bounds = (u.min() - m, u.max() + m, v.min() - m, v.max() + m)
            anchors = anchors + spread_in_view(
                anchors, cam_x, cam_y,
                min_sep=(spec.label_sep_x, spec.label_sep_y), bounds=bounds)
            # ...then forward, so a label naming a patch on the far side is
            # not drawn behind the near surface. Separately from the spread:
            # that one works in the camera plane and must not touch depth,
            # this one touches nothing else.
            lift = lift_to_front(anchors, frame.toward_viewer, both,
                                 spec.label_lift)
            moved = int((np.linalg.norm(lift, axis=1) > 1e-6).sum())
            if moved:
                report("%d label(s) lifted in front of the surface "
                       "(they named far-side patches)" % moved)
            anchors = anchors + lift
        # Draw labels IN FRONT of everything. A label sits at its patch's own
        # depth, so one naming a patch on the far side of the enzyme renders
        # behind the near surface and washes out to unreadable -- which is
        # what happened to two of the FNR names. Separating them in the camera
        # plane cannot fix that: the collision is in depth, and the separation
        # deliberately never changes depth. This does, and only for labels.
        cmd.set("float_labels", 1)   # helps in the viewport; a ray trace ignores it
        cmd.set_color("c_lab", [float(c) for c in parse_color(pal.label)])
        for (lab, _, text), pos in zip(pending_labels, anchors):
            cmd.pseudoatom(lab, pos=[float(v) for v in pos], label=text)
            cmd.set("label_size", 18, lab)
            cmd.set("label_color", "c_lab", lab)
            label_objs.append(lab)
        report("%d label(s) placed, separated in the camera plane" % len(label_objs))
    cmd.deselect()

    # The earlier "pull the two enzymes apart" workaround is gone. It existed
    # because both enzymes were painted alike on one surface each and could not
    # be told apart; now every class is its own object, named per enzyme
    # (GDH_channel, FNR_channel, ...), so they are distinguishable and
    # individually switchable without moving anything.

    # ---- scenes
    # Orthographic. Set as a SETTING as well as in the view numbers: zoom()
    # rewrites the view, and a figure comparing patches should not render a
    # near patch larger than an equal far one. It also makes the label lift
    # free -- under this camera, moving a label along the view axis does not
    # move it on screen at all, so it clears the surface without sliding away
    # from the patch it names.
    cmd.set("orthoscopic", 1)
    cmd.zoom("%s or %s" % (up_obj, dn_obj), buffer=8, complete=1)
    cmd.clip("slab", 3000)   # zoom keeps the rotation set above and fits the pair
    # up_obj / dn_obj are no longer drawn themselves: every visible element is
    # its own object (one per class, one per patch, ligands, mouths, tunnels,
    # envelopes), so any of them can be switched off independently.
    objs.base[:] = [base["up"], base["dn"]]
    note = placement.describe() if placement is not None else ""
    ic_objs = objs.annotation
    if placement is not None and "up" in exits and "dn" in exits and "p1" in arms:
        # The SAME annotation the pose sessions carry, from the same
        # function -- not a second copy. The reference placement is a real
        # iMARS geometry (O_IC held at 180, D_IC raised from its minimum
        # until nothing clashes), so the line, the arm, the arc and the
        # labels all mean here exactly what they mean there, and carry the
        # same object names. At O_IC 180 the three points are collinear and
        # build_annotation omits the arc rather than inventing a bend.
        #
        # Built, then left OFF: the map's subject is the epitopes, and these
        # numbers describe the placement rather than any patch.
        from .annotation import build_annotation

        objs.annotation[:] = build_annotation(
            cmd,
            {"p1": arms["p1"], "p2": exits["up"], "p3": exits["dn"],
             "d_ic": placement.d_ic, "o_ic": placement.o_ic, "do": placement.do},
            label="reference placement", frame=frame)
        # Black, and pushed clear of both enzymes along the camera's up axis,
        # exactly as the pose session does it -- measured from the drawn
        # extent rather than a fixed offset, so no label lands on a protein
        # whatever the pose or the zoom.
        if frame is not None:
            body = cmd.get_coords("%s or %s" % (up_obj, dn_obj))
            if body is not None and len(body):
                up_v = np.asarray(frame.up, float)
                mid = (np.asarray(exits["up"], float)
                       + np.asarray(exits["dn"], float)) / 2.0
                reach = (body - mid) @ up_v
                for nm, sign, step in (("IC_label_O", 1.0, 0.0),
                                       ("IC_label_D", 1.0, -5.0),
                                       ("IC_label_DO", -1.0, 0.0),
                                       ("IC_label_tag", -1.0, 5.0)):
                    if nm not in cmd.get_names("objects"):
                        continue
                    far = (float(reach.max()) + 6.0 if sign > 0
                           else float(reach.min()) - 6.0)
                    pos = mid + up_v * (far + sign * step)
                    cmd.alter_state(1, nm, "(x, y, z) = (%f, %f, %f)"
                                    % tuple(float(v) for v in pos))
                cmd.set_color("c_black", [0.0, 0.0, 0.0])
                for nm in ic_objs:
                    if nm.startswith("IC_label"):
                        cmd.set("label_color", "c_black", nm)
                        cmd.set("label_size", 20, nm)
        info["objects"].append(
            "IC annotation (%s), shown in both scenes: D_IC %.1f A at "
            "O_IC %.0f deg, DO %.4f. At 180 the arc's plane is undefined, so "
            "it is drawn as a half turn in the VIEW plane -- a presentational "
            "choice, not a measurement"
            % (", ".join(ic_objs), placement.d_ic, placement.o_ic, placement.do))

    # Labels and the iMARS annotation are BUILT but start off: they are in
    # the session for anyone who wants them, and out of the way by default.
    # ONE scene (user, 2026-10-25). `map`, `classes` and `envelopes` were
    # the same object set with different things switched on, and that is a
    # toggle, not a scene -- the same reason the `patches` scene went. What a
    # scene really stores that an object cannot is the CAMERA, and there is
    # one camera here.
    # ---- the assembly: the SAME pair, copied onto every other protomer.
    # Every object is replicated -- classes, patches, ligands, tunnels,
    # mouths and the binder envelopes -- so the two scenes differ only in
    # how many copies they show, which is what makes them comparable.
    copies: List[str] = []
    for ch, R, t in symmetry:
        mat = _pymol_matrix(R, t)
        for src in (objs.base + objs.classes + objs.patches + objs.tunnels
                    + objs.mouths + objs.ligands + objs.envelopes):
            dup = "%s_%s" % (src, ch)
            cmd.create(dup, src)
            cmd.transform_object(dup, mat)
            copies.append(dup)
    if copies:
        info["objects"].append(
            "assembly: %d object(s) copied onto %d further protomer(s) (%s)"
            % (len(copies), len(symmetry), ", ".join(c for c, _, _ in symmetry)))

    # Two scenes, and the ASSEMBLY one is stored last so the session opens
    # on it: its job is to show where binders would sit on every protomer at
    # once, and whether neighbouring ones collide -- which is a question
    # about the assembly, not about one pair.
    pair = objs.scene(envelopes=True, labels=True) + objs.annotation
    scenes = {
        "classes_pair": pair,
        "classes_assembly": pair + copies,
    }
    if ic_objs:
        info["objects"].append(
            "both scenes carry the annotation and the labels, in black and "
            "clear of the complexes")
    for key, objs in scenes.items():
        cmd.disable("all")
        for o in objs:
            cmd.enable(o)
        cmd.scene(key, "store", message=note)
        info["objects"].append("scene %-10s %d object(s)" % (key, len(objs)))
    opening = next(iter(scenes))
    cmd.disable("all")
    for o in scenes[opening]:
        cmd.enable(o)
    cmd.scene(opening, "recall")
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
        up_name=args.get("up_name", "Eup"),
        dn_name=args.get("dn_name", "Edn"),
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

        # every element is its own OBJECT, so any of them can be switched off
        check("E1_patch_01" in objs, "each patch is its own object")
        check(any(o.endswith("_channel") or o.endswith("_allowed") for o in objs),
              "each residue class is its own object (%s)"
              % sorted(o for o in objs if o.startswith("E1_"))[:4])
        check("E1_patch_01_binder" in objs, "every patch gets an envelope")
        check(len(info["envelopes"]) == 3, "three envelopes (%d)" % len(info["envelopes"]))
        blobs = {k: v["blobs"] for k, v in info["envelopes"].items()}
        check(set(blobs.values()) == {4},
              "a 6 A binder at 2 A steps is 4 blobs everywhere (%s)" % blobs)

        # ONE scene. `map`, `patches` and `envelopes` were this one with
        # different objects switched on, which is a toggle, not a scene.
        # What a scene stores that an object cannot is the camera, and there
        # is one camera here.
        stored = cmd.get_scene_list() or []
        check(stored == ["classes_pair", "classes_assembly"],
              "classes_pair then classes_assembly, as in the pose session (%s)"
              % stored)
        cmd.scene("classes_pair", "recall")
        on = set(cmd.get_names("objects", enabled_only=1))
        check(any(o.endswith("_binder") for o in on),
              "the pair scene carries its envelopes (%s)" % sorted(on))
        # The patches now ride ALONG with the classification (user,
        # 2026-10-25): an epitope only means something read against the
        # classes it sits on, and keeping them in separate scenes forced you
        # to answer that by flicking between the two.
        check(any("patch" in o for o in on),
              "the classes scene carries the patches too, so an epitope can "
              "be read against the classes under it")
        check(any(o.endswith("_base") for o in on),
              "...over the translucent whole-enzyme base")
        check("E1_tunnel" in on and "E1_mouth" in on,
              "...and keeps tunnels and mouths, which are objects too")
        check(sum(1 for o in on if o.startswith("E1_") and "_" in o) >= 3,
              "several class objects are shown (%s)" % sorted(on))
        cmd.scene("classes_pair", "recall")
        pair_on = set(cmd.get_names("objects", enabled_only=1))
        cmd.scene("classes_assembly", "recall")
        asm_on = set(cmd.get_names("objects", enabled_only=1))
        check(any(o.endswith("_binder") for o in pair_on),
              "the PAIR scene already carries the binder envelopes")
        check(pair_on <= asm_on,
              "the assembly shows everything the pair does, and more")
        check(all(o.startswith("IC_") for o in pair_on
                  if o.startswith("IC_")),
              "the annotation is in both")
        every = cmd.get_names("objects")
        check(any(o.startswith("IC_") for o in every),
              "the iMARS annotation exists as objects too (%s)"
              % sorted(o for o in every if o.startswith("IC_")))
        check(any(o.startswith("IC_") for o in asm_on),
              "...and is SHOWN, not hidden: the user asked for the numbers")

        t = float(cmd.get("transparency", "E1"))
        # Solid, on purpose. This session is made to be toggled: you switch a
        # patch off and on to see which surface it covers, and through a
        # translucent surface you cannot tell the patch in front from the
        # class behind it. Seeing inside is what hiding an object is for.
        check(t == 0.0, "the enzymes are solid, so a toggled patch reads (%s)" % t)
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
