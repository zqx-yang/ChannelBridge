"""Colours and render settings for channel figures, all overridable.

The defaults are the look the p5a/p5b/p5c figure sessions settled on:
ray-traced with a black outline on a flat background, surfaces lightened so
the outline reads against them, cofactors as unlightened magenta sticks.
Every one of them can be replaced from a pair file or the command line --
:data:`DEFAULT_STYLE` is a starting point, not a house rule.

Two defaults are decisions rather than taste, and are worth knowing before
overriding them:

* **Lighter hues.** Each surface colour is moved toward white by
  ``c + (1 - c) * L``. With a saturated surface the black outline of
  ``ray_trace_mode 1`` disappears into it. If you set a dark colour, raise
  its ``lighten`` too.
* **Both enzymes at 65 % transparency, nothing opaque.** The p5b/p5c rule was
  "designed chain opaque, everything else 65 %", but a posed enzyme pair has
  no designed binder yet. Leaving both transparent makes the channel geometry
  the subject and removes the occlusion-clipping step, which exists only
  because an opaque chain hides the 2D overlay. Making a chain opaque again
  means reinstating that step.

The flat settings exist because ``ray_trace_mode 1`` outlines a bumpy tunnel
surface as dots and creases; the silhouette is traced in 2D from a flat-lit
render instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

__all__ = [
    "BACKGROUNDS",
    "DEFAULT_STYLE",
    "DESIGNED_BINDER_HEX",
    "FLAT_SETTINGS",
    "NAMED_COLORS",
    "RAY_SETTINGS",
    "ColorSpec",
    "Palette",
    "RenderStyle",
    "Transparency",
    "hex_to_rgb",
    "lighten",
    "parse_color",
    "rgb_to_hex",
]

#: The designed binder's blue, used throughout this project's figures. Still
#: only a default: pass ``binder`` in a colour override to change it.
DESIGNED_BINDER_HEX = "#0065bd"

#: Names accepted wherever a colour is expected. Deliberately short: anything
#: else must be given as hex or as three numbers, so a typo is an error rather
#: than a silently wrong colour.
NAMED_COLORS: Dict[str, Tuple[float, float, float]] = {
    "white": (1.0, 1.0, 1.0),
    "black": (0.0, 0.0, 0.0),
    "red": (1.0, 0.0, 0.0),
    "crimson": (0.9, 0.1, 0.2),
    "orange": (1.0, 0.5, 0.0),
    "yellow": (1.0, 0.9, 0.0),
    "green": (0.0, 1.0, 0.0),
    "teal": (0.0, 0.6, 0.55),
    "cyan": (0.0, 1.0, 1.0),
    "blue": (0.0, 0.0, 1.0),
    "purple": (0.55, 0.1, 0.9),
    "magenta": (1.0, 0.0, 1.0),
    "grey": (0.7, 0.7, 0.7),
    "gray": (0.7, 0.7, 0.7),
    # PyMOL's own names for the residue-class palette, with PyMOL's values, so
    # that a figure built here and one built in the viewer agree. The numbered
    # greys follow PyMOL's rule: greyNN is NN per cent white.
    "wheat": (0.99, 0.82, 0.65),
    "pink": (1.0, 0.65, 0.85),
    "palegreen": (0.65, 0.9, 0.65),
    "deepteal": (0.1, 0.6, 0.6),
    "grey20": (0.2, 0.2, 0.2),
    "grey30": (0.3, 0.3, 0.3),
    "grey40": (0.4, 0.4, 0.4),
    "grey50": (0.5, 0.5, 0.5),
    "grey60": (0.6, 0.6, 0.6),
    "grey70": (0.7, 0.7, 0.7),
    "grey80": (0.8, 0.8, 0.8),
}
# the American spellings of the numbered greys, so neither has to be remembered
NAMED_COLORS.update({
    k.replace("grey", "gray"): v for k, v in NAMED_COLORS.items() if k.startswith("grey")
})

BACKGROUNDS = ("white", "E6F2F9", "CFE6F3", "BADBEE")


def hex_to_rgb(value: str) -> Tuple[float, float, float]:
    """``"#0065bd"`` to ``(0.0, 0.396..., 0.741...)``."""
    s = str(value).strip().lstrip("#")
    if len(s) != 6:
        raise ValueError("expected a 6-digit hex colour, got %r" % (value,))
    try:
        return tuple(int(s[i : i + 2], 16) / 255.0 for i in (0, 2, 4))
    except ValueError:
        raise ValueError("%r is not a hex colour" % (value,)) from None


def rgb_to_hex(rgb) -> str:
    """``(0.9, 0.05, 0.05)`` to ``"#e60d0d"``.

    The inverse of :func:`hex_to_rgb`, needed because PyMOL wants floats and
    SVG wants hex. Values are clamped rather than rejected: a colour that
    arrives marginally out of range from arithmetic (a lighten of 1.0 on a
    channel already at 1.0) should still produce a figure.
    """
    vals = [float(c) for c in rgb]
    if len(vals) != 3:
        raise ValueError("expected three channels, got %r" % (rgb,))
    return "#%02x%02x%02x" % tuple(round(min(1.0, max(0.0, float(c))) * 255) for c in vals)


def parse_color(spec) -> Tuple[float, float, float]:
    """A colour from a name, a hex string, or three numbers.

    Three numbers are read as 0-1 floats unless any of them exceeds 1, in
    which case they are read as 0-255. That ambiguity is resolved this way
    rather than guessed per channel, so ``[1, 0, 0]`` is red and
    ``[255, 0, 0]`` is also red.
    """
    if isinstance(spec, str):
        key = spec.strip().lower()
        if key in NAMED_COLORS:
            return NAMED_COLORS[key]
        return hex_to_rgb(spec)
    if isinstance(spec, (list, tuple)):
        if len(spec) != 3:
            raise ValueError("a colour needs 3 numbers, got %r" % (spec,))
        vals = [float(v) for v in spec]
        if any(v > 1.0 for v in vals):
            if any(v > 255.0 or v < 0.0 for v in vals):
                raise ValueError("colour channels out of range in %r" % (spec,))
            return tuple(v / 255.0 for v in vals)
        if any(v < 0.0 for v in vals):
            raise ValueError("colour channels out of range in %r" % (spec,))
        return tuple(vals)
    raise ValueError(
        "cannot read a colour from %r; use a name (%s), hex, or three numbers"
        % (spec, ", ".join(sorted(NAMED_COLORS)[:4]) + ", ...")
    )


def lighten(rgb, amount: float) -> Tuple[float, float, float]:
    """Move *rgb* toward white by *amount* in ``[0, 1]``: ``c + (1-c)*L``."""
    if not 0.0 <= amount <= 1.0:
        raise ValueError("lighten amount must be in [0, 1], got %r" % (amount,))
    out = []
    for c in rgb:
        if not 0.0 <= c <= 1.0:
            raise ValueError("colour channels must be in [0, 1], got %r" % (rgb,))
        out.append(c + (1.0 - c) * amount)
    return tuple(out)


@dataclass(frozen=True)
class ColorSpec:
    """One role's colour: a base hue plus how far it is moved toward white."""

    rgb: Tuple[float, float, float]
    lighten: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "rgb", parse_color(self.rgb))
        lighten(self.rgb, self.lighten)  # validates the amount

    @property
    def final(self) -> Tuple[float, float, float]:
        """The colour as it is handed to PyMOL."""
        return lighten(self.rgb, self.lighten)

    @classmethod
    def from_spec(cls, spec, default: Optional["ColorSpec"] = None) -> "ColorSpec":
        """Read an override.

        A bare colour (``"#4caf50"``, ``"green"``, ``[0, 1, 0]``) keeps the
        role's existing ``lighten``, so changing a hue does not silently
        change how washed out it is. A mapping may set either or both::

            upstream: "#4caf50"
            downstream: {color: orange, lighten: 0.4}
            ligand: {lighten: 0.2}
        """
        keep = default.lighten if default is not None else 0.0
        if isinstance(spec, Mapping):
            unknown = sorted(set(spec) - {"color", "colour", "lighten"})
            if unknown:
                raise ValueError(
                    "unknown colour key(s) %s; use 'color' and/or 'lighten'"
                    % ", ".join(repr(u) for u in unknown)
                )
            if "color" in spec and "colour" in spec:
                raise ValueError("give either 'color' or 'colour', not both")
            raw = spec.get("color", spec.get("colour"))
            amount = float(spec.get("lighten", keep))
            if raw is None:
                if default is None:
                    raise ValueError("a colour override needs a 'color' here")
                return replace(default, lighten=amount)
            return cls(parse_color(raw), amount)
        return cls(parse_color(spec), keep)


def _default_entries() -> Dict[str, ColorSpec]:
    return {
        # upstream is drawn on the LEFT, downstream on the right
        "upstream": ColorSpec((0.0, 1.0, 0.0), 0.55),
        "downstream": ColorSpec((1.0, 0.5, 0.0), 0.50),
        "binder": ColorSpec(hex_to_rgb(DESIGNED_BINDER_HEX), 0.50),
        # cofactors are not lightened: they must read through two surfaces
        "ligand": ColorSpec((1.0, 0.0, 1.0), 0.0),
        "copy": ColorSpec((0.7, 0.7, 0.7), 0.45),
        # the channel volume, never lightened: the 2D fill is drawn at 30 %
        # opacity and has to stay recognisably red
        "tunnel": ColorSpec((0.9, 0.05, 0.05), 0.0),
        "annotation": ColorSpec((0.9, 0.05, 0.05), 0.0),
    }


@dataclass(frozen=True)
class Palette:
    """Role to colour. Immutable; :meth:`with_overrides` returns a new one."""

    entries: Dict[str, ColorSpec] = field(default_factory=_default_entries)

    def __post_init__(self):
        for role, spec in self.entries.items():
            if not isinstance(spec, ColorSpec):
                raise TypeError("entry %r is not a ColorSpec" % (role,))

    @property
    def roles(self) -> List[str]:
        return sorted(self.entries)

    def spec(self, role: str) -> ColorSpec:
        if role not in self.entries:
            raise KeyError(
                "no colour for role %r; known: %s" % (role, ", ".join(self.roles))
            )
        return self.entries[role]

    def rgb(self, role: str) -> Tuple[float, float, float]:
        """The final, lightened colour for *role*."""
        return self.spec(role).final

    def base(self, role: str) -> Tuple[float, float, float]:
        """The unlightened hue for *role*."""
        return self.spec(role).rgb

    def with_overrides(self, overrides: Optional[Mapping]) -> "Palette":
        """A copy with some roles replaced.

        Unknown roles are an error: a figure silently keeping its default
        because ``upsteam`` was misspelt is worse than a failed run.
        """
        if not overrides:
            return self
        unknown = sorted(set(overrides) - set(self.entries))
        if unknown:
            raise ValueError(
                "unknown colour role(s) %s; known: %s"
                % (", ".join(repr(u) for u in unknown), ", ".join(self.roles))
            )
        new = dict(self.entries)
        for role, spec in overrides.items():
            new[role] = ColorSpec.from_spec(spec, default=self.entries[role])
        return Palette(new)

    def indistinguishable(self, roles: Sequence[str] = (), min_sep: float = 0.30) -> List[str]:
        """Pairs of roles whose final colours are too close to tell apart.

        A warning, not an error -- someone may want two greys on purpose.
        Compared on the sum of absolute channel differences of the *final*
        colours, because that is what ends up in the figure.
        """
        check = list(roles) or [r for r in self.roles if r != "annotation"]
        out = []
        for i, a in enumerate(check):
            for b in check[i + 1 :]:
                ca, cb = self.rgb(a), self.rgb(b)
                if sum(abs(x - y) for x, y in zip(ca, cb)) < min_sep:
                    out.append("%s/%s" % (a, b))
        return out


def _default_transparency() -> Dict[str, float]:
    # Only the tunnel and the inspection copy differ from the default; every
    # protein role inherits it.
    #
    # "inspect" is not a palette role: it is how transparent the surfaces are
    # in the scenes meant for LOOKING AT a structure (classes, occupancy,
    # the hotspot map) rather than for illustrating it. They were solid at
    # first, which hid the cofactors the patches are positioned against; 0.55
    # is the value the p5d map script uses for the same job.
    return {"tunnel": 0.75, "inspect": 0.55}


@dataclass(frozen=True)
class Transparency:
    """How see-through each role is drawn. ``0`` opaque, ``1`` invisible.

    A single ``default`` covers every role, with per-role exceptions, so a
    palette can gain a role without anyone having to list it here.
    """

    default: float = 0.65
    per_role: Dict[str, float] = field(default_factory=_default_transparency)

    def __post_init__(self):
        self._check("default", self.default)
        for role, v in self.per_role.items():
            self._check(role, v)

    @staticmethod
    def _check(name, v):
        if not 0.0 <= float(v) <= 1.0:
            raise ValueError(
                "transparency %s must be in [0, 1] (0 opaque, 1 invisible), got %r"
                % (name, v)
            )

    def for_role(self, role: str) -> float:
        return float(self.per_role.get(role, self.default))

    def is_opaque(self, role: str) -> bool:
        return self.for_role(role) == 0.0

    @property
    def opaque_roles(self) -> List[str]:
        return sorted(r for r, v in self.per_role.items() if float(v) == 0.0)

    def with_overrides(self, spec, roles: Sequence[str] = ()) -> "Transparency":
        """Apply an override.

        A bare number sets the default for every role::

            transparency: 0.4

        A mapping sets the default and/or individual roles::

            transparency:
              default: 0.65
              tunnel: 0.75
              binder: 0.0        # opaque
        """
        if spec is None:
            return self
        if isinstance(spec, (int, float)) and not isinstance(spec, bool):
            return replace(self, default=float(spec))
        if not isinstance(spec, Mapping):
            raise ValueError(
                "transparency takes a number or a mapping of roles, got %r" % (spec,)
            )
        known = set(roles) | set(self.per_role) | {"default"}
        unknown = sorted(set(spec) - known)
        if unknown:
            raise ValueError(
                "transparency: unknown role(s) %s; known: %s"
                % (", ".join(repr(u) for u in unknown), ", ".join(sorted(known)))
            )
        new = dict(self.per_role)
        default = float(spec.get("default", self.default))
        for role, v in spec.items():
            if role != "default":
                new[role] = float(v)
        return Transparency(default, new)


@dataclass(frozen=True)
class RenderStyle:
    """Everything a figure's appearance depends on. All of it overridable."""

    palette: Palette = field(default_factory=Palette)
    transparency: Transparency = field(default_factory=Transparency)
    background: str = "white"
    stick_radius: float = 0.45
    #: The pair file's ``render.svg:`` block, kept raw. The SVG assembler
    #: (:mod:`channelbridge.render.svg`) validates and consumes it; holding it
    #: here means a pair file carries one ``render:`` section rather than two,
    #: and an unknown key inside it is still caught -- by the assembler, which
    #: knows what the keys mean.
    svg: Mapping = field(default_factory=dict)
    #: The pair file's ``render.hotspots:`` block, kept raw for the same
    #: reason, and consumed by :class:`channelbridge.render.hotspots.ClassPalette`.
    hotspots: Mapping = field(default_factory=dict)

    _ALLOWED = (
        "colors", "colours", "palette",
        "transparency", "background", "stick_radius", "svg", "hotspots",
    )

    def __post_init__(self):
        if float(self.stick_radius) <= 0:
            raise ValueError("stick_radius must be positive, got %r" % (self.stick_radius,))
        if self.background != "white":
            parse_color(self.background)

    @classmethod
    def from_dict(cls, raw: Optional[Mapping], where: str = "render") -> "RenderStyle":
        """Build from a pair file's ``render:`` section."""
        raw = dict(raw or {})
        unknown = sorted(set(raw) - set(cls._ALLOWED))
        if unknown:
            raise ValueError(
                "%s: unknown key(s) %s; allowed: %s"
                % (where, ", ".join(repr(u) for u in unknown), ", ".join(sorted(cls._ALLOWED)))
            )
        given = [k for k in ("colors", "colours", "palette") if k in raw]
        if len(given) > 1:
            raise ValueError("%s: give only one of %s" % (where, ", ".join(given)))
        palette = Palette().with_overrides(raw.get(given[0]) if given else None)
        return cls(
            palette=palette,
            transparency=Transparency().with_overrides(
                raw.get("transparency"), roles=palette.roles
            ),
            background=str(raw.get("background", "white")),
            stick_radius=float(raw.get("stick_radius", 0.45)),
            svg=dict(raw.get("svg") or {}),
            hotspots=dict(raw.get("hotspots") or {}),
        )

    def transparency_for(self, role: str) -> float:
        """How transparent *role* is drawn."""
        return self.transparency.for_role(role)

    @property
    def needs_occlusion_clipping(self) -> bool:
        """True when an opaque role would hide the 2D overlay.

        With nothing opaque the tunnel outline can be drawn straight over the
        protein raster. An opaque chain in front of it has to be subtracted
        from the overlay (the ``occ=`` step of the p5c SVG script), so the
        caller needs to know.
        """
        return any(self.transparency.is_opaque(r) for r in self.palette.roles)

    def summary(self) -> str:
        t = self.transparency
        lines = [
            "background      %s" % self.background,
            "transparency    default %.2f" % t.default,
        ]
        for role in sorted(t.per_role):
            lines.append("  %-11s %.2f%s" % (role, t.for_role(role),
                                             "  (opaque)" if t.is_opaque(role) else ""))
        if self.needs_occlusion_clipping:
            lines.append("NOTE  an opaque role is set, so the 2D overlay must be "
                         "occlusion-clipped against it")
        else:
            lines.append("opaque          none -- the tunnels are the subject")
        lines.append("colours")
        for role in self.palette.roles:
            s = self.palette.spec(role)
            r, g, b = s.final
            lines.append("  %-11s #%02x%02x%02x  (base #%02x%02x%02x, lighten %.2f, "
                         "transparency %.2f)"
                         % (role, int(r * 255 + .5), int(g * 255 + .5), int(b * 255 + .5),
                            int(s.rgb[0] * 255 + .5), int(s.rgb[1] * 255 + .5),
                            int(s.rgb[2] * 255 + .5), s.lighten, self.transparency_for(role)))
        close = self.palette.indistinguishable()
        if close:
            lines.append("  WARNING  hard to tell apart: %s" % ", ".join(close))
        return "\n".join(lines)


#: The project default. Not a constraint: override any of it.
DEFAULT_STYLE = RenderStyle()

#: Ray-traced look: black outline, flat background, no shadows or fog.
RAY_SETTINGS: Dict[str, object] = {
    "ray_opaque_background": 1,
    "surface_quality": 1,
    "ray_trace_mode": 1,
    "ray_trace_color": "black",
    "ray_trace_gain": 0.8,
    "ray_trace_disco_factor": 0.5,
    "ray_shadows": 0,
    "ray_trace_fog": 0,
    "depth_cue": 0,
    "ambient": 0.55,
    "direct": 0.35,
    "reflect": 0.3,
    "specular": 0.0,
    "spec_reflect": 0.0,
    "antialias": 2,
    "transparency_mode": 3,
}

#: Flat-lit look for the mask / mouth / footprint helper renders: pure black
#: on white, no outline, so a silhouette can be thresholded from the image.
FLAT_SETTINGS: Dict[str, object] = {
    "ray_trace_mode": 0,
    "ambient": 1.0,
    "direct": 0.0,
    "reflect": 0.0,
    "specular": 0.0,
    "spec_reflect": 0.0,
    "ray_shadows": 0,
    "ray_trace_fog": 0,
    "depth_cue": 0,
    "antialias": 2,
}
