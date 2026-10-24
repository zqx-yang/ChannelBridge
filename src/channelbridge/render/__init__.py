"""Figures: channel-aware .pse sessions and editable SVGs.

Follows the conventions the p5a/p5b/p5c figure sessions settled on: one .pse
per camera frame with every variant as an object plus a PyMOL scene, and ONE
editable SVG per view rather than a set of PNGs. The camera moves, never the
coordinates of a posed complex.

Appearance is data, not policy: :data:`~channelbridge.render.style.DEFAULT_STYLE`
holds the project's defaults and every colour and transparency in it can be
replaced from a pair file's ``render:`` section or from the command line.
"""

from __future__ import annotations

from .frame import Frame, angle_frame, parse_selection, pymol_view, rotated, select_poses
from .style import (
    DEFAULT_STYLE,
    ColorSpec,
    Palette,
    RenderStyle,
    Transparency,
    parse_color,
)

__all__ = [
    "DEFAULT_STYLE",
    "ColorSpec",
    "Frame",
    "Palette",
    "RenderStyle",
    "Transparency",
    "angle_frame",
    "parse_color",
    "parse_selection",
    "pymol_view",
    "rotated",
    "select_poses",
]
