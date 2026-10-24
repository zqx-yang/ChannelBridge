"""ChannelBridge -- channel-aware pose search for bispecific binder design.

Given two cascade enzymes and the tunnel each uses for the channelled
species, ChannelBridge finds rigid-body poses in which the upstream product
channel opens toward the downstream substrate channel, scores them with the
iMARS DO metric, and writes the inputs a binder-design run needs.

It wraps CAVER (Chovancova et al. 2012, caver.cz) for tunnel detection and
uses the DO score of the iMARS work. Neither is distributed with it and
neither endorses it; both must be cited by anyone publishing results.
"""

from __future__ import annotations

__version__ = "0.1.0.dev0"

from .geometry import (
    DO_PUBLISHED,
    DO_ST2A,
    DOCoefficients,
    d_ic,
    do_score,
    o_ic,
)
from .tunnel import (
    BURIAL_RADIUS,
    BURIAL_SURFACE_MIN,
    Tunnel,
    burial_verdict,
)

__all__ = [
    "BURIAL_RADIUS",
    "BURIAL_SURFACE_MIN",
    "DO_PUBLISHED",
    "DO_ST2A",
    "DOCoefficients",
    "Tunnel",
    "__version__",
    "burial_verdict",
    "d_ic",
    "do_score",
    "o_ic",
]
