"""The iMARS channelling geometry: D_IC, O_IC and the DO score.

Definitions, for an upstream enzyme A whose product leaves through tunnel A
and a downstream enzyme B whose substrate enters through tunnel B::

    P1 = a point on tunnel A, one O_IC arm inward from its exit
    P2 = the exit of tunnel A            (the upstream exit)
    P3 = the exit of tunnel B            (the downstream exit)

    D_IC = |P2 - P3|                     inter-channel distance, A
    O_IC = angle P1-P2-P3                inter-channel orientation, degrees
                                         (180 = tunnel A aimed straight at P3)
    DO   = a * D_IC + b * O_IC + c

O_IC is **one-sided**: it is measured at the upstream exit only, so which
enzyme is upstream changes the score. A short downstream tunnel is therefore
not a problem -- only tunnel A's outer direction orients anything.

The trade-off the coefficients encode: 1 A of D_IC is worth
``abs(a / b)`` degrees of O_IC, which is 6.6 deg for the published values.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

__all__ = [
    "DO_PUBLISHED",
    "DO_ST2A",
    "DOCoefficients",
    "angle_deg",
    "d_ic",
    "dic_equivalent_of_oic",
    "do_score",
    "o_ic",
]


@dataclass(frozen=True)
class DOCoefficients:
    """``DO = a * D_IC + b * O_IC + c``, with provenance attached.

    Two published sets of constants disagree slightly; carrying the name
    along means a reported DO can always be traced to the one used.
    """

    a: float
    b: float
    c: float
    name: str = "unnamed"

    def __call__(self, d: float, o: float) -> float:
        return self.a * d + self.b * o + self.c

    @property
    def degrees_per_angstrom(self) -> float:
        """How many degrees of O_IC buy one angstrom of D_IC."""
        return abs(self.a / self.b)


# ---------------------------------------------------------------------------
# Coefficient sets.
#
# iMARS (Wang et al., Cell 188:1349-1362, 2025) states the DO equation three
# times -- Figure 1H, the Methods, and the Quantification section -- as
#
#     DO = -0.0079 * D_IC + 0.0012 * O_IC + 0.7865
#
# and that is the default here. It is also exactly what the paper's own
# derivation gives: their yield regression over 100 4CL-STS fusions is
# Pr = -0.71*D_IC + 0.107*O_IC + 55.10 (R^2 0.735, reproduced from their
# Table S2A to four decimals), confined to [0, 1] by putting 0 at the worst
# corner (D_IC = 100 A, O_IC = 0 deg) and 1 at the best (D_IC = 0,
# O_IC = 180 deg). That division by 90.26 yields DO_EXACT below, which rounds
# to the stated constants.
#
# Their supplementary tables are not self-consistent; see
# docs/imars-verification.md for the full audit. In short: Table S2F was
# computed with DO_EXACT, Table S2A with DO_ST2A, and Table S2G with a
# mixture plus a clamp at 0.
# ---------------------------------------------------------------------------

#: The published equation, as stated in iMARS Figure 1H and Methods. Default.
DO_PUBLISHED = DOCoefficients(-0.0079, 0.0012, 0.7865, name="imars-fig2c")

#: The same derivation carried out without rounding: -0.71/90.26, 0.107/90.26
#: and 71.0/90.26. Reproduces the DO column of iMARS Table S2F -- their
#: independent 30-linker validation set -- to 5e-05, where the rounded
#: constants are off by up to 2.9e-03. Use this to reproduce their reported
#: numbers; use :data:`DO_PUBLISHED` to match what the paper instructs.
DO_EXACT = DOCoefficients(
    -0.71 / 90.26, 0.107 / 90.26, 71.0 / 90.26, name="imars-exact"
)

#: The constants that reproduce the DO column of iMARS Table S2A (to 4e-07),
#: the 100-fusion table the model was fit on. They do **not** derive from any
#: regression on that table: they imply 5.63 degrees of O_IC per angstrom of
#: D_IC where the table's own yield regression gives 6.63, and they predict
#: the measured yields very slightly worse (r^2 0.7323 vs 0.7350). Provided
#: for reproducing that table, not recommended for scoring.
DO_ST2A = DOCoefficients(-0.00801133, 0.00142190, 0.77969380, name="imars-st2a")

#: The yield regression the DO score is derived from, in mg/L of resveratrol
#: for the 4CL-STS system: ``Pr = -0.71*D_IC + 0.107*O_IC + 55.10``. Kept
#: because it, not DO, is the quantity with units and an R^2.
PR_4CL_STS = DOCoefficients(-0.71, 0.107, 55.10, name="imars-pr-4cl-sts")

#: iMARS confines DO to [0, 1]; their Table S2G reports 0.000000 for a pose
#: at D_IC 116.5 A / O_IC 68.5 deg, where the unclamped value is -0.0488.
DO_RANGE = (0.0, 1.0)

#: The D_IC range the DO model was **fit** on (iMARS Table S2A), and the range
#: over which it was validated against measured yields (Table S2F). Scoring
#: far outside these is extrapolation and should be said so.
DO_FIT_DIC_RANGE = (45.1, 91.5)
DO_VALIDATED_DIC_RANGE = (17.4, 90.3)


def _as_point(p, what: str) -> np.ndarray:
    arr = np.asarray(p, dtype=float).reshape(-1)
    if arr.shape != (3,):
        raise ValueError("%s must be a 3-vector, got shape %r" % (what, np.shape(p)))
    if not np.all(np.isfinite(arr)):
        raise ValueError("%s contains a non-finite coordinate: %r" % (what, p))
    return arr


def d_ic(p2, p3) -> float:
    """Inter-channel distance: the straight-line gap between the two exits."""
    return float(np.linalg.norm(_as_point(p3, "p3") - _as_point(p2, "p2")))


def angle_deg(a, b, c, *, degenerate="raise") -> float:
    """The angle at *b* in the path a-b-c, in degrees.

    Parameters
    ----------
    degenerate:
        What to do when an arm has zero length, so no angle is defined:
        ``"raise"`` (default) or ``"nan"``.

    A zero-length arm is the failure mode that motivated
    :meth:`channelbridge.tunnel.Tunnel.axis_point`: measuring O_IC between
    two adjacent spheres of a densely sampled CAVER tunnel gives an arm of a
    fraction of an angstrom, where sphere-placement jitter dominates the
    direction. Raising here keeps that from passing silently as a number.
    """
    a = _as_point(a, "a")
    b = _as_point(b, "b")
    c = _as_point(c, "c")
    u, v = a - b, c - b
    nu, nv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    if nu == 0.0 or nv == 0.0:
        if degenerate == "nan":
            return float("nan")
        raise ValueError(
            "angle undefined: arm length |a-b|=%.4g, |c-b|=%.4g (a zero-length "
            "arm usually means two coincident tunnel spheres)" % (nu, nv)
        )
    cos = float(np.dot(u, v) / (nu * nv))
    return float(np.degrees(np.arccos(min(1.0, max(-1.0, cos)))))


def o_ic(p1, p2, p3, *, degenerate="raise") -> float:
    """Inter-channel orientation: the angle P1-P2-P3 at the upstream exit."""
    return angle_deg(p1, p2, p3, degenerate=degenerate)


def do_score(d: float, o: float, coeffs: DOCoefficients = DO_PUBLISHED) -> float:
    """The DO score for a given D_IC (A) and O_IC (degrees).

    No clamping: DO is a linear fit, and a value outside ``[0, 1]`` means the
    pose is outside the range the fit was made on. That is information, not
    something to hide.
    """
    d = float(d)
    o = float(o)
    if not np.isfinite(d) or not np.isfinite(o):
        raise ValueError("DO needs finite D_IC and O_IC, got %r and %r" % (d, o))
    if d < 0:
        raise ValueError("D_IC cannot be negative, got %r" % (d,))
    if not (0.0 <= o <= 180.0):
        raise ValueError("O_IC must be in [0, 180] degrees, got %r" % (o,))
    return float(coeffs(d, o))


def dic_equivalent_of_oic(delta_oic: float, coeffs: DOCoefficients = DO_PUBLISHED) -> float:
    """Angstroms of D_IC worth *delta_oic* degrees of O_IC, at equal DO."""
    return float(abs(coeffs.b / coeffs.a) * float(delta_oic))


def do_clamped(d: float, o: float, coeffs: DOCoefficients = DO_PUBLISHED) -> float:
    """:func:`do_score` confined to ``[0, 1]``, as iMARS reports it.

    Report both: the raw value says how far outside the fitted range a pose
    sits, and the clamped one is what is comparable with published DO values.
    """
    return float(min(DO_RANGE[1], max(DO_RANGE[0], do_score(d, o, coeffs))))


def extrapolation_note(d: float) -> Optional[str]:
    """Warn when a D_IC sits outside the range DO was fitted or validated on.

    The DO regression was fitted on fusion proteins whose channel exits were
    45-91 A apart, and validated against measured yields down to 17.4 A. A
    binder-held enzyme pair can be far closer than anything in that data, and
    DO is a linear fit with no support there: it keeps rising as D_IC falls,
    which is an assumption, not a measurement.
    """
    lo_fit, hi_fit = DO_FIT_DIC_RANGE
    lo_val, hi_val = DO_VALIDATED_DIC_RANGE
    if lo_val <= d <= hi_val:
        if d < lo_fit:
            return (
                "D_IC %.1f A is below the fitted range (%.1f-%.1f A) but within "
                "the range validated against measured yields (%.1f-%.1f A)"
                % (d, lo_fit, hi_fit, lo_val, hi_val)
            )
        return None
    return (
        "D_IC %.1f A is outside the range DO was validated on (%.1f-%.1f A); "
        "DO here is an extrapolation of a linear fit" % (d, lo_val, hi_val)
    )


def geometry_from_exits(
    upstream_arm_point, upstream_exit, downstream_exit, coeffs: DOCoefficients = DO_PUBLISHED
) -> Tuple[float, float, float]:
    """Return ``(D_IC, O_IC, DO)`` for one pose.

    *upstream_arm_point* is P1 and should come from
    :meth:`channelbridge.tunnel.Tunnel.axis_point`, which measures it over a
    fixed physical arc length rather than one sphere step.
    """
    d = d_ic(upstream_exit, downstream_exit)
    o = o_ic(upstream_arm_point, upstream_exit, downstream_exit)
    return d, o, do_score(d, o, coeffs)
