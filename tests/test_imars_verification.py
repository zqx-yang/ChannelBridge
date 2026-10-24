"""Verification against the iMARS paper's own published data.

Wang J, Ouyang X, Meng S, ..., Shi T, Zhao Y-L, Ni J (2025). Rational
multienzyme architecture design with iMARS. *Cell* 188:1349-1362.

These tests read the paper's supplementary workbook (Table S2) and check that
our DO implementation reproduces the paper's own numbers. They are the only
end-to-end check of the scoring layer that does not need a structure, because
iMARS publishes D_IC, O_IC, DO and the measured yield for 130 fusion enzymes
but no coordinates.

The workbook is copyrighted and is not redistributed here. Point
``CHANNELBRIDGE_IMARS_ST2`` at ``1-s2.0-S0092867424014715-mmc3.xlsx`` to run
these; they skip otherwise.

What is verified, and what is not:

* **verified** -- the yield regression, the derivation of DO from it, which
  constant set is correct, and that our `do_score` matches theirs.
* **not verified here** -- that our *geometry extraction* produces the same
  D_IC and O_IC from a structure. That needs coordinates, which the paper
  does not provide. See docs/imars-verification.md.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from channelbridge.geometry import (
    DO_EXACT,
    DO_PUBLISHED,
    DO_ST2A,
    PR_4CL_STS,
    do_clamped,
    do_score,
)

ST2 = Path(os.environ.get("CHANNELBRIDGE_IMARS_ST2", "")) if os.environ.get(
    "CHANNELBRIDGE_IMARS_ST2"
) else Path(
    "/home/hyyang/programs/ChannelBridge/claude/iMARS/"
    "ScienceDirect_files_25Sep2026_13-27-44.168/1-s2.0-S0092867424014715-mmc3.xlsx"
)

openpyxl = pytest.importorskip("openpyxl", reason="openpyxl is needed to read Table S2")

needs_st2 = pytest.mark.skipif(
    not ST2.is_file(),
    reason="iMARS Table S2 workbook not available (set CHANNELBRIDGE_IMARS_ST2)",
)


def _sheet(name):
    wb = openpyxl.load_workbook(ST2, read_only=True, data_only=True)
    try:
        rows = list(wb[name].iter_rows(values_only=True))
    finally:
        wb.close()
    out, seen_header = [], False
    for r in rows:
        if r and r[0] == "No.":
            seen_header = True
            continue
        if seen_header and r and isinstance(r[0], (int, float)):
            out.append(r)
    return out


@pytest.fixture(scope="module")
def st2a():
    """100 4CL-STS fusions: D_IC, O_IC, yield, their DO."""
    d = _sheet("ST2A")
    return (
        np.array([float(r[4]) for r in d]),
        np.array([float(r[5]) for r in d]),
        np.array([float(r[6]) for r in d]),
        np.array([float(r[7]) for r in d]),
    )


@pytest.fixture(scope="module")
def st2f():
    """The independent 30-linker validation set."""
    d = _sheet("ST2F")
    return (
        np.array([float(r[4]) for r in d]),
        np.array([float(r[5]) for r in d]),
        np.array([float(r[6]) for r in d]),
        np.array([float(r[7]) for r in d]),
    )


@needs_st2
class TestTheirYieldRegression:
    """The quantity with units, from which DO is derived."""

    def test_table_has_the_expected_shape(self, st2a):
        DIC, OIC, _Y, _DO = st2a
        assert len(DIC) == 100
        assert 45.0 < DIC.min() and DIC.max() < 92.0
        assert 0.0 <= OIC.min() and OIC.max() <= 180.0

    def test_we_reproduce_their_published_regression(self, st2a):
        # Paper, Figure 1E and Quantification:
        #   Pr = -0.71*D_IC + 0.107*O_IC + 55.10, R^2 = 0.73
        DIC, OIC, Y, _ = st2a
        M = np.column_stack([DIC, OIC, np.ones_like(DIC)])
        (a, b, c), *_ = np.linalg.lstsq(M, Y, rcond=None)
        assert a == pytest.approx(PR_4CL_STS.a, abs=5e-3)
        assert b == pytest.approx(PR_4CL_STS.b, abs=5e-4)
        assert c == pytest.approx(PR_4CL_STS.c, abs=5e-2)

    def test_their_quoted_r_squared_holds(self, st2a):
        DIC, OIC, Y, _ = st2a
        M = np.column_stack([DIC, OIC, np.ones_like(DIC)])
        coef, *_ = np.linalg.lstsq(M, Y, rcond=None)
        resid = Y - M @ coef
        r2 = 1 - (resid**2).sum() / ((Y - Y.mean()) ** 2).sum()
        assert r2 == pytest.approx(0.735, abs=0.01)  # paper says 0.73


@needs_st2
class TestTheDODerivation:
    """DO is the yield regression confined to [0, 1]."""

    def test_confining_the_regression_gives_the_published_constants(self):
        # Paper: "the minimum value for the DO score to 0 (D_IC = 100 A,
        # O_IC = 0 deg)"; the maximum is the opposite corner.
        lo = PR_4CL_STS(100.0, 0.0)
        hi = PR_4CL_STS(0.0, 180.0)
        rng = hi - lo
        assert lo == pytest.approx(-15.90, abs=0.01)
        assert hi == pytest.approx(74.36, abs=0.01)
        a, b, c = PR_4CL_STS.a / rng, PR_4CL_STS.b / rng, (PR_4CL_STS.c - lo) / rng
        # ... which is DO_EXACT.
        assert (a, b, c) == pytest.approx((DO_EXACT.a, DO_EXACT.b, DO_EXACT.c))
        # The slope terms round exactly to the paper's stated constants.
        assert round(a, 4) == DO_PUBLISHED.a
        assert round(b, 4) == DO_PUBLISHED.b
        # The intercept does not: the derivation gives 0.78662, which rounds to
        # 0.7866, while the paper states 0.7865. A 1e-4 discrepancy, recorded
        # here so it is not rediscovered as a bug in our arithmetic. It shifts
        # every DO by the same 0.0001 and so changes no ranking.
        assert round(c, 4) == 0.7866
        assert abs(c - DO_PUBLISHED.c) == pytest.approx(1.0e-4, abs=2e-5)

    def test_one_angstrom_is_worth_about_six_point_six_degrees(self):
        assert DO_PUBLISHED.degrees_per_angstrom == pytest.approx(6.58, abs=0.01)
        assert DO_EXACT.degrees_per_angstrom == pytest.approx(6.64, abs=0.01)


@needs_st2
class TestOurScoreAgainstTheirTables:
    def test_we_reproduce_their_validation_set_DO_column(self, st2f):
        # Table S2F is their independent 30-linker set. Its DO column was
        # computed with the unrounded derivation.
        DIC, OIC, _, DO = st2f
        ours = np.array([do_score(d, o, DO_EXACT) for d, o in zip(DIC, OIC)])
        assert np.abs(ours - DO).max() < 1e-4

    def test_the_rounded_constants_cost_about_three_thousandths(self, st2f):
        DIC, OIC, _, DO = st2f
        ours = np.array([do_score(d, o, DO_PUBLISHED) for d, o in zip(DIC, OIC)])
        err = np.abs(ours - DO).max()
        assert 1e-4 < err < 5e-3

    def test_their_fit_table_uses_a_different_constant_set(self, st2a):
        # Table S2A's DO column is reproduced only by DO_ST2A, not by the
        # constants the paper states. Recorded so the discrepancy cannot be
        # mistaken for a bug in our implementation.
        DIC, OIC, _, DO = st2a
        st2a_err = np.abs([do_score(d, o, DO_ST2A) for d, o in zip(DIC, OIC)] - DO).max()
        stated_err = np.abs([do_score(d, o, DO_PUBLISHED) for d, o in zip(DIC, OIC)] - DO).max()
        assert st2a_err < 1e-5
        assert stated_err > 2e-2

    def test_the_anomalous_set_does_not_derive_from_their_own_data(self, st2a):
        # Its implied trade-off is 5.63 deg/A; every regression on this table
        # gives 6.6 or more, so it is not a refit of these 100 points.
        DIC, OIC, Y, _ = st2a
        M = np.column_stack([DIC, OIC, np.ones_like(DIC)])
        coef, *_ = np.linalg.lstsq(M, Y, rcond=None)
        assert abs(coef[0] / coef[1]) == pytest.approx(6.63, abs=0.05)
        assert DO_ST2A.degrees_per_angstrom == pytest.approx(5.63, abs=0.02)

    def test_the_stated_constants_predict_yield_at_least_as_well(self, st2a):
        # The arbiter: DO only means anything as a predictor of yield.
        DIC, OIC, Y, _ = st2a
        def r2(coeffs):
            s = np.array([do_score(d, o, coeffs) for d, o in zip(DIC, OIC)])
            return np.corrcoef(s, Y)[0, 1] ** 2
        assert r2(DO_PUBLISHED) >= r2(DO_ST2A)
        assert r2(DO_PUBLISHED) == pytest.approx(0.735, abs=0.002)


@needs_st2
class TestClamping:
    def test_iMARS_clamps_at_zero(self):
        # Their Table S2G reports exactly 0.000000 for this pose.
        raw = do_score(116.525, 68.515, DO_EXACT)
        assert raw < 0.0
        assert do_clamped(116.525, 68.515, DO_EXACT) == 0.0

    def test_clamping_is_a_no_op_inside_the_range(self):
        assert do_clamped(50.0, 120.0) == pytest.approx(do_score(50.0, 120.0))


@needs_st2
class TestExtrapolationHonesty:
    """Our working range sits at the edge of theirs."""

    def test_their_fit_range_is_far_from_our_pose_grid(self, st2a):
        DIC, _, _, _ = st2a
        assert DIC.min() > 45.0  # nothing in the fit below 45 A

    def test_their_validation_set_reaches_down_to_17_angstroms(self, st2f):
        DIC, _OIC, Y, _ = st2f
        assert DIC.min() == pytest.approx(17.386, abs=0.01)
        # and that point is their best fusion, 61.3 mg/L -- so the model is
        # anchored by experiment at the distance our best pose sits at.
        i = int(np.argmin(DIC))
        assert Y[i] == pytest.approx(61.3, abs=0.1)

    def test_DO_still_tracks_yield_on_the_validation_set(self, st2f):
        DIC, OIC, Y, _ = st2f
        s = np.array([do_score(d, o, DO_PUBLISHED) for d, o in zip(DIC, OIC)])
        assert np.corrcoef(s, Y)[0, 1] ** 2 == pytest.approx(0.738, abs=0.01)

    def test_extrapolation_note_fires_below_the_validated_floor(self):
        from channelbridge.geometry import extrapolation_note

        assert extrapolation_note(60.0) is None
        assert "extrapolation" in (extrapolation_note(10.0) or "")
        assert extrapolation_note(17.0) is not None
