"""Cross-build regression: the normalisation must reconcile two CAVER builds.

The claim this package rests on is that the exit, span and O_IC arm of a
CAVER tunnel are artefacts of the build's sphere sampling, and that
normalising the chain removes the difference. These tests check that against
the two real builds available here, so the claim cannot rot silently.

They are skipped when the reference data is absent, which is the normal case
for anyone who is not on this cluster. Point ``CHANNELBRIDGE_REF_P5`` and
``CHANNELBRIDGE_REF_P6`` at the two ``out_<tag>`` folders to run them.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from channelbridge.io.pdb import read_pdb
from channelbridge.tunnel import BURIAL_SURFACE_MIN, Tunnel

P5_ROOT = Path(
    os.environ.get(
        "CHANNELBRIDGE_REF_P5",
        "/home/hyyang/projects/P5_BBD/p5d/input-prep/caver/output",
    )
)
P6_ROOT = Path(
    os.environ.get(
        "CHANNELBRIDGE_REF_P6",
        "/home/hyyang/projects/P6_CBD/p6a/input-prep/caver/output",
    )
)
STRUCTURE = Path(
    os.environ.get(
        "CHANNELBRIDGE_REF_STRUCTURE",
        "/home/hyyang/projects/P6_CBD/p6a/input-prep/caver/in/GDH/p5d_t1_GDH_BGC-NAP.pdb",
    )
)

# The GDH primary route at shell_radius 14, cluster 1: the tunnel the P5
# acceptance table is built on. Excluded in that config: HOH WAT BGC NAP.
TUN = "out_GDH_primary_sr14/data/clusters/tun_cl_001.pdb"
EXCLUDED = ("HOH", "WAT", "BGC", "NAP")

needs_refs = pytest.mark.skipif(
    not ((P5_ROOT / TUN).is_file() and (P6_ROOT / TUN).is_file() and STRUCTURE.is_file()),
    reason="reference CAVER outputs for both builds are not available here",
)


@pytest.fixture(scope="module")
def protein():
    return read_pdb(STRUCTURE).exclude_resnames(EXCLUDED).xyz


@pytest.fixture(scope="module")
def raw():
    return {
        "p5": Tunnel.from_caver_cluster(P5_ROOT / TUN),
        "p6": Tunnel.from_caver_cluster(P6_ROOT / TUN),
    }


@needs_refs
class TestTheProblem:
    """The raw chains disagree. If this stops failing, the premise changed."""

    def test_the_builds_sample_the_path_at_very_different_densities(self, raw):
        assert len(raw["p5"]) == 16
        assert len(raw["p6"]) == 145
        assert raw["p5"].mean_spacing == pytest.approx(0.98, abs=0.05)
        assert raw["p6"].mean_spacing == pytest.approx(0.17, abs=0.05)

    def test_the_bottleneck_agrees_because_it_is_interior(self, raw):
        # The shell parameters have no influence deep inside, so this is the
        # one geometric number that is quotable across builds untreated.
        assert raw["p5"].bottleneck == pytest.approx(raw["p6"].bottleneck, abs=0.01)

    def test_the_raw_exits_disagree_by_many_angstroms(self, raw):
        assert np.linalg.norm(raw["p5"].exit - raw["p6"].exit) > 7.0

    def test_the_raw_oic_arm_is_four_times_shorter_on_the_dense_chain(self, raw):
        arm5 = float(np.linalg.norm(raw["p5"].xyz[-1] - raw["p5"].xyz[-2]))
        arm6 = float(np.linalg.norm(raw["p6"].xyz[-1] - raw["p6"].xyz[-2]))
        assert arm5 == pytest.approx(0.91, abs=0.05)
        assert arm6 == pytest.approx(0.26, abs=0.05)
        assert arm5 / arm6 > 3.0

    def test_it_is_nonetheless_the_same_tunnel(self, raw, protein):
        # The dense chain passes through the sparse chain's exit, with the
        # same burial there. The routes agree; only the truncation differs.
        d = np.linalg.norm(raw["p6"].xyz - raw["p5"].exit, axis=1)
        assert d.min() < 0.15
        nearest = raw["p6"].xyz[int(d.argmin())]
        from channelbridge.tunnel import burial_at

        assert int(burial_at(nearest, protein)[0]) == pytest.approx(
            raw["p5"].burial_exit(protein), abs=1
        )


@needs_refs
class TestTheFix:
    """Normalised, the two builds agree to within a sphere."""

    @staticmethod
    @pytest.fixture(scope="class")
    def norm(raw, protein):
        return {k: v.normalised(protein, spacing=1.0) for k, v in raw.items()}

    def test_the_recorded_p5_burial_is_reproduced(self, raw, protein):
        # The P5 acceptance table records burial_exit 21 for this cluster.
        assert raw["p5"].burial_exit(protein) == pytest.approx(21, abs=1)

    def test_the_normalised_exits_agree_within_two_angstroms(self, norm):
        assert np.linalg.norm(norm["p5"].exit - norm["p6"].exit) < 2.0

    def test_the_normalised_exits_sit_at_a_real_surface_mouth(self, norm, protein):
        for k, t in norm.items():
            b = t.burial_exit(protein)
            assert b >= BURIAL_SURFACE_MIN, "%s exit is past the surface (burial %d)" % (k, b)
            assert b < 40, "%s exit is recessed (burial %d)" % (k, b)

    def test_the_normalised_spans_agree_within_two_angstroms(self, norm):
        assert abs(norm["p5"].span - norm["p6"].span) < 2.0

    def test_both_give_a_full_length_oic_arm(self, norm):
        # Measured as arc length; the chord is ~1% shorter on a curved mouth.
        for t in norm.values():
            assert t.arm_chord(1.0) == pytest.approx(1.0, rel=0.03)
            # Chord <= arc length, give or take floating point on a straight mouth.
            assert t.arm_chord(1.0) <= 1.0 + 1e-9

    def test_the_mouth_direction_keeps_a_real_cross_build_uncertainty(self, norm):
        # Normalisation fixes D_IC but does NOT make O_IC build-exact: this
        # tunnel curves as it opens, so the direction depends on the length
        # scale. Recorded here so the number cannot drift unnoticed. At the
        # published coefficients 15.5 deg is worth ~2.3 A of D_IC, ~0.019 DO.
        a5, a6 = norm["p5"].axis(1.0), norm["p6"].axis(1.0)
        cos = float(np.clip(np.dot(a5, a6), -1.0, 1.0))
        assert np.degrees(np.arccos(cos)) == pytest.approx(15.5, abs=2.0)

    def test_a_longer_arm_agrees_better_on_this_curved_mouth(self, norm):
        def between(arm):
            cos = float(np.clip(np.dot(norm["p5"].axis(arm), norm["p6"].axis(arm)), -1.0, 1.0))
            return float(np.degrees(np.arccos(cos)))

        assert between(8.0) < between(1.0)
        # And the mid range is worse than either end, so "longer is better" is
        # not a rule to apply blindly -- hence axis_spread as a quality flag.
        assert between(2.0) > between(1.0)

    def test_the_arm_spread_flag_separates_a_curved_mouth_from_a_straight_one(self, norm):
        # This is what drives the cross-build O_IC difference: the sparse
        # build's chain curves over its last 2 A, the dense build's does not.
        # The flag says which tunnel's O_IC is scale-sensitive, so a reader
        # knows whether to trust the angle to a degree or to fifteen.
        curved = norm["p5"].axis_spread((0.5, 1.0, 2.0))
        straight = norm["p6"].axis_spread((0.5, 1.0, 2.0))
        assert curved > 10.0
        assert straight < 1.0
        assert norm["p5"].describe()["oic_arm_spread_deg"] == pytest.approx(curved, abs=0.01)

    def test_normalisation_leaves_the_bottleneck_alone(self, raw, norm):
        # Resampling must not invent a narrower or wider bottleneck.
        for k in raw:
            assert norm[k].bottleneck == pytest.approx(raw[k].bottleneck, abs=0.15)

    def test_the_dic_difference_between_builds_collapses(self, raw, norm):
        # What this is all for: the same downstream exit scored against both
        # builds must give the same D_IC. A fixed probe point stands in for it.
        probe = np.array([20.0, 30.0, 10.0])
        raw_gap = abs(
            float(np.linalg.norm(raw["p5"].exit - probe) - np.linalg.norm(raw["p6"].exit - probe))
        )
        norm_gap = abs(
            float(np.linalg.norm(norm["p5"].exit - probe) - np.linalg.norm(norm["p6"].exit - probe))
        )
        assert raw_gap > 3.0
        assert norm_gap < 1.5
        assert norm_gap < raw_gap / 2.0
