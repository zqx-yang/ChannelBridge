"""The DO equation and the angles it is built on."""

from __future__ import annotations

import math

import numpy as np
import pytest

from channelbridge.geometry import (
    DO_PUBLISHED,
    DO_ST2A,
    angle_deg,
    d_ic,
    dic_equivalent_of_oic,
    do_score,
    geometry_from_exits,
    o_ic,
)


class TestDOEquation:
    def test_published_coefficients(self):
        assert DO_PUBLISHED.a == -0.0079
        assert DO_PUBLISHED.b == 0.0012
        assert DO_PUBLISHED.c == 0.7865

    def test_reproduces_the_recorded_p5_geom001(self):
        # The earlier reference: D_IC 16.0 A, O_IC 150 deg -> DO 0.8401.
        assert do_score(16.0, 150.0) == pytest.approx(0.8401, abs=5e-5)

    def test_reproduces_the_recorded_runner_up(self):
        # cl1 FNR: D_IC 22, O_IC 180 -> DO 0.829.
        assert do_score(22.0, 180.0) == pytest.approx(0.829, abs=5e-4)

    def test_one_angstrom_is_worth_six_point_six_degrees(self):
        # The trade-off quoted throughout the project documentation.
        assert DO_PUBLISHED.degrees_per_angstrom == pytest.approx(6.58, abs=0.01)
        assert dic_equivalent_of_oic(6.58) == pytest.approx(1.0, abs=0.01)

    def test_the_two_published_constant_sets_diverge_measurably(self):
        # The iMARS Figure 2C/ST2G and ST2A constants disagree. Over the grid
        # this pipeline actually scans (D_IC 6-36 A, O_IC 90-180 deg) the worst
        # case is 0.032 at the short-and-straight corner, which is larger than
        # the 0.024 quoted in the earlier notes -- so the choice of constants is
        # recorded with every DO rather than assumed negligible.
        per_pose = {
            (d, o): abs(do_score(d, o, DO_PUBLISHED) - do_score(d, o, DO_ST2A))
            for d in range(6, 37)
            for o in range(90, 181, 5)
        }
        worst_at, worst = max(per_pose.items(), key=lambda kv: kv[1])
        assert worst == pytest.approx(0.0321, abs=5e-4)
        assert worst_at == (6, 180)
        # And it is never negligible at the DO resolution the project uses
        # (differences of 0.001 are treated as meaningful).
        assert min(per_pose.values()) > 0.001

    def test_shorter_and_straighter_scores_higher(self):
        assert do_score(10.0, 180.0) > do_score(20.0, 180.0)
        assert do_score(10.0, 180.0) > do_score(10.0, 90.0)

    @pytest.mark.parametrize("d,o", [(-1.0, 150.0), (10.0, 181.0), (10.0, -1.0)])
    def test_rejects_impossible_input(self, d, o):
        with pytest.raises(ValueError):
            do_score(d, o)

    def test_rejects_nan(self):
        with pytest.raises(ValueError):
            do_score(float("nan"), 150.0)

    def test_does_not_clamp_out_of_range_results(self):
        # DO is a linear fit; a value outside [0, 1] is information.
        assert do_score(200.0, 0.0) < 0.0


class TestDistanceAndAngle:
    def test_d_ic_is_a_plain_distance(self):
        assert d_ic([0, 0, 0], [3, 4, 0]) == pytest.approx(5.0)

    def test_angle_at_the_vertex_not_at_the_ends(self):
        # a-b-c with a straight line through b is 180 deg, not 0.
        assert angle_deg([-1, 0, 0], [0, 0, 0], [1, 0, 0]) == pytest.approx(180.0)
        assert angle_deg([1, 0, 0], [0, 0, 0], [0, 1, 0]) == pytest.approx(90.0)

    def test_o_ic_180_means_aimed_straight_at_the_downstream_exit(self):
        # P1 behind P2, P3 straight ahead.
        assert o_ic([-1, 0, 0], [0, 0, 0], [10, 0, 0]) == pytest.approx(180.0)

    def test_o_ic_is_one_sided(self):
        # Swapping which enzyme is upstream changes the angle, by design.
        p1, p2, p3 = [-1, 0, 0], [0, 0, 0], [5, 5, 0]
        forward = o_ic(p1, p2, p3)
        backward = o_ic(p1, p3, p2)
        assert not math.isclose(forward, backward, abs_tol=1.0)

    def test_zero_length_arm_raises_rather_than_returning_a_number(self):
        # This is the densely-sampled-tunnel failure: two coincident spheres.
        with pytest.raises(ValueError, match="arm"):
            angle_deg([0, 0, 0], [0, 0, 0], [1, 0, 0])

    def test_zero_length_arm_can_be_asked_for_nan_instead(self):
        assert math.isnan(angle_deg([0, 0, 0], [0, 0, 0], [1, 0, 0], degenerate="nan"))

    def test_numerically_safe_at_exactly_antiparallel(self):
        # acos of a value a hair outside [-1, 1] would raise without clamping.
        v = np.array([1.0, 0.0, 0.0])
        assert angle_deg(-v * 1e8, [0, 0, 0], v * 1e8) == pytest.approx(180.0)

    @pytest.mark.parametrize("bad", [[0, 0], [0, 0, 0, 0], [[0, 0, 0], [1, 1, 1]]])
    def test_rejects_non_3_vectors(self, bad):
        with pytest.raises(ValueError):
            d_ic(bad, [0, 0, 0])


class TestGeometryFromExits:
    def test_returns_the_triple_consistently(self):
        d, o, do = geometry_from_exits([-1, 0, 0], [0, 0, 0], [16, 0, 0])
        assert d == pytest.approx(16.0)
        assert o == pytest.approx(180.0)
        assert do == pytest.approx(do_score(16.0, 180.0))


class TestOicIsOneSided:
    """O_IC measures the DONOR channel's aim and nothing else.

    Measured on the worked GDH-FNR pose: O_IC 165 deg at the upstream exit, but
    the same angle at the downstream exit is 90 deg -- the receiving mouth
    opens sideways to the line joining the two exits. DO cannot see that, so
    two geometries that differ entirely in whether the acceptor faces the
    donor score identically. The map reports it; it must never enter the
    score, whose published coefficients cover D_IC and O_IC only.
    """

    def test_the_score_ignores_the_downstream_orientation(self):
        from channelbridge.geometry import do_score

        # identical D_IC and O_IC; the acceptor's own aim is not an argument
        # the function even accepts, which is the point
        assert do_score(17.0, 165.0) == do_score(17.0, 165.0)

    def test_the_angle_is_computed_at_the_upstream_exit(self):
        import numpy as np

        from channelbridge.geometry import o_ic

        # P1 -> P2 heading +x, P3 straight ahead: perfectly aimed, 180 deg
        p1, p2, p3 = np.array([-1.0, 0, 0]), np.zeros(3), np.array([17.0, 0, 0])
        assert abs(o_ic(p1, p2, p3) - 180.0) < 1e-6
        # move P3 sideways and the angle falls, because the DONOR now misses
        p3b = np.array([0.0, 17.0, 0.0])
        assert abs(o_ic(p1, p2, p3b) - 90.0) < 1e-6


class TestYourOwnScore:
    """The published DO is the default, not the only option.

    Coefficients were always replaceable. EXTRA TERMS were not, which left
    no way to score on anything the published fit ignores -- and it ignores
    the acceptor's orientation entirely (see TestOicIsOneSided).
    """

    def _mine(self):
        from channelbridge.geometry import DOCoefficients

        return DOCoefficients(-0.0079, 0.0012, 0.7865, name="DO+acceptor",
                              extra={"o_ic_rev": 0.0012})

    def test_the_default_is_the_published_score(self):
        from channelbridge.geometry import DO_PUBLISHED, do_score

        assert DO_PUBLISHED.is_published
        assert do_score(17.0, 165.0) == pytest.approx(0.8502, abs=5e-5)

    def test_an_extra_term_is_applied(self):
        from channelbridge.geometry import do_score

        m = self._mine()
        facing = do_score(17.0, 165.0, m, o_ic_rev=180.0)
        sideways = do_score(17.0, 165.0, m, o_ic_rev=90.0)
        assert facing > sideways, "an acceptor facing the donor should score higher"
        assert facing - sideways == pytest.approx(0.0012 * 90.0)

    def test_adding_a_term_means_it_is_no_longer_iMARS(self):
        # the flag exists so a four-term number cannot be reported as "DO"
        # without someone having to notice
        assert not self._mine().is_published

    def test_it_says_which_features_it_needs(self):
        assert self._mine().features == ("o_ic_rev",)

    def test_a_missing_feature_is_refused_not_defaulted(self):
        from channelbridge.geometry import do_score

        with pytest.raises(KeyError, match="o_ic_rev"):
            do_score(17.0, 165.0, self._mine())

    def test_features_on_a_two_term_set_are_refused(self):
        # silently ignoring them would let someone believe a term was applied
        from channelbridge.geometry import do_score

        with pytest.raises(TypeError, match="two-term"):
            do_score(17.0, 165.0, o_ic_rev=180.0)

    def test_the_published_sets_are_unchanged_by_all_this(self):
        from channelbridge.geometry import DO_EXACT, DO_PUBLISHED

        for c in (DO_PUBLISHED, DO_EXACT):
            assert c.extra == {} and c.is_published


class TestSuperpose:
    """Shared by the occupancy copies and the partner footprint."""

    def _cloud(self, n=40, seed=0):
        return np.random.default_rng(seed).normal(scale=10.0, size=(n, 3))

    def test_recovers_a_known_transform(self):
        from channelbridge.geometry import superpose
        from channelbridge.reference_pose import rotation_about

        p = self._cloud()
        rot = rotation_about([0.3, 0.5, 0.8], 47.0)
        trans = np.array([3.0, -2.0, 7.0])
        back, shift, rmsd = superpose(p, p @ rot.T + trans)
        assert back == pytest.approx(rot, abs=1e-9)
        assert shift == pytest.approx(trans, abs=1e-9)
        assert rmsd == pytest.approx(0.0, abs=1e-9)

    def test_never_returns_a_reflection(self):
        """A mirror image is not a placement anything can be built from."""
        from channelbridge.geometry import superpose

        p = self._cloud(seed=1)
        mirrored = p * np.array([1.0, 1.0, -1.0])
        rot, _t, _rmsd = superpose(p, mirrored)
        assert float(np.linalg.det(rot)) == pytest.approx(1.0, abs=1e-9)

    def test_rmsd_is_reported_not_hidden(self):
        from channelbridge.geometry import superpose

        p = self._cloud(seed=2)
        noisy = p + np.random.default_rng(3).normal(scale=0.5, size=p.shape)
        _r, _t, rmsd = superpose(p, noisy)
        assert 0.1 < rmsd < 2.0

    def test_too_few_points_to_fix_a_rotation(self):
        from channelbridge.geometry import superpose

        with pytest.raises(ValueError, match="at least 3"):
            superpose(np.zeros((2, 3)), np.ones((2, 3)))

    def test_mismatched_shapes_are_refused(self):
        from channelbridge.geometry import superpose

        with pytest.raises(ValueError, match="matching shapes"):
            superpose(np.zeros((5, 3)), np.ones((6, 3)))
