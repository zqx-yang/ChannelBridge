"""Tunnels fetched from ChannelsDB 2.0 instead of computed with CAVER."""

from __future__ import annotations

import json

import numpy as np
import pytest

from channelbridge.channelsdb import (
    CATEGORIES,
    DbTunnel,
    compare,
    read_report,
    to_tunnel,
    tunnels,
)
from channelbridge.tunnel import Tunnel


def profile(n=5, r=2.0):
    return [{"X": float(i), "Y": 0.0, "Z": 0.0, "Radius": r + 0.1 * i,
             "FreeRadius": r + 0.2 + 0.1 * i, "Distance": float(i)}
            for i in range(n)]


def report(**channels):
    base = {k + "_MOLE": [] for k in ("CSATunnels", "ReviewedChannels",
                                      "CofactorTunnels", "ProcognateTunnels")}
    base.update({k + "_Caver": [] for k in ("CSATunnels", "ReviewedChannels",
                                            "CofactorTunnels", "ProcognateTunnels")})
    base.update(channels)
    return {"Annotations": [], "Channels": base}


class TestParsing:
    def test_it_reads_a_tunnel(self, tmp_path):
        p = tmp_path / "r.json"
        p.write_text(json.dumps(report(CofactorTunnels_MOLE=[
            {"Id": "15.89", "Properties": {"Charge": 1}, "Profile": profile()}])))
        got = tunnels(read_report(str(p)))
        assert len(got) == 1
        assert got[0].method == "MOLE"
        assert got[0].category == "CofactorTunnels"
        assert got[0].label == "Cofactor_MOLE_15.89"

    def test_the_substrate_channel_comes_first(self):
        # Procognate = the cognate ligand's site, which is usually what this
        # pipeline wants; CSA tunnels last
        got = tunnels(report(
            CSATunnels_MOLE=[{"Id": "a", "Profile": profile()}],
            ProcognateTunnels_Caver=[{"Id": "b", "Profile": profile()}]))
        assert [t.category for t in got] == ["ProcognateTunnels", "CSATunnels"]

    def test_both_methods_by_default(self):
        got = tunnels(report(
            CofactorTunnels_MOLE=[{"Id": "a", "Profile": profile()}],
            CofactorTunnels_Caver=[{"Id": "b", "Profile": profile()}]))
        assert {t.method for t in got} == {"MOLE", "Caver"}

    def test_a_method_filter_works(self):
        got = tunnels(report(
            CofactorTunnels_MOLE=[{"Id": "a", "Profile": profile()}],
            CofactorTunnels_Caver=[{"Id": "b", "Profile": profile()}]),
            method="Caver")
        assert [t.method for t in got] == ["Caver"]

    def test_a_one_point_profile_is_skipped_not_crashed_on(self):
        assert tunnels(report(
            CofactorTunnels_MOLE=[{"Id": "a", "Profile": profile(1)}])) == []

    def test_an_empty_report_gives_nothing(self):
        assert tunnels(report()) == []

    def test_every_category_is_known(self):
        assert set(CATEGORIES) == {"ProcognateTunnels", "CofactorTunnels",
                                   "ReviewedChannels", "CSATunnels"}


class TestGeometry:
    def test_length_and_bottleneck(self):
        t = tunnels(report(CofactorTunnels_MOLE=[
            {"Id": "a", "Profile": profile(n=5, r=2.0)}]))[0]
        assert t.length == pytest.approx(4.0)
        assert t.bottleneck == pytest.approx(2.0)

    def test_free_radius_is_kept_separately(self):
        # the two differ by ~0.1-0.2 A and only one is what CAVER's
        # bottleneck means; keeping both lets a caller check
        t = tunnels(report(CofactorTunnels_MOLE=[
            {"Id": "a", "Profile": profile()}]))[0]
        assert float(t.free_radii.min()) > t.bottleneck

    def test_conversion_keeps_the_coordinates(self):
        db = tunnels(report(CofactorTunnels_MOLE=[
            {"Id": "a", "Profile": profile()}]))[0]
        tun = to_tunnel(db)
        assert np.allclose(tun.xyz, db.xyz)
        assert np.allclose(tun.radii, db.radii)
        assert "ChannelsDB" in tun.source

    def test_free_radius_can_be_used_instead(self):
        db = tunnels(report(CofactorTunnels_MOLE=[
            {"Id": "a", "Profile": profile()}]))[0]
        assert np.allclose(to_tunnel(db, use_free_radius=True).radii, db.free_radii)


class TestMatching:
    """Picking WHICH of our clusters corresponds to a database tunnel.

    Closest approach cannot do this. Every CAVER cluster of one run starts
    from the same starting point, so they all share their first spheres and
    the minimum distance between any of them and a database tunnel passing
    that region is the same number. Measured on the reproduction set: clusters 1, 2 and 3 all
    reported 0.6 A to the same tunnel. What separates clusters is the EXIT.
    """

    def _db(self, pts):
        from channelbridge.channelsdb import DbTunnel

        pts = np.asarray(pts, float)
        return DbTunnel(category="CofactorTunnels", method="Caver", tunnel_id="t",
                        xyz=pts, radii=np.full(len(pts), 2.0),
                        free_radii=np.full(len(pts), 2.1),
                        distance=np.arange(len(pts), dtype=float), properties={})

    def _ours(self, pts):
        from channelbridge.tunnel import Tunnel

        pts = np.asarray(pts, float)
        return Tunnel(xyz=pts, radii=np.full(len(pts), 2.0), source="test")

    def test_the_exit_is_what_discriminates(self):
        from channelbridge.channelsdb import compare

        # two of ours sharing an origin and parting at the end -- the case
        # that defeated closest-approach
        common = [[0.0, 0, 0], [1, 0, 0], [2, 0, 0]]
        a = self._ours([*common, [3.0, 0, 0], [4, 0, 0]])
        b = self._ours([*common, [2.5, 3, 0], [3, 8, 0]])
        theirs = self._db([*common, [3.0, 0, 0], [4, 0, 0]])
        ma, mb = compare(a, theirs), compare(b, theirs)
        assert ma.exit_gap < 0.01, "the one that follows it should match"
        assert mb.exit_gap > 5.0, "the one that diverges should not"

    def test_coverage_catches_a_different_route_to_the_same_exit(self):
        from channelbridge.channelsdb import compare

        # long enough that the shared endpoints cannot cover half the chain
        # on their own -- at 11 points they did, which is the test being
        # wrong rather than the code
        straight = self._ours([[float(i), 0, 0] for i in range(21)])
        detour = self._db([[0.0, 0, 0], *[[float(i), 9.0, 0] for i in range(1, 20)],
                           [20.0, 0, 0]])
        m = compare(straight, detour, tol=2.0)
        assert m.exit_gap < 0.01
        assert m.coverage < 0.5
        assert m.verdict() == "same exit, different route"

    def test_a_clean_match_is_called_one(self):
        from channelbridge.channelsdb import compare

        pts = [[float(i), 0, 0] for i in range(8)]
        m = compare(self._ours(pts), self._db(pts))
        assert m.exit_gap == pytest.approx(0.0)
        assert m.coverage == pytest.approx(1.0)
        assert m.verdict() == "same tunnel"

    def test_a_transform_is_applied_to_THEIRS(self):
        from channelbridge.channelsdb import compare

        pts = [[float(i), 0, 0] for i in range(8)]
        shifted = [[p[0] + 100.0, p[1], p[2]] for p in pts]
        R, t = np.eye(3), np.array([-100.0, 0.0, 0.0])
        m = compare(self._ours(pts), self._db(shifted), R, t)
        assert m.exit_gap == pytest.approx(0.0)

    def test_without_a_transform_it_is_nonsense_and_says_so(self):
        # an unsuperposed comparison is not detectable from the numbers
        # alone; it just looks like disagreement. This pins the behaviour.
        from channelbridge.channelsdb import compare

        pts = [[float(i), 0, 0] for i in range(8)]
        shifted = [[p[0] + 100.0, p[1], p[2]] for p in pts]
        assert compare(self._ours(pts), self._db(shifted)).exit_gap > 50

    def test_best_match_picks_by_exit(self):
        from channelbridge.channelsdb import best_match

        ours = self._ours([[float(i), 0, 0] for i in range(6)])
        near = self._db([[float(i), 0, 0] for i in range(6)])
        far = self._db([[0.0, 0, 0], [1, 0, 0], [2, 4, 0], [3, 20, 0]])
        assert best_match(ours, [far, near]).exit_gap < 0.01

    def test_best_match_of_nothing_is_none(self):
        from channelbridge.channelsdb import best_match

        assert best_match(self._ours([[0.0, 0, 0], [1, 0, 0]]), []) is None


class TestLengthIsMeasuredNotTrusted:
    """``Distance`` means different things per method, so length is measured.

    MOLE writes cumulative arc length. The CAVER-derived entries write a
    non-monotonic series that saturates near 1, and reading that as a length
    reported every one of them as about 1 A long -- a plausible number that
    nothing contradicted.
    """

    @staticmethod
    def _tunnel(distance, method="Caver"):
        # a straight chain of ten spheres 1 A apart: the length is 9.0
        pts = np.array([[float(i), 0.0, 0.0] for i in range(10)])
        return DbTunnel(category="ProcognateTunnels", method=method,
                        tunnel_id="t", xyz=pts, radii=np.full(10, 2.0),
                        free_radii=np.full(10, 2.0),
                        distance=np.asarray(distance, dtype=float),
                        properties={})

    def test_length_comes_from_the_coordinates(self):
        # a Distance series claiming 1 A does not change the geometry
        caver_like = [0.0, 0.83, 0.73, 0.91, 0.88, 0.95, 1.0, 1.0, 1.0, 1.0002]
        assert self._tunnel(caver_like).length == pytest.approx(9.0)

    def test_the_reported_series_is_still_available(self):
        caver_like = [0.0, 0.83, 0.73, 0.91, 0.88, 0.95, 1.0, 1.0, 1.0, 1.0002]
        assert self._tunnel(caver_like).reported_length == pytest.approx(1.0002)

    def test_a_non_monotonic_series_is_not_arc_length(self):
        caver_like = [0.0, 0.83, 0.73, 0.91, 0.88, 0.95, 1.0, 1.0, 1.0, 1.0002]
        assert self._tunnel(caver_like).distance_is_arc_length is False

    def test_a_cumulative_series_is_arc_length(self):
        mole_like = [float(i) for i in range(10)]
        assert self._tunnel(mole_like, method="MOLE").distance_is_arc_length is True

    def test_monotonic_but_wrong_scale_is_not_arc_length(self):
        # rises steadily and still disagrees with the coordinates: the
        # monotonicity test alone would pass this
        assert self._tunnel([i * 0.1 for i in range(10)]).distance_is_arc_length is False

    def test_a_single_sphere_has_no_length(self):
        pts = np.array([[0.0, 0.0, 0.0]])
        one = DbTunnel("ProcognateTunnels", "Caver", "t", pts, np.full(1, 2.0),
                       np.full(1, 2.0), np.zeros(1), {})
        assert one.length == 0.0
        assert one.distance_is_arc_length is False


class TestRecoveryVsVerdict:
    """The two questions a comparison can answer, and why one metric is not
    enough for both.

    `verdict` asks "are these the same tunnel", exit-led, for chains of
    comparable extent. `recovery` asks "did our chain find their route",
    for the case where ours is deliberately longer because theirs stopped
    before the surface. Judging the second by the first called recovered
    routes "different tunnel".
    """

    @staticmethod
    def _pair(ours_n, theirs_n):
        # both along the x axis from the origin, 1 A apart: theirs is a
        # truncated prefix of ours, which is exactly the seeded-CAVER case
        ours = Tunnel(xyz=np.array([[float(i), 0.0, 0.0] for i in range(ours_n)]),
                      radii=np.full(ours_n, 2.0), source="ours")
        theirs = DbTunnel("ProcognateTunnels", "MOLE", "t",
                          np.array([[float(i), 0.0, 0.0] for i in range(theirs_n)]),
                          np.full(theirs_n, 2.0), np.full(theirs_n, 2.0),
                          np.arange(float(theirs_n)), {})
        return ours, theirs

    def test_a_truncated_prefix_is_fully_recovered(self):
        m = compare(*self._pair(40, 10))
        assert m.db_coverage == pytest.approx(1.0)
        assert m.recovery() == "route recovered"

    def test_and_the_exit_led_verdict_calls_it_different(self):
        # the point of having both: ours runs 30 A further, so the exits
        # cannot coincide and most of OUR spheres are past their chain
        m = compare(*self._pair(40, 10))
        assert m.exit_gap > 3.0
        assert m.coverage < 0.5
        assert m.verdict() == "different tunnel"

    def test_the_two_coverages_are_the_two_directions(self):
        m = compare(*self._pair(40, 10))
        # their chain ends at x=9, and tol is 2 A, so our spheres at x=0..11
        # count: 12 of 40. All 10 of theirs lie on ours.
        assert m.coverage == pytest.approx(12 / 40)
        assert m.db_coverage == pytest.approx(1.0)

    def test_a_genuinely_different_route_is_not_recovered(self):
        ours, theirs = self._pair(20, 10)
        theirs.xyz = theirs.xyz + np.array([0.0, 30.0, 0.0])
        m = compare(ours, theirs)
        assert m.db_coverage == pytest.approx(0.0)
        assert m.recovery() == "different route"

    def test_half_the_chain_is_partial(self):
        ours, theirs = self._pair(20, 10)
        # move the outer half of their chain away
        theirs.xyz = theirs.xyz.copy()
        theirs.xyz[5:] += np.array([0.0, 30.0, 0.0])
        m = compare(ours, theirs)
        assert m.db_coverage == pytest.approx(0.5)
        assert m.recovery() == "route partly recovered"

    def test_equal_chains_agree_both_ways(self):
        m = compare(*self._pair(10, 10))
        assert m.verdict() == "same tunnel"
        assert m.recovery() == "route recovered"
