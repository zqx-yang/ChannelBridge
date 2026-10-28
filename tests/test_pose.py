"""Stage 5: the pose search.

The tests that matter here are the invariants. A pose search can be wrong in
a way that looks entirely plausible in its output -- every column populated,
DO falling smoothly with D_IC -- while the geometry it built is not the
geometry it recorded. So the first group checks that what comes out is what
was asked for:

* the realised D_IC is the requested one, exactly;
* the realised O_IC is the requested one, exactly (this is what fails if the
  pose is built along a fitted axis and scored along the O_IC arm, the
  9-12 degree trap documented on ``Tunnel.fitted_axis``);
* ``tilt = 0`` gives ``O_IC_rev = 180``, which is the same trap seen from the
  downstream end;
* the transform is a rotation, not a reflection or a scaling.

The second group checks the filters on geometry whose answer is known by
construction, and the third that the cheap proximity test agrees with the
exact distance -- including inside the band where it refuses to answer.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from channelbridge.geometry import DO_PUBLISHED, do_score, o_ic
from channelbridge.pose import (
    MIN_DIC,
    PatchSite,
    PoseParams,
    PoseSide,
    _Proximity,
    equivalent_chains,
    path_is_clear,
    pose_grid,
    pose_transform,
    read_patch_sites,
    score_pose,
    search,
    segment_cloud_distance,
    segment_segment_distance,
    symmetry_transforms,
    write_patch_pair_best,
    write_patch_pairs,
    write_poses,
)


def ball(centre, radius=8.0, n=400, seed=0):
    """A blob of atoms, so a side has something to clash with."""
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    r = radius * rng.random(n) ** (1 / 3.0)
    return np.asarray(centre, float) + v * r[:, None]


def side(name, centre, exit_point, inner_point, *, radius=8.0, patches=(), seed=0):
    return PoseSide(name=name, xyz=ball(centre, radius, seed=seed),
                    exit=exit_point, axis_point=inner_point,
                    patches=patches, centroid=centre)


def simple_pair(gap=0.0):
    """Two blobs whose channels both point along x, facing each other."""
    up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0], seed=1)
    dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0], seed=2)
    return up, dn


#: ``tilt = 180`` turns the downstream channel to point AWAY from the upstream
#: exit, which swings its body back across the gap and into the other enzyme.
#: That is how these tests produce a real clash: at tilt 0 the two bodies lie
#: on opposite sides of their exits, so even D_IC 1 A leaves them tangent
#: rather than overlapping, which is a fact about the parameterisation and not
#: something to assert around.
BACKWARDS = 180.0


class TestRequestedIsRealised:
    """What the grid asks for has to be what comes out."""

    @pytest.mark.parametrize("d", [6.0, 12.5, 30.0])
    @pytest.mark.parametrize("angle", [180.0, 150.0, 120.0, 90.0])
    def test_d_ic_and_o_ic_are_exact(self, d, angle):
        from channelbridge.pose import GridPoint

        up, dn = simple_pair()
        rot, trans = pose_transform(up, dn, GridPoint(d, angle, 45.0, 30.0, 0.0, 70.0))
        row = score_pose(up, dn, rot, trans)
        assert row["D_IC"] == pytest.approx(d, abs=1e-9)
        assert row["O_IC"] == pytest.approx(angle, abs=1e-6)

    def test_o_ic_rev_is_180_at_tilt_zero(self):
        """Aimed straight back is the definition of tilt 0."""
        from channelbridge.pose import GridPoint

        up, dn = simple_pair()
        for angle in (180.0, 135.0, 95.0):
            for spin in (0.0, 90.0, 250.0):
                rot, trans = pose_transform(
                    up, dn, GridPoint(14.0, angle, 20.0, 0.0, 0.0, spin))
                row = score_pose(up, dn, rot, trans)
                assert row["O_IC_rev"] == pytest.approx(180.0, abs=1e-6)

    def test_tilt_moves_o_ic_rev_by_exactly_the_tilt(self):
        from channelbridge.pose import GridPoint

        up, dn = simple_pair()
        for tilt in (15.0, 45.0, 90.0):
            rot, trans = pose_transform(
                up, dn, GridPoint(20.0, 180.0, 0.0, tilt, 0.0, 0.0))
            row = score_pose(up, dn, rot, trans)
            assert row["O_IC_rev"] == pytest.approx(180.0 - tilt, abs=1e-6)

    def test_spin_does_not_change_the_score(self):
        """Spin picks the epitope, not the geometry the metric sees."""
        from channelbridge.pose import GridPoint

        up, dn = simple_pair()
        base = None
        for spin in (0.0, 30.0, 120.0, 330.0):
            rot, trans = pose_transform(
                up, dn, GridPoint(18.0, 150.0, 60.0, 30.0, 0.0, spin))
            row = score_pose(up, dn, rot, trans)
            now = (row["D_IC"], row["O_IC"], row["O_IC_rev"], row["DO"])
            if base is None:
                base = now
            else:
                assert now == pytest.approx(base, abs=1e-6)

    def test_transform_is_a_rotation(self):
        from channelbridge.pose import GridPoint

        up, dn = simple_pair()
        rot, _t = pose_transform(up, dn, GridPoint(11.0, 140.0, 30.0, 60.0, 90.0, 210.0))
        assert rot @ rot.T == pytest.approx(np.eye(3), abs=1e-9)
        assert float(np.linalg.det(rot)) == pytest.approx(1.0, abs=1e-9)

    def test_tilt_azimuth_moves_the_enzyme_but_not_the_score(self):
        """The dimension the default grid leaves out still has to work."""
        from channelbridge.pose import GridPoint

        up, dn = simple_pair()
        seen = []
        for taz in (0.0, 90.0, 180.0):
            rot, trans = pose_transform(
                up, dn, GridPoint(20.0, 180.0, 0.0, 45.0, taz, 0.0))
            row = score_pose(up, dn, rot, trans)
            assert row["D_IC"] == pytest.approx(20.0, abs=1e-9)
            assert row["O_IC"] == pytest.approx(180.0, abs=1e-6)
            assert row["O_IC_rev"] == pytest.approx(135.0, abs=1e-6)
            seen.append(rot @ np.array([1.0, 0.0, 0.0]))
        # same angles, genuinely different placements
        assert np.linalg.norm(seen[0] - seen[1]) > 0.1
        assert np.linalg.norm(seen[0] - seen[2]) > 0.1


class TestGrid:
    def test_ordered_best_do_first(self):
        grid = pose_grid(PoseParams(dic_min=6, dic_max=12, dic_step=2,
                                    oic=(180, 120), azimuth=(0, 180),
                                    tilt=(0,), spin=(0, 90)))
        dos = [do_score(g.d_ic, g.o_ic) for g in grid]
        assert dos == sorted(dos, reverse=True)

    def test_azimuth_collapses_at_oic_180(self):
        """The O_IC cone is a point there, so twelve azimuths are one pose."""
        params = PoseParams(dic_min=10, dic_max=10, oic=(180.0,),
                            azimuth=(0, 90, 180, 270), tilt=(0,), spin=(0,))
        assert len(pose_grid(params)) == 1

    def test_tilt_azimuth_collapses_at_tilt_zero(self):
        params = PoseParams(dic_min=10, dic_max=10, oic=(150.0,), azimuth=(0,),
                            tilt=(0.0,), tilt_azimuth=(0, 90, 180), spin=(0,))
        assert len(pose_grid(params)) == 1
        params = PoseParams(dic_min=10, dic_max=10, oic=(150.0,), azimuth=(0,),
                            tilt=(0.0, 30.0), tilt_azimuth=(0, 90, 180), spin=(0,))
        assert len(pose_grid(params)) == 1 + 3

    def test_degenerate_dic_is_dropped(self):
        params = PoseParams(dic_min=0.0, dic_max=2.0, dic_step=1.0,
                            oic=(180.0,), azimuth=(0,), tilt=(0,), spin=(0,))
        assert [g.d_ic for g in pose_grid(params)] == [1.0, 2.0]

    def test_all_degenerate_dic_is_refused(self):
        params = PoseParams(dic_min=0.0, dic_max=0.4, dic_step=0.2,
                            oic=(180.0,), azimuth=(0,), tilt=(0,), spin=(0,))
        with pytest.raises(ValueError, match="degenerate"):
            pose_grid(params)

    def test_default_grid_size_is_the_reference_grid(self):
        """31 x 7 x 12 x 4 x 12, less the azimuth collapse at O_IC 180."""
        params = PoseParams()
        assert len(params.dic_values) == 31
        assert len(pose_grid(params)) == 31 * (12 * 4 * 12 * 6 + 1 * 4 * 12)

    def test_bad_params_are_refused_where_they_are_made(self):
        with pytest.raises(ValueError, match="dic_step"):
            PoseParams(dic_step=0.0)
        with pytest.raises(ValueError, match="above dic_max"):
            PoseParams(dic_min=40.0, dic_max=10.0)
        with pytest.raises(ValueError, match=r"\[0, 180\]"):
            PoseParams(oic=(200.0,))
        with pytest.raises(ValueError, match="occupancy"):
            PoseParams(occupancy="yes")
        with pytest.raises(ValueError, match="bridge_lo"):
            PoseParams(bridge_lo=60.0, bridge_hi=50.0)


class TestSegments:
    def test_point_to_segment_is_exact_not_sampled(self):
        """A thin obstruction halfway along must not be stepped over."""
        p, q = np.array([0.0, 0.0, 0.0]), np.array([20.0, 0.0, 0.0])
        obstruction = np.array([[10.5, 0.4, 0.0]])
        assert segment_cloud_distance(p, q, obstruction) == pytest.approx(0.4, abs=1e-9)

    def test_beyond_the_ends_measures_to_the_endpoint(self):
        p, q = np.array([0.0, 0.0, 0.0]), np.array([10.0, 0.0, 0.0])
        far = np.array([[-3.0, 4.0, 0.0]])
        assert segment_cloud_distance(p, q, far) == pytest.approx(5.0, abs=1e-9)

    def test_parameter_range_restricts_what_is_tested(self):
        p, q = np.array([0.0, 0.0, 0.0]), np.array([10.0, 0.0, 0.0])
        at_start = np.array([[0.0, 0.5, 0.0]])
        assert segment_cloud_distance(p, q, at_start) == pytest.approx(0.5)
        # skip the first half and the point is 5 A away along the segment
        assert segment_cloud_distance(p, q, at_start, t_lo=0.5) == pytest.approx(
            math.hypot(5.0, 0.5))

    def test_segment_segment_is_segments_not_lines(self):
        # two parallel, offset segments: lines would touch, segments do not
        a = (np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]))
        b = (np.array([10.0, 0.0, 0.0]), np.array([11.0, 0.0, 0.0]))
        assert segment_segment_distance(a[0], a[1], b[0], b[1]) == pytest.approx(9.0)

    def test_crossing_segments_touch(self):
        a = (np.array([-5.0, 0.0, 0.0]), np.array([5.0, 0.0, 0.0]))
        b = (np.array([0.0, -5.0, 0.0]), np.array([0.0, 5.0, 0.0]))
        assert segment_segment_distance(a[0], a[1], b[0], b[1]) == pytest.approx(0.0)


class TestPathIsClear:
    def test_own_residues_do_not_block_the_path(self):
        """The endpoints are ON the surface, so the margins must be skipped."""
        p, q = np.array([0.0, 0.0, 0.0]), np.array([40.0, 0.0, 0.0])
        near_ends = np.array([[1.0, 0.5, 0.0], [39.0, 0.5, 0.0]])
        assert path_is_clear(p, q, near_ends, clearance=3.5, margin=7.0)

    def test_an_obstruction_in_the_middle_blocks(self):
        p, q = np.array([0.0, 0.0, 0.0]), np.array([40.0, 0.0, 0.0])
        middle = np.array([[20.0, 1.0, 0.0]])
        assert not path_is_clear(p, q, middle, clearance=3.5, margin=7.0)

    def test_short_segment_has_nothing_left_to_test(self):
        p, q = np.array([0.0, 0.0, 0.0]), np.array([5.0, 0.0, 0.0])
        blocked = np.array([[2.5, 0.0, 0.0]])
        assert path_is_clear(p, q, blocked, clearance=3.5, margin=7.0)

    def test_zero_length_is_not_a_path(self):
        p = np.array([1.0, 2.0, 3.0])
        assert not path_is_clear(p, p, np.zeros((1, 3)))


class TestProximity:
    """The cheap test has to agree with the exact one, or refuse to answer."""

    @pytest.mark.parametrize("spacing", [1.0, 0.6])
    def test_never_contradicts_the_exact_distance(self, spacing):
        rng = np.random.default_rng(7)
        cloud = ball([0.0, 0.0, 0.0], 9.0, n=250, seed=3)
        prox = _Proximity(cloud, 3.0, spacing)
        from channelbridge.reference_pose import min_interprotein_distance

        decided = 0
        for _ in range(60):
            probe = cloud[:40] + rng.normal(scale=6.0, size=3)
            verdict = prox.test(probe)
            exact = min_interprotein_distance(cloud, probe)
            if verdict is True:
                assert exact < 3.0
                decided += 1
            elif verdict is False:
                assert exact >= 3.0
                decided += 1
        assert decided > 0, "the test refused to answer anything at all"

    def test_undecided_only_near_the_boundary(self):
        cloud = ball([0.0, 0.0, 0.0], 9.0, n=250, seed=4)
        prox = _Proximity(cloud, 3.0, 1.0)
        from channelbridge.reference_pose import min_interprotein_distance

        slop = 0.5 * math.sqrt(3.0)
        for shift in np.linspace(0.0, 30.0, 61):
            probe = cloud[:40] + np.array([shift, 0.0, 0.0])
            if prox.test(probe) is None:
                exact = min_interprotein_distance(cloud, probe)
                assert abs(exact - 3.0) <= slop + 1e-6


class TestNarrowedDistance:
    """The narrowed distance must equal the full one, or it is just wrong."""

    def test_agrees_with_the_full_computation(self):
        from channelbridge.pose import NEAR_RADIUS, _clearance_grid, _min_distance
        from channelbridge.reference_pose import min_interprotein_distance

        rng = np.random.default_rng(21)
        up = ball([0.0, 0.0, 0.0], 10.0, n=600, seed=5)
        near = _clearance_grid(up, NEAR_RADIUS, 1.0)
        checked = 0
        for _ in range(40):
            moved = ball([0.0, 0.0, 0.0], 8.0, n=300, seed=6) + rng.normal(scale=12.0, size=3)
            fast = _min_distance(up, near, moved)
            slow = min_interprotein_distance(up, moved)
            assert fast == pytest.approx(slow, abs=1e-9)
            checked += 1
        assert checked == 40

    def test_exact_even_when_nothing_is_near(self):
        """Far apart, the narrowing has nothing to give and must not guess."""
        from channelbridge.pose import NEAR_RADIUS, _clearance_grid, _min_distance
        from channelbridge.reference_pose import min_interprotein_distance

        up = ball([0.0, 0.0, 0.0], 6.0, n=200, seed=7)
        near = _clearance_grid(up, NEAR_RADIUS, 1.0)
        moved = ball([0.0, 0.0, 0.0], 6.0, n=200, seed=8) + np.array([90.0, 0.0, 0.0])
        assert _min_distance(up, near, moved) == pytest.approx(
            min_interprotein_distance(up, moved), abs=1e-9)
        assert _min_distance(up, near, moved) > NEAR_RADIUS


class TestSearch:
    def test_two_separated_blobs_find_feasible_poses(self):
        up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0])),
                           PatchSite(2, np.array([0.0, -8.0, 0.0]))], seed=1)
        dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=2)
        params = PoseParams(dic_min=20, dic_max=30, dic_step=5, oic=(180.0,),
                            azimuth=(0.0,), tilt=(0.0,), spin=(0.0, 90.0, 180.0),
                            max_ok=5, bridge_lo=5.0, bridge_hi=60.0)
        result = search(up, dn, params)
        assert result.feasible, result.summary()
        best = result.best()
        assert best.status == "ok"
        assert best.geom.startswith("geom")
        # DO falls with D_IC, so the best feasible pose is the tightest one
        assert best.d_ic == min(p.d_ic for p in result.feasible)

    def test_overlapping_blobs_are_all_clash(self):
        up, dn = simple_pair()
        params = PoseParams(dic_min=1.0, dic_max=2.0, dic_step=1.0, oic=(180.0,),
                            azimuth=(0.0,), tilt=(BACKWARDS,), spin=(0.0,))
        result = search(up, dn, params)
        assert {p.status for p in result.poses} == {"clash"}
        assert result.best() is None
        assert "no feasible pose" in result.summary()

    def test_no_patches_means_no_bridgeable_pair(self):
        """A geometry that scores well and cannot be built is the point."""
        up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0], seed=1)
        dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0], seed=2)
        params = PoseParams(dic_min=25, dic_max=25, oic=(180.0,), azimuth=(0.0,),
                            tilt=(0.0,), spin=(0.0,))
        result = search(up, dn, params)
        assert [p.status for p in result.poses] == ["no_bridgeable_pair"]

    def test_max_ok_stops_the_scan(self):
        up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=1)
        dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=2)
        params = PoseParams(dic_min=22, dic_max=34, dic_step=1.0, oic=(180.0,),
                            azimuth=(0.0,), tilt=(0.0,),
                            spin=tuple(range(0, 360, 30)), max_ok=3,
                            bridge_lo=5.0, bridge_hi=60.0)
        result = search(up, dn, params)
        assert len(result.feasible) == 3
        assert result.truncated
        assert len(result.poses) < result.n_grid

    def test_only_feasible_poses_get_a_geom_id(self):
        """Numbering rejections alongside successes makes the id meaningless."""
        up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=1)
        dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=2)
        params = PoseParams(dic_min=6, dic_max=30, dic_step=6, oic=(180.0,),
                            azimuth=(0.0,), tilt=(0.0, BACKWARDS), spin=(0.0,),
                            tilt_azimuth=(0.0,), bridge_lo=5.0, bridge_hi=60.0,
                            max_ok=0)
        result = search(up, dn, params)
        assert any(p.status == "clash" for p in result.poses)
        assert all(p.geom == "" for p in result.poses if not p.feasible)
        assert [p.geom for p in result.feasible] == \
            ["geom%03d" % i for i in range(1, len(result.feasible) + 1)]

    def test_occupancy_on_without_protomers_is_refused(self):
        up, dn = simple_pair()
        with pytest.raises(ValueError, match="no equivalent protomers"):
            search(up, dn, PoseParams(occupancy="on", dic_min=20, dic_max=20,
                                      oic=(180.0,), azimuth=(0.0,), tilt=(0.0,),
                                      spin=(0.0,)))

    def test_measured_distance_only_where_it_is_needed(self):
        """A pose rejected by the cheap test reports no distance, not a guess."""
        up, dn = simple_pair()
        params = PoseParams(dic_min=1.0, dic_max=1.0, oic=(180.0,), azimuth=(0.0,),
                            tilt=(BACKWARDS,), spin=(0.0,))
        result = search(up, dn, params)
        assert result.poses[0].status == "clash"
        assert math.isnan(result.poses[0].min_distance)


class TestTiltAzimuthNote:
    """The advice has to be derived from the result, not boilerplate."""

    def _run(self, tilts, taz=(0.0,)):
        up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=1)
        dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=2)
        # the full spin sweep: with fewer, no patch pair of these blobs comes
        # into bridgeable range and every pose is no_bridgeable_pair
        return search(up, dn, PoseParams(
            dic_min=24, dic_max=30, dic_step=3, oic=(180.0,), azimuth=(0.0,),
            tilt=tilts, tilt_azimuth=taz, spin=tuple(range(0, 360, 30)),
            bridge_lo=5.0, bridge_hi=70.0, max_ok=0))

    def test_no_tilt_used_means_sweeping_cannot_help(self):
        note = self._run((0.0,)).tilt_azimuth_note()
        assert note and "cannot change this answer" in note[0]

    def test_tilt_used_means_one_direction_was_tried(self):
        result = self._run((0.0, 30.0))
        tilted = [p for p in result.feasible if p.requested.tilt]
        if not tilted:
            pytest.skip("no tilted pose was feasible on this synthetic pair")
        note = result.tilt_azimuth_note()
        assert note and "was not swept" in note[0]

    def test_silent_when_the_axis_was_swept(self):
        assert self._run((0.0, 30.0), taz=(0.0, 90.0, 180.0)).tilt_azimuth_note() == []

    def test_silent_when_nothing_is_feasible(self):
        up, dn = simple_pair()
        result = search(up, dn, PoseParams(
            dic_min=1.0, dic_max=1.0, oic=(180.0,), azimuth=(0.0,),
            tilt=(BACKWARDS,), spin=(0.0,)))
        assert result.tilt_azimuth_note() == []


class TestTagNames:
    def test_a_separator_in_an_enzyme_name_does_not_split_the_tag(self):
        """The figure stage parses `<up>p<N>-<dn>p<M>` out of a file name."""
        from channelbridge.figures import _from_stem

        up = side("enz_4CL", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(4, np.array([0.0, 8.0, 0.0]))], seed=1)
        dn = side("enz_STS", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(6, np.array([0.0, 8.0, 0.0]))], seed=2)
        assert up.tag_name == "enz4CL" and dn.tag_name == "enzSTS"
        result = search(up, dn, PoseParams(
            dic_min=24, dic_max=27, dic_step=3, oic=(180.0,), azimuth=(0.0,),
            tilt=(0.0,), spin=tuple(range(0, 360, 30)), bridge_lo=5.0,
            bridge_hi=70.0, max_ok=0))
        best = result.best()
        if best is None:
            pytest.skip("no feasible pose on this synthetic pair")
        assert best.best_pair_tag == "enz4CLp4-enzSTSp6"
        # the human column keeps the real names
        assert best.best_pair == "enz_4CL_p4 + enz_STS_p6"
        parsed = _from_stem(best.tag)
        assert (parsed.get("patch_up"), parsed.get("patch_dn")) == ("4", "6")

    def test_a_name_with_nothing_usable_still_makes_a_tag(self):
        s = side("___", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0])
        assert s.tag_name == "x"


class TestOccupancy:
    def test_binders_that_overlap_are_rejected(self):
        from channelbridge.pose import occupancy_check

        # two protomers related by a half turn about z, 30 A apart
        rot = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
        syms = [("A", np.eye(3), np.zeros(3)), ("B", rot, np.array([0.0, 0.0, 0.0]))]
        posed = ball([25.0, 0.0, 0.0], 6.0, n=100)
        # a binder axis straight through the origin: its mirror image coincides
        a, b = np.array([5.0, 0.0, 0.0]), np.array([-5.0, 0.0, 0.0])
        ok, why, metrics = occupancy_check(
            posed, a, b, syms, ball([0.0, 0.0, 0.0], 8.0), PoseParams())
        assert not ok
        assert "overlap" in why or "runs through" in why
        assert metrics["min_binder_binder"] < 2 * PoseParams().binder_radius

    def test_single_protomer_is_always_fine(self):
        from channelbridge.pose import occupancy_check

        ok, why, metrics = occupancy_check(
            ball([20.0, 0.0, 0.0]), np.zeros(3), np.ones(3),
            [("A", np.eye(3), np.zeros(3))], ball([0.0, 0.0, 0.0]), PoseParams())
        assert ok and why == "" and metrics == {}


class TestOccupancyInsideTheSearch:
    """The filter has to actually run, and only after the cheaper ones.

    `occupancy_check` being correct is not the same as the search reaching
    it. Here the upstream enzyme has a second protomer only 20 degrees away
    about z, so the two downstream copies sit close together and overlap at
    small D_IC and separate at large -- which puts `occupancy_clash` and
    `ok` in the same scan.
    """

    def _pair(self):
        from channelbridge.pose import rotation_about

        up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=1)
        dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]))], seed=2)
        up.symmetry = [("A", np.eye(3), np.zeros(3)),
                       ("B", rotation_about([0.0, 0.0, 1.0], 20.0), np.zeros(3))]
        return up, dn

    def _scan(self, **kw):
        up, dn = self._pair()
        params = PoseParams(dic_min=14, dic_max=70, dic_step=7, oic=(180.0,),
                            azimuth=(0.0,), tilt=(0.0,),
                            spin=tuple(range(0, 360, 30)), bridge_lo=5.0,
                            bridge_hi=90.0, max_ok=0, **kw)
        return search(up, dn, params)

    def test_the_filter_rejects_poses_the_others_passed(self):
        result = self._scan()
        rejected = [p for p in result.poses if p.status == "occupancy_clash"]
        assert rejected, "the occupancy filter never fired: " + result.summary()
        for p in rejected:
            assert p.occupancy_note
            assert p.min_distance >= result.params.binder_gap
            assert p.n_bridgeable > 0

    def test_turning_it_off_makes_those_poses_feasible(self):
        """Which is the only way to know the filter changed the answer."""
        on, off = self._scan(), self._scan(occupancy="off")
        assert len(off.feasible) > len(on.feasible)
        n_rejected = sum(1 for p in on.poses if p.status == "occupancy_clash")
        assert len(off.feasible) == len(on.feasible) + n_rejected

    def test_the_metrics_are_recorded_for_every_pose_that_reached_it(self):
        result = self._scan()
        reached = [p for p in result.poses
                   if p.status in ("ok", "occupancy_clash")]
        assert reached
        for p in reached:
            assert set(p.occupancy) >= {"min_copy_copy", "min_binder_binder",
                                        "min_binder_copy", "min_copy_up"}

    def test_overlapping_copies_are_why(self):
        rejected = [p for p in self._scan().poses
                    if p.status == "occupancy_clash"]
        assert any("copies overlap" in p.occupancy_note for p in rejected)


class TestSymmetry:
    def test_transforms_recover_a_known_rotation(self):
        from channelbridge.io.pdb import Atom, AtomArray

        rng = np.random.default_rng(11)
        base = rng.normal(scale=10.0, size=(30, 3))
        half = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
        shift = np.array([0.0, 0.0, 7.0])
        atoms = []
        serial = 1
        for chain, xyz in (("A", base), ("B", base @ half.T + shift)):
            for i, p in enumerate(xyz, 1):
                atoms.append(Atom(serial=serial, name=" CA ", altloc="",
                                  resname="ALA", chain=chain, resseq=i, icode=" ",
                                  xyz=tuple(p), occupancy=1.0, bfactor=0.0,
                                  element="C", hetatm=False))
                serial += 1
        arr = AtomArray(atoms)
        assert equivalent_chains(arr, "A") == ["A", "B"]
        syms = symmetry_transforms(arr, "A", ["A", "B"])
        assert len(syms) == 2
        chain, rot, trans = syms[1]
        assert chain == "B"
        assert rot == pytest.approx(half, abs=1e-6)
        assert trans == pytest.approx(shift, abs=1e-6)

    def test_a_different_protein_is_not_a_protomer(self):
        from channelbridge.io.pdb import Atom, AtomArray

        rng = np.random.default_rng(12)
        atoms = []
        serial = 1
        for chain, resname in (("A", "ALA"), ("B", "TRP")):
            for i in range(1, 31):
                atoms.append(Atom(serial=serial, name=" CA ", altloc="",
                                  resname=resname, chain=chain, resseq=i, icode=" ",
                                  xyz=tuple(rng.normal(scale=10.0, size=3)),
                                  occupancy=1.0, bfactor=0.0, element="C",
                                  hetatm=False))
                serial += 1
        assert equivalent_chains(AtomArray(atoms), "A") == ["A"]


class TestReadPatchSites:
    def _structure(self):
        from channelbridge.io.pdb import Atom, AtomArray

        atoms = []
        serial = 1
        for resi in (10, 11, 12):
            for name, offset in ((" CA ", 0.0), (" CB ", 1.5)):
                atoms.append(Atom(serial=serial, name=name, altloc="",
                                  resname="LEU", chain="A", resseq=resi, icode=" ",
                                  xyz=(float(resi), offset, 0.0), occupancy=1.0,
                                  bfactor=0.0, element="C", hetatm=False))
                serial += 1
        # a glycine, which has no CB and must fall back to CA
        atoms.append(Atom(serial=serial, name=" CA ", altloc="", resname="GLY",
                          chain="A", resseq=13, icode=" ", xyz=(13.0, 0.0, 0.0),
                          occupancy=1.0, bfactor=0.0, element="C", hetatm=False))
        return AtomArray(atoms)

    def test_centroid_is_the_mean_cb(self, tmp_path):
        table = tmp_path / "patches.tsv"
        table.write_text("rank\tscore\thotspot_residues\n1\t4.2\t10+11+12\n")
        sites = read_patch_sites(str(table), self._structure())
        assert len(sites) == 1
        assert sites[0].centroid == pytest.approx([11.0, 1.5, 0.0])
        assert sites[0].score == pytest.approx(4.2)
        assert sites[0].label == "p1"

    def test_glycine_falls_back_to_ca(self, tmp_path):
        table = tmp_path / "patches.tsv"
        table.write_text("rank\thotspot_residues\n1\t13\n")
        sites = read_patch_sites(str(table), self._structure())
        assert sites[0].centroid == pytest.approx([13.0, 0.0, 0.0])

    def test_a_residue_not_in_the_structure_is_an_error(self, tmp_path):
        """The table and the structure must come from the same prep."""
        table = tmp_path / "patches.tsv"
        table.write_text("rank\thotspot_residues\n1\t10+999\n")
        with pytest.raises(ValueError, match="not in the structure"):
            read_patch_sites(str(table), self._structure())

    def test_the_wrong_table_is_refused_by_its_columns(self, tmp_path):
        table = tmp_path / "residues.tsv"
        table.write_text("chain\tresi\tresn\n A\t10\tLEU\n")
        with pytest.raises(ValueError, match="not a hotspot patch table"):
            read_patch_sites(str(table), self._structure())

    def test_max_patches_truncates_by_rank(self, tmp_path):
        table = tmp_path / "patches.tsv"
        table.write_text("rank\thotspot_residues\n2\t11\n1\t10\n3\t12\n")
        sites = read_patch_sites(str(table), self._structure(), max_patches=2)
        assert [s.rank for s in sites] == [1, 2]


class TestTables:
    def _result(self):
        up = side("UP", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]), ("10", "11"))],
                  seed=1)
        dn = side("DN", [0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [7.0, 0.0, 0.0],
                  patches=[PatchSite(1, np.array([0.0, 8.0, 0.0]), ("20",))], seed=2)
        params = PoseParams(dic_min=6, dic_max=30, dic_step=8, oic=(180.0,),
                            azimuth=(0.0,), tilt=(0.0, BACKWARDS), spin=(0.0,),
                            tilt_azimuth=(0.0,), bridge_lo=5.0, bridge_hi=60.0,
                            max_ok=0)
        return search(up, dn, params)

    def test_poses_table_has_a_row_per_pose(self, tmp_path):
        result = self._result()
        path = write_poses(tmp_path / "poses.tsv", result)
        lines = open(path).read().strip().split("\n")
        assert len(lines) == 1 + len(result.poses)
        header = lines[0].split("\t")
        assert header[0] == "geom_id" and "O_IC_rev" in header
        assert "tilt_azimuth" in header

    def test_unmeasured_distance_is_a_dash_not_a_zero(self, tmp_path):
        result = self._result()
        path = write_poses(tmp_path / "poses.tsv", result)
        rows = [line.split("\t") for line in open(path).read().strip().split("\n")[1:]]
        header = open(path).readline().rstrip("\n").split("\t")
        col = header.index("min_interprotein_dist")
        clashes = [r for r in rows if r[header.index("status")] == "clash"]
        assert clashes, "expected the tight end of the grid to clash"
        assert all(r[col] == "-" for r in clashes)

    def test_patch_pairs_and_best_agree(self, tmp_path):
        result = self._result()
        write_patch_pairs(tmp_path / "patch_pairs.tsv", result)
        best = write_patch_pair_best(tmp_path / "patch_pair_best.tsv", result)
        head = open(best).readline().rstrip("\n").split("\t")
        assert head[:2] == ["UP_patch", "DN_patch"]
        rows = open(best).read().strip().split("\n")[1:]
        if result.feasible:
            assert rows, "a feasible pose must put its clear pair in the table"
            for row in rows:
                assert row.split("\t")[2].startswith("geom")


class TestWriteComplex:
    def test_chain_collision_is_refused(self, tmp_path):
        from channelbridge.io.pdb import Atom, AtomArray
        from channelbridge.pose import GridPoint, write_complex

        def one(chain):
            return AtomArray([Atom(serial=1, name=" CA ", altloc="", resname="ALA",
                                   chain=chain, resseq=1, icode=" ",
                                   xyz=(0.0, 0.0, 0.0), occupancy=1.0, bfactor=0.0,
                                   element="C", hetatm=False)])

        up, dn = simple_pair()
        up.atoms, dn.atoms = one("A"), one("A")
        rot, trans = pose_transform(up, dn, GridPoint(20, 180, 0, 0, 0, 0))
        from channelbridge.pose import Pose

        pose = Pose(pose=1, status="ok", d_ic=20.0, o_ic=180.0, o_ic_rev=180.0,
                    do=0.9, do_clamped=0.9, phi=0.0, theta_up=0.0, theta_dn=0.0,
                    requested=GridPoint(20, 180, 0, 0, 0, 0), rotation=rot,
                    translation=trans)
        with pytest.raises(ValueError, match="merge into one chain"):
            write_complex(tmp_path / "x.pdb", up, dn, pose)

    def test_complex_holds_both_enzymes_moved(self, tmp_path):
        from channelbridge.io.pdb import Atom, AtomArray, read_pdb
        from channelbridge.pose import GridPoint, Pose, write_complex

        def one(chain, xyz):
            return AtomArray([Atom(serial=1, name=" CA ", altloc="", resname="ALA",
                                   chain=chain, resseq=1, icode=" ", xyz=xyz,
                                   occupancy=1.0, bfactor=0.0, element="C",
                                   hetatm=False)])

        up, dn = simple_pair()
        up.atoms = one("A", (0.0, 0.0, 0.0))
        dn.atoms = one("B", (8.0, 0.0, 0.0))
        gp = GridPoint(20.0, 180.0, 0.0, 0.0, 0.0, 0.0)
        rot, trans = pose_transform(up, dn, gp)
        pose = Pose(pose=1, status="ok", d_ic=20.0, o_ic=180.0, o_ic_rev=180.0,
                    do=0.9, do_clamped=0.9, phi=0.0, theta_up=0.0, theta_dn=0.0,
                    requested=gp, rotation=rot, translation=trans)
        path = write_complex(tmp_path / "complex.pdb", up, dn, pose)
        back = read_pdb(path)
        assert {a.chain for a in back} == {"A", "B"}
        moved = next(a for a in back if a.chain == "B")
        # the downstream exit was at (8,0,0) and goes to the upstream exit
        # plus 20 A along +x, which is where its one atom now sits
        assert np.asarray(moved.xyz) == pytest.approx([28.0, 0.0, 0.0], abs=1e-6)


class TestParamsFromSpec:
    def test_pair_file_pose_section_drives_the_search(self):
        from channelbridge.config import PoseSpec

        spec = PoseSpec.from_dict("pose", {"dic_min": 8, "dic_max": 10,
                                           "oic": [180], "max_ok": 4})
        params = PoseParams.from_spec(spec, tilt=(0.0,), spin=(0.0,))
        assert params.dic_min == 8 and params.max_ok == 4
        assert params.oic == (180.0,)
        assert params.tilt == (0.0,)
        assert len(pose_grid(params)) == 3

    def test_do_coefficients_are_carried_through(self):
        up, dn = simple_pair()
        from channelbridge.pose import GridPoint

        rot, trans = pose_transform(up, dn, GridPoint(20, 180, 0, 0, 0, 0))
        row = score_pose(up, dn, rot, trans, DO_PUBLISHED)
        assert row["DO"] == pytest.approx(do_score(20.0, 180.0, DO_PUBLISHED))


def test_min_dic_is_documented_not_magic():
    """The floor exists because the O_IC vertex degenerates, not by taste."""
    assert MIN_DIC == 0.5
    with pytest.raises(ValueError):
        o_ic([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
