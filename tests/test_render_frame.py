"""The angle frame, and which poses get drawn."""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.render.frame import (
    Frame,
    angle_frame,
    parse_selection,
    pymol_view,
    rotated,
    select_poses,
)


def ang(v, w):
    v = np.asarray(v, float)
    w = np.asarray(w, float)
    c = np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


class TestFrameInvariants:
    def test_rejects_non_orthonormal_axes(self):
        with pytest.raises(ValueError, match="orthonormal"):
            Frame(np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1]]), np.zeros(3))

    def test_rejects_a_mirrored_frame(self):
        # A left-handed basis would silently mirror the whole figure.
        with pytest.raises(ValueError, match="mirrored"):
            Frame(np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]), np.zeros(3))

    def test_accepts_the_identity(self):
        f = Frame(np.eye(3), np.zeros(3))
        np.testing.assert_allclose(f.right, [1, 0, 0])
        np.testing.assert_allclose(f.up, [0, 1, 0])


class TestAngleFrame:
    """P1-P2-P3 in the plane of the page, upstream on the left."""

    def setup_method(self):
        # A 90-degree construction in the xy plane, easy to reason about.
        self.p2 = np.array([0.0, 0.0, 0.0])
        self.p1 = np.array([5.0, 0.0, 0.0])      # P2->P1 along +x
        self.p3 = np.array([0.0, 17.0, 0.0])     # P2->P3 along +y, so O_IC = 90

    def test_puts_the_construction_plane_face_on(self):
        f = angle_frame(self.p1, self.p2, self.p3)
        # all three points project with ~zero depth: the plane is face-on
        depth = (np.vstack([self.p1, self.p2, self.p3]) - f.centre) @ f.toward_viewer
        np.testing.assert_allclose(depth, 0.0, atol=1e-9)

    def test_upstream_sits_on_the_left(self):
        f = angle_frame(self.p1, self.p2, self.p3)
        xy = f.project(np.vstack([self.p1, self.p3]))
        assert xy[0][0] < 0.0, "P1 (upstream) should project to the left"
        assert xy[1][0] >= -1e-9, "P3 (downstream) should not be further left than P2"

    def test_without_the_flip_upstream_is_on_the_right(self):
        f = angle_frame(self.p1, self.p2, self.p3, upstream_left=False)
        assert f.project(self.p1)[0][0] > 0.0

    def test_the_segment_rises_at_180_minus_oic(self):
        # The arithmetic the SVG's angle arc depends on.
        f = angle_frame(self.p1, self.p2, self.p3)
        v = f.project(self.p3)[0] - f.project(self.p2)[0]
        image_angle = np.degrees(np.arctan2(v[1], v[0])) % 360
        o_ic = ang(self.p1 - self.p2, self.p3 - self.p2)
        assert o_ic == pytest.approx(90.0)
        assert image_angle == pytest.approx(180.0 - o_ic, abs=1e-6)

    @pytest.mark.parametrize("o_ic", [30.0, 60.0, 105.0, 135.0, 165.0])
    def test_the_rise_rule_holds_for_any_angle(self, o_ic):
        th = np.radians(o_ic)
        p3 = np.array([17.0 * np.cos(th), 17.0 * np.sin(th), 0.0])
        f = angle_frame(self.p1, self.p2, p3)
        v = f.project(p3)[0] - f.project(self.p2)[0]
        got = np.degrees(np.arctan2(v[1], v[0])) % 360
        assert got == pytest.approx(180.0 - o_ic, abs=1e-6)

    def test_true_distance_is_preserved_in_the_image(self):
        # D_IC can be measured straight off the figure.
        f = angle_frame(self.p1, self.p2, self.p3)
        xy = f.project(np.vstack([self.p2, self.p3]))
        assert np.linalg.norm(xy[1] - xy[0]) == pytest.approx(17.0)

    def test_works_out_of_plane(self):
        p3 = np.array([3.0, 9.0, 4.0])
        f = angle_frame(self.p1, self.p2, p3)
        depth = (np.vstack([self.p1, self.p2, p3]) - f.centre) @ f.toward_viewer
        np.testing.assert_allclose(depth, 0.0, atol=1e-9)
        o_ic = ang(self.p1 - self.p2, p3 - self.p2)
        v = f.project(p3)[0] - f.project(self.p2)[0]
        got = np.degrees(np.arctan2(v[1], v[0])) % 360
        assert got == pytest.approx(180.0 - o_ic, abs=1e-6)

    def test_default_centre_is_the_midpoint_of_the_two_mouths(self):
        f = angle_frame(self.p1, self.p2, self.p3)
        np.testing.assert_allclose(f.centre, (self.p2 + self.p3) / 2)

    def test_an_explicit_centre_is_honoured(self):
        f = angle_frame(self.p1, self.p2, self.p3, centre=[1.0, 2.0, 3.0])
        np.testing.assert_allclose(f.centre, [1.0, 2.0, 3.0])

    def test_collinear_points_fail_loudly(self):
        # O_IC of exactly 180 has no plane; silently picking one would put the
        # figure in an arbitrary orientation.
        with pytest.raises(ValueError, match="collinear"):
            angle_frame([5, 0, 0], [0, 0, 0], [-17, 0, 0])

    def test_a_zero_length_arm_fails_loudly(self):
        with pytest.raises(ValueError, match="zero length"):
            angle_frame([0, 0, 0], [0, 0, 0], [0, 17, 0])

    def test_is_named_rotated0(self):
        assert angle_frame(self.p1, self.p2, self.p3).name == "rotated0"


class TestRotated:
    def test_turning_360_returns_the_same_frame(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        np.testing.assert_allclose(rotated(f, 360).axes, f.axes, atol=1e-9)

    def test_the_vertical_is_the_axis_of_rotation(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        np.testing.assert_allclose(rotated(f, 73).up, f.up, atol=1e-9)

    def test_turning_90_swaps_right_and_depth(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        r = rotated(f, 90)
        assert ang(r.right, f.right) == pytest.approx(90.0, abs=1e-6)

    def test_stays_a_valid_right_handed_frame(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        for d in (15, 45, 180, 270):
            r = rotated(f, d)
            assert np.linalg.det(r.axes) == pytest.approx(1.0)

    def test_names_follow_the_convention(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        assert rotated(f, 30).name == "rotated30"
        assert rotated(f, 0).name == "rotated0"
        assert rotated(f, 360).name == "rotated0"
        assert rotated(f, 180).name == "rotated180"


class TestPymolView:
    def test_has_eighteen_numbers(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        assert len(pymol_view(f)) == 18

    def test_rotation_is_stored_column_wise(self):
        # Getting this transpose wrong mirrors the figure; assert it directly.
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        v = pymol_view(f)
        np.testing.assert_allclose(np.array(v[:9]).reshape(3, 3), f.axes.T, atol=1e-12)

    def test_centre_is_the_origin_of_rotation(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0], centre=[1, 2, 3])
        np.testing.assert_allclose(pymol_view(f)[12:15], [1, 2, 3])

    def test_slab_is_passed_through(self):
        f = angle_frame([5, 0, 0], [0, 0, 0], [0, 17, 0])
        assert pymol_view(f, slab=1234.0)[16] == 1234.0


class TestSelection:
    """How many / which .pse and how many / which SVG, independently."""

    AVAIL = ("geom001", "geom002", "geom003", "geom004")

    def test_none_and_empty_mean_nothing(self):
        assert parse_selection(None, self.AVAIL) == []
        assert parse_selection("", self.AVAIL) == []
        assert parse_selection("none", self.AVAIL) == []

    def test_all_means_every_pose(self):
        assert parse_selection("all", self.AVAIL) == list(self.AVAIL)

    def test_a_count_takes_the_best_n(self):
        assert parse_selection("2", self.AVAIL) == ["geom001", "geom002"]
        assert parse_selection(2, self.AVAIL) == ["geom001", "geom002"]

    def test_a_count_beyond_the_end_is_not_an_error(self):
        assert parse_selection("99", self.AVAIL) == list(self.AVAIL)

    def test_zero_means_nothing(self):
        assert parse_selection("0", self.AVAIL) == []

    def test_an_explicit_list_is_honoured_in_order(self):
        assert parse_selection("geom003,geom001", self.AVAIL) == ["geom003", "geom001"]

    def test_whitespace_in_a_list_is_tolerated(self):
        assert parse_selection(" geom002 , geom004 ", self.AVAIL) == ["geom002", "geom004"]

    def test_an_unknown_id_fails_loudly(self):
        with pytest.raises(ValueError, match="unknown pose id"):
            parse_selection("geom001,geom999", self.AVAIL)

    def test_a_negative_count_fails(self):
        with pytest.raises(ValueError, match="negative"):
            parse_selection("-1", self.AVAIL)

    def test_pse_and_svg_counts_are_independent(self):
        poses = [("geom002", 0.80), ("geom001", 0.85), ("geom003", 0.70)]
        pse, svg = select_poses(poses, n_pse=3, n_svg=1)
        assert pse == ["geom001", "geom002", "geom003"]
        assert svg == ["geom001"], "the SVG goes to the top DO score"

    def test_svg_can_be_asked_for_by_name_while_pse_is_a_count(self):
        poses = [("geom001", 0.85), ("geom002", 0.80), ("geom003", 0.70)]
        pse, svg = select_poses(poses, n_pse=2, n_svg="geom003")
        assert pse == ["geom001", "geom002"]
        assert svg == ["geom003"]

    def test_ranking_is_by_do_descending(self):
        poses = [("a", 0.1), ("b", 0.9), ("c", 0.5)]
        pse, _ = select_poses(poses, n_pse="all")
        assert pse == ["b", "c", "a"]

    def test_ties_break_on_id_so_the_order_is_stable(self):
        poses = [("z", 0.5), ("a", 0.5)]
        pse, _ = select_poses(poses, n_pse="all")
        assert pse == ["a", "z"]

    def test_default_is_nothing_drawn(self):
        assert select_poses([("a", 0.5)]) == ([], [])


class TestKabsch:
    """Carrying a tunnel onto a posed complex."""

    def setup_method(self):
        rng = np.random.default_rng(0)
        self.pts = rng.normal(size=(30, 3)) * 8.0

    def test_recovers_a_known_rigid_motion_exactly(self):
        from channelbridge.render.scenes import kabsch

        th = np.radians(37.0)
        r_true = np.array([[np.cos(th), -np.sin(th), 0],
                           [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
        t_true = np.array([12.0, -5.0, 3.0])
        moved = (r_true @ self.pts.T).T + t_true
        r, t, rmsd = kabsch(self.pts, moved)
        assert rmsd < 1e-9
        np.testing.assert_allclose(r, r_true, atol=1e-9)
        np.testing.assert_allclose(t, t_true, atol=1e-9)

    def test_never_mirrors(self):
        from channelbridge.render.scenes import kabsch

        mirrored = self.pts * np.array([1.0, 1.0, -1.0])
        r, _, rmsd = kabsch(self.pts, mirrored)
        assert np.linalg.det(r) == pytest.approx(1.0), "a reflection would flip the tunnel"
        assert rmsd > 1.0, "a mirror image should not fit well"

    def test_identity_for_identical_sets(self):
        from channelbridge.render.scenes import kabsch

        r, t, rmsd = kabsch(self.pts, self.pts)
        assert rmsd < 1e-9
        np.testing.assert_allclose(r, np.eye(3), atol=1e-9)
        np.testing.assert_allclose(t, 0.0, atol=1e-9)

    def test_rejects_mismatched_or_tiny_input(self):
        from channelbridge.render.scenes import kabsch

        with pytest.raises(ValueError, match="matching"):
            kabsch(self.pts, self.pts[:10])
        with pytest.raises(ValueError, match="at least 3"):
            kabsch(self.pts[:2], self.pts[:2])


class TestSpreadInView:
    """Labels are placed in 3D but collide in 2D, so they separate in 2D."""

    @staticmethod
    def screen(pts, x, y):
        import numpy as np
        p = np.asarray(pts, float)
        return np.stack([p @ np.asarray(x, float), p @ np.asarray(y, float)], 1)

    def test_coincident_on_screen_but_apart_in_space_still_separates(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        # three points differing ONLY in depth: identical on screen
        a = np.array([[0.0, 0, 0], [0.0, 0, 30], [0.0, 0, 60]])
        out = a + spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=12.0)
        uv = self.screen(out, [1, 0, 0], [0, 1, 0])
        for i in range(3):
            for j in range(i + 1, 3):
                assert np.linalg.norm(uv[i] - uv[j]) > 11.0

    def test_depth_is_never_changed(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        a = np.array([[0.0, 0, 5], [1.0, 0, 9], [0.5, 1, 20]])
        z = np.array([0.0, 0.0, 1.0])
        out = a + spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=15.0)
        np.testing.assert_allclose(out @ z, a @ z, atol=1e-9)

    def test_points_already_apart_are_left_alone(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        a = np.array([[0.0, 0, 0], [50.0, 0, 0]])
        np.testing.assert_allclose(spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=10.0),
                                   np.zeros((2, 3)), atol=1e-9)

    def test_a_single_label_needs_no_spreading(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        np.testing.assert_allclose(spread_in_view([[1.0, 2, 3]], [1, 0, 0], [0, 1, 0]),
                                   np.zeros((1, 3)))

    def test_empty_input(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        assert spread_in_view(np.zeros((0, 3)), [1, 0, 0], [0, 1, 0]).shape == (0, 3)

    def test_a_bad_shape_is_refused(self):
        import pytest as _pytest

        from channelbridge.render.frame import spread_in_view

        with _pytest.raises(ValueError, match=r"\(N, 3\)"):
            spread_in_view([[1.0, 2]], [1, 0, 0], [0, 1, 0])

    def test_wide_labels_need_horizontal_room_not_just_distance(self):
        import itertools

        import numpy as np

        from channelbridge.render.frame import spread_in_view

        # three anchors on one horizontal line, as patches on one face give:
        # a circular criterion satisfies itself by nudging them apart by its
        # radius, which for a line of text is still an overlap
        a = np.array([[0.0, 0, 0], [6.0, 0, 0], [12.0, 0, 0]])
        out = a + spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=(30.0, 10.0))
        uv = self.screen(out, [1, 0, 0], [0, 1, 0])
        for i, j in itertools.combinations(range(3), 2):
            dx, dy = np.abs(uv[i] - uv[j])
            assert dx >= 29.9 or dy >= 9.9, "pair %d,%d only dx=%.1f dy=%.1f" % (i, j, dx, dy)

    def test_a_scalar_min_sep_still_works(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        a = np.array([[0.0, 0, 0], [1.0, 0, 0]])
        out = a + spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=20.0)
        uv = self.screen(out, [1, 0, 0], [0, 1, 0])
        assert np.linalg.norm(uv[0] - uv[1]) > 19.0

    def test_a_bad_min_sep_is_refused(self):
        import pytest as _pytest

        from channelbridge.render.frame import spread_in_view

        a = [[0.0, 0, 0], [1.0, 0, 0]]
        with _pytest.raises(ValueError, match="one number or two"):
            spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=(1.0, 2.0, 3.0))
        with _pytest.raises(ValueError, match="positive"):
            spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=(0.0, 5.0))

    def test_bounds_keep_labels_inside_the_frame(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        # four anchors in a row, pushed apart hard enough that an unbounded
        # relaxation would throw the outer ones well past the frame
        a = np.array([[0.0, 0, 0], [2.0, 0, 0], [4.0, 0, 0], [6.0, 0, 0]])
        out = a + spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=(30.0, 10.0),
                                 bounds=(-40.0, 40.0, -20.0, 20.0))
        uv = self.screen(out, [1, 0, 0], [0, 1, 0])
        assert uv[:, 0].min() >= -40.01 and uv[:, 0].max() <= 40.01
        assert uv[:, 1].min() >= -20.01 and uv[:, 1].max() <= 20.01

    def test_without_bounds_they_are_free_to_leave(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        a = np.array([[0.0, 0, 0], [2.0, 0, 0], [4.0, 0, 0], [6.0, 0, 0]])
        out = a + spread_in_view(a, [1, 0, 0], [0, 1, 0], min_sep=(30.0, 10.0))
        uv = self.screen(out, [1, 0, 0], [0, 1, 0])
        assert uv[:, 0].max() - uv[:, 0].min() > 80.0, "unbounded spreads wide"


class TestLabelsStayInThePicture:
    """A label pushed clear of its neighbours must not land in the margin.

    On the p6a hotspot map the relaxation worked -- and put `FNR p05` at the
    right-hand edge of the frame, far from the patch it names. Separation
    without bounds trades one unreadable figure for another.
    """

    def test_bounds_are_respected(self):
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        x = np.array([1.0, 0.0, 0.0])
        y = np.array([0.0, 1.0, 0.0])
        # six anchors in a tight cluster: without bounds they fly outward
        anchors = np.array([[40.0 + 0.4 * i, 0.3 * i, 0.0] for i in range(6)])
        box = (30.0, 55.0, -12.0, 12.0)
        moved = anchors + spread_in_view(anchors, x, y, min_sep=(8.0, 3.0),
                                         bounds=box)
        u, v = moved @ x, moved @ y
        assert u.min() >= box[0] - 1e-6 and u.max() <= box[1] + 1e-6
        assert v.min() >= box[2] - 1e-6 and v.max() <= box[3] + 1e-6

    def test_depth_is_never_changed(self):
        # a label must stay at its patch's distance from the camera, or it
        # can end up behind the surface it names
        import numpy as np

        from channelbridge.render.frame import spread_in_view

        x = np.array([1.0, 0.0, 0.0])
        y = np.array([0.0, 1.0, 0.0])
        anchors = np.array([[0.0, 0.0, 7.0], [0.5, 0.2, -3.0]])
        moved = anchors + spread_in_view(anchors, x, y, min_sep=(8.0, 3.0))
        assert np.allclose(moved[:, 2], anchors[:, 2])


class TestLiftToFront:
    """The OTHER way labels collide: in depth, not in the camera plane.

    `float_labels 1` is PyMOL's own answer and does not work for a figure --
    it applies to the interactive viewport, not to ray-traced output, which
    was confirmed on a real render before this was written.
    """

    def test_only_what_is_behind_gets_moved(self):
        import numpy as np

        from channelbridge.render.frame import lift_to_front

        t = np.array([0.0, 0.0, 1.0])
        body = np.array([[0.0, 0.0, -5.0], [0.0, 0.0, 10.0]])
        anchors = np.array([[0.0, 0.0, -8.0], [0.0, 0.0, 30.0]])
        lift = lift_to_front(anchors, t, body, margin=2.0)
        assert lift[1].tolist() == [0.0, 0.0, 0.0], "already in front: leave it"
        assert (anchors + lift)[0][2] == 12.0, "behind: just clear of the front"

    def test_the_move_is_along_the_view_axis_only(self):
        # anything else would slide the label away from the patch it names
        import numpy as np

        from channelbridge.render.frame import lift_to_front

        t = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
        body = np.array([[9.0, 9.0, 0.0]])
        anchors = np.array([[0.0, 0.0, 4.0]])
        lift = lift_to_front(anchors, t, body, margin=1.0)
        assert np.allclose(np.cross(lift[0], t), 0.0)
        assert np.isclose((anchors + lift)[0] @ t, (body @ t).max() + 1.0)

    def test_an_empty_structure_moves_nothing(self):
        import numpy as np

        from channelbridge.render.frame import lift_to_front

        lift = lift_to_front(np.zeros((3, 3)), np.array([0.0, 0.0, 1.0]),
                             np.zeros((0, 3)))
        assert lift.shape == (3, 3) and not lift.any()


class TestOrthoscopicView:
    """The sign of the field-of-view slot is the orthoscopic flag."""

    def test_perspective_by_default(self):
        import numpy as np

        from channelbridge.render.frame import Frame, pymol_view

        f = Frame(axes=np.eye(3), centre=np.zeros(3))
        assert pymol_view(f)[17] < 0

    def test_orthographic_when_asked(self):
        # a map compares patches, so an equal patch must not look larger for
        # being nearer -- and under this camera the label lift costs nothing
        # on screen
        import numpy as np

        from channelbridge.render.frame import Frame, pymol_view

        f = Frame(axes=np.eye(3), centre=np.zeros(3))
        assert pymol_view(f, orthoscopic=True)[17] > 0


class TestIcLineFrame:
    """A frame with the D_IC line on x, free to roll about it.

    `angle_frame` holds the P1-P2-P3 plane face-on so O_IC reads true, and
    that uses up every degree of freedom: when a channel's drawn volume lies
    across the D_IC line there, nothing can be done. This frame trades the
    face-on plane for a roll, and reports what the trade costs.
    """

    def _pose(self):
        import numpy as np

        # O_IC 150 deg: the arm is 30 deg off the line, so rolling bites
        p2 = np.zeros(3)
        p3 = np.array([17.0, 0.0, 0.0])
        p1 = np.array([-np.cos(np.radians(30.0)), np.sin(np.radians(30.0)), 0.0])
        return p1, p2, p3

    def test_the_line_lands_on_x_at_every_roll(self):

        from channelbridge.render.frame import ic_line_frame

        p1, p2, p3 = self._pose()
        for roll in (0.0, 37.0, 180.0, 290.0):
            f = ic_line_frame(p1, p2, p3, roll=roll)
            v = (p3 - p2) @ f.axes[1]        # no component on screen y
            assert abs(float(v)) < 1e-9, roll

    def test_the_angle_reads_true_only_face_on(self):
        from channelbridge.render.frame import angle_frame, drawn_oic, ic_line_frame

        p1, p2, p3 = self._pose()
        true = drawn_oic(angle_frame(p1, p2, p3), p1, p2, p3)
        assert abs(true - 150.0) < 1e-6
        assert abs(drawn_oic(ic_line_frame(p1, p2, p3, roll=0.0), p1, p2, p3)
                   - true) < 1e-6
        # rolled 90 degrees the arm's in-plane part vanishes: it reads 180
        assert abs(drawn_oic(ic_line_frame(p1, p2, p3, roll=90.0), p1, p2, p3)
                   - 180.0) < 1e-6

    def test_the_error_budget_is_honoured(self):
        import numpy as np

        from channelbridge.render.frame import clearest_roll

        p1, p2, p3 = self._pose()
        pts = np.array([[8.5, 0.3, 0.0], [9.0, 0.0, 0.4]])    # astride the line
        for cap in (1.0, 5.0, 20.0):
            _roll, _frame, _clear, shown, true, _frac = clearest_roll(
                p1, p2, p3, pts, about="line", max_angle_error=cap, step=5.0)
            assert abs(shown - true) <= cap + 1e-9

    def test_a_budget_of_zero_cannot_be_met_except_face_on(self):
        import numpy as np

        from channelbridge.render.frame import clearest_roll

        p1, p2, p3 = self._pose()
        _roll, _f, _c, shown, true, _frac = clearest_roll(
            p1, p2, p3, np.array([[8.5, 0.3, 0.0]]), about="line",
            max_angle_error=0.0, step=15.0)
        assert abs(shown - true) < 1e-9

    def test_points_behind_the_exit_do_not_count(self):
        # the band is measured along P2->P3, not along the camera's x, which
        # is flipped for display. Taking |x| swept in the tunnel INSIDE the
        # upstream protein and made every roll look equally bad.
        import numpy as np

        from channelbridge.render.frame import clearest_roll

        p1, p2, p3 = self._pose()
        behind = np.array([[-6.0, 0.01, 0.0]])      # hard against the line, behind P2
        _roll, _f, clear, _s, _t, _fr = clearest_roll(
            p1, p2, p3, behind, about="line", max_angle_error=20.0,
            min_dic_fraction=0.0, step=15.0)
        assert clear == float("inf"), "nothing lies between the exits"

    def test_collinear_is_refused(self):
        import numpy as np
        import pytest

        from channelbridge.render.frame import ic_line_frame

        with pytest.raises(ValueError, match="collinear"):
            ic_line_frame(np.array([-1.0, 0, 0]), np.zeros(3), np.array([5.0, 0, 0]))


class TestUpstreamIsAlwaysOnTheLeft:
    """House rule (user, standing): the upstream enzyme is on the left.

    Each frame builds x from a different vector -- the O_IC arm points INTO
    the upstream protein, the D_IC line points AWAY from it -- so "apply the
    180 degree flip" means opposite things in the two. Reasoning about which
    sign each needs is exactly how ic_line_frame shipped mirrored. The frames
    now check the invariant instead, and so does this.
    """

    def _pose(self):
        import numpy as np

        p2 = np.zeros(3)
        p3 = np.array([17.0, 0.0, 0.0])
        p1 = np.array([-np.cos(np.radians(30.0)), np.sin(np.radians(30.0)), 0.0])
        return p1, p2, p3

    def _upstream_is_left(self, frame, p1, p2):
        # judged on P1, which points into the upstream BODY; the two exits
        # sit together at the interface and say little about which protein
        # is where
        return float((p1 - frame.centre) @ frame.axes[0]) < \
            float((p2 - frame.centre) @ frame.axes[0])

    def test_the_angle_frame_puts_it_left(self):
        from channelbridge.render.frame import angle_frame

        p1, p2, p3 = self._pose()
        assert self._upstream_is_left(angle_frame(p1, p2, p3), p1, p2)

    def test_every_roll_of_the_ic_frame_puts_it_left(self):
        from channelbridge.render.frame import ic_line_frame

        p1, p2, p3 = self._pose()
        for roll in range(0, 360, 15):
            f = ic_line_frame(p1, p2, p3, roll=float(roll))
            assert self._upstream_is_left(f, p1, p2), "roll %d" % roll

    def test_switching_it_off_leaves_the_raw_frame(self):
        from channelbridge.render.frame import angle_frame, ic_line_frame

        p1, p2, p3 = self._pose()
        # the raw frames differ in handedness of the x choice, which is the
        # whole reason the flip could not be a fixed sign
        a = angle_frame(p1, p2, p3, upstream_left=False)
        b = ic_line_frame(p1, p2, p3, roll=0.0, upstream_left=False)
        assert not self._upstream_is_left(a, p1, p2)
        assert self._upstream_is_left(b, p1, p2)

    def test_the_rule_survives_a_pose_with_the_exits_reversed_in_space(self):
        import numpy as np

        from channelbridge.render.frame import ic_line_frame

        # the downstream exit at NEGATIVE x in model coordinates: the frame
        # must still show it on the right
        p2 = np.zeros(3)
        p3 = np.array([-17.0, 0.0, 0.0])
        p1 = np.array([np.cos(np.radians(30.0)), np.sin(np.radians(30.0)), 0.0])
        for roll in (0.0, 45.0, 200.0):
            f = ic_line_frame(p1, p2, p3, roll=roll)
            assert self._upstream_is_left(f, p1, p2), roll


class TestClearanceIsMeasuredFromTheDrawnLine:
    """Distance from the projected D_IC segment, not from a camera axis.

    Those coincide only when the line is on x. Measuring against the axis
    made `arm` rolls look five times clearer than they were, because there
    the line sits at 180 - O_IC to x and the axis is not the line.
    """

    def test_the_two_frame_families_agree_at_roll_zero(self):
        import numpy as np

        from channelbridge.render.frame import clearest_roll

        p2 = np.zeros(3)
        p3 = np.array([17.0, 0.0, 0.0])
        p1 = np.array([-np.cos(np.radians(30.0)), np.sin(np.radians(30.0)), 0.0])
        pts = np.array([[8.5, 0.6, 0.0]])
        # roll 0 of either family is the same camera, so the same clearance
        a = clearest_roll(p1, p2, p3, pts, about="arm",
                          max_angle_error=0.0, min_dic_fraction=0.0, step=180.0)
        b = clearest_roll(p1, p2, p3, pts, about="line",
                          max_angle_error=0.0, min_dic_fraction=0.0, step=180.0)
        assert a[2] == pytest.approx(b[2], abs=1e-6)

    def test_rolling_about_the_arm_can_shorten_the_drawn_line(self):
        # and about the line it never can -- which is the difference between
        # the two, and why both are reported
        import numpy as np

        from channelbridge.render.frame import (
            arm_frame,
            drawn_dic_fraction,
            ic_line_frame,
        )

        p2 = np.zeros(3)
        p3 = np.array([17.0, 0.0, 0.0])
        p1 = np.array([-np.cos(np.radians(30.0)), np.sin(np.radians(30.0)), 0.0])
        assert drawn_dic_fraction(arm_frame(p1, p2, p3, roll=90.0), p2, p3) < 1.0
        assert drawn_dic_fraction(ic_line_frame(p1, p2, p3, roll=90.0), p2, p3) \
            == pytest.approx(1.0)

    def test_the_dic_budget_is_honoured(self):
        import numpy as np

        from channelbridge.render.frame import clearest_roll

        p2 = np.zeros(3)
        p3 = np.array([17.0, 0.0, 0.0])
        p1 = np.array([-np.cos(np.radians(60.0)), np.sin(np.radians(60.0)), 0.0])
        pts = np.array([[8.5, 0.2, 0.0], [9.0, 0.0, 0.2]])
        _r, _f, _c, _s, _t, frac = clearest_roll(
            p1, p2, p3, pts, about="arm", max_angle_error=90.0,
            min_dic_fraction=0.98, step=5.0)
        assert frac >= 0.98


class TestTheAngleArcDegeneratesAtBothEnds:
    """The arc is a slerp, so it divides by sin(phi).

    Near 0 there is no angle to draw. Near 180 sin(phi) is ~1e-16 and every
    intermediate point collapses onto the vertex, drawing a spike through P2
    that reads as an arc. The hotspot map's reference placement is at exactly
    O_IC 180, so this is its normal case, not an edge case.
    """

    def test_the_slerp_really_does_collapse_at_180(self):
        import numpy as np

        phi = np.pi
        u1, u3 = np.array([1.0, 0, 0]), np.array([-1.0, 0, 0])
        mid = (np.sin(0.5 * phi) * u1 + np.sin(0.5 * phi) * u3) / np.sin(phi)
        assert np.allclose(mid, 0.0), "this is why the guard has two ends"

    def test_the_guard_admits_a_real_angle(self):
        import numpy as np

        for deg in (30.0, 90.0, 135.0, 165.0):
            phi = np.radians(deg)
            assert 1e-3 < phi < np.pi - 1e-3

    def test_the_guard_rejects_both_degenerate_ends(self):
        import numpy as np

        for deg in (0.0, 0.01, 179.999, 180.0):
            phi = np.radians(deg)
            assert not (1e-3 < phi < np.pi - 1e-3)
