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
