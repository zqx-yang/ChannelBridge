"""Cross-section geometry and the binder envelope."""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.render.cut import (
    capsule_points,
    cut_mask,
    enzyme_section_point,
    parse_hotspot_residues,
    patch_centroid,
    section_plane,
)


class TestSectionPlane:
    def test_passes_through_the_midpoint_of_the_two_exits(self):
        pt, _n = section_plane([0, 0, 0], [10, 0, 0], [0, 0, 1])
        np.testing.assert_allclose(pt, [5, 0, 0])

    def test_faces_the_camera(self):
        _, n = section_plane([0, 0, 0], [10, 0, 0], [0, 0, 3])
        np.testing.assert_allclose(n, [0, 0, 1])

    def test_normal_is_unit_length(self):
        _, n = section_plane([0, 0, 0], [1, 2, 3], [4, 4, 4])
        assert np.linalg.norm(n) == pytest.approx(1.0)

    def test_a_zero_view_direction_fails_loudly(self):
        with pytest.raises(ValueError, match="zero length"):
            section_plane([0, 0, 0], [1, 0, 0], [0, 0, 0])


class TestCutMask:
    def setup_method(self):
        # a slab of atoms along z, the camera looking down +z
        self.xyz = np.array([[0, 0, z] for z in range(-10, 11)], float)
        self.point = np.array([0.0, 0.0, 0.0])
        self.normal = np.array([0.0, 0.0, 1.0])

    def test_removes_the_camera_side(self):
        m = cut_mask(self.xyz, self.point, self.normal)
        kept = self.xyz[~m]
        assert kept[:, 2].max() <= 0, "everything in front of the plane should go"
        assert m.sum() > 0 and (~m).sum() > 0, "the cut should actually split"

    def test_positive_offset_cuts_deeper(self):
        shallow = cut_mask(self.xyz, self.point, self.normal, offset=0.0).sum()
        deep = cut_mask(self.xyz, self.point, self.normal, offset=5.0).sum()
        assert deep > shallow, "a positive offset must remove more, not less"

    def test_negative_offset_leaves_more_in_front(self):
        base = cut_mask(self.xyz, self.point, self.normal, offset=0.0).sum()
        less = cut_mask(self.xyz, self.point, self.normal, offset=-5.0).sum()
        assert less < base

    def test_the_cut_follows_the_normal_direction(self):
        m = cut_mask(self.xyz, self.point, -self.normal)
        kept = self.xyz[~m]
        assert kept[:, 2].min() >= 0, "flipping the normal cuts the other half"

    def test_an_oblique_plane(self):
        n = np.array([1.0, 1.0, 0.0])
        pts = np.array([[5, 5, 0], [-5, -5, 0], [5, -5, 0]], float)
        m = cut_mask(pts, [0, 0, 0], n)
        assert m[0] and not m[1], "the corner along +n goes, the one along -n stays"

    def test_normalises_the_normal(self):
        a = cut_mask(self.xyz, self.point, [0, 0, 1])
        b = cut_mask(self.xyz, self.point, [0, 0, 7])
        np.testing.assert_array_equal(a, b)

    def test_empty_input(self):
        assert cut_mask(np.zeros((0, 3)), self.point, self.normal).shape == (0,)

    def test_bad_shape_and_zero_normal(self):
        with pytest.raises(ValueError, match=r"\(N, 3\)"):
            cut_mask(np.zeros((4, 2)), self.point, self.normal)
        with pytest.raises(ValueError, match="zero length"):
            cut_mask(self.xyz, self.point, [0, 0, 0])


class TestCapsule:
    def test_spans_the_two_endpoints(self):
        a, b = np.array([0.0, 0, 0]), np.array([20.0, 0, 0])
        pts, _ = capsule_points(a, b, 11.0)
        assert pts[:, 0].min() < 0.5, "the cap should reach past the start"
        assert pts[:, 0].max() > 19.5, "and past the end"

    def test_radius_is_respected_on_the_shaft(self):
        a, b = np.array([0.0, 0, 0]), np.array([20.0, 0, 0])
        pts, _ = capsule_points(a, b, 11.0)
        mid = pts[np.abs(pts[:, 0] - 10.0) < 0.5]
        r = np.linalg.norm(mid[:, 1:], axis=1)
        assert r.max() == pytest.approx(11.0, abs=0.01)

    def test_works_for_any_orientation(self):
        for axis in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]):
            a = np.zeros(3)
            b = np.asarray(axis, float) * 15.0
            pts, _ = capsule_points(a, b, 8.0)
            # every shaft point is ~8 A from the axis line
            u = b / np.linalg.norm(b)
            d = pts - a
            perp = d - np.outer(d @ u, u)
            along = d @ u
            shaft = (along > 2.0) & (along < np.linalg.norm(b) - 2.0)
            if shaft.any():
                assert np.linalg.norm(perp[shaft], axis=1).max() == pytest.approx(8.0, abs=0.01)

    def test_degenerate_capsule_is_a_sphere(self):
        pts, rad = capsule_points([1, 2, 3], [1, 2, 3], 5.0)
        assert len(pts) == 1
        np.testing.assert_allclose(pts[0], [1, 2, 3])
        assert rad[0] == 5.0

    def test_rejects_a_nonpositive_radius(self):
        with pytest.raises(ValueError, match="positive"):
            capsule_points([0, 0, 0], [1, 0, 0], 0.0)

    def test_returns_matching_points_and_radii(self):
        pts, rad = capsule_points([0, 0, 0], [10, 0, 0], 6.0)
        assert len(pts) == len(rad)
        assert (rad > 0).all()


class TestHotspotParsing:
    def test_the_plus_separated_form_the_pose_stage_writes(self):
        assert parse_hotspot_residues("33+49+58+59+60") == [33, 49, 58, 59, 60]

    def test_commas_and_spaces(self):
        assert parse_hotspot_residues("33, 49 ; 58") == [33, 49, 58]

    def test_duplicates_collapse_and_order_is_kept(self):
        assert parse_hotspot_residues("58+33+58+49") == [58, 33, 49]

    def test_empty_and_none(self):
        assert parse_hotspot_residues("") == []
        assert parse_hotspot_residues(None) == []
        assert parse_hotspot_residues("   ") == []

    def test_stray_text_is_ignored(self):
        assert parse_hotspot_residues("resi 33+49") == [33, 49]


class TestPatchCentroid:
    def test_averages_the_residues_present(self):
        coords = {1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0)}
        np.testing.assert_allclose(patch_centroid(coords, [1, 3]), [2, 0, 0])

    def test_missing_residues_are_skipped(self):
        coords = {1: (0, 0, 0), 2: (2, 0, 0)}
        np.testing.assert_allclose(patch_centroid(coords, [1, 2, 99]), [1, 0, 0])

    def test_none_when_no_residue_is_present(self):
        # an epitope outside the rendered chains must vanish, not collapse to 0
        assert patch_centroid({1: (5, 5, 5)}, [77, 88]) is None

    def test_none_for_an_empty_patch(self):
        assert patch_centroid({1: (0, 0, 0)}, []) is None


class TestBinderClashDistance:
    """Closest approach between two binder axes -- the occupancy clash test."""

    @staticmethod
    def seg(*a):
        from channelbridge.render.scenes import _segment_distance
        return _segment_distance(*a)

    def test_parallel_segments(self):
        d = self.seg([0, 0, 0], [10, 0, 0], [0, 7, 0], [10, 7, 0])
        assert d == pytest.approx(7.0)

    def test_crossing_segments_touch(self):
        d = self.seg([-5, 0, 0], [5, 0, 0], [0, -5, 0], [0, 5, 0])
        assert d == pytest.approx(0.0, abs=1e-9)

    def test_skew_segments(self):
        # offset in z by 3, crossing in xy
        d = self.seg([-5, 0, 0], [5, 0, 0], [0, -5, 3], [0, 5, 3])
        assert d == pytest.approx(3.0)

    def test_endpoints_when_segments_do_not_overlap(self):
        d = self.seg([0, 0, 0], [1, 0, 0], [5, 0, 0], [6, 0, 0])
        assert d == pytest.approx(4.0)

    def test_degenerate_points(self):
        assert self.seg([0, 0, 0], [0, 0, 0], [3, 4, 0], [3, 4, 0]) == pytest.approx(5.0)

    def test_one_degenerate_one_segment(self):
        assert self.seg([0, 0, 0], [0, 0, 0], [-5, 2, 0], [5, 2, 0]) == pytest.approx(2.0)

    def test_is_symmetric(self):
        a = self.seg([0, 0, 0], [10, 1, 2], [3, 7, 1], [8, 9, 4])
        b = self.seg([3, 7, 1], [8, 9, 4], [0, 0, 0], [10, 1, 2])
        assert a == pytest.approx(b)

    def test_two_binders_closer_than_two_radii_would_clash(self):
        # capsules of radius 11 touch when their axes are within 22 A
        d = self.seg([0, 0, 0], [20, 0, 0], [0, 18, 0], [20, 18, 0])
        assert d < 2 * 11.0, "this pair should register as a clash"


class TestEnzymeSectionPoint:
    """Each enzyme is cut at its own channel, not at a shared midpoint.

    A single plane through the midpoint of the two exits cuts both enzymes
    equally only if they sit at the same depth along the view axis. On the
    real the worked example's GDH-FNR pose they do not: the shared plane kept 1296/1971 atoms
    of GDH but only 234/2352 of FNR -- 90 % of the downstream enzyme deleted,
    which reads as no section at all.
    """

    @staticmethod
    def blob(centre, n=2000, spread=8.0, seed=0):
        rng = np.random.default_rng(seed)
        return rng.normal(size=(n, 3)) * spread + np.asarray(centre, float)

    def test_it_rescues_an_enzyme_whose_bulk_is_in_front_of_its_exit(self):
        xyz = self.blob([0, 0, 20])
        exit_xyz = np.array([0.0, 0.0, -5.0])
        normal = np.array([0.0, 0.0, 1.0])
        naive = (~cut_mask(xyz, exit_xyz, normal)).sum()
        pt = enzyme_section_point(exit_xyz, xyz, normal)
        rescued = (~cut_mask(xyz, pt, normal)).sum()
        assert naive < 0.05 * len(xyz), "the naive cut should be the disaster case"
        # the plane slides forward to the deepest the band allows, so what
        # survives is keep_min -- not more, because the cut should still pass
        # as close to the channel as it can
        assert rescued >= 0.45 * len(xyz) - 1, "kept %d of %d" % (rescued, len(xyz))
        assert rescued <= 0.50 * len(xyz), "should not retreat further than needed"

    def test_it_cuts_at_the_exit_when_that_falls_inside_the_band(self):
        xyz = self.blob([0, 0, 0])
        exit_xyz = np.array([0.0, 0.0, 1.0])     # near the middle in depth
        normal = np.array([0.0, 0.0, 1.0])
        pt = enzyme_section_point(exit_xyz, xyz, normal)
        np.testing.assert_allclose(pt, exit_xyz, atol=1e-9)

    def test_an_exit_beyond_the_enzyme_still_removes_something(self):
        # the failure the selftest caught: cutting at an exit that lies past
        # every atom removes nothing, and that is not a section
        xyz = self.blob([0, 0, 0])
        exit_xyz = np.array([0.0, 0.0, 60.0])
        normal = np.array([0.0, 0.0, 1.0])
        pt = enzyme_section_point(exit_xyz, xyz, normal)
        kept = (~cut_mask(xyz, pt, normal)).sum()
        assert 0.3 * len(xyz) < kept < 0.85 * len(xyz), "kept %d of %d" % (kept, len(xyz))

    def test_it_never_cuts_shallower_than_the_exit(self):
        # cutting in front of the exit would leave the channel mouth buried
        xyz = self.blob([0, 0, 0])
        exit_xyz = np.array([0.0, 0.0, 2.0])
        normal = np.array([0.0, 0.0, 1.0])
        pt = enzyme_section_point(exit_xyz, xyz, normal, 0.9, 0.95)
        assert float(pt @ normal) >= float(exit_xyz @ normal) - 1e-9

    def test_keep_controls_how_much_survives(self):
        xyz = self.blob([0, 0, 20])
        exit_xyz = np.array([0.0, 0.0, -5.0])
        normal = np.array([0.0, 0.0, 1.0])
        few = (~cut_mask(xyz, enzyme_section_point(exit_xyz, xyz, normal, 0.3, 0.35),
                         normal)).sum()
        many = (~cut_mask(xyz, enzyme_section_point(exit_xyz, xyz, normal, 0.8, 0.85),
                          normal)).sum()
        assert few < many

    def test_bad_inputs_are_refused(self):
        xyz = self.blob([0, 0, 0])
        with pytest.raises(ValueError, match="fraction"):
            enzyme_section_point([0, 0, 0], xyz, [0, 0, 1], keep_min=0.0)
        with pytest.raises(ValueError, match="above"):
            enzyme_section_point([0, 0, 0], xyz, [0, 0, 1], 0.8, 0.5)
        with pytest.raises(ValueError, match="zero length"):
            enzyme_section_point([0, 0, 0], xyz, [0, 0, 0])
        with pytest.raises(ValueError, match="no coordinates"):
            enzyme_section_point([0, 0, 0], np.zeros((0, 3)), [0, 0, 1])
