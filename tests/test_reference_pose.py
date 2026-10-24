"""The reference placement the hotspot map is drawn in.

It is computed before any pose search exists, so nothing here may depend on
a poses table. The cases are built from geometry whose answer is known: two
slabs with a known gap, a channel pointing a known way, a binder of a known
length.
"""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.reference_pose import (
    binder_blobs,
    min_interprotein_distance,
    occupancy_grid,
    patch_outward_normal,
    place_at_oic180,
    rotation_about,
    rotation_taking,
)


class TestRotations:
    def test_taking_one_vector_onto_another(self):
        R = rotation_taking([1, 0, 0], [0, 1, 0])
        np.testing.assert_allclose(R @ [1.0, 0, 0], [0, 1, 0], atol=1e-9)

    def test_it_is_a_proper_rotation(self):
        R = rotation_taking([1, 2, 3], [-3, 1, 1])
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-9)
        assert np.linalg.det(R) == pytest.approx(1.0)

    def test_identical_vectors_give_the_identity(self):
        np.testing.assert_allclose(rotation_taking([0, 0, 1], [0, 0, 2]), np.eye(3), atol=1e-12)

    def test_antiparallel_vectors_are_a_half_turn_not_a_crash(self):
        # the cross product vanishes here; the naive formula divides by zero
        R = rotation_taking([0, 0, 1], [0, 0, -1])
        np.testing.assert_allclose(R @ [0.0, 0, 1], [0, 0, -1], atol=1e-9)
        assert np.linalg.det(R) == pytest.approx(1.0)

    def test_a_zero_vector_is_refused(self):
        with pytest.raises(ValueError, match="zero-length"):
            rotation_taking([0, 0, 0], [1, 0, 0])

    def test_rotation_about_an_axis_leaves_the_axis_alone(self):
        R = rotation_about([0, 0, 1], 37.0)
        np.testing.assert_allclose(R @ [0.0, 0, 1], [0, 0, 1], atol=1e-9)

    def test_rotation_about_is_right_handed(self):
        R = rotation_about([0, 0, 1], 90.0)
        np.testing.assert_allclose(R @ [1.0, 0, 0], [0, 1, 0], atol=1e-9)

    def test_a_zero_axis_is_refused(self):
        with pytest.raises(ValueError, match="zero length"):
            rotation_about([0, 0, 0], 10.0)


class TestOccupancyGrid:
    def test_a_point_is_occupied_and_so_is_its_clearance_shell(self):
        grid, origin, sp = occupancy_grid(np.zeros((1, 3)), clearance=3.0, spacing=1.0)
        # the grid spans the atoms plus the clearance and no more, so anything
        # further out is simply off the array -- which counts as unoccupied
        idx = np.rint((np.array([[0.0, 0, 0], [2.0, 0, 0], [4.0, 0, 0]]) - origin) / sp).astype(int)
        assert grid[idx[0, 0], idx[0, 1], idx[0, 2]]
        assert grid[idx[1, 0], idx[1, 1], idx[1, 2]], "2 A away is inside a 3 A shell"
        assert not grid[idx[2, 0], idx[2, 1], idx[2, 2]], "4 A away is not"

    def test_the_shell_is_a_ball_not_a_cube(self):
        grid, origin, sp = occupancy_grid(np.zeros((1, 3)), clearance=3.0, spacing=1.0)
        corner = np.rint((np.array([[3.0, 3.0, 3.0]]) - origin) / sp).astype(int)[0]
        assert not grid[corner[0], corner[1], corner[2]], "the cube corner is 5.2 A away"

    def test_it_rejects_bad_input(self):
        with pytest.raises(ValueError, match=r"\(N, 3\)"):
            occupancy_grid(np.zeros((4, 2)), 3.0)
        with pytest.raises(ValueError, match="positive"):
            occupancy_grid(np.zeros((4, 3)), 3.0, spacing=0.0)


class TestMinDistance:
    def test_two_points(self):
        assert min_interprotein_distance([[0.0, 0, 0]], [[3.0, 4, 0]]) == pytest.approx(5.0)

    def test_it_finds_the_closest_pair_not_the_first(self):
        a = [[0.0, 0, 0], [10.0, 0, 0]]
        b = [[20.0, 0, 0], [11.0, 0, 0]]
        assert min_interprotein_distance(a, b) == pytest.approx(1.0)

    def test_an_empty_set_is_infinitely_far(self):
        assert min_interprotein_distance(np.zeros((0, 3)), [[0.0, 0, 0]]) == float("inf")

    def test_chunking_does_not_change_the_answer(self):
        rng = np.random.default_rng(0)
        a, b = rng.normal(size=(300, 3)), rng.normal(size=(200, 3)) + 5
        assert min_interprotein_distance(a, b, chunk=17) == pytest.approx(
            min_interprotein_distance(a, b, chunk=4096))


def _slab(centre, half=(6.0, 6.0, 6.0), n=7):
    """A small block of atoms, as a stand-in for an enzyme."""
    g = [np.linspace(-h, h, n) for h in half]
    pts = np.stack(np.meshgrid(*g, indexing="ij"), -1).reshape(-1, 3)
    return pts + np.asarray(centre, float)


class TestPlaceAtOIC180:
    """Two blocks, each with a channel along x, pushed apart until they fit."""

    def setup_method(self):
        self.up = _slab([0.0, 0, 0])
        self.up_exit = np.array([6.0, 0.0, 0.0])          # on the +x face
        self.up_axis_point = np.array([5.0, 0.0, 0.0])    # one arm inside
        self.dn = _slab([0.0, 0, 0])
        self.dn_exit = np.array([-6.0, 0.0, 0.0])         # on the -x face
        self.dn_axis_point = np.array([-5.0, 0.0, 0.0])

    def place(self, **kw):
        return place_at_oic180(
            self.up, self.up_exit, self.up_axis_point,
            self.dn, self.dn_exit, self.dn_axis_point, **kw)

    def test_it_finds_a_placement_and_reports_it_as_oic_180(self):
        p = self.place(d_min=1.0, d_max=20.0, d_step=1.0, clash_cut=2.0)
        assert p is not None
        assert p.o_ic == 180.0

    def test_the_downstream_exit_lands_on_the_upstream_channel_ray(self):
        p = self.place(d_min=1.0, d_max=20.0, clash_cut=2.0)
        moved_exit = p.apply(self.dn_exit[None, :])[0]
        u = (self.up_exit - self.up_axis_point) / np.linalg.norm(self.up_exit - self.up_axis_point)
        np.testing.assert_allclose(moved_exit, self.up_exit + p.d_ic * u, atol=1e-6)

    def test_the_result_is_actually_clash_free(self):
        p = self.place(d_min=1.0, d_max=20.0, clash_cut=2.0)
        assert min_interprotein_distance(self.up, p.apply(self.dn)) >= 2.0
        assert p.min_distance == pytest.approx(
            min_interprotein_distance(self.up, p.apply(self.dn)))

    def test_it_returns_the_smallest_clash_free_distance(self):
        p = self.place(d_min=1.0, d_max=20.0, d_step=1.0, clash_cut=2.0)
        # One step closer must not CLEAR the cut. "<=" rather than "<": the
        # voxel grid is conservative to within one voxel, so a placement that
        # exactly grazes the cut can be passed over. That costs at most one
        # d_step of D_IC and never reports a placement that is not clash-free.
        u = np.array([1.0, 0, 0])
        closer = self.up_exit + (p.d_ic - 1.0) * u
        shifted = p.apply(self.dn) - (p.apply(self.dn_exit[None, :])[0] - closer)
        assert min_interprotein_distance(self.up, shifted) <= 2.0

    def test_the_downstream_channel_points_back_at_the_upstream_exit(self):
        p = self.place(d_min=1.0, d_max=20.0, clash_cut=2.0)
        e = p.apply(self.dn_exit[None, :])[0]
        a = p.apply(self.dn_axis_point[None, :])[0]
        out = (e - a) / np.linalg.norm(e - a)
        np.testing.assert_allclose(out, [-1.0, 0, 0], atol=1e-6)

    def test_do_falls_as_the_distance_grows(self):
        near = self.place(d_min=1.0, d_max=30.0, clash_cut=2.0)
        far = self.place(d_min=near.d_ic + 5.0, d_max=30.0, clash_cut=2.0)
        assert far.do < near.do, "DO is monotonic in D_IC, so smallest clash-free is best"

    def test_none_when_nothing_fits_in_the_range(self):
        assert self.place(d_min=1.0, d_max=2.0, clash_cut=2.0) is None

    def test_a_channel_with_no_direction_is_refused(self):
        with pytest.raises(ValueError, match="upstream channel"):
            place_at_oic180(self.up, self.up_exit, self.up_exit,
                            self.dn, self.dn_exit, self.dn_axis_point)
        with pytest.raises(ValueError, match="downstream channel"):
            place_at_oic180(self.up, self.up_exit, self.up_axis_point,
                            self.dn, self.dn_exit, self.dn_exit)

    def test_the_description_names_the_numbers(self):
        p = self.place(d_min=1.0, d_max=20.0, clash_cut=2.0)
        text = p.describe()
        assert "D_IC" in text and "O_IC 180" in text and "DO" in text


class TestPatchNormal:
    def test_it_points_away_from_the_protein(self):
        xyz = _slab([0.0, 0, 0])
        n = patch_outward_normal(xyz, [6.0, 0.0, 0.0])
        assert n[0] > 0.9, "a patch on the +x face points along +x"

    def test_it_is_a_unit_vector(self):
        xyz = _slab([0.0, 0, 0])
        assert np.linalg.norm(patch_outward_normal(xyz, [0.0, 6.0, 0.0])) == pytest.approx(1.0)

    def test_it_uses_the_local_neighbourhood(self):
        # an L-shaped protein: the whole-protein centroid and the local mass
        # point different ways, and the local one is what a binder sees
        arm_a = _slab([0.0, 0, 0], half=(4, 4, 4), n=5)
        arm_b = _slab([20.0, 0, 0], half=(4, 4, 4), n=5)
        xyz = np.vstack([arm_a, arm_b])
        local = patch_outward_normal(xyz, [-4.0, 0.0, 0.0], radius=6.0)
        assert local[0] < -0.5, "away from the near arm"

    def test_a_centroid_at_the_centre_of_mass_still_returns_a_direction(self):
        xyz = _slab([0.0, 0, 0])
        n = patch_outward_normal(xyz, [0.0, 0.0, 0.0], radius=100.0)
        assert np.isfinite(n).all()
        assert np.linalg.norm(n) == pytest.approx(1.0)


class TestBinderBlobs:
    def test_the_minimum_binder_is_four_blobs(self):
        # the pose stage's own rule: a span of `binder_gap` (6 A) at 2 A steps
        pts, rad = binder_blobs([0, 0, 0], [1, 0, 0], length=6.0, step=2.0)
        assert len(pts) == 4
        assert rad[0] == pytest.approx(11.0)

    def test_they_span_exactly_the_binder_length(self):
        pts, _ = binder_blobs([0, 0, 0], [1, 0, 0], length=6.0)
        assert np.linalg.norm(pts[-1] - pts[0]) == pytest.approx(6.0)

    def test_the_first_blob_sits_on_the_patch(self):
        pts, _ = binder_blobs([3, 4, 5], [0, 0, 1], length=6.0)
        np.testing.assert_allclose(pts[0], [3, 4, 5])

    def test_they_follow_the_direction_given(self):
        pts, _ = binder_blobs([0, 0, 0], [0, -2, 0], length=6.0)
        np.testing.assert_allclose(pts[-1], [0, -6, 0], atol=1e-9)

    def test_every_patch_gets_the_same_number_so_they_compare(self):
        a = binder_blobs([0, 0, 0], [1, 0, 0])[0]
        b = binder_blobs([50, 10, 3], [0.3, -0.4, 0.9])[0]
        assert len(a) == len(b)

    def test_a_longer_binder_gets_more_blobs(self):
        assert len(binder_blobs([0, 0, 0], [1, 0, 0], length=20.0)[0]) == 11

    def test_never_fewer_than_two(self):
        assert len(binder_blobs([0, 0, 0], [1, 0, 0], length=0.0)[0]) == 2

    def test_a_zero_direction_is_refused(self):
        with pytest.raises(ValueError, match="zero length"):
            binder_blobs([0, 0, 0], [0, 0, 0])
