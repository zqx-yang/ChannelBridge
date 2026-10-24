"""Tunnel normalisation: the part that makes geometry build-independent."""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.tunnel import (
    BURIAL_SURFACE_MIN,
    Tunnel,
    burial_at,
    burial_verdict,
)


def straight_tunnel(n=16, length=15.0, radius=2.0, start=(0.0, 0.0, 0.0), axis=(1.0, 0.0, 0.0)):
    """A straight chain of *n* spheres spanning *length*, inner end first."""
    t = np.linspace(0.0, length, n)
    a = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    xyz = np.asarray(start, dtype=float)[None, :] + t[:, None] * a[None, :]
    return Tunnel(xyz, np.full(n, radius), source="synthetic")


def wall_of_atoms(centre=(0.0, 0.0, 0.0), n=12, spacing=1.5):
    """A cube of atoms around *centre*, standing in for buried protein."""
    g = (np.arange(n) - (n - 1) / 2.0) * spacing
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    return pts + np.asarray(centre, dtype=float)


class TestConstruction:
    def test_needs_at_least_two_spheres_to_have_a_direction(self):
        with pytest.raises(ValueError, match="at least 2"):
            Tunnel(np.zeros((1, 3)), np.ones(1))

    def test_rejects_mismatched_radii(self):
        with pytest.raises(ValueError, match="radii"):
            Tunnel(np.zeros((4, 3)), np.ones(3))

    def test_rejects_non_finite(self):
        xyz = np.zeros((3, 3))
        xyz[1, 0] = np.nan
        with pytest.raises(ValueError, match="non-finite"):
            Tunnel(xyz, np.ones(3))


class TestBasicGeometry:
    def test_span_path_length_and_spacing(self):
        t = straight_tunnel(n=16, length=15.0)
        assert t.span == pytest.approx(15.0)
        assert t.path_length == pytest.approx(15.0)
        assert t.mean_spacing == pytest.approx(1.0)

    def test_path_length_exceeds_span_on_a_curved_tunnel(self):
        ang = np.linspace(0.0, np.pi / 2, 20)
        xyz = np.column_stack([10 * np.cos(ang), 10 * np.sin(ang), np.zeros_like(ang)])
        t = Tunnel(xyz, np.full(len(ang), 2.0))
        assert t.path_length > t.span

    def test_bottleneck_is_the_minimum_radius(self):
        t = straight_tunnel(n=5)
        t = Tunnel(t.xyz, np.array([3.0, 2.5, 1.2, 2.8, 3.1]))
        assert t.bottleneck == pytest.approx(1.2)


class TestOrientation:
    def test_burial_puts_the_buried_end_first(self):
        # Protein sits at the origin; the tunnel runs out of it along +x.
        prot = wall_of_atoms(centre=(0, 0, 0))
        t = straight_tunnel(n=10, length=20.0, start=(0, 0, 0), axis=(1, 0, 0))
        # Hand it in backwards.
        reversed_t = Tunnel(t.xyz[::-1].copy(), t.radii[::-1].copy())
        out = reversed_t.oriented(prot)
        assert out.exit[0] > out.inner[0]
        assert "reversed" in out.history[-1]

    def test_already_correct_order_is_left_alone(self):
        prot = wall_of_atoms()
        t = straight_tunnel(n=10, length=20.0)
        out = t.oriented(prot)
        np.testing.assert_allclose(out.xyz, t.xyz)
        assert "kept" in out.history[-1]

    def test_falls_back_to_centroid_without_protein(self):
        t = straight_tunnel(n=10, length=20.0)
        out = t.oriented()
        assert "centroid-distance" in out.history[-1]


class TestTruncateAtSurface:
    def test_drops_the_spheres_that_left_the_surface(self):
        prot = wall_of_atoms(centre=(0, 0, 0), n=14, spacing=1.5)
        # Runs from inside the cube out to 30 A, well past it.
        t = straight_tunnel(n=61, length=30.0, start=(0, 0, 0), axis=(1, 0, 0))
        out = t.truncate_at_surface(prot)
        assert len(out) < len(t)
        assert out.burial_exit(prot) >= BURIAL_SURFACE_MIN
        # And the sphere just outside the kept one is below the threshold.
        nxt = t.xyz[len(out)]
        assert burial_at(nxt, prot)[0] < BURIAL_SURFACE_MIN

    def test_keeps_a_tunnel_that_ends_at_the_surface(self):
        prot = wall_of_atoms(n=14, spacing=1.5)
        t = straight_tunnel(n=8, length=7.0, start=(0, 0, 0), axis=(1, 0, 0))
        out = t.truncate_at_surface(prot)
        assert len(out) == len(t)
        assert "nothing to drop" in out.history[-1]

    def test_warns_and_keeps_all_when_nothing_is_buried(self):
        prot = wall_of_atoms(centre=(500, 500, 500))
        t = straight_tunnel(n=10, length=10.0)
        with pytest.warns(UserWarning, match="never reaches burial"):
            out = t.truncate_at_surface(prot)
        assert len(out) == len(t)

    def test_a_blob_beyond_the_mouth_does_not_extend_the_exit(self):
        # A neighbouring subunit or additive sitting out in the solvent also
        # raises burial. Scanning inward from the far end would accept it and
        # put the exit ~40 A out; scanning outward from the buried inner end
        # stops at the real mouth.
        prot = np.vstack(
            [
                wall_of_atoms(centre=(0, 0, 0), n=14, spacing=1.5),
                wall_of_atoms(centre=(40, 0, 0), n=8, spacing=1.0),
            ]
        )
        t = straight_tunnel(n=101, length=50.0, start=(0, 0, 0), axis=(1, 0, 0))
        out = t.truncate_at_surface(prot)
        assert out.exit[0] < 20.0
        # The naive inward scan really would have been fooled here.
        b = t.burial_profile(prot)
        assert b[-1] >= BURIAL_SURFACE_MIN

    def test_a_single_noisy_dip_does_not_truncate_early(self):
        # One sphere below threshold in the middle of a buried stretch is
        # noise, not the mouth.
        prot = wall_of_atoms(centre=(0, 0, 0), n=20, spacing=1.5)
        t = straight_tunnel(n=61, length=30.0, start=(0, 0, 0), axis=(1, 0, 0))
        b = t.burial_profile(prot)
        deep = int(np.argmax(b < BURIAL_SURFACE_MIN)) // 2  # safely inside
        out = t.truncate_at_surface(prot, persist=2.0)
        # The real crossing is far outside `deep`, so the exit is well past it.
        assert np.linalg.norm(out.exit - t.inner) > t._arc()[deep] + 2.0

    def test_warns_when_the_inner_end_is_not_buried(self):
        # Usually means the sphere order is wrong, or the exclude list removed
        # the protein instead of the ligand.
        prot = wall_of_atoms(centre=(30, 0, 0), n=14, spacing=1.5)
        t = straight_tunnel(n=31, length=15.0, start=(0, 0, 0), axis=(1, 0, 0))
        with pytest.warns(UserWarning, match="inner end is not"):
            out = t.truncate_at_surface(prot)
        assert len(out) == len(t)


class TestResample:
    def test_hits_the_requested_spacing(self):
        t = straight_tunnel(n=145, length=19.0)
        out = t.resample(1.0)
        assert out.mean_spacing == pytest.approx(1.0, abs=0.05)

    def test_preserves_both_endpoints_exactly(self):
        t = straight_tunnel(n=145, length=19.069)
        out = t.resample(1.0)
        np.testing.assert_allclose(out.inner, t.inner, atol=0)
        np.testing.assert_allclose(out.exit, t.exit, atol=0)

    def test_preserves_path_length_on_a_straight_tunnel(self):
        t = straight_tunnel(n=145, length=19.0)
        assert t.resample(1.0).path_length == pytest.approx(19.0, abs=1e-6)

    def test_interpolates_radii(self):
        t = Tunnel(np.column_stack([np.linspace(0, 10, 11), np.zeros(11), np.zeros(11)]),
                   np.linspace(1.0, 3.0, 11))
        out = t.resample(2.0)
        assert out.radii[0] == pytest.approx(1.0)
        assert out.radii[-1] == pytest.approx(3.0)
        assert np.all(np.diff(out.radii) > 0)

    def test_upsamples_a_sparse_tunnel_too(self):
        t = straight_tunnel(n=3, length=12.0)
        assert len(t.resample(1.0)) > len(t)

    def test_rejects_nonpositive_spacing(self):
        with pytest.raises(ValueError, match="positive"):
            straight_tunnel().resample(0.0)


class TestAxisPoint:
    """P1 over a fixed arc length -- the fix for the O_IC arm."""

    def test_arm_length_is_independent_of_sphere_count(self):
        sparse = straight_tunnel(n=16, length=15.0)
        dense = straight_tunnel(n=145, length=15.0)
        for t in (sparse, dense):
            p1 = t.axis_point(1.0)
            assert np.linalg.norm(t.exit - p1) == pytest.approx(1.0, abs=1e-6)

    def test_dense_and_sparse_chains_give_the_same_axis(self):
        sparse = straight_tunnel(n=16, length=15.0, axis=(1, 1, 0))
        dense = straight_tunnel(n=145, length=15.0, axis=(1, 1, 0))
        np.testing.assert_allclose(sparse.axis(1.0), dense.axis(1.0), atol=1e-9)

    def test_the_old_second_to_last_sphere_convention_is_what_it_fixes(self):
        # On a 145-sphere chain the last step is ~0.1 A; over a 1 A arm the
        # same chain gives a ten times longer baseline for the same direction.
        dense = straight_tunnel(n=145, length=15.0)
        old_arm = float(np.linalg.norm(dense.xyz[-1] - dense.xyz[-2]))
        new_arm = float(np.linalg.norm(dense.exit - dense.axis_point(1.0)))
        assert old_arm < 0.2
        assert new_arm == pytest.approx(1.0, abs=1e-6)

    def test_axis_is_a_unit_vector_pointing_outward(self):
        t = straight_tunnel(n=20, length=15.0, axis=(0, 0, 1))
        a = t.axis(1.0)
        assert np.linalg.norm(a) == pytest.approx(1.0)
        assert a[2] > 0

    def test_refuses_an_arm_longer_than_the_tunnel(self):
        t = straight_tunnel(n=5, length=2.0)
        with pytest.raises(ValueError, match="shorter than the requested"):
            t.axis_point(5.0)

    def test_can_be_asked_to_fall_back_instead(self):
        t = straight_tunnel(n=5, length=2.0)
        np.testing.assert_allclose(t.axis_point(5.0, strict=False), t.inner)

    def test_follows_a_curve_rather_than_cutting_across_it(self):
        ang = np.linspace(0.0, np.pi / 2, 50)
        xyz = np.column_stack([10 * np.cos(ang), 10 * np.sin(ang), np.zeros_like(ang)])
        t = Tunnel(xyz, np.full(len(ang), 2.0))
        p1 = t.axis_point(2.0)
        # On the arc, so ~10 A from the circle centre, not chopping the chord.
        assert np.linalg.norm(p1) == pytest.approx(10.0, abs=0.05)


class TestBurial:
    def test_bands_match_the_calibration(self):
        assert burial_verdict(0) == "past the surface"
        assert burial_verdict(14) == "past the surface"
        assert burial_verdict(15) == "surface mouth"
        assert burial_verdict(39) == "surface mouth"
        assert burial_verdict(40) == "recessed"
        assert burial_verdict(59) == "recessed"
        assert burial_verdict(60) == "buried"

    def test_counts_only_atoms_within_the_radius(self):
        prot = np.array([[0.0, 0, 0], [7.9, 0, 0], [8.1, 0, 0], [100.0, 0, 0]])
        assert burial_at([0, 0, 0], prot)[0] == 2

    def test_empty_protein_is_zero_not_an_error(self):
        assert burial_at([0, 0, 0], np.zeros((0, 3)))[0] == 0

    def test_excluded_residues_must_not_be_counted(self):
        # Counting a cofactor CAVER never searched inflates burial: the P5 bug.
        from channelbridge.io.pdb import Atom, AtomArray

        def atom(resname, x):
            return Atom(1, "C", "", resname, "A", 1, "", (x, 0.0, 0.0), 1.0, 0.0, "C", True)

        arr = AtomArray([atom("ALA", 1.0), atom("NAP", 2.0), atom("HOH", 3.0)])
        kept = arr.exclude_resnames(["NAP", "HOH", "WAT"])
        assert len(kept) == 1
        assert burial_at([0, 0, 0], kept.xyz)[0] == 1


class TestNormalisedPipeline:
    def test_normalised_applies_all_three_steps_in_order(self):
        prot = wall_of_atoms(centre=(0, 0, 0), n=14, spacing=1.5)
        t = straight_tunnel(n=151, length=30.0, start=(0, 0, 0), axis=(1, 0, 0))
        out = t.normalised(prot, spacing=1.0)
        joined = "; ".join(out.history)
        assert "oriented" in joined
        assert "truncate_at_surface" in joined
        assert "resample" in joined
        assert out.mean_spacing == pytest.approx(1.0, abs=0.1)
        assert out.burial_exit(prot) >= BURIAL_SURFACE_MIN

    def test_history_is_carried_into_the_report(self):
        prot = wall_of_atoms()
        out = straight_tunnel(n=40, length=20.0).normalised(prot)
        d = out.describe(prot)
        assert d["history"] != "raw"
        assert "burial_exit" in d and "verdict" in d
        assert d["spheres"] == len(out)


class TestWriteCaverPdb:
    """Round-tripping a normalised chain through CAVER's own format."""

    def test_roundtrip_preserves_centres_and_radii(self, tmp_path):
        t = straight_tunnel(n=16, length=15.0, radius=2.5)
        p = tmp_path / "tun_cl_001.pdb"
        assert t.write_caver_pdb(p) == 16
        back = Tunnel.from_caver_cluster(p)
        np.testing.assert_allclose(back.xyz, t.xyz, atol=1e-3)
        np.testing.assert_allclose(back.radii, t.radii, atol=1e-2)

    def test_radius_goes_in_the_b_factor_column(self, tmp_path):
        t = straight_tunnel(n=3, radius=3.14)
        p = tmp_path / "t.pdb"
        t.write_caver_pdb(p)
        line = next(ln for ln in p.read_text().splitlines() if ln.startswith("ATOM"))
        assert line[54:60].strip() == ""          # occupancy blank, as CAVER writes it
        assert float(line[60:66]) == pytest.approx(3.14)

    def test_written_with_lf_endings(self, tmp_path):
        from channelbridge.io.pdb import has_crlf

        p = tmp_path / "t.pdb"
        straight_tunnel().write_caver_pdb(p)
        assert has_crlf(p) is False

    def test_order_is_preserved_so_the_exit_stays_the_exit(self, tmp_path):
        prot = wall_of_atoms(centre=(0, 0, 0), n=14, spacing=1.5)
        t = straight_tunnel(n=61, length=30.0, start=(0, 0, 0), axis=(1, 0, 0))
        norm = t.normalised(prot, spacing=1.0)
        p = tmp_path / "t.pdb"
        norm.write_caver_pdb(p)
        back = Tunnel.from_caver_cluster(p)
        np.testing.assert_allclose(back.exit, norm.exit, atol=1e-3)
        np.testing.assert_allclose(back.inner, norm.inner, atol=1e-3)

    def test_resampling_to_1A_makes_the_last_sphere_step_a_1A_arm(self, tmp_path):
        # Why this writer exists: a reader implementing the published
        # "second-to-last sphere" construction then gets a 1 A arm for free.
        prot = wall_of_atoms(centre=(0, 0, 0), n=14, spacing=1.5)
        t = straight_tunnel(n=151, length=30.0, start=(0, 0, 0), axis=(1, 0, 0))
        norm = t.normalised(prot, spacing=1.0)
        p = tmp_path / "t.pdb"
        norm.write_caver_pdb(p)
        back = Tunnel.from_caver_cluster(p)
        last_step = float(np.linalg.norm(back.xyz[-1] - back.xyz[-2]))
        assert last_step == pytest.approx(1.0, abs=0.1)
