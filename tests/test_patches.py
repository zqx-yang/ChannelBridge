"""Stage 3: residue classification and hotspot patch selection.

The measurements come from PyMOL; everything here is the arithmetic done on
them, which is where the decisions live and where they can be tested without
a renderer.
"""

from __future__ import annotations

import pytest

from channelbridge.patches import (
    PatchParams,
    Residue,
    classify,
    find_patches,
    geometry_score,
    grow_avoid,
    patch_score,
    write_patch_table,
)

P = PatchParams()


def res(resi, resn="LEU", **kw):
    kw.setdefault("sasa", 40.0)
    kw.setdefault("perp", 20.0)
    kw.setdefault("dist_mouth", 25.0)
    kw.setdefault("axial", -30.0)
    return Residue(resi=str(resi), resn=resn, **kw)


class TestClassesPartition:
    """Every residue gets exactly one class, and the order of the tests is
    the priority. Several of these pairs are deliberately contradictory --
    a residue can be in the channel AND in the shadow -- so the order is the
    whole specification."""

    def test_the_cofactor_wins_over_everything(self):
        r = res(1, sasa=0.0, perp=0.0, axial=0.0)
        assert classify(r, in_core=True, in_channel=True, in_avoid=True,
                        params=P) == "LIGAND_CORE"

    def test_channel_beats_shadow(self):
        # a channel-lining residue is usually also in the shadow cylinder;
        # calling it SHADOW would lose the fact that it lines the tunnel
        r = res(2, perp=0.0, axial=0.0)
        assert classify(r, in_core=False, in_channel=True, in_avoid=False,
                        params=P) == "CHANNEL"

    def test_shadow_beats_avoid(self):
        r = res(3, perp=0.0, axial=0.0)
        assert classify(r, in_core=False, in_channel=False, in_avoid=True,
                        params=P) == "SHADOW"

    def test_channel_and_shadow_both_need_exposure(self):
        # buried residues cannot be bound, so labelling them by geometry
        # would paint the inside of the protein
        r = res(4, sasa=1.0, perp=0.0, axial=0.0)
        assert classify(r, in_core=False, in_channel=True, in_avoid=False,
                        params=P) == "BURIED"

    def test_interface_before_buried(self):
        r = res(5, sasa=20.0, dsasa=40.0)
        assert classify(r, in_core=False, in_channel=False, in_avoid=False,
                        params=P) == "IFACE"

    def test_the_rest_is_allowed(self):
        assert classify(res(6), in_core=False, in_channel=False,
                        in_avoid=False, params=P) == "ALLOWED"

    def test_the_shadow_reaches_back_inside_the_mouth(self):
        # a binder just inside the opening blocks it as surely as one outside
        inside = res(7, perp=5.0, axial=-5.0)      # within shadow_back
        behind = res(8, perp=5.0, axial=-50.0)     # far up the channel
        assert classify(inside, in_core=False, in_channel=False,
                        in_avoid=False, params=P) == "SHADOW"
        assert classify(behind, in_core=False, in_channel=False,
                        in_avoid=False, params=P) == "ALLOWED"


class TestAvoidGrowth:
    """A hotspot is the CENTRE of a binder footprint, so keeping a binder off
    a residue means keeping hotspots further away than the residue itself."""

    def test_an_exposed_residue_spreads(self):
        rs = {"1": res(1, sasa=60.0, cb=(0.0, 0.0, 0.0)),
              "2": res(2, cb=(5.0, 0.0, 0.0)),
              "3": res(3, cb=(50.0, 0.0, 0.0))}
        got = grow_avoid(["1"], rs, P)
        assert got == {"1", "2"}, "the far one is untouched"

    def test_a_buried_residue_protects_only_itself(self):
        # growing a 12 A ball from a catalytic residue under the surface
        # deletes usable epitope to protect something nothing can reach
        rs = {"1": res(1, sasa=2.0, cb=(0.0, 0.0, 0.0)),
              "2": res(2, cb=(5.0, 0.0, 0.0))}
        assert grow_avoid(["1"], rs, P) == {"1"}

    def test_a_per_residue_radius_overrides(self):
        rs = {"1": res(1, sasa=60.0, cb=(0.0, 0.0, 0.0)),
              "2": res(2, cb=(20.0, 0.0, 0.0))}
        assert grow_avoid(["1"], rs, P) == {"1"}
        assert grow_avoid(["1"], rs, P, radius={"1": 25.0}) == {"1", "2"}

    def test_an_unknown_name_is_ignored_not_fatal(self):
        assert grow_avoid(["99"], {"1": res(1)}, P) == set()


class TestGeometryScore:
    def test_zero_inside_the_shadow_cylinder(self):
        # no amount of good chemistry redeems capping the channel
        assert geometry_score(25.0, P.shadow_r - 0.1, P) == 0.0

    def test_one_in_both_plateaus(self):
        assert geometry_score(25.0, 24.0, P) == pytest.approx(1.0)

    def test_it_falls_off_outside_the_bands(self):
        near = geometry_score(25.0, P.axis_hi + 6.0, P)
        far = geometry_score(25.0, P.axis_hi + 11.0, P)
        assert 0.0 < far < near < 1.0

    def test_the_bands_are_ranges_not_optima(self):
        # "far enough to clear the mouth, near enough to bridge" is a band;
        # a single best distance would be a claim we cannot support
        lo = geometry_score(P.reach_lo, 24.0, P)
        hi = geometry_score(P.reach_hi, 24.0, P)
        assert lo == hi == pytest.approx(1.0)


class TestPatchScore:
    def test_it_reproduces_the_published_p6a_row(self):
        # GDH rank 1 of the p6a run: the table records score 4.54 beside the
        # feature values below. This is the regression that matters -- if it
        # ever drifts, patch tables from before and after stop comparing.
        assert patch_score(anchors=3, ss_frac=1.00, phob_frac=0.50, geom=1.00,
                           mean_sasa=63, mean_n_chi=2.10) == pytest.approx(4.54, abs=0.005)

    def test_anchors_saturate_at_three(self):
        kw = dict(ss_frac=0.5, phob_frac=0.5, geom=0.5, mean_sasa=50,
                  mean_n_chi=2.0)
        assert patch_score(anchors=3, **kw) == patch_score(anchors=9, **kw)

    def test_side_chain_entropy_is_the_only_penalty(self):
        kw = dict(anchors=2, ss_frac=0.5, phob_frac=0.5, geom=0.5, mean_sasa=50)
        assert patch_score(mean_n_chi=4.0, **kw) < patch_score(mean_n_chi=1.0, **kw)


class TestFindPatches:
    def _surface(self, n=12, spacing=4.0):
        rs = {}
        for i in range(n):
            resn = "PHE" if i % 3 == 0 else "ALA"
            rs[str(i + 1)] = res(i + 1, resn=resn, ss="H",
                                 cb=(i * spacing, 0.0, 0.0), sasa=60.0)
        return rs, {str(i + 1): "ALLOWED" for i in range(n)}

    def test_only_allowed_residues_are_eligible(self):
        rs, cls = self._surface()
        cls["1"] = "CHANNEL"
        patches = find_patches(rs, cls, P)
        assert all("1" not in p.members for p in patches)

    def test_a_patch_needs_an_anchor(self):
        rs, cls = self._surface()
        for r in rs.values():
            r.resn = "ALA"          # no anchors anywhere
        assert find_patches(rs, cls, P) == []

    def test_patches_are_kept_apart(self):
        rs, cls = self._surface(n=24, spacing=3.0)
        patches = find_patches(rs, cls, P)
        for i, a in enumerate(patches):
            for b in patches[i + 1:]:
                d = sum((x - y) ** 2 for x, y in zip(a.centroid, b.centroid)) ** 0.5
                assert d > P.patch_sep

    def test_ranks_are_dense_and_one_based(self):
        rs, cls = self._surface(n=24, spacing=3.0)
        patches = find_patches(rs, cls, P)
        assert [p.rank for p in patches] == list(range(1, len(patches) + 1))
        assert patches[0].label == "patch_01"

    def test_no_more_than_asked_for(self):
        rs, cls = self._surface(n=60, spacing=3.0)
        few = PatchParams(n_patches=2)
        assert len(find_patches(rs, cls, few)) <= 2

    def test_they_come_back_best_first(self):
        rs, cls = self._surface(n=24, spacing=3.0)
        scores = [p.score for p in find_patches(rs, cls, P)]
        assert scores == sorted(scores, reverse=True)


class TestTableFormat:
    def test_the_columns_match_the_original(self, tmp_path):
        # the map session and the generator-input stage both read this file;
        # a renamed column breaks them silently
        rs, cls = TestFindPatches()._surface(n=24, spacing=3.0)
        out = write_patch_table(str(tmp_path / "p.tsv"), find_patches(rs, cls, P))
        head = open(out, encoding="utf-8").readline().rstrip("\n").split("\t")
        assert head == ["rank", "score", "geom_score", "seed", "n_res",
                        "anchors", "ss_pct", "hydrophobic_pct", "mean_n_chi",
                        "mean_sasa", "dist_mouth", "angle_mouth",
                        "perp_offset", "hotspot_residues", "sequence"]

    def test_residues_are_plus_joined_as_the_readers_expect(self, tmp_path):
        rs, cls = TestFindPatches()._surface(n=24, spacing=3.0)
        out = write_patch_table(str(tmp_path / "p.tsv"), find_patches(rs, cls, P))
        row = open(out, encoding="utf-8").read().splitlines()[1].split("\t")
        assert "+" in row[13] and all(x.isdigit() for x in row[13].split("+"))


class TestEveryMouthCastsAShadow:
    """SHADOW is tested at ALL openings, not just the primary.

    An enzyme may open more than once, and a binder capping a secondary
    mouth blocks the channel as surely as one capping the primary. The
    residue table reports the primary's geometry, which is why this is
    invisible in the table: p6a's GDH residues 100 and 103 sit 22 A off the
    primary axis -- far outside its 13 A cylinder -- and are SHADOW on the
    strength of a second mouth.
    """

    def test_a_second_opening_can_make_it_shadow(self):
        r = res(1, perp=22.0, axial=-2.5, openings=[(-2.5, 22.0), (0.3, 6.0)])
        assert classify(r, in_core=False, in_channel=False, in_avoid=False,
                        params=P) == "SHADOW"

    def test_the_primary_alone_would_not_have(self):
        r = res(2, perp=22.0, axial=-2.5, openings=[(-2.5, 22.0)])
        assert classify(r, in_core=False, in_channel=False, in_avoid=False,
                        params=P) == "ALLOWED"

    def test_no_openings_falls_back_to_the_reported_geometry(self):
        # a one-tunnel run leaves the list empty; the columns still decide
        r = res(3, perp=6.0, axial=0.0)
        assert classify(r, in_core=False, in_channel=False, in_avoid=False,
                        params=P) == "SHADOW"
