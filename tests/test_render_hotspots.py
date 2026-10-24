"""The overall hotspot map: tables, gradient, centroids, palette.

The PyMOL half is covered by ``hotspots.py --selftest``; everything here is
the part that decides *what* gets drawn, which is where a wrong answer is
invisible in the finished picture.
"""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.render.hotspots import (
    CLASS_ORDER,
    ClassPalette,
    HotspotMapSpec,
    Patch,
    patch_centroids,
    patch_gradient,
    read_patch_table,
    read_residue_classes,
)

CLASS_TSV = (
    "chain\tresi\tresn\tss\tclass\tsasa\n"
    "B\t1\tMET\tL\tIFACE\t81.1\n"
    "B\t19\tHIS\tL\tAVOID\t233.2\n"
    "B\t33\tLYS\tH\tALLOWED\t120.0\n"
    "B\t49\tLEU\tH\tALLOWED\t60.0\n"
    "B\t58\tGLU\tH\tCHANNEL\t12.0\n"
    "B\t60\tILE\tH\tSHADOW\t44.0\n"
)

PATCH_TSV = (
    "rank\tscore\tgeom_score\tn_res\thotspot_residues\tsequence\n"
    "1\t4.54\t1.00\t3\t33+49+58\tLYS33 LEU49 GLU58\n"
    "2\t4.10\t0.90\t2\t60+61\tILE60 ALA61\n"
    "3\t3.80\t0.80\t2\t78+81\tSER78 LYS81\n"
)


@pytest.fixture()
def classes_file(tmp_path):
    p = tmp_path / "GDH_residue_classes.tsv"
    p.write_text(CLASS_TSV)
    return str(p)


@pytest.fixture()
def patch_file(tmp_path):
    p = tmp_path / "GDH_hotspot_patches.tsv"
    p.write_text(PATCH_TSV)
    return str(p)


class TestReadResidueClasses:
    def test_it_groups_by_class(self, classes_file):
        t = read_residue_classes(classes_file)
        assert set(t) == {"IFACE", "AVOID", "ALLOWED", "CHANNEL", "SHADOW"}
        assert sorted(t["ALLOWED"]) == [33, 49]

    def test_it_keeps_the_residue_name(self, classes_file):
        assert read_residue_classes(classes_file)["CHANNEL"][58] == "GLU"

    def test_every_class_it_produces_has_a_default_colour(self, classes_file):
        pal = ClassPalette()
        for cls in read_residue_classes(classes_file):
            assert cls in pal.classes, "%s would be painted with nothing" % cls

    def test_a_missing_file_is_an_error_not_an_empty_map(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_residue_classes(str(tmp_path / "nope.tsv"))

    def test_the_wrong_table_is_named(self, tmp_path):
        p = tmp_path / "wrong.tsv"
        p.write_text("a\tb\n1\t2\n")
        with pytest.raises(ValueError, match="residue-class table"):
            read_residue_classes(str(p))

    def test_rows_with_an_unparseable_residue_number_are_skipped(self, tmp_path):
        p = tmp_path / "t.tsv"
        p.write_text("chain\tresi\tresn\tclass\nB\t12A\tALA\tAVOID\nB\t13\tALA\tAVOID\n")
        assert sorted(read_residue_classes(str(p))["AVOID"]) == [13]


class TestReadPatchTable:
    def test_rank_order_and_residues(self, patch_file):
        ps = read_patch_table(patch_file)
        assert [p.rank for p in ps] == ["1", "2", "3"]
        assert ps[0].residues == [33, 49, 58]

    def test_scores_and_sequences_come_along(self, patch_file):
        ps = read_patch_table(patch_file)
        assert ps[0].score == pytest.approx(4.54)
        assert ps[0].sequence.startswith("LYS33")

    def test_the_label_is_zero_padded_so_sorting_keeps_the_rank_order(self):
        labels = [Patch(str(r), []).label for r in (1, 2, 10)]
        assert labels == ["patch_01", "patch_02", "patch_10"]
        assert sorted(labels) == labels

    def test_a_table_written_out_of_order_comes_back_ranked(self, tmp_path):
        p = tmp_path / "t.tsv"
        p.write_text("rank\thotspot_residues\n3\t1\n1\t2\n2\t3\n")
        assert [x.rank for x in read_patch_table(str(p))] == ["1", "2", "3"]

    def test_the_wrong_table_is_named(self, tmp_path):
        p = tmp_path / "wrong.tsv"
        p.write_text("rank\tscore\n1\t2\n")
        with pytest.raises(ValueError, match="hotspot-patch table"):
            read_patch_table(str(p))

    def test_a_missing_file_is_an_error(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_patch_table(str(tmp_path / "nope.tsv"))


class TestPatchGradient:
    def test_one_shade_per_patch(self):
        assert len(patch_gradient(4)) == 4

    def test_the_shades_deepen_down_the_list(self):
        shades = patch_gradient(4, "red")
        lightness = [sum(s) for s in shades]
        assert lightness == sorted(lightness, reverse=True), \
            "patch_01 palest, patch_04 deepest -- an increasing shade"

    def test_the_order_can_be_flipped_for_best_is_darkest(self):
        shades = patch_gradient(4, "red", best_is_darkest=True)
        lightness = [sum(s) for s in shades]
        assert lightness == sorted(lightness), "rank 1 deepest instead"

    def test_every_shade_is_still_the_same_hue(self):
        # mixed with white, not scaled to black: red stays the strongest channel
        for s in patch_gradient(5, "red"):
            assert s[0] >= s[1] and s[0] >= s[2]

    def test_a_single_patch_gets_the_deep_end_not_the_palest(self):
        (only,) = patch_gradient(1, "red")
        assert only == pytest.approx(patch_gradient(4, "red")[-1])

    def test_none_for_no_patches(self):
        assert patch_gradient(0) == []

    def test_any_base_hue_works(self):
        blues = patch_gradient(3, "#0000ff")
        for s in blues:
            assert s[2] >= s[0] and s[2] >= s[1]

    def test_shades_are_distinguishable_from_each_other(self):
        shades = patch_gradient(6, "red")
        gaps = [abs(sum(a) - sum(b)) for a, b in zip(shades, shades[1:])]
        assert min(gaps) > 0.05, "two patches must not come out the same colour"


class TestClassPalette:
    def test_defaults_follow_the_p5d_map(self):
        p = ClassPalette()
        assert p.classes["CHANNEL"] == "magenta"
        assert p.classes["AVOID"] == "grey60"
        assert p.classes["SHADOW"] == "pink"
        assert p.classes["ALLOWED"] == "cyan"

    def test_it_does_not_split_hydrophobes_out_of_the_surface(self):
        assert not any("HYDROPH" in c for c in CLASS_ORDER)

    def test_an_override_keeps_the_classes_it_does_not_name(self):
        p = ClassPalette.from_dict({"classes": {"CHANNEL": "yellow"}})
        assert p.classes["CHANNEL"] == "yellow"
        assert p.classes["AVOID"] == "grey60", "the rest must survive"

    def test_an_unknown_class_is_an_error_not_a_silent_no_op(self):
        with pytest.raises(ValueError, match="unknown residue class"):
            ClassPalette.from_dict({"classes": {"CHANEL": "yellow"}})

    def test_an_unknown_key_is_an_error(self):
        with pytest.raises(ValueError, match="unknown key"):
            ClassPalette.from_dict({"patch_bass": "red"})

    def test_a_bad_colour_is_caught_at_construction(self):
        with pytest.raises(ValueError):
            ClassPalette.from_dict({"patch_base": "octarine"})

    def test_the_gradient_ends_must_be_fractions(self):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            ClassPalette(patch_light=1.5)

    def test_the_flip_is_settable_from_a_pair_file(self):
        assert ClassPalette.from_dict({"best_is_darkest": True}).best_is_darkest


class TestPatchCentroids:
    def test_the_centroid_of_the_residues_present(self):
        coords = {1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0)}
        got = patch_centroids(coords, [Patch("1", [1, 3])])
        np.testing.assert_allclose(got["patch_01"], [2, 0, 0])

    def test_missing_residues_are_skipped(self):
        coords = {1: (0, 0, 0), 2: (2, 0, 0)}
        got = patch_centroids(coords, [Patch("1", [1, 2, 99])])
        np.testing.assert_allclose(got["patch_01"], [1, 0, 0])

    def test_a_patch_with_no_residue_present_is_left_out(self):
        # better no envelope than a green ball at the origin
        assert patch_centroids({1: (5, 5, 5)}, [Patch("1", [77])]) == {}

    def test_one_entry_per_patch(self):
        coords = {r: (float(r), 0.0, 0.0) for r in range(1, 11)}
        got = patch_centroids(coords, [Patch("1", [1, 2]), Patch("2", [9, 10])])
        assert set(got) == {"patch_01", "patch_02"}


class TestSpec:
    def test_chains_are_cleaned_and_named(self):
        s = HotspotMapSpec("c.pdb", up_name="GDH", dn_name="FNR",
                           up_chains=["B", "", "C"], dn_chains=["N"])
        assert s.chains("up") == ["B", "C"]
        assert s.name("dn") == "FNR"

    def test_the_binder_radius_matches_the_pose_stage_default(self):
        assert HotspotMapSpec("c.pdb").binder_radius == 11.0


class TestColourSeparation:
    """The drawn colours, not the named ones: lightening changes them."""

    def test_the_patches_are_tellable_from_every_class(self):
        from channelbridge.render.hotspots import ClassPalette

        clashes = ClassPalette().indistinguishable(6)
        mixed = [c for c in clashes if c[0].startswith("patch_") or c[1].startswith("patch_")]
        assert mixed == [], "a patch must not look like a residue class: %s" % mixed

    def test_a_gradient_neighbour_is_not_reported_as_a_clash(self):
        from channelbridge.render.hotspots import ClassPalette

        clashes = ClassPalette().indistinguishable(6)
        assert not any(a.startswith("patch_") and b.startswith("patch_")
                       for a, b, _ in clashes), "the gradient is deliberate"

    def test_only_the_background_classes_are_lightened(self):
        from channelbridge.render.hotspots import ClassPalette
        from channelbridge.render.style import parse_color

        p = ClassPalette()
        # a constraint the user must see keeps its full colour
        assert p.final_class_rgb("CHANNEL") == parse_color("magenta")
        # one that means "nothing special here" does not
        assert p.final_class_rgb("AVOID") != parse_color("grey60")

    def test_lightening_everything_pushes_the_pale_classes_together(self):
        # The reason background_classes exists. Not "more clashing pairs" --
        # uniform lightening actually reports fewer -- but a WORSE worst pair:
        # wheat, grey60 and pink all converge on pale, and BURIED/SHADOW go
        # from 0.25 apart to 0.17, which on a shaded surface is one colour.
        from channelbridge.render.hotspots import CLASS_ORDER, ClassPalette

        uniform = ClassPalette(background_classes=tuple(CLASS_ORDER)).indistinguishable(6)
        selective = ClassPalette().indistinguishable(6)
        assert uniform[0][2] < selective[0][2], (
            "the closest pair must be further apart when only the background "
            "classes are lightened (%s vs %s)" % (selective[0], uniform[0]))

    def test_an_unknown_background_class_is_an_error(self):
        from channelbridge.render.hotspots import ClassPalette

        with pytest.raises(ValueError, match="unknown background class"):
            ClassPalette(background_classes=("BURIED", "CHANEL"))


class TestNeutralTwin:
    """The patch scenes paint everything that is not a patch one flat colour.

    The classification is seven colours; the patches are what the user is
    judging. Shown together the rainbow wins -- so `patches` and `envelopes`
    use a neutral copy while `classes` and `map` keep the classification.
    """

    def test_the_neutral_colour_is_a_real_colour(self):
        from channelbridge.render.hotspots import ClassPalette
        from channelbridge.render.style import parse_color

        assert parse_color(ClassPalette().neutral)

    def test_it_is_overridable(self):
        from channelbridge.render.hotspots import ClassPalette

        assert ClassPalette.from_dict({"neutral": "white"}).neutral == "white"

    def test_a_bad_neutral_is_caught_at_construction(self):
        from channelbridge.render.hotspots import ClassPalette

        with pytest.raises(ValueError):
            ClassPalette(neutral="greige")

    def test_every_patch_shade_stands_out_against_the_neutral(self):
        import math

        from channelbridge.render.hotspots import ClassPalette, patch_gradient
        from channelbridge.render.style import parse_color

        p = ClassPalette()
        bg = parse_color(p.neutral)
        for i, shade in enumerate(patch_gradient(6, p.patch_base, p.patch_light,
                                                 p.patch_dark), 1):
            d = math.dist(bg, shade)
            assert d > 0.40, "patch_%02d is only %.2f from the neutral surface" % (i, d)

    def test_background_classes_survive_from_dict(self):
        # regression: `neutral` and `background_classes` were both briefly in
        # the float-coercion loop, where float() rejects a colour and a list
        from channelbridge.render.hotspots import ClassPalette

        p = ClassPalette.from_dict({"background_classes": "BURIED, AVOID",
                                    "neutral": "white"})
        assert p.background_classes == ("BURIED", "AVOID")
        assert p.neutral == "white"
