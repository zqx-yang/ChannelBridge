"""Colours, transparency and render settings -- defaults and overrides."""

from __future__ import annotations

import pytest

from channelbridge.render.style import (
    BACKGROUNDS,
    DEFAULT_STYLE,
    DESIGNED_BINDER_HEX,
    FLAT_SETTINGS,
    NAMED_COLORS,
    RAY_SETTINGS,
    ColorSpec,
    Palette,
    RenderStyle,
    Transparency,
    hex_to_rgb,
    lighten,
    parse_color,
)


class TestColorParsing:
    def test_hex_with_and_without_hash(self):
        assert parse_color("#0065bd") == parse_color("0065bd")

    def test_the_binder_blue(self):
        r, g, b = hex_to_rgb(DESIGNED_BINDER_HEX)
        assert (round(r, 4), round(g, 4), round(b, 4)) == (0.0, 0.3961, 0.7412)

    def test_names(self):
        assert parse_color("green") == (0.0, 1.0, 0.0)
        assert parse_color(" GREEN ") == (0.0, 1.0, 0.0)
        assert parse_color("gray") == parse_color("grey")

    def test_floats_0_to_1(self):
        assert parse_color([0.0, 0.5, 1.0]) == (0.0, 0.5, 1.0)

    def test_ints_0_to_255(self):
        # Resolved per-colour, not per-channel: [255,0,0] is red, [1,0,0] too.
        assert parse_color([255, 0, 0]) == (1.0, 0.0, 0.0)
        assert parse_color([1, 0, 0]) == (1.0, 0.0, 0.0)

    @pytest.mark.parametrize("bad", ["#fff", "nosuchcolour", [1, 2], [0, 0, 300],
                                     [-1, 0, 0], 7, None])
    def test_rejects_malformed(self, bad):
        with pytest.raises(ValueError):
            parse_color(bad)

    def test_an_unknown_name_is_an_error_not_a_silent_default(self):
        with pytest.raises(ValueError):
            parse_color("chartreuse")


class TestLighten:
    def test_zero_is_a_no_op_and_one_is_white(self):
        assert lighten((0.2, 0.4, 0.6), 0.0) == (0.2, 0.4, 0.6)
        assert lighten((0.2, 0.4, 0.6), 1.0) == (1.0, 1.0, 1.0)

    def test_matches_the_documented_formula(self):
        assert lighten((0.0, 1.0, 0.0), 0.55) == pytest.approx((0.55, 1.0, 0.55))

    @pytest.mark.parametrize("amount", [-0.1, 1.1])
    def test_rejects_out_of_range_amount(self, amount):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            lighten((0.2, 0.2, 0.2), amount)


class TestDefaults:
    """The current look, which must stay the default."""

    def test_upstream_is_green_downstream_orange(self):
        up = DEFAULT_STYLE.palette.rgb("upstream")
        dn = DEFAULT_STYLE.palette.rgb("downstream")
        assert up[1] > up[0] and up[1] > up[2]
        assert dn[0] > dn[2] and dn[1] > dn[2]

    def test_binder_is_the_project_blue(self):
        assert DEFAULT_STYLE.palette.base("binder") == hex_to_rgb(DESIGNED_BINDER_HEX)

    def test_ligands_are_not_lightened(self):
        """Cofactors read through two translucent surfaces, so they keep
        their full colour. Purple since 2026-10-27: the D_IC/O_IC annotation
        took pure magenta and the two were then the same colour."""
        assert DEFAULT_STYLE.palette.spec("ligand").lighten == 0.0
        assert DEFAULT_STYLE.palette.rgb("ligand") == (0.55, 0.1, 0.9)

    def test_tunnel_is_magenta_and_unlightened(self):
        """Magenta since 2026-10-27: everything to do with the channel is
        one colour, and lightening a 30 %-opacity fill loses it."""
        r, g, b = DEFAULT_STYLE.palette.rgb("tunnel")
        assert r > 0.8 and g < 0.2 and b > 0.8
        assert DEFAULT_STYLE.palette.spec("tunnel").lighten == 0.0

    def test_both_enzymes_are_transparent_and_nothing_is_opaque(self):
        assert DEFAULT_STYLE.transparency_for("upstream") == pytest.approx(0.65)
        assert DEFAULT_STYLE.transparency_for("downstream") == pytest.approx(0.65)
        assert DEFAULT_STYLE.transparency.opaque_roles == []
        assert DEFAULT_STYLE.needs_occlusion_clipping is False

    def test_the_tunnel_is_more_transparent_than_the_proteins(self):
        assert DEFAULT_STYLE.transparency_for("tunnel") == pytest.approx(0.75)

    def test_all_default_roles_are_distinguishable(self):
        assert DEFAULT_STYLE.palette.indistinguishable() == []

    def test_background_is_white(self):
        assert DEFAULT_STYLE.background == "white"

    def test_summary_mentions_the_key_facts(self):
        s = DEFAULT_STYLE.summary()
        assert "upstream" in s and "tunnel" in s
        assert "tunnels are the subject" in s


class TestColorOverrides:
    def test_a_bare_colour_keeps_the_role_lighten(self):
        # Changing a hue must not silently change how washed out it is.
        p = Palette().with_overrides({"upstream": "#4caf50"})
        assert p.spec("upstream").lighten == Palette().spec("upstream").lighten
        assert p.base("upstream") == hex_to_rgb("#4caf50")

    def test_a_mapping_can_set_both(self):
        p = Palette().with_overrides({"downstream": {"color": "teal", "lighten": 0.3}})
        assert p.base("downstream") == NAMED_COLORS["teal"]
        assert p.spec("downstream").lighten == pytest.approx(0.3)

    def test_a_mapping_can_set_lighten_alone(self):
        p = Palette().with_overrides({"ligand": {"lighten": 0.4}})
        assert p.base("ligand") == (0.55, 0.1, 0.9)
        assert p.spec("ligand").lighten == pytest.approx(0.4)

    def test_british_spelling_accepted(self):
        p = Palette().with_overrides({"upstream": {"colour": "purple"}})
        assert p.base("upstream") == NAMED_COLORS["purple"]

    def test_both_spellings_at_once_is_an_error(self):
        with pytest.raises(ValueError, match="either"):
            Palette().with_overrides({"upstream": {"color": "red", "colour": "red"}})

    def test_a_misspelt_role_is_an_error(self):
        # Silently keeping the default would be worse than failing.
        with pytest.raises(ValueError, match="unknown colour role"):
            Palette().with_overrides({"upsteam": "green"})

    def test_an_unknown_colour_key_is_an_error(self):
        with pytest.raises(ValueError, match="unknown colour key"):
            Palette().with_overrides({"upstream": {"hue": "green"}})

    def test_overrides_do_not_mutate_the_default(self):
        before = Palette().rgb("upstream")
        Palette().with_overrides({"upstream": "red"})
        assert Palette().rgb("upstream") == before
        assert DEFAULT_STYLE.palette.rgb("upstream") == before

    def test_empty_overrides_return_the_same_palette(self):
        p = Palette()
        assert p.with_overrides(None) is p
        assert p.with_overrides({}) is p

    def test_indistinguishable_flags_a_clash_the_user_created(self):
        p = Palette().with_overrides({"downstream": {"color": "green", "lighten": 0.55}})
        assert "downstream/upstream" in p.indistinguishable()

    def test_indistinguishable_is_a_warning_not_an_error(self):
        # Two identical roles are allowed; the caller decides what to do.
        p = Palette().with_overrides({"copy": {"color": "green", "lighten": 0.55}})
        assert p.indistinguishable() != []
        assert p.rgb("copy") == p.rgb("upstream")


class TestTransparencyOverrides:
    def test_default_covers_every_role(self):
        t = Transparency()
        assert t.for_role("upstream") == pytest.approx(0.65)
        assert t.for_role("a_role_added_later") == pytest.approx(0.65)

    def test_a_bare_number_sets_the_default(self):
        t = Transparency().with_overrides(0.3)
        assert t.for_role("upstream") == pytest.approx(0.3)
        # the per-role tunnel exception survives
        assert t.for_role("tunnel") == pytest.approx(0.75)

    def test_a_mapping_sets_default_and_roles(self):
        t = Transparency().with_overrides(
            {"default": 0.5, "upstream": 0.2, "tunnel": 0.9}, roles=["upstream", "tunnel"])
        assert t.for_role("downstream") == pytest.approx(0.5)
        assert t.for_role("upstream") == pytest.approx(0.2)
        assert t.for_role("tunnel") == pytest.approx(0.9)

    def test_zero_means_opaque(self):
        t = Transparency().with_overrides({"binder": 0.0}, roles=["binder"])
        assert t.is_opaque("binder")
        assert t.opaque_roles == ["binder"]

    def test_one_means_invisible_and_is_allowed(self):
        t = Transparency().with_overrides({"copy": 1.0}, roles=["copy"])
        assert t.for_role("copy") == 1.0

    @pytest.mark.parametrize("bad", [-0.1, 1.5])
    def test_out_of_range_is_an_error(self, bad):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            Transparency().with_overrides({"upstream": bad}, roles=["upstream"])

    def test_a_misspelt_role_is_an_error(self):
        with pytest.raises(ValueError, match="unknown role"):
            Transparency().with_overrides({"upsteam": 0.2}, roles=["upstream"])

    def test_a_bad_type_is_an_error(self):
        with pytest.raises(ValueError, match="number or a mapping"):
            Transparency().with_overrides("very see-through")

    def test_none_returns_the_same_object(self):
        t = Transparency()
        assert t.with_overrides(None) is t

    def test_overrides_do_not_mutate_the_default(self):
        Transparency().with_overrides({"upstream": 0.1}, roles=["upstream"])
        assert Transparency().for_role("upstream") == pytest.approx(0.65)


class TestRenderStyleFromDict:
    def test_empty_gives_the_defaults(self):
        s = RenderStyle.from_dict(None)
        assert s.palette.rgb("upstream") == DEFAULT_STYLE.palette.rgb("upstream")
        assert s.transparency_for("tunnel") == pytest.approx(0.75)

    def test_a_full_override(self):
        s = RenderStyle.from_dict({
            "colors": {"upstream": "#4caf50", "downstream": {"color": "teal", "lighten": 0.3}},
            "transparency": {"default": 0.5, "binder": 0.0},
            "background": "E6F2F9",
            "stick_radius": 0.6,
        })
        assert s.palette.base("upstream") == hex_to_rgb("#4caf50")
        assert s.transparency_for("downstream") == pytest.approx(0.5)
        assert s.transparency_for("binder") == 0.0
        assert s.background == "E6F2F9"
        assert s.stick_radius == pytest.approx(0.6)

    def test_making_a_role_opaque_flags_the_occlusion_step(self):
        s = RenderStyle.from_dict({"transparency": {"binder": 0.0}})
        assert s.needs_occlusion_clipping is True
        assert "occlusion-clipped" in s.summary()

    def test_british_and_american_palette_keys_both_work(self):
        for key in ("colors", "colours", "palette"):
            s = RenderStyle.from_dict({key: {"upstream": "red"}})
            assert s.palette.base("upstream") == (1.0, 0.0, 0.0)

    def test_two_palette_keys_at_once_is_an_error(self):
        with pytest.raises(ValueError, match="only one of"):
            RenderStyle.from_dict({"colors": {}, "colours": {}})

    def test_an_unknown_key_is_an_error(self):
        with pytest.raises(ValueError, match="unknown key"):
            RenderStyle.from_dict({"colour_scheme": "dark"})

    def test_a_bad_background_is_an_error(self):
        with pytest.raises(ValueError):
            RenderStyle.from_dict({"background": "nosuchcolour"})

    def test_each_agreed_background_is_accepted(self):
        for bg in BACKGROUNDS:
            assert RenderStyle.from_dict({"background": bg}).background == bg

    def test_a_nonpositive_stick_radius_is_an_error(self):
        with pytest.raises(ValueError, match="positive"):
            RenderStyle.from_dict({"stick_radius": 0})


class TestSettings:
    def test_ray_mode_draws_the_black_outline(self):
        assert RAY_SETTINGS["ray_trace_mode"] == 1
        assert RAY_SETTINGS["ray_trace_color"] == "black"

    def test_flat_mode_has_no_outline_and_is_fully_ambient(self):
        # ray_trace_mode 1 would outline the bumpy tunnel surface as creases.
        assert FLAT_SETTINGS["ray_trace_mode"] == 0
        assert FLAT_SETTINGS["ambient"] == 1.0
        assert FLAT_SETTINGS["direct"] == 0.0

    def test_neither_mode_casts_shadows_or_fog(self):
        for s in (RAY_SETTINGS, FLAT_SETTINGS):
            assert s["ray_shadows"] == 0 and s["ray_trace_fog"] == 0 and s["depth_cue"] == 0

    def test_the_two_modes_differ_only_in_the_lighting_keys(self):
        shared = set(RAY_SETTINGS) & set(FLAT_SETTINGS)
        differing = {k for k in shared if RAY_SETTINGS[k] != FLAT_SETTINGS[k]}
        assert differing == {"ray_trace_mode", "ambient", "direct", "reflect"}


class TestColorSpec:
    def test_validates_its_colour_on_construction(self):
        with pytest.raises(ValueError):
            ColorSpec("nosuchcolour")

    def test_validates_its_lighten_on_construction(self):
        with pytest.raises(ValueError):
            ColorSpec("red", 1.5)

    def test_final_applies_the_lightening(self):
        assert ColorSpec("green", 0.55).final == pytest.approx((0.55, 1.0, 0.55))

    def test_accepts_a_name_a_hex_or_numbers(self):
        assert ColorSpec("red").rgb == ColorSpec("#ff0000").rgb == ColorSpec([1, 0, 0]).rgb


class TestTheChannelIsMagentaAndThePatchIsRed:
    """Everything to do with the channel is magenta; the epitope is red.

    One colour was doing three jobs: the channel, the D_IC/O_IC marks and
    the hotspot patch all read `annotation`. They are three roles now, which
    is what lets the channel and its measurement share a colour deliberately
    while the patch takes a different one. User's defaults, 2026-10-27;
    every role is still overridable.
    """

    def test_the_channel_is_magenta(self):
        assert DEFAULT_STYLE.palette.rgb("tunnel") == (1.0, 0.0, 1.0)

    def test_the_measurement_is_the_same_magenta(self):
        assert (DEFAULT_STYLE.palette.rgb("annotation")
                == DEFAULT_STYLE.palette.rgb("tunnel"))

    def test_it_is_still_its_own_role(self):
        """Sharing a value is not sharing a role: either can move alone."""
        p = Palette().with_overrides({"annotation": (0.0, 0.0, 0.0)})
        assert p.rgb("annotation") == (0.0, 0.0, 0.0)
        assert p.rgb("tunnel") == (1.0, 0.0, 1.0)

    def test_the_pose_patch_is_red_and_has_its_own_role(self):
        assert DEFAULT_STYLE.palette.rgb("patch") == (0.9, 0.05, 0.05)
        assert "patch" in DEFAULT_STYLE.palette.roles

    def test_the_patch_is_not_the_channel_colour(self):
        """The epitope and the channel must not be mistaken for each other."""
        assert (DEFAULT_STYLE.palette.rgb("patch")
                != DEFAULT_STYLE.palette.rgb("tunnel"))

    def test_the_cofactor_is_not_the_channel_colour_either(self):
        """It sits inside the channel, so it cannot share its hue."""
        assert (DEFAULT_STYLE.palette.rgb("ligand")
                != DEFAULT_STYLE.palette.rgb("tunnel"))

    def test_the_defaults_do_not_warn_about_themselves(self):
        """A default configuration that reports its own colours as
        indistinguishable trains people to ignore the warning."""
        assert DEFAULT_STYLE.palette.indistinguishable() == []

    def test_the_deliberate_pair_is_skipped_and_nothing_else_is(self):
        from channelbridge.render.style import SAME_BY_DESIGN

        assert frozenset(("tunnel", "annotation")) in SAME_BY_DESIGN
        p = Palette().with_overrides({"downstream": {"color": (0.0, 1.0, 0.0),
                                                     "lighten": 0.55}})
        assert "downstream/upstream" in p.indistinguishable()

    def test_every_role_is_still_overridable(self):
        p = Palette().with_overrides({"annotation": (0.0, 0.0, 0.0),
                                      "patch": (1.0, 1.0, 0.0)})
        assert p.rgb("annotation") == (0.0, 0.0, 0.0)
        assert p.rgb("patch") == (1.0, 1.0, 0.0)
