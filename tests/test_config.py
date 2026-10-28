"""Pair-file loading and validation."""

from __future__ import annotations

import json

import pytest

from channelbridge.config import ConfigError, PairFile, PoseSpec, load_pair_file

MINIMAL = {
    "name": "gdh_fnr",
    "channelled": "NADPH",
    "upstream": 1,
    "targets": [
        {
            "id": 1,
            "name": "GDH",
            "file": "gdh.pdb",
            "chain": "B",
            "site": "resn NAP and name C4N",
            "vent": "resn NAP and name N1A",
        },
        {
            "id": 2,
            "name": "FNR",
            "file": "fnr.pdb",
            "chain": "N",
            "site": "resn FDA and name N5",
            "vent": "resn A2P",
        },
    ],
}


def write_json(tmp_path, data, name="pair.json"):
    p = tmp_path / name
    p.write_text(json.dumps(data), encoding="utf-8")
    return p


def load(tmp_path, data):
    return load_pair_file(write_json(tmp_path, data))


class TestMinimal:
    def test_loads(self, tmp_path):
        pair = load(tmp_path, MINIMAL)
        assert pair.name == "gdh_fnr"
        assert pair.channelled == "NADPH"
        assert pair.upstream_target.name == "GDH"
        assert pair.downstream_target.name == "FNR"

    def test_relative_paths_resolve_against_the_pair_file(self, tmp_path):
        pair = load(tmp_path, MINIMAL)
        assert pair.target(1).file == (tmp_path / "gdh.pdb").resolve()

    def test_defaults_are_the_normalising_ones(self, tmp_path):
        pair = load(tmp_path, MINIMAL)
        assert pair.tunnel.normalise is True
        assert pair.tunnel.spacing == 1.0
        assert pair.tunnel.burial_min == 15
        assert pair.tunnel.oic_arm == 1.0

    def test_default_do_set_is_the_published_one(self, tmp_path):
        assert load(tmp_path, MINIMAL).do_coefficients.name == "imars-fig2c"

    def test_summary_mentions_the_direction_and_the_normalisation(self, tmp_path):
        s = load(tmp_path, MINIMAL).summary()
        assert "upstream" in s and "GDH" in s
        assert "normalised" in s

    def test_summary_warns_when_normalisation_is_off(self, tmp_path):
        d = dict(MINIMAL, tunnel={"normalise": False})
        assert "RAW" in load(tmp_path, d).summary()


class TestRejections:
    def test_a_typo_is_an_error_not_a_silent_default(self, tmp_path):
        # channel_ligand misspelt would silently disable the channel.
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][0]["channel_ligands"] = "NAP"
        with pytest.raises(ConfigError, match="unknown key"):
            load(tmp_path, bad)

    def test_unknown_top_level_key(self, tmp_path):
        with pytest.raises(ConfigError, match="unknown key"):
            load(tmp_path, dict(MINIMAL, oops=1))

    def test_missing_required_key_names_it(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        del bad["channelled"]
        with pytest.raises(ConfigError, match="channelled"):
            load(tmp_path, bad)

    def test_needs_exactly_two_targets(self, tmp_path):
        with pytest.raises(ConfigError, match="exactly 2 targets"):
            load(tmp_path, dict(MINIMAL, targets=MINIMAL["targets"][:1]))

    def test_upstream_must_be_1_or_2(self, tmp_path):
        with pytest.raises(ConfigError, match="one-sided"):
            load(tmp_path, dict(MINIMAL, upstream=3))

    def test_duplicate_chain_is_rejected(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][1]["chain"] = "B"
        with pytest.raises(ConfigError, match="chain"):
            load(tmp_path, bad)

    def test_duplicate_name_is_rejected(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][1]["name"] = "GDH"
        with pytest.raises(ConfigError, match="same name"):
            load(tmp_path, bad)

    def test_multi_character_chain_is_rejected(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][0]["chain"] = "BC"
        with pytest.raises(ConfigError, match="single PDB chain"):
            load(tmp_path, bad)

    def test_target_name_must_be_filename_safe(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][0]["name"] = "GDH/1"
        with pytest.raises(ConfigError, match="file names"):
            load(tmp_path, bad)

    def test_target_ids_must_be_one_and_two(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][0]["id"] = 3
        with pytest.raises(ConfigError, match="id must be 1 or 2"):
            load(tmp_path, bad)

    def test_unknown_do_set(self, tmp_path):
        with pytest.raises(ConfigError, match="unknown set"):
            load(tmp_path, dict(MINIMAL, do_coefficients="made-up"))

    def test_bad_oic_angle(self, tmp_path):
        with pytest.raises(ConfigError, match=r"\[0, 180\]"):
            load(tmp_path, dict(MINIMAL, pose={"oic": [200]}))

    def test_inverted_dic_range(self, tmp_path):
        with pytest.raises(ConfigError, match="above"):
            load(tmp_path, dict(MINIMAL, pose={"dic_min": 30, "dic_max": 10}))

    def test_nonpositive_spacing(self, tmp_path):
        with pytest.raises(ConfigError, match="spacing must be positive"):
            load(tmp_path, dict(MINIMAL, tunnel={"spacing": 0}))

    def test_cluster_is_one_based(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][0]["caver"] = {"cluster": 0}
        with pytest.raises(ConfigError, match="1-based"):
            load(tmp_path, bad)

    def test_bad_start_point(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][0]["caver"] = {"start": [1.0, 2.0]}
        with pytest.raises(ConfigError, match=r"\[x, y, z\]"):
            load(tmp_path, bad)

    def test_missing_file(self, tmp_path):
        with pytest.raises(ConfigError, match="not found"):
            load_pair_file(tmp_path / "nope.json")

    def test_unsupported_extension(self, tmp_path):
        p = tmp_path / "pair.toml"
        p.write_text("x = 1", encoding="utf-8")
        with pytest.raises(ConfigError, match=r"expected a \.yaml"):
            load_pair_file(p)

    def test_empty_file(self, tmp_path):
        p = tmp_path / "pair.json"
        p.write_text("null", encoding="utf-8")
        with pytest.raises(ConfigError, match="empty"):
            load_pair_file(p)

    def test_invalid_json_says_so(self, tmp_path):
        p = tmp_path / "pair.json"
        p.write_text("{not json", encoding="utf-8")
        with pytest.raises(ConfigError, match="not valid JSON"):
            load_pair_file(p)


class TestDiskValidation:
    def test_reports_a_missing_target_file(self, tmp_path):
        pair = load(tmp_path, MINIMAL)
        problems = pair.validate_against_disk()
        assert len(problems) == 2
        assert "file not found" in problems[0]

    def test_rejects_mmcif_because_caver_cannot_read_it(self, tmp_path):
        bad = json.loads(json.dumps(MINIMAL))
        bad["targets"][0]["file"] = "gdh.cif"
        bad["targets"][1]["file"] = "fnr.pdb"
        (tmp_path / "gdh.cif").write_text("data_x", encoding="utf-8")
        (tmp_path / "fnr.pdb").write_text("END", encoding="utf-8")
        problems = load(tmp_path, bad).validate_against_disk()
        assert any("mmCIF" in p for p in problems)

    def test_clean_when_files_exist(self, tmp_path):
        (tmp_path / "gdh.pdb").write_text("END", encoding="utf-8")
        (tmp_path / "fnr.pdb").write_text("END", encoding="utf-8")
        assert load(tmp_path, MINIMAL).validate_against_disk() == []


class TestPoseGrid:
    def test_default_grid_size_matches_the_documented_scan(self):
        # D_IC 6-36 step 1 (31) x 12 spin x 4 tilt, and per O_IC: 12 azimuths
        # for each of the six below 180, but ONE at 180, where the cone the
        # azimuth sweeps is a single point. The naive product (31*7*12*12*4 =
        # 124,992) is what the search does not scan.
        assert PoseSpec().n_poses == 31 * (6 * 12 + 1) * 12 * 4
        assert PoseSpec().n_poses == 108_624

    def test_the_estimate_equals_what_the_search_actually_scans(self):
        """Two numbers for one grid is a bug report waiting to happen."""
        from channelbridge.pose import PoseParams, pose_grid

        for raw in ({}, {"oic": [180]}, {"oic": [150]}, {"oic": [180, 90]},
                    {"dic_min": 8, "dic_max": 20, "dic_step": 4},
                    {"tilt": [0, 45], "spin": [0, 90], "azimuth": [0, 120, 240]}):
            spec = PoseSpec.from_dict("pose", raw)
            assert spec.n_poses == len(pose_grid(PoseParams.from_spec(spec))), raw

    def test_a_smaller_grid_is_reported_smaller(self):
        small = PoseSpec.from_dict("pose", {"tilt": [0], "azimuth": [0], "spin": [0]})
        assert small.n_poses == 31 * 7


class TestYaml:
    def test_yaml_loads_when_pyyaml_is_available(self, tmp_path):
        pytest.importorskip("yaml")
        import yaml

        p = tmp_path / "pair.yaml"
        p.write_text(yaml.safe_dump(MINIMAL), encoding="utf-8")
        assert load_pair_file(p).name == "gdh_fnr"

    def test_from_dict_rejects_a_non_mapping(self, tmp_path):
        with pytest.raises(ConfigError, match="mapping"):
            PairFile.from_dict([1, 2], base=tmp_path)


class TestRenderSection:
    """Figure appearance travels in the pair file, with the defaults built in."""

    def test_absent_render_section_gives_the_defaults(self, tmp_path):
        from channelbridge.render.style import DEFAULT_STYLE

        pair = load(tmp_path, MINIMAL)
        assert pair.render.palette.rgb("upstream") == DEFAULT_STYLE.palette.rgb("upstream")
        assert pair.render.transparency_for("tunnel") == pytest.approx(0.75)
        assert pair.render.needs_occlusion_clipping is False

    def test_colours_and_transparency_are_overridable(self, tmp_path):
        d = dict(MINIMAL, render={
            "colors": {"upstream": "#4caf50", "downstream": {"color": "teal", "lighten": 0.3}},
            "transparency": {"default": 0.5, "binder": 0.0, "tunnel": 0.6},
            "background": "CFE6F3",
        })
        pair = load(tmp_path, d)
        assert pair.render.palette.base("upstream") == pytest.approx((0.298, 0.6863, 0.3137), abs=1e-3)
        assert pair.render.palette.spec("downstream").lighten == pytest.approx(0.3)
        assert pair.render.transparency_for("downstream") == pytest.approx(0.5)
        assert pair.render.transparency_for("tunnel") == pytest.approx(0.6)
        assert pair.render.transparency.is_opaque("binder")
        assert pair.render.background == "CFE6F3"

    def test_a_misspelt_colour_role_is_rejected(self, tmp_path):
        d = dict(MINIMAL, render={"colors": {"upsteam": "green"}})
        with pytest.raises(ConfigError, match="unknown colour role"):
            load(tmp_path, d)

    def test_a_misspelt_transparency_role_is_rejected(self, tmp_path):
        d = dict(MINIMAL, render={"transparency": {"upsteam": 0.2}})
        with pytest.raises(ConfigError, match="unknown role"):
            load(tmp_path, d)

    def test_an_unknown_render_key_is_rejected(self, tmp_path):
        with pytest.raises(ConfigError, match="unknown key"):
            load(tmp_path, dict(MINIMAL, render={"colour_scheme": "dark"}))

    def test_an_out_of_range_transparency_is_rejected(self, tmp_path):
        with pytest.raises(ConfigError, match=r"\[0, 1\]"):
            load(tmp_path, dict(MINIMAL, render={"transparency": {"tunnel": 1.5}}))

    def test_the_summary_says_whether_anything_is_opaque(self, tmp_path):
        assert "nothing opaque" in load(tmp_path, MINIMAL).summary()
        d = dict(MINIMAL, render={"transparency": {"binder": 0.0}})
        assert "occlusion clipping" in load(tmp_path, d).summary()
