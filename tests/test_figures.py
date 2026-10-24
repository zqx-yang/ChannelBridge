"""Finding the poses in an output folder, and planning the work for them."""

from __future__ import annotations

import os

import pytest

from channelbridge.figures import (
    GeomFiles,
    Runner,
    build,
    collect_poses,
    plan_scene_names,
    scan_output,
)
from channelbridge.render.select import resolve

STEMS = [
    ("geom001", "geom001_DO0.850_d017.0_o165_GDHp4-FNRp1_az000_ti090_sp240_ok"),
    ("geom002", "geom002_DO0.846_d013.0_o135_GDHp4-FNRp1_az330_ti090_sp240_ok"),
    ("geom003", "geom003_DO0.842_d018.0_o165_GDHp2-FNRp3_az330_ti060_sp270_ok"),
]


@pytest.fixture()
def out_dir(tmp_path):
    for _pid, stem in STEMS:
        (tmp_path / (stem + "_complex.pdb")).write_text("END\n")
        (tmp_path / (stem + "_occupancy.pdb")).write_text("END\n")
        (tmp_path / (stem + ".pml")).write_text("# pml\n")
    (tmp_path / (STEMS[0][1] + ".pse")).write_bytes(b"\x00")
    return str(tmp_path)


class TestScanOutput:
    def test_it_finds_every_pose(self, out_dir):
        assert [g.pose_id for g in scan_output(out_dir)] == ["geom001", "geom002", "geom003"]

    def test_it_reads_the_numbers_out_of_the_file_name(self, out_dir):
        g = scan_output(out_dir)[0]
        assert g.do == pytest.approx(0.850)
        assert g.d_ic == pytest.approx(17.0)
        assert g.o_ic == pytest.approx(165.0)

    def test_it_reads_which_patch_each_enzyme_contributes(self, out_dir):
        g = scan_output(out_dir)[2]
        assert (g.patch_up, g.patch_dn) == ("2", "3")

    def test_it_pairs_the_occupancy_and_session_files(self, out_dir):
        g = scan_output(out_dir)[0]
        assert g.occupancy_pdb and g.occupancy_pdb.endswith("_occupancy.pdb")
        assert g.pse and g.pse.endswith(".pse")

    def test_a_pose_without_a_complex_is_not_a_pose(self, tmp_path):
        (tmp_path / "geom009_DO0.5_whatever.pml").write_text("x")
        assert scan_output(str(tmp_path)) == []

    def test_a_missing_folder_says_so(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            scan_output(str(tmp_path / "nope"))

    def test_the_values_dict_only_holds_what_was_found(self):
        g = GeomFiles("geom001", "stem", ".", do=0.8)
        assert g.values == {"DO": 0.8}

    def test_poses_come_out_ready_for_the_selection_language(self, out_dir):
        poses = collect_poses(scan_output(out_dir))
        assert poses[0] == ("geom001", pytest.approx(0.850))


class TestPlanSceneNames:
    def test_monomer_only_by_default(self):
        plan = resolve([("geom001", 0.8)], svg="all")
        assert plan_scene_names(plan) == ["render_monomer"]

    def test_both_gives_both(self):
        plan = resolve([("geom001", 0.8)], svg="all", render="both")
        assert plan_scene_names(plan) == ["render_monomer", "render_multimer"]

    def test_the_cross_section_is_added_when_asked_for(self):
        plan = resolve([("geom001", 0.8)], svg="all", cross_section="monomer")
        assert "render_cs" in plan_scene_names(plan)

    def test_a_render_of_none_still_leaves_something_to_assemble(self):
        # asking for an SVG and no variant is a contradiction; assembling the
        # monomer is the least surprising reading of it
        plan = resolve([("geom001", 0.8)], svg="all", render="none")
        assert plan_scene_names(plan) == ["render_monomer"]


class TestBuild:
    """Dry runs: what would be done, without a PyMOL on the machine."""

    def test_a_pdb_selection_copies_the_complex(self, out_dir, tmp_path):
        geoms = scan_output(out_dir)
        plan = resolve(collect_poses(geoms), pdb="top2")
        made = build(plan, geoms, str(tmp_path / "fig"),
                     runner=Runner(dry_run=True, report=lambda *_: None),
                     report=lambda *_: None)
        assert len(made["pdb"]) == 2
        assert all(p.endswith("_complex.pdb") for p in made["pdb"])

    def test_an_svg_selection_plans_a_session_for_a_pose_that_has_none(self, out_dir, tmp_path):
        geoms = scan_output(out_dir)
        plan = resolve(collect_poses(geoms), svg="top2",
                       existing_pse=[g.pose_id for g in geoms if g.pse])
        assert "geom002" in plan.pse, "the second pose has no session and needs one"
        assert "geom001" in plan.pse_reused, "the first already has one"

    def test_a_pose_not_in_the_folder_is_reported_not_crashed_on(self, out_dir, tmp_path):
        geoms = scan_output(out_dir)
        plan = resolve(collect_poses(geoms), pdb="all")
        plan.pdb.append("geom404")
        made = build(plan, geoms, str(tmp_path / "fig"),
                     runner=Runner(dry_run=True, report=lambda *_: None),
                     report=lambda *_: None)
        assert any("geom404" in f for f in made["failed"])
        assert len(made["pdb"]) == 3, "the real poses still got done"

    def test_a_dry_run_writes_nothing(self, out_dir, tmp_path):
        geoms = scan_output(out_dir)
        plan = resolve(collect_poses(geoms), pdb="top1")
        target = tmp_path / "fig"
        build(plan, geoms, str(target), runner=Runner(dry_run=True, report=lambda *_: None),
              report=lambda *_: None)
        assert not target.exists()

    def test_a_real_copy_lands_where_it_says(self, out_dir, tmp_path):
        geoms = scan_output(out_dir)
        plan = resolve(collect_poses(geoms), pdb="rank1")
        target = str(tmp_path / "fig")
        made = build(plan, geoms, target, runner=Runner(dry_run=False, report=lambda *_: None),
                     report=lambda *_: None)
        assert os.path.isfile(made["pdb"][0])
        assert os.path.dirname(made["pdb"][0]) == target


class TestRunner:
    def test_it_prints_the_command_it_would_run(self, tmp_path):
        seen = []
        r = Runner(pymol="/nope/pymol", dry_run=True, report=seen.append)
        assert r.run("script.py", ["a=1"]) == 0
        assert any("script.py" in s and "a=1" in s for s in seen)

    def test_a_missing_pymol_is_reported_rather_than_raised(self, tmp_path):
        seen = []
        r = Runner(pymol=str(tmp_path / "nope"), dry_run=False, report=seen.append)
        assert r.run("script.py", []) == 127
        assert any("not found" in s for s in seen)

    def test_an_unset_pymol_flag_falls_back_to_the_default(self):
        seen = []
        Runner(pymol=None, dry_run=True, report=seen.append).run("s.py", [])
        assert any("pymol" in s for s in seen), "a None binary must not reach the command"

    def test_it_finds_the_render_scripts_in_the_package(self):
        assert os.path.isfile(Runner().script("run_pse.py"))
        assert os.path.isfile(Runner().script("run_png.py"))
        assert os.path.isfile(Runner().script("run_hotspots.py"))


class TestSvgRasterScene:
    """An SVG is built on an annotation-free raster, not on the render scene.

    Building it on the render scene prints D_IC/O_IC/DO twice -- once baked
    into the raster by PyMOL and once as the SVG's own editable text -- and
    leaves the 3D channel's solvent-side smear under the traced outline.
    """

    def test_each_render_variant_maps_to_its_base(self):
        from channelbridge.figures import _svg_raster_scene

        assert _svg_raster_scene("render_monomer") == "base_monomer"
        assert _svg_raster_scene("render_multimer") == "base_multimer"
        assert _svg_raster_scene("render_cs") == "base_cs"

    def test_the_session_builder_knows_those_scene_names(self):
        from channelbridge.figures import _svg_raster_scene, plan_scene_names
        from channelbridge.render.scenes import SCENE_NAMES

        plan = resolve([("geom001", 0.8)], svg="all", render="both",
                       cross_section="monomer")
        for scene in plan_scene_names(plan):
            assert _svg_raster_scene(scene) in SCENE_NAMES, scene


class TestPymolResolution:
    """Where PyMOL is found. A hard-coded cluster path made the tool
    unusable anywhere else, including in its own dev container."""

    def test_the_environment_variable_wins(self, monkeypatch):
        from channelbridge.figures import default_pymol

        monkeypatch.setenv("CHBRIDGE_PYMOL", "/somewhere/else/pymol")
        assert default_pymol() == "/somewhere/else/pymol"

    def test_otherwise_whatever_is_on_path(self, monkeypatch, tmp_path):
        import os
        import stat

        from channelbridge.figures import default_pymol

        fake = tmp_path / "pymol"
        fake.write_text("#!/bin/sh\n")
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
        monkeypatch.delenv("CHBRIDGE_PYMOL", raising=False)
        monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
        assert default_pymol() == str(fake)

    def test_the_runner_asks_at_call_time_not_import_time(self, monkeypatch):
        # a container sets CHBRIDGE_PYMOL after the module is imported
        from channelbridge.figures import Runner

        monkeypatch.setenv("CHBRIDGE_PYMOL", "/late/pymol")
        assert Runner().binary == "/late/pymol"
