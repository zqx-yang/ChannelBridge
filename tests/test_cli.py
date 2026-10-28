"""The ``chbridge`` command line, driven in-process."""

from __future__ import annotations

import json
import os

import pytest

from channelbridge.cli import _unique_labels, main

CAVER_SPHERE = "ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f      %6.2f"


def make_caver_out(root, name, n=12, step=1.0, radius=2.5):
    """A minimal CAVER output folder: data/clusters/tun_cl_001.pdb."""
    d = root / name / "data" / "clusters"
    d.mkdir(parents=True)
    lines = [
        CAVER_SPHERE % (i + 1, i * step, 0.0, 0.0, radius) for i in range(n)
    ]
    (d / "tun_cl_001.pdb").write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return root / name


def make_structure(path, n=12, spacing=1.5):
    """A cube of CA atoms around the origin, so burial is non-zero inside."""
    lines = []
    k = 0
    g = [(i - (n - 1) / 2.0) * spacing for i in range(n)]
    for x in g:
        for y in g:
            for z in g:
                k += 1
                lines.append(
                    "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C"
                    % (k, 1, x, y, z)
                )
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path


class TestEntryPoints:
    def test_no_command_prints_help_and_succeeds(self, capsys):
        assert main([]) == 0
        assert "chbridge" in capsys.readouterr().out

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as e:
            main(["--version"])
        assert e.value.code == 0
        assert "chbridge" in capsys.readouterr().out


class TestDo:
    def test_reproduces_the_recorded_geom001(self, capsys):
        assert main(["do", "--dic", "16", "--oic", "150"]) == 0
        out = capsys.readouterr().out
        assert "0.8401" in out
        assert "6.6 deg" in out

    def test_all_sets_shows_the_divergence(self, capsys):
        assert main(["do", "--dic", "16", "--oic", "150", "--all-sets"]) == 0
        out = capsys.readouterr().out
        # The paper's stated constants are the default; the other two sets are
        # shown with their offset so the choice is never silent.
        assert "imars-st2a" in out
        assert "imars-exact" in out
        assert "+0.024" in out          # st2a sits ~0.025 above
        assert "6.6 deg" in out

    def test_impossible_angle_fails_cleanly(self, capsys):
        assert main(["do", "--dic", "16", "--oic", "200"]) == 2
        assert "error:" in capsys.readouterr().err


class TestValidate:
    def test_reports_a_valid_pair_file(self, tmp_path, capsys):
        (tmp_path / "a.pdb").write_text("END\n", encoding="utf-8")
        (tmp_path / "b.pdb").write_text("END\n", encoding="utf-8")
        pair = {
            "name": "t",
            "channelled": "X",
            "upstream": 1,
            "targets": [
                {"id": 1, "name": "A", "file": "a.pdb", "chain": "B", "site": "s", "vent": "v"},
                {"id": 2, "name": "B", "file": "b.pdb", "chain": "N", "site": "s", "vent": "v"},
            ],
        }
        p = tmp_path / "pair.json"
        p.write_text(json.dumps(pair), encoding="utf-8")
        assert main(["validate", str(p)]) == 0
        assert "ok: pair file is valid" in capsys.readouterr().out

    def test_missing_target_file_is_a_nonzero_exit(self, tmp_path, capsys):
        pair = {
            "name": "t",
            "channelled": "X",
            "upstream": 1,
            "targets": [
                {"id": 1, "name": "A", "file": "a.pdb", "chain": "B", "site": "s", "vent": "v"},
                {"id": 2, "name": "B", "file": "b.pdb", "chain": "N", "site": "s", "vent": "v"},
            ],
        }
        p = tmp_path / "pair.json"
        p.write_text(json.dumps(pair), encoding="utf-8")
        assert main(["validate", str(p)]) == 1
        assert "PROBLEM" in capsys.readouterr().err

    def test_malformed_pair_file_exits_two(self, tmp_path, capsys):
        p = tmp_path / "pair.json"
        p.write_text("{", encoding="utf-8")
        assert main(["validate", str(p)]) == 2
        assert "error:" in capsys.readouterr().err


class TestTunnel:
    def test_raw_only_without_a_structure(self, tmp_path, capsys):
        out = make_caver_out(tmp_path, "out_X")
        assert main(["tunnel", str(out)]) == 0
        cap = capsys.readouterr()
        assert "raw" in cap.out
        assert "without --structure" in cap.err

    def test_normalises_when_given_a_structure(self, tmp_path, capsys):
        out = make_caver_out(tmp_path, "out_X", n=40, step=0.5)
        s = make_structure(tmp_path / "s.pdb")
        assert main(["tunnel", str(out), "--structure", str(s)]) == 0
        o = capsys.readouterr().out
        assert "normalised" in o
        assert "resample" in o

    def test_compare_reports_exit_separation(self, tmp_path, capsys):
        a = make_caver_out(tmp_path / "A", "out_X", n=40, step=0.5)
        b = make_caver_out(tmp_path / "B", "out_X", n=12, step=1.0)
        s = make_structure(tmp_path / "s.pdb")
        assert main(["tunnel", str(a), str(b), "--structure", str(s), "--compare"]) == 0
        o = capsys.readouterr().out
        assert "largest exit separation" in o

    def test_identical_folder_names_are_disambiguated(self, tmp_path, capsys):
        a = make_caver_out(tmp_path / "A", "out_X")
        b = make_caver_out(tmp_path / "B", "out_X")
        assert main(["tunnel", str(a), str(b)]) == 0
        o = capsys.readouterr().out
        assert "A/out_X" in o and "B/out_X" in o

    def test_a_folder_with_no_clusters_explains_itself(self, tmp_path, capsys):
        d = tmp_path / "out_empty" / "data"
        d.mkdir(parents=True)
        assert main(["tunnel", str(tmp_path / "out_empty")]) == 2
        assert "no pathway" in capsys.readouterr().err

    def test_missing_structure_fails_cleanly(self, tmp_path, capsys):
        out = make_caver_out(tmp_path, "out_X")
        assert main(["tunnel", str(out), "--structure", str(tmp_path / "nope.pdb")]) == 2
        assert "structure not found" in capsys.readouterr().err

    def test_excluding_everything_is_caught(self, tmp_path, capsys):
        out = make_caver_out(tmp_path, "out_X")
        s = make_structure(tmp_path / "s.pdb")
        rc = main(["tunnel", str(out), "--structure", str(s), "--exclude", "ALA"])
        assert rc == 2
        assert "every atom" in capsys.readouterr().err

    def test_a_named_cluster_that_does_not_exist_lists_what_does(self, tmp_path, capsys):
        out = make_caver_out(tmp_path, "out_X")
        assert main(["tunnel", str(out), "--cluster", "7"]) == 2
        assert "this run has clusters" in capsys.readouterr().err

    def test_a_tun_cl_file_can_be_passed_directly(self, tmp_path, capsys):
        out = make_caver_out(tmp_path, "out_X")
        f = out / "data" / "clusters" / "tun_cl_001.pdb"
        assert main(["tunnel", str(f)]) == 0
        assert "raw" in capsys.readouterr().out


class TestUniqueLabels:
    def test_distinct_basenames_stay_short(self, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        assert _unique_labels([str(tmp_path / "a"), str(tmp_path / "b")]) == ["a", "b"]

    def test_colliding_basenames_grow_leftwards(self, tmp_path):
        (tmp_path / "A").mkdir()
        (tmp_path / "B").mkdir()
        (tmp_path / "A" / "out_X").mkdir()
        (tmp_path / "B" / "out_X").mkdir()
        labels = _unique_labels([str(tmp_path / "A" / "out_X"), str(tmp_path / "B" / "out_X")])
        assert labels == ["A/out_X", "B/out_X"]

    def test_single_path_is_its_basename(self, tmp_path):
        assert _unique_labels([str(tmp_path / "only")]) == ["only"]


class TestSessionAndRenderAreSeparateSteps:
    """``session`` builds the .pse, ``render`` draws from one.

    Kept apart on the user's instruction (2026-10-25) so that "open it,
    rotate it to the frame you want, save it, render THAT" is a first-class
    path -- and so ``render`` accepts a session the user built themselves.
    """

    def test_session_names_its_output_from_the_prefix(self, tmp_path, capsys):
        out = tmp_path / "out"
        out.mkdir()
        (out / "geom001_complex.pdb").write_text("END\n", encoding="utf-8")
        assert main(["session", str(out), "--prefix", "run1", "--dry-run"]) == 0
        printed = capsys.readouterr().out
        assert "run1_geom001_render.pse" in printed

    def test_render_refuses_a_session_that_is_not_there(self, tmp_path, capsys):
        rc = main(["render", str(tmp_path / "nope.pse"), "--dry-run"])
        assert rc != 0
        assert "no such session" in capsys.readouterr().err

    def test_render_takes_any_session_and_names_both_variants(self, tmp_path, capsys):
        pse = tmp_path / "hand_rotated.pse"
        pse.write_bytes(b"not really a session")
        assert main(["render", str(pse), "--prefix", "run1",
                     "--pose", "geom007", "--dry-run"]) == 0
        printed = capsys.readouterr().out
        assert "run1_geom007_render_pair.svg" in printed
        assert "run1_geom007_render_assembly.svg" in printed


class TestPseArgReachesTheSessionBuilder:
    """``--pse-arg`` is passed on verbatim, one token per occurrence.

    It was dropped in silence. `cmd_session` handed an ``action="append"``
    list to a splitter written for a single comma-separated string, whose
    ``if not value`` guard made the empty list falsy -- so every ``--pse-arg``
    vanished and the session builder failed with "no chains given for the
    upstream enzyme", naming neither the option nor the reason. Comma
    splitting is wrong too: these are ``key=value`` tokens whose values
    contain commas.
    """

    def _cmd(self, tmp_path, capsys, extra):
        out = tmp_path / "out"
        out.mkdir()
        (out / "geom001_complex.pdb").write_text("END\n", encoding="utf-8")
        assert main(["session", str(out), "--prefix", "p", "--dry-run",
                     *extra]) == 0
        return capsys.readouterr().out

    def test_one_arg_arrives(self, tmp_path, capsys):
        printed = self._cmd(tmp_path, capsys, ["--pse-arg", "up_chains=A"])
        assert "up_chains=A" in printed

    def test_a_value_may_contain_commas(self, tmp_path, capsys):
        printed = self._cmd(tmp_path, capsys, ["--pse-arg", "up_chains=B,C,D,E"])
        assert "up_chains=B,C,D,E" in printed

    def test_the_option_is_repeatable(self, tmp_path, capsys):
        printed = self._cmd(tmp_path, capsys,
                            ["--pse-arg", "up_chains=A",
                             "--pse-arg", "dn_chains=B",
                             "--pse-arg", "ligands=FPP,GST"])
        for want in ("up_chains=A", "dn_chains=B", "ligands=FPP,GST"):
            assert want in printed

    def test_case_is_not_changed(self, tmp_path, capsys):
        printed = self._cmd(tmp_path, capsys, ["--pse-arg", "label=enz_GPPS"])
        assert "label=enz_GPPS" in printed

    def test_none_given_is_still_fine(self, tmp_path, capsys):
        printed = self._cmd(tmp_path, capsys, [])
        assert "p_geom001_render.pse" in printed


class TestRenderLabelsComeFromTheSession:
    """`render` must not label 0.0 just because nobody told it the numbers.

    `--dic/--oic/--do` default to 0.0, which is a real-looking value, so a
    render driven from a session alone published "D_IC = 0.0 A, O_IC = 0.0
    deg, DO = 0.000" over an otherwise correct figure with nothing in the
    output saying the labels were unset. `chbridge session` is the only place
    those numbers exist, so it writes them beside the .pse and they are read
    back here.
    """

    def _session(self, tmp_path, **geom):
        pse = tmp_path / "s.pse"
        pse.write_bytes(b"not really a session")
        if geom:
            (tmp_path / "s_geometry.json").write_text(json.dumps(geom),
                                                      encoding="utf-8")
        return pse

    def test_the_sidecar_is_read(self, tmp_path, capsys):
        pse = self._session(tmp_path, d_ic=6.0, o_ic=165.0, do=0.9371)
        assert main(["render", str(pse), "--variant", "pair",
                     "--dry-run"]) == 0
        out = capsys.readouterr().out
        assert "D_IC 6.0 A" in out and "O_IC 165.0 deg" in out
        assert "DO 0.9371" in out

    def test_an_explicit_flag_wins(self, tmp_path, capsys):
        pse = self._session(tmp_path, d_ic=6.0, o_ic=165.0, do=0.9371)
        assert main(["render", str(pse), "--variant", "pair", "--dic", "12",
                     "--dry-run"]) == 0
        out = capsys.readouterr().out
        assert "D_IC 12.0 A" in out
        assert "O_IC 165.0 deg" in out        # the others still come from it

    def test_no_numbers_anywhere_warns_rather_than_labelling_zero(
            self, tmp_path, capsys):
        pse = self._session(tmp_path)
        assert main(["render", str(pse), "--variant", "pair",
                     "--dry-run"]) == 0
        err = capsys.readouterr().err
        assert "reads as a measurement" in err

    def test_a_corrupt_sidecar_is_reported_not_swallowed(self, tmp_path, capsys):
        pse = self._session(tmp_path)
        (tmp_path / "s_geometry.json").write_text("{not json", encoding="utf-8")
        assert main(["render", str(pse), "--variant", "pair",
                     "--dry-run"]) == 0
        assert "ignoring" in capsys.readouterr().err


class TestRenderHelpersAreNotCachedOnExistence:
    """A corrected session must not render with the old session's helpers.

    The flat helper renders -- mask, mouth markers, footprint, transparent
    base -- live in `intermediate/` and were reused whenever the file existed.
    Measured: a session rebuilt with its tunnels correctly carried onto the
    pose produced a figure with corrected D_IC/O_IC/DO labels over exit
    markers at byte-identical pixel coordinates from the old, wrong frame --
    a picture that disagrees with its own caption, and nothing said so.
    """

    def test_a_newer_session_reruns_the_helpers(self, tmp_path):
        from channelbridge.cli import _older_than

        inter = tmp_path / "intermediate"
        inter.mkdir()
        helper = inter / "p_geom001_mask.png"
        helper.write_bytes(b"old")
        pse = tmp_path / "s.pse"
        pse.write_bytes(b"session")
        os.utime(helper, (1, 1))
        os.utime(pse, (100, 100))
        assert _older_than(str(helper), str(pse)) is True

    def test_an_up_to_date_helper_is_kept(self, tmp_path):
        from channelbridge.cli import _older_than

        helper = tmp_path / "h.png"
        pse = tmp_path / "s.pse"
        helper.write_bytes(b"x")
        pse.write_bytes(b"y")
        os.utime(pse, (1, 1))
        os.utime(helper, (100, 100))
        assert _older_than(str(helper), str(pse)) is False

    def test_a_missing_helper_is_stale(self, tmp_path):
        from channelbridge.cli import _older_than

        assert _older_than(str(tmp_path / "nope.png")) is True
