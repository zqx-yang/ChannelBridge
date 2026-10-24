"""The ``chbridge`` command line, driven in-process."""

from __future__ import annotations

import json

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
