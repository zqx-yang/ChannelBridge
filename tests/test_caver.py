"""Stage 2: CAVER configuration and which tunnels to carry forward.

The selection rule is not ours. It is Vavra et al. (2024), measured on
~13,158 experimental structures and ~29,752 CaverDock energy profiles:
the top-ranked tunnel is the energetically most favourable in 75 % of cases,
the second in a further 19 %. We implement it rather than re-deriving a
worse version of it from three domains.
"""

from __future__ import annotations

import pytest

from channelbridge.caver import (
    KEEP_TOP,
    PRIORITY_FLOOR,
    CaverParams,
    Cluster,
    caverdock_note,
    cavity_start_point,
    parse_residue_id,
    read_log,
    read_summary,
    select,
    write_config,
)


def cl(i, priority, length=20.0, br=2.0):
    return Cluster(cluster_id=i, n_tunnels=1, avg_bottleneck=br,
                   max_bottleneck=br, avg_length=length, priority=priority)


class TestDefaults:
    """Ours, not CAVER's, and the difference is deliberate."""

    def test_the_probe_is_a_water_probe(self):
        # CAVER defaults to 0.9, which threads gaps no water could pass. On
        # our six pilot domains that gave 51 clusters against 33 at 1.4,
        # with cluster 1 identical in span, bottleneck and exit.
        assert CaverParams().probe_radius == 1.4

    def test_the_seed_is_fixed(self):
        # CAVER's search is stochastic; without this the cluster NUMBERING
        # can change between runs, and every downstream reference is by
        # number
        assert CaverParams().seed == 1

    def test_waters_are_excluded(self):
        assert "HOH" in CaverParams().exclude

    def test_every_value_reaches_the_config(self):
        cfg = CaverParams(probe_radius=0.9, shell_radius=3.0).as_config()
        assert cfg["probe_radius"] == "0.90"
        assert cfg["shell_radius"] == "3.0"


class TestWriteConfig:
    def test_the_start_point_is_written_and_echoed(self, tmp_path):
        # CAVER silently MOVES a start point that does not fit its probe,
        # so the one we asked for has to be recoverable from the file
        p = write_config(str(tmp_path / "c.txt"), start_point=(1.0, 2.0, 3.0))
        text = open(p, encoding="utf-8").read()
        assert "starting_point_coordinates 1.000 2.000 3.000" in text
        assert "# starting point as given: 1.000 2.000 3.000" in text

    def test_a_bad_start_point_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="three numbers"):
            write_config(str(tmp_path / "c.txt"), start_point=(1.0, 2.0))

    def test_extras_are_passed_through(self, tmp_path):
        p = write_config(str(tmp_path / "c.txt"), start_point=(1.0, 2.0, 3.0),
                         extra={"max_tunnels": "50"})
        assert "max_tunnels 50" in open(p, encoding="utf-8").read()

    def test_no_starting_point_at_all_is_refused(self, tmp_path):
        # CAVER has no default and refuses to run; failing here names the
        # cause, failing in CAVER does not.
        with pytest.raises(ValueError, match="no starting point"):
            write_config(str(tmp_path / "c.txt"))

    def test_start_residues_are_written(self, tmp_path):
        params = CaverParams(start_residues=("B:128", "B:303", "20A"))
        text = open(write_config(str(tmp_path / "c.txt"), params),
                    encoding="utf-8").read()
        assert "starting_point_residue B:128 B:303 20A" in text

    def test_start_atoms_are_written(self, tmp_path):
        params = CaverParams(start_atoms=(578, 745))
        text = open(write_config(str(tmp_path / "c.txt"), params),
                    encoding="utf-8").read()
        assert "starting_point_atom 578 745" in text

    def test_the_three_start_routes_combine(self, tmp_path):
        # CAVER takes the centre of gravity of everything named, weighting
        # each entity equally -- they do not override one another, so a
        # config may legitimately carry all three.
        params = CaverParams(start_residues=("B:128",), start_atoms=(578,))
        text = open(write_config(str(tmp_path / "c.txt"), params,
                                 start_point=(1.0, 2.0, 3.0)),
                    encoding="utf-8").read()
        assert "starting_point_coordinates 1.000 2.000 3.000" in text
        assert "starting_point_residue B:128" in text
        assert "starting_point_atom 578" in text


class TestReadSummary:
    HEADER = ("  ID      No   No_snaps   Avg_BR       SD   Max_BR    Avg_L"
              "      SD   Avg_C      SD    Priority")

    def _summary(self, tmp_path, rows):
        p = tmp_path / "summary.txt"
        p.write_text("preamble\n" + self.HEADER + "\n" + "\n".join(rows) + "\n",
                     encoding="utf-8")
        return str(p)

    def test_it_reads_a_table(self, tmp_path):
        path = self._summary(tmp_path, [
            "  1  10  1  2.910  0.1  2.9  24.206  0.1  1.2  0.1  0.86052",
            "  2   8  1  2.101  0.1  2.1  21.010  0.1  1.1  0.1  0.81603",
        ])
        got = read_summary(path)
        assert [c.cluster_id for c in got] == [1, 2]
        assert got[0].priority == pytest.approx(0.86052)
        assert got[0].avg_length == pytest.approx(24.206)

    def test_columns_are_found_by_NAME_not_position(self, tmp_path):
        # a CAVER version that adds a column must not shift Priority into
        # another field and be believed
        header = self.HEADER.replace("Avg_BR", "NewCol       SD   Avg_BR")
        p = tmp_path / "summary.txt"
        p.write_text(header + "\n  1  10  1  9.9  0.0  2.910  0.1  2.9  "
                     "24.206  0.1  1.2  0.1  0.86052\n", encoding="utf-8")
        got = read_summary(str(p))
        assert got[0].priority == pytest.approx(0.86052)

    def test_a_directory_is_accepted(self, tmp_path):
        self._summary(tmp_path, ["  1  1  1  2.0  0  2.0  20.0  0  1  0  0.9"])
        assert read_summary(str(tmp_path))[0].cluster_id == 1

    def test_a_file_that_is_not_a_summary_is_refused(self, tmp_path):
        p = tmp_path / "summary.txt"
        p.write_text("this is not a CAVER summary\n", encoding="utf-8")
        with pytest.raises(ValueError, match="no cluster table"):
            read_summary(str(p))

    def test_they_come_back_best_first(self, tmp_path):
        path = self._summary(tmp_path, [
            "  1  1  1  2.0  0  2.0  20.0  0  1  0  0.40",
            "  2  1  1  2.0  0  2.0  20.0  0  1  0  0.90",
        ])
        assert [c.cluster_id for c in read_summary(path)] == [2, 1]


class TestSelect:
    def test_it_keeps_the_top_three_above_the_floor(self):
        got = select([cl(i, 0.9 - 0.05 * i) for i in range(1, 8)])
        assert len(got.kept) == KEEP_TOP

    def test_the_floor_is_the_published_one(self):
        assert PRIORITY_FLOOR == 0.55

    def test_a_low_priority_cluster_is_dropped_even_in_the_top_three(self):
        got = select([cl(1, 0.90), cl(2, 0.40), cl(3, 0.30)])
        assert [c.cluster_id for c in got.kept] == [1]

    def test_everything_below_the_floor_is_called_out(self):
        got = select([cl(1, 0.40), cl(2, 0.30)])
        assert not got.kept
        assert any("below priority" in n for n in got.notes)

    def test_no_clusters_at_all_is_not_a_crash(self):
        got = select([])
        assert got.best is None
        assert any("no tunnel clusters" in n for n in got.notes)

    def test_a_short_top_tunnel_is_flagged(self):
        # "If the first tunnel is shorter than 5 A, the pocket may be on the
        # surface and tunnel analysis may not be relevant" -- which is
        # exactly what the worked example's FNR cluster 1 is, at 3.2 A
        got = select([cl(1, 0.96, length=3.2), cl(2, 0.93, length=9.3)])
        assert any("may be on the surface" in n for n in got.notes)

    def test_a_near_tie_at_the_top_is_flagged(self):
        # cluster 2's priority was 91.6-99.9 % of cluster 1's across our
        # pilot domains; a seed change could swap them, and the swap moves
        # D_IC by tens of angstroms
        got = select([cl(1, 0.861), cl(2, 0.816)])
        assert any("near-tie" in n for n in got.notes)

    def test_a_clear_winner_is_not_flagged_as_a_tie(self):
        got = select([cl(1, 0.95), cl(2, 0.60)])
        assert not any("near-tie" in n for n in got.notes)

    def test_best_is_the_highest_priority_not_the_lowest_id(self):
        got = select([cl(7, 0.95), cl(1, 0.60)])
        assert got.best.cluster_id == 7


class TestCaverdockNote:
    """The note exists so that 'just run CaverDock' is not the default
    reading. Our own seed study: five identical runs gave E_a 2.00, 3.60,
    3.30, 3.30, 4.80 -- sd 1.00 against margins of 0.30-0.70."""

    def test_it_states_the_measured_noise(self):
        text = caverdock_note(select([cl(1, 0.9), cl(2, 0.8)]))
        assert "sd ~1.0" in text

    def test_it_says_caverdock_cannot_find_tunnels(self):
        text = caverdock_note(select([cl(1, 0.9)]))
        assert "cannot find one" in text

    def test_it_offers_the_published_alternative(self):
        assert "75 %" in caverdock_note(select([cl(1, 0.9)]))


class TestParseResidueId:
    @pytest.mark.parametrize("text,want", [
        ("157", (None, 157, "")),
        ("A:157", ("A", 157, "")),
        (" B:20 ", ("B", 20, "")),
        ("20A", (None, 20, "A")),
        ("B:20A", ("B", 20, "A")),
    ])
    def test_forms(self, text, want):
        assert parse_residue_id(text) == want

    @pytest.mark.parametrize("bad", ["", "   ", "abc", "A:", "A:xy"])
    def test_rejects_nonsense(self, bad):
        with pytest.raises(ValueError):
            parse_residue_id(bad)


class TestCavityStartPoint:
    """The derived start point, on a synthetic cavity of known geometry.

    A shell of atoms with a hollow centre: the answer is knowable, so this
    tests the rule rather than agreeing with whatever the code does.
    """

    @staticmethod
    def _shell(radius=8.0, n=14, hole=None):
        from channelbridge.io.pdb import Atom, AtomArray

        atoms = []
        g = [-1.0, 0.0, 1.0]
        serial = 0
        for i in g:
            for j in g:
                for k in g:
                    if i == j == k == 0:
                        continue
                    v = (i, j, k)
                    norm = sum(c * c for c in v) ** 0.5
                    xyz = tuple(radius * c / norm for c in v)
                    serial += 1
                    atoms.append(Atom(serial, "CA", "", "ALA", "A", serial, "",
                                      xyz, 1.0, 0.0, "C", False))
        # three residues around the hollow centre, pointing at each other:
        # their centroid is the origin, and each sits close to it
        for n_, xyz in enumerate([(2.0, 0, 0), (-1.0, 1.7, 0), (-1.0, -1.7, 0)], 1):
            serial += 1
            atoms.append(Atom(serial, "CB", "", "CYS", "A", 900 + n_, "",
                              xyz, 1.0, 0.0, "C", False))
        return AtomArray(atoms)

    def test_finds_a_point_off_the_blocked_centroid(self):
        arr = self._shell()
        got = cavity_start_point(arr, ["A:901", "A:902", "A:903"],
                                 probe_radius=1.0, reach=6.0, step=0.5)
        # the centroid itself is ~1.4 A from the CB atoms, too close for a
        # 1.0 A probe needing 2.6 A of clearance, so it must have moved
        assert got.offset > 0.0
        assert got.free >= 1.0 + 1.6

    def test_reports_the_evidence_it_is_a_pocket(self):
        got = cavity_start_point(self._shell(), ["A:901"], probe_radius=1.0,
                                 reach=6.0, step=0.5)
        assert got.burial > 0
        assert "burial" in got.summary()
        assert "free sphere" in got.summary()

    def test_prefers_the_nearest_qualifying_point_not_the_roomiest(self):
        # The roomiest point in a shell is the geometric centre; what we
        # want is the nearest one that fits, so the offset must be small.
        arr = self._shell(radius=12.0)
        got = cavity_start_point(arr, ["A:901", "A:902", "A:903"],
                                 probe_radius=1.0, reach=8.0, step=0.5)
        assert got.offset < 4.0

    def test_no_room_is_an_error_not_a_silent_bad_point(self):
        with pytest.raises(ValueError, match="can hold a"):
            cavity_start_point(self._shell(), ["A:901", "A:902", "A:903"],
                               probe_radius=20.0, reach=4.0, step=0.5)

    def test_unknown_residue_is_an_error(self):
        with pytest.raises(ValueError, match="none of the residues"):
            cavity_start_point(self._shell(), ["A:4242"], probe_radius=1.0)

    def test_no_residues_is_an_error(self):
        with pytest.raises(ValueError, match="no residues given"):
            cavity_start_point(self._shell(), [])


class TestReadLog:
    """What CAVER loaded, which is not what the file holds.

    The gap is silent in both directions and was measured the hard way: a
    cargo not named in `exclude_residue_names` can still be missing from the
    measurement (a clashing AlphaFill transplant is dropped), and a residue
    holding the starting point is removed without a word. A pipeline that
    wants to claim the cargo was out of the channel has to read this count.
    """

    LOG = """Caver computation started.
Mapping snapshots...
...finished.

Basic parameters:
starting_point_coordinates [4.811, 9.895, -12.308]
probe_radius 1.4

* Processing 1 *
Residues: 597
Atoms: 4911
Constructing additively weighted voronoi diagram for 4911 spheres.
Voronoi vertices: 33299
"""

    def test_it_reads_both_counts(self, tmp_path):
        p = tmp_path / "log.txt"
        p.write_text(self.LOG, encoding="utf-8")
        got = read_log(str(p))
        assert (got.residues, got.atoms) == (597, 4911)

    def test_a_run_directory_works_too(self, tmp_path):
        (tmp_path / "log.txt").write_text(self.LOG, encoding="utf-8")
        assert read_log(str(tmp_path)).atoms == 4911

    def test_the_first_snapshot_wins(self, tmp_path):
        p = tmp_path / "log.txt"
        p.write_text(self.LOG + "\n* Processing 2 *\nResidues: 1\nAtoms: 2\n",
                     encoding="utf-8")
        assert read_log(str(p)).atoms == 4911

    def test_a_number_elsewhere_in_the_line_is_not_taken(self, tmp_path):
        p = tmp_path / "log.txt"
        p.write_text("Atoms: 7 of them\nResidues: 3\nAtoms: 11\n",
                     encoding="utf-8")
        # "Atoms: 7 of them" is four fields, not two: not a count line.
        assert read_log(str(p)).atoms == 11

    def test_no_count_is_an_error_not_a_zero(self, tmp_path):
        """A check that cannot run must not report a pass."""
        p = tmp_path / "log.txt"
        p.write_text("Caver computation started.\n", encoding="utf-8")
        with pytest.raises(ValueError, match="Atoms:"):
            read_log(str(p))

    def test_a_missing_log_is_an_error(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_log(str(tmp_path / "nope"))
