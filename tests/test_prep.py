"""Stage 1: selecting the atoms the rest of the pipeline measures."""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.io.pdb import read_pdb
from channelbridge.prep import FURNITURE, PrepParams, prepare, write_pdb

ATOM = ("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s")


def structure(tmp_path, rows):
    lines = []
    for i, (rec, name, alt, resn, ch, seq, xyz, occ) in enumerate(rows, 1):
        lines.append(ATOM % (rec, i, " %-3s" % name, alt, resn, ch, seq, " ",
                             xyz[0], xyz[1], xyz[2], occ, 20.0, " C"))
    p = tmp_path / "s.pdb"
    p.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return str(p)


def protein(ch="A", seq=1, n=4):
    return [("ATOM  ", nm, "", "ALA", ch, seq + i, (float(i), 0.0, 0.0), 1.0)
            for i in range(n) for nm in ("CA",)]


class TestChains:
    def test_the_first_polymer_chain_by_default(self, tmp_path):
        s = structure(tmp_path, protein("A") + protein("B", 10))
        _atoms, rep = prepare(s)
        assert rep.chains_kept == ["A"] and rep.chains_dropped == ["B"]

    def test_and_it_says_so_when_there_was_a_choice(self, tmp_path):
        # the tunnel, the residue classes and the pose all refer to ONE
        # protomer, so a silent pick is a trap
        s = structure(tmp_path, protein("A") + protein("B", 10))
        _atoms, rep = prepare(s)
        assert any("ONE protomer" in n for n in rep.notes)

    def test_no_note_when_there_was_no_choice(self, tmp_path):
        _atoms, rep = prepare(structure(tmp_path, protein("A")))
        assert not rep.notes

    def test_an_absent_chain_is_an_error_not_an_empty_file(self, tmp_path):
        s = structure(tmp_path, protein("A"))
        with pytest.raises(ValueError, match="not in"):
            prepare(s, PrepParams(chains=["Z"]))

    def test_keeping_everything_leaves_nothing_is_refused(self, tmp_path):
        s = structure(tmp_path, [("HETATM", "O", "", "HOH", "A", 1,
                                  (0.0, 0.0, 0.0), 1.0)])
        with pytest.raises(ValueError, match="no polymer atoms"):
            prepare(s)


class TestLigands:
    def _with(self, tmp_path, resn):
        return structure(tmp_path, [*protein("A"),
                                    ("HETATM", "C1", "", resn, "A", 99,
                                     (9.0, 0, 0), 1.0)])

    def test_waters_go(self, tmp_path):
        _a, rep = prepare(self._with(tmp_path, "HOH"))
        assert rep.waters_dropped == 1

    def test_waters_can_be_kept_if_you_insist(self, tmp_path):
        _a, rep = prepare(self._with(tmp_path, "HOH"),
                          PrepParams(keep_waters=True))
        assert rep.waters_dropped == 0

    def test_crystallisation_furniture_goes(self, tmp_path):
        _a, rep = prepare(self._with(tmp_path, "GOL"))
        assert "GOL" in rep.ligands_dropped

    def test_an_unlisted_ligand_goes_too(self, tmp_path):
        # default is to keep nothing het: a ligand left in the channel blocks
        # the route whose geometry we are measuring
        _a, rep = prepare(self._with(tmp_path, "NAP"))
        assert "NAP" in rep.ligands_dropped

    def test_a_named_ligand_stays(self, tmp_path):
        _a, rep = prepare(self._with(tmp_path, "NAP"),
                          PrepParams(keep_ligands=["NAP"]))
        assert rep.ligands_kept == {"NAP": 1}

    def test_asking_for_a_ligand_that_is_not_there_is_reported(self, tmp_path):
        _a, rep = prepare(self._with(tmp_path, "NAP"),
                          PrepParams(keep_ligands=["NAP", "FAD"]))
        assert any("FAD" in n for n in rep.notes)

    def test_drop_beats_keep(self, tmp_path):
        _a, rep = prepare(self._with(tmp_path, "NAP"),
                          PrepParams(keep_ligands=["NAP"], drop_ligands=["NAP"]))
        assert "NAP" in rep.ligands_dropped

    def test_water_is_furniture(self):
        assert {"HOH", "GOL", "SO4"} <= FURNITURE


class TestAltlocs:
    def test_the_highest_occupancy_wins(self, tmp_path):
        s = structure(tmp_path, [*protein("A"),
                                 ("ATOM  ", "CB", "A", "SER", "A", 2,
                                  (0.0, 1.0, 0.0), 0.3),
                                 ("ATOM  ", "CB", "B", "SER", "A", 2,
                                  (0.0, 2.0, 0.0), 0.7)])
        atoms, rep = prepare(s)
        cb = [a for a in atoms.atoms if a.name == "CB"]
        assert len(cb) == 1 and cb[0].xyz[1] == pytest.approx(2.0)
        assert rep.altlocs_resolved == 1

    def test_both_can_be_kept(self, tmp_path):
        s = structure(tmp_path, [*protein("A"),
                                 ("ATOM  ", "CB", "A", "SER", "A", 2,
                                  (0.0, 1.0, 0.0), 0.3),
                                 ("ATOM  ", "CB", "B", "SER", "A", 2,
                                  (0.0, 2.0, 0.0), 0.7)])
        atoms, _rep = prepare(s, PrepParams(keep_altlocs=True))
        assert len([a for a in atoms.atoms if a.name == "CB"]) == 2


class TestCoordinatesAreNeverMoved:
    """The guarantee the ChannelsDB branch depends on.

    A tunnel fetched from a database is in the deposited entry's frame.
    Moving the protein invalidates it, and the result looks like a 20 A
    disagreement rather than like the mistake it is.
    """

    def test_prep_moves_nothing(self, tmp_path):
        s = structure(tmp_path, [*protein("A"),
                                 ("HETATM", "C1", "", "NAP", "A", 99,
                                  (9.0, 1.0, 2.0), 1.0)])
        before = {(a.chain, a.resseq, a.name): np.array(a.xyz)
                  for a in read_pdb(s).atoms}
        atoms, _rep = prepare(s, PrepParams(keep_ligands=["NAP"]))
        for a in atoms.atoms:
            assert np.allclose(a.xyz, before[(a.chain, a.resseq, a.name)])

    def test_it_survives_a_write_and_read(self, tmp_path):
        s = structure(tmp_path, [*protein("A"),
                                 ("HETATM", "C1", "", "NAP", "A", 99,
                                  (9.123, 1.456, 2.789), 1.0)])
        atoms, _rep = prepare(s, PrepParams(keep_ligands=["NAP"]))
        out = write_pdb(str(tmp_path / "out.pdb"), atoms)
        back = {(a.chain, a.resseq, a.name): np.array(a.xyz)
                for a in read_pdb(out).atoms}
        for a in atoms.atoms:
            assert np.allclose(a.xyz, back[(a.chain, a.resseq, a.name)], atol=1e-3)

    def test_the_file_records_that_it_was_prepared(self, tmp_path):
        atoms, _rep = prepare(structure(tmp_path, protein("A")))
        out = write_pdb(str(tmp_path / "o.pdb"), atoms, header=["chains A"])
        assert "REMARK" in open(out, encoding="utf-8").read()

    def test_it_is_written_LF_only(self, tmp_path):
        # a CR sneaks into a field a downstream tool then fails to parse,
        # and the failure is a long way from the cause
        atoms, _rep = prepare(structure(tmp_path, protein("A")))
        out = write_pdb(str(tmp_path / "o.pdb"), atoms)
        assert b"\r" not in open(out, "rb").read()


class TestRelabel:
    """Renaming the kept chain.

    Needed more often than it looks: most deposited entries of a given
    enzyme have only chain A, so a pair built from two of them has both
    enzymes on A and `chbridge validate` rejects it -- the pose stage cannot
    tell them apart.
    """

    def test_the_kept_chain_is_renamed(self, tmp_path):
        s = structure(tmp_path, protein("A"))
        atoms, rep = prepare(s, PrepParams(chains=["A"], relabel="B"))
        assert {a.chain for a in atoms.atoms} == {"B"}
        assert rep.relabelled_to == "B"
        assert "renamed to B" in rep.summary()

    def test_ligands_come_with_it(self, tmp_path):
        s = structure(tmp_path, [
            *protein("A"),
            ("HETATM", "C1", "", "NAP", "A", 50, (1.0, 1.0, 1.0), 1.0)])
        atoms, _rep = prepare(s, PrepParams(chains=["A"], keep_ligands=["NAP"],
                                            relabel="B"))
        assert {a.chain for a in atoms.atoms} == {"B"}
        assert any(a.resname == "NAP" for a in atoms.atoms)

    def test_renaming_does_not_move_anything(self, tmp_path):
        # the whole point: a renamed chain is still in the deposited frame,
        # so a ChannelsDB tunnel fetched for that entry still applies
        s = structure(tmp_path, protein("A"))
        before = np.array([a.xyz for a in read_pdb(s).atoms])
        atoms, _rep = prepare(s, PrepParams(chains=["A"], relabel="B"))
        after = np.array([a.xyz for a in atoms.atoms])
        assert np.abs(after - before).max() == 0.0

    def test_merging_two_chains_onto_one_id_is_refused(self, tmp_path):
        # it would fuse two protomers into one object
        s = structure(tmp_path, protein("A") + protein("B", 10))
        with pytest.raises(ValueError, match="separate protomers"):
            prepare(s, PrepParams(chains=["A", "B"], relabel="C"))

    def test_a_multi_character_id_is_refused(self, tmp_path):
        s = structure(tmp_path, protein("A"))
        with pytest.raises(ValueError, match="single chain id"):
            prepare(s, PrepParams(chains=["A"], relabel="BB"))

    def test_no_relabel_leaves_the_id_alone(self, tmp_path):
        atoms, rep = prepare(structure(tmp_path, protein("A")))
        assert {a.chain for a in atoms.atoms} == {"A"}
        assert rep.relabelled_to is None
        assert "renamed" not in rep.summary()
