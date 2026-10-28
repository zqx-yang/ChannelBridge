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


def ligand(resname, ch, seq, xyz, n=2):
    """A het residue of *n* atoms starting at *xyz*, 1 A apart along x."""
    return [("HETATM", "C%d" % i, "", resname, ch, seq,
             (xyz[0] + i, xyz[1], xyz[2]), 1.0) for i in range(n)]


class TestLigandsFromAnotherChain:
    """A modelled structure does not put a ligand in its protein's chain.

    AlphaFill gives every transplant its own chain id, so chain-based
    selection alone keeps the protein and loses the cofactor -- and the report
    then says "ligands kept: none" while the channel being measured runs
    through the site whose occupant was the point of the measurement. Contact
    decides instead.
    """

    def protein_and_ligands(self, tmp_path, copies):
        rows = protein("A", 1, 4)
        for ch, xyz in copies:
            rows += ligand("FPP", ch, 1, xyz)
        return structure(tmp_path, rows)

    def test_a_touching_ligand_is_adopted(self, tmp_path):
        s = self.protein_and_ligands(tmp_path, [("S", (1.0, 2.0, 0.0))])
        _, rep = prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"]))
        assert rep.ligands_kept == {"FPP": 2}
        assert [r[:3] for r in rep.ligand_copies] == [("FPP", "S", 1)]
        assert rep.ligand_copies[0][4] is True

    def test_a_distant_copy_is_not_adopted(self, tmp_path):
        s = self.protein_and_ligands(tmp_path, [("S", (0.0, 50.0, 0.0))])
        _, rep = prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"]))
        assert rep.ligands_kept == {}
        # it is still reported, so "none kept" is explained rather than silent
        assert rep.ligand_copies[0][3] == 0 and rep.ligand_copies[0][4] is False

    def test_only_a_named_ligand_is_adopted(self, tmp_path):
        # A neighbouring protomer's cofactor is not pulled in: nobody asked
        # for it, and it belongs to the chain that was dropped.
        rows = protein("A", 1, 4) + ligand("HEM", "S", 1, (1.0, 2.0, 0.0))
        _, rep = prepare(structure(tmp_path, rows),
                         PrepParams(chains=["A"], keep_ligands=["FPP"]))
        assert rep.ligands_kept == {} and rep.ligand_copies == []

    def test_contact_distance_is_a_flag(self, tmp_path):
        s = self.protein_and_ligands(tmp_path, [("S", (1.0, 6.0, 0.0))])
        kw = dict(chains=["A"], keep_ligands=["FPP"])
        assert prepare(s, PrepParams(ligand_contact=4.0, **kw))[1].ligands_kept == {}
        assert prepare(s, PrepParams(ligand_contact=8.0, **kw))[1].ligands_kept == {"FPP": 2}

    def test_all_copies_by_default(self, tmp_path):
        s = self.protein_and_ligands(
            tmp_path, [("S", (1.0, 2.0, 0.0)), ("T", (2.0, 2.5, 0.0))])
        _, rep = prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"]))
        assert rep.ligands_kept == {"FPP": 4}

    def test_best_keeps_the_most_contacted_copy(self, tmp_path):
        # Four donors give four copies of one ligand; the spare ones are walls
        # in somebody's channel, so 'best' exists.
        s = self.protein_and_ligands(
            tmp_path, [("S", (0.0, 3.9, 0.0)), ("T", (1.0, 1.0, 0.0))])
        _, rep = prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"],
                                       ligand_copies="best"))
        assert rep.ligands_kept == {"FPP": 2}
        kept = [r for r in rep.ligand_copies if r[4]]
        assert len(kept) == 1 and kept[0][1] == "T"
        assert max(r[3] for r in rep.ligand_copies) == kept[0][3]

    def test_an_unknown_copy_policy_is_refused(self, tmp_path):
        s = self.protein_and_ligands(tmp_path, [("S", (1.0, 2.0, 0.0))])
        with pytest.raises(ValueError, match=r"all.*best"):
            prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"],
                                  ligand_copies="first"))

    def test_relabel_renumbers_a_colliding_ligand(self, tmp_path):
        # AlphaFill numbers every transplant 1, so moving it into the
        # protein's chain would give it the address of residue 1. Downstream a
        # residue is selected by chain and number, so this is not cosmetic.
        s = self.protein_and_ligands(tmp_path, [("S", (1.0, 2.0, 0.0))])
        atoms, rep = prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"],
                                           relabel="U"))
        assert rep.ligands_renumbered == [("FPP", "S", 1, 5)]
        het = [a for a in atoms.atoms if a.hetatm]
        assert {a.chain for a in het} == {"U"}
        assert {a.resseq for a in het} == {5}
        poly = {(a.chain, a.resseq) for a in atoms.atoms if not a.hetatm}
        assert ("U", 5) not in poly

    def test_a_ligand_that_does_not_collide_keeps_its_number(self, tmp_path):
        rows = protein("A", 1, 4) + ligand("FPP", "S", 900, (1.0, 2.0, 0.0))
        atoms, rep = prepare(structure(tmp_path, rows),
                             PrepParams(chains=["A"], keep_ligands=["FPP"],
                                        relabel="U"))
        assert rep.ligands_renumbered == []
        assert {a.resseq for a in atoms.atoms if a.hetatm} == {900}

    def test_adoption_moves_nothing(self, tmp_path):
        s = self.protein_and_ligands(tmp_path, [("S", (1.0, 2.0, 0.0))])
        before = read_pdb(s)
        atoms, _ = prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"],
                                         relabel="U"))
        want = {tuple(round(v, 3) for v in a.xyz)
                for a in before.atoms if a.hetatm}
        got = {tuple(round(v, 3) for v in a.xyz)
               for a in atoms.atoms if a.hetatm}
        assert got == want

    def test_the_report_names_each_copy_and_its_fate(self, tmp_path):
        s = self.protein_and_ligands(
            tmp_path, [("S", (0.0, 3.9, 0.0)), ("T", (1.0, 1.0, 0.0))])
        _, rep = prepare(s, PrepParams(chains=["A"], keep_ligands=["FPP"],
                                       ligand_copies="best"))
        text = rep.summary()
        assert "from another chain" in text
        assert "redundant copy" in text


class TestRelabelMapping:
    """A kept dimer needs each chain renamed, which is not a merge.

    The pose stage refuses a complex whose two sides share a chain id, so a
    dimeric side has to be given ids the other side does not use -- and a
    single `--relabel` id cannot express that without collapsing the two
    protomers into one chain.
    """

    def two_chains(self, tmp_path):
        return structure(tmp_path, protein("A", 1, 3) + protein("B", 1, 3))

    def test_each_chain_gets_its_own_id(self, tmp_path):
        atoms, rep = prepare(self.two_chains(tmp_path),
                             PrepParams(chains=["A", "B"], relabel="A:U,B:V"))
        got = {(a.chain, a.resseq) for a in atoms.atoms}
        assert {c for c, _ in got} == {"U", "V"}
        assert len(got) == 6
        assert "A->U" in rep.relabelled_to and "B->V" in rep.relabelled_to

    def test_a_single_id_still_works_for_one_chain(self, tmp_path):
        _, rep = prepare(structure(tmp_path, protein("A", 1, 3)),
                         PrepParams(chains=["A"], relabel="U"))
        assert rep.relabelled_to == "U"

    def test_a_single_id_for_two_chains_is_still_refused(self, tmp_path):
        with pytest.raises(ValueError, match="merge"):
            prepare(self.two_chains(tmp_path),
                    PrepParams(chains=["A", "B"], relabel="U"))

    def test_the_refusal_suggests_the_mapping(self, tmp_path):
        with pytest.raises(ValueError, match=r"A:A,B:B"):
            prepare(self.two_chains(tmp_path),
                    PrepParams(chains=["A", "B"], relabel="U"))

    def test_mapping_two_chains_onto_one_id_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="merge"):
            prepare(self.two_chains(tmp_path),
                    PrepParams(chains=["A", "B"], relabel="A:U,B:U"))

    def test_a_kept_chain_left_out_of_the_mapping_is_refused(self, tmp_path):
        # Leaving one out would let it keep an id the mapping just handed to
        # another chain.
        with pytest.raises(ValueError, match="does not name kept chain"):
            prepare(self.two_chains(tmp_path),
                    PrepParams(chains=["A", "B"], relabel="A:U"))

    def test_naming_a_chain_that_is_not_kept_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="not kept"):
            prepare(self.two_chains(tmp_path),
                    PrepParams(chains=["A"], relabel="A:U,B:V"))

    def test_an_adopted_ligand_follows_the_first_kept_chain(self, tmp_path):
        rows = (protein("A", 1, 3) + protein("B", 1, 3)
                + ligand("FPP", "S", 1, (1.0, 2.0, 0.0)))
        atoms, _ = prepare(structure(tmp_path, rows),
                           PrepParams(chains=["A", "B"], keep_ligands=["FPP"],
                                      relabel="A:U,B:V"))
        het = [a for a in atoms.atoms if a.hetatm]
        assert het and {a.chain for a in het} == {"U"}

    def test_a_malformed_mapping_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="not old:new"):
            prepare(self.two_chains(tmp_path),
                    PrepParams(chains=["A", "B"], relabel="A:U,B"))


class TestDropResidues:
    """A disordered terminus is attached, so it cannot be dropped by name."""

    def test_a_range_goes(self, tmp_path):
        s = structure(tmp_path, protein("A", 1, 8))
        atoms, rep = prepare(s, PrepParams(chains=["A"], drop_residues=[(1, 3)]))
        assert {a.resseq for a in atoms.atoms} == {4, 5, 6, 7, 8}
        assert rep.residues_dropped == [(1, 3)]
        assert rep.atoms_dropped_by_range == 3

    def test_several_ranges(self, tmp_path):
        s = structure(tmp_path, protein("A", 1, 8))
        atoms, _ = prepare(s, PrepParams(chains=["A"],
                                       drop_residues=[(1, 2), (7, 8)]))
        assert {a.resseq for a in atoms.atoms} == {3, 4, 5, 6}

    def test_a_kept_ligand_is_not_dropped_by_a_range(self, tmp_path):
        # Ranges are for polymer termini; a ligand's number is unrelated and
        # may well fall inside one.
        rows = protein("A", 1, 5) + ligand("FPP", "A", 2, (1.0, 1.0, 0.0))
        atoms, _ = prepare(structure(tmp_path, rows),
                           PrepParams(chains=["A"], keep_ligands=["FPP"],
                                      drop_residues=[(1, 3)]))
        assert [a.resname for a in atoms.atoms if a.hetatm] == ["FPP", "FPP"]

    def test_a_backwards_range_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="backwards"):
            prepare(structure(tmp_path, protein("A", 1, 8)),
                    PrepParams(chains=["A"], drop_residues=[(8, 1)]))

    def test_the_report_names_what_went(self, tmp_path):
        _, rep = prepare(structure(tmp_path, protein("A", 1, 8)),
                         PrepParams(chains=["A"], drop_residues=[(1, 3)]))
        assert "1-3" in rep.summary() and "residues dropped" in rep.summary()

    def test_dropping_every_polymer_residue_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="nothing is left"):
            prepare(structure(tmp_path, protein("A", 1, 4)),
                    PrepParams(chains=["A"], drop_residues=[(1, 4)]))
