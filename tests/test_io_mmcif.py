"""mmCIF reading: the auth/label choice, quoting, models, PDB column limits.

The reader exists so the prep stage can convert; every test here is about a
way the conversion can be wrong while still producing a file.
"""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.io import read_structure
from channelbridge.io.mmcif import is_mmcif, pdb_column_problems, read_mmcif
from channelbridge.io.pdb import Atom

HEAD = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num"""


def row(**kw):
    f = dict(group="ATOM", id=1, el="C", atom="CA", alt=".", comp="ALA",
             asym="A", seq="1", ins="?", x=1.0, y=2.0, z=3.0, occ=1.0, b=20.0,
             auth_seq="101", auth_comp="ALA", auth_asym="P", auth_atom="CA",
             model="1")
    f.update(kw)
    return ("%(group)s %(id)s %(el)s %(atom)s %(alt)s %(comp)s %(asym)s "
            "%(seq)s %(ins)s %(x).3f %(y).3f %(z).3f %(occ).2f %(b).2f "
            "%(auth_seq)s %(auth_comp)s %(auth_asym)s %(auth_atom)s "
            "%(model)s" % f)


def cif(tmp_path, rows, name="s.cif"):
    p = tmp_path / name
    p.write_text(HEAD + "\n" + "\n".join(rows) + "\n#\n", encoding="utf-8")
    return str(p)


class TestAuthNumberingWins:
    """The decision the module exists to get right.

    `label_*` numbers a polymer from 1 and gives every non-polymer entity its
    own asym id; `auth_*` carries the depositor's chain letters and residue
    numbers -- the ones in the paper, in ChannelsDB and in the `--chain` a
    user types. Reading `label_*` turns a kept cofactor into a dropped one.
    """

    def test_residue_number_and_chain_come_from_auth(self, tmp_path):
        a = read_mmcif(cif(tmp_path, [row()])).atoms[0]
        assert (a.chain, a.resseq) == ("P", 101)

    def test_label_is_the_fallback_not_the_source(self, tmp_path):
        # A writer that omits the auth columns entirely: then label_* is all
        # there is, and dropping the atom would be worse than using it.
        head = "\n".join(line for line in HEAD.splitlines()
                         if not line.startswith("_atom_site.auth_")
                         and not line.endswith("pdbx_PDB_model_num"))
        p = tmp_path / "noauth.cif"
        body = "ATOM 1 C CA . ALA A 7 ? 1.000 2.000 3.000 1.00 20.00"
        p.write_text(head + "\n" + body + "\n#\n", encoding="utf-8")
        a = read_mmcif(str(p)).atoms[0]
        assert (a.chain, a.resseq, a.name) == ("A", 7, "CA")

    def test_a_ligand_in_its_own_label_asym_keeps_the_auth_chain(self, tmp_path):
        # The failure this prevents: in a deposited entry the haem is
        # label_asym B but auth_asym A, so reading label_* puts it in a chain
        # the user did not keep and the site measures as empty.
        rows = [row(), row(group="HETATM", id=2, comp="HEM", asym="B",
                           seq=".", auth_comp="HEM", auth_seq="400",
                           auth_asym="P", atom="FE", auth_atom="FE", el="FE")]
        atoms = read_mmcif(cif(tmp_path, rows)).atoms
        assert [x.chain for x in atoms] == ["P", "P"]
        assert atoms[1].hetatm and atoms[1].resseq == 400

    def test_a_null_label_seq_id_does_not_crash_on_int(self, tmp_path):
        # label_seq_id is "." for everything that is not a polymer residue.
        rows = [row(group="HETATM", seq=".", auth_seq="500", comp="FPP",
                    auth_comp="FPP")]
        assert read_mmcif(cif(tmp_path, rows)).atoms[0].resseq == 500


class TestNulls:
    def test_dot_and_question_become_empty_not_literal(self, tmp_path):
        a = read_mmcif(cif(tmp_path, [row(alt=".", ins="?")])).atoms[0]
        assert a.altloc == "" and a.icode == ""

    def test_a_real_altloc_and_icode_survive(self, tmp_path):
        a = read_mmcif(cif(tmp_path, [row(alt="B", ins="A")])).atoms[0]
        assert a.altloc == "B" and a.icode == "A"


class TestQuoting:
    def test_a_quoted_atom_name_does_not_shift_the_columns(self, tmp_path):
        # Unquoted, a name with a space would push Cartn_x into Cartn_y and
        # the structure would read as a different shape, not as an error.
        a = read_mmcif(cif(tmp_path, [row(atom='"O5\'"', auth_atom='"O5\'"')])).atoms[0]
        assert a.name == "O5'"
        np.testing.assert_allclose(a.xyz, [1.0, 2.0, 3.0])

    def test_an_apostrophe_inside_a_single_quoted_value(self, tmp_path):
        # mmCIF ends a quoted value at quote-then-whitespace, so the
        # apostrophe in O5' is part of the name.
        a = read_mmcif(cif(tmp_path, [row(atom="'O5'A'", auth_atom="'O5'A'")])).atoms[0]
        assert a.name == "O5'A"


class TestModels:
    def test_only_the_first_model_by_default(self, tmp_path):
        rows = [row(id=1, model="1"), row(id=2, model="1"),
                row(id=3, model="2", x=99.0)]
        atoms = read_mmcif(cif(tmp_path, rows))
        assert len(atoms) == 2
        assert atoms.xyz[:, 0].max() < 50.0

    def test_all_models_on_request(self, tmp_path):
        rows = [row(id=1, model="1"), row(id=2, model="2")]
        assert len(read_mmcif(cif(tmp_path, rows), models="all")) == 2

    def test_models_argument_is_checked(self, tmp_path):
        with pytest.raises(ValueError, match=r"first.*all"):
            read_mmcif(cif(tmp_path, [row()]), models="last")


class TestFailuresAreLoud:
    def test_a_file_with_no_atom_site_is_an_error(self, tmp_path):
        p = tmp_path / "meta.cif"
        p.write_text("data_x\n_entry.id X\n#\n", encoding="utf-8")
        with pytest.raises(ValueError, match="no atom_site"):
            read_mmcif(str(p))

    def test_a_non_numeric_coordinate_is_an_error_not_a_skip(self, tmp_path):
        # Skipping it quietly would delete part of a protein and the channel
        # would run through the hole.
        with pytest.raises(ValueError, match="cannot parse coordinates"):
            read_mmcif(cif(tmp_path, [row().replace("1.000", "whoops")]))

    def test_a_truncated_last_row_is_an_error(self, tmp_path):
        p = tmp_path / "short.cif"
        p.write_text(HEAD + "\n" + row() + "\nATOM 2 C CA .\n#\n",
                     encoding="utf-8")
        with pytest.raises(ValueError, match="truncated"):
            read_mmcif(str(p))


class TestDispatch:
    def test_is_mmcif_by_extension(self):
        assert is_mmcif("x.cif") and is_mmcif("X.MMCIF")
        assert not is_mmcif("x.pdb") and not is_mmcif("x.cif.gz")

    def test_read_structure_picks_the_reader(self, tmp_path):
        c = read_structure(cif(tmp_path, [row()]))
        pdb = tmp_path / "s.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA P 101       1.000   2.000   3.000  1.00 20.00"
            "           C\nEND\n", encoding="utf-8")
        p = read_structure(str(pdb))
        assert len(c) == len(p) == 1
        for a, b in ((c.atoms[0], p.atoms[0]),):
            assert (a.chain, a.resseq, a.resname, a.name) == \
                   (b.chain, b.resseq, b.resname, b.name)
            np.testing.assert_allclose(a.xyz, b.xyz)


def atom(**kw):
    f = dict(serial=1, name="CA", altloc="", resname="ALA", chain="A",
             resseq=1, icode="", xyz=(0.0, 0.0, 0.0), occupancy=1.0,
             bfactor=0.0, element="C", hetatm=False)
    f.update(kw)
    return Atom(**f)


class TestPdbColumnProblems:
    """mmCIF has no column widths and PDB does, so the limits are reachable."""

    def test_nothing_wrong_with_an_ordinary_atom(self):
        assert pdb_column_problems([atom()]) == []

    def test_a_two_character_chain_id_is_reported(self):
        # AlphaFill really does emit chains AA..AW: each transplant gets its
        # own asym id and they run past Z.
        (msg,) = pdb_column_problems([atom(chain="AA")])
        assert "AA" in msg and "--relabel" in msg

    def test_a_residue_number_past_the_field_is_reported(self):
        (msg,) = pdb_column_problems([atom(resseq=12345)])
        assert "-999..9999" in msg

    def test_a_long_atom_name_is_reported(self):
        (msg,) = pdb_column_problems([atom(name="HD123")])
        assert "HD123" in msg

    def test_a_serial_past_the_field_is_reported(self):
        (msg,) = pdb_column_problems([atom(serial=100001)])
        assert "99999" in msg
