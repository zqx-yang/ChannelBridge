"""PDB reading: fixed columns, and CAVER's radius-in-B-factor convention."""

from __future__ import annotations

import numpy as np
import pytest

from channelbridge.io.pdb import has_crlf, read_pdb, write_pdb

# A real CAVER tunnel cluster line. Occupancy is blank; 3.20 is the radius.
CAVER_LINE = "ATOM      1  H   FIL T   1      42.446  12.914  -6.843        3.20"

# A protein line with a four-character atom name, which fills columns 13-16.
WIDE_NAME = "ATOM     12 HD13 LEU B  49      10.000  20.000  30.000  1.00 25.00           H"


def write(tmp_path, name, lines, newline="\n"):
    p = tmp_path / name
    with open(p, "w", newline=newline) as fh:
        fh.write("\n".join(lines) + "\n")
    return p


class TestCaverRadius:
    def test_radius_comes_from_the_b_factor_column(self):
        import pathlib
        import tempfile

        with tempfile.TemporaryDirectory() as d:
            p = write(pathlib.Path(d), "t.pdb", [CAVER_LINE])
            atoms = read_pdb(p)
            assert len(atoms) == 1
            assert atoms.atoms[0].radius == pytest.approx(3.20)
            assert atoms.atoms[0].occupancy is None
            np.testing.assert_allclose(atoms.xyz[0], [42.446, 12.914, -6.843])

    def test_radii_property_matches(self, tmp_path):
        p = write(tmp_path, "t.pdb", [CAVER_LINE, CAVER_LINE.replace("3.20", "2.91")])
        assert read_pdb(p).radii.tolist() == pytest.approx([3.20, 2.91])

    def test_a_tunnel_of_zero_radii_is_rejected(self, tmp_path):
        from channelbridge.tunnel import Tunnel

        p = write(tmp_path, "t.pdb", [CAVER_LINE.replace("3.20", "0.00")] * 3)
        with pytest.raises(ValueError, match="B-factor"):
            Tunnel.from_caver_cluster(p)

    def test_a_single_sphere_cluster_is_rejected_with_a_useful_message(self, tmp_path):
        from channelbridge.tunnel import Tunnel

        p = write(tmp_path, "t.pdb", [CAVER_LINE])
        with pytest.raises(ValueError, match="no pathway"):
            Tunnel.from_caver_cluster(p)


class TestFixedColumns:
    def test_four_character_atom_name_does_not_shift_the_residue(self, tmp_path):
        p = write(tmp_path, "t.pdb", [WIDE_NAME])
        a = read_pdb(p).atoms[0]
        assert a.name == "HD13"
        assert a.resname == "LEU"
        assert a.chain == "B"
        assert a.resseq == 49
        assert a.occupancy == pytest.approx(1.00)
        assert a.bfactor == pytest.approx(25.00)

    def test_negative_coordinates_parse(self, tmp_path):
        p = write(tmp_path, "t.pdb", [CAVER_LINE])
        assert read_pdb(p).xyz[0][2] == pytest.approx(-6.843)

    def test_malformed_coordinates_raise_rather_than_being_skipped(self, tmp_path):
        bad = "ATOM      1  H   FIL T   1      4x.446  12.914  -6.843        3.20"
        p = write(tmp_path, "t.pdb", [bad])
        with pytest.raises(ValueError, match="cannot parse coordinates"):
            read_pdb(p)

    def test_non_coordinate_records_are_ignored(self, tmp_path):
        p = write(tmp_path, "t.pdb", ["REMARK  hello", "TER", CAVER_LINE, "END"])
        assert len(read_pdb(p)) == 1


class TestModels:
    def test_first_model_only_by_default(self, tmp_path):
        p = write(
            tmp_path,
            "t.pdb",
            ["MODEL        5", CAVER_LINE, "ENDMDL", "MODEL        6", CAVER_LINE, "ENDMDL"],
        )
        assert len(read_pdb(p)) == 1

    def test_all_models_when_asked(self, tmp_path):
        p = write(
            tmp_path,
            "t.pdb",
            ["MODEL        5", CAVER_LINE, "ENDMDL", "MODEL        6", CAVER_LINE, "ENDMDL"],
        )
        assert len(read_pdb(p, models="all")) == 2

    def test_a_model_numbered_other_than_one_is_still_read(self, tmp_path):
        # CAVER labels the kept tunnel by its index in the snapshot, e.g. MODEL 5.
        p = write(tmp_path, "t.pdb", ["MODEL        5", CAVER_LINE, "ENDMDL"])
        assert len(read_pdb(p)) == 1

    def test_rejects_a_bad_models_argument(self, tmp_path):
        p = write(tmp_path, "t.pdb", [CAVER_LINE])
        with pytest.raises(ValueError, match="models must be"):
            read_pdb(p, models="every")


class TestLineEndings:
    def test_crlf_is_detected(self, tmp_path):
        # CRLF makes CAVER report "No PDB files found" and exit 0.
        p = write(tmp_path, "crlf.pdb", [CAVER_LINE], newline="\r\n")
        assert has_crlf(p) is True

    def test_lf_is_clean(self, tmp_path):
        p = write(tmp_path, "lf.pdb", [CAVER_LINE], newline="\n")
        assert has_crlf(p) is False

    def test_write_pdb_always_writes_lf(self, tmp_path):
        src = write(tmp_path, "in.pdb", [WIDE_NAME])
        atoms = read_pdb(src)
        out = tmp_path / "out.pdb"
        assert write_pdb(out, atoms) == 1
        assert has_crlf(out) is False

    def test_roundtrip_preserves_the_fields_we_use(self, tmp_path):
        src = write(tmp_path, "in.pdb", [WIDE_NAME])
        a0 = read_pdb(src).atoms[0]
        out = tmp_path / "out.pdb"
        write_pdb(out, [a0])
        a1 = read_pdb(out).atoms[0]
        assert (a1.name, a1.resname, a1.chain, a1.resseq) == (
            a0.name,
            a0.resname,
            a0.chain,
            a0.resseq,
        )
        np.testing.assert_allclose(a1.xyz, a0.xyz, atol=1e-3)


class TestSelection:
    def test_exclude_resnames_is_case_insensitive_and_trims(self, tmp_path):
        lines = [
            "HETATM    1  C1  BGC B 262      10.000  10.000  10.000  1.00 20.00           C",
            "HETATM    2  PA  NAP B 601      11.000  10.000  10.000  1.00 20.00           P",
            "ATOM      3  CA  ALA B   1      12.000  10.000  10.000  1.00 20.00           C",
        ]
        p = write(tmp_path, "t.pdb", lines)
        arr = read_pdb(p)
        assert len(arr) == 3
        assert len(arr.exclude_resnames([" bgc ", "nap"])) == 1

    def test_select_by_chain_resname_and_record(self, tmp_path):
        lines = [
            "HETATM    1  C1  BGC B 262      10.000  10.000  10.000  1.00 20.00           C",
            "ATOM      2  CA  ALA C   1      12.000  10.000  10.000  1.00 20.00           C",
        ]
        p = write(tmp_path, "t.pdb", lines)
        arr = read_pdb(p)
        assert len(arr.select(chain="B")) == 1
        assert len(arr.select(resname="bgc")) == 1
        assert len(arr.select(hetatm=False)) == 1

    def test_centroid_of_empty_raises(self):
        from channelbridge.io.pdb import AtomArray

        with pytest.raises(ValueError, match="empty"):
            AtomArray([]).centroid()
