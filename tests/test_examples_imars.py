"""The decision helpers inside ``examples/iMARS/run_examples.py``.

An example script is normally not worth a test file. These four functions
are, because each of them is a place where the pipeline gave a wrong answer
without failing, and two of them gave the SAME wrong answer twice:

* `_cargo_seed` -- a start point averaged over several copies of the cargo
  lands between them and in none of them. Found first on a kept dimer, where
  it was loud (CAVER returned nothing), then again on two copies in one
  chain, where it was silent: CAVER ran from a point 4.1 A outside either
  site and reported tunnels of it.
* `_caver_exclude` -- the cargo was excluded only when it also supplied the
  start point, so four enzymes were measured with their own cargo asked to
  sit in the channel.
* `_check_loaded` -- the bound that catches an exclusion that did not take.
* `shared_entries` -- keyed on the superseded table, it reported four shared
  structures that no longer existed and cost three pairs of the comparison.

The script is loaded by path: it is not a package, and importing it the way
a user runs it is part of what is being checked.
"""

from __future__ import annotations

import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "examples", "iMARS",
                      "run_examples.py")


@pytest.fixture(scope="module")
def ex():
    spec = importlib.util.spec_from_file_location("imars_run_examples", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def Atom(chain, resseq, resname, xyz, hetatm=True, occupancy=1.0,
         name="C", element="C"):
    """One real atom record, with the fields these helpers do not care about
    filled in. The package's own dataclass, so a field added to it shows up
    here as a signature error rather than as a stub that silently diverges.
    """
    from channelbridge.io.pdb import Atom as Rec

    return Rec(serial=-1, name=name, altloc="", resname=resname, chain=chain,
               resseq=resseq, icode="", xyz=tuple(float(v) for v in xyz),
               occupancy=occupancy, bfactor=0.0, element=element,
               hetatm=hetatm)


class Struct:
    def __init__(self, atoms):
        self.atoms = atoms


def cargo(chain, resseq, centre):
    """Three atoms whose centroid is *centre*."""
    x, y, z = centre
    return [Atom(chain, resseq, "LIG", (x - 1, y, z)),
            Atom(chain, resseq, "LIG", (x + 1, y, z)),
            Atom(chain, resseq, "LIG", (x, y, z))]


class TestCargoSeed:
    def test_one_copy_gives_its_own_centroid(self, ex):
        st = Struct(cargo("A", 1, (10.0, 0.0, 0.0)))
        xyz, how = ex._cargo_seed(st, "LIG", "A")
        assert xyz == pytest.approx([10.0, 0.0, 0.0])
        assert "chain A" in how

    def test_two_copies_do_not_give_the_midpoint(self, ex):
        """The whole bug, in one assertion.

        Two copies 20 A apart. The pooled centroid is at x = 10, which is in
        neither site; the answer has to be one of the two.
        """
        st = Struct([*cargo("A", 1, (0.0, 0.0, 0.0)),
                     *cargo("B", 1, (20.0, 0.0, 0.0))])
        xyz, _how = ex._cargo_seed(st, "LIG", "A")
        assert xyz == pytest.approx([0.0, 0.0, 0.0])
        assert xyz[0] != pytest.approx(10.0)

    def test_the_chain_being_measured_wins(self, ex):
        st = Struct([*cargo("A", 1, (0.0, 0.0, 0.0)),
                     *cargo("B", 1, (20.0, 0.0, 0.0))])
        xyz, how = ex._cargo_seed(st, "LIG", "B")
        assert xyz == pytest.approx([20.0, 0.0, 0.0])
        assert "chain B" in how

    def test_two_copies_in_ONE_chain_pick_one_and_admit_it(self, ex):
        """Taxadiene synthase: 3P5P has two A3C copies in one chain.

        Nothing in the file prefers either, so the lower residue number is
        taken and the sentence says the choice is arbitrary -- rather than
        averaging them, which is what produced a tunnel of the space in
        between.
        """
        st = Struct([*cargo("B", 911, (0.0, 0.0, 0.0)),
                     *cargo("B", 912, (8.2, 0.0, 0.0))])
        xyz, how = ex._cargo_seed(st, "LIG", "B")
        assert xyz == pytest.approx([0.0, 0.0, 0.0])
        assert "arbitrary" in how
        assert "911" in how

    def test_an_atom_record_counts_as_cargo(self, ex):
        """AlphaFill writes a transplant as ATOM, not HETATM."""
        atoms = [Atom("A", 1, "LIG", (5.0, 0.0, 0.0), hetatm=False)]
        xyz, _how = ex._cargo_seed(Struct(atoms), "LIG", "A")
        assert xyz == pytest.approx([5.0, 0.0, 0.0])

    def test_absent_cargo_is_None_not_a_guess(self, ex):
        assert ex._cargo_seed(Struct([]), "LIG", "A") is None


class TestCaverExclude:
    def test_waters_always(self, ex):
        assert ex._caver_exclude("") == ["HOH", "WAT"]

    def test_the_cargo_whenever_there_is_one(self, ex):
        """Not only when it also gives the start point."""
        assert ex._caver_exclude("GST") == ["HOH", "WAT", "GST"]

    def test_a_cargo_that_is_water_is_not_listed_twice(self, ex):
        assert ex._caver_exclude("HOH") == ["HOH", "WAT"]


class TestCheckLoaded:
    def _run(self, ex, tmp_path, atoms, loaded, exclude):
        from channelbridge.io.pdb import write_pdb
        struct = str(tmp_path / "s.pdb")
        write_pdb(struct, atoms)
        run_dir = tmp_path / "caver"
        run_dir.mkdir()
        (run_dir / "log.txt").write_text(
            "* Processing 1 *\nResidues: 1\nAtoms: %d\n" % loaded,
            encoding="utf-8")
        return ex._check_loaded(str(run_dir), struct, exclude)

    def test_a_deficit_is_reported_not_judged(self, ex, tmp_path):
        """CAVER drops atoms of its own accord and does not say so."""
        atoms = [*cargo("A", 1, (0.0, 0.0, 0.0)),
                 *cargo("A", 2, (9.0, 0.0, 0.0))]
        got = self._run(ex, tmp_path, atoms, loaded=3, exclude=["HOH", "WAT"])
        assert "3 of 6" in got
        assert "dropped by CAVER itself" in got

    def test_the_excluded_cargo_is_out_of_the_ceiling(self, ex, tmp_path):
        atoms = [*cargo("A", 1, (0.0, 0.0, 0.0)),
                 Atom("A", 2, "ALA", (9.0, 0.0, 0.0), hetatm=False)]
        got = self._run(ex, tmp_path, atoms, loaded=1,
                        exclude=["HOH", "WAT", "LIG"])
        assert "1 of 1" in got
        assert "3 excluded by name" in got

    def test_a_surplus_is_an_error(self, ex, tmp_path):
        """An exclusion that did not take is the failure this exists for."""
        atoms = cargo("A", 1, (0.0, 0.0, 0.0))
        with pytest.raises(RuntimeError, match="meant to be kept out"):
            self._run(ex, tmp_path, atoms, loaded=3,
                      exclude=["HOH", "WAT", "LIG"])

    def test_low_occupancy_does_not_shrink_the_ceiling(self, ex, tmp_path):
        """Measured: CAVER loads half-occupancy altloc atoms.

        Treating occupancy as a cutoff made the bound too tight and failed a
        run that was correct -- so the bound counts every atom in the file.
        """
        atoms = [Atom("A", 1, "LIG", (float(i), 0.0, 0.0), occupancy=0.5)
                 for i in range(3)]
        got = self._run(ex, tmp_path, atoms, loaded=3, exclude=["HOH", "WAT"])
        assert "3 of 3" in got


class TestSharedEntries:
    def test_the_real_set_has_no_duplicate_computations(self, ex):
        """19 pairs over 24 distinct inputs, not 16 over 20.

        Keyed on the superseded `enzymes.tsv` this reported four shared
        structures; `inputs.tsv` gave every one of those enzymes a model of
        its own, so the answer is none.
        """
        assert ex.shared_entries() == {}
        assert ex.duplicate_pairs() == []

    def test_it_keys_on_the_whole_computation(self, ex):
        """Two enzymes reading one file are not one computation unless every
        decision about it also matches."""
        rows = {
            "X": dict(file="f.cif", chains="A", cargo="", keep_also="",
                      trim="", start_mode="cavity", start_residues="1,2"),
            "Y": dict(file="f.cif", chains="A", cargo="", keep_also="",
                      trim="", start_mode="cavity", start_residues="1,2"),
            "Z": dict(file="f.cif", chains="A+B", cargo="", keep_also="",
                      trim="", start_mode="cavity", start_residues="1,2"),
        }
        real, ex.inputs = ex.inputs, lambda: rows
        try:
            got = ex.shared_entries()
        finally:
            ex.inputs = real
        assert list(got.values()) == [["X", "Y"]]
