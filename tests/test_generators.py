"""Stage 7: the generator input for a chosen pose.

The tests that matter here are the ones about the FRAME and about NUMBERING.
Both failures are silent: a spec that frees the inter-enzyme placement, or one
whose hotspots are offset, runs to completion and produces designs against a
geometry or an epitope nobody chose.
"""

from __future__ import annotations

import json
import os
from typing import ClassVar

import pytest

from channelbridge.generators import (
    DEFAULT_BINDER_LENGTH,
    NO_TIP,
    TIP_ATOMS,
    TOOLS,
    PoseSpec,
    Target,
    bindcraft_settings,
    contig,
    crop_runs,
    design_jobs,
    job_atoms,
    rfd3_spec,
    rfd4_conditions,
    tools_table,
    write_bindcraft,
    write_flexcraft_args,
    write_for,
    write_rfd3,
    write_rfd4,
    write_rfd4_cif,
    write_rfdiffusion_args,
)


class Atom:
    """Just enough of the io atom record for the mmCIF writer."""

    def __init__(self, name, resname, chain, resseq, xyz, hetatm=False,
                 element=""):
        self.name = name
        self.resname = resname
        self.chain = chain
        self.resseq = resseq
        self.xyz = xyz
        self.hetatm = hetatm
        self.element = element
        self.occupancy = 1.0
        self.bfactor = 0.0


def make_spec(**kw):
    up = Target("A", 1, 381, [117, 118, 121],
                {117: "TYR", 118: "GLU", 121: "ARG"}, "upstream GPPS patch 2")
    dn = Target("B", 1, 598, [441, 512],
                {441: "ILE", 512: "TYR"}, "downstream MYS patch 6")
    args = dict(name="iMARS_02_GPPS-MYS_geom001", complex_pdb="/x/cplx.pdb",
                upstream=up, downstream=dn, geom_id="geom001",
                d_ic=6.0, o_ic=165.0, do=0.9371)
    args.update(kw)
    return PoseSpec(**args)


class TestTipAtoms:
    def test_every_standard_residue_has_tip_atoms(self):
        assert len(TIP_ATOMS) == 20

    def test_the_three_residues_without_a_tip_keyword_are_covered(self):
        # RFdiffusion3 rejects TIP for these, so the table must name atoms.
        for r in NO_TIP:
            assert r in TIP_ATOMS
            assert TIP_ATOMS[r]

    def test_tip_atoms_are_side_chain_not_backbone(self):
        # Gly is the sole exception: it has no side chain, so CA is the tip.
        backbone = {"N", "C", "O"}
        for res, atoms in TIP_ATOMS.items():
            assert not (set(atoms) & backbone), res
            if res != "GLY":
                assert "CA" not in atoms, res


class TestTargetRefuses:
    def test_a_hotspot_outside_the_span_is_an_error_not_a_warning(self):
        with pytest.raises(ValueError, match="outside the chain's span"):
            Target("A", 10, 20, [25], {25: "ALA"})

    def test_a_reversed_span_is_an_error(self):
        with pytest.raises(ValueError, match="before first"):
            Target("A", 20, 10, [])

    def test_a_hotspot_without_a_residue_name_is_an_error(self):
        with pytest.raises(ValueError, match="no residue name"):
            Target("A", 1, 50, [10, 11], {10: "ALA"})

    def test_an_unknown_residue_type_refuses_rather_than_selecting_everything(self):
        t = Target("A", 1, 50, [10], {10: "XYZ"})
        with pytest.raises(ValueError, match="no tip-atom definition"):
            t.tip_atoms(10)


class TestPoseSpecRefuses:
    def test_both_enzymes_on_one_chain_is_an_error(self):
        up = Target("A", 1, 100, [10], {10: "ALA"})
        dn = Target("A", 1, 100, [20], {20: "ALA"})
        with pytest.raises(ValueError, match="both enzymes are on chain A"):
            PoseSpec(name="n", complex_pdb="c.pdb", upstream=up, downstream=dn)

    def test_a_bad_binder_range_is_an_error(self):
        with pytest.raises(ValueError, match="0 < lo <= hi"):
            make_spec(binder_length=(150, 100))


class TestOffsetAndContig:
    def test_a_tool_that_keeps_the_input_numbering_gets_offset_zero(self):
        t = Target("A", 185, 385, [209], {209: "LEU"})
        assert t.offset(False) == 0

    def test_a_tool_that_renumbers_from_one_gets_first_minus_one(self):
        # The exact case that cost another campaign a whole batch: a crop
        # starting at 185 means UniProt 209 is residue 25 after renumbering.
        t = Target("A", 185, 385, [209], {209: "LEU"})
        assert t.offset(True) == 184
        assert 209 - t.offset(True) == 25

    def test_contig_puts_the_binder_first_then_each_target_after_a_break(self):
        s = make_spec()
        assert contig(s, renumbers=False) == "100-150, /0, A1-381, /0, B1-598"

    def test_contig_and_hotspots_are_renumbered_TOGETHER(self):
        # Renumbering the span but not the hotspots (or the reverse) is the
        # failure this pairing exists to prevent.
        up = Target("A", 185, 385, [209], {209: "LEU"})
        dn = Target("B", 50, 250, [100], {100: "TYR"})
        s = make_spec(upstream=up, downstream=dn)
        assert contig(s, renumbers=True) == "100-150, /0, A1-201, /0, B1-201"
        hot = [sel for _k, _v, sels in rfd4_conditions(s)
               for sel in sels if "/" in sel]
        assert "sel('A/*/25/[CD1,CD2]')" in hot      # 209 - 184
        assert "sel('B/*/51/[CZ,OH]')" in hot        # 100 - 49

    def test_the_default_binder_range_is_the_long_class(self):
        assert DEFAULT_BINDER_LENGTH == (100, 150)


class TestRFd3:
    def test_one_key_per_design_task_with_both_spans_in_one_contig(self):
        s = make_spec()
        spec = rfd3_spec(s)
        assert list(spec) == [s.name]
        e = spec[s.name]
        assert e["dialect"] == 2
        assert e["contig"] == "100-150,/0,A1-381,/0,B1-598"
        assert e["infer_ori_strategy"] == "hotspots"
        assert e["is_non_loopy"] is True

    def test_every_hotspot_uses_RFd3s_own_TIP_keyword(self):
        """Never our atom table.

        `TIP` is resolved against RFdiffusion3's own per-residue-type
        definition, so writing atoms out would mean keeping a second copy of
        that table in sync with theirs. The only residues where the two
        could disagree are Ala, Gly and Pro, which have no tip -- and those
        are not hotspots, because `pick_hotspots` chooses from a patch's
        exposed F/Y/W/I/L/M.
        """
        s = make_spec()
        hot = rfd3_spec(s)[s.name]["select_hotspots"]
        assert hot == {"A117": "TIP", "A118": "TIP", "A121": "TIP",
                       "B441": "TIP", "B512": "TIP"}
        assert set(hot.values()) == {"TIP"}

    def test_the_json_carries_the_pose_and_the_frame_warning(self):
        s = make_spec()
        d = json.loads(write_rfd3(s))
        assert "geom001" in d["_comment"]
        assert "DO 0.9371" in d["_comment"]
        assert "FIXED RELATIVE TO EACH OTHER" in d["_comment"]


class TestRFD4:
    def test_the_frame_is_pinned_with_C_CRD_not_C_CRDGRP(self):
        # C_CRDGRP keeps a group rigid internally and lets its global pose
        # float. With two chains that frees the inter-enzyme placement, which
        # IS the result of the pose search.
        keys = [k for k, _v, _s in rfd4_conditions(make_spec())]
        assert "C_CRD" in keys
        assert "C_CRDGRP" not in keys
        # It appears in the file only as the comment saying why it is not
        # used, never as a condition key.
        body = [ln for ln in write_rfd4(make_spec()).splitlines()
                if not ln.startswith("#")]
        assert not [ln for ln in body if "C_CRDGRP" in ln]
        assert [ln for ln in write_rfd4(make_spec()).splitlines()
                if ln.startswith("#") and "C_CRDGRP" in ln]

    def test_C_CRD_and_C_SEQ_cover_both_chains_whole(self):
        conds = {k: sels for k, _v, sels in rfd4_conditions(make_spec())}
        assert conds["C_CRD"] == ["sel('A')", "sel('B')"]
        assert conds["C_SEQ"] == ["sel('A')", "sel('B')"]

    def test_C_HOT_names_tip_atoms_of_every_patch_residue_on_both_chains(self):
        conds = {k: sels for k, _v, sels in rfd4_conditions(make_spec())}
        assert conds["C_HOT"] == [
            "sel('A/*/117/[CZ,OH]')", "sel('A/*/118/[CD,OE1,OE2]')",
            "sel('A/*/121/[CZ,NH1,NH2]')",
            "sel('B/*/441/[CD1]')", "sel('B/*/512/[CZ,OH]')"]

    def test_the_yaml_is_the_documented_shape(self):
        text = write_rfd4(make_spec())
        assert "iMARS_02_GPPS-MYS_geom001:" in text
        assert "  input: iMARS_02_GPPS-MYS_geom001.cif" in text
        assert '  contig: "100-150, /0, A1-381, /0, B1-598"' in text
        assert "    C_CRD: {true: [" in text
        assert text.startswith("# iMARS_02_GPPS-MYS_geom001")

    def test_the_yaml_parses_and_round_trips_if_pyyaml_is_available(self):
        yaml = pytest.importorskip("yaml")
        d = yaml.safe_load(write_rfd4(make_spec()))
        e = d["iMARS_02_GPPS-MYS_geom001"]
        assert e["contig"] == "100-150, /0, A1-381, /0, B1-598"
        assert e["conditions"]["C_CRD"][True] == [
            {"select": "sel('A')"}, {"select": "sel('B')"}]
        assert len(e["conditions"]["C_HOT"][True]) == 5

    def test_the_cif_keeps_the_author_numbering_in_label_seq_id(self):
        # A 1..N counter here would put the patch tables' numbers and the
        # file's numbers one renumbering apart, in opposite directions.
        atoms = [Atom("CA", "TYR", "A", 117, (1.0, 2.0, 3.0), element="C"),
                 Atom("CZ", "TYR", "A", 117, (4.0, 5.0, 6.0), element="C"),
                 Atom("CA", "ILE", "B", 441, (7.0, 8.0, 9.0), element="C")]
        text = write_rfd4_cif(make_spec(), atoms)
        rows = [ln.split() for ln in text.splitlines()
                if ln.startswith(("ATOM", "HETATM"))]
        assert len(rows) == 3
        label_seq = [r[8] for r in rows]
        auth_seq = [r[15] for r in rows]
        assert label_seq == ["117", "117", "441"]
        assert auth_seq == label_seq

    def test_each_chain_gets_its_own_entity_id(self):
        atoms = [Atom("CA", "TYR", "A", 117, (0.0, 0.0, 0.0)),
                 Atom("CA", "ILE", "B", 441, (0.0, 0.0, 0.0))]
        rows = [ln.split() for ln in write_rfd4_cif(make_spec(), atoms)
                .splitlines() if ln.startswith("ATOM")]
        assert rows[0][7] != rows[1][7]

    def test_a_ligand_is_written_as_HETATM(self):
        atoms = [Atom("C1", "GST", "B", 601, (0.0, 0.0, 0.0), hetatm=True)]
        text = write_rfd4_cif(make_spec(), atoms)
        assert text.splitlines()[-2].startswith("HETATM")


class TestBindCraft:
    def test_both_chains_are_named_so_neither_enzyme_is_dropped(self):
        assert bindcraft_settings(make_spec())["chains"] == "A,B"

    def test_hotspots_are_chain_qualified_and_in_the_input_numbering(self):
        s = bindcraft_settings(make_spec())
        assert s["target_hotspot_residues"] == "A117,A118,A121,B441,B512"

    def test_design_path_is_absent_because_the_job_wrapper_writes_it(self):
        assert "design_path" not in bindcraft_settings(make_spec())

    def test_lengths_and_count(self):
        s = bindcraft_settings(make_spec(n_designs=42))
        assert s["lengths"] == [100, 150]
        assert s["number_of_final_designs"] == 42

    def test_the_json_is_loadable(self):
        d = json.loads(write_bindcraft(make_spec()))
        assert d["starting_pdb"] == "/x/cplx.pdb"


class TestJobScriptArgs:
    def test_flexcraft_writes_one_length_and_says_so(self):
        text = write_flexcraft_args(make_spec())
        assert "--num_aa 125" in text          # midpoint of 100-150
        assert "our range is 100-150" in text
        assert '--hotspots "A117,A118,A121,B441,B512"' in text
        assert 'TARGET_PDB="/x/cplx.pdb"' in text

    def test_rfdiffusion_v1_says_it_cannot_express_tip_atoms(self):
        text = write_rfdiffusion_args(make_spec())
        assert "Whole-residue hotspots only" in text
        assert "'ppi.hotspot_res=[A117,A118,A121,B441,B512]'" in text
        assert "contigmap.contigs=[100-150/0 A1-381,B1-598]" in text

    def test_both_carry_the_pose_provenance(self):
        for text in (write_flexcraft_args(make_spec()),
                     write_rfdiffusion_args(make_spec())):
            assert "geom001" in text
            assert "DO 0.9371" in text


class TestToolSurvey:
    def test_the_tools_that_need_a_spec_file_are_the_ones_recorded(self):
        need = {n for n, t in TOOLS.items() if t.spec_file}
        assert need == {"RFdiffusion3", "RFD4-Proteina", "BindCraft",
                        "protpardelle-1c"}

    def test_the_tools_whose_spec_is_job_script_arguments(self):
        args = {n for n, t in TOOLS.items() if not t.spec_file}
        assert args == {"flexcraft", "RFdiffusion", "RFDpoly"}

    def test_every_refusal_says_why(self):
        for t in TOOLS.values():
            if not t.supported:
                assert t.note.startswith("Refused.")
                assert len(t.note) > 60

    def test_every_supported_tool_has_a_writer(self):
        from channelbridge.generators import WRITERS
        assert set(WRITERS) == {n for n, t in TOOLS.items() if t.supported}

    def test_the_table_lists_every_tool(self):
        table = tools_table()
        for n in TOOLS:
            assert n in table


class TestWriteFor:
    def test_writes_the_expected_files(self, tmp_path):
        s = make_spec()
        out = str(tmp_path)
        assert [os.path.basename(p) for p in write_for("RFdiffusion3", s, out)] \
            == ["iMARS_02_GPPS-MYS_geom001.json"]
        atoms = [Atom("CA", "TYR", "A", 117, (0.0, 0.0, 0.0))]
        got = [os.path.basename(p)
               for p in write_for("RFD4-Proteina", s, out, atoms)]
        assert got == ["iMARS_02_GPPS-MYS_geom001.cif",
                       "iMARS_02_GPPS-MYS_geom001.yaml"]

    def test_rfd4_without_atoms_refuses_rather_than_writing_half_an_input(
            self, tmp_path):
        with pytest.raises(ValueError, match="needs the posed complex"):
            write_for("RFD4-Proteina", make_spec(), str(tmp_path))

    def test_an_unsupported_tool_refuses_with_its_reason(self, tmp_path):
        with pytest.raises(ValueError, match="Refused"):
            write_for("protpardelle-1c", make_spec(), str(tmp_path))

    def test_an_unknown_tool_lists_the_known_ones(self, tmp_path):
        with pytest.raises(ValueError, match="known:"):
            write_for("nosuchtool", make_spec(), str(tmp_path))


class Pt:
    """An atom with coordinates, for the cropping tests."""

    def __init__(self, chain, resseq, xyz, hetatm=False, name="CA",
                 resname="ALA"):
        self.chain = chain
        self.resseq = resseq
        self.xyz = xyz
        self.hetatm = hetatm
        self.name = name
        self.resname = resname
        self.element = "C"
        self.occupancy = 1.0
        self.bfactor = 0.0
        self.serial = -1
        self.altloc = ""
        self.icode = ""


class TestDesignModes:
    """One pose, three ways to ask for a binder.

    The three are a ladder: monospecific is easiest and contains none of the
    pose's geometry, context adds the partner as a wall, bispecific asks for
    the bridge directly. A campaign can walk up it.
    """

    def test_bispecific_is_the_pose_itself(self):
        s = make_spec()
        jobs = design_jobs(s, "bispecific")
        assert jobs == [s]
        assert jobs[0].mode == "bispecific"
        assert len(jobs[0].designed) == 2

    def test_monospecific_is_one_job_per_enzyme_with_the_other_absent(self):
        jobs = design_jobs(make_spec(), "monospecific")
        assert len(jobs) == 2
        for j in jobs:
            assert j.mode == "monospecific"
            assert len(j.targets) == 1 and not j.context
        assert {j.targets[0].chain for j in jobs} == {"A", "B"}

    def test_each_monospecific_job_keeps_its_own_hotspots(self):
        up, dn = design_jobs(make_spec(), "monospecific")
        assert list(up.targets[0].hotspots) == [117, 118, 121]
        assert list(dn.targets[0].hotspots) == [441, 512]

    def test_context_needs_the_atoms_to_crop_with(self):
        with pytest.raises(ValueError, match="needs the posed complex"):
            design_jobs(make_spec(), "context")

    def test_an_unknown_mode_is_refused(self):
        with pytest.raises(ValueError, match="unknown mode"):
            design_jobs(make_spec(), "trispecific")

    def test_a_single_target_spec_cannot_be_split(self):
        s = make_spec()
        mono = design_jobs(s, "monospecific")[0]
        with pytest.raises(ValueError, match="single target already"):
            design_jobs(mono, "monospecific")


class TestContextCropping:
    def _atoms(self):
        # chain A hotspot at the origin; chain B runs away along x
        out = [Pt("A", 117, (0.0, 0.0, 0.0))]
        out += [Pt("B", 400 + i, (float(i), 0.0, 0.0)) for i in range(60)]
        return out

    def test_only_what_is_near_the_epitope_is_kept(self):
        runs = crop_runs(self._atoms(), "B", [(0.0, 0.0, 0.0)], cut=20.0)
        assert runs == [(400, 420)]

    def test_the_cut_is_a_distance_not_a_count(self):
        near = crop_runs(self._atoms(), "B", [(0.0, 0.0, 0.0)], cut=10.0)
        far = crop_runs(self._atoms(), "B", [(0.0, 0.0, 0.0)], cut=30.0)
        assert near == [(400, 410)]
        assert far == [(400, 430)]

    def test_a_short_isolated_run_is_dropped(self):
        """BindCraft: an isolated fragment with no core is rearranged
        randomly, so it is not worth keeping."""
        atoms = [Pt("B", 1, (0.0, 0.0, 0.0)), Pt("B", 2, (1.0, 0.0, 0.0)),
                 Pt("B", 50, (0.0, 0.0, 0.0)), Pt("B", 51, (1.0, 0.0, 0.0)),
                 Pt("B", 52, (2.0, 0.0, 0.0)), Pt("B", 53, (3.0, 0.0, 0.0))]
        runs = crop_runs(atoms, "B", [(0.0, 0.0, 0.0)], cut=20.0, min_run=4)
        assert runs == [(50, 53)]

    def test_nothing_near_is_empty_not_an_error(self):
        assert crop_runs(self._atoms(), "B", [(500.0, 0.0, 0.0)],
                         cut=5.0) == []

    def test_context_mode_crops_and_marks_the_partner(self):
        jobs = design_jobs(make_spec(), "context", atoms=self._atoms(),
                           context_cut=20.0)
        job = next(j for j in jobs if j.upstream.chain == "A")
        assert job.mode == "context"
        assert job.downstream.is_context is True
        assert list(job.downstream.hotspots) == []
        assert job.downstream.keep == [(400, 420)]

    def test_a_far_partner_falls_back_to_monospecific_rather_than_failing(self):
        """A measurement, not an error: nothing of the partner is near this
        epitope, so there is no wall to add."""
        atoms = [Pt("A", 117, (0.0, 0.0, 0.0)),
                 Pt("B", 441, (900.0, 0.0, 0.0)),
                 Pt("B", 512, (901.0, 0.0, 0.0))]
        jobs = design_jobs(make_spec(), "context", atoms=atoms,
                           context_cut=20.0)
        job = next(j for j in jobs if j.upstream.chain == "A")
        assert job.downstream is None
        assert job.mode == "monospecific"


class TestCroppedNumbering:
    """A crop is not a constant shift, and the two tools differ."""

    def _cropped(self):
        ctx = Target("B", 400, 600, (), {}, "ctx", keep=[(400, 410),
                                                         (500, 520)],
                     is_context=True)
        return make_spec(downstream=ctx)

    def test_a_renumbering_tool_sees_one_contiguous_span(self):
        """The gaps a crop made do not exist once the preserved residues are
        counted 1..N, so emitting them separately would describe chain
        breaks the tool is not going to make."""
        assert contig(self._cropped(), renumbers=True).endswith("B1-32")

    def test_a_tool_that_keeps_the_input_numbering_sees_the_real_segments(self):
        got = contig(self._cropped(), renumbers=False)
        assert got.endswith("B400-410, B500-520")

    def test_a_cropped_chain_refuses_to_give_a_constant_offset(self):
        ctx = self._cropped().downstream
        with pytest.raises(ValueError, match="cropped into 2 segments"):
            ctx.offset(True)

    def test_a_context_chain_may_not_carry_hotspots(self):
        with pytest.raises(ValueError, match="context and carries"):
            Target("B", 1, 100, [50], {50: "PHE"}, is_context=True)

    def test_cropping_away_the_epitope_is_refused(self):
        with pytest.raises(ValueError, match="cropped-away part"):
            Target("B", 1, 100, [90], {90: "PHE"}, keep=[(1, 50)])

    def test_the_cif_is_written_cropped(self):
        """BindCraft: delete entire residues, leave no residual parts."""
        atoms = [Pt("A", 117, (0.0, 0.0, 0.0)),
                 Pt("B", 405, (1.0, 0.0, 0.0)),
                 Pt("B", 450, (2.0, 0.0, 0.0))]      # outside the kept runs
        text = write_rfd4_cif(self._cropped(), atoms)
        assert " 405 " in text
        assert " 450 " not in text


class TestWhatAJobsFileContains:
    """The written structure must agree with the job, for every tool.

    It did not. A monospecific job's partner is not in `spec.targets` at
    all, so nothing removed it: the file held the enzyme the mode exists to
    take away and only the contig left it out. A cropped context chain was
    filtered in the mmCIF and nowhere else, so BindCraft -- the tool whose
    memory guidance motivated cropping -- got the whole partner and a
    `memory_note` describing a file that was never written.
    """

    ATOMS: ClassVar[list] = [
        Pt("A", 5, (0.0, 0.0, 0.0)),
        Pt("B", 405, (1.0, 0.0, 0.0)),
        Pt("B", 450, (2.0, 0.0, 0.0)),
        Pt("C", 900, (3.0, 0.0, 0.0), hetatm=True, resname="GST"),
    ]

    def _cropped(self):
        ctx = Target("B", 400, 600, (), {}, "ctx", keep=[(400, 410)],
                     is_context=True)
        return make_spec(downstream=ctx)

    def test_a_dropped_chain_is_not_in_the_file(self):
        mono = design_jobs(make_spec(), "monospecific")[0]
        chains = {a.chain for a in job_atoms(mono, self.ATOMS)}
        assert "B" not in chains

    def test_a_cropped_chain_is_trimmed_in_the_file(self):
        kept = job_atoms(self._cropped(), self.ATOMS)
        got = {(a.chain, a.resseq) for a in kept}
        assert ("B", 405) in got
        assert ("B", 450) not in got

    def test_cargo_on_its_own_chain_survives(self):
        """It is outside every contig anyway; dropping it here would be a
        second, silent decision about the cargo."""
        kept = job_atoms(self._cropped(), self.ATOMS)
        assert any(a.chain == "C" for a in kept)

    def test_an_uncropped_target_keeps_its_cargo_past_the_span(self):
        atoms = [Pt("A", 5, (0.0, 0.0, 0.0)),
                 Pt("A", 999, (1.0, 0.0, 0.0), hetatm=True, resname="GST")]
        kept = job_atoms(make_spec(), atoms)
        assert len(kept) == 2

    def test_bindcraft_is_pointed_at_the_trimmed_file(self):
        s = bindcraft_settings(self._cropped(), "/tmp/job_target.pdb")
        assert s["starting_pdb"] == "/tmp/job_target.pdb"

    def test_bindcraft_falls_back_to_the_complex_when_nothing_was_trimmed(self):
        assert bindcraft_settings(make_spec())["starting_pdb"] \
            == make_spec().complex_pdb

    def test_a_trimmed_job_needs_the_atoms(self, tmp_path):
        mono = design_jobs(make_spec(), "monospecific")[0]
        with pytest.raises(ValueError, match="needs the job's own structure"):
            write_for("BindCraft", mono, str(tmp_path), atoms=None)

    def test_a_trimmed_job_writes_a_target_pdb(self, tmp_path):
        mono = design_jobs(make_spec(), "monospecific")[0]
        made = write_for("BindCraft", mono, str(tmp_path), atoms=self.ATOMS)
        assert any(p.endswith("_target.pdb") for p in made)


class TestTheFrameWarningIsOnlyWhereItApplies:
    def test_a_two_chain_job_carries_it(self):
        assert "FIXED RELATIVE" in write_rfd3(make_spec())

    def test_a_single_chain_job_does_not(self):
        """Printing it anyway is noise that trains a reader to skip it where
        it matters -- and makes a one-chain file read as deliberately two."""
        mono = design_jobs(make_spec(), "monospecific")[0]
        assert "FIXED RELATIVE" not in write_rfd3(mono)
        assert "FIXED RELATIVE" not in write_bindcraft(mono)


class TestAHotspotWithoutATipIsRefused:
    def test_an_alanine_hotspot_fails_here_not_at_load_time(self):
        up = Target("A", 1, 100, [10], {10: "ALA"})
        with pytest.raises(ValueError, match="no TIP definition"):
            rfd3_spec(make_spec(upstream=up))

    def test_the_message_points_at_the_selection_rule(self):
        up = Target("A", 1, 100, [10], {10: "GLY"})
        with pytest.raises(ValueError, match="pick_hotspots"):
            rfd3_spec(make_spec(upstream=up))


class TestJobNamesAreDistinct:
    def test_labels_sharing_a_first_word_still_give_distinct_names(self):
        up = Target("A", 1, 100, [10], {10: "PHE"}, "enzyme one")
        dn = Target("B", 1, 100, [20], {20: "LEU"}, "enzyme two")
        jobs = design_jobs(make_spec(upstream=up, downstream=dn),
                           "monospecific")
        assert len({j.name for j in jobs}) == 2

    def test_a_label_with_no_word_characters_falls_back_to_the_chain(self):
        up = Target("A", 1, 100, [10], {10: "PHE"}, "   ")
        dn = Target("B", 1, 100, [20], {20: "LEU"}, "   ")
        jobs = design_jobs(make_spec(upstream=up, downstream=dn),
                           "monospecific")
        assert len({j.name for j in jobs}) == 2


class TestGappedChainNumbering:
    """A chain with unmodelled residues does not shift by a constant.

    The same class of error as the 184-residue offset, reached from the
    other direction: a tool that renumbers the PRESERVED residues 1..N
    counts the missing ones as absent, and an offset built from the span
    counts them as present.
    """

    def _gapped(self):
        names = {i: "PHE" for i in list(range(1, 50)) + list(range(61, 101))}
        return Target("A", 1, 100, [70], names, "up", present=sorted(names))

    def test_n_residues_counts_what_is_there(self):
        assert self._gapped().n_residues() == 89

    def test_the_output_index_skips_the_gap(self):
        assert self._gapped().output_index(70, True) == 59

    def test_the_input_numbering_is_unchanged(self):
        assert self._gapped().output_index(70, False) == 70

    def test_a_constant_offset_is_refused(self):
        with pytest.raises(ValueError, match="unmodelled gaps"):
            self._gapped().offset(True)

    def test_a_dense_chain_still_shifts_by_a_constant(self):
        t = Target("A", 185, 385, [209], {i: "PHE" for i in range(185, 386)})
        assert t.offset(True) == 184
        assert t.output_index(209, True) == 25

    def test_a_residue_that_is_not_preserved_has_no_index(self):
        with pytest.raises(ValueError, match="not among the"):
            self._gapped().output_index(55, True)
