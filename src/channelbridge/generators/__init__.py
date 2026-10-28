"""Stage 7: the input a binder generator needs in order to hold a chosen pose.

Stages 1-6 end with a *geometry*: two enzymes placed so a channelled
intermediate can pass between them, and a patch on each surface that a binder
should grip. Stage 7 turns that into the input file a design tool actually
reads. It writes no structures of its own beyond a copy of the posed complex:
every number in it comes from the pose table and the patch tables.

**The one thing that is easy to get wrong, and silently.** A binder design
tool is normally pointed at *one* target, and the conditioning vocabulary of
every tool reflects that: it keeps the target rigid *internally* and lets the
whole thing float, because where the target sits in the box does not matter
when there is only one. Here there are two, and the entire result of stages
2-5 is their RELATIVE placement. A spec that makes each enzyme separately
rigid, and the pair free, throws away D_IC and O_IC while looking completely
reasonable -- the run succeeds, the designs score well, and they hold a
geometry nobody chose. So:

* for RFD4-Proteina we condition with ``C_CRD`` (per-atom coordinate, pinned
  where the input puts it) on **both** chains, not ``C_CRDGRP`` per chain.
  ``C_CRDGRP`` preserves a group's internal geometry and leaves its global
  pose to the model, which is right for one target and wrong for two.
* for RFdiffusion3 both target spans go in one ``contig`` with the binder
  segment first, which is what fixes them in the input frame.
* for BindCraft both chains go in ``chains`` -- it treats the starting PDB as
  given, so the pose survives by construction.

:func:`frame_warning` states this in one line and every writer puts it in the
file it writes, because a spec that is copied into another campaign loses this
module but keeps its comments.

**Which tools need a file at all.** The user's rule: a generator that takes
its whole specification as command-line arguments inside the SLURM script does
not need stage 7 to write anything -- it needs the posed structure and the
argument text. :data:`TOOLS` records which is which, with what we checked:

=================  =============  ====================================
tool               spec file      what stage 7 emits
=================  =============  ====================================
RFdiffusion3       JSON or YAML   ``<name>.json``
RFD4-Proteina      YAML or JSON   ``<name>.yaml`` + ``<name>.cif``
BindCraft          JSON           ``<name>_settings.json``
protpardelle-1c    YAML           refused: see :data:`TOOLS`
flexcraft          none           ``<name>_flexcraft_args.txt``
RFdiffusion        none           ``<name>_rfdiffusion_args.txt``
RFDpoly            none           refused: see :data:`TOOLS`
=================  =============  ====================================

The argument files are not inputs to the tool. They are the lines to paste
into the sample job script, which is where those tools take their spec.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field, replace
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: The functional end of each side chain: the atoms a hotspot condition should
#: name. This is RFdiffusion3's ``TIP`` convention written out, because the
#: tools that take explicit atom lists (RFD4-Proteina's ``sel()``) have no
#: ``TIP`` keyword and would otherwise get every atom of the residue, which
#: asks the binder to contact the backbone as well as the side chain.
#:
#: Gly has no side chain and Ala's is one carbon, so they are CA and CB; Pro's
#: ring is its only surface. RFdiffusion3 itself rejects ``TIP`` for exactly
#: those three, which is the check on this table.
TIP_ATOMS: Dict[str, Tuple[str, ...]] = {
    "ALA": ("CB",),
    "ARG": ("CZ", "NH1", "NH2"),
    "ASN": ("CG", "OD1", "ND2"),
    "ASP": ("CG", "OD1", "OD2"),
    "CYS": ("SG",),
    "GLN": ("CD", "OE1", "NE2"),
    "GLU": ("CD", "OE1", "OE2"),
    "GLY": ("CA",),
    "HIS": ("ND1", "CE1", "NE2"),
    "ILE": ("CD1",),
    "LEU": ("CD1", "CD2"),
    "LYS": ("CE", "NZ"),
    "MET": ("SD", "CE"),
    "PHE": ("CZ", "CE1", "CE2"),
    "PRO": ("CG", "CD"),
    "SER": ("OG",),
    "THR": ("OG1", "CG2"),
    "TRP": ("CZ2", "CH2"),
    "TYR": ("CZ", "OH"),
    "VAL": ("CG1", "CG2"),
}

#: Residue types for which RFdiffusion3's ``TIP`` keyword is undefined. A
#: hotspot on one of these fails the whole batch at load time, so stage 7
#: refuses it here with a message naming the residue -- `pick_hotspots`
#: never selects one, and a hand-edited patch table or an API-built
#: :class:`Target` is the only way to reach it.
NO_TIP = ("ALA", "GLY", "PRO")

#: RFdiffusion3's atom-selection keyword for a hotspot residue. ``TIP`` fixes
#: that residue's tip atoms using RFd3's own per-residue-type definition;
#: ``ALL`` would fix every atom, asking the binder to contact the backbone,
#: and ``BKBN`` only the backbone. Always ``TIP``, so we are not keeping a
#: copy of their table in sync with theirs.
HOTSPOT_ATOM_SEL = "TIP"

#: Binder length range, residues. Two enzymes have to be bridged, so the short
#: class that a single-target campaign would use is not offered: the default is
#: the long class of the house convention. It is a DEFAULT, not a measurement
#: -- nothing here has measured how long a binder must be to span a given
#: patch separation, and the pose table's ``centroid_dist`` is the quantity
#: that would inform it.
DEFAULT_BINDER_LENGTH = (100, 150)


@dataclass(frozen=True)
class ToolSpec:
    """What one generator takes, and how we know."""

    #: Name as the cluster's module and sample script spell it.
    name: str
    #: ``True`` if the tool reads a structured specification file that stage 7
    #: must write; ``False`` if its whole specification is command-line
    #: arguments inside the job script.
    spec_file: bool
    #: ``json``, ``yaml``, ``args`` (text to paste into the job script), or
    #: ``""`` when the tool is refused.
    fmt: str
    #: Why, in one sentence. For a refusal, why it cannot be used here.
    note: str
    #: ``False`` when stage 7 will not emit for this tool at all.
    supported: bool = True


#: Every generator we checked, and what it takes. Checked against the sample
#: job scripts under ``/home/programs/scripts/<tool>/`` on this cluster and
#: against each tool's own documentation; the two disagreed nowhere.
TOOLS: Dict[str, ToolSpec] = {
    "RFdiffusion3": ToolSpec(
        "RFdiffusion3", True, "json",
        "INPUT_JSON in the sample script: a .json/.yaml design spec, one "
        "entry per design task, whose file stem becomes the design-name "
        "prefix."),
    "RFD4-Proteina": ToolSpec(
        "RFD4-Proteina", True, "yaml",
        "INPUT_PATH in the sample script: a YAML/JSON condition spec. Its "
        "`input` must be .cif/.cif.gz/.bcif, so stage 7 writes the posed "
        "complex as mmCIF beside the spec."),
    "BindCraft": ToolSpec(
        "BindCraft", True, "json",
        "SETTINGS_JSON in the sample script: target settings (starting_pdb, "
        "chains, target_hotspot_residues, lengths, number_of_final_designs). "
        "design_path is NOT written -- the job wrapper sets it."),
    "flexcraft": ToolSpec(
        "flexcraft", False, "args",
        "No spec file: the target is TARGET_PDB and everything else is "
        "FLEXCRAFT_ARGS in the job script (--num_aa, --hotspots "
        "\"A952,A953\", --radius ...). Stage 7 emits those lines."),
    "RFdiffusion": ToolSpec(
        "RFdiffusion", False, "args",
        "No spec file: Hydra key=value arguments in RFD_ARGS "
        "(contigmap.contigs, ppi.hotspot_res, inference.input_pdb). Stage 7 "
        "emits those lines."),
    "protpardelle-1c": ToolSpec(
        "protpardelle-1c", True, "", supported=False,
        note="Refused. Its shipped sampling configs are unconditional, "
             "partial diffusion, motif scaffolding and structure prediction; "
             "none is two-chain binder design against a held target, and "
             "--motif-pdb conditions on a motif to REBUILD around, not on a "
             "target to bind. Writing a spec would be inventing a mode the "
             "tool does not have."),
    "RFDpoly": ToolSpec(
        "RFDpoly", False, "", supported=False,
        note="Refused. It is the multi-polymer (DNA/RNA/protein) model; its "
             "contigs name polymer chain types and the task here is a protein "
             "binder against two protein chains, which RFdiffusion3 and "
             "RFD4-Proteina already cover."),
}


@dataclass
class Target:
    """One enzyme as the posed complex holds it.

    ``first``/``last`` are the residue numbers IN THE POSED COMPLEX FILE, and
    ``hotspots`` are in that same numbering. Nothing here renumbers anything;
    the per-tool offset (:meth:`offset`) is computed from ``first`` and
    asserted against ``resnames``.
    """

    chain: str
    first: int
    last: int
    hotspots: Sequence[int]
    #: residue number -> three-letter name, for every hotspot at least. Used
    #: to pick tip atoms and to assert that a renumbered selection still
    #: lands on the residue we meant.
    resnames: Mapping[int, str] = field(default_factory=dict)
    #: free-text label for comments, e.g. "upstream, GPPS, patch 2"
    label: str = ""
    #: Residue runs of this chain that the input keeps, ``[(lo, hi), ...]``.
    #: ``None`` means the whole ``first..last`` span. A cropped CONTEXT chain
    #: has several: BindCraft's guidance is to trim a target "only to the
    #: parts necessary to facilitate correct binding", and the memory it
    #: saves is the difference between a job that fits on a card and one that
    #: does not.
    keep: Optional[Sequence[Tuple[int, int]]] = None
    #: True when this chain is present so the binder sees it, but is NOT
    #: being designed against: no hotspots, and nothing asks the binder to
    #: contact it. It is there to occupy the space it occupies in the pose.
    is_context: bool = False
    #: The residue numbers this chain ACTUALLY has, sorted. ``None`` means
    #: "assume ``first..last`` with no gaps", which is only safe for a chain
    #: you know is dense. A crystal structure with an unmodelled loop is not,
    #: and a tool that renumbers the preserved residues 1..N then counts the
    #: missing ones as present: measured on a 1-100 chain missing 50-60, a
    #: hotspot at 70 was emitted as 70 where the tool calls it 59. That is
    #: the same class of error as the 184-residue offset in :meth:`offset`,
    #: arrived at from the other direction.
    present: Optional[Sequence[int]] = None

    def __post_init__(self) -> None:
        if self.is_context and self.hotspots:
            raise ValueError(
                "chain %s is context and carries %d hotspot(s). A context "
                "chain is there to be avoided, not bound; a hotspot on it "
                "asks the binder for a second interface, which is a "
                "different design job." % (self.chain, len(self.hotspots)))
        if self.keep is not None:
            runs = [(int(a), int(b)) for a, b in self.keep]
            if not runs:
                raise ValueError("chain %s: keep is empty; pass None for the "
                                 "whole chain" % self.chain)
            for a, b in runs:
                if b < a:
                    raise ValueError("chain %s: kept run %d-%d is reversed"
                                     % (self.chain, a, b))
            object.__setattr__(self, "keep", sorted(runs))
            out = [r for r in self.hotspots
                   if not any(a <= r <= b for a, b in runs)]
            if out:
                raise ValueError(
                    "chain %s: hotspot(s) %s fall in a cropped-away part of "
                    "the chain (kept: %s). Cropping that removes the epitope "
                    "is not a smaller job, it is a different one."
                    % (self.chain, ",".join(str(r) for r in out),
                       ", ".join("%d-%d" % r for r in runs)))
        if self.last < self.first:
            raise ValueError("%s: last (%d) is before first (%d)"
                             % (self.chain, self.last, self.first))
        bad = [r for r in self.hotspots
               if not self.first <= r <= self.last]
        if bad:
            raise ValueError(
                "chain %s: hotspot(s) %s are outside the chain's span %d-%d. "
                "A hotspot the contig does not preserve is silently dropped "
                "by every tool here, so this is an error, not a warning."
                % (self.chain, ",".join(str(r) for r in bad),
                   self.first, self.last))
        missing = [r for r in self.hotspots if r not in self.resnames]
        if self.resnames and missing:
            raise ValueError(
                "chain %s: no residue name for hotspot(s) %s, so their tip "
                "atoms cannot be named"
                % (self.chain, ",".join(str(r) for r in missing)))

    def runs(self) -> List[Tuple[int, int]]:
        """The residue runs this chain contributes, in order."""
        return list(self.keep) if self.keep else [(self.first, self.last)]

    def preserved(self) -> List[int]:
        """Every residue number the input hands the tool, in order.

        The kept runs, intersected with the residues that actually exist.
        This -- not the span -- is what a tool renumbering from 1 counts.
        """
        runs = self.runs()
        if self.present is not None:
            have = sorted(int(r) for r in self.present)
        else:
            have = list(range(self.first, self.last + 1))
        return [r for r in have if any(a <= r <= b for a, b in runs)]

    def n_residues(self) -> int:
        """How many residues of this chain reach the generator."""
        return len(self.preserved())

    def output_index(self, resi: int, renumbers: bool) -> int:
        """Where *resi* lands in the tool's numbering.

        The one place a residue number is converted, so the conversion can
        be right for a dense chain, a chain with unmodelled gaps and a
        cropped chain at once. :meth:`offset` is the special case of this
        for a dense, uncropped chain, and refuses the others.
        """
        if not renumbers:
            return int(resi)
        order = self.preserved()
        try:
            return order.index(int(resi)) + 1
        except ValueError:
            raise ValueError(
                "chain %s: residue %d is not among the %d residues this "
                "input preserves, so it has no position in the tool's "
                "numbering" % (self.chain, resi, len(order))) from None

    def output_runs(self, renumbers: bool) -> List[Tuple[int, int]]:
        """The kept runs as the TOOL will number them.

        With ``renumbers`` false they are the input's own numbers. With it
        true the preserved residues of this chain are counted 1..N in order
        -- "Output res_id restarts at 1 per chain" -- so neither a crop nor
        an unmodelled gap shifts them by a constant: the whole chain becomes
        one span of however many residues it actually contributes.
        """
        if not renumbers:
            return self.runs()
        return [(1, self.n_residues())]

    def offset(self, renumbers: bool) -> int:
        """How much to subtract from a residue number for one tool.

        A tool that renumbers each preserved chain from 1 (RFD4-Proteina does:
        "Output res_id restarts at 1 per chain") needs ``first - 1`` taken
        off; a tool that keeps the input numbering needs 0. This was found the
        hard way on another campaign, where a hotspot set 184 residues off
        still produced designs -- against the wrong epitope.
        """
        if renumbers:
            if self.keep is not None and len(self.runs()) > 1:
                raise ValueError(
                    "chain %s is cropped into %d segments, so a tool that "
                    "renumbers from 1 does not shift its residues by a "
                    "constant and this offset would be wrong for every "
                    "residue after the first gap. Use output_index()."
                    % (self.chain, len(self.runs())))
            if self.n_residues() != self.last - self.first + 1:
                raise ValueError(
                    "chain %s spans %d-%d but has only %d residues, so it "
                    "has unmodelled gaps and does not shift by a constant "
                    "under renumbering. Use output_index()."
                    % (self.chain, self.first, self.last, self.n_residues()))
        return (self.first - 1) if renumbers else 0

    def tip_atoms(self, resi: int) -> Tuple[str, ...]:
        name = (self.resnames.get(resi) or "").upper()
        if name not in TIP_ATOMS:
            raise ValueError(
                "chain %s residue %d is %r, which has no tip-atom definition. "
                "Add it to TIP_ATOMS rather than letting the whole residue be "
                "selected: a hotspot on every atom asks the binder to contact "
                "the backbone too." % (self.chain, resi, name or "?"))
        return TIP_ATOMS[name]


@dataclass
class PoseSpec:
    """A chosen pose, ready to be written out for any generator."""

    #: prefix for everything written, per the house convention
    name: str
    #: the posed complex as stage 5 wrote it (PDB)
    complex_pdb: str
    upstream: Target
    #: ``None`` only in a single-target job -- see :func:`design_jobs`.
    downstream: Optional[Target]
    #: pose id from the pose table, e.g. ``geom001``
    geom_id: str = ""
    d_ic: Optional[float] = None
    o_ic: Optional[float] = None
    do: Optional[float] = None
    binder_length: Tuple[int, int] = DEFAULT_BINDER_LENGTH
    #: how many designs to ask for, where the tool takes a count
    n_designs: int = 100

    def __post_init__(self) -> None:
        if self.downstream is not None and \
                self.upstream.chain == self.downstream.chain:
            raise ValueError(
                "both enzymes are on chain %s. The pose stage writes upstream "
                "and downstream on different chains (CHAIN_IDS upstream ACEG, "
                "downstream BDFH) precisely so a generator can tell them "
                "apart; one chain means the prep relabelling was skipped."
                % self.upstream.chain)
        if not self.designed:
            raise ValueError(
                "every chain in this job is context, so nothing is being "
                "designed against. At least one target needs hotspots.")
        lo, hi = self.binder_length
        if not 0 < lo <= hi:
            raise ValueError("binder_length must be 0 < lo <= hi, got %r"
                             % (self.binder_length,))

    @property
    def targets(self) -> Tuple[Target, ...]:
        """Every chain in the input, designed or context, in chain order."""
        if self.downstream is None:
            return (self.upstream,)
        return (self.upstream, self.downstream)

    @property
    def designed(self) -> Tuple[Target, ...]:
        """The chains the binder is being designed against."""
        return tuple(t for t in self.targets if not t.is_context)

    @property
    def context(self) -> Tuple[Target, ...]:
        """The chains that are present but not targeted."""
        return tuple(t for t in self.targets if t.is_context)

    @property
    def mode(self) -> str:
        if len(self.designed) > 1:
            return "bispecific"
        return "context" if self.context else "monospecific"

    def n_residues(self) -> int:
        """Target residues reaching the generator, binder excluded.

        The number BindCraft's memory guidance is about: roughly 550
        residues of target plus binder on a 32 GB card, 950 on an 80 GB one.
        """
        return sum(t.n_residues() for t in self.targets)

    def provenance(self) -> str:
        """One line naming the pose these numbers came from."""
        bits = [self.name]
        if self.geom_id:
            bits.append(self.geom_id)
        for k, v in (("D_IC", self.d_ic), ("O_IC", self.o_ic), ("DO", self.do)):
            if v is not None:
                bits.append("%s %.4g" % (k, v))
        return "  ".join(bits)


def frame_warning(spec: Optional["PoseSpec"] = None) -> str:
    """The sentence a generated spec carries about the inter-chain frame.

    Empty for a single-chain job: there is no relative placement to preserve
    and printing the warning anyway would be noise that trains a reader to
    skip it where it matters.
    """
    if spec is not None and len(spec.targets) < 2:
        return ""
    return ("the two target chains must stay FIXED RELATIVE TO EACH OTHER: "
            "their relative placement IS the result of stages 2-5, and a "
            "spec that only keeps each chain rigid on its own discards D_IC "
            "and O_IC without failing")


def tools_table() -> str:
    """The survey, as a plain-text table. ``chbridge generate --list-tools``."""
    rows = [("tool", "spec file", "format", "stage 7")]
    for t in TOOLS.values():
        rows.append((t.name, "yes" if t.spec_file else "no (job script)",
                     t.fmt or "-", "writes it" if t.supported else "REFUSED"))
    w = [max(len(r[i]) for r in rows) for i in range(4)]
    out = []
    for i, r in enumerate(rows):
        out.append("  ".join(r[j].ljust(w[j]) for j in range(4)).rstrip())
        if i == 0:
            out.append("  ".join("-" * w[j] for j in range(4)))
    out.append("")
    for t in TOOLS.values():
        out.append("%s: %s" % (t.name, t.note))
    return "\n".join(out)


# --------------------------------------------------------------------------
# contigs
# --------------------------------------------------------------------------

def contig(spec: PoseSpec, *, renumbers: bool, join: str = ", ") -> str:
    """``<lo>-<hi>, /0, A<a>-<b>, /0, B<c>-<d>``.

    The binder segment comes first and each target span follows after its own
    chain break. Both tools read the same grammar here: a bare range is a
    designed segment of sampled length, ``/0`` starts a new chain, and
    ``<chain><lo>-<hi>`` preserves that span of the input.

    ``renumbers`` says whether the tool renumbers preserved chains from 1. It
    changes the SPAN WRITTEN IN THE CONTIG as well as the hotspot numbers, and
    getting only one of the two right is the failure this argument exists to
    prevent.
    """
    lo, hi = spec.binder_length
    parts = ["%d-%d" % (lo, hi)]
    for t in spec.targets:
        parts.append("/0")
        # A cropped chain contributes several preserved segments. They are
        # listed without a chain break between them: `/0` starts a new CHAIN,
        # and these are discontinuous pieces of one.
        parts += ["%s%d-%d" % (t.chain, a, b)
                  for a, b in t.output_runs(renumbers)]
    return join.join(parts)


# --------------------------------------------------------------------------
# design modes
# --------------------------------------------------------------------------

#: The three ways to ask for a binder against a posed pair.
MODES = ("bispecific", "monospecific", "context")

#: How far from the epitope a context residue is kept, in angstrom.
#:
#: Measured, not inherited: over the 38 enzyme sides of ``examples/iMARS``
#: the closest partner atom sits 8.7 to 34.6 A from the epitope (median
#: 15.2). The worked example's 20 A reaches the partner on 28 of those 38;
#: 30 reaches it on 37 and 40 on all of them, while keeping more of a
#: partner that the memory budget wants gone. 30 is the compromise, and the
#: quantity it has to cover is a binder's reach from the epitope -- roughly
#: twice the ~9 A binder radius the pose stage uses.
CONTEXT_CUT = 30.0

#: BindCraft's memory guidance, target + binder residues per card. Quoted
#: rather than derived: "a 32 GB card fits about 550 residues, an 80 GB one
#: about 950". Used only to warn.
CARD_RESIDUES = ((550, "32 GB"), (950, "80 GB"))

#: Runs shorter than this are dropped from a crop. An isolated two-residue
#: fragment has no structural core, and BindCraft is explicit that such
#: pieces "will be rearranged randomly".
MIN_CONTEXT_RUN = 4


def crop_runs(atoms: Sequence[object], chain: str, near_xyz: Sequence[Sequence[float]],
              *, cut: float = CONTEXT_CUT,
              min_run: int = MIN_CONTEXT_RUN) -> List[Tuple[int, int]]:
    """Residue runs of *chain* within *cut* of any point in *near_xyz*.

    The crop BindCraft asks for -- "trim the input only to the parts
    necessary to facilitate correct binding" -- expressed as a measurement
    rather than a hand-drawn boundary: keep what is near the epitope the
    binder is being designed against, drop the rest of the partner.

    Runs shorter than *min_run* are dropped. Two caveats, both from
    BindCraft's own guidance, and neither of which this can check for you:
    a trim should fall at a domain boundary or a hinge, and it should not
    expose a hydrophobic core. A distance cut does neither on purpose, so
    read what comes out before trusting it on a new target.
    """
    import numpy as np

    pts = np.asarray(near_xyz, dtype=float).reshape(-1, 3)
    if not len(pts):
        return []
    by_res: Dict[int, List[Sequence[float]]] = {}
    for a in atoms:
        if getattr(a, "chain", None) != chain or getattr(a, "hetatm", False):
            continue
        by_res.setdefault(int(a.resseq), []).append(a.xyz)
    keep = []
    for resi, xyz in by_res.items():
        q = np.asarray(xyz, dtype=float)
        d = np.linalg.norm(q[:, None, :] - pts[None, :, :], axis=2).min()
        if d <= float(cut):
            keep.append(resi)
    keep.sort()
    runs: List[Tuple[int, int]] = []
    for r in keep:
        if runs and r == runs[-1][1] + 1:
            runs[-1] = (runs[-1][0], r)
        else:
            runs.append((r, r))
    return [t for t in runs if t[1] - t[0] + 1 >= int(min_run)]


def memory_note(spec: PoseSpec) -> str:
    """What BindCraft's memory guidance says about a job this size."""
    n = spec.n_residues() + spec.binder_length[1]
    for limit, card in CARD_RESIDUES:
        if n <= limit:
            return ("%d residues (targets %d + binder up to %d): fits "
                    "BindCraft's guidance for a %s card"
                    % (n, spec.n_residues(), spec.binder_length[1], card))
    return ("%d residues (targets %d + binder up to %d): OVER BindCraft's "
            "guidance of about %d for an %s card. Crop the context "
            "(--context-cut) or design against one target at a time."
            % (n, spec.n_residues(), spec.binder_length[1],
               CARD_RESIDUES[-1][0], CARD_RESIDUES[-1][1]))


def design_jobs(spec: PoseSpec, mode: str = "bispecific", *,
                atoms: Optional[Sequence[object]] = None,
                context_cut: float = CONTEXT_CUT) -> List[PoseSpec]:
    """The job(s) one pose becomes, for one design mode.

    * ``bispecific`` -- one job: a binder against BOTH enzymes at once, which
      is the thing the pose is for. Hardest to design, and the only mode that
      can produce the bridge directly.
    * ``monospecific`` -- one job per enzyme: a binder against that enzyme
      alone, the other absent. Easiest, and the pose's geometry is nowhere in
      it: two binders designed this way have no reason to be compatible, and
      joining them is a separate problem.
    * ``context`` -- one job per enzyme: a binder against that enzyme with
      the OTHER present but not targeted, cropped to what sits within
      ``context_cut`` of the epitope. The partner is there to occupy the
      space it occupies in the pose, so a binder is not designed into it.
      Needs *atoms* -- the posed complex -- to do the cropping.

    The three are a ladder, not alternatives: a campaign can take the
    monospecific answer as a starting point, re-run it in context, and keep
    the bispecific job for the designs that survive.
    """
    if mode not in MODES:
        raise ValueError("unknown mode %r; known: %s"
                         % (mode, ", ".join(MODES)))
    if mode == "bispecific":
        return [spec]
    if spec.downstream is None:
        raise ValueError("%s mode splits a two-target pose; this one has "
                         "a single target already" % mode)
    if mode == "context" and atoms is None:
        raise ValueError(
            "context mode crops the partner chain to what is near the "
            "epitope, which needs the posed complex's atoms")

    out: List[PoseSpec] = []
    pairs = ((spec.upstream, spec.downstream), (spec.downstream, spec.upstream))
    for designed, partner in pairs:
        other: Optional[Target] = None
        if mode == "context":
            anchor = [a.xyz for a in (atoms or ())
                      if getattr(a, "chain", None) == designed.chain
                      and int(getattr(a, "resseq", -1)) in set(designed.hotspots)]
            runs = crop_runs(atoms or (), partner.chain, anchor,
                             cut=context_cut)
            if not runs:
                # Not an error: it is a measurement. No part of the partner
                # is near this epitope, so there is no wall for the binder to
                # avoid and context mode has nothing to add -- the job IS the
                # monospecific one. Failing here would stop the other side's
                # job too, over a fact about this one.
                other = None
            else:
                other = replace(partner, hotspots=(), is_context=True,
                                keep=runs)
        # The tag has to survive being a file name AND be unique. A label
        # is free text ("upstream, GPPS, patch 2"), so it is sanitised and
        # the chain id appended: two labels sharing a first word would
        # otherwise give two jobs one name, and `write_for` would overwrite
        # the first without a word.
        raw = (designed.label or "").strip().split()[:1] or [designed.chain]
        tag = re.sub(r"\W+", "_", raw[0]).strip("_") or designed.chain
        name = "%s_vs_%s_%s%s" % (spec.name, tag, designed.chain,
                                  "_in_context" if mode == "context" else "")
        out.append(replace(spec, name=name, upstream=designed, downstream=other))
    names = [j.name for j in out]
    if len(set(names)) != len(names):
        raise ValueError("design_jobs produced duplicate names (%s); the "
                         "second would overwrite the first" % ", ".join(names))
    return out


def job_atoms(spec: PoseSpec, atoms: Sequence[object]) -> List[object]:
    """The atoms a job's input file should contain, and no others.

    Two filters, and both were missing. A **monospecific** job's partner is
    not in ``spec.targets`` at all, so nothing removed it: the file still
    held the enzyme the mode exists to take away, and only the contig left
    it out. A **cropped** context chain was filtered in the mmCIF and
    nowhere else, so BindCraft -- the tool whose memory guidance motivated
    cropping in the first place -- was handed the whole partner.

    Chains the spec does not describe at all (a cargo written on its own
    chain, say) are kept: this decides what the JOB is, not what the cargo
    is, and the contig already excludes het records.
    """
    # Only a CROPPED chain is trimmed by residue. An uncropped target is
    # kept whole, cargo and all: a kept cofactor sits past the polymer span
    # by construction, the contig does not preserve it either way, and
    # dropping it here would be a second, silent decision about the cargo.
    cropped = {t.chain: t.runs() for t in spec.targets if t.keep is not None}
    known = {t.chain for t in spec.targets}
    out = []
    for a in atoms:
        ch = getattr(a, "chain", None)
        if ch in cropped:
            r = int(getattr(a, "resseq", 0))
            if not any(lo <= r <= hi for lo, hi in cropped[ch]):
                continue
        elif ch not in known and _is_polymer_chain(ch, atoms, known):
            continue          # a target chain this job does not include
        out.append(a)
    return out


def _is_polymer_chain(chain, atoms, known) -> bool:
    """True if *chain* is a protein chain the spec left out."""
    if chain in known:
        return False
    return any(getattr(a, "chain", None) == chain and not getattr(a, "hetatm", False)
               for a in atoms)


def write_job_pdb(spec: PoseSpec, atoms: Sequence[object], path: str) -> str:
    """The job's structure as a PDB, for the tools that take one.

    BindCraft is explicit that a trimmed target must be trimmed in the FILE
    -- "delete entire residues, and not to leave any residual atoms or
    residue parts" -- and it never moves what it is given, so this file is
    also what preserves the pose.
    """
    from ..io.pdb import write_pdb

    write_pdb(path, job_atoms(spec, atoms))
    return path


# --------------------------------------------------------------------------
# RFdiffusion3
# --------------------------------------------------------------------------

def rfd3_spec(spec: PoseSpec, *, key: Optional[str] = None,
              non_loopy: bool = True) -> Dict[str, object]:
    """The RFdiffusion3 input JSON, as a dict (one key, one design task).

    ``select_hotspots`` maps ``<chain><resi>`` to **``TIP``**, always.

    That is RFdiffusion3's own keyword for "all of this residue's tip atoms",
    resolved against its own per-residue-type table. Writing the atoms out
    instead means keeping a second copy of that table in sync with theirs,
    for no gain -- and the only residues where the two could differ are Ala,
    Gly and Pro, which have no tip and are not hotspots: `pick_hotspots`
    chooses from the exposed F/Y/W/I/L/M of a patch. An earlier version
    handed over whole patches and therefore had to emit `CA` and `CB` for
    those three, which is the table this no longer needs.

    RFdiffusion3 keeps the input numbering, so nothing is renumbered here.
    """
    hot: Dict[str, str] = {}
    for t in spec.designed:
        for r in t.hotspots:
            name = (t.resnames.get(r) or "").upper()
            if name in NO_TIP:
                raise ValueError(
                    "chain %s residue %d is %s, and RFdiffusion3 has no TIP "
                    "definition for it -- the batch would fail at load. %s "
                    "has no side chain a binder can grip; it should not be a "
                    "hotspot. patches.pick_hotspots selects only exposed "
                    "F/Y/W/I/L/M, so this came from a hand-edited patch "
                    "table or a Target built by hand."
                    % (t.chain, r, name, name))
            hot["%s%d" % (t.chain, r)] = HOTSPOT_ATOM_SEL
    return {
        key or spec.name: {
            "dialect": 2,
            "input": spec.complex_pdb,
            "contig": contig(spec, renumbers=False, join=","),
            "infer_ori_strategy": "hotspots",
            "select_hotspots": hot,
            "is_non_loopy": bool(non_loopy),
        }
    }


def write_rfd3(spec: PoseSpec, **kw) -> str:
    """The JSON text. JSON has no comments, so the provenance goes in a key."""
    body = rfd3_spec(spec, **kw)
    out = {
        "_comment": "%s | %s | %s" % (
            spec.provenance(), frame_warning(spec),
            "written by channelbridge.generators (stage 7)"),
    }
    out.update(body)
    return json.dumps(out, indent=2) + "\n"


# --------------------------------------------------------------------------
# RFD4-Proteina
# --------------------------------------------------------------------------

def _sel(chain: str, resi: Optional[int] = None,
         atoms: Sequence[str] = ()) -> str:
    """One atomworks path selection: ``A/*/25/[CZ,OH]``.

    Trailing fields default to ``*``, so a bare chain is just ``A``.
    """
    if resi is None:
        return "sel('%s')" % chain
    if not atoms:
        return "sel('%s/*/%d')" % (chain, resi)
    return "sel('%s/*/%d/[%s]')" % (chain, resi, ",".join(atoms))


def rfd4_conditions(spec: PoseSpec) -> List[Tuple[str, str, List[str]]]:
    """``[(key, value, [select, ...]), ...]`` for the ``conditions`` mapping.

    Three conditions, and the first is the one that matters:

    ``C_CRD: true`` on every atom of both target chains pins them WHERE THE
    INPUT PUTS THEM. The documented alternative, ``C_CRDGRP``, keeps a group's
    internal geometry and leaves its global pose to the model: correct for a
    single target, and here it would free the inter-enzyme placement that is
    the whole result of the pose search. See :func:`frame_warning`.

    ``C_SEQ: true`` gives the model both target sequences, so it designs
    against the real surfaces rather than redesigning them.

    ``C_HOT: true`` on the tip atoms of each patch residue asks for contact
    there. Residue numbers are renumbered because RFD4-Proteina restarts
    ``res_id`` at 1 in each preserved chain and resolves ``conditions`` after
    ``contig``.
    """
    chains = [t.chain for t in spec.targets]
    hot: List[str] = []
    # Only the DESIGNED chains. A context chain has no hotspots by
    # construction, and asking it for an offset would raise where it is
    # cropped -- rightly, since a cropped chain's residues do not shift by a
    # constant. Iterating `designed` says that, instead of relying on the
    # inner loop being empty.
    for t in spec.designed:
        for r in t.hotspots:
            hot.append(_sel(t.chain, t.output_index(r, True), t.tip_atoms(r)))
    return [
        ("C_CRD", "true", [_sel(c) for c in chains]),
        ("C_SEQ", "true", [_sel(c) for c in chains]),
        ("C_HOT", "true", hot),
    ]


def write_rfd4(spec: PoseSpec, *, cif_name: Optional[str] = None,
               key: Optional[str] = None) -> str:
    """The YAML text. ``input`` is resolved relative to the spec file."""
    name = cif_name or (spec.name + ".cif")
    lines = [
        "# %s" % spec.provenance(),
        "# RFD4-Proteina condition spec, written by channelbridge.generators",
        "# (stage 7). `input` is resolved relative to THIS file.",
        "# C_CRD, not C_CRDGRP, on every target chain because %s."
        % frame_warning(spec),
        "%s:" % (key or spec.name),
        "  input: %s" % name,
        '  contig: "%s"' % contig(spec, renumbers=True),
        "  conditions:",
    ]
    for cond, value, sels in rfd4_conditions(spec):
        body = ", ".join('{select: "%s"}' % s for s in sels)
        lines.append("    %s: {%s: [%s]}" % (cond, value, body))
    return "\n".join(lines) + "\n"


def write_rfd4_cif(spec: PoseSpec, atoms: Sequence[object]) -> str:
    """The posed complex as a minimal mmCIF, which is what RFD4 will read.

    RFD4-Proteina refuses a PDB (``input`` must be ``.cif``/``.cif.gz``/
    ``.bcif``), so the complex is rewritten rather than converted by a third
    tool. Only ``_atom_site`` is emitted: nothing downstream of here reads a
    header, and a header we invented would be a header somebody later trusts.

    ``label_seq_id`` is written as the author numbering, NOT as a 1..N
    counter. A reader that prefers ``label_seq_id`` then sees the same numbers
    the patch tables use, and the contig written by :func:`contig` with
    ``renumbers=True`` accounts for the renumbering RFD4 does on top.

    A CROPPED chain is written cropped: only the residues the contig
    preserves are emitted. BindCraft puts the same rule more bluntly for its
    own input -- "delete entire residues, and not to leave any residual atoms
    or residue parts" -- and writing the whole chain and trusting the contig
    to hide the rest would leave the file disagreeing with the spec beside
    it.
    """
    atoms = job_atoms(spec, atoms)
    cols = ("group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
            "label_comp_id", "label_asym_id", "label_entity_id",
            "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y",
            "Cartn_z", "occupancy", "B_iso_or_equiv", "auth_seq_id",
            "auth_asym_id", "pdbx_PDB_model_num")
    ent = {}
    for a in atoms:
        ent.setdefault(a.chain, str(len(ent) + 1))
    out = ["data_%s" % spec.name, "#", "loop_"]
    out += ["_atom_site.%s" % c for c in cols]
    for i, a in enumerate(atoms, 1):
        el = (getattr(a, "element", "") or a.name[0]).strip() or "C"
        out.append(" ".join((
            "HETATM" if a.hetatm else "ATOM", str(i), el, a.name, ".",
            a.resname, a.chain, ent[a.chain], str(a.resseq), "?",
            "%.3f" % a.xyz[0], "%.3f" % a.xyz[1], "%.3f" % a.xyz[2],
            "%.2f" % getattr(a, "occupancy", 1.0),
            "%.2f" % getattr(a, "bfactor", 0.0),
            str(a.resseq), a.chain, "1")))
    return "\n".join(out) + "\n#\n"


# --------------------------------------------------------------------------
# BindCraft
# --------------------------------------------------------------------------

def bindcraft_settings(spec: PoseSpec,
                       starting_pdb: Optional[str] = None) -> Dict[str, object]:
    """BindCraft target settings.

    ``chains`` holds every chain of the job: BindCraft takes ``starting_pdb``
    as given and never moves the target, so the pose survives by construction
    -- this is the one tool here where the frame needs no condition at all.
    A context chain is listed too, and simply carries no hotspot: BindCraft
    will then avoid it rather than bind it.

    ``target_hotspot_residues`` uses BindCraft's chain-qualified spelling
    (``A117,A118,B441``). Residue numbers are the input's: BindCraft does not
    renumber. ``design_path`` is deliberately absent -- the cluster's job
    wrapper writes it into a patched copy, and a value here would be
    overwritten, so writing one would only mislead.
    """
    hot: List[str] = []
    for t in spec.designed:
        hot += ["%s%d" % (t.chain, r) for r in t.hotspots]
    lo, hi = spec.binder_length
    return {
        "binder_name": spec.name,
        # The job's own structure where one was written: a monospecific or
        # cropped job must not be handed the whole posed complex, or the
        # trim that the memory note describes never reaches the tool.
        "starting_pdb": starting_pdb or spec.complex_pdb,
        "chains": ",".join(t.chain for t in spec.targets),
        "target_hotspot_residues": ",".join(hot),
        "lengths": [lo, hi],
        "number_of_final_designs": int(spec.n_designs),
    }


def write_bindcraft(spec: PoseSpec, starting_pdb: Optional[str] = None) -> str:
    out = {"_comment": "%s | %s" % (spec.provenance(), frame_warning(spec))}
    out.update(bindcraft_settings(spec, starting_pdb))
    return json.dumps(out, indent=2) + "\n"


# --------------------------------------------------------------------------
# tools whose whole spec is job-script arguments
# --------------------------------------------------------------------------

def write_flexcraft_args(spec: PoseSpec,
                         starting_pdb: Optional[str] = None) -> str:
    """The lines to paste into ``FC_sample_sh_1gpu.sh``.

    flexcraft takes one ``--hotspots`` string of chain-qualified residues and
    checks it against the target PDB, so a wrong number fails loudly rather
    than quietly -- the only tool here with that property.

    ``--num_aa`` is a single length, not a range, so the midpoint of the
    binder range is written and the range is stated in the comment: a silent
    choice of one end would bias every design one way.
    """
    hot = []
    for t in spec.designed:
        hot += ["%s%d" % (t.chain, r) for r in t.hotspots]
    lo, hi = spec.binder_length
    mid = (lo + hi) // 2
    return "\n".join([
        "# %s" % spec.provenance(),
        "# flexcraft takes no spec file: these are job-script lines.",
        "# %s -- flexcraft never moves TARGET_PDB, so the pose is safe." %
        frame_warning(spec),
        'TARGET_PDB="%s"' % (starting_pdb or spec.complex_pdb),
        "FLEXCRAFT_ARGS=(",
        "    --num_aa %d                  # one length only; our range is %d-%d"
        % (mid, lo, hi),
        '    --hotspots "%s"' % ",".join(hot),
        "    --num_designs %d" % spec.n_designs,
        ")",
    ]) + "\n"


def write_rfdiffusion_args(spec: PoseSpec,
                           starting_pdb: Optional[str] = None) -> str:
    """The lines to paste into ``RFd_sample_sh_1gpu.sh``.

    RFdiffusion 1 has no per-atom hotspot vocabulary: ``ppi.hotspot_res``
    takes whole residues, so the tip-atom choice in :data:`TIP_ATOMS` is not
    expressible and is simply not used. That is a real loss of specificity
    relative to RFdiffusion3 and RFD4-Proteina, and it is why this is emitted
    as a fallback rather than recommended.
    """
    hot = []
    for t in spec.designed:
        hot += ["%s%d" % (t.chain, r) for r in t.hotspots]
    lo, hi = spec.binder_length
    # RFdiffusion keeps the input numbering, so the contig names the kept
    # runs at their own numbers -- a cropped chain is several segments.
    spans = ",".join("%s%d-%d" % (t.chain, a, b)
                     for t in spec.targets
                     for a, b in t.output_runs(False))
    return "\n".join([
        "# %s" % spec.provenance(),
        "# RFdiffusion (v1) takes no spec file: these are job-script lines.",
        "# %s." % frame_warning(spec),
        "# Whole-residue hotspots only: the tip atoms stage 7 picked for",
        "# RFdiffusion3/RFD4-Proteina cannot be expressed here.",
        "RFD_ARGS=(",
        "    'contigmap.contigs=[%d-%d/0 %s]'" % (lo, hi, spans),
        "    'ppi.hotspot_res=[%s]'" % ",".join(hot),
        "    inference.num_designs=%d" % spec.n_designs,
        "    inference.input_pdb=%s" % (starting_pdb or spec.complex_pdb),
        ")",
    ]) + "\n"


#: Tools whose input is a structure FILE they read as given, rather than a
#: spec that names spans. A trimmed job has to hand these a trimmed file.
STRUCTURE_TOOLS = ("BindCraft", "flexcraft", "RFdiffusion")

#: ``tool -> (writer, file suffix)`` for every supported tool.
WRITERS = {
    "RFdiffusion3": (write_rfd3, ".json"),
    "RFD4-Proteina": (write_rfd4, ".yaml"),
    "BindCraft": (write_bindcraft, "_settings.json"),
    "flexcraft": (write_flexcraft_args, "_flexcraft_args.txt"),
    "RFdiffusion": (write_rfdiffusion_args, "_rfdiffusion_args.txt"),
}


def write_for(tool: str, spec: PoseSpec, out_dir: str,
              atoms: Optional[Sequence[object]] = None) -> List[str]:
    """Write every file *tool* needs for *spec*. Returns the paths written."""
    t = TOOLS.get(tool)
    if t is None:
        raise ValueError("unknown tool %r; known: %s"
                         % (tool, ", ".join(sorted(TOOLS))))
    if not t.supported:
        raise ValueError("%s is not supported by stage 7. %s"
                         % (tool, t.note))
    writer, suffix = WRITERS[tool]
    os.makedirs(out_dir, exist_ok=True)
    written: List[str] = []
    # Does this job's input differ from the posed complex? It does whenever a
    # chain was dropped (monospecific) or cropped (context), and then every
    # tool that takes a structure PATH must be given the trimmed file --
    # otherwise the crop exists only in the contig, and BindCraft, which
    # reads no contig at all, sees the whole partner.
    trimmed = (spec.downstream is None
               or any(t.keep is not None for t in spec.targets))
    if atoms is None and (trimmed or tool == "RFD4-Proteina"):
        raise ValueError(
            "%s needs the job's own structure written out (a chain is "
            "dropped or cropped), and so needs the posed complex's atoms; "
            "read %s and pass them in" % (tool, spec.complex_pdb))
    if tool == "RFD4-Proteina":
        cif = os.path.join(out_dir, spec.name + ".cif")
        with open(cif, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(write_rfd4_cif(spec, atoms))
        written.append(cif)
    start_pdb = None
    if trimmed and tool in STRUCTURE_TOOLS:
        start_pdb = write_job_pdb(
            spec, atoms, os.path.join(out_dir, spec.name + "_target.pdb"))
        written.append(start_pdb)
    path = os.path.join(out_dir, spec.name + suffix)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(writer(spec, start_pdb) if tool in STRUCTURE_TOOLS
                 else writer(spec))
    written.append(path)
    return written
