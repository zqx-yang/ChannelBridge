"""Stage 1: turn a deposited structure into one this pipeline can measure.

Small, deliberately. Preparation for channel work is not the same as
preparation for docking or simulation: no protonation, no minimisation, no
rotamer fixing. What matters here is which atoms are present, because every
later stage measures distances to them -- a water left in the channel is a
wall, and a ligand left in it blocks the route whose geometry you are trying
to measure.

Four decisions, and each is a flag:

**Which chains.** A crystal entry often holds several copies. The tunnel,
the residue classification and the pose all refer to one protomer, so pick
it; the others come back later as symmetry copies in the assembly views.

**Which ligands to keep.** Keep the ones that are part of the protein for
this purpose -- a prosthetic group, a structural metal -- and drop the rest.
The CHANNELLED species is a special case: keep it in the file, because the
start point is derived from where it sits, and exclude it by name in the
CAVER configuration so it does not block its own channel.

**Waters: dropped.** Always, unless you ask otherwise. A water in a channel
is indistinguishable from protein to a sphere-rolling algorithm.

**Alternate locations.** Keep the highest-occupancy one. Leaving both puts
two copies of a side chain in the structure, which narrows a channel that is
not actually narrow.

**Coordinates are never moved.** Not re-centred, not symmetry-expanded, not
superposed. If you took a tunnel from ChannelsDB it is in the deposited
entry's frame, and moving the protein silently invalidates it -- the result
looks like a 20 A disagreement rather than like the mistake it is. Expanding
a biological assembly is a separate step, after measuring.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..io.pdb import AtomArray, read_pdb

__all__ = ["PrepParams", "PrepReport", "prepare", "write_pdb"]

#: Dropped unless kept explicitly. Waters, cryoprotectants, buffer and the
#: other furniture of crystallisation: present because of how the structure
#: was obtained, not because of what the protein is.
FURNITURE = frozenset("""
HOH WAT DOD
GOL EDO PEG PG4 PGE 1PE 2PE P6G MPD DMS
SO4 PO4 ACT ACY FMT CIT TRS EPE MES BTB
IOD BR CL NA K CS RB
""".split())


@dataclass
class PrepParams:
    """What to keep. Everything here is a flag on `chbridge prep`."""

    #: Chains to keep. Empty keeps the first polymer chain found, and says so.
    chains: Sequence[str] = ()
    #: Ligand residue names to KEEP (cofactors, prosthetic groups, the
    #: channelled species). Everything else het is dropped.
    keep_ligands: Sequence[str] = ()
    #: Extra residue names to drop even if they look like ligands.
    drop_ligands: Sequence[str] = ()
    #: Keep waters. Off, and rarely right: a water in a channel is a wall.
    keep_waters: bool = False
    #: Keep every altloc rather than the highest-occupancy one.
    keep_altlocs: bool = False
    #: Rename the kept chain to this id. Needed more often than it looks: a
    #: pair built from two single-chain entries has both enzymes on chain A,
    #: and the pose stage cannot tell them apart. Renaming a chain is not
    #: moving it, so this is safe on the ChannelsDB branch -- a database
    #: sphere chain carries no chain id, only coordinates.
    #:
    #: Only one chain may be kept when renaming: collapsing several chains
    #: onto one id would merge two protomers into one object.
    relabel: Optional[str] = None


@dataclass
class PrepReport:
    """What prepare() did, in enough detail to argue with."""

    chains_kept: List[str] = field(default_factory=list)
    chains_dropped: List[str] = field(default_factory=list)
    relabelled_to: Optional[str] = None
    ligands_kept: Dict[str, int] = field(default_factory=dict)
    ligands_dropped: Dict[str, int] = field(default_factory=dict)
    waters_dropped: int = 0
    altlocs_resolved: int = 0
    atoms_in: int = 0
    atoms_out: int = 0
    notes: List[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            "chains   kept %s%s%s" % (
                "+".join(self.chains_kept) or "(none)",
                ", dropped %s" % "+".join(self.chains_dropped)
                if self.chains_dropped else "",
                ", renamed to %s" % self.relabelled_to
                if self.relabelled_to else ""),
        ]
        if self.ligands_kept:
            lines.append("ligands  kept %s" % ", ".join(
                "%s x%d" % (k, v) for k, v in sorted(self.ligands_kept.items())))
        if self.ligands_dropped:
            lines.append("         dropped %s" % ", ".join(
                "%s x%d" % (k, v) for k, v in sorted(self.ligands_dropped.items())))
        if self.waters_dropped:
            lines.append("waters   dropped %d atom(s)" % self.waters_dropped)
        if self.altlocs_resolved:
            lines.append("altlocs  %d atom(s) dropped, highest occupancy kept"
                         % self.altlocs_resolved)
        lines.append("atoms    %d -> %d" % (self.atoms_in, self.atoms_out))
        lines.extend("note     " + n for n in self.notes)
        return "\n".join(lines)


def prepare(structure, params: Optional[PrepParams] = None) -> Tuple[AtomArray, PrepReport]:
    """Select atoms from a structure. Nothing is moved and nothing is added."""
    params = params or PrepParams()
    atoms = read_pdb(structure) if isinstance(structure, str) else structure
    rep = PrepReport(atoms_in=len(atoms.atoms))

    polymer_chains = []
    for a in atoms.atoms:
        if not a.hetatm and a.chain not in polymer_chains:
            polymer_chains.append(a.chain)
    want = [c for c in params.chains if c]
    if not want:
        want = polymer_chains[:1]
        if len(polymer_chains) > 1:
            rep.notes.append(
                "no --chain given; kept %r of %s. The tunnel, the residue "
                "classes and the pose all refer to ONE protomer, so this is a "
                "choice worth making explicitly"
                % (want[0], "+".join(polymer_chains)))
    missing = [c for c in want if c not in polymer_chains]
    if missing:
        raise ValueError(
            "chain(s) %s are not in %s (it has %s)"
            % (", ".join(missing), getattr(atoms, "source", "the structure"),
               "+".join(polymer_chains) or "none"))
    rep.chains_kept = list(want)
    rep.chains_dropped = [c for c in polymer_chains if c not in want]

    keep = {r.upper() for r in params.keep_ligands}
    drop = {r.upper() for r in params.drop_ligands}

    kept = []
    for a in atoms.atoms:
        if a.chain not in want:
            continue
        resn = (a.resname or "").upper()
        if a.hetatm:
            if resn in drop:
                rep.ligands_dropped[resn] = rep.ligands_dropped.get(resn, 0) + 1
                continue
            if resn in ("HOH", "WAT", "DOD") and not params.keep_waters:
                rep.waters_dropped += 1
                continue
            if resn not in keep and resn in FURNITURE:
                rep.ligands_dropped[resn] = rep.ligands_dropped.get(resn, 0) + 1
                continue
            if resn not in keep:
                rep.ligands_dropped[resn] = rep.ligands_dropped.get(resn, 0) + 1
                continue
            rep.ligands_kept[resn] = rep.ligands_kept.get(resn, 0) + 1
        kept.append(a)

    if not params.keep_altlocs:
        kept, n = _resolve_altlocs(kept)
        rep.altlocs_resolved = n

    if params.relabel is not None:
        new = str(params.relabel)
        if len(new) != 1:
            raise ValueError(
                "--relabel takes a single chain id, got %r" % params.relabel)
        if len(want) > 1:
            raise ValueError(
                "cannot rename %d chains (%s) to one id %r: that would merge "
                "separate protomers into one chain. Keep one chain, or leave "
                "the ids alone." % (len(want), "+".join(want), new))
        import dataclasses

        kept = [dataclasses.replace(a, chain=new) for a in kept]
        rep.relabelled_to = new

    rep.atoms_out = len(kept)
    asked = keep - set(rep.ligands_kept)
    if asked:
        rep.notes.append(
            "asked to keep %s, which %s not in the kept chain(s)"
            % (", ".join(sorted(asked)), "is" if len(asked) == 1 else "are"))
    if not any(not a.hetatm for a in kept):
        raise ValueError("nothing is left after preparation: no polymer atoms kept")
    return AtomArray(kept), rep


def _resolve_altlocs(atoms):
    """Keep one position per atom: the highest occupancy, first on a tie.

    Two copies of a side chain narrow a channel that is not narrow, and a
    sphere-rolling algorithm cannot tell they are alternatives.
    """
    best: Dict[Tuple, object] = {}
    order: List[Tuple] = []
    plain = []
    for a in atoms:
        if not (a.altloc or "").strip():
            plain.append(a)
            continue
        key = (a.chain, a.resseq, a.icode, a.name)
        if key not in best:
            best[key] = a
            order.append(key)
        elif float(a.occupancy or 0) > float(best[key].occupancy or 0):
            best[key] = a
    chosen = list(best[k] for k in order)
    dropped = sum(1 for a in atoms if (a.altloc or "").strip()) - len(chosen)
    out = plain + chosen
    out.sort(key=lambda a: a.serial)
    return out, dropped


_FMT = ("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f"
        "          %2s")


def write_pdb(path: str, atoms: AtomArray, header: Sequence[str] = ()) -> str:
    """Write a PDB, LF-only, coordinates untouched.

    LF because a CR sneaks into a field that a downstream tool then fails to
    parse, and the failure is a long way from the cause.
    """
    lines = ["REMARK   1 %s" % h for h in header]
    for a in atoms.atoms:
        name = a.name if len(a.name) >= 4 else " %-3s" % a.name
        lines.append(_FMT % (
            "HETATM" if a.hetatm else "ATOM  ", a.serial, name,
            a.altloc or " ", a.resname, a.chain, a.resseq, a.icode or " ",
            a.xyz[0], a.xyz[1], a.xyz[2],
            a.occupancy if a.occupancy is not None else 1.0,
            a.bfactor if a.bfactor is not None else 0.0,
            (a.element or "").rjust(2)))
    lines.append("END")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return path
