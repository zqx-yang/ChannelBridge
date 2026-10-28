"""Stage 1: turn a deposited structure into one this pipeline can measure.

Reads PDB or mmCIF and writes PDB. This is the only stage that accepts
mmCIF: CAVER cannot read it, so a later stage taking one would show up as
an empty tunnel rather than as a wrong file type. `channelbridge.io.mmcif`
explains the auth_* vs label_* numbering trap that conversion walks into.

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

from ..io import AtomArray, pdb_column_problems, read_structure

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
    #: How close a het residue in ANOTHER chain must come to a kept chain to
    #: count as belonging to it.
    #:
    #: Needed because a predicted or transplanted model does not put a ligand
    #: in its protein's chain. AlphaFill gives every transplant its own chain
    #: id (one model here has the protein in A and ligands in B..Z and AA..AW),
    #: so `--chain A --keep FPP` keeps the protein and none of the ligands, the
    #: report says "ligands kept: none", and the channel then measured runs
    #: through the empty site whose occupant was the point of the exercise.
    #: Contact, not chain id, is what makes a ligand this protein's.
    ligand_contact: float = 4.0
    #: Residue-number ranges to drop, as ``[(first, last)]``.
    #:
    #: For a disordered terminus, which is a wall to a sphere-rolling
    #: algorithm exactly as a water is -- and unlike a water it is attached, so
    #: it cannot be dropped by residue name. GGPPS's model carries an
    #: 18-residue N-terminal tail that swings 9 A across AlphaFold3's five
    #: samples, around a core of 326 of 358, sitting where a channel mouth
    #: would be. Dropping atoms changes the answer, so each range is reported.
    drop_residues: Sequence[Tuple[int, int]] = ()
    #: What to do when several copies of one kept ligand touch the kept chain.
    #: ``"all"`` keeps them (faithful; correct for a protein with two real
    #: cofactor sites). ``"best"`` keeps the copy with the most protein
    #: contacts, which is what a transplanted model needs: AlphaFill will
    #: place four copies of GPP from four donors, three of them redundant or
    #: half on the surface, and each extra copy is a wall in somebody's
    #: channel. Both choices are reported per copy.
    ligand_copies: str = "all"


@dataclass
class PrepReport:
    """What prepare() did, in enough detail to argue with."""

    chains_kept: List[str] = field(default_factory=list)
    chains_dropped: List[str] = field(default_factory=list)
    relabelled_to: Optional[str] = None
    ligands_kept: Dict[str, int] = field(default_factory=dict)
    ligands_dropped: Dict[str, int] = field(default_factory=dict)
    #: One row per copy of a kept ligand that came from another chain:
    #: ``(resname, chain, resseq, contacts, kept)``. Printed in full, because
    #: "which copy" is a decision and a silent one is a wrong one.
    ligand_copies: List[Tuple[str, str, int, int, bool]] = field(default_factory=list)
    #: ``(resname, chain, old resseq, new resseq)`` for a ligand renumbered to
    #: avoid colliding with a protein residue after `--relabel`.
    ligands_renumbered: List[Tuple[str, str, int, int]] = field(default_factory=list)
    #: ``[(first, last)]`` residue ranges dropped by `--drop-resi`.
    residues_dropped: List[Tuple[int, int]] = field(default_factory=list)
    atoms_dropped_by_range: int = 0
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
        for resname, chain, resseq, contacts, kept in self.ligand_copies:
            lines.append("         %s %s/%d from another chain, %d contact(s): %s"
                         % (resname, chain, resseq, contacts,
                            "kept" if kept else "dropped as a redundant copy"))
        for resname, chain, old, new in self.ligands_renumbered:
            lines.append("         %s %s/%d renumbered to %d (it collided with "
                         "a protein residue)" % (resname, chain, old, new))
        if self.residues_dropped:
            lines.append("residues dropped %s (%d atom(s))" % (
                ", ".join("%d-%d" % r for r in self.residues_dropped),
                self.atoms_dropped_by_range))
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
    atoms = read_structure(structure) if isinstance(structure, str) else structure
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
    if params.ligand_copies not in ("all", "best"):
        raise ValueError("ligand_copies must be 'all' or 'best', got %r"
                         % (params.ligand_copies,))

    adopted = _adopt_ligands(atoms, want, keep, params, rep)

    kept = []
    for a in atoms.atoms:
        if a.chain not in want and _res_key(a) not in adopted:
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

    if params.drop_residues:
        ranges = [(int(a), int(b)) for a, b in params.drop_residues]
        for a, b in ranges:
            if a > b:
                raise ValueError(
                    "--drop-resi range %d-%d runs backwards" % (a, b))
        before = len(kept)
        kept = [x for x in kept
                if x.hetatm or not any(a <= x.resseq <= b for a, b in ranges)]
        rep.residues_dropped = ranges
        rep.atoms_dropped_by_range = before - len(kept)

    if not params.keep_altlocs:
        kept, n = _resolve_altlocs(kept)
        rep.altlocs_resolved = n

    if params.relabel is not None:
        import dataclasses

        mapping = _relabel_map(str(params.relabel), want)
        if len(mapping) == 1:
            new = next(iter(mapping.values()))
            if adopted:
                kept = _renumber_collisions(kept, new, adopted, rep)
            kept = [dataclasses.replace(a, chain=new) for a in kept]
            rep.relabelled_to = new
        else:
            # A mapping, one new id per kept chain. An adopted ligand follows
            # the chain it was adopted INTO, which is the first kept chain:
            # it was selected by contact with the kept protein as a whole and
            # a ligand cannot belong to two protomers at once.
            first = mapping[want[0]]
            if adopted:
                kept = _renumber_collisions(kept, first, adopted, rep)
            kept = [dataclasses.replace(
                a, chain=mapping.get(a.chain, first)
                if _res_key(a) not in adopted else first) for a in kept]
            rep.relabelled_to = ",".join(
                "%s->%s" % (k, mapping[k]) for k in want)

    rep.atoms_out = len(kept)
    asked = keep - set(rep.ligands_kept)
    if asked:
        rep.notes.append(
            "asked to keep %s, which %s not in the kept chain(s)"
            % (", ".join(sorted(asked)), "is" if len(asked) == 1 else "are"))
    if not any(not a.hetatm for a in kept):
        raise ValueError("nothing is left after preparation: no polymer atoms kept")
    return AtomArray(kept), rep


def _relabel_map(spec: str, want: Sequence[str]) -> Dict[str, str]:
    """``{old chain: new chain}`` from a `--relabel` value.

    Two forms, because two things are wanted and only one of them used to be
    possible. A bare id (``U``) renames the single kept chain, and refuses
    more than one, because collapsing two protomers onto one id merges them.
    A mapping (``A:U,B:V``) renames each kept chain separately, which is NOT
    a merge and is what a kept dimer needs: the pose stage refuses a complex
    whose two sides share a chain id, so a dimeric side has to be given ids
    that the other side does not use.

    The mapping must be injective and must name every kept chain, since a
    chain left out would keep an id the mapping may have just handed to
    another chain.
    """
    spec = spec.strip()
    if ":" not in spec:
        if len(spec) != 1:
            raise ValueError(
                "--relabel takes a single chain id or a mapping like "
                "'A:U,B:V', got %r" % spec)
        if len(want) > 1:
            raise ValueError(
                "cannot rename %d chains (%s) to one id %r: that would merge "
                "separate protomers into one chain. Use a mapping "
                "(--relabel %s) or leave the ids alone."
                % (len(want), "+".join(want), spec,
                   ",".join("%s:%s" % (c, c) for c in want)))
        return {want[0]: spec}
    out: Dict[str, str] = {}
    for part in spec.split(","):
        if part.count(":") != 1:
            raise ValueError("--relabel mapping entry %r is not old:new" % part)
        old_id, new_id = (x.strip() for x in part.split(":"))
        if len(new_id) != 1:
            raise ValueError("chain id %r is not a single character" % new_id)
        if old_id in out:
            raise ValueError("--relabel names chain %r twice" % old_id)
        out[old_id] = new_id
    missing = [c for c in want if c not in out]
    if missing:
        raise ValueError(
            "--relabel does not name kept chain(s) %s; every kept chain needs "
            "a new id or it may collide with one this mapping hands out"
            % ", ".join(missing))
    extra = [c for c in out if c not in want]
    if extra:
        raise ValueError("--relabel names chain(s) %s, which are not kept"
                         % ", ".join(extra))
    if len(set(out.values())) != len(out):
        raise ValueError("--relabel maps two chains onto one id, which would "
                         "merge them: %r" % spec)
    return out


def _res_key(a) -> Tuple[str, int, str, str]:
    """One het residue, identified the way a PDB file identifies it."""
    return (a.chain, a.resseq, a.icode, a.resname)


def _adopt_ligands(atoms: AtomArray, want: Sequence[str], keep, params,
                   rep: "PrepReport"):
    """Which het residues OUTSIDE the kept chains belong to them anyway.

    Only residues whose name was explicitly asked for are considered: a
    neighbouring protomer's cofactor is not adopted because nobody named it.
    Contact is counted against the kept chains' atoms, so a copy sitting on
    the surface loses to one in the site.
    """
    if not keep:
        return set()
    import numpy as np

    poly = np.asarray([a.xyz for a in atoms.atoms
                       if a.chain in want and not a.hetatm], dtype=float)
    if not len(poly):
        return set()

    candidates: Dict[Tuple, list] = {}
    for a in atoms.atoms:
        if a.hetatm and a.chain not in want and (a.resname or "").upper() in keep:
            candidates.setdefault(_res_key(a), []).append(a.xyz)
    if not candidates:
        return set()

    contacts = {}
    for key, xyz in candidates.items():
        d = np.linalg.norm(poly[:, None, :] - np.asarray(xyz)[None, :, :], axis=2)
        contacts[key] = int((d.min(axis=1) <= params.ligand_contact).sum())

    touching = {k for k, n in contacts.items() if n > 0}
    chosen = set(touching)
    if params.ligand_copies == "best":
        by_name: Dict[str, list] = {}
        for k in touching:
            by_name.setdefault(k[3].upper(), []).append(k)
        chosen = set()
        for keys in by_name.values():
            keys.sort(key=lambda k: (-contacts[k], k[0], k[1]))
            chosen.add(keys[0])
    for key in sorted(candidates, key=lambda k: (k[3], -contacts[k], k[0], k[1])):
        rep.ligand_copies.append(
            (key[3], key[0], key[1], contacts[key], key in chosen))
    return chosen


def _renumber_collisions(kept, new_chain, adopted, rep: "PrepReport"):
    """Give an adopted ligand a residue number nothing else in the chain uses.

    Renaming the chain is what makes this necessary: AlphaFill numbers every
    transplant 1, so a ligand moved into the protein's chain becomes A/1 --
    the same address as the protein's first residue. Downstream, a residue is
    selected by chain and number, so the collision is not cosmetic. Numbers
    are identity, not geometry: nothing is moved here either.
    """
    import dataclasses

    # Every non-adopted atom is about to become *new_chain*, so the numbers
    # already taken are taken THERE. Comparing against the chain they still
    # carry finds no collision and lets the ligand keep its number: the whole
    # point of the renumber.
    taken = {(new_chain, a.resseq, a.icode) for a in kept
             if _res_key(a) not in adopted}
    nxt = max((n for _, n, _ in taken), default=0) + 1
    assigned: Dict[Tuple, int] = {}
    out = []
    for a in kept:
        key = _res_key(a)
        if key not in adopted:
            out.append(a)
            continue
        if key not in assigned:
            if (new_chain, a.resseq, a.icode) in taken:
                while (new_chain, nxt, "") in taken:
                    nxt += 1
                assigned[key] = nxt
                taken.add((new_chain, nxt, ""))
                rep.ligands_renumbered.append(
                    (a.resname, key[0], a.resseq, nxt))
                nxt += 1
            else:
                assigned[key] = a.resseq
        new = assigned[key]
        out.append(a if new == a.resseq
                   else dataclasses.replace(a, resseq=new, icode=""))
    return out


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

    Raises if anything kept cannot fit a fixed-column PDB field. mmCIF input
    makes that reachable: it has no column widths, so a two-letter chain id or
    a five-character atom name writes into the field beside it and is then
    parsed, without complaint, as something else.
    """
    problems = pdb_column_problems(atoms.atoms)
    if problems:
        raise ValueError(
            "%s cannot be written as PDB:\n  - %s" % (path, "\n  - ".join(problems)))
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
