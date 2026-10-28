"""Read coordinates from mmCIF, so the prep stage can write the PDB.

CAVER reads PDB only, and so does everything downstream of it here. That is
not a reason to refuse mmCIF at the *door*: most structures worth running
arrive as mmCIF and nothing else. An AlphaFold DB model, an AlphaFill model
and an AlphaFold3 prediction are mmCIF with no PDB offered; a deposited entry
large enough to have no PDB format at all is mmCIF or nothing. So the
conversion belongs in the prep stage -- which is the one stage allowed to
decide which atoms exist -- and nowhere else.

**auth_* ids, not label_*.** This is the one decision in the file and it is
not cosmetic. mmCIF carries two numbering schemes. `label_*` is the
structure-model scheme: every distinct chemical entity gets its own
`label_asym_id`, so in a deposited entry the protein is `A` and *each* of its
cofactors, metals and waters becomes a separate asym id (`B`, `C`, `D`...).
`auth_*` is the depositor's scheme -- the chain letters and residue numbers
that appear in the paper, in ChannelsDB, in a reviewer's email and in every
`--chain A` a user will type.

Reading `label_*` therefore breaks the pipeline in a way that looks like
biology. `chbridge prep --chain A --keep HEM` would keep chain A's protein and
silently drop the haem, because the haem is in chain `B` of the label scheme:
the report says "ligands kept: none" and the channel you then measure runs
straight through the empty cofactor site. Residue numbers diverge too --
`label_seq_id` counts from 1 along the entity and is null for anything that is
not a polymer, so a Tyr314 named in a paper is some other number here, or
none. This reader prefers `auth_*` for the atom name, residue name, residue
number and chain, and falls back to `label_*` per field when a file omits it
(AlphaFold3 output omits nothing relevant; some minimal writers omit
`auth_atom_id`).

Nothing is moved, nothing is added, and hydrogens are kept: the same
contract as :func:`channelbridge.io.pdb.read_pdb`.
"""

from __future__ import annotations

import os
from typing import Dict, Iterator, List, Optional, Sequence

from .pdb import Atom, AtomArray, _maybe_float

__all__ = ["is_mmcif", "pdb_column_problems", "read_mmcif"]

#: mmCIF's two spellings of "no value". Both mean the field is absent, and
#: both have to become an empty string rather than a literal "." in a PDB
#: column -- a "." altloc or insertion code is a character PDB readers treat
#: as real.
NULLS = (".", "?")

#: The fields taken from the file, in the order they are tried. The auth_*
#: spelling wins where both exist; see the module docstring.
_FIELDS = {
    "name": ("auth_atom_id", "label_atom_id"),
    "resname": ("auth_comp_id", "label_comp_id"),
    "chain": ("auth_asym_id", "label_asym_id"),
    "resseq": ("auth_seq_id", "label_seq_id"),
}


def is_mmcif(path) -> bool:
    """True if *path* looks like mmCIF, by extension."""
    return os.path.splitext(os.fspath(path))[1].lower() in (".cif", ".mmcif")


def _tokenise(line: str) -> List[str]:
    """Split one mmCIF data line into values, respecting quotes.

    An atom name may be quoted (``"O5'"``, ``'CA'``) and a quoted value may
    contain spaces, so a bare ``str.split()`` shifts every later column by one
    and the coordinates read as a residue name. Quoting is rare in
    ``atom_site`` and routine in the rest of a file, which is exactly the kind
    of asymmetry that makes a split-based reader work until it doesn't.
    """
    out: List[str] = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c in "'\"":
            # mmCIF ends a quoted value at quote-then-whitespace, so an
            # apostrophe inside a name (5'-phosphate) does not end it.
            j = i + 1
            while j < n:
                if line[j] == c and (j + 1 >= n or line[j + 1] in " \t"):
                    break
                j += 1
            out.append(line[i + 1:j])
            i = j + 1
            continue
        j = i
        while j < n and line[j] not in " \t":
            j += 1
        out.append(line[i:j])
        i = j
    return out


def _atom_site_rows(path) -> Iterator[Dict[str, str]]:
    """Yield one dict per ``atom_site`` row, keyed by the short tag name.

    Handles both forms: the ``loop_`` used by every real structure, and the
    single-row key/value form a one-atom file may use.
    """
    cols: List[str] = []
    single: Dict[str, str] = {}
    pending: List[str] = []
    in_loop = False
    in_rows = False
    text: Optional[List[str]] = None

    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            stripped = line.strip()

            # A semicolon-delimited text field spans lines. Not expected in
            # atom_site, but a value that swallowed the rest of the file would
            # be a silent truncation rather than an error.
            if text is not None:
                if stripped == ";":
                    pending.append("\n".join(text))
                    text = None
                else:
                    text.append(line)
                continue
            if line.startswith(";") and in_rows:
                text = [line[1:]]
                continue

            if stripped == "loop_":
                if in_rows:
                    break
                in_loop, cols, pending = True, [], []
                continue
            if stripped.startswith("_atom_site."):
                tag = stripped.split(".", 1)[1].split()[0]
                if in_loop and not in_rows:
                    cols.append(tag)
                else:
                    parts = _tokenise(stripped)
                    if len(parts) > 1:
                        single[tag] = parts[1]
                continue
            if stripped.startswith("_") or stripped.startswith("data_"):
                # Another category. If rows have been read, atom_site is done.
                if in_rows:
                    break
                in_loop, cols = False, []
                continue
            if not cols:
                continue
            if not stripped or stripped == "#" or stripped.startswith("#"):
                if in_rows:
                    break
                continue
            if stripped in ("stop_",) or stripped.startswith("save_"):
                break

            in_rows = True
            pending.extend(_tokenise(line))
            while len(pending) >= len(cols):
                row = pending[:len(cols)]
                del pending[:len(cols)]
                yield dict(zip(cols, row))

    if pending and cols:
        raise ValueError(
            "%s: the atom_site loop ends with %d leftover value(s) for %d "
            "column(s), so a row is truncated"
            % (os.fspath(path), len(pending), len(cols)))
    if single:
        yield single


def _pick(row: Dict[str, str], field: str) -> str:
    for tag in _FIELDS[field]:
        v = row.get(tag)
        if v is not None and v not in NULLS:
            return v
    return ""


def read_mmcif(path, *, models="first") -> AtomArray:
    """Read the ATOM/HETATM records of an mmCIF file.

    Parameters
    ----------
    path:
        A ``.cif`` / ``.mmcif`` file.
    models:
        ``"first"`` (default) keeps the first ``pdbx_PDB_model_num`` seen --
        an NMR ensemble or a multi-sample prediction otherwise arrives as
        several overlapping copies of the protein, which reads as a structure
        with no free space anywhere. ``"all"`` keeps every model.

    Notes
    -----
    Residue names, numbers, chain ids and atom names come from the ``auth_*``
    columns where present. The module docstring says why that matters.
    """
    if models not in ("first", "all"):
        raise ValueError("models must be 'first' or 'all', got %r" % (models,))

    atoms: List[Atom] = []
    first_model: Optional[str] = None
    for row in _atom_site_rows(path):
        model = row.get("pdbx_PDB_model_num", "")
        if models == "first":
            if first_model is None:
                first_model = model
            elif model != first_model:
                continue
        try:
            xyz = (float(row["Cartn_x"]), float(row["Cartn_y"]),
                   float(row["Cartn_z"]))
        except (KeyError, ValueError):
            raise ValueError(
                "%s: cannot parse coordinates from atom_site row %r"
                % (os.fspath(path), row)) from None
        serial_raw = (row.get("id") or "").strip()
        try:
            serial = int(serial_raw)
        except ValueError:
            serial = -1
        seq_raw = _pick(row, "resseq")
        altloc = row.get("label_alt_id", "")
        icode = row.get("pdbx_PDB_ins_code", "")
        atoms.append(Atom(
            serial=serial,
            name=_pick(row, "name"),
            altloc="" if altloc in NULLS else altloc.strip(),
            resname=_pick(row, "resname"),
            chain=_pick(row, "chain"),
            resseq=int(seq_raw) if seq_raw.lstrip("-").isdigit() else 0,
            icode="" if icode in NULLS else icode.strip(),
            xyz=xyz,
            occupancy=_maybe_float(row.get("occupancy", "")),
            bfactor=_maybe_float(row.get("B_iso_or_equiv", "")),
            element=(row.get("type_symbol", "") or "").strip(),
            hetatm=(row.get("group_PDB", "ATOM").strip().upper() == "HETATM"),
        ))
    if not atoms:
        raise ValueError(
            "%s: no atom_site records found. An mmCIF with only metadata, or "
            "a file that is not mmCIF despite its name" % os.fspath(path))
    return AtomArray(atoms)


def pdb_column_problems(atoms: Sequence[Atom]) -> List[str]:
    """What in *atoms* a fixed-column PDB file cannot represent.

    mmCIF has no column widths; PDB has. A two-character chain id, a residue
    number past 9999 or an atom name of five characters does not fail to
    write -- it writes into the neighbouring column and corrupts the field
    beside it, which a reader then parses without complaint. Checked before
    writing rather than discovered afterwards.
    """
    problems = []
    wide = sorted({a.chain for a in atoms if len(a.chain) > 1})
    if wide:
        problems.append(
            "chain id(s) %s need more than the 1 column PDB gives them; keep "
            "one chain (--chain) and rename it (--relabel)" % ", ".join(wide))
    big = [a for a in atoms if not -999 <= a.resseq <= 9999]
    if big:
        problems.append(
            "%d residue number(s) are outside the PDB range -999..9999 (first: "
            "%s%d)" % (len(big), big[0].resname, big[0].resseq))
    long_names = sorted({a.name for a in atoms if len(a.name) > 4})
    if long_names:
        problems.append(
            "atom name(s) %s are longer than the 4 columns PDB gives them"
            % ", ".join(long_names[:5]))
    over = [a for a in atoms if a.serial > 99999]
    if over:
        problems.append(
            "%d atom serial(s) are past 99999, which overruns the field and "
            "shifts every column after it" % len(over))
    return problems
