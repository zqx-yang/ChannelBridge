"""Minimal, dependency-free PDB reading and writing.

Deliberately not a general PDB library. It exists so that the geometry layer
(:mod:`channelbridge.tunnel`, :mod:`channelbridge.geometry`) can be used and
tested without PyMOL, which is only needed for the prep and pose stages.

Two things here are easy to get wrong and are the reason this module exists:

* **Column slices, not whitespace splitting.** PDB is a fixed-column format.
  ``line.split()`` silently merges fields when a coordinate reaches five
  digits or an atom name fills all four columns.
* **CAVER puts the sphere radius in the B-factor column.** A tunnel cluster
  file (``data/clusters/tun_cl_*.pdb``) has an empty occupancy field and the
  probe radius in columns 61-66::

      ATOM      1  H   FIL T   1      42.446  12.914  -6.843        3.20

  Reading that as an occupancy, or splitting on whitespace, loses the radius.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Iterable, Iterator, Optional, Sequence

import numpy as np

__all__ = ["Atom", "AtomArray", "has_crlf", "iter_atom_lines", "read_pdb", "write_pdb"]

# 0-based slices for the fields we use (PDB v3.3 fixed columns).
_REC = slice(0, 6)
_SERIAL = slice(6, 11)
_NAME = slice(12, 16)
_ALTLOC = slice(16, 17)
_RESNAME = slice(17, 20)
_CHAIN = slice(21, 22)
_RESSEQ = slice(22, 26)
_ICODE = slice(26, 27)
_X = slice(30, 38)
_Y = slice(38, 46)
_Z = slice(46, 54)
_OCC = slice(54, 60)
_BFAC = slice(60, 66)
_ELEMENT = slice(76, 78)


@dataclass(frozen=True)
class Atom:
    """One ATOM/HETATM record, with the fields this package needs."""

    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    xyz: tuple
    occupancy: Optional[float]
    bfactor: Optional[float]
    element: str
    hetatm: bool

    @property
    def radius(self) -> Optional[float]:
        """Alias for :attr:`bfactor`.

        CAVER writes tunnel probe radii into the B-factor column; naming it
        explicitly keeps call sites honest about what they are reading.
        """
        return self.bfactor


class AtomArray:
    """A column store over a list of :class:`Atom`, for vectorised geometry."""

    __slots__ = ("_xyz", "atoms")

    def __init__(self, atoms: Sequence[Atom]):
        self.atoms = list(atoms)
        self._xyz = (
            np.asarray([a.xyz for a in self.atoms], dtype=float).reshape(-1, 3)
            if self.atoms
            else np.zeros((0, 3), dtype=float)
        )

    def __len__(self) -> int:
        return len(self.atoms)

    def __iter__(self) -> Iterator[Atom]:
        return iter(self.atoms)

    @property
    def xyz(self) -> np.ndarray:
        """``(N, 3)`` coordinates. A view; do not mutate in place."""
        return self._xyz

    @property
    def radii(self) -> np.ndarray:
        """``(N,)`` B-factor column, which for CAVER tunnels is the radius."""
        return np.asarray([0.0 if a.bfactor is None else a.bfactor for a in self.atoms], dtype=float)

    def centroid(self) -> np.ndarray:
        if not len(self):
            raise ValueError("centroid of an empty AtomArray")
        return self._xyz.mean(axis=0)

    def exclude_resnames(self, resnames: Iterable[str]) -> "AtomArray":
        """Drop atoms whose residue name is in *resnames* (case-insensitive).

        This is the CAVER ``exclude_residue_names`` filter. Burial has to be
        counted against the atoms CAVER actually searched, so the same
        exclusion must be applied here or a cofactor CAVER never saw inflates
        the number -- the bug recorded in CAVER_HOWTO section 4.
        """
        drop = {r.strip().upper() for r in resnames}
        return AtomArray([a for a in self.atoms if a.resname.upper() not in drop])

    def select(self, *, chain=None, resname=None, hetatm=None) -> "AtomArray":
        keep = []
        for a in self.atoms:
            if chain is not None and a.chain != chain:
                continue
            if resname is not None and a.resname.upper() != resname.strip().upper():
                continue
            if hetatm is not None and a.hetatm is not hetatm:
                continue
            keep.append(a)
        return AtomArray(keep)


def _maybe_float(raw: str) -> Optional[float]:
    raw = raw.strip()
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def iter_atom_lines(path) -> Iterator[str]:
    """Yield the ATOM/HETATM lines of *path*, ignoring everything else."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "ATOM ", "HETATM")):
                yield line


def read_pdb(path, *, models="first") -> AtomArray:
    """Read ATOM/HETATM records from *path*.

    Parameters
    ----------
    path:
        A ``.pdb`` file. mmCIF is not supported: CAVER cannot read it either,
        so the pipeline writes PDB throughout.
    models:
        ``"first"`` (default) keeps only the first MODEL, which is what a
        single-snapshot CAVER run and a prepared target both contain.
        ``"all"`` concatenates every model -- correct for a CAVER tunnel
        cluster file, where each MODEL is one tunnel of the cluster.

    Notes
    -----
    Hydrogens are kept. CAVER tunnel files name every sphere ``H``/``FIL``,
    so filtering by element here would silently empty a tunnel.
    """
    if models not in ("first", "all"):
        raise ValueError("models must be 'first' or 'all', got %r" % (models,))

    atoms = []
    seen_model = False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = line[_REC].strip()
            if rec == "MODEL":
                if seen_model and models == "first":
                    break
                seen_model = True
                continue
            if rec == "ENDMDL" and models == "first" and seen_model:
                # Keep reading only if no coordinates arrived yet (some writers
                # emit an empty leading model).
                if atoms:
                    break
                continue
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                xyz = (float(line[_X]), float(line[_Y]), float(line[_Z]))
            except ValueError:
                # A truncated or malformed coordinate line. Skipping it
                # quietly would corrupt a tunnel, so make it loud.
                raise ValueError(
                    "%s: cannot parse coordinates from line: %r" % (os.fspath(path), line.rstrip())
                ) from None
            serial_raw = line[_SERIAL].strip()
            try:
                serial = int(serial_raw)
            except ValueError:
                # CAVER 3.01 allows letters as extra digits in the serial for
                # very large files. The serial is not used for geometry, so a
                # fallback is safe here.
                serial = -1
            atoms.append(
                Atom(
                    serial=serial,
                    name=line[_NAME].strip(),
                    altloc=line[_ALTLOC].strip(),
                    resname=line[_RESNAME].strip(),
                    chain=line[_CHAIN].strip(),
                    resseq=int(line[_RESSEQ]) if line[_RESSEQ].strip().lstrip("-").isdigit() else 0,
                    icode=line[_ICODE].strip(),
                    xyz=xyz,
                    occupancy=_maybe_float(line[_OCC]),
                    bfactor=_maybe_float(line[_BFAC]),
                    element=line[_ELEMENT].strip(),
                    hetatm=(rec == "HETATM"),
                )
            )
    return AtomArray(atoms)


def write_pdb(path, atoms: Iterable[Atom], *, newline="\n") -> int:
    """Write *atoms* to *path* with LF line endings, returning the count.

    ``newline="\\n"`` is not a default to argue with. A PDB with CRLF endings
    makes CAVER report "No PDB files found" and **exit 0**, so a runner
    believes the run succeeded (CAVER_HOWTO, trap 4).
    """
    n = 0
    with open(path, "w", encoding="ascii", newline=newline) as fh:
        for a in atoms:
            rec = "HETATM" if a.hetatm else "ATOM  "
            name = a.name if len(a.name) >= 4 else " %-3s" % a.name
            occ = 1.00 if a.occupancy is None else a.occupancy
            bfac = 0.00 if a.bfactor is None else a.bfactor
            fh.write(
                "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                % (
                    rec,
                    a.serial if a.serial >= 0 else (n + 1),
                    name,
                    a.altloc or " ",
                    a.resname,
                    a.chain or " ",
                    a.resseq,
                    a.icode or " ",
                    a.xyz[0],
                    a.xyz[1],
                    a.xyz[2],
                    occ,
                    bfac,
                    a.element,
                )
            )
            n += 1
        fh.write("END\n")
    return n


def has_crlf(path) -> bool:
    """True if *path* contains a CRLF line ending.

    Worth calling on anything about to be handed to CAVER; see
    :func:`write_pdb`.
    """
    with open(path, "rb") as fh:
        return b"\r\n" in fh.read()
