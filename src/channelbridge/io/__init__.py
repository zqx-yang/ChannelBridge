"""Readers and writers.

Everything downstream of the prep stage is PDB, because CAVER reads PDB only.
The prep stage itself reads either format: see :mod:`channelbridge.io.mmcif`
for why mmCIF has to be accepted at the door, and which of its two numbering
schemes is the right one.
"""

from __future__ import annotations

from .mmcif import is_mmcif, pdb_column_problems, read_mmcif
from .pdb import Atom, AtomArray, has_crlf, read_pdb, write_pdb

__all__ = ["Atom", "AtomArray", "has_crlf", "is_mmcif", "pdb_column_problems",
           "read_mmcif", "read_pdb", "read_structure", "write_pdb"]


def read_structure(path, *, models="first") -> AtomArray:
    """Read a structure, PDB or mmCIF, chosen by the file extension.

    Only the prep stage should call this. A later stage that accepts mmCIF
    would be accepting a file CAVER cannot read, and the mismatch would
    surface as an empty tunnel rather than as a wrong file type.
    """
    if is_mmcif(path):
        return read_mmcif(path, models=models)
    return read_pdb(path, models=models)
