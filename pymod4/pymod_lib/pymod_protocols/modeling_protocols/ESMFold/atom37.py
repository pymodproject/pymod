# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Small PDB <-> atom37 helper, vendored from the standalone prototype.

This avoids depending on the ESM SDK (and even on Biopython) for the
coordinate conversion the Forge ``/inverse_fold`` endpoint expects. It is
intentionally simple:

- reads ATOM records from a PDB file;
- selects one chain;
- maps standard amino-acid residues into AlphaFold/ESM-style atom37 order;
- missing atoms are represented as None, which becomes JSON null.

For inverse folding, Forge expects coordinates as N x 37 x 3.

Known limitations (see the prototype README): insertion codes, multiple
models, nonstandard residues, missing backbone atoms, alternate conformations
beyond blank/A, and mmCIF input are handled only minimally. For PyMod the input
is normally a single, clean chain exported from PyMOL, which this covers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Iterable


ATOM37_NAMES = [
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1",
    "SG", "CD", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD",
    "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2",
    "CH2", "NH1", "NH2", "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT",
]
ATOM37_INDEX = {name: i for i, name in enumerate(ATOM37_NAMES)}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # common protonation / special names mapped conservatively
    "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C",
    "MSE": "M",
}
AA1_TO_3 = {v: k for k, v in AA3_TO_1.items()
            if len(k) == 3 and k not in {"HID", "HIE", "HIP", "CYX", "CYM", "MSE"}}
AA1_TO_3["M"] = "MET"


@dataclass
class ResidueAtom37:
    chain_id: str
    resseq: int
    icode: str
    resname: str
    aa: str
    coords: List[Optional[List[float]]]


def _parse_pdb_atom_line(line: str):
    record = line[0:6].strip()
    atom = line[12:16].strip()
    altloc = line[16].strip()
    resname = line[17:20].strip()
    chain = line[21].strip() or "_"
    resseq = int(line[22:26])
    icode = line[26].strip()
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    return record, atom, altloc, resname, chain, resseq, icode, x, y, z


def pdb_to_atom37(
    pdb_path: str,
    chain_id: Optional[str] = None,
    *,
    include_hetatm: bool = False,
) -> Tuple[str, List[List[Optional[List[float]]]], List[ResidueAtom37]]:
    """
    Convert a PDB chain to sequence + N x 37 x 3 coordinates.

    Missing atoms are None. Only standard amino-acid residues are retained.
    For alternate locations, blank or A is preferred. If ``chain_id`` is None
    and several chains are present, the first chain is used.
    """
    residues: Dict[Tuple[str, int, str], ResidueAtom37] = {}
    chain_seen_order: List[str] = []

    with open(pdb_path, "r") as handle:
        for line in handle:
            if not (line.startswith("ATOM  ") or (include_hetatm and line.startswith("HETATM"))):
                continue

            try:
                record, atom, altloc, resname, chain, resseq, icode, x, y, z = _parse_pdb_atom_line(line)
            except Exception:
                continue

            if chain not in chain_seen_order:
                chain_seen_order.append(chain)

            if chain_id is not None and chain != chain_id:
                continue

            if altloc not in ("", "A"):
                continue

            resname = resname.upper()
            if resname not in AA3_TO_1:
                continue

            if atom not in ATOM37_INDEX:
                continue

            key = (chain, resseq, icode)
            if key not in residues:
                residues[key] = ResidueAtom37(
                    chain_id=chain,
                    resseq=resseq,
                    icode=icode,
                    resname=resname,
                    aa=AA3_TO_1[resname],
                    coords=[None for _ in ATOM37_NAMES],
                )
            residues[key].coords[ATOM37_INDEX[atom]] = [x, y, z]

    if not residues and chain_id is None and chain_seen_order:
        raise ValueError("No standard amino-acid residues found in %s" % pdb_path)

    ordered = [
        residues[k]
        for k in sorted(residues.keys(), key=lambda t: (t[0], t[1], t[2]))
    ]

    # If no chain_id was passed and multiple chains are present, keep the first.
    if chain_id is None and ordered:
        first_chain = ordered[0].chain_id
        ordered = [r for r in ordered if r.chain_id == first_chain]

    if not ordered:
        raise ValueError(
            "No residues found for chain %r in %s. Available chains seen: %s"
            % (chain_id, pdb_path, chain_seen_order)
        )

    sequence = "".join(r.aa for r in ordered)
    coordinates = [r.coords for r in ordered]
    return sequence, coordinates, ordered


def coordinates_with_null_atoms(coordinates):
    """
    Return a copy of an atom37 coordinate list in the shape the Forge REST API
    expects: missing atoms (represented internally as ``None``) become
    ``[None, None, None]`` -- a valid 3-element list whose entries are null.

    The API's own fold output uses this convention. A bare ``None`` serializes
    to JSON ``null``, which the structure-input validator rejects with
    "Input should be a valid list"; filling missing atoms with a sentinel
    coordinate (e.g. the origin) instead badly distorts the structure the model
    sees (much worse inverse-folding agreement). ``[None, None, None]`` is both
    accepted and faithful.
    """
    return [[([None, None, None] if atom is None else atom) for atom in residue]
            for residue in coordinates]


def write_fasta(sequence: str, out_path: str, header: str = "designed_sequence") -> None:
    with open(out_path, "w") as handle:
        handle.write(">%s\n" % header)
        for i in range(0, len(sequence), 80):
            handle.write(sequence[i:i + 80] + "\n")


def read_fasta(path: str) -> str:
    seq_parts = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts)


def _write_chain_atoms(handle, sequence, coordinates, chain_id, b_list, serial):
    """
    Write the ATOM records for a single chain and return the next serial number.

    Residues are numbered from 1 within the chain; the serial number is shared
    across the whole file and so is passed in / returned. Missing atoms (null
    coordinates) are skipped. ``b_list`` (e.g. pLDDT) goes in the B-factor
    column.
    """
    for i, aa in enumerate(sequence):
        resname = AA1_TO_3.get(aa.upper(), "UNK")
        resseq = i + 1
        residue_coords = coordinates[i]

        for atom_name, xyz in zip(ATOM37_NAMES, residue_coords):
            if xyz is None:
                continue
            if len(xyz) != 3 or any(v is None for v in xyz):
                continue

            x, y, z = xyz
            element = atom_name[0]

            bfac = float(b_list[i]) if i < len(b_list) else 0.0
            handle.write(
                "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                % (serial, _format_atom_name(atom_name), resname, chain_id,
                   resseq, x, y, z, 1.00, bfac, element)
            )
            serial += 1
    return serial


def write_pdb_from_atom37(
    sequence: str,
    coordinates: List[List[Optional[List[float]]]],
    out_path: str,
    chain_id: str = "A",
    b_factors: Optional[Iterable[float]] = None,
) -> None:
    """
    Write a simple single-chain PDB from atom37 coordinates.

    Intentionally minimal and mostly for visualization in PyMOL. All non-null
    atom37 atoms are written. If pLDDT is passed as ``b_factors``, it is stored
    in the B-factor column.
    """
    b_list = list(b_factors) if b_factors is not None else [0.0] * len(sequence)
    with open(out_path, "w") as handle:
        _write_chain_atoms(handle, sequence, coordinates, chain_id, b_list, 1)
        handle.write("TER\nEND\n")


# Chain IDs assigned, in order, to the chains of a multimer model.
PDB_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def write_multichain_pdb_from_atom37(chains, out_path: str) -> None:
    """
    Write a multi-chain PDB from a list of per-chain atom37 records.

    ``chains`` is a list of ``(sequence, coordinates, b_factors)`` tuples (one
    per chain, in order). Chains get consecutive PDB chain IDs (A, B, C, ...),
    each numbered from residue 1, separated by ``TER`` records. The serial
    number is continuous across the whole file. This is the format PyMod's
    structure parser splits back into one element/object per chain.
    """
    serial = 1
    with open(out_path, "w") as handle:
        for chain_number, (sequence, coordinates, b_factors) in enumerate(chains):
            chain_id = PDB_CHAIN_IDS[chain_number % len(PDB_CHAIN_IDS)]
            b_list = (list(b_factors) if b_factors is not None
                      else [0.0] * len(sequence))
            serial = _write_chain_atoms(handle, sequence, coordinates,
                                        chain_id, b_list, serial)
            handle.write("TER\n")
        handle.write("END\n")


def _format_atom_name(atom_name: str) -> str:
    """
    PDB atom-name column formatting: names shorter than 4 chars are padded so
    the element starts in column 14 (one leading space), matching how PyMOL and
    most parsers expect single-letter-element atoms.
    """
    if len(atom_name) >= 4:
        return atom_name
    return " " + atom_name.ljust(3)
