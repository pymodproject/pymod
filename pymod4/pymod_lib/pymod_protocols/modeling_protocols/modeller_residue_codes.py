# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
The alphabet MODELLER uses for the residues it reads, and what it means for an
alignment.

A module of its own, with no imports, so that the rule can be read and tested on
its own: it decides how the target row of a PIR alignment is built, and getting it
wrong is silent -- the rows still come out the same length, MODELLER still builds a
model, and the model is wrong.
"""


def is_modelled_residue_code(code):
    """
    Tells whether a residue as MODELLER reads it is a residue the target sequence
    has to provide, or a block residue the model just inherits from the template.

    MODELLER writes 'residue.code' with an alphabet where the case carries the
    meaning: standard amino acids are uppercase one-letter codes, while nucleic
    acids and everything read as a heteroresidue are lowercase or '.'. The DNA of a
    zinc finger complex reads as 'eljltllljlt', a zinc as 'z', a water as 'w'.

    So only an uppercase code is a position the target aligns to. Anything else is a
    block residue, and the target must carry a '.' in front of it, which is how a
    PIR alignment says "keep this residue as it is in the template".
    """

    return code.isupper()


# The residue names of the bases, as they are written in a PDB file. Read from the
# name and not from the record type, since a base can be an ATOM or a HETATM
# depending on who wrote the file.
DNA_BASE_NAMES = frozenset(("DA", "DT", "DG", "DC", "DU"))
RNA_BASE_NAMES = frozenset(("A", "U", "G", "C", "T"))


def is_nucleic_residue_name(residue_name):
    """
    Whether a residue name from a PDB file is a nucleic acid.

    Used to leave the nucleic acids out of the file given to MODELLER, which is the
    only way to leave them out of a model: MODELLER reads them whatever env.io.hetatm
    says.
    """

    name = (residue_name or "").strip().upper()
    return (name in DNA_BASE_NAMES
            or name in RNA_BASE_NAMES
            or ("D" + name) in DNA_BASE_NAMES)


def file_has_nucleic_residues(file_path):
    """
    Whether a structure file carries any nucleic acid at all. Read line by line and
    stopped at the first base, so it costs nothing on a file which has none.
    """

    try:
        with open(file_path) as structure_file:
            for line in structure_file:
                if (line.startswith(("ATOM", "HETATM")) and len(line) >= 27
                        and is_nucleic_residue_name(line[17:20])):
                    return True
    except (OSError, UnicodeDecodeError):
        return False
    return False


def nucleic_residues_in_order(lines):
    """
    The nucleic acid residues of a PDB file, in the order they appear, as a list of
    the (residue number, insertion code) that identifies each of them.
    """

    residues = []
    seen = set()
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
            continue
        if not is_nucleic_residue_name(line[17:20]):
            continue
        key = (line[22:26], line[26])
        if key not in seen:
            seen.add(key)
            residues.append(key)
    return residues


def split_nucleic_acids_into_chains(lines, chain_lengths):
    """
    Gives the nucleic acids of a model their own chains back.

    To be read by MODELLER a template has its nucleic acids merged into the protein
    chain, as heteroresidues renumbered after the protein. The model comes back the
    same way, which is not how anyone looks at DNA: in PyMOL a double helix is two
    chains, one per filament.

    'chain_lengths' is the [(chain id, number of residues)] of the filaments as they
    were in the structure the user opened, in the order they appeared there, which is
    the order the merge kept. The bases of the model are handed back to those chains
    in that order.

    Returns the new lines, or None when the count does not add up -- which means the
    model is not the shape this assumes, and the safe thing is to change nothing.
    """

    residues = nucleic_residues_in_order(lines)
    expected = sum(count for _, count in chain_lengths)
    if not residues or expected != len(residues):
        return None

    # Which chain each residue goes back to.
    chain_of_residue = {}
    position = 0
    for chain_id, count in chain_lengths:
        for key in residues[position:position + count]:
            chain_of_residue[key] = chain_id
        position += count

    new_lines = []
    for line in lines:
        if (line.startswith(("ATOM", "HETATM")) and len(line) >= 27
                and is_nucleic_residue_name(line[17:20])):
            chain_id = chain_of_residue.get((line[22:26], line[26]))
            if chain_id:
                line = line[:21] + chain_id[0] + line[22:]
        # A TER carries a chain too, and a stale one would claim the bases still
        # belong to the protein chain.
        elif line.startswith("TER") and len(line) >= 27:
            chain_id = chain_of_residue.get((line[22:26], line[26]))
            if chain_id:
                line = line[:21] + chain_id[0] + line[22:]
        new_lines.append(line)
    return new_lines
