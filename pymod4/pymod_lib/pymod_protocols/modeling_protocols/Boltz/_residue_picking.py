# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Naming residues of a chain: the epitope of a design, the pocket of a prediction.

Both ask the user the same question -- which residues of which chain -- and both
have the same trap in them. The API counts the residues of a chain from zero
("All residue indices in constraints are 0-indexed"), while the numbers a user
reads in PyMod and in PyMOL are those of the structure, which start where the PDB
file starts and may skip. Chain H of 3DD2 begins at 104, so its residue 110 is
index 6.

The numbers are therefore looked up among the residues of the element rather than
offset by a constant, and a number the chain does not have is refused by name.
"""

from pymol import cmd

from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_design as bd


def compress(numbers):
    """[45, 46, 47, 60] as "45-47,60", which is what a field accepts back."""

    pieces, start, previous = [], None, None
    for number in numbers:
        if start is None:
            start = previous = number
            continue
        if number == previous + 1:
            previous = number
            continue
        pieces.append(str(start) if start == previous else "%s-%s" % (start, previous))
        start = previous = number
    if start is not None:
        pieces.append(str(start) if start == previous else "%s-%s" % (start, previous))
    return ",".join(pieces)


def residue_positions(element):
    """
    {the number the user reads: the position in the sequence which is sent}.

    The position is the index of the residue among the polymer residues of the
    element, which is exactly the sequence a request carries once its gaps are
    removed.
    """

    positions = {}
    try:
        residues = element.get_polymer_residues()
    except Exception:
        residues = []
    for position, residue in enumerate(residues):
        number = getattr(residue, "db_index", None)
        if number is None:
            continue
        try:
            number = int(number)
        except (TypeError, ValueError):
            continue
        positions.setdefault(number, position)
    return positions


def element_name(element):
    return (getattr(element, "my_header", None)
            or getattr(element, "compact_header", None) or "?")


class Residue_picking_mixin:
    """
    Turns residue numbers typed by a user, or selected in PyMOL, into the indices
    a Boltz request carries. Mixed into a protocol, which supplies the chains:
    'assignment' is [(element, chain_id)], as the request numbers them.
    """

    def residues_to_indices(self, text, assignment):
        """{chain_id: [0-based indices]} from what the user typed, or {}."""

        chain_ids = [chain_id for _, chain_id in assignment]
        numbers_per_chain = bd.parse_residue_numbers(text, chain_ids)
        if not numbers_per_chain:
            return {}

        indices = {}
        for element, chain_id in assignment:
            numbers = numbers_per_chain.get(chain_id)
            if not numbers:
                continue

            positions = residue_positions(element)
            chosen, missing = [], []
            for number in numbers:
                if number in positions:
                    chosen.append(positions[number])
                else:
                    missing.append(number)

            if missing:
                if not positions:
                    raise bd.BoltzDesignError(
                        "The residues of %s are not known to PyMod, so residue %s"
                        " cannot be found in it."
                        % (element_name(element),
                           ", ".join(str(number) for number in missing[:8])))
                raise bd.BoltzDesignError(
                    "Chain %s (%s) has no residue %s. Its residues are numbered"
                    " %s to %s."
                    % (chain_id, element_name(element),
                       ", ".join(str(number) for number in missing[:8]),
                       min(positions), max(positions)))
            indices[chain_id] = sorted(chosen)

        return indices


    def residue_text_from_pymol(self, assignment, selection="sele"):
        """
        The residues of a PyMOL selection which belong to these chains, written as
        a field accepts them: "H:45-52,60, L:88".

        Returns (text, message): the text for the field, and what to tell the user
        when nothing came of it.
        """

        # No selection at all and an empty one are the same thing to a user.
        try:
            selected_atoms = cmd.count_atoms(selection)
        except Exception:
            selected_atoms = 0
        if not selected_atoms:
            return "", ("There is nothing selected in PyMOL. Select the residues"
                        " in the viewer, then press this button again.")

        pieces, without_structure = [], []
        for element, chain_id in assignment:
            if not element.has_structure():
                without_structure.append(element_name(element))
                continue

            found = set()
            try:
                cmd.iterate("(%s) and (%s) and polymer and name CA"
                            % (selection, element.get_pymol_selector()),
                            "found.add(int(resv))", space={"found": found})
            except Exception:
                continue
            if found:
                pieces.append("%s:%s" % (chain_id, compress(sorted(found))))

        if pieces:
            return ", ".join(pieces), ""

        if without_structure:
            return "", ("The selected sequences have no structure in PyMOL (%s),"
                        " so a selection cannot say which residues they are."
                        " Type the residue numbers instead."
                        % ", ".join(without_structure))
        return "", ("Nothing in the selection belongs to the sequences of this"
                    " run.")
