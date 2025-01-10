# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Methods to compute conservation values in an alignment column.
"""

def get_conservation_symbol(column):
    """
    Calculate the conservation symbol for an alignment based on the positively scoring groups that
    occur in the Gonnet Pam250 matrix (see: http://www.clustal.org/download/clustalx_help.html).
    Takes as an input an alignment column represented by a list.

    "*" indicates positions which have a single, fully conserved residue.

    ":" indicates that one of the following 'strong' groups is fully conserved:
    STA, NEQK, NHQK, NDEQ, QHRK, MILV, MILF, HY, FYW

    "." indicates that one of the following 'weaker' groups is fully conserved:
    CSA, ATV, SAG, STNK, STPA, SGND, SNDEQK, NDEQHK, NEQHRK, FVLIM, HFY
    """

    symbol = "-"
    # If there is a gap in that position of the alignment, then the symbol is automatically a
    # "-".
    if "-" in column:
        symbol = "-"
    else:
        # If it is really slow try to use frozenset.
        residues = set(column)
        # All residues in the column are identical: full conservation.
        if len(residues) == 1:
            symbol = "*"
        else:
            # Strong conservation.
            if residues.issubset("STA"):
                symbol = ":"
            elif residues.issubset("NEQK"):
                symbol = ":"
            elif residues.issubset("NHQK"):
                symbol = ":"
            elif residues.issubset("NDEQ"):
                symbol = ":"
            elif residues.issubset("QHRK"):
                symbol = ":"
            elif residues.issubset("MILV"):
                symbol = ":"
            elif residues.issubset("MILF"):
                symbol = ":"
            elif residues.issubset("HY"):
                symbol = ":"
            elif residues.issubset("FYW"):
                symbol = ":"
            # Weak conservation.
            elif residues.issubset("CSA"):
                symbol = "."
            elif residues.issubset("ATV"):
                symbol = "."
            elif residues.issubset("SAG"):
                symbol = "."
            elif residues.issubset("STNK"):
                symbol = "."
            elif residues.issubset("STPA"):
                symbol = "."
            elif residues.issubset("SGND"):
                symbol = "."
            elif residues.issubset("SNDEQK"):
                symbol = "."
            elif residues.issubset("NDEQHK"):
                symbol = "."
            elif residues.issubset("NEQHRK"):
                symbol = "."
            elif residues.issubset("FVLIM"):
                symbol = "."
            elif residues.issubset("HFY"):
                symbol = "."
    return symbol


def compute_stars(elements, adjust_elements=False):
    """
    Computes the "stars" of an alignment.
    """
    if adjust_elements:
        for e in elements:
            e.my_sequence = e.my_sequence.rstrip("-")
        max_length = max([len(e.my_sequence) for e in elements])
        for e in elements:
            e.my_sequence = e.my_sequence.ljust(max_length, "-")

    sequences = [e.my_sequence for e in elements] # str(e.my_sequence)
    if not adjust_elements:
        minimum_length = max([len(s) for s in sequences]) # min([len(s) for s in sequences])
    else:
        minimum_length = max_length

    stars = ""
    for i in range(0, minimum_length):
        column = [_get_ali_char(s, i) for s in sequences] # column = [s[i] for s in sequences]
        stars += get_conservation_symbol(column)

    return stars

def _get_ali_char(seq, i):
    if len(seq) > i:
        return seq[i]
    else:
        return "-"
