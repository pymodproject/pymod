# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Turns an alignment built in PyMod into the a3m the Boltz API expects.

This is what makes the integration worth more than the web interface of the
service: instead of letting the server search for homologues, a prediction can be
run on an alignment the user has built and curated in PyMod, with whichever tool
and database they chose (BLAST, the HMMER family, HH-suite, or one loaded from a
file).

The API accepts a custom alignment inline, base64 encoded, so nothing has to be
published anywhere; 'boltz_payload.msa_custom' does that encoding. This module only
produces the a3m text.

**The a3m format.** The first sequence is the query, written WITHOUT gaps: its
residues define the "match" columns of the alignment. Every other sequence is
written against those columns:

  - in a match column (the query has a residue) the character is upper case, or
    '-' when that sequence has a gap there;
  - in an insert column (the query has a gap) the character is lower case, and is
    left out entirely when that sequence also has a gap there.

So every sequence after the first has exactly as many upper case characters and
dashes as the query has residues, and any number of lower case ones in between.
That is what distinguishes a3m from a plain aligned FASTA, in which the insertions
of the query would be kept as gap columns for everyone.
"""


# Characters which count as a gap in a PyMod alignment.
GAP_CHARS = "-."


class BoltzMSAError(ValueError):
    """
    Raised when an alignment can not be converted. The message is meant to be
    shown to the user as it is.
    """
    pass


def _is_gap(character):
    return character in GAP_CHARS


def alignment_to_a3m(query_aligned, homologs_aligned,
                     query_name="query", homolog_names=None):
    """
    Converts an alignment to a3m and returns the text.

    Parameters
    ----------
    query_aligned:
        The row of the sequence being predicted, with its gaps.
    homologs_aligned:
        The other rows, all of the same length as the query row.
    query_name, homolog_names:
        Names for the FASTA headers. Missing names are numbered.

    Raises BoltzMSAError when the rows do not all have the same length, or when
    the query has no residue at all.
    """

    if query_aligned is None:
        raise BoltzMSAError("The alignment has no query sequence.")

    query_aligned = "".join(str(query_aligned).split())
    homologs_aligned = ["".join(str(row).split()) for row in (homologs_aligned or [])]

    length = len(query_aligned)
    for index, row in enumerate(homologs_aligned):
        if len(row) != length:
            raise BoltzMSAError(
                "The sequences of the alignment do not all have the same length:"
                " the query has %s columns while sequence %s has %s. The"
                " alignment can not be used."
                % (length, index + 1, len(row)))

    # The columns in which the query has a residue are the match columns.
    match_columns = [i for i, character in enumerate(query_aligned)
                     if not _is_gap(character)]
    if not match_columns:
        raise BoltzMSAError("The query sequence of the alignment is empty.")

    if homolog_names is None:
        homolog_names = []
    names = list(homolog_names) + ["sequence_%s" % (i + 1)
                                   for i in range(len(homolog_names),
                                                  len(homologs_aligned))]

    lines = [">%s" % _clean_name(query_name),
             "".join(query_aligned[i] for i in match_columns).upper()]

    match_column_set = set(match_columns)
    for row, name in zip(homologs_aligned, names):
        characters = []
        for column, character in enumerate(row):
            if column in match_column_set:
                # Aligned to a residue of the query, or a gap in front of one.
                characters.append("-" if _is_gap(character) else character.upper())
            elif not _is_gap(character):
                # Inserted with respect to the query.
                characters.append(character.lower())
        lines.append(">%s" % _clean_name(name))
        lines.append("".join(characters))

    return "\n".join(lines) + "\n"


def _clean_name(name):
    """
    A FASTA header must be a single line, and the service reads the names, so
    anything which would break the file is removed.
    """
    text = " ".join(str(name or "sequence").split())
    return text or "sequence"


def count_sequences(a3m_text):
    """Number of sequences in an a3m, the query included."""
    return sum(1 for line in (a3m_text or "").splitlines()
               if line.startswith(">"))


###################################################################################################
# Reading an alignment out of PyMod.                                                              #
###################################################################################################

def alignment_of_element(element):
    """
    Returns the other sequences aligned with 'element' in PyMod, or an empty list
    when it is not part of an alignment.

    An element belongs to an alignment when it is the child of a cluster; the
    other children of that cluster are its homologues, and they all carry the
    aligned sequence, gaps included.
    """

    try:
        if not element.is_child():
            return []
        mother = element.mother
        if mother is None or not hasattr(mother, "get_children"):
            return []
        return [child for child in mother.get_children() if child is not element]
    except Exception:
        # An element which is not in a cluster, or a cluster shaped differently
        # from what is expected, simply has no alignment to offer.
        return []


def a3m_for_element(element, max_sequences=None):
    """
    Builds the a3m of a PyMod element from the alignment it belongs to, or returns
    None when it belongs to none, or to one in which it is the only sequence.

    'max_sequences' caps how many homologues are sent, the query excluded. The
    request travels base64 encoded inside the JSON body, so a very deep alignment
    would make it large; None means no limit.
    """

    homologs = alignment_of_element(element)
    if not homologs:
        return None

    if max_sequences is not None and len(homologs) > max_sequences:
        homologs = homologs[:max_sequences]

    return alignment_to_a3m(
        element.my_sequence,
        [child.my_sequence for child in homologs],
        query_name=_header_of(element),
        homolog_names=[_header_of(child) for child in homologs])


def _header_of(element):
    for attribute in ("compact_header", "my_header"):
        value = getattr(element, attribute, None)
        if value:
            return value
    return "sequence"
