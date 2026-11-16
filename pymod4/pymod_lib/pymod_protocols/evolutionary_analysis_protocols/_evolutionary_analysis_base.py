# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol


# What a nucleotide alignment is made of: the four bases of DNA, the uracil of
# RNA, and 'N' for a base which was not read.
NUCLEOTIDE_LETTERS = set("ACGTUN")
# Characters which say nothing about the kind of sequence.
UNINFORMATIVE_LETTERS = set("-.?*Xx ")


def looks_like_nucleotides(sequences):
    """
    Whether these sequences are made of nucleotides. Every residue has to be
    one: a protein of twenty amino acids of which none is a leucine, a serine
    or any of the other sixteen does not occur.
    """
    letters = set("".join(sequences).upper()) - UNINFORMATIVE_LETTERS
    return bool(letters) and letters <= NUCLEOTIDE_LETTERS


def is_a_nucleic_acid_alignment(cluster_element):
    """
    Whether an alignment is one of nucleic acids.
    """

    # What PyMod worked out when the sequences were loaded.
    if getattr(cluster_element, "polymer_type", "protein") in ("dna", "rna"):
        return True

    # A sequence which was already aligned when it was loaded is read as a
    # protein, since its gaps are not nucleotides, so the sequences are asked
    # as well.
    try:
        sequences = [child.my_sequence for child in cluster_element.get_children()
                     if not child.is_cluster()]
    except AttributeError:
        return False
    return looks_like_nucleotides(sequences)


class Evolutionary_analysis_protocol(PyMod_protocol):

    protocol_name = "evolutionary_analysis"

    def __init__(self, pymod, input_cluster_element, *args):
        self.input_cluster_element = input_cluster_element
        PyMod_protocol.__init__(self, pymod, *args)


    def refuse_a_nucleic_acid_alignment(self, title):
        """
        Says that this analysis cannot be run on an alignment of nucleic acids,
        and whether that is the case. The conservation scores of PyMod are
        computed with amino acid substitution matrices and amino acid
        frequencies, which have nothing to say about nucleotides: the matrices
        have no entry for a pair of bases, they have one for the amino acids
        which happen to be written with the same letters, so the analysis runs
        to the end and returns scores which mean nothing.
        """

        if not is_a_nucleic_acid_alignment(self.input_cluster_element):
            return False
        self.pymod.main_window.show_error_message(
            title,
            "This analysis is not available for nucleic acid alignments: its"
            " scores come from amino acid substitution matrices and residue"
            " frequencies.")
        return True
