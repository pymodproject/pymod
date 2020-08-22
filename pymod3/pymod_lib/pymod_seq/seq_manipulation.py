# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import re

import Bio

from pymod_lib import pymod_vars


###################################################################################################
# Some useful functions used in the plugin.                                                       #
###################################################################################################

def compute_sequence_identity(seq1, seq2, toss_modres = False, return_matches=False):
    """
    Computes the identity [(identical residues)/(total matches in the alignment)] between two
    aligned sequences. If 'toss_modres' is set to 'True', 'X' in sequences will be replaced with
    '-' so that the mdofied residues in PDB chain sequences will not be used to compute the
    sequence identity.
    """
    sequence_identity = 0.0
    identical_positions = 0.0
    aligned_positions = 0.0

    if toss_modres:
        seq1 = seq1.replace("X","-")
        seq2 = seq2.replace("X","-")

    if str(seq1) != str(seq2):
        # Computes the identity the two sequences. Range only for the shortest sequence.
        for i in range( min(len(seq1),len(seq2)) ):
            if seq1[i] != "-" and seq2[i] != "-":
                aligned_positions += 1
                if seq1[i] == seq2[i]:
                    identical_positions += 1
        try:
            sequence_identity = (identical_positions/aligned_positions)*100
        except ZeroDivisionError:
            sequence_identity = 0.0

    # If the two sequence are the same, just return 1 (or 100).
    else:
        sequence_identity = (1.0/1.0)*100

    if not return_matches:
        return round(sequence_identity,2)
    else:
        return round(sequence_identity,2), identical_positions


def get_residue_id_in_aligned_sequence(aligned_sequence, real_id):
    """
    aligned_sequence: a sequence with indels.
    real_id: the id (position -1) of the residue in the gapless sequence.
    returns: the id (position -1) of the same residue in the alignment (the sequence + indels).
    """
    assert( len(str(aligned_sequence).replace("-","")) >= real_id )
    real_counter = 0
    for i,p in enumerate(aligned_sequence):
        if p != "-":
            if real_counter == real_id:
                return i
            else:
                real_counter += 1


def get_residue_id_in_gapless_sequence(aligned_sequence, id_in_alignment):
    """
    This does the opposite of the function above.
    aligned_sequence: a sequence with indels.
    id_in_alignment: the id (position -1) of the same residue in the alignment (the sequence + indels).
    returns: the id (position -1) of the residue in the gapless sequence.
    """
    assert(len(aligned_sequence) >= id_in_alignment)
    residue_counter = 0
    real_id = None
    for (k, r) in enumerate(aligned_sequence):
        if k != id_in_alignment:
            if r != "-":
                residue_counter += 1
        else:
            if r != "-":
                real_id = residue_counter
            break
    return real_id


def find_residue_conservation(template_seq, target_seq, real_id):
    """
    Takes two aligned sequences, a "template" and a "target", and sees if a
    residue with the id "i" of the template is conserved in the target.
    """
    # seq_counter: real id of the residue
    # k: id of the residue in the alignment
    alignment_position = get_residue_id_in_aligned_sequence(template_seq,real_id)
    return template_seq[alignment_position] == target_seq[alignment_position]


def get_leading_gaps(sequence, index):

    gap_counter = 0
    for i, position in enumerate(sequence[index:]):
        if i != 0:
            if position == "-":
                gap_counter += 1
            else:
                break
    return gap_counter


def global_pairwise_alignment(seq1, seq2, toss_modres=False):
    """
    If seq1 contains gaps, also aseq1 will maintain these gaps.
    """

    import Bio.pairwise2

    # The original version was: globalms(seq1, seq2, 2, -1, -.5, -.1)
    ali = Bio.pairwise2.align.globalms(seq1, seq2,
                                       2,  # match score.
                                       -2, # mismatch score.
                                       -5, # gap opening.
                                       -.5 # gap extension.
                                       )

    aseq1 = ali[0][0]
    aseq2 = ali[0][1]
    seqid, matches = compute_sequence_identity(aseq1, aseq2, return_matches=True, toss_modres=toss_modres)
    return {"seq1": aseq1, "seq2": aseq2, "id": seqid, "matches": matches}


def get_limit_residues_ids(seq1, seq2):
    c1 = 0
    c2 = 0
    s1a = None # first match residue id of sequence 1
    e1a = None # last match residue id of sequence 2
    s2a = None
    e2a = None
    last_match_position = 0
    for p1, p2 in zip(seq1, seq2):
        if p1 != "-" and p2 != "-":
            if s1a != None:
                s1a = c1
                s2a = c2
            else:
                e1a = c1
                e2a = c2
        if p1 != "-":
            c1 += 1
        if p2 != "-":
            c2 += 1
    return {"s1a": s1a, "e1a": e1a,"s2a": s2a, "e2a": e2a}


#######################
# Alignments.         #
#######################

def adjust_aligned_elements_length(elements, remove_right_indels=True):
    if len(set([len(e.my_sequence) for e in elements])) == 1:
        return False
    # First remove indels at the end of the sequences.
    if remove_right_indels:
        for e in elements:
            e.my_sequence = str(e.my_sequence).rstrip("-")
    # Then pad each sequence with the right number of indels to make them of the same length as
    # the longest sequence.
    max_length = max([len(e.my_sequence) for e in elements])
    for e in elements:
        e.my_sequence = str(e.my_sequence).ljust(max_length,"-")
        # e.set_sequence(str(e.my_sequence).ljust(max_length,"-"), permissive=False)


def all_positions_are_gaps(column):
    for i in column:
        if i != "-":
            return False
    return True


def remove_gap_only_columns(cluster):
    """
    Remove the columns containing only gaps in the child elements of a PyMod cluster element.
    """
    children = cluster.get_children()
    list_seqs = [list(s.my_sequence) for s in children]
    columns_to_remove = []
    for i in range(0, min([len(s) for s in list_seqs])):
        if all_positions_are_gaps([s[i] for s in list_seqs]):
            columns_to_remove.append(i)

    for columns_removed, i in enumerate(columns_to_remove):
        for s in list_seqs:
            s.pop(i-columns_removed)

    for child, list_seq in zip(children, list_seqs):
        child.my_sequence = "".join(list_seq)


#######################
# Clean up sequences. #
#######################

def clean_sequence_from_input(sequence):
    #     sequence = sequence.replace("Z","X") # ambiguity, E or Q
    #     sequence = sequence.replace("B","X") # ambiguity, D or N
    #     sequence = sequence.replace("J","X") # ambiguity, I or L
    #     sequence = sequence.replace("O","X") # pyrrolysine
    #     sequence = sequence.replace("U","X") # selenocysteine
    #     sequence = sequence.replace(".","X") # selenocysteine
    cleaned_seq = re.sub(pymod_vars.non_prot_standard_regex, pymod_vars.modified_residue_one_letter, clean_white_spaces_from_input(sequence))
    return cleaned_seq.upper()

def clean_white_spaces_from_input(string):
    return string.replace('\n','').replace('\r','').replace(' ','').replace('\t','')

def check_correct_sequence(sequence, remove_indels=True):
    """
    Check if string contains any characters not belonging to the standard protein residues alphabet
    (plus the 'X' characters for heteroresidues.)
    """
    non_protein_character_found = False
    if remove_indels:
        sequence = sequence.replace("-", "")
    for c in sequence:
        if not c in pymod_vars.prot_standard_and_x_one_letter:
            non_protein_character_found = True
    if non_protein_character_found:
        return False
    else:
        return True

def get_invalid_characters_list(sequence, remove_indels=True):
    """
    Check if string contains any characters not belonging to the standard protein residues alphabet
    (plus the 'X' characters for heteroresidues.)
    """
    non_protein_character_found = False
    invalid_characters_list = []
    if remove_indels:
        sequence = sequence.replace("-","")
    for c in sequence:
        if not c in pymod_vars.prot_standard_and_x_one_letter:
            invalid_characters_list.append(c)
    return invalid_characters_list


###################################################################################################
# Classes for standard bioinformatics operations.                                                 #
###################################################################################################

def one2three(letter):
    """
    Returns a three letter code for a residue corresponding to a one letter symbol.
    """
    if letter in pymod_vars.prot_one_to_three_code:
        return pymod_vars.prot_one_to_three_code[letter]
    else:
        return "???"

def three2one(res, force_standard_parent=False):
    """
    Returns a one letter symbol corresponding to a three letter code for a residue. If
    'force_standard_parent' is set to 'True', if the three letter code of a modified residue
    is supplied, the method will attempt to return the code of the corresponding unmodified
    residue.
    """
    # Try to set modified residues of the protein with the letter belonging to the unmodified
    # residue
    if not force_standard_parent:
        # If it is any kind of known modified amminoacid set it as its original non
        # modified amminoacid
        if res in pymod_vars.code_standard:
            return pymod_vars.code_standard[res]
        else:
            return "X"
    else:
        pass
