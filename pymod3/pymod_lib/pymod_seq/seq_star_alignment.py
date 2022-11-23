# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing functions to merge alignments and simplified center star alignments.
"""

import re

from pymod_lib.pymod_vars import dict_of_matrices

import Bio.pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# matrix = matlist.blosum62.copy()

matrix = dict_of_matrices["BLOSUM62"]
matrix.update({("X", "X"): 5})
[matrix.update({("-", i): -10}) for i in "QWERTYIPASDFGHKLXCVNM"]
[matrix.update({(i, "-"): -10}) for i in "QWERTYIPASDFGHKLXCVNM"]


def _get_matrix_val(ca, cb):
    if (ca, cb) in matrix:
        return matrix[(ca, cb)]
    elif (cb, ca) in matrix:
        return matrix[(cb, ca)]
    else:
        return 0

def global_pairwise_alignment_cs(seq1, seq2, gop=10, gep=0.2):

    # Function for scoring aligned residues.
    sf = _get_matrix_val
    # Preparing this object for affine gap penalty scoring greatly speeds up the alignment execution.
    pf = Bio.pairwise2.affine_penalty(open=-gop, extend=-gep)
    # Actually performs the alignment.
    ali = Bio.pairwise2.align.globalcc(seq1, seq2, sf, pf, pf)

    return ali[0][0], ali[0][1]


def build_center_star_alignment(seqs, pairwise_alis=None, verbose=False):

    n_tot_alis = int(len(seqs)*(len(seqs)-1)/2)

    # Use the caller-supplied pairwise alignments.
    if pairwise_alis != None:
        if not len(pairwise_alis) == n_tot_alis:
            raise ValueError("Provided %s parwise alignments, but %s were needed." % (len(pairwise_alis), n_tot_alis))
    # Build from scratch the pairwise alignments.
    else:
        pairwise_alis = []
        for seq_i_id, seq_i in enumerate(seqs):
            for seq_j_id, seq_j in enumerate(seqs[seq_i_id:]):
                if seq_j_id == 0:
                    continue
                if verbose:
                    print("*", seq_i, seq_j)
                aliseq_i, aliseq_j = global_pairwise_alignment_cs(seq_i, seq_j)
                pairwise_alis.append((aliseq_i, aliseq_j))


    # Prepares the center star blocks.
    cs_seq = seqs[0].replace("-", "")
    cs_blocks = ["."]
    for aa in cs_seq:
        cs_blocks.append(aa)
        cs_blocks.append(".")

    blocks_list = [cs_blocks]


    #----------------------------------------------------------------------------------------------
    # Build the blocks.                                                                           -
    #----------------------------------------------------------------------------------------------

    n_cs_alis = len(seqs)-1
    for pairwise_ali_i, pairwise_ali in enumerate(pairwise_alis[:n_cs_alis]):

        os_seq = pairwise_ali[1].replace("-", "")
        os_seq_p = seqs[pairwise_ali_i+1].replace("-", "")
        if not os_seq == os_seq_p:
            raise ValueError("Sequences mismatch: %s/%s." % (os_seq, os_seq_p))

        matches_dict = {}
        os_insertions_dict = {}
        cs_res_count = 0
        os_res_count = 0
        for cs_p, os_p in zip(*pairwise_ali):

            if verbose:
                print(cs_p, os_p)

            # Match.
            if cs_p != "-" and os_p != "-":
                matches_dict[cs_res_count] = os_res_count

            # Insertion in the center star.
            elif cs_p != "-" and os_p == "-":
                pass

            # Insertion in the other sequence.
            elif cs_p == "-" and os_p != "-":
                if cs_res_count in os_insertions_dict.keys():
                    os_insertions_dict[cs_res_count].append(os_res_count)
                else:
                    os_insertions_dict[cs_res_count] = [os_res_count]

            # Double gap.
            elif cs_p == "-" and os_p == "-":
                raise ValueError("Gap only column encountered in alignment: %s" % repr(pairwise_ali))

            if cs_p != "-":
                cs_res_count += 1
            if os_p != "-":
                os_res_count += 1

        if verbose:
            print(matches_dict)
            print(os_insertions_dict)

        # Recontructs the other sequence blocks.
        os_blocks = ["."]
        for aa in cs_seq:
            os_blocks.append("-")
            os_blocks.append(".")

        cs_block_count = 0
        for aa_i, aa in enumerate(cs_seq+"-"): # Adds a gap character to take into C-terminal insertions.
            if aa_i in os_insertions_dict:
                os_blocks[cs_block_count] = "".join([os_seq[i] for i in os_insertions_dict[aa_i]])
            cs_block_count += 1

            if aa_i in matches_dict:
                os_blocks[cs_block_count] = os_seq[matches_dict[aa_i]]
            cs_block_count += 1

        if verbose:
            print(os_blocks)
            print(cs_blocks)
            print("")
        blocks_list.append(os_blocks)


    if verbose:
        print("\n# Results")


    #----------------------------------------------------------------------------------------------
    # Build the final alignment.                                                                  -
    #----------------------------------------------------------------------------------------------

    ali_seqs = [[] for i in range(0, len(seqs))]
    for block_index in range(0, len(cs_blocks)):

        # Add the matches.
        if verbose:
            print(cs_blocks[block_index], [b[block_index] for b in blocks_list])

        # Add insertions.
        if cs_blocks[block_index] == ".":
            all_insertions_dict = []
            res_insertions_ids = []

            for b_i, b in enumerate(blocks_list):
                if b[block_index] != ".":
                    all_insertions_dict.append(b[block_index])
                    res_insertions_ids.append(b_i)
                else:
                    all_insertions_dict.append("-")

            if len(res_insertions_ids) > 0:
                max_insertion_len = max([len(ins) for ins in all_insertions_dict])

                # Insertions are short, no need to refine them.
                if max_insertion_len <= 3 or len(res_insertions_ids) == 1:
                    pass

                # Refine an insertion block.
                else:
                    refine_seqs = [all_insertions_dict[ins_id] for ins_id in res_insertions_ids]
                    if verbose:
                        print("\n= Refining:", refine_seqs)

                    # Performs a center star alignment with the inserted fragments in order to
                    # align them.
                    refine_seqs = build_center_star_alignment(refine_seqs)
                    for ins_id_i, ins_id in enumerate(res_insertions_ids):
                        all_insertions_dict[ins_id] = refine_seqs[ins_id_i]
                    max_insertion_len = max([len(ins) for ins in all_insertions_dict])

                # Actually include the insertions in the alignment.
                for s, ins in zip(ali_seqs, all_insertions_dict):
                    s.append(ins.ljust(max_insertion_len, "-"))

        else:
            for s, b in zip(ali_seqs, blocks_list):
                s.append(b[block_index])

    ali_seqs = ["".join(s) for s in ali_seqs]
    if verbose:
        print("")
        for ali_seq in ali_seqs:
            print(ali_seq)

    return ali_seqs


def _all_positions_are_gaps(pos_list):
    for p in pos_list:
        if p != "-":
            return False
    return True

def save_cstar_alignment(seqs, all_ids, pairwise_alis, output_filepath=None, verbose=False):
    """
    Performs a center star alignment between a list of sequences.

    seqs:
        must contain a list of gapless sequences. The first one is assumed to be the center star.
    all_ids:
        must contain a list of ids associated to the sequences in 'seqs'.
    pairwise_alis:
        must contain a list of n*(n-1)/2 elements, where n=len(seqs). Eache element
        represents an alignment between two sequences. For example, if we have
            seqs = ["A", "B", "C", "D", "E"]
        'pairwise_alis' must contain the following tuples, each corresponding to a
        pairwise alignment between two sequences:
            pairwise_alis = [("A", "B"), ("A", "C"), ("A", "D"), ("A", "E"),
                             ("B", "C"), ("B", "D"), ("B", "E"),
                             ("C", "D"), ("C", "E"),
                             ("D", "E")]
        Instead of providing tuples, a 'None' element can be provided for alignments in which the
        center star is not present (in the examples above, all the tuples where "A" is not
        present).
    output_filepath:
        alignment filepath where the center star alignment will be written (in clustal
        format).
    """

    # Remove gap-only columns from pairwise alignments, if necessary.
    if pairwise_alis != None:
        for pairwise_ali_id, pairwise_ali in enumerate(pairwise_alis[:]):
            if pairwise_ali != None:
                if len(pairwise_ali[0]) != len(pairwise_ali[1]):
                    raise ValueError("Sequences are not aligned: \n%s\n%s" % (pairwise_ali[0], pairwise_ali[1]))
                columns_to_remove = []
                for i in range(0, min([len(s) for s in pairwise_ali])):
                    if _all_positions_are_gaps([s[i] for s in pairwise_ali]):
                        columns_to_remove.append(i)
                if len(columns_to_remove) > 0:
                    list_seqs = [list(s) for s in pairwise_ali]
                    for columns_removed, i in enumerate(columns_to_remove):
                        for s in list_seqs:
                            s.pop(i-columns_removed)
                    pairwise_alis[pairwise_ali_id] = ("".join(list_seqs[0]), "".join(list_seqs[1]))

    # Actually builds the center star alignment.
    msa_0 = build_center_star_alignment(seqs=seqs, pairwise_alis=pairwise_alis, verbose=verbose)
    max_length = max([len(si) for si in msa_0])
    msa_0 = [si.ljust(max_length, "-") for si in msa_0]

    # Saves the alignment.
    seq_records = []
    for aliseq, rec_id in zip(msa_0, all_ids):
        seq_records.append(SeqRecord(Seq(str(aliseq)), id=rec_id))
    if output_filepath != None:
        SeqIO.write(seq_records, output_filepath, "clustal")

    return seq_records


if __name__ == "__main__":

    pairwise_alis = [("AAAA-",
                      "BBB-B"),
                     ("AA----AA",
                      "-CCCCCCC"),
                     ("AA---AA-",
                      "--CCADDD"),
                    #  ("A-A-A-A",
                    #   "EEEEEE-"),
                     None,
                     None,
                     None,
                    #  None,
                    #  None,
                    #  None
                     ]

    build_center_star_alignment(
                                # seqs=["AAAA", "BBBB", "CCCCC", "DDD", "EEEEEE"],
                                seqs=["AAAA", "BBBB", "CCCCCCC", "CCADDD", ],
                                pairwise_alis=pairwise_alis,
                                verbose=True,
                                )


###################################################################################################
# Join alignments.                                                                                #
###################################################################################################

def _get_ali_dicts(msa, verbose=False):
    msa_1_seq_1 = msa[0]
    match_cols_dict_1 = {}
    insertion_cols_dict_1 = {}
    res_1_count = 0

    in_insertion = False
    for ai in range(0, len(msa_1_seq_1)):

        if verbose:
            print([s[ai] for s in msa])

        if msa_1_seq_1[ai] != "-":
            match_cols_dict_1[res_1_count] = [s[ai] for s in msa[1:]]
            res_1_count += 1
            in_insertion = False

        if msa_1_seq_1[ai] == "-":
            if in_insertion:
                insertion_cols_dict_1[res_1_count].append([s[ai] for s in msa])
            else:
                insertion_cols_dict_1[res_1_count] = [[s[ai] for s in msa]]
                in_insertion = True

    return match_cols_dict_1, insertion_cols_dict_1


def _get_complete_block(j_seq, msa, match_cols_dict, insertion_cols_dict):
    anchor_block = ["."]
    p_count = 0

    for p1 in j_seq:

        if p1 != "-":
            if p_count in insertion_cols_dict:
                anchor_block[-1] = insertion_cols_dict[p_count]
            anchor_block.append([p1] + match_cols_dict[p_count])
            p_count += 1
        else:
            anchor_block.append(["-"]*len(msa))

        anchor_block.append(".")

    if p_count in insertion_cols_dict:
        anchor_block[-1] = insertion_cols_dict[p_count]

    return anchor_block


def join_alignments(j_msa, msa_l, verbose=False):
    """
    Join two multiple sequence alignments provided in 'msa_l' using an "anchor" pairwise alignment
    provided in 'j_msa'.

    j_msa:
        must contain a pairwise alignment between two "anchor" sequences. For example:
            j_msa = ["AAA-AA", "DD-DDD"]
    msa_l:
        must contain list with two multiple sequence alignments in which the first sequences must
        be the anchor sequences. For example:
            msa_l = [
                     ["AAAA-A", "-BBBBB", "CCC-CC"],
                     ["-DD-DDD", "E-EEEE-", "F-FF-FF"],
                    ]
        Note how in the first multiple sequence alignment, the first sequence is the first anchor
        in 'j_msa' (that is, "AAAAA"). In the second multiple sequence alignment, the first
        sequence is the second anchor in 'j_msa' (that is, "DDDDD").

    return new_msa:
        Returns a list with all the sequences present in 'msa_l' aligned according to the new
        alignment.
    """

    if len(j_msa) > 2:
        raise ValueError("The anchor alignment can have only two sequences (%s provided)." % len(j_msa))

    for a_i, a in enumerate(j_msa):
        a_seq = a.replace("-", "")
        msa_seq = msa_l[a_i][0].replace("-", "")
        if not a_seq == msa_seq:
            raise ValueError("The anchor sequences provided in 'j_msa' and 'msa_l' are not the same: \n%s\n%s." % (a_seq, msa_seq))


    match_cols_dict_1, insertion_cols_dict_1 = _get_ali_dicts(msa_l[0], verbose=False)
    match_cols_dict_2, insertion_cols_dict_2 = _get_ali_dicts(msa_l[1], verbose=False)
    if verbose:
        print(match_cols_dict_1, insertion_cols_dict_1)
        print(match_cols_dict_2, insertion_cols_dict_2)
        print("")
        print(j_msa[0])
        print(j_msa[1])
        print("")


    # Prepares the anchor blocks.
    anchor_block_1 = _get_complete_block(j_msa[0], msa_l[0], match_cols_dict_1, insertion_cols_dict_1)
    # print (anchor_block_1, len(anchor_block_1))
    anchor_block_2 = _get_complete_block(j_msa[1], msa_l[1], match_cols_dict_2, insertion_cols_dict_2)
    # print (anchor_block_2, len(anchor_block_2))

    if verbose:
        print("\n###")
        for a in msa_l[0]:
            print(a)
            print("")
        print("\n###")
        for a in msa_l[1]:
            print(a)
        print("")


    # Actually rebuilts the full alignment.
    alignment_matrix = []
    for b1, b2 in zip(anchor_block_1, anchor_block_2):

        if verbose:
            print(b1, b2)

        if b1 == "." and b2 == ".":
            pass

        elif b1 == "." and b2 != ".":
            for si in b2:
                alignment_column = ["-"]*len(msa_l[0]) + si
                alignment_matrix.append(alignment_column)

        elif b1 != "." and b2 == ".":
            for si in b1:
                alignment_column = si + ["-"]*len(msa_l[1])
                alignment_matrix.append(alignment_column)

        elif b1 != "-" and b2 != "-":
            if type(b1[0]) != list and type(b2[0]) != list:
                alignment_column = b1 + b2
                alignment_matrix.append(alignment_column)
            else:
                if type(b1[0]) == list:
                    for si in b1:
                        alignment_column = si + ["-"]*len(msa_l[1])
                        alignment_matrix.append(alignment_column)
                if type(b2[0]) == list:
                    for si in b2:
                        alignment_column = ["-"]*len(msa_l[0]) + si
                        alignment_matrix.append(alignment_column)


    new_msa = [[] for i in range(0, len(msa_l[0]) + len(msa_l[1]) )]

    for c in alignment_matrix:
        for ci, s in zip(c, new_msa):
            s.append(ci)

    if verbose:
        print("\n@ Results")
        print(alignment_matrix)
        for s in new_msa:
            print("".join(s))

    new_msa = ["".join(s) for s in new_msa]

    return new_msa


if __name__ == "__main__":

    j_msa = ("AAAAAA-A",
             "-DDDD-DD")

    msa_l = [("AAAAA-AA",
              "B-BBBBB-",
              "CCCCCCC-"
              ),
             ("-DD-D--DDD",
              "EEEE-EE---")]

    results = join_alignments(j_msa=j_msa, msa_l=msa_l, verbose=True)

    print("")
    for s in results:
        print(s)
