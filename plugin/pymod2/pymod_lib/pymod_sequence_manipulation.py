import Bio
from Bio.pairwise2 import format_alignment
import pymod_vars as pmdt

###################################################################################################
# Some useful global functions used by all classes.                                               #
###################################################################################################

def check_correct_sequence(sequence, remove_indels=True):
    """
    Check if string contains any characters not belonging to the standard protein residues alphabet
    (plus the 'X' characters for heteroresidues.)
    """
    non_protein_character_found = False
    if remove_indels:
        sequence = sequence.replace("-","")
    for c in sequence:
        if not c in pmdt.protein_residues_set:
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
        if not c in pmdt.protein_residues_set:
            invalid_characters_list.append(c)
    return invalid_characters_list


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


def get_residue_id_in_aligned_sequence(aligned_sequence,real_id):
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


def get_residue_id_in_gapless_sequence(aligned_sequence,id_in_alignment):
    """
    This does the opposite of the function above.
    aligned_sequence: a sequence with indels.
    id_in_alignment: the id (position -1) of the same residue in the alignment (the sequence + indels).
    returns: the id (position -1) of the residue in the gapless sequence.
    """
    assert(len(aligned_sequence) >= id_in_alignment)
    residue_counter = 0
    real_id = 0
    for (k,r) in enumerate(aligned_sequence):
        if k != id_in_alignment:
            if r != "-":
                residue_counter += 1
        else:
            real_id = residue_counter
    return real_id


def find_residue_conservation(template_seq,target_seq,real_id):
    """
    Takes two aligned sequences, a "template" and a "target", and sees if a
    residue with the id "i" of the template is conserved in the target.
    """
    conservation = False
    # seq_counter: real id of the residue
    # k: id of the residue in the alignment
    alignment_position = get_residue_id_in_aligned_sequence(template_seq,real_id)
    if template_seq[alignment_position] == target_seq[alignment_position]:
        conservation = True
    else:
        conservation = False
    return conservation


def get_leading_gaps(sequence, index):
    gap_counter = 0
    for i, position in enumerate(sequence[index:]):
        if i != 0:
            if position == "-":
                gap_counter += 1
            else:
                break
    return gap_counter


def get_starting_gaps(sequence):
    """
    Gets the number of starting gaps in an aligned sequence. For '---ALKNMFG-R'
    returns 3.
    """
    gc = 0
    for p in sequence:
        if p == "-":
            gc += 1
        else:
            break
    return gc


def count_indels_to_next_residue(seq, start_real_id):
    """
    Takes an aligned sequence, and counts how many indels there are from one residues (specified
    by "start_real_id") to the next. For example for the sequence "A-TLSKRC" start_real_id=0,
    this will return a dictionary with:
        "indels-to-next-residue" : 1 (there is one indel from the first A to the second T)
        "start-aligned-id" : 0 (id of the first A in the aligned sequence)
        "end-aligned-id" : 2 (id of the second T in the aligned sequence)
    """
    start_aligned_id = get_residue_id_in_aligned_sequence(seq,start_real_id)
    for i,p in enumerate(seq[start_aligned_id:]):
        if i != 0:
            if p != "-":
                results = {"indels-to-next-residue": i-1,
                           "start-aligned-id": start_aligned_id,
                           "end-aligned-id": start_aligned_id+i}
                return results
    # This will be returned for the last position.
    return {"indels-to-next-residue": 0,
             "start-aligned-id": start_aligned_id,
             "end-aligned-id": start_aligned_id+1}


def global_pairwise_alignment(seq1, seq2, toss_modres=False):
    ali = Bio.pairwise2.align.globalms(seq1, seq2, 2, -1, -.5, -.1)
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

###################################################################################################
# Classes needed to join alignments and to perform a center star alignment.                       #
###################################################################################################

class Star_alignment:
    """
    A class used to build star alignments. This will be used when building alignments comprising
    the query and the retrieved hits. Or also to append new sequences to an already existing
    alignment by using as guide (that is, a 'star') a sequence in the alignment.
    """

    def __init__(self, query):
        """
        The query is just the 'star' of star alignment. This is called 'query' because this class
        will mostly be used to build alignments in which a (PSI)-BLAST query is the star.
        """
        self.query_str = query
        self.query = list(query)
        self.original_query = self.query[:]
        self.query_length = len(filter(lambda x: x != "-", self.query))
        self.query_starting_indels = get_starting_gaps(self.original_query)

        # A list comprising the query (as its first element) and all the other sequences of the
        # alignment.
        self.aligned_sequences = []
        self.aligned_sequences.append(self.query)


    def extend_with_aligned_sequences(self, aligned_sequences):
        """
        Builds a list of sequences already aligned to the query to populate the 'aligned_sequences'
        list. This is used to reconstruct an already existing cluster before appending new
        sequences to it.
        """
        for seq in aligned_sequences:
            self.aligned_sequences.append(list(seq))


    def build_blast_local_alignment_list(self, hsp_list):
        """
        Takes as input a list of HSP class objects from Biopython. In this way a star alignment can
        be generated through the 'generate_blast_pseudo_alignment()' method.
        """
        self.hsp_list = hsp_list
        # This will contain a list of 'Pairwise_alignment' objects.
        self.blast_local_alignments_list = []
        for hsp in self.hsp_list:
            local_alignment = Pairwise_alignment(str(hsp.query),str(hsp.sbjct),query_start= hsp.query_start, full_query=self.query_str)
            self.blast_local_alignments_list.append(local_alignment)


    def generate_blast_pseudo_alignment(self):
        """
        Builds a star alignment using as a center self.query, stored inside
        self.aligned_sequences[0].
        """
        for j,local_alignment in enumerate(self.blast_local_alignments_list):
            self.append_new_sequence(local_alignment)


    def append_new_sequence(self, pairwise_alignment):
        """
        Takes as input a 'Pairwise_alignment' object and appends it to the cluster of sequences
        present in 'aligned_sequences' using as a guide the sequence in 'aligned_sequences[0]' and
        by keeping the other sequences in frame.
        """
        # Adds a new line to the multiple alignment.
        self.add_new_sequence()

        if self.query_starting_indels != 0:
            self.add_residues_to_new_sequence(["-"]*(self.query_starting_indels))

        counter = 0
        while (counter < self.query_length):

            fq_results = count_indels_to_next_residue(self.query, counter)
            aq_results = count_indels_to_next_residue(pairwise_alignment.aligned_query, counter)
            aq_start = aq_results["start-aligned-id"]
            aq_end = aq_results["end-aligned-id"]
            fq_start = fq_results["start-aligned-id"]
            fq_end = fq_results["end-aligned-id"]

            if aq_results["indels-to-next-residue"] == fq_results["indels-to-next-residue"]:
                residues_to_add = pairwise_alignment.aligned_sequence[aq_start:aq_end]
                self.add_residues_to_new_sequence(residues_to_add)

            elif aq_results["indels-to-next-residue"] > fq_results["indels-to-next-residue"]:
                indels_to_add = aq_results["indels-to-next-residue"] - fq_results["indels-to-next-residue"]
                # Is the order important here?
                self.add_gaps_to_alignment(fq_start, indels_to_add)
                residues_to_add = pairwise_alignment.aligned_sequence[aq_start:aq_end]
                self.add_residues_to_new_sequence(residues_to_add)

            elif aq_results["indels-to-next-residue"] < fq_results["indels-to-next-residue"]:
                indels_to_add = fq_results["indels-to-next-residue"] - aq_results["indels-to-next-residue"]
                residues_to_add = pairwise_alignment.aligned_sequence[aq_start:aq_end]
                self.add_residues_to_new_sequence(residues_to_add)
                self.add_residues_to_new_sequence(["-"]*indels_to_add)

            counter += 1

        # print "# fq:", self.query
        # print "# aq:", pairwise_alignment.aligned_query
        # print "# as:", pairwise_alignment.aligned_sequence
        # print "\n"


    def add_new_sequence(self):
        self.aligned_sequences.append([])


    def add_residues_to_new_sequence(self,residues):
        self.aligned_sequences[-1].extend(residues)


    def add_gaps_to_alignment(self,gap_index, number_of_gaps):
        for seq in self.aligned_sequences[0:-1]:
            seq[gap_index+1:gap_index+1] = ["-"]*(number_of_gaps)


    def update_pymod_elements(self, elements_to_update):
        """
        Updates the sequences of a list of 'PyMod_elements' objects provided in the
        'elements_to_update' argument with the results obtained while building the alignment.
        """
        for updated_sequence, element in zip(self.aligned_sequences, elements_to_update):
            element.my_sequence = "".join(updated_sequence)


class Pairwise_alignment:
    """
    A class to store information about a pairwise alignments. Objects of this class will be used
    inside the "Star_alignment" class in order to build star alignments and append new sequences
    to existing clusters using as guides 'star' sequences.
    """

    def __init__(self, aligned_query, aligned_sequence, query_start=1, full_query=None):
        # Actually the so-called center of the star alignment.
        self.aligned_query = list(aligned_query)
        # The sequence to be added to the star alignment.
        self.aligned_sequence = list(aligned_sequence)
        # The position (real_id+1) of the residue of the full query.
        self.query_start = query_start
        '''
        if full_query != None:
            self.full_query = list(full_query) if full_query != None else self.aligned_query
            self.complete_sequences()
        '''
        if full_query != None:
            self.full_query = list(full_query)
            self.complete_sequences()
        else:
            self.full_query = self.aligned_query

    def complete_sequences(self):
        """
        Pad with indels the two aligned sequences to adjust their lenght to the lenght of the
        full query.
        """
        # Pad on the left.
        aligned_start_id = get_residue_id_in_aligned_sequence(self.full_query, self.query_start-1)
        # Pads with residues from the full query.
        self.aligned_query[0:0] = self.full_query[0:aligned_start_id]
        # Pads with indels.
        self.aligned_sequence[0:0] = ["-"]*(aligned_start_id)

        # Pad on the right.
        gapless_full_query = filter(lambda p:p!="-",self.full_query)
        gapless_aligned_query = filter(lambda p:p!="-",self.aligned_query)
        length_difference = len(gapless_full_query) - len(gapless_aligned_query)
        if length_difference != 0:
            extended_sequence = gapless_full_query[len(gapless_aligned_query):]
            self.aligned_query.extend(extended_sequence)
            self.aligned_sequence.extend(["-"]*len(extended_sequence))
