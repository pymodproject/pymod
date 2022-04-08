# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Protocol to compute Shannon's entropy or relative entropy of a multiple sequence alignment loaded
in PyMod.
"""


import os
import math

import numpy as np
from Bio import AlignIO

from pymod_lib.pymod_vars import yesno_dict
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_protocol_window_qt, PyMod_radioselect_qt
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols._evolutionary_analysis_base import Evolutionary_analysis_protocol


# Amino acid frequencies found in the UniProt/SwissProt database (retrieved from:
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
# on 31/12/2019). Used to compute the relative entropy of a multiple sequence alignment.
swissprot_freq_dict = {'Q': 0.03932447341942903, 'W': 0.010987376799787335, 'E': 0.06731857081478208,
                       'R': 0.055339639013144266, 'T': 0.053554782580885245, 'Y': 0.029203273379595142,
                       'I': 0.05922496839490527, 'P': 0.047360775631143305, 'A': 0.08258683854220726,
                       'S': 0.06632123826472652, 'D': 0.054624431023512963, 'F': 0.03865735720637375,
                       'G': 0.07081179428789602, 'H': 0.02275979219540371, 'K': 0.05814545801909458,
                       'L': 0.09654911798449017, 'C': 0.013828227451923757, 'V': 0.06864307992926043,
                       'N': 0.04060773272101554, 'M': 0.024151072340423636}


class Entropy_analysis(Evolutionary_analysis_protocol):
    """
    Class implementing a PyMod protocol to score the entropy of a multiple sequence alignment.
    """

    def launch_from_gui(self):
        self.build_entropy_scorer_window()


    entropy_func_labels = ["Shannon's Entropy", "Relative Entropy"]
    entropy_func_names = ["entropy", "relative_entropy"]
    entropy_func_dict = dict(zip(entropy_func_labels, entropy_func_names))

    def build_entropy_scorer_window(self):
        """
        Builds a window with options for the entropy scorer algorithm.
        """

        self.entropy_scorer_window = Entropy_scorer_window_qt(parent=self.pymod.main_window,
            protocol=self,
            title="Entropy analysis options",
            upper_frame_title="Here you can modify options for Entropy analysis",
            submit_command=self.entropy_scorer_state)
        self.entropy_scorer_window.show()


    def entropy_scorer_state(self):
        """
        Called when the "SUBMIT" button is pressed on the Entropy analysis window. Contains
        the code to compute Shannon's entropies.
        """

        # Saves a .fasta file for the alignment.
        aligned_sequences = self.input_cluster_element.get_children()
        self.pymod.save_alignment_fasta_file("temp", aligned_sequences)
        input_file_shortcut = os.path.join(self.pymod.alignments_dirpath, "temp.fasta")

        try:

            # Get the options from the GUI.
            toss_gaps = yesno_dict[self.entropy_scorer_window.entropy_scorer_exclude_gaps_rds.getvalue()]
            toss_gap_threshold = 0.25
            entropy_func = self.entropy_func_dict[self.entropy_scorer_window.entropy_scorer_function_rds.getvalue()]

            # Gets the list of entropy score. There are as many values as positions in the alignment.
            entropy_scores = self.score_entropy(input_file_shortcut,
                                                toss_gaps=toss_gaps,
                                                toss_gap_threshold=toss_gap_threshold,
                                                entropy_func=entropy_func)

            # Filter the entropy scores.
            valid_entropy_scores = [i for i in entropy_scores if i is not None]
            if len(valid_entropy_scores) == 0:
                raise ValueError("All the colums have high levels of gaps (> %s)." % toss_gap_threshold)

            # Assigns bins to the entropy values.
            entropy_items = self.get_entropy_items_list(entropy_scores, entropy_func=entropy_func)

            # Assigns entropy scores to each one of the aligned sequences.
            for seq in aligned_sequences:
                residues = seq.get_polymer_residues()
                rc = 0
                for (r, v) in zip(seq.my_sequence, entropy_items):
                    if r != "-":
                        residues[rc].entropy_score = v
                        rc += 1
                seq._has_entropy_scores = True
                self.pymod.main_window.color_element_by_entropy_scores(seq)

        except Exception as e:
            message = "Could not compute Entropy scores because of the following error: '%s'." % e
            self.pymod.main_window.show_error_message("Entropy Analysis Error", message)

        # Removes the temporary alignment file.
        os.remove(input_file_shortcut)
        self.entropy_scorer_window.destroy()


    amino_acids = tuple("QWERTYIPASDFGHKLCVNM")
    # amino_acids_and_gap = tuple("QWERTYIPASDFGHKLCVNM-")

    def score_entropy(self, alignment_filepath,
                      toss_gaps=False, toss_gap_threshold=0.25,
                      entropy_func="entropy",
                      only_aa=True):
        """
        Reads the alignment provided in the 'alignment_filepath' argument and computes the Shannon's
        entropy of each column in the alignment. The alignment must be in the FASTA format.
        """

        if not entropy_func in self.entropy_func_names:
            raise KeyError("Unknown 'entropy_func': %s" % entropy_func)

        alignment = list(AlignIO.parse(alignment_filepath, "fasta"))[0]
        if len(set([len(seq.seq) for seq in alignment])) != 1:
            raise ValueError("Not all sequences in the alignment have the same length.")

        ali_len = len(alignment[0].seq)

        # Iterate through every column of the alignment.
        entropy_scores = []
        for i in range(0, ali_len):

            ali_col_i = list(alignment[:,i])

            # Exclude columns with too many gaps.
            if toss_gaps:
                gap_freq = ali_col_i.count("-")/float(len(ali_col_i))
                if gap_freq >= toss_gap_threshold: # There are too many gaps.
                    entropy_scores.append(None)
                    continue

            # Exclude gap characters from the entropy calculations.
            if only_aa:
                ali_col_i = [p_ij for p_ij in ali_col_i if p_ij != "-"]

            # There is a gap-only column.
            if len(ali_col_i) == 0:
                entropy_scores.append(None)
                continue

            if entropy_func == "entropy":
                entropy_i = sum([self.get_information(aa, ali_col_i) for aa in self.amino_acids])
            elif entropy_func == "relative_entropy":
                entropy_i = sum([self.get_relative_information(aa, ali_col_i) for aa in self.amino_acids])
            entropy_scores.append(entropy_i)

        return entropy_scores


    def get_information(self, aa, ali_col, pseudo_count=0):
        aa_counts = ali_col.count(aa)
        if aa_counts == 0:
            return 0
        else:
            aa_freq = (aa_counts + pseudo_count)/float(len(ali_col))
            return -aa_freq*math.log(aa_freq)

    def get_relative_information(self, aa, ali_col, pseudo_count=0):
        aa_counts = ali_col.count(aa)
        if aa_counts == 0:
            return 0
        else:
            aa_freq = (aa_counts + pseudo_count)/float(len(ali_col))
            return aa_freq*math.log(aa_freq/swissprot_freq_dict[aa])


    n_bins = 9

    def get_entropy_items_list(self, entropy_scores, entropy_func="entropy"):
        """
        Prepares dictionaries which will store the entropy score and bin for each residue.
        """

        # Tosses alignment positions with too many gaps.
        clist = np.array([i for i in entropy_scores if i is not None])
        if entropy_func == "relative_entropy":
            clist = -clist

        bins = np.linspace(min(clist), max(clist), num=self.n_bins+1)

        list_of_entropy_items = []
        for pos_idx, val in enumerate(entropy_scores):
            # Adds 'None' values for position tossed out because of their high gap content.
            if val is None:
                list_of_entropy_items.append({"entropy-score": None, "interval": None})
            else:
                if entropy_func == "entropy":
                    b = self._get_bin(np.digitize(val, bins, right=True))
                elif entropy_func == "relative_entropy":
                    b = self._get_bin(np.digitize(-val, bins, right=True))
                list_of_entropy_items.append({"entropy-score": round(val, 3),
                                              "interval": b})

        return list_of_entropy_items

    def _get_bin(self, b):
        if b <= 0:
            return 1
        else:
            return b


class Entropy_scorer_window_qt(PyMod_protocol_window_qt):
    """
    Options window for MSA entropy calculations.
    """

    def add_middle_frame_widgets(self):
        # Type of entropy measure.
        self.entropy_scorer_function_rds = PyMod_radioselect_qt(label_text="Entropy Measure",
                                                                buttons=self.protocol.entropy_func_labels)
        self.entropy_scorer_function_rds.setvalue(self.protocol.entropy_func_labels[0])
        self.middle_formlayout.add_widget_to_align(self.entropy_scorer_function_rds)

        # Toss gaps.
        self.entropy_scorer_exclude_gaps_rds = PyMod_radioselect_qt(label_text="Toss gaps",
                                                                    buttons=('Yes', 'No'))
        self.entropy_scorer_exclude_gaps_rds.setvalue("Yes")
        self.middle_formlayout.add_widget_to_align(self.entropy_scorer_exclude_gaps_rds)

        self.middle_formlayout.set_input_widgets_width("auto", padding=10)
