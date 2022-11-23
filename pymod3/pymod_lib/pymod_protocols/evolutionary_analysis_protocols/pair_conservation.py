# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Protocol to compute the conservation scores between a reference sequence in a multiple
sequence alignment and the rest of the sequences in that alignment.
"""


import os
import math

from pymod_lib.pymod_vars import dict_of_matrices

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_radioselect_qt, PyMod_combobox_qt)
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols._evolutionary_analysis_base import Evolutionary_analysis_protocol


class Pair_conservation_analysis(Evolutionary_analysis_protocol):
    """
    Class implementing a PyMod protocol to score the conservation between pair of
    aligned residues of a multiple sequence alignment.
    """

    def launch_from_gui(self):
        self.build_pc_window()


    def build_pc_window(self):
        """
        Builds a window with options for the pair conservation scorer algorithm.
        """

        self.pc_scorer_window = Pair_conservation_window_qt(parent=self.pymod.main_window,
            protocol=self,
            title="Pair Conservation options",
            upper_frame_title="Here you can modify options for Pair Conservation analysis",
            submit_command=self.pc_scorer_state)
        self.pc_scorer_window.show()


    def pc_scorer_state(self):
        """
        Called when the "SUBMIT" button is pressed options window. Contains the
        the code to compute the pairwise conservation scores.
        """

        try:
            # Get the options from the GUI.
            self.reference_id = self.pc_scorer_window.get_reference_id()
            self.conservation_mode = "blosum62"

            # Gets the list of conservation scores. There are as many values as
            # positions in the alignment.
            all_pc_scores = self.pair_conservation_scorer()

            # Assigns conservation scores to each one of the aligned sequences.
            for seq_pc_scores, seq in zip(all_pc_scores, self.input_cluster_element.get_children()):
                residues = seq.get_polymer_residues()
                rc = 0
                for (r, v) in zip(seq.my_sequence, seq_pc_scores):
                    if r != "-":
                        residues[rc].pc_score = v
                        rc += 1
                seq._has_pc_scores = True
                self.pymod.main_window.color_element_by_pc_scores(seq)

        except Exception as e:
            message = "Could not compute Pair Conservation scores because of the following error: '%s'." % e
            self.pymod.main_window.show_error_message("Pair Conservation Error", message)

        # Removes the temporary alignment file.
        self.pc_scorer_window.close()


    # amino_acids = tuple("QWERTYIPASDFGHKLCVNM")
    # amino_acids_and_gap = tuple("QWERTYIPASDFGHKLCVNM-")

    def pair_conservation_scorer(self):
        """
        Computes the pair conservation scores. Get the reference sequence, and
        for each of its residues, check if the corresponding residue in another
        sequence is conserved (according to different criteria).
        """

        # Prepare the scoring matrix.
        if self.conservation_mode == "blosum62":
            #sub_matrix = MatrixInfo.blosum62.copy()
            sub_matrix = dict_of_matrices["BLOSUM62"]
            for pair in list(sub_matrix.keys()):
                value = sub_matrix[pair]
                sub_matrix.update({(pair[1], pair[0]): value})
        else:
            raise KeyError("Unknown 'conservation_mode': " % self.conservation_mode)

        # Get the elements and check their sequences.
        aligned_elements = self.input_cluster_element.get_children()
        reference_element = aligned_elements[self.reference_id]
        if len(set([len(seq.my_sequence) for seq in aligned_elements])) != 1:
            raise ValueError("Not all sequences in the alignment have the same length.")

        # Actually assign the scores.
        ali_len = len(aligned_elements[0].my_sequence)
        all_pc_scores = []

        for aligned_seq in aligned_elements:

            seq_pc_scores = []

            for i in range(0, ali_len):

                ref_pos = reference_element.my_sequence[i]
                sel_pos = aligned_seq.my_sequence[i]

                if ref_pos != "-" and sel_pos != "-":
                    # Modify this to include other conservation scores.
                    if ref_pos == sel_pos:
                        seq_pc_scores.append(2)
                    else:
                        score = sub_matrix.get((ref_pos, sel_pos), 0)
                        if score > 0:
                            seq_pc_scores.append(1)
                        else:
                            seq_pc_scores.append(0)
                else:
                    seq_pc_scores.append(0)

            all_pc_scores.append(seq_pc_scores)

        return all_pc_scores


class Pair_conservation_window_qt(PyMod_protocol_window_qt):
    """
    Options window for MSA pair conservation calculations.
    """

    def add_middle_frame_widgets(self):

        # Type of conservation measure.
        # TODO.
        pass

        # Reference sequence.
        sequences_list = [element.my_header for element in self.protocol.input_cluster_element.get_children()]
        self.reference_combobox = PyMod_combobox_qt(label_text="Reference Sequence",
                                                    items=sequences_list)
        self.reference_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.reference_combobox)

        self.middle_formlayout.set_input_widgets_width(175)


    def get_reference_id(self):
        try:
            return self.reference_combobox.get_index()
        except:
            print("- Warning: could not obtain the reference sequence id.")
            return 0
