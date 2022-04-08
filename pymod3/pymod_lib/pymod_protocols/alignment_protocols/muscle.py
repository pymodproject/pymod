# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
MUSCLE.
"""

import os

from Bio.Align.Applications import MuscleCommandline

# Protocols.
from ._base_alignment._base_regular_alignment import Regular_sequence_alignment

# GUI.
from ._base_alignment._gui import Regular_alignment_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_radioselect_qt


class MUSCLE_alignment:

    alignment_program = "muscle"

    def additional_initialization(self):
        self.tool = self.pymod.muscle

    def get_options_from_gui(self):
        self.params_from_gui = {}
        self.params_from_gui["muscle_mode"] = self.alignment_window.get_muscle_mode()
        return True


class MUSCLE_regular_alignment(MUSCLE_alignment, Regular_sequence_alignment):

    protocol_name = "muscle"

    def get_alignment_window_class_qt(self):
        return MUSCLE_regular_window_qt

    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_muscle(sequences_to_align, output_file_name=output_file_name,
                        muscle_mode=self.params_from_gui["muscle_mode"])

    def run_muscle(self, sequences_to_align, output_file_name, muscle_mode):
        """
        This method allows to interact with the local MUSCLE.
        """
        # TODO: to insert the following options:
        #           - guide tree from:
        #               - none
        #               - first iteration
        #               - second iteration
        self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)
        # Input FASTA for MUSCLE.
        infasta=os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
        # Output FASTA from MUSCLE, in tree order.
        outfasta_tree=os.path.join(self.pymod.alignments_dirpath, output_file_name + ".out_fasta")
        # Output ALN.
        outaln=os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
        muscle_exec = self.tool["exe_file_path"].get_value()
        if muscle_mode == "highest_accuracy":
            cline = MuscleCommandline(muscle_exec, input=infasta, out=outfasta_tree, clwout=outaln)
        elif muscle_mode == "large_datasets":
            cline = MuscleCommandline(muscle_exec, input=infasta, out=outfasta_tree, clwout=outaln,
                                      maxiters=2)
        elif muscle_mode == "fastest":
            cline = MuscleCommandline(muscle_exec, input=infasta, out=outfasta_tree, clwout=outaln,
                                      maxiters=1, diags=True, sv=True, distance1="kbit20_3")
        else:
            raise KeyError(muscle_mode)
        self.pymod.execute_subprocess(str(cline))


class MUSCLE_regular_window_qt(Regular_alignment_window_qt):

    def build_algorithm_options_widgets(self):

        # MUSCLE mode radioselect (for more information see: https://www.drive5.com/muscle/manual/).
        self.muscle_modes = ["Highest Accuracy", "Large Datasets", "Fastest"]
        self.muscle_modes_short = ["highest_accuracy", "large_datasets", "fastest"]
        self.muscle_modes_dict = dict((k, v) for (k, v) in zip(self.muscle_modes, self.muscle_modes_short))

        self.muscle_mode_rds = PyMod_radioselect_qt(label_text="MUSCLE Mode",
                                                    buttons=self.muscle_modes)
        self.muscle_mode_rds.setvalue(self.muscle_modes[0])
        self.middle_formlayout.add_widget_to_align(self.muscle_mode_rds)
        self.middle_formlayout.set_input_widgets_width("auto")

    def get_muscle_mode(self):
        return self.muscle_modes_dict[self.muscle_mode_rds.getvalue()]
