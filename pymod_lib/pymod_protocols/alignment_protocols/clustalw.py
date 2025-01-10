# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
ClustalW.
"""

import os

from Bio.Align.Applications import ClustalwCommandline

# Protocols.
from ._clustal_common import Clustal_regular_alignment, Clustal_profile_alignment

# GUI.
from ._base_alignment._gui import Regular_alignment_window_qt, Profile_alignment_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_radioselect_qt, PyMod_entryfield_qt


class Clustalw_alignment:
    """
    General Clustal W alignments.
    """

    protocol_name = "clustalw"

    def additional_initialization(self):
        self.tool = self.pymod.clustalw

    def get_options_from_gui(self):

        self.params_from_gui = {}
        error_from_gui = False

        try:
            self.params_from_gui["selected_matrix"] = self.alignment_window.get_matrix_value()
        except Exception as e:
            error_from_gui = True
            error_message = "Invalid Matrix."
        try:
            self.params_from_gui["gapopen_value"] = float(self.alignment_window.get_gapopen_value())
        except Exception as e:
            error_from_gui = True
            error_message = str(e) # "Invalid Gap Open Value."
        try:
            self.params_from_gui["gapextension_value"] = float(self.alignment_window.get_gapextension_value())
        except Exception as e:
            error_from_gui = True
            error_message = str(e) # "Invalid Gap Extension Value."

        if error_from_gui:
            self.pymod.main_window.show_error_message("Parameters Error", error_message)
            return False
        else:
            return True


class Clustalw_regular_alignment(Clustalw_alignment, Clustal_regular_alignment):

    def get_alignment_window_class_qt(self):
        return Clustalw_regular_window_qt


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_clustalw(sequences_to_align,
                          output_file_name=output_file_name,
                          matrix=self.params_from_gui["selected_matrix"],
                          gapopen=self.params_from_gui["gapopen_value"],
                          gapext=self.params_from_gui["gapextension_value"])


    def run_clustalw(self, sequences_to_align, output_file_name, matrix="blosum", gapopen=10, gapext=0.2):
        """
        This method allows to interact with the local ClustalW.
        """

        # First build an input FASTA file containing the sequences to be aligned.
        self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)
        # Sets the full paths of input and output files.
        input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
        output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
        # Run an alignment with all the sequences using ClustalW command line, through Biopython.
        cline = ClustalwCommandline(self.pymod.clustalw["exe_file_path"].get_value(),
                infile=input_file_path, outfile=output_file_path, outorder="INPUT",
                matrix=matrix, gapopen=gapopen, gapext=gapext)
        self.pymod.execute_subprocess(str(cline))


class Clustalw_profile_alignment(Clustalw_alignment, Clustal_profile_alignment):

    def get_alignment_window_class(self):
        return Clustalw_profile_window

    def get_alignment_window_class_qt(self):
        return Clustalw_profile_window_qt


    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        clustalw_path = self.tool["exe_file_path"].get_value()
        cline='"'         +clustalw_path+'"'+ \
            ' -PROFILE1="'+profile_file_shortcut+'"'+ \
            ' -PROFILE2="'+sequences_to_add_file_shortcut+'" -SEQUENCES -OUTORDER=INPUT'+ \
            ' -MATRIX='   +self.params_from_gui["selected_matrix"] + \
            ' -GAPOPEN='  +str(self.params_from_gui["gapopen_value"]) + \
            ' -GAPEXT='   +str(self.params_from_gui["gapextension_value"]) + \
            ' -OUTFILE="' +output_file_shortcut+'.aln"'
        return cline


    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalw_path = self.tool["exe_file_path"].get_value()
        cline='"'          +clustalw_path+'"' \
            ' -PROFILE1="' +profile1+'"'+ \
            ' -PROFILE2="' +profile2+'" -OUTORDER=INPUT' \
            ' -MATRIX='    +self.params_from_gui["selected_matrix"]+ \
            ' -GAPOPEN='   +str(self.params_from_gui["gapopen_value"])+ \
            ' -GAPEXT='    +str(self.params_from_gui["gapextension_value"])+ \
            ' -OUTFILE="'  +output_file_shortcut+'.aln"'
        return cline


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Clustalw_base_window_qt:
    """
    Base class for ClustalW protocols.
    """

    def build_algorithm_options_widgets(self):

        # Scoring matrix radioselect.
        self.clustal_matrices = ["Blosum", "Pam", "Gonnet", "Id"]
        self.clustal_matrices_dict = {"Blosum": "blosum", "Pam": "pam", "Gonnet": "gonnet", "Id": "id"}
        self.matrix_rds = PyMod_radioselect_qt(label_text="Scoring Matrix Selection",
                                               buttons=self.clustal_matrices)
        self.matrix_rds.setvalue("Blosum")
        self.middle_formlayout.add_widget_to_align(self.matrix_rds)

        # Gap open entryfield.
        self.gapopen_enf = PyMod_entryfield_qt(label_text="Gap Opening Penalty",
                                               value="10.0",
                                               validate={'validator': 'real',
                                                         'min': 0, 'max': 1000})
        self.middle_formlayout.add_widget_to_align(self.gapopen_enf)

        # Gap extension entryfield.
        self.gapextension_enf = PyMod_entryfield_qt(label_text="Gap Extension Penalty",
                                                    value="0.2",
                                                    validate={'validator': 'real',
                                                              'min': 0, 'max': 1000})
        self.middle_formlayout.add_widget_to_align(self.gapextension_enf)

        self.middle_formlayout.set_input_widgets_width("auto")


    def get_matrix_value(self):
        return self.clustal_matrices_dict[self.matrix_rds.getvalue()]

    def get_gapopen_value(self):
        return self.gapopen_enf.getvalue(validate=True)

    def get_gapextension_value(self):
        return self.gapextension_enf.getvalue(validate=True)


class Clustalw_regular_window_qt(Clustalw_base_window_qt, Regular_alignment_window_qt):
    pass

class Clustalw_profile_window_qt(Clustalw_base_window_qt, Profile_alignment_window_qt):
    pass
