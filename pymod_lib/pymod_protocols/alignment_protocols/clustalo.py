# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Clustal Omega.
"""

import os

from Bio.Align.Applications import ClustalOmegaCommandline

# Protocols.
from ._clustal_common import Clustal_regular_alignment, Clustal_profile_alignment

# GUI.
from ._base_alignment._gui import Regular_alignment_window_qt, Profile_alignment_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_entryfield_qt, PyMod_radioselect_qt


class Clustalomega_alignment:
    """
    General Clustal Omega alignments.
    """

    alignment_program = "clustalo"
    protocol_name = "clustalo"

    def additional_initialization(self):
        self.tool = self.pymod.clustalo


    def get_options_from_gui(self):
        self.params_from_gui = {}
        error_from_gui = False

        self.params_from_gui["use_full_dm"] = self.alignment_window.get_use_full_dm_value()

        try:
            self.params_from_gui["iterations"] = self.alignment_window.get_iterations_value()
        except Exception as e:
            error_from_gui = True
            error_message = str(e) # "Invalid Combined Iterations Value."

        if error_from_gui:
            self.pymod.main_window.show_error_message("Parameters Error", error_message)
            return False
        else:
            return True


class Clustalomega_regular_alignment(Clustalomega_alignment, Clustal_regular_alignment):
    """
    Regular alignments using Clustal Omega.
    """

    def get_alignment_window_class_qt(self):
        return Clustalo_regular_window_qt


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_clustalo(sequences_to_align,
                          output_file_name=output_file_name,
                          extraoption="",
                          iterations=self.params_from_gui["iterations"],
                          use_full_dm=self.params_from_gui["use_full_dm"])


    def run_clustalo(self, sequences_to_align, output_file_name=None, extraoption="", iterations=0, use_full_dm=False):

        self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)

        input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
        output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
        guidetree_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".dnd")

        cline = ClustalOmegaCommandline(
            self.tool["exe_file_path"].get_value(),
            infile= input_file_path,
            outfile= output_file_path,
            guidetree_out=guidetree_file_path,
            force=True, outfmt="clustal")

        cline = str(cline)

        if iterations != 0:
            cline = "%s --iter=%s" % (cline, iterations)
        if use_full_dm:
            cline = "%s --full --full-iter" % (cline)

        # Run MSA with all sequences using CLustalO command line.
        self.pymod.execute_subprocess(cline)


class Clustalomega_profile_alignment(Clustalomega_alignment, Clustal_profile_alignment):
    """
    Profile alignments for Clustal Omega.
    """

    def get_alignment_window_class_qt(self):
        return Clustalo_profile_window_qt


    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        clustalo_path = self.tool["exe_file_path"].get_value()
        cline='"'           +clustalo_path+'"'+ \
            ' --profile1="' +profile_file_shortcut+'"'+ \
            ' --outfile="'  +output_file_shortcut+'.aln"'+ \
            ' --outfmt=clustal --force'+ \
            ' '#  +self.alignment_window.get_extraoption_value()
        if len(self.elements_to_add)>1:
            cline+=' --infile="'  +sequences_to_add_file_shortcut+'"'
        else:
            cline+=' --profile2="'+sequences_to_add_file_shortcut+'"'

        if self.params_from_gui["iterations"] != 0:
            cline = "%s --iter=%s" % (cline, self.params_from_gui["iterations"])
        if self.params_from_gui["use_full_dm"]:
            cline = "%s --full --full-iter" % (cline)

        return cline


    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalo_path = self.tool["exe_file_path"].get_value()
        cline='"'           +clustalo_path+'"' \
            ' --profile1="' +profile1+'"'+ \
            ' --profile2="' +profile2+'"'+ \
            ' --outfile="'  +output_file_shortcut+'.aln"' \
            ' --outfmt=clustal --force' \
            ' '#  +self.alignment_window.get_extraoption_value()

        if self.params_from_gui["iterations"] != 0:
            cline = "%s --iter=%s" % (cline, self.params_from_gui["iterations"])
        if self.params_from_gui["use_full_dm"]:
            cline = "%s --full --full-iter" % (cline)
        return cline


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Clustalo_base_window_qt:
    """
    Base class for ClustalOmega protocols.
    """

    def build_algorithm_options_widgets(self):

        # Use full distance matrix.
        self.use_full_dm_rds = PyMod_radioselect_qt(label_text="Use Full Distance Matrix",
                                                    buttons=('Yes', 'No'))
        self.use_full_dm_rds.setvalue("No")
        self.middle_formlayout.add_widget_to_align(self.use_full_dm_rds)

        # Number of (combined guide-tree/HMM) iterations.
        self.clustalo_iterations_enf = PyMod_entryfield_qt(label_text="Combined Iterations",
                                                           value='0',
                                                           validate={'validator': 'integer',
                                                                     'min': 0, 'max': 5})
        self.middle_formlayout.add_widget_to_align(self.clustalo_iterations_enf)

        self.middle_formlayout.set_input_widgets_width("auto")

    def get_iterations_value(self):
        return self.clustalo_iterations_enf.getvalue(validate=True)

    def get_use_full_dm_value(self):
        return self.use_full_dm_rds.getvalue() == "Yes"


class Clustalo_regular_window_qt(Clustalo_base_window_qt, Regular_alignment_window_qt):
    pass

class Clustalo_profile_window_qt(Clustalo_base_window_qt, Profile_alignment_window_qt):
    pass
