# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Launch BLAST using an executable and databases found on the user's system.
The classes in this model build up on the ones for PSI-BLAST.
"""

import os

from pymod_lib.pymod_protocols.similarity_searches_protocols.psiblast import PSI_BLAST_search, PSI_BLAST_options_window_qt
from pymod_lib.pymod_threading import Protocol_exec_dialog


###################################################################################################
# LOCAL BLAST.                                                                                    #
###################################################################################################

class LOC_BLAST_search(PSI_BLAST_search):

    blast_version = "blastp"
    protocol_name = blast_version
    exe_filename = "blastp"
    min_inclusion_eval_default = 0.005

    def get_blast_window_class_qt(self):
        return Local_BLAST_options_window_qt

    def run_psiblast(self):
        """
        Launches a standalone version of BLAST installed locally when using the BLAST
        option in the plugin main menu.
        """
        # Builds a temporary file with the sequence of the query needed by psiblast.
        query_file_name = "query"
        self.pymod.build_sequence_file([self.blast_query_element], query_file_name,
                                       file_format="fasta", remove_indels=True,
                                       new_directory=self.output_directory)

        # Sets some parameters in needed to run PSI-BLAST.
        ncbi_dir = self.pymod.blast_plus["exe_dir_path"].get_value()

        args = {"ncbi_dir": ncbi_dir,
                "db_path": self.db_path,
                "query": os.path.join(self.output_directory, query_file_name+".fasta"),
                "inclusion_ethresh": None,
                "outfmt": 5,
                "out": os.path.join(self.output_directory, self.xml_blast_output_file_name),
                "num_iterations": None,
                "evalue": self.evalue_cutoff,
                "max_target_seqs": self.max_hsp_num,
                "blast_version": "blastp"}

        if not self.pymod.use_protocol_threads:
            self.execute_psiblast(**args)

        else:
            label_text = "Running BLASTp. Please wait for the process to complete..."
            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.execute_psiblast,
                                            args=args,
                                            title="Running BLASTp",
                                            label_text=label_text)
            p_dialog.exec_()

        # If everything went ok, return 'True', so that the results window can be opened.
        return True


    def get_options_from_gui_specific(self):
        self.get_options_from_gui_blast()


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Local_BLAST_options_window_qt(PSI_BLAST_options_window_qt):

    def build_algorithm_advanced_options_widgets(self):
        pass
