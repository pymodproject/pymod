# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for performing hmmsearch (profile-vs-sequence) searches in PyMod. It mainly
builds up on the code used to perform BLAST searches and phmmer searches.
"""

import os

from Bio import SeqIO

from pymod_lib.pymod_protocols.similarity_searches_protocols.phmmer import PHMMER_search, Phmmer_options_window_qt
from pymod_lib.pymod_os_specific import get_exe_file_name
from pymod_lib.pymod_seq.seq_manipulation import remove_gap_only_columns


################################################################################################
# HMMSEARCH.                                                                                   #
################################################################################################

class Hmmsearch_search(PHMMER_search):

    blast_version = "hmmsearch"
    protocol_name = blast_version
    exe_filename = None
    all_exe_filenames = ["hmmbuild", "hmmsearch"]


    def launch_from_gui(self):

        # Check if a correct selection is provided.
        selected_clusters_list = self.pymod.get_selected_clusters()

        if len(selected_clusters_list) == 0:
            title = "Selection Error"
            message = "Please select an entire alignment object to perform a %s search" % (self.blast_version_full_name)
            self.pymod.main_window.show_error_message(title, message)
            return None

        if len(selected_clusters_list) > 1:
            title = "Selection Error"
            message = "Please select only one alignment object to perform a %s search" % (self.blast_version_full_name)
            self.pymod.main_window.show_error_message(title, message)
            return None

        if not self.check_blast_program():
            return None

        # Gets the selected sequence. The main index will be used later to build the cluster.
        self.blast_query_element = selected_clusters_list[0]

        # Opens the options window.
        self.build_blast_window()


    def get_blast_window_class_qt(self):
        return Hmmsearch_options_window_qt


    def execute_hmmer_program(self, query_filepath, out_filepath, db_filepath,
                              exe_filepath, hmmbuild_exe_filepath,
                              report_evalue=10e-6):
        """
        Execute the locally installed hmmbuild and hmmsearch.
        """

        # Launch hmmbuild to buil a profile HMM from the alignment file saved from
        # PyMod.
        hmm_filepath = query_filepath + ".hmm"
        command_ls = [hmmbuild_exe_filepath,
                      hmm_filepath,
                      query_filepath]
        self.pymod.new_execute_subprocess(command_ls)

        # Launche hmmsearch.
        self.hmmscan_ali_filepath = out_filepath + ".sth"
        command_ls = [exe_filepath,
                      "-o", out_filepath,
                      "-A", self.hmmscan_ali_filepath,
                      "--domE", str(report_evalue),
                      hmm_filepath, db_filepath]
        self.pymod.new_execute_subprocess(command_ls)


    def import_results_in_pymod(self):
        """
        Builds a cluster with the hit sequences.
        """

        # The list of elements whose sequences will be updated according to the star alignment.
        hsp_elements = []
        use_hmmer_pdb = self.db_filename.startswith("pdbaa")

        #------------------------------------------------------------
        # Builds a new cluster with the query and all the new hits. -
        #------------------------------------------------------------

        query_original_index = self.pymod.get_pymod_element_index_in_container(self.blast_query_element)

        # Creates PyMod elements for all the imported hits and add them to the cluster.
        for h in self.hsp_imported_from_blast:
            cs = self.pymod.build_pymod_element_from_hsp(self.get_hsp_subj_seq(h["hsp"]),
                                                         self.get_hsp_element_title(h, use_hmmer_pdb))
            self.pymod.add_element_to_pymod(cs)
            hsp_elements.append(cs)

        # Builds the "BLAST search" cluster element.
        new_blast_cluster = self.pymod.add_new_cluster_to_pymod(
            cluster_type="profile-cluster",
            child_elements=hsp_elements,
            algorithm=self.blast_version_full_name,
            update_stars=False)

        # Move the new cluster to the same position of the original query element in PyMod main
        # window.
        self.pymod.change_pymod_element_list_index(new_blast_cluster, query_original_index)


        #---------------------------------
        # Updates the aligned sequences. -
        #---------------------------------

        if os.path.isfile(self.hmmscan_ali_filepath):

            # For each sequence imported in PyMod, searches in the output alignment
            # of hmmsearch for a sequence with the same amino acids in order to
            # import in PyMod the alignment generated by hmmsearch.
            hsp_recs = SeqIO.parse(self.hmmscan_ali_filepath, "stockholm")

            # Stores the gapless sequences of the hits imported in PyMod.
            hsp_seqs = [str(e.my_sequence).replace("-", "") for e in hsp_elements]

            # Iter through the msa sequences (which contain also the hits not
            # imported in PyMod).
            for h in hsp_recs:
                hsp_msa_aliseq = str(h.seq).upper()
                hsp_msa_seq = hsp_msa_aliseq.replace("-", "")
                # Search a matching sequence in the hits imported in PyMod.
                for hsp_idx, (hsp_seq, hsp_element) in enumerate(zip(hsp_seqs, hsp_elements)):
                    if hsp_seq == hsp_msa_seq:
                        # Updates the sequence.
                        hsp_element.my_sequence = hsp_msa_aliseq
                        hsp_seqs.pop(hsp_idx)
                        hsp_elements.pop(hsp_idx)
                        break
                if not hsp_elements:
                    break

            remove_gap_only_columns(new_blast_cluster)


        #------------
        # Finishes. -
        #------------

        # Cleans and update the PyMod main window.
        self.finish_blast_search()

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                       update_menus=True, update_elements=True)


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Hmmsearch_options_window_qt(Phmmer_options_window_qt):
    """
    Window for HMMSEARCH searches.
    """

    pass
