# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Profile alignments.
"""

from pymod_lib import pymod_vars

from pymod_lib.pymod_protocols.alignment_protocols._base_alignment import Alignment_protocol
from pymod_lib.pymod_exceptions import catch_protocol_exception


class Profile_alignment(Alignment_protocol):

    alignment_strategy = "profile-alignment"

    #################################################################
    # Start the alignment process.                                  #
    #################################################################

    def check_alignment_selection(self):
        """
        Checks if the selected elements can be used to perform a profile alignment.
        """
        # This will be set to True if there is an adequate selection in order to align two profiles.
        self.can_perform_ptp_alignment = False

        # Checks if there is at least one cluster which is entirely selected.
        number_of_selected_clusters = len(self.selected_clusters_list)
        number_of_involved_clusters = len(self.involved_clusters_list)
        number_of_root_sequences = len(self.selected_root_sequences_list)

        # No clusters are involved.
        if number_of_selected_clusters == 0:
            return False
        # If there is only one selected cluster and if there is at least one selected sequence this
        # cluster, then a sequence to profile alignment can be performed.
        if number_of_involved_clusters == 1 and number_of_selected_clusters == 1 and number_of_root_sequences > 0:
            return True
        # Two involved clusters.
        elif number_of_involved_clusters == 2:
            # If there aren't any other selected sequences a profile to profile alignment can be
            # performed.
            if number_of_selected_clusters == 2 and number_of_root_sequences == 0:
                self.can_perform_ptp_alignment = True
            return True
        # Only sequence to profile alignments can be performed.
        elif number_of_involved_clusters >= 3:
            return True
        else:
            return False


    def selection_not_valid(self):
        title = "Selection Error"
        message = ("Please select at least one entire cluster and some other sequences"
                   " in order to perform a profile alignment.")
        self.pymod.main_window.show_error_message(title, message)


    def check_sequences_level(self):
        self.clusters_are_involved = True
        return True


    #################################################################
    # Perform the alignment.                                        #
    #################################################################

    def define_alignment_mode(self):
        """
        Gets several parameters from the GUI in order to define the alignment mode.
        """
        # It can be either "sequence-to-profile" or "profile-to-profile".
        self.alignment_mode = self.alignment_window.get_alignment_mode()
        # Takes the index of the target cluster.
        self.target_cluster_index = None
        # Takes the index of the target cluster for the "keep-previous-alignment" mode.
        if self.alignment_mode == "sequence-to-profile":
            # If there is only one cluster involved its index its going to be 0.
            if len(self.selected_clusters_list) == 1:
                self.target_cluster_index = 0 # Cluster index.
            # Get the index of the cluster from the combobox.
            elif len(self.selected_clusters_list) > 1:
                if hasattr(self.alignment_window.target_profile_frame, "get_selected_cluster_index"):
                    self.target_cluster_index = self.alignment_window.target_profile_frame.get_selected_cluster_index()
                else:
                    self.target_cluster_index = self.alignment_window.get_selected_cluster_index()


    @catch_protocol_exception
    def perform_alignment_protocol(self):
        if self.alignment_mode == "sequence-to-profile":
            self.perform_sequence_to_profile_alignment()

        elif self.alignment_mode == "profile-to-profile":
            self.perform_profile_to_profile_alignment()


    ######################################################
    # Methods to perform sequence-to-profile alignments. #
    ######################################################

    def perform_sequence_to_profile_alignment(self):
        self.run_sequence_to_profile_alignment_program()


    #####################################################
    # Methods to perform profile-to-profile alignments. #
    #####################################################

    def perform_profile_to_profile_alignment(self):
        self.run_profile_to_profile_alignment_program()


    #################################################################
    # Import the updated sequences in PyMod.                        #
    #################################################################

    def create_alignment_element(self):

        #---------------------------------------------------------
        # Expand an already existing cluster with new sequences. -
        #---------------------------------------------------------

        if self.alignment_mode == "sequence-to-profile":
            # Gets the target cluster element.
            self.alignment_element = self.involved_clusters_list[self.target_cluster_index]
            # Appends new sequences to the target cluster.
            for element in self.elements_to_add:
                self.alignment_element.add_child(element)
            # Updates the alignment element with new information about the new alignment.
            self.alignment_element.algorithm = "merged"
            # alignment_description = "merged with %s" % (pymod_vars.algs_full_names_dict[self.protocol_name])
            alignment_description = "merged"
            self.alignment_element.my_header = self.pymod.set_alignment_element_name(alignment_description, self.alignment_element.cluster_id)

        #--------------------------------------
        # Join two or more existing clusters. -
        #--------------------------------------

        elif self.alignment_mode == "profile-to-profile":
            # Find the right mother index in order to build the new cluster where one of the
            # original ones was placed.
            lowest_index = min([self.pymod.get_pymod_element_index_in_root(e) for e in self.elements_to_align])
            # Orders them.
            self.elements_to_align = sorted(self.elements_to_align, key=lambda el: (self.pymod.get_pymod_element_index_in_root(el),
                                                                                    self.pymod.get_pymod_element_index_in_container(el)))

            alignment_description = "joined by using " + pymod_vars.algs_full_names_dict[self.protocol_name]
            # ali_name = "Joined " + self.pymod.set_alignment_element_name(alignment_description, self.pymod.alignment_count)

            # Builds the new "PyMod_element" object for the new alignment.
            new_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                # cluster_name=ali_name,
                                                child_elements=self.elements_to_align,
                                                algorithm=self.protocol_name,
                                                update_stars=True)
            # Moves the new element from the bottom of the list to its new position.
            self.pymod.change_pymod_element_list_index(new_cluster, lowest_index)
