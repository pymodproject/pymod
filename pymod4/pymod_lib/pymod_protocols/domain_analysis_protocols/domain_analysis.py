# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Protocols to search and assign domains to a protein sequence loaded in PyMod.
"""

import os

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.domain_analysis_protocols.hmmscan import Hmmscan_protocol
from pymod_lib.pymod_protocols.domain_analysis_protocols.split_domains import Split_into_domains_protocol

from pymod_lib.pymod_os_specific import clean_file_name
from pymod_lib.pymod_seq.seq_star_alignment import join_alignments

from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt


class Domain_Analysis_Protocol(PyMod_protocol):
    """
    This class will launch a 'Hmmscan_protocol' protocol used to identify the domains of a protein
    query. It will then store all the information useful to perform the operations flow of the domain
    analysis such as splitting and fuse.
    """

    protocol_name = "Domain Analysis"

    def __init__(self, pymod, domain_search_mode):
        """
        The 'domain_search_mode' can be:
            - "local" (a local hmmscan executable will be used to search on local databases)
            - "remote" (a remote hmmscan search will be performed on the EBI servers).
        """
        if not domain_search_mode in ("local", "remote"):
            raise KeyError("Unknown 'domain_search_mode': %s" % domain_search_mode)
        PyMod_protocol.__init__(self, pymod,
                                output_directory=pymod.domain_analysis_dirpath)
        self.domain_search_mode = domain_search_mode


    def launch_from_gui(self):

        # Check for only one selected sequence.
        selected_elements = self.pymod.get_selected_sequences()
        if len(selected_elements) != 1:
            title = "Selection Error"
            message = "Please select one and only one sequence to perform a domain search."
            self.pymod.main_window.show_error_message(title, message)
            return None

        self.pymod_element = selected_elements[0]

        # An index to distinguish each domain analysis performed in PyMod.
        self._index_of_domain_protocol = self.pymod.active_domain_analysis_count

        # Sets the name which will be used to label the query sequence in the domain analysis.
        self.query_element_name = "hmmscan_search_%s_%s" % (self.pymod_element.unique_index, self._index_of_domain_protocol)

        # Saves a FASTA file with the query sequence.
        self._pymod_element_seq_filepath = os.path.join(self.output_directory, self.query_element_name + '.fasta')
        self.pymod.build_sequence_file([self.pymod_element], self._pymod_element_seq_filepath, file_format="fasta", remove_indels=True,
                                       use_structural_information=False, add_extension=False,
                                       unique_indices_headers=False)

        # This will be storing instances of the various protocol classes for domain analysis.
        self.search_protocol = None

        # Launches an hmmscan search.
        self.run_hmmscan_search()


    ###########################################################################
    # Search domains with hmmscan.                                            #
    ###########################################################################

    def run_hmmscan_search(self):

        # Reinitializes the domains list of the query element before searching again.
        if self.pymod_element.get_domains_features():
            message = ("Would you like to perform a new domain search operation on this protein?"
                       " Its previous domain search results will be lost.")
            confirmation = askyesno_qt("Confirm Domain Search", message, parent=self.pymod.get_qt_parent())
            if not confirmation:
                return None
            self.pymod_element.clear_domains()
            self.pymod.main_window.gridder(update_clusters=True, update_menus=True)

        # Initializes the domain search protocol.
        self.search_protocol = Hmmscan_protocol(self.pymod, output_directory=self.output_directory,
                                                father_protocol=self)
        self.search_protocol.launch_from_gui()


    def evaluate_domain_search(self):
        # Updates PyMod data on domain searches.
        self.pymod_element.domain_analysis_id = self.pymod.active_domain_analysis_count
        self.pymod.active_domain_analysis_count += 1
        # Updates the menus.
        self.pymod.main_window.gridder(update_menus=True)


###################################################################################################
# Fuse split domains.                                                                             #
###################################################################################################

class Fuse_domains_protocol(PyMod_protocol):
    """
    Protocol used to join in a single alignment the fragments derived from using a
    'Split_into_domains_protocol' to split a query protein in its domains.
    """

    protocol_name = "Domain Fuse"

    def __init__(self, pymod, pymod_element, output_directory=None):
        PyMod_protocol.__init__(self, pymod, output_directory)
        self.pymod_element = pymod_element


    def launch_from_gui(self):
        if not self.check_fuse_conditions():
            return None

        self.run_fuse_protocol()


    def check_fuse_conditions(self):

        fuse_error_title = "Fuse Error"
        fuse_error_message = "Can not perform the Fuse operation."

        # Checks if all the derived domains are currently present in PyMod (that is, that no one
        # got deleted or altered).
        all_pymod_sequences = self.pymod.get_all_sequences()
        for domain_element in self.pymod_element.derived_domains_list:
            if not domain_element in all_pymod_sequences:
                message = "The '%s' sequence has been deleted. %s" % (domain_element.my_header, fuse_error_message)
                self.pymod.main_window.show_error_message(fuse_error_title, message)
                return False

        # Checks whether all the derived domains are currently aligned to other sequences.
        for domain_element in self.pymod_element.derived_domains_list:
            if domain_element.splitted_domain == None:
                message = "The sequence of the domain element '%s' has been edited, it will not be possible to map it to the original query sequence. %s" % (domain_element.my_header, fuse_error_message)
                self.pymod.main_window.show_info_message(fuse_error_title, message)
                return False

            if not domain_element.is_child():
                message = "Not all the domains derived from element '%s' are currently aligned to other sequences. %s" % (self.pymod_element.my_header, fuse_error_message)
                self.pymod.main_window.show_info_message(fuse_error_title, message)
                return False

            other_domain_elements = [d_el for d_el in self.pymod_element.derived_domains_list if not d_el is domain_element]
            if len(set(other_domain_elements) & set(domain_element.get_siblings())) != 0:
                message = "Multiple domain sequences derived from the original query element '%s' can not be present in the same alignment. %s" % (self.pymod_element.my_header, fuse_error_message)
                self.pymod.main_window.show_info_message(fuse_error_title, message)
                return False

        return True


    def run_fuse_protocol(self, verbose=False):

        if verbose:
            print('\n\n__________FUSE___________\n')

        #-----------------------------------------------------------------
        # Builds multiple sequence alignments for each cluster involved. -
        #-----------------------------------------------------------------

        query_seq_gapless = self.pymod_element.my_sequence.replace("-", "")

        # This will contain a list of multiple sequence alignments (each represented by a list of
        # aligned sequences).
        msa_l = [] # This will contain sequences.
        msa_elements_list = [] # This will contain the PyMod elements to which the sequences above belong to.


        # Adds the sequences of the domains and their siblings.
        for domain_i, domain_el in enumerate(self.pymod_element.derived_domains_list):

            # For each domain, prepare its multiple sequence alignment.
            domain_msa = [domain_el.my_sequence]
            domain_elements_msa = [domain_el]

            for sibling in domain_el.get_siblings():
                domain_msa.append(sibling.my_sequence)
                domain_elements_msa.append(sibling)


            # In the first iteration, adds as a first MSA the MSA of the original full length query.
            if domain_i == 0:

                # Prepares the first multiple sequence alignment (the one for the full length query).
                query_msa = [self.pymod_element.my_sequence]
                query_elements_msa = [self.pymod_element]
                if self.pymod_element.is_child():
                    for sibling in self.pymod_element.get_siblings():
                        query_msa.append(sibling.my_sequence)
                        query_elements_msa.append(sibling)
                msa_l = [query_msa]
                msa_elements_l = [query_elements_msa]

            # In the successive iterations, adds as a first MSA the MSA produced in the previous
            # iteration.
            else:
                msa_l = [new_msa]

                _msa_elements_l = []
                for msa_i in msa_elements_l:
                    for seq_j in msa_i:
                        _msa_elements_l.append(seq_j)
                msa_elements_l = [_msa_elements_l]

            msa_l.append(domain_msa)
            msa_elements_l.append(domain_elements_msa)

            # Prepares the "reference" alignment (an alignment containing the two anchor sequences
            # of the two alignments which will be joined).
            j_msa = [query_seq_gapless]
            domain_seq_gapless = domain_el.my_sequence.replace("-", "")
            domain_info = domain_el.splitted_domain
            domain_frag_start = max((0, domain_info.parent_start - domain_info.offset[0]))
            domain_frag_end = domain_info.parent_end + domain_info.offset[0]
            domain_seq_gapped = "-"*domain_frag_start + domain_seq_gapless
            domain_seq_gapped += "-"*(max(0, len(query_seq_gapless)-len(domain_seq_gapped)))
            j_msa.append(domain_seq_gapped)

            # Progressively join the alignments.
            new_msa = join_alignments(j_msa, msa_l)
            if verbose:
                print("\n# Temporary MSA %s" % domain_i)
                for seq in new_msa:
                    print(seq)


        #-------------------
        # Get the results. -
        #-------------------

        elements_to_add = []
        for domain_el in self.pymod_element.derived_domains_list:
            elements_to_add.append(domain_el)
            elements_to_add.extend(domain_el.get_siblings())

        # Adds the domains and their siblings to the full length query alignment.
        if self.pymod_element.is_child():
            self.pymod_element.mother.add_children(elements_to_add)

        # Builds a new alignment object.
        else:
            query_original_index = self.pymod.get_pymod_element_index_in_container(self.pymod_element)
            new_blast_cluster = self.pymod.add_new_cluster_to_pymod(
                cluster_type="alignment",
                child_elements=[self.pymod_element] + elements_to_add,
                algorithm="joined",
                update_stars=False)
            # Move the new cluster to the same position of the original query element in PyMod main
            # window.
            self.pymod.change_pymod_element_list_index(new_blast_cluster, query_original_index)


        # Updates the sequence of the elements involved.
        new_msa_elements = []
        for msa_i in msa_elements_l:
            for seq_j in msa_i:
                new_msa_elements.append(seq_j)
        for element, updated_seq in zip(new_msa_elements, new_msa):
            element.set_sequence(updated_seq)


        # Shows the results in the main window.
        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True, update_elements=True)
