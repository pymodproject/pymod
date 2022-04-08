# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module with a base class for similarity searches protocols. It is used to build classes for every
similarity search protocol in PyMod, that is, (PSI-)BLAST protocols and phmmer/jackhmmer protocols.
"""

import os

from pymod_lib.pymod_os_specific import get_exe_file_name, clean_file_name
from pymod_lib.pymod_vars import algs_full_names_dict
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_seq.seq_star_alignment import save_cstar_alignment

from pymol.Qt import QtWidgets, QtCore
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          small_font_style, highlight_color)


###############################################################################################
# BLAST ALGORITHMS.                                                                           #
###############################################################################################

hmmer_protocols_names = ("phmmer", "jackhmmer", "hmmsearch")

results_header_options = {'background': 'black', 'fg': 'red', 'height': 1, 'pady': 10, 'font': "comic 12"}
results_row_options = {'background': 'black', 'fg': 'white', 'height': 1, 'highlightbackground': 'black'}


class Generic_BLAST_search(PyMod_protocol):
    """
    Base class, which is used to build up classes for specific similarity searches programs.
    """

    #################################################################
    # Initialize the launching of BLAST programs.                   #
    #################################################################

    # These attributes will be defined in each subclass for implementing a different similarity
    # search protocol.
    blast_version = None
    protocol_name = None
    exe_filename = None
    # Common attributes.
    xml_blast_output_file_name = "blast_out.xml"
    default_max_number_of_hits = 100
    e_value_threshold_default = 10.0

    def additional_initialization(self):
        self.blast_version_full_name = algs_full_names_dict[self.blast_version]

    def launch_from_gui(self):

        # Check if a correct selection is provided.
        selected_sequences = self.pymod.get_selected_sequences()
        if len(selected_sequences) == 0:
            title = "Selection Error"
            message = "Please select a sequence to perform a %s search" % (self.blast_version_full_name)
            self.pymod.main_window.show_error_message(title, message)
            return None
        if len(selected_sequences) > 1:
            title = "Selection Error"
            message = "Please select only one sequence to perform a %s search" % (self.blast_version_full_name)
            self.pymod.main_window.show_error_message(title, message)
            return None

        if not self.check_blast_program():
            return None

        # Gets the selected sequence. The main index will be used later to build the cluster.
        self.blast_query_element = selected_sequences[0]

        # Opens the options window.
        self.build_blast_window()


    def check_blast_program(self):
        """
        Default method to check for an executable and locally installed databases.
        Called as one of the first things when performing a similarity search. Actually used in the
        PSI-BLAST, local BLAST and phmmer protocols. May be overridden in child classes.
        """

        # Checks for the tool's directory.
        if not self.tool.tool_dir_exists():
            self.tool.tool_dir_not_found()
            return False

        # Checks for the tool's executable file.
        exe_filepaths = self.get_similariry_search_exe_path()
        for exe_filepath in exe_filepaths:
            if not os.path.isfile(exe_filepath):
                exe_dirpath = self.tool["exe_dir_path"].get_value()
                exe_filename = os.path.basename(exe_filepath)
                title = "%s is Missing" % (self.blast_version_full_name)
                message = ("The %s executable file ('%s') does not exists in the executables"
                           " directory specified in the PyMod Options Window ('%s')."
                           " If you want to use this function, please specify in the Options"
                           " Window an executables directory where a '%s' file is"
                           " found." % (self.blast_version_full_name, exe_filename, exe_dirpath, exe_filename))
                self.pymod.main_window.show_error_message(title, message)
                return False

        # Checks if a database directory is defined in the PyMod options.
        db_dirpath = self.tool["database_dir_path"].get_value()
        if db_dirpath.replace(" ", "") == "":
            title = "No Database Directory"
            message = ("No databases directory is defined for %s. Please define"
                       " a database directory in the PyMod options window in order"
                       " to perform a %s search." % (self.blast_version_full_name, self.blast_version_full_name))
            self.pymod.main_window.show_error_message(title, message)
            return False

        # Check if the database directory actually exits.
        if not os.path.isdir(db_dirpath):
            title = "Database Directory not Found"
            message = ("The database directory '%s' does not exists. Please define"
                       " an existing database directory in the PyMod options window"
                       " in order to perform a %s search." % (db_dirpath, self.blast_version_full_name))
            self.pymod.main_window.show_error_message(title, message)
            return False

        # Checks if any databases are present.
        self.databases_directories_list = self.build_blast_db_list()
        if len(self.databases_directories_list) == 0:
            title = "No Databases"
            message = ("No %s databases were found at path '%s'. Please install"
                       " some databases in order to perform a %s"
                       " search." % (self.blast_version_full_name, db_dirpath, self.blast_version_full_name))
            self.pymod.main_window.show_error_message(title, message)
            return False

        return True


    def get_similariry_search_exe_path(self, return_dict=False):
        """
        Get the absolute filepaths of the executables files needed to perform the
        similarity search.
        """
        if self.protocol_name != "hmmsearch": # Returns only one executable filepath.
            exe_filepath_list = [os.path.join(self.tool["exe_dir_path"].get_value(),
                                              get_exe_file_name(self.exe_filename))]
            exe_filepath_dict = {self.exe_filename: exe_filepath_list[0]}

        else: # Returns two executables filepath.
            exe_filepath_list = []
            exe_filepath_dict = {}
            for exe_filename in self.all_exe_filenames:
                exe_filepath = os.path.join(self.tool["exe_dir_path"].get_value(),
                                            get_exe_file_name(exe_filename))
                exe_filepath_list.append(exe_filepath)
                exe_filepath_dict[exe_filename] = exe_filepath

        if return_dict:
            return exe_filepath_dict
        else:
            return exe_filepath_list


    #################################################################
    # BLAST programs options window.                                #
    #################################################################

    def build_blast_window(self):
        """
        Builds a window containing the widget necessary to define the options for BLAST and
        PSI-BLAST searches.
        """

        blast_window_class_qt = self.get_blast_window_class_qt()
        self.blast_options_window = blast_window_class_qt(self.pymod.main_window, protocol=self,
            title="%s Search Options" % (self.blast_version_full_name),
            upper_frame_title="Here you can modify options for %s" % (self.blast_version_full_name),
            submit_command=self.blast_window_state)
        self.blast_options_window.show()


    #################################################################
    # Actually Run BLAST programs.                                  #
    #################################################################

    def blast_window_state(self):
        """
        This function is called when the 'SUBMIT' button in the BLAST options window is pressed.
        """

        # Do not proceed if users have not provided a correct set of input parameters through
        # the GUI.
        if not self.check_blast_input_parameters():
            return False

        self.get_options_from_gui()
        self.blast_options_window.destroy()

        # Performs a similarity search with the appropriate program.
        blast_status = self.run_blast_program()

        # Displays the window with results.
        if blast_status and self.check_blast_results():
            self.show_blast_output_window()
        # Removes temp files at the end of the whole process, both if hits were found or not.
        else:
            self.finish_blast_search()


    def check_blast_results(self):
        """
        Checks if at least one hsp was identified in the search and stores the results in
        'self.search_record'.
        """

        # Parse the output file built by the search program.
        self.search_record = self.get_search_record()

        # Obtains a list of HSPs objects from Biopython.
        self.hsp_list = self.build_hsp_list()

        # Filters the BLAST record according to the advanced options.
        if self.blast_options_window.showing_advanced_widgets:
            self.filter_results_with_advanced_options()

        # Only take the top hits defined by the "max hits" parameter in the GUI.
        self.hsp_list = self.hsp_list[:self.max_hsp_num]

        # Exit the whole process if no hits were found.
        if len(self.hsp_list) == 0:
            self.pymod.main_window.show_warning_message("%s Message" % self.blast_version_full_name,
                "No hits were found by %s." % self.blast_version_full_name)
            return False
        # Returns 'True' if some hits were found.
        else:
            return True


    def check_min_max_seqid_input(self):
        if self.blast_options_window.showing_advanced_widgets:
            min_id = float(self.blast_options_window.min_id_enf.getvalue())
            max_id = float(self.blast_options_window.max_id_enf.getvalue())
            if min_id >= max_id:
                title = "Input Error"
                message = "The 'Min Id%%' value (%s) must be lower than the 'Max Id%%' one (%s)." % (min_id, max_id)
                self.pymod.main_window.show_error_message(title, message)
                return False
        return True

    def get_options_from_gui(self):
        """
        Get the options and parameters from the options window.
        """
        self.max_hsp_num = int(self.blast_options_window.max_hits_enf.getvalue())
        if self.blast_query_element.is_child():
            self.new_sequences_import_mode = self.blast_options_window.get_import_mode()
        else:
            self.new_sequences_import_mode = "build-new"
        self.get_options_from_gui_specific()

    def get_options_from_gui_specific(self):
        """
        To be implemented in each child classes.
        """
        pass

    def get_options_from_gui_blast(self):
        """
        Get from the GUI those options which are common to all BLAST protocols.
        """
        self.get_db_from_gui_blast()
        self.evalue_cutoff = self.blast_options_window.e_value_threshold_enf.getvalue()
        if self.blast_options_window.showing_advanced_widgets:
            self.min_id = self.blast_options_window.min_id_enf.getvalue()
            self.max_id = self.blast_options_window.max_id_enf.getvalue()
            self.min_coverage = self.blast_options_window.min_coverage_enf.getvalue()

    def get_db_from_gui_blast(self):
        """
        Gets the database value for BLAST protocols from the GUI.
        """
        self.db_path = self.get_database_from_gui()


    def get_search_record(self):
        # An attribute where is going to be stored a Biopython "Blast" record class object.
        result_handle = open(os.path.join(self.output_directory, self.xml_blast_output_file_name), "r")
        # The 'get_blast_record' is overridden in child classes (each class represents a specific
        # BLAST vesion).
        search_record = self.get_blast_record(result_handle)
        result_handle.close()
        return search_record

    def filter_results_with_advanced_options(self):
        min_id = float(self.min_id)
        max_id = float(self.max_id)
        min_coverage = float(self.min_coverage)
        for hsp in self.hsp_list[:]:
            # if ((int(self.max_id) <= hsp.pymod_info["seqid"]*100 < int(self.min_id)) or
            #     hsp.pymod_info["query_span"]*100 < int(self.min_coverage)):
            #     self.hsp_list.remove(hsp)
            remove_hsp = False
            seqid = hsp.pymod_info["seqid"]*100
            if seqid < min_id:
                remove_hsp = True
            if seqid > max_id:
                remove_hsp = True
            if hsp.pymod_info["query_span"]*100 < min_coverage:
                remove_hsp = True
            if remove_hsp:
                self.hsp_list.remove(hsp)


    def get_hsp_info(self, hsp, full_query_sequence = None):
        """
        Gets a Biopython HSP object and computes additional information on it and returns it as a
        dictionary.
        """

        if self.protocol_name in hmmer_protocols_names:
            offset = 1
        else:
            offset = 0

        # Gets the query span.
        if self.protocol_name != "hmmsearch":
            qt = len(str(self.blast_query_element.my_sequence).replace("-",""))
        else:
            qt = self.profile_length
        qs = hsp.query_start + offset
        qe = len(str(self.get_hsp_query_seq(hsp)).replace("-", "")) + qs
        query_span = float(qe - qs)/float(qt)

        # Gets the id% of the hsp.
        matches, identities = self.get_hsp_matches(hsp)
        idp = float(identities)/matches

        # Gets the subject span.
        hs = self.get_subj_start(hsp) + offset
        he = len(self.get_hsp_subj_seq(hsp).replace("-","")) + hs

        # Returns the information.
        additional_infor_dict = {"identities": identities, "seqid": idp,
                                 "query_span": query_span, "matches": matches,
                                 "query_start": qs, "query_end": qe,
                                 "sbjct_start": hs, "sbjct_end": he}

        return additional_infor_dict


    def get_hsp_matches(self, hsp):
        q = self.get_hsp_query_seq(hsp)
        h = self.get_hsp_subj_seq(hsp)
        matches_count = 0
        identities_count = 0
        for qi, hi in zip(q, h):
            if qi != "-" and hi != "-":
                matches_count += 1
                if qi == hi:
                    identities_count += 1
        return matches_count, identities_count


    def get_blast_output_basename(self):
        basename = (clean_file_name(self.blast_query_element.compact_header) + "_" +
                    self.blast_version_full_name + "_" +
                    "search_%s" % (self.pymod.blast_cluster_counter) )
        return basename


    def finish_blast_search(self):
        output_filename = self.get_blast_output_basename() + ".xml"
        try:
            os.rename(os.path.join(self.output_directory, self.xml_blast_output_file_name),
                      os.path.join(self.output_directory, output_filename))
            for f in os.listdir(self.output_directory):
                if not os.path.splitext(f)[-1] == ".xml":
                    os.remove(os.path.join(self.output_directory, f))
        except IOError:
            pass

    def quit_protocol(self):
        self.finish_blast_search()
        PyMod_protocol.quit_protocol(self)


    #################################################################
    # Show BLAST programs output.                                   #
    #################################################################

    def show_blast_output_window(self):
        """
        Displays the window with results from BLAST in a new window.
        """
        self.blast_output_window = Similarity_searches_results_window_qt(parent=self.pymod.main_window,
                                                           protocol=self)
        self.blast_output_window.show()


    def build_hsp_list(self):
        hsp_list = []
        for alignment in self.search_record.alignments:
            for hsp in alignment.hsps:
                hsp._title = alignment.title
                # Gets additional information on the hsp (such as its the id % and query span).
                hsp.pymod_info = self.get_hsp_info(hsp)
                hsp_list.append(hsp)
        # Sort by evalue.
        hsp_list = list(sorted(hsp_list, key=lambda h: h.expect))
        return hsp_list

    def get_subject_name(self, hsp):
        full_title = hsp._title
        if len(full_title) > 100:
            return full_title[0:100] + "..."
        else:
            return full_title

    def get_hsp_evalue(self, hsp):
        return hsp.expect

    def get_hsp_identities(self, hsp):
        return hsp.identities

    def get_subj_start(self, hsp):
        return hsp.sbjct_start


    #################################################################
    # Import BLAST results in PyMod.                                #
    #################################################################

    def blast_results_state(self):
        """
        Called when the 'SUBMIT' button is pressed on some BLAST results window.
        """
        # For each hsp takes the state of its check-box.
        self.my_blast_map = [int(chk.isChecked()) for chk in self.blast_output_window.sbjct_checkbuttons_list]

        # If the user selected at least one HSP.
        if 1 in self.my_blast_map:
            self.build_hits_to_import_list()
            # This will actually import the sequences inside Pymod.
            self.import_results_in_pymod()

        self.blast_output_window.destroy()


    def build_hits_to_import_list(self):
        """
        Builds a list containing those hits that were selected by the user in the BLAST results
        window.
        """

        # This will be used to build PyMod elements out of the subjects of the HSP identified by
        # BLAST.
        self.hsp_imported_from_blast = []
        self.total_hsp_counter = 0 # Counts the total number of hsp.
        self.total_fetched_hsp_counter = 0 # Counts the total number of fetched hsp.

        for hsp in self.hsp_list:
            fetch_hsp = False
            if self.my_blast_map[self.total_hsp_counter] == 1:
                # Appends the hits (subjects).
                self.hsp_imported_from_blast.append({"hsp": hsp, "title": self.get_subject_name(hsp)})
                self.total_fetched_hsp_counter += 1
            self.total_hsp_counter+=1


    def import_results_in_pymod(self):
        """
        Builds a cluster with the query sequence as a mother and retrieved hits as children.
        """

        # The list of elements whose sequences will be updated according to the star alignment.
        elements_to_update = [self.blast_query_element]
        hsp_elements = []
        use_hmmer_pdb = self.protocol_name in hmmer_protocols_names and self.db_filename.startswith("pdbaa")

        #------------------------------------------------------------
        # Builds a new cluster with the query and all the new hits. -
        #------------------------------------------------------------

        if self.new_sequences_import_mode == "build-new":

            # Gets the original index of the query element in its container.
            query_container = self.blast_query_element.mother
            query_original_index = self.pymod.get_pymod_element_index_in_container(self.blast_query_element)

            # Creates PyMod elements for all the imported hits and add them to the cluster.
            for h in self.hsp_imported_from_blast:
                # Gives them the query mother_index, to make them its children.
                cs = self.pymod.build_pymod_element_from_hsp(self.get_hsp_subj_seq(h["hsp"]),
                                                             self.get_hsp_element_title(h, use_hmmer_pdb))
                self.pymod.add_element_to_pymod(cs)
                elements_to_update.append(cs)
                hsp_elements.append(cs)

            # Builds the "BLAST search" cluster element.
            new_blast_cluster = self.pymod.add_new_cluster_to_pymod(
                cluster_type="blast-cluster",
                query=self.blast_query_element,
                child_elements=elements_to_update,
                algorithm=self.blast_version_full_name,
                update_stars=False)

            # Move the new cluster to the same position of the original query element in PyMod main
            # window.
            keep_in_mother_cluster = False
            if keep_in_mother_cluster:
                if not query_container.is_root():
                    query_container.add_child(new_blast_cluster)
                self.pymod.change_pymod_element_list_index(new_blast_cluster, query_original_index)
            else:
                if query_container.is_root():
                    self.pymod.change_pymod_element_list_index(new_blast_cluster, query_original_index)
                else:
                    self.pymod.change_pymod_element_list_index(new_blast_cluster, self.pymod.get_pymod_element_index_in_container(query_container)+1)

            sibling_elements = []


        #----------------------------------------------------------------------------
        # Expand the original cluster of the query by appending to it the new hits. -
        #----------------------------------------------------------------------------

        elif self.new_sequences_import_mode == "expand":

            # The list of elements whose sequences will be updated according to the star alignment.
            elements_to_update = []
            # Begins with the query element.
            elements_to_update.append(self.blast_query_element)
            sibling_elements = self.blast_query_element.get_siblings(sequences_only=True)
            elements_to_update.extend(sibling_elements)

            new_blast_cluster = self.blast_query_element.mother
            # Creates PyMod elements for all the imported hits and add them to the cluster.
            for h in self.hsp_imported_from_blast:
                # Gives them the query mother_index, to make them its children.
                cs = self.pymod.build_pymod_element_from_hsp(self.get_hsp_subj_seq(h["hsp"]),
                                                             self.get_hsp_element_title(h, use_hmmer_pdb))
                self.pymod.add_element_to_pymod(cs)
                elements_to_update.append(cs)
                hsp_elements.append(cs)
                new_blast_cluster.add_child(cs)

            # Sets the query elements as the lead of its cluster.
            self.blast_query_element.set_as_query()


        #------------------------------------------------------------------------
        # Builds a center star alignment in which the query is the center star. -
        #------------------------------------------------------------------------

        aligned_pairs = []
        seqs = []
        all_ids = []

        query_seq = self.blast_query_element.my_sequence.replace("-", "")
        seqs_to_update_dict = {self.blast_query_element.get_unique_index_header(): self.blast_query_element}

        for hsp_i, hsp in enumerate(self.hsp_imported_from_blast):

            hsp_query_seq = self.get_hsp_query_seq(hsp["hsp"])
            hsp_subj_seq = self.get_hsp_subj_seq(hsp["hsp"])

            if self.protocol_name in hmmer_protocols_names:
                offset = 1
            else:
                offset = 0

            hsp_q_res_len = len(hsp_query_seq.replace("-", ""))
            hsp_q_start_index = hsp["hsp"].query_start - 1 + offset

            # Substitute the original query sequence with the portion found in the hsp alignment.
            updated_query_seq = query_seq[:hsp_q_start_index] + hsp_query_seq + query_seq[hsp_q_start_index+hsp_q_res_len:]

            # Add gap extensions to the subject N and C-term.
            updated_subj_seq = len(query_seq[:hsp_q_start_index])*"-" + hsp_subj_seq + len(query_seq[hsp_q_start_index+hsp_q_res_len:])*"-"

            if query_seq.replace("-", "") != updated_query_seq.replace("-", ""):
                raise Exception("The query sequence returned by '%s' is not the same of the original one." % self.blast_version)

            # Add the (query, subject) tuple to 'aligned_pairs'.
            aligned_pairs.append((updated_query_seq, updated_subj_seq))
            # Add the gapless sequences and unique ids to the 'seqs' and 'all_ids'.
            if hsp_i == 0:
                seqs.append(query_seq)
                all_ids.append(self.blast_query_element.get_unique_index_header())
            seqs.append(hsp_subj_seq.replace("-", ""))
            all_ids.append(hsp_elements[hsp_i].get_unique_index_header())

            seqs_to_update_dict[hsp_elements[hsp_i].get_unique_index_header()] = hsp_elements[hsp_i]

        # Add sibling elements if necessary.
        if self.new_sequences_import_mode == "expand":
            for sibling in sibling_elements:
                aligned_pairs.append((self.blast_query_element.my_sequence, sibling.my_sequence))
                seqs.append(sibling.my_sequence.replace("-", ""))
                all_ids.append(sibling.get_unique_index_header())
                seqs_to_update_dict[sibling.get_unique_index_header()] = sibling

        # Adds empty (None) elements to the 'aligned_pairs' lists. In this way, the subjects will
        # be pairwise aligned in the 'save_cstar_alignment' method.
        aligned_pairs.extend([None]*int((len(seqs)-1)*(len(seqs)-2)/2))

        # Joins the alignments and get the updated Biopython records.
        recs = save_cstar_alignment(seqs=seqs, all_ids=all_ids, pairwise_alis=aligned_pairs)


        #------------
        # Finishes. -
        #------------

        # Updates the PyMod sequences.
        for r in recs:
            seqs_to_update_dict[r.id].my_sequence = str(r.seq)

        # Cleans and update the PyMod main window.
        self.finish_blast_search()

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True, update_elements=True)


    def get_hsp_subj_seq(self, hsp_obj):
        return str(hsp_obj.sbjct)


    def get_hsp_query_seq(self, hsp_obj):
        return str(hsp_obj.query)


    def get_hsp_element_title(self, hsp, use_hmmer_pdb=False):
        hsp_header = hsp["title"]
        if use_hmmer_pdb:
            # Changes the name of the element header in a format similar to those of BLAST hits from
            # the NCBI databases.
            try:
                pdb_fields = hsp["title"].split("_")
                pdb_id = pdb_fields[0]
                pdb_chain = pdb_fields[1]
                hsp_header = "pdb|%s|%s" % (pdb_id, pdb_chain)
            except:
                pass
        return hsp_header


###############################################################################################
# GUI for BLAST options.                                                                      #
###############################################################################################

class BLAST_base_options_window_qt(PyMod_protocol_window_qt):
    """
    Base class for windows used in the similarity searches protocols.
    """

    input_widget_width = 10
    # geometry_string = "450x550"

    def build_protocol_middle_frame(self):

        #------------------
        # Simple options. -
        #------------------

        # Let users decide how to import new sequences when the query is a child element (that
        # is, it is already present in a cluster).
        if self.protocol.blast_query_element.is_child():

            self.import_mode_vbox = QtWidgets.QVBoxLayout()
            self.import_mode_label = QtWidgets.QLabel("Hit Import Mode")
            self.import_mode_vbox.addWidget(self.import_mode_label)

            self.import_mode_button_group = QtWidgets.QButtonGroup()
            self.middle_formlayout.addRow(self.import_mode_vbox)

            # Build new alignment.
            new_alignment_radiobutton = QtWidgets.QRadioButton("Build a new alignment with the query and the new hit sequences")
            new_alignment_radiobutton._value = "build-new"
            new_alignment_radiobutton.setChecked(True)
            self.import_mode_vbox.addWidget(new_alignment_radiobutton)
            self.import_mode_button_group.addButton(new_alignment_radiobutton)

            # Expand alignment.
            expand_alignment_radiobutton = QtWidgets.QRadioButton("Expand the existing alignment by appending the new hit sequences")
            # "Expand the already existing cluster by appending to it the new hit sequences"
            expand_alignment_radiobutton._value = "expand"
            self.import_mode_vbox.addWidget(expand_alignment_radiobutton)
            self.import_mode_button_group.addButton(expand_alignment_radiobutton)


        # Each algorithm will have its own standard widgets.
        self.build_algorithm_standard_options_widgets()


        # E-value selection.
        if self.protocol.protocol_name in hmmer_protocols_names:
            e_value_threshold_enf_text = "c-Evalue Threshold"
        else:
            e_value_threshold_enf_text = "E-value Threshold"

        self.e_value_threshold_enf = PyMod_entryfield_qt(label_text=e_value_threshold_enf_text,
                                                         value=str(self.protocol.e_value_threshold_default),
                                                         validate={'validator': 'real',
                                                                   'min': 0.0, 'max': 1000.0})
        self.middle_formlayout.add_widget_to_align(self.e_value_threshold_enf, validate=True)


        # Max hit number selection.
        self.max_hits_enf = PyMod_entryfield_qt(label_text="Max Number of Hits",
                                                value=str(self.protocol.default_max_number_of_hits),
                                                validate={'validator': 'integer',
                                                          'min': 1, 'max': 5000})
        self.middle_formlayout.add_widget_to_align(self.max_hits_enf, validate=True)


        # -------------------
        # Advanced options. -
        # -------------------

        self.show_advanced_button()

        # Minimum id% on with query.
        self.min_id_enf = PyMod_entryfield_qt(label_text="Min Id% Threshold",
                                              value="0",
                                              validate={'validator': 'integer',
                                                        'min': 0, 'max': 100})
        self.middle_formlayout.add_widget_to_align(self.min_id_enf, advanced_option=True, validate=True)


        # Maximum id% on with query.
        self.max_id_enf = PyMod_entryfield_qt(label_text="Max Id% Threshold",
                                              value="100",
                                              validate={'validator': 'integer',
                                                        'min': 0, 'max': 100})
        self.middle_formlayout.add_widget_to_align(self.max_id_enf, advanced_option=True, validate=True)


        # Minimum coverage on the query.
        self.min_coverage_enf = PyMod_entryfield_qt(label_text="Min Coverage% Threshold",
                                                    value="0",
                                                    validate={'validator': 'integer',
                                                              'min': 0, 'max': 100})
        self.middle_formlayout.add_widget_to_align(self.min_coverage_enf, advanced_option=True, validate=True)


        # Advanced options for a specific algorithm.
        self.build_algorithm_advanced_options_widgets()


        self.middle_formlayout.set_input_widgets_width("auto")


    def get_import_mode(self):
        for radiobutton in self.import_mode_button_group.buttons():
            if radiobutton.isChecked():
                return radiobutton._value
        raise ValueError("No import mode was selected.")


    # Override in child classes.
    def build_algorithm_standard_options_widgets(self):
        pass

    def build_algorithm_advanced_options_widgets(self):
        pass


class Similarity_searches_results_window_qt(QtWidgets.QMainWindow):
    """
    Window for showing similarity searches results.
    """

    is_pymod_window = True

    def __init__(self, parent, protocol):

        super(Similarity_searches_results_window_qt, self).__init__(parent)
        self.protocol = protocol

        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle(self._get_window_title())

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        ################
        # Upper frame. #
        ################

        title_text = self._get_upper_frame_title()

        self.upper_frame_title = QtWidgets.QLabel(title_text)
        self.main_vbox.addWidget(self.upper_frame_title)


        #################
        # Middle frame. #
        #################

        # Scroll area which contains the widgets, set as the centralWidget.
        self.middle_scroll = QtWidgets.QScrollArea()
        self.main_vbox.addWidget(self.middle_scroll)
        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # Scroll area properties.
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        # QFormLayout in the middle frame.
        self.middle_formlayout = QtWidgets.QFormLayout()
        self.middle_widget.setLayout(self.middle_formlayout)


        #-----------------
        # Buttons frame. -
        #-----------------

        # Set the frame and its layout.
        self.buttons_frame = QtWidgets.QFrame()
        self.middle_formlayout.addRow(self.buttons_frame)
        self.buttons_hbox = QtWidgets.QHBoxLayout()
        self.buttons_frame.setLayout(self.buttons_hbox)

        # Build the control buttons.
        self.blast_select_all_button = QtWidgets.QPushButton(text="Select All")
        self.blast_select_all_button.clicked.connect(self.blast_select_all)
        self.blast_select_none_button = QtWidgets.QPushButton(text="Select None")
        self.blast_select_none_button.clicked.connect(self.blast_select_none)
        self.blast_select_n_button = QtWidgets.QPushButton(text="Select Top:")
        self.blast_select_n_button.clicked.connect(self.blast_select_n)
        for button in [self.blast_select_all_button, self.blast_select_none_button, self.blast_select_n_button]:
            self.buttons_hbox.addWidget(button)

        # Build the line-edit for selecting only top entries.
        self.blast_select_n_enf = PyMod_entryfield_qt(label_text="", value="10",
                                                      validate={'validator': 'integer',
                                                                'min': 1, 'max': 5000})
        self.blast_select_n_enf.entry.setFixedWidth(70)
        self.buttons_hbox.addWidget(self.blast_select_n_enf.entry)

        # Align to the left all these widgets.
        self.buttons_hbox.setAlignment(QtCore.Qt.AlignLeft)
        for button in [self.blast_select_all_button, self.blast_select_none_button, self.blast_select_n_button]:
            button.setFixedWidth(button.sizeHint().width()+30)

        #-----------------
        # Results frame. -
        #-----------------

        # Set the frame and its layout.
        self.results_frame = QtWidgets.QFrame()
        self.middle_formlayout.addRow(self.results_frame)
        self.results_grid = QtWidgets.QGridLayout()
        self.results_frame.setLayout(self.results_grid)

        # Calls a method which actually displays the similarity searches results.
        self.display_blast_hits()

        # Align the gridded widgets to the left.
        self.results_grid.setAlignment(QtCore.Qt.AlignLeft)
        self.results_grid.setHorizontalSpacing(30)


        #################
        # Bottom frame. #
        #################

        self.main_button = QtWidgets.QPushButton("Submit")
        self.main_button.clicked.connect(lambda a=None: self.protocol.blast_results_state())
        self.main_vbox.addWidget(self.main_button)
        self.main_button.setFixedWidth(self.main_button.sizeHint().width())


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


    def _get_window_title(self):
        return "%s Results" % self.protocol.blast_version_full_name


    def _get_upper_frame_title(self):
        title_text = ("%s Output for: %s\nFound %s sequences\nPlease Select the Sequences to"
                      " Import" % (self.protocol.blast_version_full_name,
                                   self.protocol.blast_query_element.compact_header,
                                   len(self.protocol.hsp_list)))
        return title_text


    def display_blast_hits(self):
        """
        This is used to display in the BLAST results window information for each
        hit and a checkbutton to select it for importing it inside PyMod.
        """

        # Shows the headers of the columns.
        headers_font_style = "%s; color: %s; font-weight: bold" % (small_font_style, highlight_color)
        self.blast_seq_label = QtWidgets.QLabel("Name")
        self.blast_seq_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.blast_seq_label, 0, 0)

        if not self.protocol.blast_version in hmmer_protocols_names:
            evalue_header_text = "E-Value"
        else:
            evalue_header_text = "E-Value (c-Evalue)"
        self.blast_e_val_label = QtWidgets.QLabel(evalue_header_text)
        self.blast_e_val_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.blast_e_val_label, 0, 1)

        if self.protocol.protocol_name != "hmmsearch":
            self.blast_iden_label = QtWidgets.QLabel("Identity")
            self.blast_iden_label.setStyleSheet(headers_font_style)
            self.results_grid.addWidget(self.blast_iden_label, 0, 2)

        self.query_span_label = QtWidgets.QLabel("Query span")
        self.query_span_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.query_span_label, 0, 3)

        self.subject_span_label = QtWidgets.QLabel("Subject span")
        self.subject_span_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.subject_span_label, 0, 4)


        # Displays in the results window the hsps found in the output file of the
        # similarity search program.
        self.blast_output_row = 1

        self.sbjct_checkbuttons_list = [] # List containing the checkbutton widgets.

        for hsp in self.protocol.hsp_list:

            # Hit name and checkbox.
            subject_name = self.protocol.get_subject_name(hsp)
            hsp_checkbox = QtWidgets.QCheckBox(subject_name)
            hsp_checkbox.setStyleSheet(small_font_style)
            self.sbjct_checkbuttons_list.append(hsp_checkbox)
            self.results_grid.addWidget(hsp_checkbox, self.blast_output_row, 0)

            # E-value info.
            if not self.protocol.blast_version in hmmer_protocols_names:
                evalue_label_text = "%.2e" % (self.protocol.get_hsp_evalue(hsp))
            else:
                evalue_tuple = self.protocol.get_hsp_evalue(hsp)
                evalue_label_text = "%.1e (%.1e)" % (evalue_tuple[0], evalue_tuple[1])
                # evalue_label_text = "%.1e" % (evalue_tuple[0])
            evalue_label = QtWidgets.QLabel(evalue_label_text)
            evalue_label.setStyleSheet(small_font_style)
            self.results_grid.addWidget(evalue_label, self.blast_output_row, 1)

            # HSP seqid info.
            if self.protocol.protocol_name != "hmmsearch":
                id_text = "%s/%s (%.1f%%)" % (hsp.pymod_info["identities"], int(hsp.pymod_info["matches"]), hsp.pymod_info["seqid"]*100)
                identities_label = QtWidgets.QLabel(id_text)
                identities_label.setStyleSheet(small_font_style)
                self.results_grid.addWidget(identities_label, self.blast_output_row, 2)

            # Query span info.
            span_info_text = "%s-%s (%.1f%%)" % (hsp.pymod_info["query_start"], hsp.pymod_info["query_end"], hsp.pymod_info["query_span"]*100)
            span_info_label = QtWidgets.QLabel(span_info_text)
            span_info_label.setStyleSheet(small_font_style)
            self.results_grid.addWidget(span_info_label, self.blast_output_row, 3)

            # Subject span info.
            hspan_info_text = "%s-%s" % (hsp.pymod_info["sbjct_start"], hsp.pymod_info["sbjct_end"])
            hspan_info_label = QtWidgets.QLabel(hspan_info_text)
            hspan_info_label.setStyleSheet(small_font_style)
            self.results_grid.addWidget(hspan_info_label, self.blast_output_row, 4)

            self.blast_output_row += 1


    def blast_select_all(self):
        for chk in self.sbjct_checkbuttons_list:
            chk.setChecked(True)


    def blast_select_none(self):
        for chk in self.sbjct_checkbuttons_list:
            chk.setChecked(False)


    def blast_select_n(self):
        try:
            select_top = int(self.blast_select_n_enf.getvalue())
            self.blast_select_none()
            count = 0
            for chk in self.sbjct_checkbuttons_list:
                chk.setChecked(True)
                count += 1
                if count == select_top:
                    break
        except ValueError: # Can not convert the input value in an integer.
            pass
