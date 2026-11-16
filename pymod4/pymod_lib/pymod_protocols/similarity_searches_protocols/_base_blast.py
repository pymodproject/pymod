# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module with a base class for similarity searches protocols. It is used to build classes for every
similarity search protocol in PyMod, that is, (PSI-)BLAST protocols and phmmer/jackhmmer protocols.
"""

import os
import re

from pymod_lib.pymod_os_specific import get_exe_file_name, clean_file_name
from pymod_lib.pymod_vars import algs_full_names_dict
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.base_protocols import AFDB_common
from pymod_lib.pymod_protocols.alignment_protocols.clustalw import Clustalw_regular_alignment, Clustalw_profile_alignment
from pymod_lib.pymod_protocols.modeling_protocols.modeling_protocols_common import Fetch_models_file
from pymod_lib.pymod_seq.seq_star_alignment import save_cstar_alignment
from pymod_lib.pymod_threading import Protocol_exec_dialog

from pymol.Qt import QtWidgets, QtCore
from pymod_lib.pymod_gui.hits_table_qt import (Hits_table_qt,
                                              span_from_fraction)
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_entryfield_qt,
                                                          small_font_style,
                                                          highlight_color,
                                                          small_font_style_colored)


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
    # Whether this is an AlphaFold DB search. Set to True via launch_from_gui only
    # for AFDB searches; the class-level default keeps subclasses that don't set it
    # (e.g. hmmsearch) from raising AttributeError in check_blast_results().
    AF = False

    # Which fetcher downloads the structure of each hit when AF is on. A
    # search which takes its structures from another database of predicted
    # models -- the ESM Atlas -- subclasses this and swaps the fetcher.
    models_fetcher_class = Fetch_models_file

    # What to call the database of predicted structures in front of the user. A
    # search which takes its structures elsewhere sets its own.
    models_database_name = "AFDB"

    # The same database written as it belongs in a sentence: the column of the
    # results window, the exported tables and the name of the cluster the hits end
    # up in all read "<this> model".
    models_label = "AlphaFold"

    def additional_initialization(self):
        self.blast_version_full_name = algs_full_names_dict[self.blast_version]


    def launch_from_gui(self, AF = False):

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
        self.AF = AF
        if AF:
            self.build_blast_window_for_AFDB_search()
        else:
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


    def build_blast_window_for_AFDB_search(self):
        """
        Builds a window containing the widget necessary to define the options for BLAST and
        PSI-BLAST searches.
        """

        blast_window_class_qt = self.get_blast_window_class_qt(AF = True)
        self.blast_options_window = blast_window_class_qt(self.pymod.main_window, protocol=self,
            title= "Local BLAST Search - Options",
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

        # Check results
        if blast_status and self.check_blast_results():

            # Displays the window with results.
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
        if self.AF:
            label_text = ("Running %s search. Please wait for the process to complete..."
                          % self.models_database_name)
            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.build_hsp_list,
                                            args=(),
                                            title="Running %s search"
                                                  % self.models_database_name,
                                            label_text=label_text)
            p_dialog.exec_()

        else:
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

        # The model of a database of predicted structures which corresponds to this
        # hit. Which database is asked depends on the search, so the name is
        # resolved by a method: a search which takes its structures from elsewhere
        # then does not have to ask the AlphaFold database about every hit.
        af_mod = self.resolve_model_for_hit(hsp) if self.AF else None
        # Returns the information.
        additional_infor_dict = {"identities": identities, "seqid": idp,
                                 "query_span": query_span, "matches": matches,
                                 "query_start": qs, "query_end": qe,
                                 "sbjct_start": hs, "sbjct_end": he,
                                 "af_model": af_mod}

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

    def resolve_model_for_hit(self, hsp):
        """
        The name of the predicted model which corresponds to a hit, or None.

        By default the AlphaFold database is asked, by the UniProt accession of the hit,
        which is what the "Search By: Sequence" mode of the AlphaFold retrieval needs. A
        search which takes its structures from another database overrides this.
        """

        resolver = AFDB_common(self, hsp)
        resolver.init_run_AFDB_uniprot_id_search()
        return getattr(resolver, "af_model_id", None)


    def show_blast_output_window(self):
        """
        Displays the window with results from BLAST in a new window.
        """
        self.blast_output_window = Similarity_searches_results_window_qt(parent=self.pymod.main_window,
                                                           protocol=self,
                                                           AF = self.AF)
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
        self.hsp_list = list(sorted(hsp_list, key=lambda h: h.expect))
        return self.hsp_list

    def get_subject_name(self, hsp):
        full_title = hsp._title
        if len(full_title) > 100:
            return full_title[0:100] + "..."
        else:
            return full_title


    # A hit is named by an accession and described by what follows it. Which of
    # the two a search puts in the name of the hit depends on the program and on
    # the database: a BLAST alignment carries both in one string, HMMER keeps
    # them apart, and a database of the PDB gives an accession with nothing after
    # it at all. The window asks for them separately and shows each in its own
    # column, so that a list of accessions is not the only thing on the screen.

    def get_subject_id(self, hsp):
        """The accession of a hit, without what describes it."""

        name = str(getattr(hsp, "_title", "") or "")
        if getattr(hsp, "_description", None):
            # The search kept the two apart, so the name is already the accession.
            return name.strip()
        return name.split(None, 1)[0] if name.strip() else name


    def get_subject_description(self, hsp):
        """What a search says about a hit besides its accession, or ''."""

        described = getattr(hsp, "_description", None)
        if described:
            return str(described).strip()
        name = str(getattr(hsp, "_title", "") or "")
        pieces = name.split(None, 1)
        return pieces[1].strip() if len(pieces) > 1 else ""

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
                if self.AF:
                    self.hsp_imported_from_blast.append({"hsp": hsp, "title": hsp.pymod_info["af_model"]})
                else:
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

            if self.AF:
                cluster_type = "af-cluster"
            else:
                cluster_type = "blast-cluster"
            # Builds the "BLAST search" cluster element.
            new_blast_cluster = self.pymod.add_new_cluster_to_pymod(
                cluster_type= cluster_type,
                cluster_label=self.models_label,
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

            # Only perform sequence comparison check for protein sequences, not nucleic acids
            if self.blast_query_element.polymer_type == "protein":
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

        if self.AF:

            # ##### MODIFIED #####
            # self.pymod.af_modeling_session_list[self.pymod.current_af_session]["import_mode_choice"] = self.blast_output_window.import_mode_cb.currentText()
            # ##### END #####

            # The user is asked whether to import the whole structure or just the hit fragment
            if self.blast_output_window.import_mode_cb.currentText() == "Whole Structure":
                import_mode = "multiple-chains"
            else:
                import_mode = "single-chain"

            self.blast_output_window.destroy()

            for pymod_elements in hsp_elements:
                # Download AF model
                fp = self.models_fetcher_class(self.pymod)
                fp.initialize_from_gui("single", pymod_elements)
                fp.import_mode = import_mode
                fp.fetch_models_files()

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

    def __init__(self, parent, protocol, AF):

        super(Similarity_searches_results_window_qt, self).__init__(parent)
        self.protocol = protocol
        self._described = None

        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle(self._get_window_title())
        ##### MODIFIED #####
        self.setGeometry(300, 200, 1300, 600)
        ##### END #####

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

        if AF:

            # The download of the 3D-structure is automatically triggered after having checked the BLAST outputs of interest
            self.import_mode_label = QtWidgets.QLabel("Structure import mode: ")
            self.import_mode_cb = QtWidgets.QComboBox()
            self.import_mode_cb.addItems(["Whole Structure", "Only hit fragment"])

            self.buttons_hbox.addWidget(self.import_mode_label)
            self.buttons_hbox.addWidget(self.import_mode_cb)

        ##### MODIFIED #####
        # The search bar sits over the table, in the table's own widget.

        # Export buttons
        self.export_xlsx_button = QtWidgets.QPushButton("Export XLSX")
        self.export_xlsx_button.clicked.connect(self.export_to_xlsx)
        self.export_csv_button = QtWidgets.QPushButton("Export CSV")
        self.export_csv_button.clicked.connect(self.export_to_csv)

        # Layout
        self.buttons_hbox.addWidget(self.export_xlsx_button)
        self.buttons_hbox.addWidget(self.export_csv_button)
        ##### END #####

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

        # The hits, in the table which every search of PyMod shows them in: a
        # search bar over a table with a box to tick on each row, and the whole
        # of a long name written out underneath.
        self.hits_table = Hits_table_qt(
            parent=self, columns=self.result_columns(),
            detail_column="Name",
            search_placeholder="part of a hit name or of its description...")
        self.middle_formlayout.addRow(self.hits_table)

        # Calls a method which actually displays the similarity searches results.
        self.display_blast_hits()


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


    def result_columns(self):
        """
        The columns of this search. They differ a little between the programs:
        hmmsearch reports no identity, and a search of a model database has a
        column for the model of each hit.
        """

        columns = [("", 0), ("Name", 260 if self.hits_are_described() else 420)]
        if self.hits_are_described():
            columns.append(("Description", 340))
        if self.protocol.blast_version not in hmmer_protocols_names:
            columns.append(("E-Value", 0))
        else:
            columns.append(("E-Value (c-Evalue)", 0))
        if self.protocol.protocol_name != "hmmsearch":
            columns.append(("Identity", 0))
        columns.append(("Query span", 150))
        columns.append(("Subject span", 0))
        if self.protocol.AF:
            columns.append(("%s model" % self.protocol.models_label, 0))
        return tuple(columns)


    # A protocol which does not answer these -- one which does not descend from
    # Generic_BLAST_search -- gets the reading which works off the hit itself,
    # rather than an error.

    def subject_id(self, hsp):
        reader = getattr(self.protocol, "get_subject_id", None)
        if reader is not None:
            return reader(hsp)
        return Generic_BLAST_search.get_subject_id(self.protocol, hsp)

    def subject_description(self, hsp):
        reader = getattr(self.protocol, "get_subject_description", None)
        if reader is not None:
            return reader(hsp)
        return Generic_BLAST_search.get_subject_description(self.protocol, hsp)


    def hits_are_described(self):
        """
        Whether this search said anything about its hits besides their accession.

        A search of a database of the PDB gives "5FYK_A" and nothing else, and a
        column of empty cells beside a column of accessions helps nobody; a
        search of UniProt gives a name for every hit, and that is the column
        worth having. Asked of the hits themselves rather than of the database,
        which is what the answer really depends on.
        """

        if self._described is None:
            self._described = any(self.subject_description(hsp)
                                  for hsp in getattr(self.protocol, "hsp_list", []))
        return self._described


    def display_blast_hits(self):
        """
        Fills the table with one row per hit: the figures which say how good it
        is, and a box to tick to import it into PyMod.
        """

        # What the protocols read to know which hits were ticked. They read it as
        # a list of checkboxes in the order of the hits, which is what it still
        # is; the table keeps them in that order.
        self.sbjct_checkbuttons_list = []

        for hsp in self.protocol.hsp_list:

            subject_name = self.protocol.get_subject_name(hsp)
            values = [self.subject_id(hsp)]
            if self.hits_are_described():
                values.append(self.subject_description(hsp))

            if not self.protocol.blast_version in hmmer_protocols_names:
                values.append("%.2e" % (self.protocol.get_hsp_evalue(hsp)))
            else:
                evalue_tuple = self.protocol.get_hsp_evalue(hsp)
                values.append("%.1e (%.1e)" % (evalue_tuple[0], evalue_tuple[1]))

            if self.protocol.protocol_name != "hmmsearch":
                values.append("%s/%s (%.1f%%)"
                              % (hsp.pymod_info["identities"],
                                 int(hsp.pymod_info["matches"]),
                                 hsp.pymod_info["seqid"]*100))

            values.append("%s-%s (%.1f%%)" % (hsp.pymod_info["query_start"],
                                              hsp.pymod_info["query_end"],
                                              hsp.pymod_info["query_span"]*100))
            values.append("%s-%s" % (hsp.pymod_info["sbjct_start"],
                                     hsp.pymod_info["sbjct_end"]))
            if self.protocol.AF:
                values.append("%s" % (hsp.pymod_info["af_model"]))

            # The span of the query is drawn as a bar as well as written: which
            # part of the query a hit covers is what says whether it answers the
            # whole question or one domain of it. The search reports the covered
            # fraction rather than the length of the query, which is enough.
            spans = {"Query span": span_from_fraction(hsp.pymod_info["query_start"],
                                                      hsp.pymod_info["query_end"],
                                                      hsp.pymod_info["query_span"])}

            # The strip underneath and the search bar both take the whole of it:
            # a hit is looked for by its accession as often as by a word of what
            # describes it, and the accession alone is what the column shows.
            description = self.subject_description(hsp)
            whole = ("%s \u2014 %s" % (self.subject_id(hsp), description)
                     if description else subject_name)

            self.hits_table.add_row(values, payload=hsp, spans=spans,
                                    detail=whole, searchable=whole)

            self.sbjct_checkbuttons_list.append(self.hits_table.rows[-1][0])

        self.hits_table.finish()

        # A hit is looked up in its database by double clicking its row: the name
        # used to be a link, and a link in a table cell is a widget per row.
        self.hits_table.table.cellDoubleClicked.connect(self.open_link_of_row)


    def open_link_of_row(self, row, column=0):
        hsp = self.hits_table.payload_of_row(row)
        if hsp is not None:
            self.open_database_link(self.protocol.get_subject_name(hsp))


    ###############################################################################################
    # The search bar and the buttons which select hits.                                           #
    ###############################################################################################

    @property
    def searchbar(self):
        return self.hits_table.searchbar

    @property
    def search_count_label(self):
        return self.hits_table.search_count_label

    def update_display(self):
        self.hits_table.update_display()

    def shown_checkboxes(self):
        return self.hits_table.shown_checkboxes()

    # These act on the hits which are shown: with an empty search bar that is
    # every hit, and with something typed in it, what was searched for.

    def blast_select_all(self):
        self.hits_table.select_all()

    def blast_select_none(self):
        self.hits_table.select_none()

    def blast_select_n(self):
        try:
            self.hits_table.select_first(int(self.blast_select_n_enf.getvalue()))
        except ValueError: # Can not convert the input value in an integer.
            pass


    def open_database_link(self, subject_name):
        """Open the appropriate database link based on the subject name format."""
        import webbrowser
        import re

        def get_ncbi_identifier(subject_name):
            # Try to find accessions like gb|AAB23349.1| or ref|NP_00125.1|
            match = re.search(r'(gb|ref|emb)\|([^|]+)\|', subject_name)
            if match:
                return match.group(2)
            
            # Try to get GI number
            match = re.search(r'gi\|(\\d+)', subject_name)
            if match:
                return match.group(1)

            # Fallback to the first part of the name
            return subject_name.split(' ')[0]

        try:
            database = ""
            if hasattr(self.protocol, 'get_ncbiblast_database'):
                database = self.protocol.get_ncbiblast_database()

            # UniProt/SwissProt
            if database == "swissprot" or subject_name.startswith('sp|') or subject_name.startswith('tr|'):
                try:
                    uniprot_id = subject_name.split('|')[1].split('.')[0]
                    url = f"https://www.uniprot.org/uniprot/{uniprot_id}"
                    webbrowser.open(url)
                    return
                except IndexError:
                    uniprot_pattern = r'([A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}|[OPQ][0-9][A-Z0-9]{3}[0-9])'
                    match = re.search(uniprot_pattern, subject_name)
                    if match:
                        uniprot_id = match.group(1).split('.')[0]
                        url = f"https://www.uniprot.org/uniprot/{uniprot_id}"
                        webbrowser.open(url)
                        return

            # PDB
            if subject_name.startswith('pdb|'):
                pdb_id = subject_name.split('|')[1]
                url = f"https://www.rcsb.org/structure/{pdb_id}"
                webbrowser.open(url)
                return

            # NCBI search for everything else (nucleotide and other protein)
            identifier = get_ncbi_identifier(subject_name)
            url = f"https://www.ncbi.nlm.nih.gov/search/all/?term={identifier}"
            webbrowser.open(url)
            return
            
        except Exception as e:
            self.protocol.pymod.main_window.show_error_message(
                "Link Error", 
                f"Error opening link: {str(e)}"
            )

    def export_to_csv(self):
        """Export BLAST results to CSV format."""
        try:
            import csv
            from pymol.Qt import QtWidgets
            
            # Get save file path
            file_path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, 
                "Export BLAST Results to CSV",
                f"{self.protocol.blast_version_full_name}_results.csv",
                "CSV files (*.csv)"
            )
            
            if not file_path:
                return
            
            with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write header
                header = ['Name']
                if self.hits_are_described():
                    header.append('Description')
                header.append('E-Value')
                if self.protocol.protocol_name != "hmmsearch":
                    header.append('Identity')
                header.extend(['Query Span', 'Subject Span'])
                if self.protocol.AF:
                    header.append('%s Model' % self.protocol.models_label)
                writer.writerow(header)
                
                # Write data
                for hsp in self.protocol.hsp_list:
                    row = [self.subject_id(hsp)]
                    if self.hits_are_described():
                        row.append(self.subject_description(hsp))
                    
                    # E-value
                    if not self.protocol.blast_version in hmmer_protocols_names:
                        row.append("%.2e" % self.protocol.get_hsp_evalue(hsp))
                    else:
                        evalue_tuple = self.protocol.get_hsp_evalue(hsp)
                        row.append("%.1e (%.1e)" % (evalue_tuple[0], evalue_tuple[1]))
                    
                    # Identity
                    if self.protocol.protocol_name != "hmmsearch":
                        row.append("%s/%s (%.1f%%)" % (
                            hsp.pymod_info["identities"], 
                            int(hsp.pymod_info["matches"]), 
                            hsp.pymod_info["seqid"]*100
                        ))
                    
                    # Query span
                    row.append("%s-%s (%.1f%%)" % (
                        hsp.pymod_info["query_start"], 
                        hsp.pymod_info["query_end"], 
                        hsp.pymod_info["query_span"]*100
                    ))
                    
                    # Subject span
                    row.append("%s-%s" % (
                        hsp.pymod_info["sbjct_start"], 
                        hsp.pymod_info["sbjct_end"]
                    ))
                    
                    # The model of this hit
                    if self.protocol.AF:
                        row.append(hsp.pymod_info["af_model"])
                    
                    writer.writerow(row)
            
            self.protocol.pymod.main_window.show_info_message(
                "Export Complete", 
                f"Results exported to {file_path}"
            )
            
        except Exception as e:
            self.protocol.pymod.main_window.show_error_message(
                "Export Error", 
                f"Error exporting CSV: {str(e)}"
            )

    def export_to_xlsx(self):
        """Export BLAST results to Excel format."""
        try:
            # Try to import openpyxl
            try:
                from openpyxl import Workbook
                from openpyxl.styles import Font, PatternFill
            except ImportError:
                self.protocol.pymod.main_window.show_error_message(
                    "Missing Dependency", 
                    "openpyxl is required for Excel export. Install it with: pip install openpyxl"
                )
                return
            
            from pymol.Qt import QtWidgets
            
            # Get save file path
            file_path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, 
                "Export BLAST Results to Excel",
                f"{self.protocol.blast_version_full_name}_results.xlsx",
                "Excel files (*.xlsx)"
            )
            
            if not file_path:
                return
            
            # Create workbook and worksheet
            wb = Workbook()
            ws = wb.active
            ws.title = f"{self.protocol.blast_version_full_name} Results"
            
            # Header style
            header_font = Font(bold=True)
            header_fill = PatternFill(start_color="366092", end_color="366092", fill_type="solid")
            
            # Write header
            header = ['Name']
            if self.hits_are_described():
                header.append('Description')
            header.append('E-Value')
            if self.protocol.protocol_name != "hmmsearch":
                header.append('Identity')
            header.extend(['Query Span', 'Subject Span'])
            if self.protocol.AF:
                header.append('%s Model' % self.protocol.models_label)
            
            for col, title in enumerate(header, 1):
                cell = ws.cell(row=1, column=col, value=title)
                cell.font = header_font
                cell.fill = header_fill
            
            # Write data
            for row_idx, hsp in enumerate(self.protocol.hsp_list, 2):
                col = 1
                
                # Name, and what describes it when the search said anything.
                ws.cell(row=row_idx, column=col, value=self.subject_id(hsp))
                col += 1
                if self.hits_are_described():
                    ws.cell(row=row_idx, column=col,
                            value=self.subject_description(hsp))
                    col += 1
                
                # E-value
                if not self.protocol.blast_version in hmmer_protocols_names:
                    ws.cell(row=row_idx, column=col, value="%.2e" % self.protocol.get_hsp_evalue(hsp))
                else:
                    evalue_tuple = self.protocol.get_hsp_evalue(hsp)
                    ws.cell(row=row_idx, column=col, value="%.1e (%.1e)" % (evalue_tuple[0], evalue_tuple[1]))
                col += 1
                
                # Identity
                if self.protocol.protocol_name != "hmmsearch":
                    ws.cell(row=row_idx, column=col, value="%s/%s (%.1f%%)" % (
                        hsp.pymod_info["identities"], 
                        int(hsp.pymod_info["matches"]), 
                        hsp.pymod_info["seqid"]*100
                    ))
                    col += 1
                
                # Query span
                ws.cell(row=row_idx, column=col, value="%s-%s (%.1f%%)" % (
                    hsp.pymod_info["query_start"], 
                    hsp.pymod_info["query_end"], 
                    hsp.pymod_info["query_span"]*100
                ))
                col += 1
                
                # Subject span
                ws.cell(row=row_idx, column=col, value="%s-%s" % (
                    hsp.pymod_info["sbjct_start"], 
                    hsp.pymod_info["sbjct_end"]
                ))
                col += 1
                
                # The model of this hit
                if self.protocol.AF:
                    ws.cell(row=row_idx, column=col, value=hsp.pymod_info["af_model"])
            
            # Auto-adjust column widths
            for column in ws.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = min(max_length + 2, 50)
                ws.column_dimensions[column_letter].width = adjusted_width
            
            # Save workbook
            wb.save(file_path)
            
            self.protocol.pymod.main_window.show_info_message(
                "Export Complete", 
                f"Results exported to {file_path}"
            )
            
        except Exception as e:
            self.protocol.pymod.main_window.show_error_message(
                "Export Error", 
                f"Error exporting Excel: {str(e)}"
            )
    ##### END #####
