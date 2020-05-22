# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for performing phmmer searches in PyMod. It mainly builds up on the code used to perform BLAST
searches.
"""

import os

from Bio import SearchIO

from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import Generic_BLAST_search, BLAST_base_options_window_qt
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import catch_protocol_exception
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_radioselect_qt


################################################################################################
# PHMMER.                                                                                      #
################################################################################################

class PHMMER_search(Generic_BLAST_search):

    blast_version = "phmmer"
    protocol_name = blast_version
    exe_filename = "phmmer"
    e_value_threshold_default = 1.0


    def additional_initialization(self):
        Generic_BLAST_search.additional_initialization(self)
        self.tool = self.pymod.hmmer_tool

    def get_blast_window_class_qt(self):
        return Phmmer_options_window_qt


    def build_blast_db_list(self):
        db_dirpath = self.pymod.hmmer_tool["database_dir_path"].get_value()
        db_list = [filename for filename in os.listdir(db_dirpath) if filename.endswith(".fasta")]
        return db_list


    def check_blast_input_parameters(self):
        # Check if a valid database for PSI-BLAST was provided.
        db_full_path = self.get_database_from_gui()
        if db_full_path == None:
            title = "Input Error"
            message = "Please choose a valid database."
            self.pymod.main_window.show_error_message(title, message)
            return False
        # Check all the other input fields.
        try:
            self.blast_options_window.check_general_input()
        except Exception as general_input_error:
            title = "Input Error"
            message = str(general_input_error)
            self.pymod.main_window.show_error_message(title, message)
            return False
        if not self.check_min_max_seqid_input():
            return False
        # Returns 'True' only if all input parameters are valid.
        return True


    @catch_protocol_exception
    def run_blast_program(self):
        return self.run_hmmer()


    def run_hmmer(self):
        """
        Launches a standalone version of PHMMER installed locally.
        """

        # Builds a temporary file with the sequence of the query needed.
        query_file_name = "query"

        if self.protocol_name != "hmmsearch": # Saves a sequence file with only the query.
            self.pymod.build_sequence_file([self.blast_query_element], query_file_name, file_format="fasta",
                                           remove_indels=True, new_directory=self.output_directory)
        else: # Saves the entire query alignment.
            self.pymod.build_sequence_file(self.blast_query_element.get_children(), query_file_name,
                                           file_format="fasta", remove_indels=False,
                                           new_directory=self.output_directory)


        # Sets some parameters needed to run PHMMER.
        if self.protocol_name != "hmmsearch": # Get only one filepath.
            exe_filepath = self.get_similariry_search_exe_path()[0]
        else: # Get all filepaths.
            exe_filepath_dict = self.get_similariry_search_exe_path(return_dict=True)

        out_filepath = os.path.join(self.output_directory, "%s_out.txt" % self.protocol_name)
        self.result_filepath = out_filepath

        query_filepath = os.path.join(self.output_directory, "query.fasta")
        db_dirpath = self.pymod.hmmer_tool["database_dir_path"].get_value()
        db_filepath = os.path.join(db_dirpath, self.db_filename)

        args = {"query_filepath": query_filepath,
                "out_filepath": out_filepath,
                "db_filepath": db_filepath,
                "report_evalue": self.report_evalue}
        args.update(self.additional_params)

        if self.protocol_name != "hmmsearch": # Adds only one exec filepath.
            args["exe_filepath"] = exe_filepath
        else: # Adds two exec filepaths.
            args["exe_filepath"] = exe_filepath_dict[self.all_exe_filenames[1]]
            args["hmmbuild_exe_filepath"] = exe_filepath_dict[self.all_exe_filenames[0]]


        # Actually runs phmmer.
        if not self.pymod.use_protocol_threads:
            self.execute_hmmer_program(**args)

        else:
            label_text = "Running %s. Please wait for the process to complete..." % self.blast_version_full_name
            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.execute_hmmer_program,
                                            args=args,
                                            title="Running %s" % self.blast_version_full_name,
                                            label_text=label_text)
            p_dialog.exec_()

        return True


    def execute_hmmer_program(self, query_filepath, out_filepath, db_filepath, exe_filepath, report_evalue=10e-6):
        """
        Execute the locally installed PHMMER. Used when running PHMMER through the 'PHMMER'.
        """

        phmmer_command_ls = [exe_filepath, "-o", out_filepath, "--domE", report_evalue,
                             query_filepath, db_filepath]
        self.pymod.new_execute_subprocess(phmmer_command_ls)


    def get_options_from_gui_specific(self):
        self.db_filename = self.get_database_from_gui()
        self.report_evalue = self.blast_options_window.e_value_threshold_enf.getvalue()
        if self.blast_options_window.showing_advanced_widgets:
            self.min_id = self.blast_options_window.min_id_enf.getvalue()
            self.max_id = self.blast_options_window.max_id_enf.getvalue()
            self.min_coverage = self.blast_options_window.min_coverage_enf.getvalue()
        self.additional_params = self.blast_options_window.get_additional_hmmer_parameters()


    def get_database_from_gui(self):
        return self.blast_options_window.get_phmmer_database()


    def get_search_record(self):
        # Parse HMMER's output.
        with open(self.result_filepath, "r") as result_handle:
            phmmer_query_result = SearchIO.read(result_handle, 'hmmer3-text')

        # Gets the hmmsearch profile length.
        if self.protocol_name == "hmmsearch":
            self.profile_length = None
            with open(self.result_filepath, "r") as result_handle:
                for line in result_handle:
                    # The length of the profile is stored in a line of the output
                    # file starting with "Query:".
                    if line.startswith("Query:"):
                        try:
                            self.profile_length = int(line.split()[-1].split("=")[1].replace("]", ""))
                        except (IndexError, TypeError):
                            pass
            if self.profile_length is None:
                self.profile_length = len(self.blast_query_element.get_children()[0].my_sequence)

        return phmmer_query_result


    def get_hsp_query_seq(self, hsp_obj):
        return str(hsp_obj.query.seq).upper().replace(".", "-")

    def get_hsp_subj_seq(self, hsp_obj):
        return str(hsp_obj.hit.seq).upper()

    def get_hsp_evalue(self, hsp):
        return hsp.evalue, hsp.evalue_cond

    def get_subj_start(self, hsp):
        return hsp.hit_start

    def build_hsp_list(self):
        hsp_list = []
        for hsp in self.search_record.hsps:
                hsp._title = hsp.hit.id
                # Gets additional information on the hsp (such as its the id % and query span).
                hsp.pymod_info = self.get_hsp_info(hsp)
                hsp_list.append(hsp)
        # Sort by evalue.
        hsp_list = list(sorted(hsp_list, key=lambda h: h.evalue))
        return hsp_list


    def finish_blast_search(self):
        output_filename = self.get_blast_output_basename() + ".txt"
        try:
            os.rename(self.result_filepath, os.path.join(self.output_directory, output_filename))
        except IOError:
            pass


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Phmmer_options_window_qt(BLAST_base_options_window_qt):
    """
    Window for PHMMER searches.
    """

    def build_algorithm_standard_options_widgets(self):

        self.hmmer_database_labels_dict = {}

        # A list containing information about the databases present in PyMod hmmer database folder.
        db_list = []
        for db in self.protocol.databases_directories_list:
            db_label = "".join(db.split(".")[0:-1])
            db_list.append(db_label)
            self.hmmer_database_labels_dict[db_label] = db

        # Makes the user choose the folder where the hmmer database files are stored locally.
        self.hmmer_database_rds = PyMod_radioselect_qt(label_text="Database Selection",
                                                       buttons=db_list)
        self.middle_formlayout.add_widget_to_align(self.hmmer_database_rds)


        self.build_additional_hmmer_widgets()


    def build_additional_hmmer_widgets(self):
        pass


    def build_algorithm_advanced_options_widgets(self):
        pass


    def get_phmmer_database(self):
        return self.hmmer_database_labels_dict.get(self.hmmer_database_rds.getvalue(), None)


    def get_additional_hmmer_parameters(self):
        return {}
