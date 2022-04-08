# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for running PSI-BLAST in PyMod.
"""

import os
import shutil

from Bio.Blast import NCBIXML

import pymod_lib.pymod_os_specific as pmos
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import Generic_BLAST_search
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import catch_protocol_exception


temp_output_dir_name = "__blast_temp__"

class PSI_BLAST_common:
    """
    A mixin class for using PSI-BLAST in other protocols (such as in psipred) other than the
    PSI-BLAST protocol itself.
    """

    def verify_valid_blast_dbdir(self, dbpath, remove_temp_files=False):
        """
        Checks if the folder specified in 'dbpath' contains a valid set of sequence database files.
        """
        if os.path.isfile(dbpath):
            return ""
        if remove_temp_files:
            for fn in os.listdir(dbpath):
                if fn == temp_output_dir_name:
                    fp = os.path.join(dbpath, fn)
                    if os.path.isfile(fp):
                        os.remove(fp)
                    elif os.path.isdir(fp):
                        shutil.rmtree(fp)

        dbprefix = self._get_blast_database_prefix(dbpath)
        return dbprefix

    def _get_blast_database_prefix(self, dbpath):
        database_files = []
        for fn in os.listdir(dbpath):
            if fn == ".DS_Store":
                pass
            elif fn.startswith("taxdb"):
                pass
            else:
                prefix = fn.split(".")[0]
                prefix.replace("_v5", "")
                database_files.append(prefix)
        return os.path.commonprefix(database_files)


    def execute_psiblast(self, ncbi_dir, db_path, query,
                               inclusion_ethresh=0.001, num_iterations=3,
                               evalue=None, max_target_seqs=None, num_alignments=None,
                               out=None, outfmt=None, out_pssm=None,
                               blast_version="psiblast"):
        """
        Execute the locally installed PSI-BLAST. Used when running PSI-BLAST through the 'PSI-BLAST'
        command on the plugin main menu or when predicting secondary structures with PSIPRED.
        """

        self.check_blast_version(blast_version)

        # Gests the prefix of the database folder.
        moved_to_db_dir = False
        try:
            dp_prefix = self.verify_valid_blast_dbdir(db_path)
            # Makes a temporary directory in the folder of the selected database.
            os.mkdir(os.path.join(db_path, temp_output_dir_name))
            # Copies the .fasta file of the query in the temporary folder.
            query_file_name = os.path.split(query)[1]
            shutil.copy(query, os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves in the database directory.
            os.chdir(db_path)
            moved_to_db_dir = True
            # Sets the input and  output file names.
            temp_query_shortcut = os.path.join(temp_output_dir_name, query_file_name)
            temp_out_file_shortcut = None
            if out != None:
                temp_out_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out)[1])
            temp_out_pssm_file_shortcut = None
            if out_pssm != None:
                temp_out_pssm_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out_pssm)[1])

            # Builds the PSI-BLAST commandline.
            psiblast_command = self.build_psiblast_commandline(
                ncbi_dir = ncbi_dir,
                db_path = dp_prefix,
                query = temp_query_shortcut,
                inclusion_ethresh = inclusion_ethresh,
                num_iterations = num_iterations,
                evalue = evalue,
                max_target_seqs = max_target_seqs,
                num_alignments = num_alignments,
                out = temp_out_file_shortcut,
                outfmt = outfmt,
                out_pssm = temp_out_pssm_file_shortcut,
                blast_version=blast_version)

            # Execute PSI-BLAST.
            self.pymod.execute_subprocess(psiblast_command)

            # Goes back to the original directory.
            os.chdir(self.pymod.current_project_dirpath)
            # Removes the query temp file.
            os.remove(os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves the temporary files in the originally specified output directory.
            for output_file in os.listdir(os.path.join(db_path, temp_output_dir_name)):
                new_output_filepath = os.path.join(os.path.split(query)[0], output_file)
                if os.path.isfile(new_output_filepath):
                    os.remove(new_output_filepath)
                shutil.move(os.path.join(db_path, temp_output_dir_name, output_file), os.path.split(query)[0])
            # Remove the temporary directory.
            shutil.rmtree(os.path.join(db_path, temp_output_dir_name))

        except Exception as e:
            # If something goes wrong while executing PSI-BLAST, go back to the project directory
            # and removes the temporary directory in the database folder, it it was built.
            if moved_to_db_dir:
                os.chdir(self.pymod.current_project_dirpath)
            if os.path.isdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
            raise e


    def build_psiblast_commandline(self, ncbi_dir, db_path, query,
                                   inclusion_ethresh=0.001, num_iterations=3,
                                   evalue=None, max_target_seqs=None, num_alignments=None,
                                   out=None, outfmt=None, out_pssm=None,
                                   blast_version="psiblast"):

        self.check_blast_version(blast_version)

        # blastdbcmd -db "\"Users\joeuser\My Documents\Downloads\mydb\"" -info
        # blastdbcmd -db ' "path with spaces/mydb" ' -info
        psiblast_path = pmos.build_commandline_path_string(os.path.join(ncbi_dir, pmos.get_exe_file_name(blast_version)))
        db_path = pmos.build_commandline_file_argument("db", db_path)
        query = pmos.build_commandline_file_argument("query", query)

        if inclusion_ethresh != None:
            inclusion_ethresh = " -inclusion_ethresh %s" % (inclusion_ethresh)
        else:
            inclusion_ethresh = ""

        if num_iterations != None:
            num_iterations = " -num_iterations %s" % (num_iterations)
        else:
            num_iterations = ""

        if evalue != None:
            evalue = " -evalue %s" % (evalue)
        else:
            evalue = ""

        if outfmt != None:
            outfmt = " -outfmt %s" % (outfmt) # 5 produces an .xml output file.
        else:
            outfmt = ""

        if out != None:
            out = pmos.build_commandline_file_argument("out", out)
        else:
            out = ""

        if max_target_seqs != None:
            max_target_seqs = " -max_target_seqs %s" % (max_target_seqs)
        else:
            max_target_seqs = ""

        if out_pssm != None:
            out_pssm = pmos.build_commandline_file_argument("out_pssm", out_pssm)
        else:
            out_pssm = ""

        if num_alignments != None:
            num_alignments = " -num_alignments %s" % (num_alignments)
        else:
            num_alignments = ""

        psiblast_command = (psiblast_path + db_path + query +
                            inclusion_ethresh + out + outfmt + out_pssm +
                            num_iterations + evalue + max_target_seqs +
                            num_alignments)

        return psiblast_command

    def check_blast_version(self, blast_version):
        if not blast_version in ("psiblast", "blastp"):
            raise KeyError("Unknown BLAST version: %s." % blast_version)


###################################################################################################
# PSI-BLAST protocol.                                                                             #
###################################################################################################

class PSI_BLAST_search(Generic_BLAST_search, PSI_BLAST_common):

    blast_version = "psi-blast"
    protocol_name = blast_version
    exe_filename = "psiblast"
    # PSI-BLAST minimum inclusion E-value.
    min_inclusion_eval_default = 0.005

    def additional_initialization(self):
        Generic_BLAST_search.additional_initialization(self)
        self.tool = self.pymod.blast_plus


    def get_blast_window_class_qt(self):
        return PSI_BLAST_options_window_qt


    def check_blast_input_parameters(self):
        # Check if a valid database for PSI-BLAST was provided.
        db_full_path = self.get_database_from_gui()
        if db_full_path == None:
            title = "Input Error"
            message = "Please choose a valid database."
            self.pymod.main_window.show_error_message(title, message)
            return False
        if not self.verify_valid_blast_dbdir(db_full_path, remove_temp_files=True):
            title = "Input Error"
            message = "The database '%s' directory does not contain a valid set of database files." % (db_full_path)
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
        return self.run_psiblast()


    def get_blast_record(self, result_handle):
        """
        Convert it to a list because when a using .parse(), Biopython returns a generator.
        """
        records = list(NCBIXML.parse(result_handle))
        return records[0]


    #################################################################
    # PSI-BLAST specific.                                           #
    #################################################################

    def run_psiblast(self):
        """
        Launches a standalone version of PSI-BLAST installed locally when using the PSI-BLAST
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
                "inclusion_ethresh": self.evalue_inclusion_cutoff,
                "outfmt": 5,
                "out": os.path.join(self.output_directory, self.xml_blast_output_file_name),
                "num_iterations": self.iterations,
                "evalue": self.evalue_cutoff,
                "max_target_seqs": self.max_hsp_num}

        if not self.pymod.use_protocol_threads:
            self.execute_psiblast(**args)

        else:
            label_text = "Running PSI-BLAST. Please wait for the process to complete..."
            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.execute_psiblast,
                                            args=args,
                                            title="Running PSI-BLAST",
                                            label_text=label_text)
            p_dialog.exec_()

        # If everything went ok, return 'True', so that the results window can be opened.
        return True


    def get_options_from_gui_specific(self):
        self.get_options_from_gui_blast()
        self.iterations = self.blast_options_window.psiblast_iterations_enf.getvalue()
        if self.blast_options_window.showing_advanced_widgets:
            self.evalue_inclusion_cutoff = self.blast_options_window.psiblast_eval_threshold_enf.getvalue()
        else:
            self.evalue_inclusion_cutoff = self.min_inclusion_eval_default


    def build_blast_db_list(self):
        """
        Generates a list of dictionaries each containing information about the sequence databases
        present in the default BLAST database directory.
        """

        blast_db_dir = self.tool["database_dir_path"].get_value()
        list_of_databases_directories = []

        # If there are multiple directories containing database files with the same prefixes, this
        # will rename their prefixes so that the database radioselect will not have multiple buttons
        # with the same name.
        def get_new_prefix(prefix, list_of_databases_directories, n=1, prefix_root=None):
            if prefix_root == None:
                prefix_root = prefix
            if prefix in [dbd["prefix"] for dbd in list_of_databases_directories]:
                new_prefix = prefix_root + "-" + str(n)
                return get_new_prefix(new_prefix, list_of_databases_directories, n+1, prefix_root)
            else:
                return prefix

        for path in os.listdir(blast_db_dir):
            full_path = os.path.join(blast_db_dir, path)
            prefix = self.verify_valid_blast_dbdir(full_path, remove_temp_files=True)
            if prefix:
                prefix = get_new_prefix(prefix, list_of_databases_directories)
                dbd = {"full-path": full_path, "prefix": prefix}
                list_of_databases_directories.append(dbd)

        return list_of_databases_directories


    def get_database_from_gui(self):
        button_name = self.blast_options_window.psiblast_database_rds.getvalue()
        for dbd in self.databases_directories_list:
            if dbd["prefix"] == button_name:
                return dbd["full-path"]
        return None


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_entryfield_qt, PyMod_radioselect_qt
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import BLAST_base_options_window_qt


class PSI_BLAST_options_window_qt(BLAST_base_options_window_qt):
    """
    Window for PSI-BLAST searches.
    """

    def build_algorithm_standard_options_widgets(self):

        # Makes the user chose the folder where the BLAST database files are stored locally.
        # A list containing information about the databases present in PyMod BLAST database
        # folder.
        db_list = [k["prefix"] for k in self.protocol.databases_directories_list]
        self.psiblast_database_rds = PyMod_radioselect_qt(label_text="Database Selection",
                                                          buttons=db_list)
        self.middle_formlayout.add_widget_to_align(self.psiblast_database_rds)


        # Number of PSI-BLAST iterations.
        if self.protocol.blast_version == "psi-blast":
            self.psiblast_iterations_enf = PyMod_entryfield_qt(label_text="PSI-BLAST Iterations",
                                                               value='3',
                                                               validate={'validator': 'integer',
                                                                         'min': 1, 'max': 10})
            self.middle_formlayout.add_widget_to_align(self.psiblast_iterations_enf, validate=True)


    def build_algorithm_advanced_options_widgets(self):

        if self.protocol.blast_version == "psi-blast":
            self.psiblast_eval_threshold_enf = PyMod_entryfield_qt(
                label_text="Inclusion E-value",
                value=str(self.protocol.min_inclusion_eval_default),
                validate={'validator': 'real', 'min': 0.0, 'max': 1000.0})
            self.middle_formlayout.add_widget_to_align(self.psiblast_eval_threshold_enf,
                                                       advanced_option=True, validate=True)

        # Use current cluster for PSI-BLAST PSSM.
        # if self.blast_query_element.is_child:
        #     self.use_current_pssm_rds = PyMod_radioselect(self.midframe, label_text = 'Use current cluster as PSSM')
        #     for text in ('Yes', 'No'):
        #         self.use_current_pssm_rds.add(text)
        #     self.use_current_pssm_rds.setvalue('No')
        #     # self.use_current_pssm_rds.pack(side = 'top', padx = 10, pady = 10, anchor="w")
        #     self.add_widget_to_align(self.use_current_pssm_rds)
        #     self.add_advanced_widget(self.use_current_pssm_rds)
