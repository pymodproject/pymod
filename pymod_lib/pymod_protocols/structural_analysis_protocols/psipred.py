# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os

from pymod_lib import pymod_os_specific as pmos
from pymod_lib.pymod_vars import psipred_output_extensions
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.similarity_searches_protocols.psiblast import PSI_BLAST_common
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import catch_protocol_exception


class PSIPRED_prediction(PyMod_protocol, PSI_BLAST_common):
    """
    Class for a psipred prediction protocol.
    """

    protocol_name = "psipred"
    error_title = "PSIPRED error"

    def __init__(self, pymod, target_sequences=None):
        PyMod_protocol.__init__(self, pymod)
        self.target_sequences = self.get_pymod_elements(target_sequences)


    @catch_protocol_exception
    def launch_from_gui(self):
        if len(self.target_sequences) == 0:
            self.pymod.main_window.show_error_message("PSIPRED Error",
                "Please select at least one sequence to be analyzed with PSIPRED.")
            return False

        if not self.check_psipred_parameters():
            return False

        # Run psipred on all the selected sequences.
        self.predicted_sequences = []
        if not self.pymod.use_protocol_threads:
            self.predict_all_sequences()
        else:
            label_text = ("Running psipred on %s. Please wait for the process to"
                          " complete..." % self.get_seq_text(self.target_sequences))

            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.predict_all_sequences,
                                            args=(),
                                            title="Running psipred",
                                            label_text=label_text)
            p_dialog.exec_()

        # Colors the sequences by predicted secondary structure.
        for sequence in self.predicted_sequences:
            sequence.predicted_secondary_structure = True
            self.pymod.main_window.color_element_by_pred_sec_str(sequence, color_structure=True)


    def check_psipred_parameters(self): # predict_secondary_structure(self, elements=None):
        """
        Checks that the files needed to run PSIPRED exists on users' machines.
        """

        # First checks for PSIPRED installation.
        if not self.pymod.psipred.tool_dir_exists(parameter_name="exe_dir_path"):
            self.pymod.psipred.tool_dir_not_found()
            return False

        psipred_exe_dirpath = self.pymod.psipred["exe_dir_path"].get_value()
        for _exe_filename in ("chkparse", "psipred", "psipass2"):
            exe_filename = pmos.get_exe_file_name(_exe_filename)
            if not os.path.isfile(os.path.join(psipred_exe_dirpath, exe_filename)):
                message = ("A PSIPRED executable file ('%s') does not exists in"
                           " the PSIPRED executables directory specified in the PyMod"
                           " Options Window ('%s'). If you want to use PSIPRED, please"
                           " specify in the Options Window a PSIPRED executables directory"
                           " where a '%s' file is found." % (exe_filename, psipred_exe_dirpath, exe_filename))
                self.pymod.main_window.show_error_message(self.error_title, message)
                return False

        # Then checks for PSIPRED datafiles.
        if not self.pymod.psipred.tool_dir_exists(parameter_name="data_dir_path"):
            message = ("PSIPRED 'data' directory not found! Please specify an existent"
                       " directory in the PSIPRED Options Window of PyMod.")
            self.pymod.main_window.show_error_message(self.error_title, message)
            return False
        psipred_data_dirpath = self.pymod.psipred["data_dir_path"].get_value()
        if len(os.listdir(psipred_data_dirpath)) == 0:
            message = ("PSIPRED 'data' directory is empty! Please specify a data"
                       " directory actually containing PSIPRED data file in the"
                       " in the Options Window of PyMod.")
            self.pymod.main_window.show_error_message(self.error_title, message)
            return False


        # Checks for PSI-BLAST on the user's system.
        if not self.pymod.blast_plus.tool_dir_exists(parameter_name="exe_dir_path"):
            self.pymod.blast_plus.tool_dir_not_found()
            return False

        psiblast_exe_dirpath = self.pymod.blast_plus["exe_dir_path"].get_value()
        if not os.path.isfile(os.path.join(psiblast_exe_dirpath, pmos.get_exe_file_name("psiblast"))):
            message = ("The PSI-BLAST executable file ('psiblast') does not exists in"
                       " the BLAST+ suite executables directory specified in the PyMod"
                       " Options Window ('%s'). If you want to use PSIPRED, please"
                       " specify in the Options Window a BLAST+ suite executables directory"
                       " where a 'psiblast' file is found." % (psiblast_exe_dirpath))
            self.pymod.main_window.show_error_message(self.error_title, message)
            return False

        # And finally checks for a BLAST database.
        if not self.pymod.psipred.tool_dir_exists(parameter_name="database_dir_path"):
            message = ("A directory containing a BLAST database was not found! Please"
                       " specify an existent directory in the PSIPRED options in the options window of PyMod.")
            self.pymod.main_window.show_error_message(self.error_title, message)
            return False

        dbpath = self.pymod.psipred["database_dir_path"].get_value()
        if not self.verify_valid_blast_dbdir(dbpath, remove_temp_files=True):
            message = ("The database '%s' directory does not contain a valid set"
                       " of database files." % (dbpath))
            self.pymod.main_window.show_error_message(self.error_title, message)
            return False

        return True


    def predict_all_sequences(self):
        for sequence in self.target_sequences:
            # Actually calls the method that launches PSIPRED.
            prediction_successful = self.run_psipred(sequence)
            if prediction_successful:
                self.predicted_sequences.append(sequence)


    def run_psipred(self, element):
        """
        Actually runs PSIPRED, collects its results and map them on the sequences in PyMod main
        window.
        """
        print_output = False
        sequence_header = element.my_header
        if print_output:
            print("- Beginning PSIPRED prediction for:", sequence_header)

        # The name of the BLAST database file.
        # If the database files are contained in a folder like this: /home/user/pymod/databases/swissprot/swissprot
        dbpath = self.pymod.psipred["database_dir_path"].get_value() # e.g.: /home/user/pymod/databases/swissprot
        if print_output:
            print("- dbpath:", dbpath)

        # Where the NCBI programs have been installed.
        ncbidir = self.pymod.blast_plus["exe_dir_path"].get_value()
        if print_output:
            print("- ncbidir:", ncbidir)

        # Where the PSIPRED V2 programs have been installed.
        execdir = self.pymod.psipred["exe_dir_path"].get_value()
        if print_output:
            print("- execdir:", execdir)

        # Where the PSIPRED V2 data files have been installed.
        datadir = self.pymod.psipred["data_dir_path"].get_value()
        if print_output:
            print("- datadir",datadir)

        # Write the temporary input fasta file, setting its basename.
        basename = "psipred_temp"
        if print_output:
            print("- basename: ", basename)
        self.pymod.build_sequence_file([element], basename, file_format="fasta", remove_indels=True, new_directory=self.pymod.psipred_dirpath)

        #---------------------
        # Execute PSI-BLAST. -
        #---------------------

        if print_output:
            print("- Running PSI-BLAST with sequence", basename ,"...")
        try:
            self.execute_psiblast(
                ncbi_dir = ncbidir,
                db_path = dbpath,
                query = os.path.join(self.pymod.psipred_dirpath, basename+".fasta"),
                inclusion_ethresh = 0.001,
                out_pssm = os.path.join(self.pymod.psipred_dirpath, basename+".chk"),
                out = os.path.join(self.pymod.psipred_dirpath, basename+".blast"),
                num_iterations = 3,
                num_alignments = 0)
            # psiblast_output = open("%s.blast" % os.path.join(self.pymod.psipred_dirpath, basename),"w")
            # self.pymod.execute_subprocess(psiblast_command, new_stdout=psiblast_output)
            # psiblast_output.close()

        except:
            if print_output:
                print("- FATAL: Error whilst running psiblast - script terminated!")
            raise Exception("There was an error while running PSI-BLAST, so PSIPRED cannot perform a prediction for %s." % (sequence_header))

        #--------------------
        # Execute chkparse. -
        #--------------------

        if print_output:
            print("- Predicting secondary structure...")
        chkdir_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("chkparse"))) + " " +
                          pmos.build_commandline_path_string("%s.chk" % os.path.join(self.pymod.psipred_dirpath, basename)))
        try:
            chkdir_output = open("%s.mtx" % os.path.join(self.pymod.psipred_dirpath, basename),"w")
            self.pymod.execute_subprocess(chkdir_command, new_stdout=chkdir_output)
            chkdir_output.close()
        except:
            if print_output:
                print("- FATAL: Error whilst running chkdir - script terminated!")
            raise Exception("No homologous sequences were found by PSI-BLAST for %s, so PSIPRED cannot perform a prediction for this sequence." % (sequence_header))

        #--------------------------
        # Execute PSIPRED pass 1. -
        #--------------------------

        psipass1_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipred"))) + " " +
                            pmos.build_commandline_path_string("%s.mtx" % os.path.join(self.pymod.psipred_dirpath, basename)) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat")) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat2")) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat3")))
        try:
            psipass1_output = open("%s.ss" % os.path.join(self.pymod.psipred_dirpath, basename),"w")
            self.pymod.execute_subprocess(psipass1_command, new_stdout=psipass1_output)
            psipass1_output.close()
        except Exception as e:
            if print_output:
                print("- FATAL: Error whilst running psipred 1 - script terminated!")
            raise Exception("There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))

        #--------------------------
        # Execute PSIPRED pass 2. -
        #--------------------------

        psipass2_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipass2"))) + " " +
                            "%s 1 1.0 1.0" % pmos.build_commandline_path_string(os.path.join(datadir,"weights_p2.dat")) + " " +
                            pmos.build_commandline_path_string("%s.ss2" % os.path.join(self.pymod.psipred_dirpath, basename)) + " " +
                            pmos.build_commandline_path_string("%s.ss" % os.path.join(self.pymod.psipred_dirpath, basename)))
        try:
            psipass2_output = open("%s.horiz" % os.path.join(self.pymod.psipred_dirpath, basename),"w")
            self.pymod.execute_subprocess(psipass2_command, new_stdout=psipass2_output)
            psipass2_output.close()
        except:
            if print_output:
                print("- FATAL: Error whilst running psipass 2 - script terminated!")
            raise Exception("There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))

        #--------------------------
        # Clean up PSIPRED files. -
        #--------------------------

        if print_output:
            print("- Cleaning up ...")

        # Remove temporary files.
        self.remove_psipred_temp_files()

        # Renames the output files.
        output_files_name = pmos.clean_file_name(element.my_header)
        for ext in psipred_output_extensions:
            os.rename(os.path.join(self.pymod.psipred_dirpath, basename+ext),
                      os.path.join(self.pymod.psipred_dirpath, output_files_name+ext))

        if print_output:
            print("- Final output files:" + output_files_name + ".ss2 " + output_files_name + ".horiz")
            print("- Finished.")

        #----------------------------------------------
        # Parses the results from .horiz output file. -
        #----------------------------------------------

        results_file = open(os.path.join(self.pymod.psipred_dirpath, output_files_name+".horiz"),"r")
        confs = "" # String for confidence scores of each residue.
        preds = "" # String for the secondary structure elements prediction of each residue.
        for l in results_file.readlines():
            if l.startswith("Conf:"):
                rl = l[6:66].replace(" ","").replace("\n","")
                confs += rl
            elif l.startswith("Pred:"):
                rl = l[6:66].replace(" ","").replace("\n","")
                preds += rl
        results_file.close()

        # Actually stores in the PyMod elements the results.
        element.psipred_elements_list = []
        for c, e, res in zip(confs, preds, element.get_polymer_residues()):
            res.psipred_result = {"confidence": int(c), "sec-str-element": e}

        return True


    def remove_psipred_temp_files(self):
        try:
            for f in os.listdir(self.pymod.psipred_dirpath):
                if not os.path.splitext(f)[1] in psipred_output_extensions:
                    os.remove(os.path.join(self.pymod.psipred_dirpath,f))
        except:
            pass

    def quit_protocol(self):
        self.remove_psipred_temp_files()
        PyMod_protocol.quit_protocol(self)
