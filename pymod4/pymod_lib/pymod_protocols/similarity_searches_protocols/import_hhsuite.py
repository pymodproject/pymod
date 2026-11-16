# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for parsing HH-suite (https://github.com/soedinglab/hh-suite) results files
and to import their content in PyMod.
"""

import os
import re
import shutil
import gzip
import urllib.request

from Bio import SeqIO

from pymol.Qt import QtWidgets, QtCore
from pymol import cmd

from pymod_lib.pymod_externals.hhsuite.hh_reader import read_result
from pymod_lib.pymod_externals.hhsuite.hhmakemodel import to_seq
from pymod_lib.pymod_externals.hhsuite.hhmakemodel import main as hhmakemodel_main
from pymod_lib.pymod_gui.hits_table_qt import coverage_span
from pymod_lib.pymod_gui.shared_gui_components_qt import (small_font_style, highlight_color,
                                                          PyMod_protocol_window_qt,
                                                          PyMod_hbox_option_qt,
                                                          askyesno_qt, askopenfile_qt)
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import (Generic_BLAST_search,
                                                                                 Similarity_searches_results_window_qt)
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import PyModInterruptedProtocol
from pymod_lib.pymod_protocols.structural_databases_protocols import Associate_structure
from pymod_lib.pymod_os_specific import check_network_connection


###############################################################################
# Classes to import HH-suite results in PyMod.                                #
###############################################################################

class Import_HHsuite_results(Generic_BLAST_search):

    def additional_initialization(self):
        pass


    def launch_from_gui(self, mode="hhr"):

        # Import template search results from HH-suite.
        if mode == "hhr":
            self.launch_from_gui_hhr()

        # Import a multiple sequence alignment from HH-suite.
        elif mode == "a3m":
            self.launch_from_gui_a3m()

        else:
            raise KeyError("Unknown 'mode': {}".format(mode))


    ###########################################################################
    # Import results from a HHR file.                                         #
    ###########################################################################

    import_error_message = "Import HH-suite Results Error"

    def launch_from_gui_hhr(self):
        """
        Import template search results from HH-suite.
        """

        # Select the file to open.
        self.hhr_filepath = askopenfile_qt("Select a HHR file to open",
                                           name_filter="*.hhr",
                                           parent=self.pymod.get_qt_parent())

        if not self.hhr_filepath:
            return None


        # Parses the hhr results file.
        try:
            hhr_results = read_result(self.hhr_filepath)
        except Exception as e:
            self.pymod.main_window.show_error_message(self.import_error_message,
                ("The HHR results file appears to be invalid and the following error"
                " was raised: {}".format(str(e))))
            return None

        # Check for empty results.
        if not hhr_results:
            self.pymod.main_window.show_warning_message(self.import_error_message,
                "Empty HHR file: this file does not contain any alignment.")
            return None

        # Store the results.
        self.hhr_results = hhr_results
        self.query_id = self.hhr_results[0].query_id
        if len(self.query_id) > 35:
            self.query_id = self.query_id[0:35] + "..."


        # Shows the results window.
        # 'AF' says whether the hits carry models of a prediction database, which
        # an HH-suite search does not: it is a search of templates. The window it
        # shares with the other searches asks for it, and not answering was a
        # TypeError which stopped this command before it opened anything.
        self.import_hhr_window = HHsuite_import_hhr_window_qt(parent=self.pymod.main_window,
                                                              protocol=self,
                                                              AF=self.AF)
        self.import_hhr_window.show()


    def get_subject_name(self, hsp, max_len=100):
        t = "{} ({})".format(hsp.template_id, hsp.template_info[1:])
        if len(t) > max_len:
            t = t[0:max_len] + "..."
        return t.replace("\n", "")

    def get_hsp_evalue(self, hsp):
        return hsp.evalue

    def get_prob_value(self, hsp):
        return hsp.probability

    def get_hsp_query_seq(self, hsp):
        return to_seq(hsp.query_ali)

    def get_hsp_subj_seq(self, hsp):
        return to_seq(hsp.template_ali)


    def blast_results_state(self):
        """
        Called when the "SUBMIT" button is pressed on the results window.
        """

        # For each hsp takes the state of its check-box.
        self.my_blast_map = [chk.isChecked() for chk in self.import_hhr_window.sbjct_checkbuttons_list]

        self.import_hhr_window.destroy()

        # If the user selected at least one template.
        if any(self.my_blast_map):

            # This will actually import the sequences inside PyMod.
            self.import_results_in_pymod()


    def import_results_in_pymod(self):
        """
        Builds a list containing those hits that were selected by the user in the BLAST results
        window.
        """

        if not check_network_connection("https://google.com", timeout=3):
            has_network = False
        else:
            has_network = True

        #------------------------------------------
        # Get the templates selected by the user. -
        #------------------------------------------

        self.imported_hsp_list = []

        for hsp_counter, hsp in enumerate(self.hhr_results):

            if self.my_blast_map[hsp_counter]:

                template_id = hsp.template_id

                hsp_dict = {"hsp": hsp, "pdb_code": None, "chain_id": None,
                            "hsp_num": str(hsp_counter + 1), "successful": False}

                if re.match("([a-zA-Z0-9]{4})_([A-Za-z])$", template_id):
                    pdb_code, chain_id = template_id.split("_")
                    hsp_dict["pdb_code"] = pdb_code
                    hsp_dict["chain_id"] = chain_id

                self.imported_hsp_list.append(hsp_dict)


        # Check if the CIF files of the hits can be fetched from the PDB.
        pdb_hsp_list = [h for h in self.imported_hsp_list if h["pdb_code"] is not None]

        if len(pdb_hsp_list) == len(self.imported_hsp_list):

            if has_network:
                message = ("Would you like to fetch from the Internet the 3D"
                           " structures of the templates that you selected?")
                fetch_cif_files = askyesno_qt("Fetch 3D structures?", message,
                                              parent=self.pymod.get_qt_parent())
            else:
                message = ("No network connection. The template 3D structures can"
                           " not be fetched now from the Internet.")
                self.pymod.main_window.show_warning_message(self.import_error_message,
                                                            message)
                fetch_cif_files = False

        elif len(pdb_hsp_list) == 0:

            fetch_cif_files = False

        else:

            if has_network:
                n_pdb_missing = len(self.imported_hsp_list)-len(pdb_hsp_list)
                message = ("You selected {} (out of {}) hit sequences which do not appear"
                           " to be valid templates from the PDB. No 3D structure"
                           " will be fetched now. You can fetch the PDB structures"
                           " of each hit sequence having a valid PDB id later at any"
                           " moment".format(n_pdb_missing, len(self.imported_hsp_list)))
                self.pymod.main_window.show_warning_message(self.import_error_message,
                                                            message)
            fetch_cif_files = False


        #--------------------------------------
        # Prepare the input and output files. -
        #--------------------------------------

        # Prepare the templates input CIF directory. This will be needed by the
        # 'hhmakemodel.py' script.
        self.tem_in_dirpath = os.path.join(self.output_directory, "hhsuite_tem_input")
        if os.path.isdir(self.tem_in_dirpath):
            shutil.rmtree(self.tem_in_dirpath)
        os.mkdir(self.tem_in_dirpath)

        # Prepare the templates output CIF directory. This will be needed by the
        # 'hhmakemodel.py' script.
        self.tem_out_dirpath = os.path.join(self.output_directory, "hhsuite_tem_output")
        if os.path.isdir(self.tem_out_dirpath):
            shutil.rmtree(self.tem_out_dirpath)
        os.mkdir(self.tem_out_dirpath)

        # Set the path of the ouput PIR file generated by the 'hhmakemodel.py' script.
        self.pir_out_filepath = os.path.join(self.output_directory, "hhsuite_alignment.pir")


        #------------------------------------
        # Actually downloads the CIF files. -
        #------------------------------------

        if fetch_cif_files:
            try:

                if not self.pymod.use_protocol_threads:
                    self.download_all_cif_files()

                else:
                    label_text = ("Connecting to the PDB to download %s. Please wait for"
                                  " completion..." % self.get_seq_text(self.imported_hsp_list, "CIF file"))
                    p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                                    function=self.download_all_cif_files,
                                                    args=(), title="Downloading CIF files",
                                                    wait_start=0.15, wait_end=0.15,
                                                    label_text=label_text)
                    p_dialog.exec_()

            except PyModInterruptedProtocol:
                return None


            # Check if there were some structures which could not be fetched.
            n_failures = len([h for h in self.imported_hsp_list if not h["successful"]])
            if n_failures != 0:

                title = "Download Error"
                message = ("Can not access the PDB database to download %s structures out of %s."
                           " These templates will not be imported in PyMod."
                           "\nPlease check your Internet connection or if the PDB ids of the"
                           " structures are valid." % (n_failures, len(self.imported_hsp_list)))
                self.pymod.main_window.show_warning_message(title, message)

                if n_failures == len(self.imported_hsp_list):
                    self.pymod.main_window.show_warning_message("Import HHsuite Results Error",
                        "No templates could be fetched. Quitting.")
                    return None

            self.selected_templates_nums = [h["hsp_num"] for h in self.imported_hsp_list if h["successful"]]

        else:

            self.selected_templates_nums = [h["hsp_num"] for h in self.imported_hsp_list]


        #------------------------------------
        # Runs the 'hhmakemodel.py' script. -
        #------------------------------------

        # Prepare the arguments for the hhsuite function.
        hhsuite_args = {"input": self.hhr_filepath,
                        "cifs": self.tem_in_dirpath,
                        "pir": self.pir_out_filepath,
                        "output": self.tem_out_dirpath,
                        "verbose": True,
                        "m": self.selected_templates_nums,
                        "e": None, "r": 0, # Do not use any filters.
                        "c": True}

        # Run the hhsuite function.
        hhmakemodel_main(hhsuite_args)


        #----------------------------------------------------------
        # Parse the output PIR file produced by the above script. -
        #----------------------------------------------------------

        elements_to_load = []

        ali_records = list(SeqIO.parse(self.pir_out_filepath, "pir"))
        for record_idx, record in enumerate(ali_records):
            if record_idx == 0:
                element_name = self.query_id
            else:
                if fetch_cif_files:
                    element_name = "__temp_template_{}__".format(record_idx-1)
                else:
                    element_name = self.imported_hsp_list[record_idx-1]["hsp"].template_id

            new_element = self.pymod.build_pymod_element_from_args(element_name,
                                                                   str(record.seq).replace("*", "-"))
            elements_to_load.append(new_element)
            self.pymod.add_element_to_pymod(new_element)

        # Add a new alignment object to PyMod in which contains the query sequence
        # and that of the selected templates.
        new_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                          child_elements=elements_to_load,
                                                          algorithm="imported")


        #-----------------------------------------------
        # Actually imports the final results in PyMod. -
        #-----------------------------------------------

        if fetch_cif_files:

            for record_idx, record in enumerate(ali_records):
                if record_idx == 0:
                    continue
                chain = record.description.split(":")[3]
                modified_cif_filepath = os.path.join(self.tem_in_dirpath, "{}.cif".format(record.id))
                modified_pdb_filepath = os.path.join(self.tem_in_dirpath, "{}.pdb".format(record.id))
                tem_temp_name = "_hhsuite_template_{}".format(record_idx)
                cmd.load(modified_cif_filepath, tem_temp_name)
                cmd.save(modified_pdb_filepath, tem_temp_name)
                cmd.delete(tem_temp_name)

                try:
                    a = Associate_structure(self.pymod, elements_to_load[record_idx])
                    a.associate(modified_pdb_filepath, chain)
                except Exception as e:
                    title = "Associate Structure Failure"
                    message = ("The structure association for %s chain %s failed because"
                               " of the following error: %s" % (record.id, chain, str(e)))
                    self.pymod.main_window.show_error_message(title, message)


        #-------------
        # Cleans up. -
        #-------------

        if os.path.isdir(self.tem_in_dirpath):
            shutil.rmtree(self.tem_in_dirpath)

        if os.path.isdir(self.tem_out_dirpath):
            shutil.rmtree(self.tem_out_dirpath)

        if os.path.isfile(self.pir_out_filepath):
            os.remove(self.pir_out_filepath)

        self.pymod.main_window.gridder(update_clusters=True, update_menus=True, update_elements=True)


    def download_all_cif_files(self):
        """
        This method will call other methods to tetch the CIF files needed to import
        the HHsuite templates in PyMod.
        """
        for hsp_dict in self.imported_hsp_list:
            self._fetch_single_element(hsp_dict)


    def _fetch_single_element(self, hsp_dict):
        """
        Download the CIF file for a single element.
        """

        pdb_code = hsp_dict["pdb_code"]
        output_name = "{}.cif".format(pdb_code)

        try:
            pdb_url = "https://files.rcsb.org/download/{}.cif.gz".format(pdb_code)
            temp_gzip_file_name = urllib.request.urlretrieve(pdb_url)[0]
            open_gzip_file = gzip.open(temp_gzip_file_name) # Uncompress the file while reading
            output_path = os.path.join(self.tem_in_dirpath, output_name)
            saved_file = open(output_path, 'wb')
            saved_file.write(open_gzip_file.read()) # Write pdb file
            open_gzip_file.close()
            saved_file.close()

            hsp_dict["successful"] = True

        except:

            pass


    ###########################################################################
    # Import a multiple sequence alignment from a A3M file built by HH-suite. #
    ###########################################################################

    def launch_from_gui_a3m(self):

        # Select the file to open.
        self.a3m_filepath = askopenfile_qt("Select a A3M file to open",
                                           name_filter="*.a3m *.oa3m",
                                           parent=self.pymod.get_qt_parent())

        if not self.a3m_filepath:
            return None


        # Parses the a3m file.
        error_title = "Import HH-suite Results Error"
        try:
            a3m_content = []
            with open(self.a3m_filepath) as a3m_handle:
                records = list(SeqIO.parse(a3m_handle, "fasta"))
            for record in records:
                a3m_content.append((record.id, str(record.seq)))

        except Exception as e:
            self.pymod.main_window.show_error_message(self.import_error_message,
                ("The A3M file appears to be invalid and the following error"
                " was raised: {}".format(str(e))))
            return None

        # Check for empty results.
        if not a3m_content:
            self.pymod.main_window.show_warning_message(self.import_error_message,
                "Empty A3M file: this file does not contain any aligned sequence.")
            return None

        # Store the results.
        self.a3m_content = a3m_content
        self.query_id = self.a3m_content[0][0]
        if len(self.query_id) > 35:
            self.query_id = self.query_id[0:35] + "..."

        # Shows the options window.
        self.import_a3m_window = HHsuite_import_a3m_window_qt(parent=self.pymod.main_window,
                                                              protocol=self,
                                                              title="A3M import options",
                                                              upper_frame_title="Otions for importing an A3M file content",
                                                              submit_command=self.a3m_import_state)
        self.import_a3m_window.show()


    def a3m_import_state(self):

        # Gets the number of sequences to import.
        try:
            n_seqs_to_import = int(self.import_a3m_window.seqs_to_import_spinbox.value())
        except Exception as e:
            self.pymod.main_window.show_error_message("Input Error",
                                                      "Non valid numeric input. Please correct.")
            return None

        # Convert the sequences in the A3M format in FASTA format.
        fasta_content = convert_items_a3m_to_fasta(self.a3m_content[0:n_seqs_to_import+1])

        # Import the results in PyMod.
        elements_to_load = []
        for item_idx, item in enumerate(fasta_content):
            new_element = self.pymod.build_pymod_element_from_args(item[0], item[1])
            elements_to_load.append(new_element)
            self.pymod.add_element_to_pymod(new_element)
        new_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                          child_elements=elements_to_load,
                                                          algorithm="imported")

        self.import_a3m_window.destroy()

        self.pymod.main_window.gridder(update_clusters=True, update_menus=True, update_elements=True)


def convert_items_a3m_to_fasta(a3m_items):
    """
    Receives a 'items' list. Each element in the list is a two-element tuple
    containg the id of a sequence and its aligned sequence in the A3M format
    (see: https://github.com/soedinglab/hh-suite/wiki#generating-a-multiple-sequence-alignment-using-hhblits).
    Returns a similar list, in which the sequences are in the FASTA format.
    """

    full_alignments = [list(item[1]) for item in a3m_items]

    for i, seq_i in enumerate(full_alignments):

        for pos_idx, pos in enumerate(seq_i):

            if pos.islower():

                for j, seq_j in enumerate(full_alignments):

                    if i != j:
                        seq_j.insert(pos_idx, "-")

    fasta_items = []
    for item_idx, a3m_item in enumerate(a3m_items):
        fasta_item = (a3m_item[0], "".join(full_alignments[item_idx]).upper())
        fasta_items.append(fasta_item)

    return fasta_items


###############################################################################
# GUI.                                                                        #
###############################################################################

class HHsuite_import_hhr_window_qt(Similarity_searches_results_window_qt):
    """
    Window for showing similarity searches results.
    """

    def _get_window_title(self):
        return "Import HH-suite Results"


    def _get_upper_frame_title(self):
        title_text = ("HHsearch Output for: %s\nFound %s sequences\nPlease Select the Sequences to"
                      " Import" % (self.protocol.query_id,
                                   len(self.protocol.hhr_results)))
        return title_text


    def result_columns(self):
        """
        The columns of an HH-suite result. There is no model database behind it
        and the identity is always reported, so the columns are fixed; the hit
        is a template, which is what its span is called.
        """

        return (("", 0), ("Name", 420), ("E-Value", 0), ("Probability", 0),
                ("Identity", 0), ("Query span", 150), ("Template span", 150))


    def display_blast_hits(self):
        """
        One row per hit of the .hhr file, in the table which every search of
        PyMod shows its hits in.
        """

        self.sbjct_checkbuttons_list = []

        for hsp in self.protocol.hhr_results:

            subject_name = self.protocol.get_subject_name(hsp)

            matches_count, identities_count = self.protocol.get_hsp_matches(hsp)
            seqid = identities_count/matches_count

            query_span_val = (hsp.end[0]-hsp.start[0])/hsp.query_length
            template_span_val = (hsp.end[1]-hsp.start[1])/hsp.template_length

            values = [subject_name,
                      "%.2e" % (self.protocol.get_hsp_evalue(hsp)),
                      str(self.protocol.get_prob_value(hsp)),
                      "{}/{} ({:.1f}%)".format(identities_count, matches_count,
                                               seqid*100),
                      "{}-{} ({:.1f}%)".format(hsp.start[0], hsp.end[0],
                                               query_span_val*100),
                      "{}-{} ({:.1f}%)".format(hsp.start[1], hsp.end[1],
                                               template_span_val*100)]

            # Both lengths are known here, so both spans are drawn.
            spans = {"Query span": coverage_span(hsp.query_length,
                                                 hsp.start[0], hsp.end[0]),
                     "Template span": coverage_span(hsp.template_length,
                                                    hsp.start[1], hsp.end[1])}

            self.hits_table.add_row(values, payload=hsp, spans=spans,
                                    detail=subject_name, searchable=subject_name)
            self.sbjct_checkbuttons_list.append(self.hits_table.rows[-1][0])

        self.hits_table.finish()


class HHsuite_import_a3m_window_qt(PyMod_protocol_window_qt):
    """
    Window with options for importing multiple sequence alignments from a3m files.
    """

    default_offset = 10

    def build_protocol_middle_frame(self):

        n_tot_sequences = len(self.protocol.a3m_content) - 1

        # Sub-frame created to select the number of sequences to import.
        range_label = "Import the first n. sequences (tot: {})".format(n_tot_sequences)
        self.range_subframe = PyMod_hbox_option_qt(label_text=range_label)

        # Spinbox to select the number of sequences to import.
        self.seqs_to_import_spinbox = QtWidgets.QSpinBox()
        self.seqs_to_import_spinbox.setRange(1, n_tot_sequences)
        self.seqs_to_import_spinbox.setValue(min((100, n_tot_sequences)))
        self.range_subframe.hbox.addWidget(self.seqs_to_import_spinbox)
        self.range_subframe.set_auto_input_widget_width()
        self.middle_formlayout.add_widget_to_align(self.range_subframe)
