# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Commands executed when interacting with the PyMod main window.
"""

import os
import sys
import shutil
import re
import webbrowser
import json

from Bio import Phylo

from pymol.Qt import QtWidgets

import pymol

import pymod_lib
from pymod_lib import pymod_os_specific as pmos
from pymod_lib import pymod_vars
from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt, asksaveasfile_qt, askopenfile_qt, askopenfiles_qt
from pymod_lib.pymod_gui.specific_gui_components_qt import Raw_sequence_window_qt, PyMod_options_window_qt

from pymod_lib.pymod_plot import pymod_plot_qt
from pymod_lib.pymod_seq.seq_manipulation import compute_sequence_identity
from pymod_lib.pymod_exceptions import PyModInvalidFile, PyModUnknownFile
from pymod_lib.pymod_gui.pymod_table import QtTableWindow

from pymod_lib.pymod_protocols.similarity_searches_protocols.ncbi_blast import NCBI_BLAST_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.local_blast import LOC_BLAST_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.psiblast import PSI_BLAST_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.phmmer import PHMMER_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.jackhmmer import Jackhmmer_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.hmmsearch import Hmmsearch_search

from pymod_lib.pymod_protocols.alignment_protocols.clustalw import Clustalw_regular_alignment, Clustalw_profile_alignment
from pymod_lib.pymod_protocols.alignment_protocols.clustalo import Clustalomega_regular_alignment, Clustalomega_profile_alignment
from pymod_lib.pymod_protocols.alignment_protocols.muscle import MUSCLE_regular_alignment
from pymod_lib.pymod_protocols.alignment_protocols.salign_seq import SALIGN_seq_regular_alignment, SALIGN_seq_profile_alignment

from pymod_lib.pymod_protocols.alignment_protocols.ce_alignment import CEalign_regular_alignment
from pymod_lib.pymod_protocols.alignment_protocols.salign_str import SALIGN_str_regular_alignment

from pymod_lib.pymod_protocols.domain_analysis_protocols.domain_analysis import Domain_Analysis_Protocol, Fuse_domains_protocol
from pymod_lib.pymod_protocols.domain_analysis_protocols.split_domains import Split_into_domains_protocol

from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.campo import CAMPO_analysis
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.scr_find import SCR_FIND_analysis
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.entropy_scorer import Entropy_analysis
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.pair_conservation import Pair_conservation_analysis
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.weblogo import WebLogo_analysis
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.espript import ESPript_analysis
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.tree_building import Tree_building

from pymod_lib.pymod_protocols.structural_analysis_protocols.ramachandran_plot import Ramachandran_plot
from pymod_lib.pymod_protocols.structural_analysis_protocols.contact_map_analysis import Contact_map_analysis
from pymod_lib.pymod_protocols.structural_analysis_protocols.structural_divergence_plot import Structural_divergence_plot
from pymod_lib.pymod_protocols.structural_analysis_protocols.dope_assessment import DOPE_assessment, show_dope_plot
from pymod_lib.pymod_protocols.structural_analysis_protocols.secondary_structure_assignment import Secondary_structure_assignment
from pymod_lib.pymod_protocols.structural_analysis_protocols.psipred import PSIPRED_prediction

from pymod_lib.pymod_protocols.modeling_protocols.homology_modeling import MODELLER_homology_modeling
from pymod_lib.pymod_protocols.modeling_protocols.loop_refinement import MODELLER_loop_refinement

from pymod_lib.pymod_protocols.updater_protocols import UpdaterProtocol


class PyMod_main_menu_commands:

    ###############################################################################################
    # FILES MENU COMMANDS.                                                                        #
    ###############################################################################################

    def open_file_from_the_main_menu(self):
        """
        This method is called when using the 'File -> Open from File...' command in PyMod main menu.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        file_paths = askopenfiles_qt("Select files to open",
                                     name_filter="*.fasta *.fa *.gp *.pdb *.ent",
                                     parent=self.get_qt_parent())

        # Loads each files in PyMod.
        loaded_sequences = False
        for single_file_path in file_paths:
            extension = os.path.splitext(single_file_path)[1].replace(".","").lower()
            try:
                if extension in ("fasta", "fa"):
                    self.open_sequence_file(single_file_path, "fasta")
                elif extension == "gp":
                    self.open_sequence_file(single_file_path, "genbank")
                elif extension in ("pdb", "ent"):
                    self.open_structure_file(single_file_path, extension)
                else:
                    raise PyModUnknownFile()
                loaded_sequences = True

            except PyModInvalidFile:
                title = "File Type Error"
                message = "The selected File ('%s') is not a valid '%s' file." % (single_file_path, pymod_vars.supported_sequence_file_types[extension])
                self.main_window.show_error_message(title, message)
            except PyModUnknownFile:
                title = "File Type Error"
                message = "The selected File ('%s') has an unknown file extension '%s'." % (single_file_path, extension)
                self.main_window.show_error_message(title, message)

        if loaded_sequences:
            self.main_window.gridder()


    def open_alignment_from_main_menu(self):
        """
        Lets users import in Pymod an alignment stored in an external file.
        """
        openfilename, extension = self.choose_alignment_file()
        if not None in (openfilename, extension):
            self.build_cluster_from_alignment_file(openfilename, extension)
        self.main_window.gridder(update_menus=True, update_elements=True)


    #################################################################
    # Add new sequences.                                            #
    #################################################################

    def show_raw_seq_input_window(self):
        """
        Launched when the user wants to add a new sequence by directly typing it into a Text entry.
        """
        self.raw_seq_window = Raw_sequence_window_qt(self.main_window,
                                                     title="Add Raw Sequence",
                                                     upper_frame_title="Type or Paste your Sequence",
                                                     submit_command=self.raw_seq_input_window_state)
        self.raw_seq_window.show()
        # self.raw_seq_window.resize(700, self.raw_seq_window.sizeHint().height())


    def raw_seq_input_window_state(self):
        """
        This is called when the SUBMIT button of the 'Add Raw Sequence' is pressed.
        """
        def special_match(strg, search=re.compile(r'[^QWERTYIPASDFGHKLXCVNM-]').search):
            return not bool(search(strg))
        def name_match(strg, search2=re.compile(r'[^a-zA-Z0-9_]').search):
            return not bool(search2(strg))

        sequence = self.raw_seq_window.get_sequence()
        sequence_name = self.raw_seq_window.get_sequence_name()

        if special_match(sequence) and len(sequence):
            if len(sequence_name) and name_match(sequence_name):
                self.add_element_to_pymod(self.build_pymod_element_from_args(sequence_name, sequence))
                self.raw_seq_window.destroy()
                self.main_window.gridder()
            else:
                title = 'Sequence Name Error'
                message = 'Please check the sequence name: only letters, numbers and "_" are allowed.'
                self.main_window.show_error_message(title, message)
        else:
            title = 'Sequence Error'
            message = 'Please check your sequence: only standard amino acid letters, "X" and "-" are allowed.'
            self.main_window.show_error_message(title, message)


    #################################################################
    # Saving files.                                                 #
    #################################################################

    def sequence_save_dialog(self, element):
        """
        Save a single sequence to a file.
        """
        # Ask to remove indels.
        if "-" in element.my_sequence:
            remove_indels_choice = self.ask_to_remove_indels()
        else:
            remove_indels_choice = False

        # Choose the file path.
        filepath = asksaveasfile_qt("Save FASTA file", name_filter="*.fasta", parent=self.get_qt_parent())
        if not filepath:
            return None

        # Actually saves the file. The file extension will not be added automatically.
        dirpath = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        self.build_sequence_file([element], filename, file_format="fasta", remove_indels=remove_indels_choice,
                                 new_directory=dirpath, add_extension=False)


    def save_selection_dialog(self, mode="selection"):
        """
        Save selection in a single file.
        """
        # Builds the selection.
        if mode == "selection":
            selection = self.get_selected_sequences()
        elif mode == "all":
            selection = self.get_all_sequences()
        else:
            raise KeyError("Unknown 'mode': %s" % mode)

        # Ask users if they want to include indels in the sequences to save.
        remove_indels_choice = False
        for e in selection:
            if "-" in e.my_sequence:
                remove_indels_choice = self.ask_to_remove_indels()
                break

        # Ask users to choose a directory where to save the file.
        filepath = asksaveasfile_qt("Save FASTA file", name_filter="*.fasta", parent=self.get_qt_parent())
        if not filepath:
            return None

        dirpath = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        self.build_sequence_file(selection, filename, file_format="fasta",
                                 remove_indels=remove_indels_choice, same_length=remove_indels_choice,
                                 new_directory=dirpath, add_extension=False)


    def ask_to_remove_indels(self):
        title = "Save File"
        message = "Would you like to remove indels from the sequence when saving it to a file?"
        remove_indels_choice = askyesno_qt(title, message, parent=self.get_qt_parent())
        return remove_indels_choice


    def alignment_save_dialog(self, alignment_element):
        """
        Lets the user choose the path to which an alignment file is going to be saved, and saves
        an alignment file there.
        """
        save_file_full_path = asksaveasfile_qt("Save an alignment file",
                                               name_filter="*.fasta;;*.aln;;*.sto",
                                               parent=self.get_qt_parent())
        if not save_file_full_path:
            return None

        alignment_file_name, extension = os.path.splitext(os.path.basename(save_file_full_path))
        extension = extension.replace(".", "")

        # The get all the aligned elements.
        aligned_elements = alignment_element.get_children()

        # Saves a file with all the sequences in the project "Alignments" directory.
        if extension == "fasta":
            self.save_alignment_fasta_file(alignment_file_name, aligned_elements)
        elif extension == "aln":
            self.build_sequence_file(aligned_elements, alignment_file_name, file_format="clustal", remove_indels=False)
        elif extension == "sto":
            self.build_sequence_file(aligned_elements, alignment_file_name, file_format="stockholm", remove_indels=False)
        # If the user didn't provide a valid extension.
        else:
            title = "Format Error"
            if extension != "":
                message = "Unknown alignment file extension: '%s'." % (extension)
            else:
                message = "No alignment file extension provided."
            message += " Please provide a valid extension. Example: .fasta (FASTA), .aln (Clustal) or .sto (Stockholm)"
            self.main_window.show_error_message(title, message)
            return

        # Moves the saved file to the path chosen by the user.
        try:
            old_path = os.path.join(self.alignments_dirpath, alignment_file_name + "." + extension)
            os.rename(old_path, save_file_full_path)
        except Exception as e:
            title = "File Error"
            message = "Could not save the alignment file to path '%s' for the following reason: %s." % (save_file_full_path, e)
            self.main_window.show_error_message(title, message)


    def save_pdb_chain_to_file_dialog(self, pymod_element):
        """
        Save a PDB single chain to a file.
        """

        # Choose the file path.
        filepath = asksaveasfile_qt("Save PDB file for this chain", name_filter="*.pdb", parent=self.get_qt_parent())

        if not filepath:
            return None

        # Actually saves the file.
        try:
            if os.path.isfile(filepath):
                os.remove(filepath)
            shutil.copy(pymod_element.structure.current_chain_file_path, filepath)
        except Exception as e:
            title = "File Error"
            message = "Could not save the PDB chain file to path '%s' for the following reason: %s." % (filepath, e)
            self.main_window.show_error_message(title, message)


    ###############################################################################################
    # SIMILARITY SEARCHES.                                                                        #
    ###############################################################################################

    def launch_blast_algorithm(self, blast_version):
        """
        Called when BLAST or PSI-BLAST is launched from the main menu.
        """
        if blast_version == "blast":
            blast_search = NCBI_BLAST_search(self, output_directory=self.similarity_searches_dirpath)
        elif blast_version == "psi-blast":
            blast_search = PSI_BLAST_search(self, output_directory=self.similarity_searches_dirpath)
        elif blast_version == "blastp":
            blast_search = LOC_BLAST_search(self, output_directory=self.similarity_searches_dirpath)
        else:
            raise KeyError(blast_version)

        blast_search.launch_from_gui()


    def launch_hmmer_algorithm(self, hmmer_version):
        if hmmer_version == "phmmer":
            hmmer_search = PHMMER_search(self, output_directory=self.similarity_searches_dirpath)
        elif hmmer_version == "jackhmmer":
            hmmer_search = Jackhmmer_search(self, output_directory=self.similarity_searches_dirpath)
        elif hmmer_version == "hmmsearch":
            hmmer_search = Hmmsearch_search(self, output_directory=self.similarity_searches_dirpath)
        else:
            raise KeyError(hmmer_version)
        hmmer_search.launch_from_gui()


    ###############################################################################################
    # DOMAIN ANALYSIS                                                                             #
    ###############################################################################################

    def launch_domain_analysis(self, mode):
        dom_an = Domain_Analysis_Protocol(self, mode)
        dom_an.launch_from_gui()


    ###############################################################################################
    # ALIGNMENT BUILDING.                                                                         #
    ###############################################################################################

    def launch_alignment_from_the_main_menu(self, program, strategy):
        """
        Launched from the 'Sequence', 'Structure Alignment' or 'Profile Alignment' from the submenus
        of the main window.
        """

        # Regular.
        if strategy == "regular":

            # Sequence alignments.
            if program == "clustalw":
                aligment_protocol_class = Clustalw_regular_alignment
            elif program == "clustalo":
                aligment_protocol_class = Clustalomega_regular_alignment
            elif program == "muscle":
                aligment_protocol_class = MUSCLE_regular_alignment
            elif program == "salign-seq":
                aligment_protocol_class = SALIGN_seq_regular_alignment

            # Structural alignments.
            elif program == "ce":
                aligment_protocol_class = CEalign_regular_alignment
            elif program == "salign-str":
                aligment_protocol_class = SALIGN_str_regular_alignment

        # Profile.
        elif strategy == "profile":

            if program == "clustalw":
                aligment_protocol_class = Clustalw_profile_alignment
            elif program == "clustalo":
                aligment_protocol_class = Clustalomega_profile_alignment
            elif program == "salign-seq":
                aligment_protocol_class = SALIGN_seq_profile_alignment

        # Actually launches the alignment protocol.
        a = aligment_protocol_class(self, output_directory=self.alignments_dirpath)
        a.launch_from_gui()


    ###############################################################################################
    # STRUCTURAL ANALYSIS TOOLS.                                                                  #
    ###############################################################################################

    def assign_secondary_structure(self, element):
        sec_str_assignment = Secondary_structure_assignment(self, element)
        sec_str_assignment.assign_secondary_structure()


    def dope_from_main_menu(self):
        """
        Called when users decide calculate DOPE of a structure loaded in PyMod.
        """
        dope_assessment = DOPE_assessment(self)
        dope_assessment.launch_from_gui()


    def ramachandran_plot_from_main_menu(self):
        """
        PROCHEK style Ramachandran Plot.
        """
        ramachandran_plot = Ramachandran_plot(self)
        ramachandran_plot.launch_from_gui()


    def contact_map_from_main_menu(self):
        """
        Contact/distance map analysis for one or more sequences.
        """
        contact_analysis = Contact_map_analysis(self)
        contact_analysis.launch_from_gui()


    def sda_from_main_menu(self):
        """
        Analyze the structural divergence between two or more aligned structures.
        """
        sda_analysis = Structural_divergence_plot(self)
        sda_analysis.launch_from_gui()


    def launch_psipred_from_main_menu(self):
        """
        Called when users decide to predict the secondary structure of a sequence using PSIPRED.
        """
        psipred_protocol = PSIPRED_prediction(self)
        psipred_protocol.launch_from_gui()


    # def superpose_from_main_menu(self):
    #     superpose_protocol = Superpose(self)
    #     superpose_protocol.launch_from_gui()


    ###############################################################################################
    # MODELING.                                                                                   #
    ###############################################################################################

    def launch_modeller_hm_from_main_menu(self):
        modeller_session = MODELLER_homology_modeling(self)
        modeller_session.launch_from_gui()


    def launch_modeller_lr_from_main_menu(self):
        modeller_session = MODELLER_loop_refinement(self)
        modeller_session.launch_from_gui()


    ###############################################################################################
    # PYMOD OPTIONS WINDOW.                                                                       #
    ###############################################################################################

    def show_pymod_options_window(self):
        """
        Builds a window that allows to modify some PyMod options.
        """

        self.pymod_options_window = PyMod_options_window_qt(self.main_window,
            pymod=self,
            title="PyMod Options",
            upper_frame_title="Here you can modify options for PyMod",
            submit_command=self.set_pymod_options_state)
        self.pymod_options_window.show()
        # self.pymod_options_window.resize(1100, 800)


    def set_pymod_options_state(self):
        """
        This function is called when the SUBMIT button is pressed in the PyMod options window.
        """
        old_projects_dir = self.pymod_plugin["pymod_dir_path"].get_value()
        new_projects_dir = self.pymod_options_window.get_value_from_gui(self.pymod_plugin["pymod_dir_path"])
        if not os.path.isdir(new_projects_dir):
            title = "Configuration Error"
            message = "The PyMod Projects Directory you specified ('%s') does not exist on your system. Please choose an existing directory." % (new_projects_dir)
            self.main_window.show_error_message(title, message)
            return False

        # Saves the changes to PyMod configuration file.
        with open(self.cfg_file_path, 'w') as  cfgfile:
            pymod_config_data = {}
            for tool in self.pymod_tools:
                new_tool_parameters = {}
                for parameter in tool.parameters:
                    if parameter.can_be_updated_from_gui():
                        new_tool_parameters.update({parameter.name: self.pymod_options_window.get_value_from_gui(parameter)})
                    else:
                        new_tool_parameters.update({parameter.name: parameter.get_value()})
                new_tool_dict = {tool.name: new_tool_parameters}
                pymod_config_data.update(new_tool_dict)
            cfgfile.write(json.dumps(pymod_config_data))

        # Then updates the values of the parameters of the tools contained in "self.pymod_tools" so
        # that they can be used in the current PyMod session.
        try:
            # Prevents the user from changing the project directory during a session.
            self.get_parameters_from_configuration_file()
            if old_projects_dir != new_projects_dir:
                title = "Configuration Updated"
                message = "You changed PyMod projects directory, the new directory will be used the next time you launch PyMod."
                self.main_window.show_warning_message(title, message)
            self.pymod_options_window.destroy()

        except Exception as e:
            self.show_configuration_file_error(e, "read")
            self.main_window.close()


    ###############################################################################################
    # ALIGNMENT MENU AND ITS BEHAVIOUR.                                                           #
    ###############################################################################################

    def save_alignment_to_file_from_ali_menu(self, alignment_element):
        self.alignment_save_dialog(alignment_element)


    def launch_campo_from_main_menu(self, pymod_cluster):
        campo = CAMPO_analysis(self, pymod_cluster)
        campo.launch_from_gui()

    def launch_scr_find_from_main_menu(self, pymod_cluster):
        scr_find = SCR_FIND_analysis(self, pymod_cluster)
        scr_find.launch_from_gui()

    def launch_entropy_scorer_from_main_menu(self, pymod_cluster):
        entropy_scorer = Entropy_analysis(self, pymod_cluster)
        entropy_scorer.launch_from_gui()

    def launch_pc_from_main_menu(self, pymod_cluster):
        entropy_scorer = Pair_conservation_analysis(self, pymod_cluster)
        entropy_scorer.launch_from_gui()

    def launch_weblogo_from_main_menu(self, pymod_cluster):
        weblogo = WebLogo_analysis(self, pymod_cluster)
        weblogo.launch_from_gui()

    def launch_espript_from_main_menu(self, pymod_cluster):
        espript = ESPript_analysis(self, pymod_cluster)
        espript.launch_from_gui()


    #################################################################
    # Build and display sequence identity and RMSD matrices of      #
    # alignments.                                                   #
    #################################################################

    def display_identity_matrix(self, alignment_element):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        # Then get all its children (the aligned elements).
        aligned_elements = alignment_element.get_children()
        n = len(aligned_elements)

        # identity_matrix = [[None]*n]*n # [] # Builds an empty (nxn) "matrix".
        identity_matrix = []
        for a in range(n):
            identity_matrix.append([None]*n)

        # Computes the identities (or anything else) and builds the matrix.
        for i in range(len(aligned_elements)):
            for j in range(len(aligned_elements)):
                if j >= i:
                    sid = compute_sequence_identity(aligned_elements[i].my_sequence, aligned_elements[j].my_sequence)
                    # This will fill "half" of the matrix.
                    identity_matrix[i][j] = sid
                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    identity_matrix[j][i] = sid

        # Build the list of sequences names.
        sequences_names = []
        for e in aligned_elements:
            sequences_names.append(e.compact_header)

        title = 'Identity matrix for ' + alignment_element.my_header
        self.show_table(sequences_names, sequences_names, identity_matrix, title)


    def display_rmsd_matrix_from_alignments_menu(self, alignment_element):
        # Checks the elements of the alignment.
        alignments_elements = alignment_element.get_children()
        aligned_structures = [e for e in alignments_elements if e.has_structure()]
        if len(aligned_structures) < 2:
            self.main_window.show_error_message("Selection Error",
                ("A RMSD matrix can only be computed on alignments containing at least"
                 " two elements having a 3D structure loaded in PyMOL."))
            return None
        if any([e.polymer_type == "nucleic_acid" for e in aligned_structures]):
            self.main_window.show_error_message("Selection Error",
                "Can not perform the analysis for nucleic acids structures.")
            return None

        if not alignment_element.algorithm in pymod_vars.structural_alignment_tools:
            message = pymod_vars.structural_alignment_warning % "RMSD matrix"
            self.main_window.show_warning_message("Alignment Warning", message)

        ali_protocol = CEalign_regular_alignment(self)
        alignment_element.rmsd_dict = ali_protocol.compute_rmsd_dict(aligned_structures)
        self.display_rmsd_matrix(alignment_element)


    def display_rmsd_matrix(self, alignment_element):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        aligned_elements = [e for e in alignment_element.get_children() if e.has_structure()]
        rmsd_dict = alignment_element.rmsd_dict
        rmsd_matrix_to_display = []
        n = len(aligned_elements)
        rmsd_matrix_to_display = [[None]*n for a in range(n)]

        for i, ei in enumerate(aligned_elements):
            for j, ej in enumerate(aligned_elements):
                if j >= i:
                    # This will fill "half" of the matrix.
                    if rmsd_dict[(ei.unique_index, ej.unique_index)] is not None:
                        rmsd_matrix_to_display[i][j] = round(rmsd_dict[(ei.unique_index, ej.unique_index)], 4)

                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    if rmsd_dict[(ej.unique_index, ei.unique_index)] is not None:
                        rmsd_matrix_to_display[j][i] = round(rmsd_dict[(ej.unique_index, ei.unique_index)], 4)

        # Build the list of sequences names.
        sequences_names = [e.compact_header for e in aligned_elements]
        title = 'RMSD matrix for ' + alignment_element.my_header
        self.show_table(sequences_names, sequences_names, rmsd_matrix_to_display, title)


    def show_table(self, column_headers=None, row_headers=None, data_array=[], title="New Table",
                   # columns_title=None, rows_title=None,
                   # rowheader_width=20, number_of_tabs=2,
                   sortable=False, row_headers_height=25,
                   width=800, height=450
                   ):
        """
        Displayes in a new window a table with data from the bidimensional 'data_array' numpy array.
        """

        # Builds a new window in which the table will be displayed.
        new_window = QtTableWindow(parent=self.main_window, title=title, data=data_array,
                                   row_labels=row_headers, row_labels_height=row_headers_height,
                                   column_labels=column_headers,
                                   sortable=sortable,)
        new_window.show()


    #################################################################
    # Show guide trees and build trees out of alignments.           #
    #################################################################

    def show_guide_tree_from_alignments_menu(self, alignment_element):
        """
        Shows the guide tree that was constructed in order to perform a multiple alignment.
        """
        # Gets the path of the .dnd file of the alignment.
        self.show_tree(alignment_element.get_tree_file_path())


    def show_tree(self, tree_file_path):
        # Reads a tree file using Phylo.
        tree = Phylo.read(tree_file_path, "newick")
        tree.ladderize() # Flip branches so deeper clades are displayed at top
        # Displays its content using pyqtgraph.
        pymod_plot_qt.draw_tree(tree=tree, pymod=self)


    def show_dendrogram_from_alignments_menu(self, alignment_element):
        """
        Shows dendrograms built by SALIGN.
        """
        pymod_plot_qt.draw_modeller_dendrogram(dendrogram_filepath=alignment_element.get_tree_file_path(),
                                               pymod=self)


    def build_tree_from_alignments_menu(self, alignment_element):
        """
        Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        menu.
        """
        tree_building = Tree_building(self, input_cluster_element=alignment_element)
        tree_building.launch_from_gui()


    ###############################################################################################
    # MODELS MENU AND ITS BEHAVIOUR.                                                              #
    ###############################################################################################

    def save_modeling_session(self, modeling_session):
        """
        Build a zip file of the modeling directory of a certain session.
        """
        archive_path = asksaveasfile_qt("Save PyMod Session file", name_filter="*.zip", parent=self.get_qt_parent())
        if not archive_path:
            return None

        try:
            pmos.zip_directory(directory_path=os.path.join(self.models_dirpath, os.path.basename(modeling_session.modeling_directory_path)),
                               zipfile_path=archive_path)
        except:
            title = "File Error"
            message = "Could not save the modeling session file to path: %s" % (archive_path)
            self.main_window.show_error_message(title, message)

    def show_session_profile(self, modeling_session):
        """
        Shows a DOPE profile of a modeling session.
        """
        show_dope_plot(dope_plot_data=modeling_session.dope_profile_data,
                       parent_window=self.main_window, pymod=self)

    def show_assessment_table(self, modeling_session):
        self.show_table(**modeling_session.assessment_table_data)


    ###############################################################################################
    # DOMAINS.                                                                                    #
    ###############################################################################################

    def launch_domain_splitting(self, pymod_element):
        protocol = Split_into_domains_protocol(self, pymod_element, output_directory=self.domain_analysis_dirpath)
        protocol.launch_from_gui()

    def launch_domain_fuse(self, pymod_element):
        protocol = Fuse_domains_protocol(self, pymod_element, output_directory=self.domain_analysis_dirpath)
        protocol.launch_from_gui()


    ###############################################################################################
    # SELECTION MENU COMMANDS.                                                                    #
    ###############################################################################################

    def select_all_from_main_menu(self):
        self.select_all_sequences()

    def select_all_sequences(self):
        for element in self.get_pymod_elements_list():
            element.widget_group.select_element(select_all=True)

    def deselect_all_from_main_menu(self):
        self.deselect_all_sequences()

    def deselect_all_sequences(self):
        for element in self.get_pymod_elements_list():
            element.widget_group.deselect_element(deselect_all=True)


    def show_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure():
                self.show_chain_in_pymol(element)

    def hide_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure():
                self.hide_chain_in_pymol(element)

    def select_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure() and not element.selected:
                element.widget_group.toggle_element()

    def deselect_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure() and element.selected:
                element.widget_group.toggle_element()

    def expand_all_clusters_from_main_menu(self):
        for element in self.get_cluster_elements():
            element.widget_group.expand_cluster()

    def collapse_all_clusters_from_main_menu(self):
        for element in self.get_cluster_elements():
            element.widget_group.collapse_cluster()


    ###############################################################################################
    # DISPLAY MENU COMMANDS.                                                                      #
    ###############################################################################################

    def change_font_from_action(self, font_type=None, font_size=None):
        self.main_window.update_font(font_type, font_size)

    def change_font_from_main_menu(self):
        """
        Opens a font selector widget from Qt.
        """
        font, font_is_valid = QtWidgets.QFontDialog.getFont()
        print("- Selected font:", font, font_is_valid)


    ###############################################################################################
    # HELP MENU COMMANDS.                                                                         #
    ###############################################################################################

    developer_email = "giacomo.janson@uniroma1.it"

    def show_about_dialog(self):

        # Initializes the message box.
        about_dialog = PyMod_QMessageBox(self.get_qt_parent())
        about_dialog.setIcon(QtWidgets.QMessageBox.Information)
        about_dialog.setWindowTitle(self.pymod_plugin_name)

        # Sets the main text.
        about_dialog.setText("Version: %s" % self.pymod_version + "." + self.pymod_revision)
        infomative_text = ('Copyright (C): 2020 Giacomo Janson, Alessandro Paiardini\n'
                           'Copyright (C): 2016 Giacomo Janson, Chengxin Zhang, Alessandro Paiardini'
                           '\n\nFor information on PyMod %s visit:\n'
                           '  http://schubert.bio.uniroma1.it/pymod/\n\n'
                           'Or send us an email at:\n %s' % (self.pymod_version, self.developer_email))
        about_dialog.setInformativeText(infomative_text)

        # Adds detailed information.
        pymod_plugin_path = os.path.dirname(os.path.dirname(pymod_lib.__file__))

        try:
            import PyQt5
            pyqt5_version = PyQt5.QtCore.PYQT_VERSION_STR
        except:
            pyqt5_version = "-"

        try:
            from pymol import Qt
            pymol_pyqt_name = Qt.PYQT_NAME
        except:
            pymol_pyqt_name = "-"

        try:
            import Bio
            biopython_version = Bio.__version__
        except:
            biopython_version = "-"

        try:
            import numpy
            numpy_version = numpy.__version__
        except:
            numpy_version = "-"

        try:
            import modeller
            modeller_version = modeller.__version__
            modeller_path = repr(modeller.__path__)
        except:
            modeller_version = "-"
            modeller_path = "-"

        try:
            import conda
            import conda.cli.python_api as conda_api
            conda_version = conda.__version__
            if self.DEVELOP:
                conda_info_dict = json.loads(conda_api.run_command(conda_api.Commands.INFO, "--json")[0])
                conda_info_text = "\n# Conda\n"
                for k in sorted(conda_info_dict.keys()):
                    conda_info_text += "- %s: %s\n" % (k, repr(conda_info_dict[k]))
        except:
            conda_version = "-"
            conda_info_dict = {}
            conda_info_text = ""

        has_pymol_conda = str(hasattr(pymol, "externing") and hasattr(pymol.externing, "conda"))

        def _get_path_string(path):
            _path = path
            if os.path.isdir(_path):
                return _path
            else:
                return _path + " (not found)"

        additional_text = ("# PyMod\n"
                           "- Version: " + self.pymod_version + "\n"
                           "- Revision: " + self.pymod_revision + "\n"
                           "- Plugin path: " + _get_path_string(pymod_plugin_path) + " \n"
                           "- Config directory: " + _get_path_string(self.cfg_directory_path) + "\n"
                           "- PyMod Directory: " + _get_path_string(self.current_pymod_dirpath) + "\n"
                           "- Current PyMod project: " + _get_path_string(self.current_project_dirpath) + "\n\n"

                           "# PyMOL\n"
                           "- Version: " + str(pymol.cmd.get_version()[0]) + "\n"
                           "- Path: " + sys.executable + "\n"
                           "- Qt: " + str(pymol_pyqt_name) + "\n"
                           "- Has Conda: " + has_pymol_conda + "\n\n"

                           "# Python\n"
                           "- Version: " + str(sys.version) + "\n"
                           "- Arch: " + pmos.get_python_architecture() + "\n"
                           "- Path: " + sys.executable + "\n\n"

                           "# Operating system\n"
                           "- Platform: " + sys.platform + "\n"
                           "- Arch: " + pmos.get_os_architecture() + "\n\n"

                           "# Python libs\n"
                           "- PyQt5: " + pyqt5_version + "\n"
                           "- Conda version: " + conda_version + "\n"
                           "- Numpy version: " + numpy_version + "\n"
                           "- Biopython version: " + biopython_version + "\n"
                           "- MODELLER version: " + modeller_version + "\n"
                           "- MODELLER path: " + modeller_path + "\n"
                          )

        if self.DEVELOP:
            additional_text += conda_info_text

        about_dialog.setDetailedText(additional_text)

        # Actually shows the message box.
        about_dialog.setModal(True)
        about_dialog.exec_()


    def open_online_documentation(self):
        webbrowser.open("http://schubert.bio.uniroma1.it/pymod/documentation.html")


    def launch_databases_update(self):
        db_updater = UpdaterProtocol(self)
        db_updater.launch_from_gui()


    ###############################################################################################
    # SESSIONS.                                                                                   #
    ###############################################################################################

    def exit_from_main_menu(self):
        self.main_window.confirm_close()


    def start_new_session_from_main_menu(self):
        title = "Begin New Project?"
        message = ("Are you really sure you want to begin a new PyMod project? If"
                   " you do not save your current project, its data will be"
                   " permanently lost.")
        answer = askyesno_qt(title, message, parent=self.get_qt_parent())
        if not answer:
            return None
        self.start_new_session()


    def save_session_from_main_menu(self):
        save_project_full_path = asksaveasfile_qt("Save PyMod Session file",
                                                  name_filter="*.pmse",
                                                  parent=self.get_qt_parent())
        if not save_project_full_path:
            return None
        self.save_pymod_session(save_project_full_path)


    def open_session_from_main_menu(self):
        project_archive_filepath = askopenfile_qt("Open a PyMod Session file",
                                                  name_filter="*.pmse",
                                                  parent=self.get_qt_parent())
        if not project_archive_filepath:
            return None
        if not os.path.isfile(project_archive_filepath):
            return None
        self.open_pymod_session(project_archive_filepath)


class PyMod_QMessageBox(QtWidgets.QMessageBox):

    def closeEvent(self, event):
        self.close()
