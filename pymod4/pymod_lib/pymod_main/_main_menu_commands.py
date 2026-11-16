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
##### MODIFIED #####
import requests
##### END #####

from Bio import Phylo
##### MODIFIED #####
from Bio import SeqIO
##### END #####
#################Add ETE code 21/05/2026-SARA###################
try:
    from ete3 import Tree, PhyloTree, AttrFace, TreeFace, TreeStyle, faces, TextFace, NodeStyle
except:
    pass
####################END ETE######################
from pymol.Qt import QtWidgets

import pymol

import pymod_lib
from pymod_lib import pymod_os_specific as pmos
from pymod_lib import pymod_vars
from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt, asksaveasfile_qt, askopenfile_qt, askopenfiles_qt
from pymod_lib.pymod_gui.specific_gui_components_qt import Raw_sequence_window_qt, PyMod_options_window_qt, Uniprot_id_window_qt, PDB_id_window_qt
##### END #####

from pymod_lib.pymod_plot import pymod_plot_qt
from pymod_lib.pymod_seq.seq_manipulation import (compute_sequence_identity,
                                                  get_invalid_typed_characters)
from pymod_lib.pymod_exceptions import PyModInvalidFile, PyModUnknownFile
from pymod_lib.pymod_gui.pymod_table import QtTableWindow
##### MODIFIED #####
from pymod_lib.pymod_os_specific import check_network_connection
from pymod_lib.pymod_os_specific import check_importable_openbabel
##### END #####

from pymod_lib.pymod_protocols.similarity_searches_protocols.ncbi_blast import NCBI_BLAST_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.local_blast import LOC_BLAST_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.psiblast import PSI_BLAST_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.phmmer import PHMMER_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.jackhmmer import Jackhmmer_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.hmmsearch import Hmmsearch_search
from pymod_lib.pymod_protocols.similarity_searches_protocols.import_hhsuite import Import_HHsuite_results

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
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.alignment_figure import Alignment_figure_analysis
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.tree_building import Tree_building
##################Add ETE code 21/05/2026-SARA#################
from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.tree_building import Tree_visualization
##########################END ETE##########################
from pymod_lib.pymod_protocols.structural_analysis_protocols.ramachandran_plot import Ramachandran_plot
from pymod_lib.pymod_protocols.structural_analysis_protocols.contact_map_analysis import Contact_map_analysis
from pymod_lib.pymod_protocols.structural_analysis_protocols.structural_divergence_plot import Structural_divergence_plot
from pymod_lib.pymod_protocols.structural_analysis_protocols.dope_assessment import DOPE_assessment, show_dope_plot
from pymod_lib.pymod_protocols.structural_analysis_protocols.secondary_structure_assignment import Secondary_structure_assignment
from pymod_lib.pymod_protocols.structural_analysis_protocols.psipred import PSIPRED_prediction

from pymod_lib.pymod_protocols.modeling_protocols.homology_modeling import MODELLER_homology_modeling
from pymod_lib.pymod_protocols.modeling_protocols.loop_refinement import MODELLER_loop_refinement
from pymod_lib.pymod_protocols.modeling_protocols.clash_and_refine import ClashAndRefine
from pymod_lib.pymod_gui.clash_refine_gui import ClashRefineWindow
from pymod_lib.pymod_protocols.modeling_protocols.molecule_optimizer import MoleculeOptimizer
from pymod_lib.pymod_protocols.modeling_protocols.AlphaFold import AFDB_search, show_af_plots
from pymod_lib.pymod_protocols.modeling_protocols.ESMFold import ESM_search, show_esm_plots
from pymod_lib.pymod_protocols.modeling_protocols.ESMFold.atlas_protocols import (
    ESM_atlas_search, build_plddt_plot_window, plot_plddt_of_model)
from pymod_lib.pymod_protocols.modeling_protocols.ESMFold.forge_protocols import ESM_forge_fold, ESM_inverse_fold
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import BOLTZ_model
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.design_protocol import BOLTZ_design
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.boltz_plots import show_boltz_plots
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.results_gui_qt import Boltz_results_window_qt

from pymod_lib.pymod_protocols.updater_protocols import UpdaterProtocol
################Add ETE code 21/05/2026-SARA################
from pymod_lib.pymod_plot.pyqtgraph.functions import makeQImage
########################END ETE###################

from pymod_lib.pymod_protocols.base_protocols import AFDB_common
from pymod_lib.pymod_residue import PyMod_residue


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
                                     name_filter="*.fasta *.fa *.fas *.gp *.gb *.genbank *.pdb *.ent *.cif *.mmcif *.sdf *.mol *.mol2 *.xyz",
                                     parent=self.get_qt_parent())

        # Loads each files in PyMod.
        loaded_sequences = False
        for single_file_path in file_paths:
            extension = os.path.splitext(single_file_path)[1].replace(".","").lower()
            try:
                # Sequence file formats
                if extension in ("fasta", "fa", "fas"):
                    self.open_sequence_file(single_file_path, "fasta")
                elif extension in ("gp", "gb", "genbank"):
                    self.open_sequence_file(single_file_path, "genbank")
                elif extension in ("embl"):
                    self.open_sequence_file(single_file_path, "embl")
                elif extension in ("swiss", "sw"):
                    self.open_sequence_file(single_file_path, "swiss")
                elif extension in ("pir"):
                    self.open_sequence_file(single_file_path, "pir")
                elif extension in ("phylip", "phy"):
                    self.open_sequence_file(single_file_path, "phylip")
                elif extension in ("nexus", "nex"):
                    self.open_sequence_file(single_file_path, "nexus")
                elif extension in ("clustal", "aln"):
                    self.open_sequence_file(single_file_path, "clustal")
                elif extension in ("stockholm", "sto"):
                    self.open_sequence_file(single_file_path, "stockholm")
                # Structure file formats
                elif extension in ("pdb", "ent"):
                    self.open_structure_file(single_file_path, "pdb")
                elif extension in ("cif", "mmcif"):
                    self.open_structure_file(single_file_path, "cif")
                elif extension in ("sdf", "mol", "mol2"):
                    self.open_structure_file_via_biopython(single_file_path, extension)
                elif extension in ("xyz"):
                    self.open_structure_file_via_biopython(single_file_path, extension)
                else:
                    raise PyModUnknownFile()
                loaded_sequences = True

            except PyModInvalidFile:
                title = "File Type Error"
                message = "The selected File ('%s') is not a valid file or cannot be processed." % single_file_path
                self.main_window.show_error_message(title, message)
            except PyModUnknownFile:
                title = "File Type Error"
                message = "The selected File ('%s') has an unknown file extension '%s'." % (single_file_path, extension)
                self.main_window.show_error_message(title, message)
            except Exception as e:
                title = "File Loading Error"
                message = "Could not load file '%s': %s" % (single_file_path, str(e))
                self.main_window.show_error_message(title, message)

        if loaded_sequences:
            self.main_window.gridder()


    def open_structure_file_via_biopython(self, file_path, file_format):
        """
        Handle structure files that require conversion or special processing via Biopython.
        For formats like SDF, MOL, MOL2, XYZ that aren't directly supported by PyMOL/PDB parsers.
        """
        try:
            # For now, these formats are not fully supported but we can provide a meaningful error
            supported_structure_formats = {
                'sdf': 'SDF (Structure Data Format)',
                'mol': 'MDL MOL',
                'mol2': 'Tripos MOL2',
                'xyz': 'XYZ coordinates'
            }
            
            title = "Format Not Yet Supported"
            format_name = supported_structure_formats.get(file_format, file_format.upper())
            message = (f"The {format_name} format is recognized but not yet fully supported in PyMod4. "
                      f"Please convert your file to PDB or CIF format for now.\n\n"
                      f"Future versions of PyMod4 will support this format directly.")
            self.main_window.show_warning_message(title, message)
            
        except Exception as e:
            title = "File Loading Error"
            message = f"Could not process {file_format.upper()} file '{file_path}': {str(e)}"
            self.main_window.show_error_message(title, message)

    #####################Add ETE code 21/05/2026-SARA######################
    def open_alignment_from_main_menu(self):
        """
        Lets users import in Pymod an alignment stored in an external file.
        """
        openfilename, extension = self.choose_alignment_file()
        if not None in (openfilename, extension):
            self.build_cluster_from_alignment_file(openfilename, extension)
        self.main_window.gridder(update_menus=True, update_elements=True)


    def open_tree_from_main_menu(self):

        file_paths = askopenfiles_qt("Select files to open",
                                     ##############MODIFICATION on 31/07/2025#########
                                     #Standardize extension recognition
                                     name_filter="*.ph *.phy *.phylip *.nw *.nwk *.newick *.tree",
                                     ####END###
                                     parent=self.get_qt_parent())

        loaded_trees = False
        for single_file_path in file_paths:
            extension = os.path.splitext(single_file_path)[1].replace(".","").lower()
            try:
                ##############MODIFICATION on 31/07/2025#########
                #Standardize extension recognition
                if extension in ("phylip", "ph", "phy", "nw", "nwk", "newick", "tree"):
                ####END###
                    ltree = Phylo.read(single_file_path, "newick")
                else:
                    raise PyModUnknownFile()
                loaded_trees = True

            except PyModInvalidFile:
                title = "File Type Error"
                message = "The selected File ('%s') is not a valid '%s' file." % (single_file_path, pymod_vars.supported_sequence_file_types[extension])
                self.main_window.show_error_message(title, message)
            except PyModUnknownFile:
                title = "File Type Error"
                message = "The selected File ('%s') has an unknown file extension '%s'." % (single_file_path, extension)
                self.main_window.show_error_message(title, message)

        if loaded_trees:
            self.loaded_trees_list.append(single_file_path)

            ############MODIFICATION SARA on 11/04/2026#####################
            self.main_window.show_info_message(
                "Importing tree",
                "The tree has been imported successfully."
            )

            ###########END MODIFICATION #######################
            self.main_window.gridder(update_menus=True, update_elements=True)
    #########################END ETE##########################



    ##################Original##################
    """
    def open_alignment_from_main_menu(self):
        
        #Lets users import in Pymod an alignment stored in an external file.
        
        openfilename, extension = self.choose_alignment_file()
        if not None in (openfilename, extension):
            self.build_cluster_from_alignment_file(openfilename, extension)
        self.main_window.gridder(update_menus=True, update_elements=True)
    """
        ###################END original###################

    #------------
    # HH-suite. -
    #------------

    def open_hhsuite_hhr_from_main_menu(self):
        import_protocol = Import_HHsuite_results(self, output_directory=self.similarity_searches_dirpath)
        import_protocol.launch_from_gui(mode="hhr")

    def open_hhsuite_a3m_from_main_menu(self):
        import_protocol = Import_HHsuite_results(self, output_directory=self.similarity_searches_dirpath)
        import_protocol.launch_from_gui(mode="a3m")


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


    ##### MODIFIED #####
    def show_uniprot_id_input_window(self):
        """
        Launched when the user wants to add a new sequence by giving the corresponding UniProt ID.
        """
        self.uniprot_id_window = Uniprot_id_window_qt(self.main_window,
                                                     title="Fetch by UniProt ID",
                                                     upper_frame_title="Type your UniProt ID",
                                                     submit_command=self.uniprot_id_input_window_state)
        self.uniprot_id_window.show()
    ##### END #####


    ################################MODIFICATION on 31/03/2025###############################################
    ##Enhanced sequence type identification for raw sequences.
    def raw_seq_input_window_state(self):
        """
        This is called when the SUBMIT button of the 'Add Raw Sequence' is pressed.
        """
        def name_match(strg, search2=re.compile(r'[^a-zA-Z0-9_]').search):
            """
            Checks if the sequence name contains only letters, numbers, and underscores.
            """
            return not bool(search2(strg))

        sequence = self.raw_seq_window.get_sequence()
        sequence_name = self.raw_seq_window.get_sequence_name()

        # What is not a residue. The check used to be a regular expression which let
        # through every letter, every digit and the underscore, so "MKT_1" was accepted
        # as a sequence and the digit and the underscore became residues of it.
        invalid_characters = get_invalid_typed_characters(sequence)

        # Validate the sequence
        if len(sequence) and not invalid_characters:
            # Validate the sequence name
            if len(sequence_name) and name_match(sequence_name):
                # Identify the polymer type (DNA, RNA, or protein)
                polymer_type = self.identify_sequence_type(sequence)

                # Create a new PyMod element and associate the polymer type
                element = self.build_pymod_element_from_args(sequence_name, sequence)
                element.polymer_type = polymer_type  # Assign polymer type to the element

                # Add the element to PyMod
                self.add_element_to_pymod(element)
                self.raw_seq_window.destroy()
                self.main_window.gridder()
            else:
                title = 'Sequence Name Error'
                message = 'Please check the sequence name: only letters, numbers, and "_" are allowed.'
                self.main_window.show_error_message(title, message)
        else:
            title = 'Sequence Error'
            if not len(sequence):
                message = 'Please type a sequence.'
            else:
                message = ('%s %s not a residue. A sequence may hold the standard amino'
                           ' acid letters, "X" for a residue which is not known, the'
                           ' bases of DNA and RNA (A, T, C, G, U) and "-" for a gap.'
                           % (", ".join('"%s"' % c for c in invalid_characters),
                              "is" if len(invalid_characters) == 1 else "are"))
            self.main_window.show_error_message(title, message)
    #####################END MODIFICATION###########################################


    ##### MODIFIED #####
    def uniprot_id_input_window_state(self):
        """
        This is called when the SUBMIT button of the 'Fetch by UniProt ID' is pressed.
        """
        if not check_network_connection("https://google.com", timeout=3):
            title = "Connection Error"
            message = ("An internet connection is not available, can not connect to the UniProt"
                       " server to find the corresponding sequence.")
            self.main_window.show_error_message(title, message)
            return None

        uniprot_id = self.uniprot_id_window.get_sequence_uniprot_id()
        ########################MODIFICATION on 28/02/2025############################################
        #####it can recognize also lower case chars for the uniprot ID search
        uniprot_id = uniprot_id.upper()  # Convert to uppercase before searching
        print(uniprot_id)

        uniprot_id_search = re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", uniprot_id)
        print(uniprot_id_search)
        
        ########################MODIFICATION on 28/02/2025############################################
        #########Supports the case where the user provides a dirty string that contains other characters besides the UNIPROT ID, otherwise extracts the UNIPROT code from that string###########
        if not uniprot_id_search:
            title = 'UniProt ID Error'
            message = 'Invalid UniProt ID. Please check that the UniProt ID inserted is correct.'
            self.main_window.show_error_message(title, message)
            return None
        
        uniprot_header = uniprot_id_search.group(0)
        path = os.path.join(self.temp_directory_dirpath, uniprot_header + ".fasta")
        #print(path)
        fastaurl = "https://rest.uniprot.org/uniprotkb/" + uniprot_header + ".fasta"
        response = requests.get(fastaurl)

        if response.status_code == 200 and response.content:
            file = open(path, "wb")
            file.write(response.content)
            file.close()
            with open(path) as fasta_handle:
                sequence_record = SeqIO.read(fasta_handle, 'fasta')
            #print(sequence_record)
            sequence = str(sequence_record.seq)
            self.add_element_to_pymod(self.build_pymod_element_from_args(uniprot_header, sequence))
            self.uniprot_id_window.destroy()
            self.main_window.gridder()

        else:
            title = 'UniProt ID Error'
            message = f'UniProt ID {uniprot_id} not found in the UniProt database.'
            self.main_window.show_error_message(title, message)

        """OLD CODE (It does not support the case where the user provides a messy string that contains other characters besides the UniProt ID.)
        if uniprot_id_search and len(uniprot_id_search.group(0)) == len(uniprot_id):
            uniprot_header = uniprot_id_search.group(0)
            path = os.path.join(self.temp_directory_dirpath, uniprot_header + ".fasta")
            #print(path)
            fastaurl = "https://rest.uniprot.org/uniprotkb/" + uniprot_header + ".fasta"
            response = requests.get(fastaurl)

            if response.status_code == 200 and response.content:
                file = open(path, "wb")
                file.write(response.content)
                file.close()
                sequence_record = SeqIO.read(path, 'fasta')
                #print(sequence_record)
                sequence = str(sequence_record.seq)
                self.add_element_to_pymod(self.build_pymod_element_from_args(uniprot_header, sequence))
                self.uniprot_id_window.destroy()
                self.main_window.gridder()

            else:
                title = 'UniProt ID Error'
                message = f'UniProt ID {uniprot_id} not found in the UniProt database.'
                self.main_window.show_error_message(title, message)

        else:
            title = 'UniProt ID Error'
            message = 'Invalid UniProt ID. Please check that the UniProt ID inserted is correct.'
            self.main_window.show_error_message(title, message)
            """
    ##### END #####


    def fetch_by_pdb_code_from_main_menu(self):
        """
        Launched when the user wants to add a new sequence by giving the corresponding PDB ID.
        """
        self.pdb_id_window = PDB_id_window_qt(self.main_window,
                                                     title="Fetch by PDB ID",
                                                     upper_frame_title="Type your PDB ID",
                                                     submit_command=self.pdb_id_input_window_state)
        self.pdb_id_window.show()

    def pdb_id_input_window_state(self):
        """
        This is called when the SUBMIT button of the 'Fetch by PDB ID' is pressed.
        """
        if not check_network_connection("https://files.rcsb.org", timeout=5):
            title = "Connection Error"
            message = ("An internet connection is not available, can not connect to the PDB"
                       " server to find the corresponding structure.")
            self.main_window.show_error_message(title, message)
            return None

        pdb_id = self.pdb_id_window.get_sequence_pdb_id().upper()

        if len(pdb_id) != 4:
            title = 'PDB ID Error'
            message = 'Invalid PDB ID. Please check that the PDB ID inserted is correct.'
            self.main_window.show_error_message(title, message)
            return None

        # Try to download as PDB format
        url_pdb = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response_pdb = requests.get(url_pdb)

        if response_pdb.status_code == 200:
            file_format = "pdb"
            file_content = response_pdb.text
            file_extension = ".pdb"
        else:
            # If PDB format fails, try CIF format
            url_cif = f"https://files.rcsb.org/download/{pdb_id}.cif"
            response_cif = requests.get(url_cif)
            if response_cif.status_code == 200:
                file_format = "cif"
                file_content = response_cif.text
                file_extension = ".cif"
            else:
                title = 'PDB ID Error'
                message = f'PDB ID {pdb_id} not found in the PDB database.'
                self.main_window.show_error_message(title, message)
                return

        # Save the downloaded content to a temporary file
        temp_file_path = os.path.join(self.temp_directory_dirpath, pdb_id + file_extension)
        with open(temp_file_path, "w") as f:
            f.write(file_content)

        try:
            self.open_structure_file(temp_file_path, file_format)
            self.new_objects_index += 1
            self.pdb_id_window.destroy()
            self.main_window.gridder()
        except Exception as e:
            print(e)
            title = 'File Loading Error'
            message = f'Could not load the structure for PDB ID {pdb_id}.'
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
                                               ################Add ETE code 21/05/2026###############
                                               name_filter="*.fasta;;*.aln;;*.sto;;*.phylip;;*.nexus",
                                               #############END ETE#################
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
        ####################Add ETE code 21/05/2026-SARA##################à
        elif extension == "phylip":
            self.build_sequence_file(aligned_elements, alignment_file_name, file_format="phylip", remove_indels=False)
        elif extension == "nexus":
            self.build_sequence_file(aligned_elements, alignment_file_name, file_format="nexus", remove_indels=False)
        #####################END ETE#####################À
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

        selected_sequences = self.get_selected_sequences()
        if len(selected_sequences) == 1:
            if selected_sequences[0].polymer_type in ["dna", "rna"] and blast_version == "blastp":
                self.main_window.show_error_message("Input Error", "Local BLAST search is not available for nucleotide sequences.")
                return

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
        # Check if any selected sequences are nucleic acids before launching the protocol
        selected_sequences = self.get_selected_sequences()
        if len(selected_sequences) > 0:
            nucleic_acid_sequences = [e for e in selected_sequences if e.polymer_type in ["dna", "rna", "nucleic_acid"]]
            if nucleic_acid_sequences:
                self.main_window.show_error_message("Selection Error", 
                    "PSIPRED secondary structure prediction is only available for protein sequences. "
                    "Please select only protein sequences.")
                return
        
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


    def launch_afdb_from_main_menu(self):
        afdb_session = AFDB_search(self)
        afdb_session.launch_from_gui()


    def launch_esm_from_main_menu(self):
        esm_session = ESM_search(self)
        esm_session.launch_from_gui()

    def launch_esm_atlas_from_main_menu(self):
        """
        Called by "Tools -> Structure Modeling -> ESMFold -> Fetch from the ESM Atlas".
        Retrieves a structure which the Atlas holds already, instead of folding one.
        """
        esm_atlas_session = ESM_atlas_search(self)
        esm_atlas_session.launch_from_gui()

    def launch_forge_fold_from_main_menu(self):
        ESM_forge_fold(self).launch_from_gui()

    def launch_inverse_fold_from_main_menu(self):
        ESM_inverse_fold(self).launch_from_gui()

    def launch_boltz_from_main_menu(self):
        boltz_session = BOLTZ_model(self)
        boltz_session.launch_from_gui()

    def import_boltz_prediction_by_id_from_main_menu(self):
        boltz_session = BOLTZ_model(self)
        boltz_session.launch_import_by_id_from_gui()

    def show_esm_atlas_figures(self, atlas_model):
        """
        The figures the ESM Atlas reports for a model which was retrieved from it.
        """

        from pymod_lib.pymod_protocols.modeling_protocols.ESMFold import atlas_client

        record = atlas_model.get("record", {})
        lines = ["Model: %s" % atlas_model.get("model_name", "?"), ""]
        if record.get("accession"):
            lines.append("Atlas entry: %s, from %s. The Atlas is keyed on the sequence, so"
                         " this is whichever entry deposited it, which is not always the"
                         " one the model is named after."
                         % (record["accession"], record.get("source", "the Atlas")))
        url = atlas_client.record_url(record)
        if url:
            lines.append("Record: %s" % url)
        if record.get("sequence_length"):
            lines.append("Length: %s residues" % record["sequence_length"])
        if record.get("mean_plddt") is not None:
            lines.append("Mean pLDDT: %.1f out of 100, the confidence of ESMFold2 in this"
                         " model, averaged over its residues"
                         % (float(record["mean_plddt"]) * 100))
        if record.get("ptm") is not None:
            lines.append("pTM: %.2f, its confidence in the fold as a whole; above about"
                         " 0.5 the fold is likely to be right" % float(record["ptm"]))
        if record.get("residues_plddt"):
            values = [float(v) * 100 for v in record["residues_plddt"]]
            lines.append("Per residue: from %.0f to %.0f; see the pLDDT plot for where the"
                         " model is confident and where it is not"
                         % (min(values), max(values)))
        if record.get("folded_on_demand"):
            lines.append("This entry was folded on demand rather than precomputed.")

        features = record.get("sae_features") or []
        if features:
            lines.append("")
            lines.append("What the model recognises in this protein, as the interpretable"
                         " features of ESMC name it:")
            for feature in features[:5]:
                label = feature.get("label", "?")
                regions = feature.get("residue_regions") or []
                where = ", ".join("%s-%s" % (r.get("start"), r.get("end"))
                                  for r in regions[:3])
                lines.append("  - %s%s" % (label, " (residues %s)" % where if where else ""))

        self.main_window.show_info_message("Figures of %s"
                                           % atlas_model.get("model_name", "the model"),
                                           "\n".join(lines))

    def show_esm_atlas_plddt(self, atlas_model):
        """
        The confidence of an Atlas model along its sequence, drawn as the other
        per-residue confidences are.
        """

        element = atlas_model.get("model_element") or atlas_model.get("element")
        model_name = atlas_model.get("model_name", "the model")
        try:
            path = element.get_structure_file(basename_only=False)
        except Exception:
            path = None
        if not path or not os.path.isfile(path):
            self.main_window.show_warning_message(
                "No plot to draw",
                "The confidence of %s is in the B-factor column of its structure, and that"
                " file could not be found." % model_name)
            return None

        cp = build_plddt_plot_window(self)
        plot_plddt_of_model(self, cp, model_name, path, element)

    def show_plddt_plot_boltz(self, boltz_session):
        show_boltz_plots(self, boltz_session, plot_type="plddt")

    def show_pae_plot_boltz(self, boltz_session):
        show_boltz_plots(self, boltz_session, plot_type="pae")

    def show_boltz_confidence_figures(self, boltz_session):
        """The table of confidence and binding figures of a Boltz prediction."""
        self.boltz_results_window = Boltz_results_window_qt(
            self.main_window,
            job=boltz_session.job,
            model_name=boltz_session.model_name,
            output_directory=boltz_session.output_directory,
            sample=getattr(boltz_session, "sample", None))
        self.boltz_results_window.show()

    def launch_foldseek_search(self):
        # Imported here rather than at the top of the module: it reaches the
        # network layer, and the menu is built before any search is ever run.
        from pymod_lib.pymod_protocols.similarity_searches_protocols.foldseek import (
            Foldseek_search)
        foldseek_session = Foldseek_search(self)
        foldseek_session.launch_from_gui()


    def launch_boltz_design_from_main_menu(self):
        boltz_design_session = BOLTZ_design(self)
        boltz_design_session.launch_from_gui()

    def optimize_from_main_menu(self):
        # Open Babel is what the minimisation runs on. When it is missing, PyMod offers to
        # install it rather than only saying that it is needed: see
        # offer_openbabel_installation() in the installer. Either way the command stops
        # here, because a module installed now is importable only after a restart.
        if not check_importable_openbabel():
            self.offer_openbabel_installation("Optimize Interface")
            return

        selected_elements = self.get_selected_elements()

        if len(selected_elements) != 2:
            self.main_window.show_warning_message(
                "Selection Error",
                "You must select exactly two objects to optimize an interface. "
                "If you want to minimize a single molecule, please use the 'optimize mol' function."
            )
            return

        elements_with_structure = [e for e in selected_elements if e.has_structure()]
        if len(elements_with_structure) != 2:
            self.main_window.show_warning_message(
                "Selection Error",
                "Both selected elements must have an associated 3D structure.")
            return

        # Lancia il tool passando l'intera istanza di Pymod e la lista di elementi.
        clash_refine = ClashAndRefine(self, selected_elements=elements_with_structure)
        clash_refine.run_clash_refine()


    def launch_molecule_optimizer(self):
        if not check_importable_openbabel():
            self.offer_openbabel_installation("Optimize Molecule")
            return

        selected_elements = self.get_selected_elements()

        if len(selected_elements) != 1:
            message = "You must select exactly one object to optimize."
            if len(selected_elements) == 2:
                # Two objects is the interface command's selection, so say so rather than
                # only saying no.
                message += ("\n\nTo relax the contact between two chains, please use the"
                            " 'Optimize Interface' tool first.")
            self.main_window.show_warning_message("Selection Error", message)
            return

        elements_with_structure = [e for e in selected_elements if e.has_structure()]
        if len(elements_with_structure) != 1:
            self.main_window.show_warning_message(
                "Selection Error",
                "The selected element must have an associated 3D structure.")
            return

        molecule_optimizer = MoleculeOptimizer(self, selected_elements=elements_with_structure)
        molecule_optimizer.run()


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

    # Kept, but not bound to any menu command anymore: the ESPript entry of the
    # "Render Alignment" submenu was removed, since that protocol depends on two
    # external servers and never returns its figure to PyMod. See 'espript.py' and
    # 'build_alignment_submenu' in 'pymod_gui/main_window/_main_menu_qt.py'.
    def launch_espript_from_main_menu(self, pymod_cluster):
        espript = ESPript_analysis(self, pymod_cluster)
        espript.launch_from_gui()

    def launch_alignment_figure_from_main_menu(self, pymod_cluster):
        alignment_figure = Alignment_figure_analysis(self, pymod_cluster)
        alignment_figure.launch_from_gui()


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
                                   sortable=sortable)
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
    ########################## Add ETE code 21/05/2026-SARA##################
    ##############MODIFICATION on 31/07/2025###########
    def build_tree_from_alignments_menu(self, alignment_element):
        """
        Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        menu.
        """
        tree_building = Tree_building(self, input_cluster_element=alignment_element)
        tree_building.launch_from_gui()
        ##############END############



    """OLD CODE
    def build_tree_from_alignments_menu(self, alignment_element):
        
        #Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        #menu.
        
        tree_building = Tree_building(self, input_cluster_element=alignment_element)
        tree_building.launch_from_gui()
    """

    def visualize_tree(self, tree):

        tree_visualization = Tree_visualization(self, input_cluster_element = tree)
        tree_visualization.launch_from_gui()

    ################MODIFICATION on 21/05/2026-SARA#################
    ### Remove Tree option (UNLOAD ONLY - DOES NOT DELETE FILE)
    def remove_tree(self, filepath):
        """
        Unload a phylogenetic tree from the PyMod session.

        This function performs the following steps:
        1. Converts the given filepath to its absolute path.
        2. Removes the tree from the global 'loaded_trees_list' if present.
        3. Removes the tree from the list of trees associated with any alignment (alignment_element.trees_list).
        4. Rebuilds the trees menu to reflect the removal.
        5. Shows an error message if the tree was not found in any list.
        
        NOTE: The original file on disk is NOT deleted.

        Parameters:
        filepath (str): The file path of the tree to be unloaded.
        """

        # Log the unload attempt
        print(f"[DEBUG] Trying to unload tree from session: {filepath}")

        # Normalize the file path to an absolute path
        target_path = os.path.abspath(filepath)

        # Flags to track if the tree was removed from internal lists
        removed_from_loaded_list = False
        removed_from_alignment = False

        # Try to remove the tree from the list of globally loaded trees
        for tree_path in list(self.loaded_trees_list):  # Make a copy to avoid modifying while iterating
            if os.path.abspath(tree_path) == target_path:
                self.loaded_trees_list.remove(tree_path)
                print(f"[DEBUG] Removed tree from loaded_trees_list: {tree_path}")
                removed_from_loaded_list = True

        # Try to remove the tree from all alignments' tree lists
        for alignment in self.get_cluster_elements():
            initial_count = len(alignment.trees_list)  # Save the original length
            # Keep only trees whose filepath doesn't match the one to remove
            alignment.trees_list = [tree for tree in alignment.trees_list if os.path.abspath(tree['filepath']) != target_path]
            # Check if anything was actually removed
            if len(alignment.trees_list) < initial_count:
                print(f"[DEBUG] Removed tree from alignment: {alignment.my_header}")
                removed_from_alignment = True

        # --- MODIFICA SARA: RIMOZIONE DELLA CANCELLAZIONE DAL DISCO ---
        # Il blocco os.remove è stato rimosso per evitare di cancellare il file originale.
        print(f"[DEBUG] File kept on disk: {target_path}")
        # -------------------------------------------------------------

        # Rebuild the "Trees" menu to reflect the changes in the GUI
        self.main_window.build_trees_submenu()

        # If the tree wasn't found in any internal list, inform the user
        if not removed_from_loaded_list and not removed_from_alignment:
            print(f"[DEBUG] Tree was not tracked in loaded_trees_list or alignments")
            self.main_window.show_error_message("Tree not found", f"The tree '{os.path.basename(filepath)}' was not found in this session.")
    ##########################END######################################
    """OLD CODE: working only with external imported tree and not with tree built from alignments
    def remove_tree(self, filepath):
        
        #Remove a tree from the loaded trees list by its filepath and update the Trees menu accordingly.
        #If the tree is found and removed, the Trees submenu is rebuilt to reflect the changes.
        
        
        # Print the current loaded trees for debugging purposes
        print("[DEBUG] Loaded trees list:")
        for i, tree_path in enumerate(self.loaded_trees_list):
            print(f"{i}: {tree_path}")
        
        # Print which tree filepath we are trying to remove
        print(f"[DEBUG] Trying to remove tree with filepath: {filepath}")

        # Get the absolute path of the target filepath for consistent comparison
        target_path = os.path.abspath(filepath)
        
        # Iterate over the loaded trees list to find the tree to remove
        for tree_path in self.loaded_trees_list:
            # Compare absolute paths to avoid issues with relative paths
            if os.path.abspath(tree_path) == target_path:
                # Remove the tree from the loaded trees list
                self.loaded_trees_list.remove(tree_path)
                print(f"[DEBUG] Removed tree: {filepath}")

                # Rebuild the Trees submenu to update the menu with the current tree list
                self.main_window.build_trees_submenu()

                # Exit the function early since the tree has been found and removed
                return

        # If the tree was not found in the list, print a debug message
        print(f"[DEBUG] Tree not found: {filepath}")
        """

    ###################END############


    def visualize_imported_tree(self, tree):

        tree_visualization = Tree_visualization(self, input_cluster_element = tree)
        tree_visualization.launch_from_gui()


    def save_tree_from_trees_menu(self, tree_file_path):
        output_filepath = asksaveasfile_qt("Save to Newick file",
                                            parent=self.get_qt_parent(),
                                            name_filter="*.nw")
        if not output_filepath:
            return None
        try:
            shutil.copy(tree_file_path,output_filepath)
        except Exception as e:
            print("- WARNING: could not write a newick file: %s" % str(e))

##########################END ETE#####################################


######################Original code#############################
    """
    def build_tree_from_alignments_menu(self, alignment_element):
        
        #Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        #menu.
        
        tree_building = Tree_building(self, input_cluster_element=alignment_element)
        tree_building.launch_from_gui()
    """
    ############### END original ##########################


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

    ##### MODIFIED #####
    ###AF
    def show_plddt_plot(self, af_session):
        show_af_plots(self, af_session, type = "plddt")

    def show_pae_plot(self, af_session):
        show_af_plots(self, af_session, type = "pae")

    ################MODIFICATION on 24/07/2025##################
    ###########ESMFold##################
    def show_plddt_plot_ESM(self, esm_session):
        show_esm_plots(self, esm_session, type = "plddt")

    def show_pae_plot_ESM(self, esm_session):
        show_esm_plots(self, esm_session, type = "pae")
    #################END###################
    



    # def show_plddt_plot(self, af_model):
    #     show_af_plots(self, af_model, type = "plddt")
    # 
    # def show_pae_plot(self, af_model):
    #     show_af_plots(self, af_model, type = "pae")
    ##### END #####



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

    # The address the user's guide gives for PyMod 4 ("For more information, please
    # contact:"). It used to be the one of the author of PyMod 3, which is where a bug
    # report of PyMod 4 went when it was sent from this box.
    developer_email = "alessandro.paiardini@uniroma1.it"

    def show_about_dialog(self):

        # Initializes the message box.
        about_dialog = PyMod_QMessageBox(self.get_qt_parent())
        about_dialog.setIcon(QtWidgets.QMessageBox.Information)
        about_dialog.setWindowTitle(self.pymod_plugin_name)

        # Sets the main text. The release, as the rest of the interface names it;
        # the exact version and revision are in the details below, which is what a
        # bug report needs them for.
        about_dialog.setText(self.pymod_plugin_name)
        # The name of the plugin is above this text and in the title of the
        # window, so it is not written again here: only the release which each
        # copyright belongs to, and only where there is more than one.
        infomative_text = ('Copyright (C): 2026 Serena Rosignoli, Sara Giordani,'
                           ' Elisa Lustrino, Alessandro Paiardini\n\n'
                           'Earlier releases\n'
                           'Copyright (C): 2020 Giacomo Janson, Alessandro Paiardini\n'
                           'Copyright (C): 2016 Giacomo Janson, Chengxin Zhang,'
                           ' Alessandro Paiardini'
                           '\n\nFor more information visit:\n'
                           '  http://schubert.bio.uniroma1.it/pymod/\n\n'
                           'Or send us an email at:\n %s' % (self.developer_email,))
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
            # The deprecated 'conda.cli.python_api' has been replaced by a
            # subprocess-based wrapper exposed as pmos.conda_api.
            conda_api = pmos.conda_api
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

        additional_text = ("# PyMod\n"                           "- Version: " + self.pymod_version + "\n"                           "- Revision: " + self.pymod_revision + "\n"                          "- Plugin path: " + _get_path_string(pymod_plugin_path) + " \n"                           "- Config directory: " + _get_path_string(self.cfg_directory_path) + "\n"                           "- PyMod Directory: " + _get_path_string(self.current_pymod_dirpath) + "\n"                           "- Current PyMod project: " + _get_path_string(self.current_project_dirpath) + "\n\n"                           
                           "# PyMOL\n"                           "- Version: " + str(pymol.cmd.get_version()[0]) + "\n"                           "- Path: " + sys.executable + "\n"                           "- Qt: " + str(pymol_pyqt_name) + "\n"                           "- Has Conda: " + has_pymol_conda + "\n\n"                           
                           "# Python\n"                           "- Version: " + str(sys.version) + "\n"                           "- Arch: " + pmos.get_python_architecture() + "\n"                           "- Path: " + sys.executable + "\n\n"                           
                           "# Operating system\n"                           "- Platform: " + sys.platform + "\n"                           "- Arch: " + pmos.get_os_architecture() + "\n\n"                           
                           "# Python libs\n"                           "- PyQt5: " + pyqt5_version + "\n"                           "- Conda version: " + conda_version + "\n"                           "- Numpy version: " + numpy_version + "\n"                           "- Biopython version: " + biopython_version + "\n"                           "- MODELLER version: " + modeller_version + "\n"                           "- MODELLER path: " + modeller_path + "\n"                          )

        if self.DEVELOP:
            additional_text += conda_info_text

        about_dialog.setDetailedText(additional_text)

        # Actually shows the message box.
        about_dialog.setModal(True)
        about_dialog.exec_()


    def open_online_documentation(self):
        webbrowser.open("https://github.com/pymodproject/pymod#how-to-use-pymod")


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