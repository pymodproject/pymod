# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for performing modeling or database search with AlphaFold in PyMod.
"""

import os
import sys
import shutil
import importlib
import pickle
import re
import urllib.request
import json
import requests
import Bio.PDB
import numpy as np
import math

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymol import cmd

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.__gui_qt import Boltz_search_window_qt
from pymod_lib.pymod_protocols.structural_databases_protocols import Associate_structure

##### MODIFIED #####
from pymol.Qt import QtWidgets, QtCore
##### END #####


###################################################################################################
# Boltz Modeling.                                                                      #
###################################################################################################


class BOLTZ_model(PyMod_protocol):

    """
    Class to represent a modeling run with Boltz-1.
    """


    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)

        self.cut_structure = False

        self.boltz_initialization()

        self.selected_sequences = self.pymod.get_selected_sequences()


    def boltz_initialization(self):
        self.boltz_modeling_window_class_qt = Boltz_search_window_qt


    def launch_from_gui(self):
        """
        This method is called when the "Database search" option is clicked in the "Tools" menu.
        """

        self.initialize_modeling_protocol_from_gui()


    ###########################################################################
    # Checks the input selection.                                             #
    ###########################################################################

    def initialize_modeling_protocol_from_gui(self):

        #---------------------------------------------------------------------
        # Get the selected sequences to see if there is a correct selection. -
        #---------------------------------------------------------------------

        # First check if at least one sequence is selected.
        if not len(self.selected_sequences) > 0:
            self.pymod.main_window.show_error_message("Selection Error", "Please select at least one target sequence to use Boltz-1 modeling tool.")
            return None

        # Checks if all the selected sequences can be used to build a model.
        if False in [s.can_be_modeled() for s in self.selected_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only sequences that do not have a structure loaded in PyMOL.")
            return None

        self.build_boltz_modeling_window()


    def build_boltz_modeling_window(self):

        """
        Builds the Boltz-1 modeling tool window.
        """

        self.boltz_modeling_window = self.boltz_modeling_window_class_qt(parent=self.pymod.main_window,
                                                             protocol=self, selected_sequences = self.selected_sequences)
        self.boltz_modeling_window.show()


    def run_boltz_modeling(self):

        # Make Session directory 
        models_dir = os.path.join(self.pymod.current_project_dirpath, self.pymod.models_dirname)
        model_subdir_name = "%s_%s_%s" % (self.pymod.models_subdirectory, self.pymod.performed_boltz_modeling_count + 1, "_boltz")
        self.boltz_output_dir_path = os.path.join(models_dir, model_subdir_name)
        self.boltz_input_fasta = os.path.join(self.boltz_output_dir_path, "inputs.fasta")
        if not os.path.isdir(self.boltz_output_dir_path):
            os.mkdir(self.boltz_output_dir_path)

        # Get inputs from Boltz GUI
        self.get_input_from_boltz_modeling_window()

        # Format selected sequences (protein)
        boltz_inputs = []

        for selseq in self.selected_sequences:

            if selseq.polymer_type == "dna":
                boltz_inputs.append(self.format_inputs(selseq.structure.chain_id, "dna", None, selseq.my_sequence))
            elif selseq.polymer_type == "rna":
                boltz_inputs.append(self.format_inputs(selseq.structure.chain_id, "rna", None, selseq.my_sequence))
            else:
                boltz_inputs.append(self.format_inputs(selseq.structure.chain_id, "protein", None, selseq.my_sequence))

            # TODO: # Format selected sequences (dna/rna)

        for selhetres in self.selected_hetres:

            #TODO: chain number
            boltz_inputs.append(self.format_inputs("A", "ccd", None, selhetres.text()))

        
        # Format inputs as a dictionary
        # boltz_inputs = [
        # self.format_inputs("A", "protein", "./examples/msa/seq1.a3m", "MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ"),
        # self.format_inputs("B", "protein", None, "MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ"),
        # self.format_inputs("C", "ccd", None, "SAH"),
        # self.format_inputs("E", "smiles", None, "N[C@@H](Cc1ccc(O)cc1)C(=O)O")
        # ]

        # Make FASTA file as Boltz input 
        self.create_fasta_boltz_input(self.boltz_input_fasta, boltz_inputs)


        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True, update_elements=True)


    def get_input_from_boltz_modeling_window(self):

        # If the selected sequence is part of a cluster, Pymod asks whether to build a new alignment or to extend the existing one
        # if self.selected_sequences[0].is_child():
        #     if self.boltz_modeling_window.al_mode_cb.currentText() == "Build new alignment":
        #         self.new_sequences_import_mode = "build-new"
        #     else:
        #         self.new_sequences_import_mode = "expand"
        # else:
        #     self.new_sequences_import_mode = "build-new"

        # Get the search mode
        # if self.boltz_modeling_window.sequence_id_rb.isChecked():
        #     self.search_choice = "by_id"
        # else:
        #     self.search_choice = "by_sequence"

        # # The user is asked whether to import the whole structure or just the hit fragment
        # if self.boltz_modeling_window.import_mode_cb.currentText() == "Whole Structure":
        #     self.new_structure_import_mode = "all_structure"
        # else:
        #     self.new_structure_import_mode = "hit_fragment"

        self.selected_hetres = [a for a in self.boltz_modeling_window.list_of_cb if a.isChecked()]


    def create_fasta_boltz_input(self, output_file, chains):
        """
        Creates a FASTA file in the specified format.

        Parameters:
        - output_file (str): Path to the output FASTA file.
        - chains (list of dict): List of chain details. Each dict should have the keys:
            - 'chain_id' (str): Unique identifier for the chain.
            - 'entity_type' (str): Type of entity (protein, dna, rna, smiles, ccd).
            - 'msa_path' (str or None): Path to MSA file, "empty", or None for non-protein chains.
            - 'sequence' (str): Sequence (amino acid, nucleotide, SMILES, or CCD code).

        Example:
        chains = [
            {"chain_id": "A", "entity_type": "protein", "msa_path": "./examples/msa/seq1.a3m", "sequence": "MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ"},
            {"chain_id": "C", "entity_type": "ccd", "msa_path": None, "sequence": "SAH"}
        ]
        create_fasta_file("output.fasta", chains)
        """
        with open(output_file, 'w') as fasta:
            for chain in chains:
                chain_id = chain['chain_id']
                entity_type = chain['entity_type']
                msa_path = chain['msa_path'] if chain['msa_path'] is not None else 'empty' if entity_type == 'protein' else ''
                sequence = chain['sequence'].replace('-', '') 

                header = f">{chain_id}|{entity_type}|{msa_path}"
                fasta.write(header + "\n")
                fasta.write(sequence + "\n")

    def format_inputs(self, chain_id, entity_type, msa_path, sequence):
        """
        Creates a dictionary for a chain entry.

        Parameters:
        - chain_id (str): Unique identifier for the chain.
        - entity_type (str): Type of entity (protein, dna, rna, smiles, ccd).
        - msa_path (str or None): Path to MSA file, "empty", or None for non-protein chains.
        - sequence (str): Sequence (amino acid, nucleotide, SMILES, or CCD code).

        Returns:
        - dict: A dictionary containing the chain details.
        """
        if entity_type == 'protein' and msa_path is None:
            msa_path = 'empty'
        return {
            "chain_id": chain_id,
            "entity_type": entity_type,
            "msa_path": msa_path,
            "sequence": sequence
        }


    def connect_to_afdb(self):
        pass


    def by_id_mode_import_results_in_pymod(self):

        """
        Builds a cluster with the query sequence as a mother and retrieved hits as children.
        """

        self.query_element = self.selected_sequences[0]

        # The list of elements whose sequences will be updated according to the star alignment.
        elements_to_update = [self.query_element]

        #------------------------------------------------------------
        # Builds a new cluster with the query and all the new hits. -
        #------------------------------------------------------------

        if self.new_sequences_import_mode == "build-new":

            # Gets the original index of the query element in its container.
            query_container = self.query_element.mother
            query_original_index = self.pymod.get_pymod_element_index_in_container(self.query_element)

            # Gets the name and sequence of the query_element
            name = self.pymod.get_header(self.query_element)
            sequence = self.query_element.my_sequence

            # Creates PyMod element
            # Gives them the query mother_index, to make them its children.
            # cs = self.pymod.build_pymod_element_from_args(str(self.name),
            #                                              sequence)

            # Initialize a PyMod element of the new sequence
            self.pymod.add_element_to_pymod(self.new_model_element)
            elements_to_update.append(self.new_model_element)


            self.query_element.selected = True
            self.new_model_element.selected = True

            selected_cluster = self.query_element.get_ancestor()

            # for seq_elements in selected_cluster.get_children():
            #     seq_elements.selected = True

            a = Clustalw_regular_alignment(self.pymod, output_directory=self.pymod.alignments_dirpath)
            a.rebuild_single_alignment_choice = False
            a.involved_clusters_list = []
            a.alignment_mode = "build-new-alignment"
            a.target_cluster_index = 0
            ###
            a.params_from_gui = {}
            a.params_from_gui["selected_matrix"] = "blosum"
            a.params_from_gui["gapopen_value"] = 10.0
            a.params_from_gui["gapextension_value"] = 0.2
            a.elements_to_align_dict = {}
            a.AF = True

            a.launch_no_gui()

            sibling_elements = []

        #----------------------------------------------------------------------------
        # Expand the original cluster of the query by appending to it the new hits. -
        #----------------------------------------------------------------------------

        elif self.new_sequences_import_mode == "expand":

            # The list of elements whose sequences will be updated according to the star alignment.
            elements_to_update = []
            # Begins with the query element.
            elements_to_update.append(self.query_element)
            sibling_elements = self.query_element.get_siblings(sequences_only=True)
            elements_to_update.extend(sibling_elements)

            new_cluster = self.query_element.mother
            # Creates PyMod elements for all the imported hits and add them to the cluster.

            name = self.pymod.get_header(self.query_element)
            sequence = self.query_element.my_sequence

            self.pymod.add_element_to_pymod(self.new_model_element)
            elements_to_update.append(self.new_model_element)

            self.query_element.selected = True
            self.new_model_element.selected = True

            selected_cluster = self.query_element.get_ancestor()

            # for seq_elements in selected_cluster.get_children():
            #     seq_elements.selected = True

            a = Clustalw_regular_alignment(self.pymod, output_directory=self.pymod.alignments_dirpath)
            a.rebuild_single_alignment_choice = False
            a.involved_clusters_list = [selected_cluster]
            a.alignment_mode = "keep-previous-alignment"
            a.target_cluster_index = 0
            ###
            a.params_from_gui = {}
            a.params_from_gui["selected_matrix"] = "blosum"
            a.params_from_gui["gapopen_value"] = 10.0
            a.params_from_gui["gapextension_value"] = 0.2
            a.elements_to_align_dict = {}
            a.AF = True
            # {'selected_matrix': 'blosum', 'gapopen_value': 10.0, 'gapextension_value': 0.2}

            a.launch_no_gui()

        if self.cut_structure:
            #-------------------------------------------------------------------------------------
            # Gets information about matching and missing residues in the two aligned sequences. -
            #-------------------------------------------------------------------------------------
            pc = 0 # Alignment position counter.
            hc = 0 # Target residue counter.
            tc = 0 # PDB structure residue counter.
            matching_positions = [] # list of matching positions.
            missing_positions = [] # list of missing residues in the pdb structure with respect to the target sequence.
            for hr, tr in zip(self.query_element.my_sequence, self.new_model_element.my_sequence):
                if hr != "-" and tr != "-" and hr == tr:
                    matching_positions.append({"pc": pc, "hc": hc, "tc": tc})
                if tr == "-" and hr != "-":
                    missing_positions.append({"pc": pc, "hc": hc, "tc": tc})
                if hr != "-":
                    hc += 1
                if tr != "-":
                    tc += 1
                pc += 1

            "{'pc': 110, 'hc': 0, 'tc': 110}, {'pc': 111, 'hc': 1, 'tc': 111}"

            starting_pos = matching_positions[0]["pc"]
            ending_pos = matching_positions[-1]["pc"]

            if self.new_sequences_import_mode == "expand":

                matching_seq = self.new_model_element.my_sequence[starting_pos:ending_pos+1]
                cropped_sequence = "-" * starting_pos + matching_seq + "-" * (ending_pos - (starting_pos + len(matching_seq)))
                self.new_model_element.my_sequence = cropped_sequence

                matching_seq = self.query_element.my_sequence[starting_pos:ending_pos+1]
                cropped_sequence = "-" * starting_pos + matching_seq + "-" * (ending_pos - (starting_pos + len(matching_seq)))
                self.query_element.my_sequence = cropped_sequence

            else:

                matching_seq = self.new_model_element.my_sequence[starting_pos:ending_pos+1]
                cropped_sequence = matching_seq + "-" * (ending_pos - (starting_pos + len(matching_seq)))
                self.new_model_element.my_sequence = cropped_sequence

                matching_seq = self.query_element.my_sequence[starting_pos:ending_pos+1]
                cropped_sequence = matching_seq + "-" * (ending_pos - (starting_pos + len(matching_seq)))
                self.query_element.my_sequence = cropped_sequence


            fp = Fetch_models_file(self.pymod)
            fp.initialize_from_gui("single", self.new_model_element)
            fp.import_mode = "single-chain"
            fp.fetch_models_files()

        else:

            fp = Fetch_models_file(self.pymod)
            fp.initialize_from_gui("single", self.new_model_element)
            fp.import_mode = "multiple-chains"
            fp.fetch_models_files()

        fp.new_element_with_structure._is_af_model = True

        #------------------------
        # Color imported models -
        #------------------------

        # List of symbols to be assigned to a color
        af_scores_dict = {}
        af_scores_dict["-"] = 1
        af_scores_dict["*"] = 3
        af_scores_dict[":"] = 3
        af_scores_dict["."] = 2
        af_scores_dict["missing"] = None

        selected_cluster = self.query_element.get_ancestor()

        # List elements to compare (the query and the new Pymod element associated to the model)
        new_cluster_sequences_list = [self.query_element, fp.new_element_with_structure]

        # List of columns for each position in the alignment
        list_of_alignment_columns = []

        # Length of alignment
        alignment_length = len(new_cluster_sequences_list[0].my_sequence)

        # Build the list of alignment columns
        for num in range(alignment_length):
            temp_list = []
            for cluster in new_cluster_sequences_list:
                temp_list.append(cluster.my_sequence[num])
            list_of_alignment_columns.append(temp_list)

        # Get the residues of the sequence to be colored
        residues = fp.new_element_with_structure.get_polymer_residues()

        # 'rc' keeps track of the residue index of the sequence to be modeled
        rc = 0

        # iterate in parallel over the sequence of the new element (the sequence as it is in the alignment) and the list of alignment columns
        for num, (seq_res, column) in enumerate(zip(fp.new_element_with_structure.my_sequence, list_of_alignment_columns)):
            if column[0] == "-" and column[1] != "-":
                symbol = get_conservation_symbol(column)
                residues[rc].af_score = af_scores_dict[symbol]
                rc += 1
            elif column[0] != "-" and column[1] != "-":
                if column[0] == column[1]:
                    residues[rc].af_score = af_scores_dict["*"]
                    rc += 1
                else:
                    residues[rc].af_score = af_scores_dict["."]
                    rc += 1

        self.new_model_element._has_af_scores = True
        self.pymod.main_window.color_element_by_af_scores(fp.new_element_with_structure)

        return "prova"


    def run_sequence_similarity_search(self):

        self.pymod.af_blast_search = LOC_BLAST_search(self.pymod, output_directory=self.pymod.similarity_searches_dirpath)
        self.pymod.af_blast_search.launch_from_gui(AF = True)

        #self.prova()

    def prova(self):

        selected_cluster = self.query_element.get_ancestor()

        # for seq_elements in selected_cluster.get_children():
        #     seq_elements.selected = True

        a = Clustalw_regular_alignment(self.pymod, output_directory=self.pymod.alignments_dirpath)
        a.rebuild_single_alignment_choice = False
        a.involved_clusters_list = [selected_cluster]
        a.alignment_mode = "keep-previous-alignment"
        a.target_cluster_index = 0
        ###
        a.params_from_gui = {}
        a.params_from_gui["selected_matrix"] = "blosum"
        a.params_from_gui["gapopen_value"] = 10.0
        a.params_from_gui["gapextension_value"] = 0.2
        a.elements_to_align_dict = {}
        a.AF = True
        # {'selected_matrix': 'blosum', 'gapopen_value': 10.0, 'gapextension_value': 0.2}

        a.launch_no_gui()


        return self.pymod.af_blast_search


##### MODIFIED #####

###############################################
# Window for selecting the PAE maps to plot #
###############################################

class Select_pae_window_qt(QtWidgets.QMainWindow):
    """
    Window with a checkbox for each model present in that AF modeling session. Select a model to obtain the corresponding PAE plot.
    """

    is_pymod_window = True

    def __init__(self, parent, af_session):

        super(Select_pae_window_qt, self).__init__(parent)
        self.af_session = af_session
        self.initUI()


    def initUI(self):

        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle("Model selection for plotting PAE")

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()

        ################
        # Upper frame. #
        ################

        self.upper_frame_title = QtWidgets.QLabel("Select the model for which you want to plot the PAE")
        self.main_vbox.addWidget(self.upper_frame_title)

        self.model_group = QtWidgets.QGroupBox()
        self.model_vbox = QtWidgets.QVBoxLayout()
        self.model_group.setLayout(self.model_vbox)

        self.models_in_session = ["AF-O14965_chain_A"]

        for model in self.models_in_session:
            self.checkbox = QtWidgets.QCheckBox(model)
            self.model_vbox.addWidget(self.checkbox)

        self.main_vbox.addWidget(self.model_group)

        self.central_widget.setLayout(self.main_vbox)

##### END #####



##### MODIFIED #####
def show_af_plots(self, af_session, type):

    # print(dir(self))
    # print(dir(af_session))

    if type == "plddt":

        """
        The predicted local distance difference test (pLDDT) score (0-100) is a per-residue confidence score,
        with values greater than 90 indicating high confidence, and values below 50 indicating low confidence.
        AlphaFold models are often shown with high confidence residues colored blue, and lower confidence in yellow,
        orange and red. This measure estimates whether the predicted residue has similar distances to neighboring
        C-alpha atoms (within 15 Angstroms) in agreement with distances in the true structure.
        """
        cp = PyMod_plot_window_qt(self.main_window)
        cp.initialize(pymod=self, title="Predicted Local Distance Difference Test (pLDDT) per Position")
        messagebar_text_on_update_qt = "Selected: __residue_name__ __residue_pdb_position__ of __plot_name__, pLDDT value: __y__"
        cp.build_plotting_area(messagebar_initial_text="Set the 'Interact on click' option to 'Yes' to click on the plot and to highlight residues in PyMOL.",
                               update_messagebar=True,
                               messagebar_text_on_update=messagebar_text_on_update_qt,
                               on_click_action=_highlight_in_pymol_from_plddt_plot_qt,
                               x_label_text="Positions",
                               y_label_text="pLDDT")

        for af_model in self.af_modeling_session_list[af_session]:
            print(af_model)
            af_file_dirpath = self.af_modeling_session_list[af_session][af_model]["af_session_dirpath"] + "/" + af_model + ".pdb"
            print(af_file_dirpath)
            model_element = self.af_modeling_session_list[af_session][af_model]["model_element"]
            print(model_element)

            structure = Bio.PDB.PDBParser(QUIET=True).get_structure(af_model, af_file_dirpath)

            pLDDT_data = []
            res_numbers = []
            res_names = []
            pymol_selectors = []

            for model in structure:

                for chain in model:
                    plot_name = af_model

                    for residue in chain:
                        #print(dir(residue))
                        pymol_selector = f"{structure.get_id()} and resi {residue.get_id()[1]}"
                        pymol_selectors.append(pymol_selector)

                        if "CA" in residue:
                            alpha_carbon = residue["CA"]
                            b_factor = alpha_carbon.get_bfactor()
                            pLDDT_data.append(b_factor)

                            res_number = residue.get_id()[1]
                            res_numbers.append(res_number)

                            res_name = residue.get_resname()
                            res_names.append(res_name)

            # print(pLDDT_data)
            # print(res_numbers)
            # print(res_names)
            # print(pymol_selectors)

            residue_additional_data = []
            for r_name, r_position, r_selector in zip(res_names, res_numbers, pymol_selectors):
                residue_additional_data.append({"residue_name": r_name,
                                                    "residue_pdb_position": r_position,
                                                    "pymol_selector": r_selector,
                                                    "export_label": "%s %s"%(r_name, r_position)})
            #print(residue_additional_data)

            plddt_scores_dict = {} # Items will contain DOPE scores of only polymer residues.
            # Retain only the pLDDT values for residues of the chain (standard and modified residues).

            plddt_scores = []
            for res, score in zip(model_element.residues, pLDDT_data):
                if res.is_polymer_residue():
                    plddt_scores.append(score)
            plddt_scores_dict.update({model_element: plddt_scores})

            # Builds a list of all pLDDT values of the residues in the selection.
            l_plddt = []

            l_plddt.extend(plddt_scores_dict[model_element])

            # Takes the min and max values among all the selected residues.
            min_value = min(l_plddt)
            max_value = max(l_plddt)
            # An array with the equally sapced limits generated with the list above.
            bins = np.array(np.linspace(min_value, max_value, num=10))

            # An array with all the DOPE values of a single chain in the selection.
            a_plddt = np.array(plddt_scores_dict[model_element])
            # An array with the id of the bins where those values reside.
            inds = np.digitize(a_plddt, bins)

            plddt_items = []
            for plddt_score, bin_id in zip(a_plddt, inds): # zip(ldope, inds):
                plddt_items.append({"score": plddt_score, "interval": bin_id})
            #print(plddt_items)
            # Actually assigns the pLDDT score to the residues of the PyMod element.
            for res, plddt_item in zip(model_element.get_polymer_residues(), plddt_items):
                res.plddt_score = plddt_item
            model_element._has_plddt_scores = True

            self.main_window.color_element_by_plddt(model_element)

            cp.add_plot(res_numbers, pLDDT_data, label = plot_name, additional_data = residue_additional_data)
            cp.show()
            print("plotting plddt")


    elif type == "pae":

        """
        Predicted aligned error (PAE) gives a distance error for every pair of residues.
        It gives AlphaFold's estimate of position error at residue x when the predicted and true structures are aligned on residue y.
        Values range from 0 - 35 Angstroms. It is usually shown as a heatmap image with residue numbers running along vertical and
        horizontal axes and color at each pixel indicating PAE value for the corresponding pair of residues.
        If the relative position of two domains is confidently predicted then the PAE values will be low (less than 5A)
        for pairs of residues with one residue in each domain.
        """

        self.select_pae_window = Select_pae_window_qt(parent=self.main_window, af_session=af_session)

        self.select_pae_window.show()



        if not check_network_connection("https://google.com", timeout=3):
            title = "Connection Error"
            message = ("An internet connection is not available, can not connect to the AFDB"
                       " to download the corresponding PAE values.")
            self.main_window.show_error_message(title, message)
            return None

        for af_model in self.af_modeling_session_list[af_session]:
            print(af_model)
            af_file_dirpath = self.af_modeling_session_list[af_session][af_model]["af_session_dirpath"] + "/" + af_model + ".pdb"
            print(af_file_dirpath)
            uniprot_header = os.path.basename(os.path.normpath(self.af_modeling_session_list[af_session][af_model]["af_session_dirpath"])).rsplit('-', 2)[1]
            print(uniprot_header)
            model_element = self.af_modeling_session_list[af_session][af_model]["model_element"]
            print(model_element)

            # Link to AFDB
            link_to_afdb_api = "https://alphafold.ebi.ac.uk/api/prediction/"

            # Make full link with Uniprot-ID of interest
            full_link = link_to_afdb_api + str(uniprot_header)
            print(full_link)

            # Open link and get data in JSON format
            response = urllib.request.urlopen(full_link)
            data = json.loads(response.read())

            # Get url, from the JSON data, to download the document containing pae values
            #TODO: SEE WHAT HAPPENS IF ENTRYid IS NOT FOUND
            url_to_pae_doc = data[0]["paeDocUrl"]
            print(url_to_pae_doc)

            # Download PAE values corresponding to AlphaFold model in json format
            response = requests.get(url_to_pae_doc)
            pae_doc_dirpath = os.path.join(self.af_session_dirpath, "AF-" + uniprot_header + ".json")
            pae_file = open(pae_doc_dirpath, "wb")
            pae_file.write(response.content)
            pae_file.close()

            # To open the file containing PAE values in reading mode
            with open(pae_doc_dirpath,'r') as file:
                data = json.load(file)

            # if the json file contains a list with one single element, that is a dictionary
            if isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict):
                # to obtain the dictionary 'pae_data' containing the values of PAE
                pae_data = data[0]
                # to get the values associated to the key 'predicted_aligned_error'
                pae_values = pae_data.get("predicted_aligned_error")
                #print(pae_values)

            # if the json file directly contains a dictionary
            else:
                # to get the values associated to the key 'predicted_aligned_error'
                pae_values = data["predicted_aligned_error"]


            structure = Bio.PDB.PDBParser(QUIET=True).get_structure(af_model, af_file_dirpath)

            res_numbers = []

            for model in structure:

                for chain in model:

                    for residue in chain:

                        if "CA" in residue:

                            res_number = residue.get_id()[1]
                            res_numbers.append(res_number)

            # Remove residues not present in the specified range
            pae_values_to_keep = remove_residues_not_in_range(pae_values, res_numbers)
            #print(pae_values_to_keep)

            model_element_list = [model_element]

            ref_residues = []
            ref_selectors = []

            # Get the coordinates of the residues.
            residues, coords_array, selectors = af_session.get_coords_array(model_element, "ca")

            ref_residues = residues
            #print(ref_residues)
            ref_selectors = selectors
            #print(ref_selectors)


            # Initializes the distance map window.
            cp = Contact_map_analysis_window_qt(self.main_window)
            cp.initialize_map(pymod=self,
                              data_array=pae_values_to_keep,
                              pymod_elements=model_element_list,
                              ref_residues=ref_residues,
                              ref_selectors=ref_selectors,
                              title="Predicted Aligned Error (PAE) for " + af_model,
                              pixel_size=5,
                              feature_type="pae",
                              threshold=35,
                              interaction_center="ca")
            cp.show()
            print("plotting pae")


def _highlight_in_pymol_from_plddt_plot_qt(point_data):

    if point_data is not None:

        cmd.select("pymod_selection", point_data["pymol_selector"])
        cmd.center("pymod_selection")


def remove_residues_not_in_range(pae_values, res_numbers):

    lower_bound, upper_bound = min(res_numbers), max(res_numbers)
    pae_lists_to_keep = [pae_list for i, pae_list in enumerate(pae_values) if i+1 >= lower_bound and i+1 <= upper_bound]
    #print(pae_lists_to_keep)
    pae_values_to_keep = [pae_list[lower_bound-1:upper_bound] for pae_list in pae_lists_to_keep]

    return pae_values_to_keep

##### END #####



# def show_af_plots(self, af_model, type):
#
#     if type == "plddt":
#
#         """
#         The predicted local distance difference test (pLDDT) score (0-100) is a per-residue confidence score,
#         with values greater than 90 indicating high confidence, and values below 50 indicating low confidence.
#         AlphaFold models are often shown with high confidence residues colored blue, and lower confidence in yellow,
#         orange and red. This measure estimates whether the predicted residue has similar distances to neighboring
#         C-alpha atoms (within 15 Angstroms) in agreement with distances in the true structure.
#         """
#
#         #print(dir(self))
#
#         print(af_model)
#         af_file_dirpath = self.af_modeling_session_list[self.current_af_session][af_model]["af_session_dirpath"] + "/" + af_model + ".pdb"
#         print(af_file_dirpath)
#         model_element = self.af_modeling_session_list[self.current_af_session][af_model]["model_element"]
#         print(model_element)
#
#         structure = Bio.PDB.PDBParser(QUIET=True).get_structure(af_model, af_file_dirpath)
#
#         pLDDT_data = []
#         res_numbers = []
#         res_names = []
#         pymol_selectors = []
#
#         for model in structure:
#
#             for chain in model:
#                 plot_name = af_model
#
#                 for residue in chain:
#                     #print(dir(residue))
#                     pymol_selector = f"{structure.get_id()} and resi {residue.get_id()[1]}"
#                     pymol_selectors.append(pymol_selector)
#
#                     if "CA" in residue:
#                         alpha_carbon = residue["CA"]
#                         b_factor = alpha_carbon.get_bfactor()
#                         pLDDT_data.append(b_factor)
#
#                         res_number = residue.get_id()[1]
#                         res_numbers.append(res_number)
#
#                         res_name = residue.get_resname()
#                         res_names.append(res_name)
#
#         # print(pLDDT_data)
#         # print(res_numbers)
#         # print(res_names)
#         #print(pymol_selectors)
#
#         residue_additional_data = []
#         for r_name, r_position, r_selector in zip(res_names, res_numbers, pymol_selectors):
#             residue_additional_data.append({"residue_name": r_name,
#                                                 "residue_pdb_position": r_position,
#                                                 "pymol_selector": r_selector,
#                                                 "export_label": "%s %s"%(r_name, r_position)})
#         #print(residue_additional_data)
#
#         plddt_scores_dict = {} # Items will contain DOPE scores of only polymer residues.
#         # Retain only the pLDDT values for residues of the chain (standard and modified residues).
#
#         plddt_scores = []
#         for res, score in zip(model_element.residues, pLDDT_data):
#             if res.is_polymer_residue():
#                 plddt_scores.append(score)
#         plddt_scores_dict.update({model_element: plddt_scores})
#
#         # Builds a list of all pLDDT values of the residues in the selection.
#         l_plddt = []
#
#         l_plddt.extend(plddt_scores_dict[model_element])
#
#         # Takes the min and max values among all the selected residues.
#         min_value = min(l_plddt)
#         max_value = max(l_plddt)
#         # An array with the equally sapced limits generated with the list above.
#         bins = np.array(np.linspace(min_value, max_value, num=10))
#
#         # An array with all the DOPE values of a single chain in the selection.
#         a_plddt = np.array(plddt_scores_dict[model_element])
#         # An array with the id of the bins where those values reside.
#         inds = np.digitize(a_plddt, bins)
#
#         plddt_items = []
#         for plddt_score, bin_id in zip(a_plddt, inds): # zip(ldope, inds):
#             plddt_items.append({"score": plddt_score, "interval": bin_id})
#         #print(plddt_items)
#         # Actually assigns the pLDDT score to the residues of the PyMod element.
#         for res, plddt_item in zip(model_element.get_polymer_residues(), plddt_items):
#             res.plddt_score = plddt_item
#         model_element._has_plddt_scores = True
#
#         self.main_window.color_element_by_plddt(model_element)
#
#
#         cp = PyMod_plot_window_qt(self.main_window)
#         cp.initialize(pymod=self, title="Predicted Local Distance Difference Test (pLDDT) per Position")
#         messagebar_text_on_update_qt = "Selected: __residue_name__ __residue_pdb_position__ of __plot_name__, pLDDT value: __y__"
#         cp.build_plotting_area(messagebar_initial_text="Set the 'Interact on click' option to 'Yes' to click on the plot and to highlight residues in PyMOL.",
#                                update_messagebar=True,
#                                messagebar_text_on_update=messagebar_text_on_update_qt,
#                                on_click_action=_highlight_in_pymol_from_plddt_plot_qt,
#                                x_label_text="Positions",
#                                y_label_text="pLDDT")
#
#         cp.add_plot(res_numbers, pLDDT_data, label = plot_name, additional_data = residue_additional_data)
#         cp.show()
#         print("plotting plddt")
#
#
#     elif type == "pae":
#
#         """
#         Predicted aligned error (PAE) gives a distance error for every pair of residues.
#         It gives AlphaFold's estimate of position error at residue x when the predicted and true structures are aligned on residue y.
#         Values range from 0 - 35 Angstroms. It is usually shown as a heatmap image with residue numbers running along vertical and
#         horizontal axes and color at each pixel indicating PAE value for the corresponding pair of residues.
#         If the relative position of two domains is confidently predicted then the PAE values will be low (less than 5A)
#         for pairs of residues with one residue in each domain.
#         """
#
#         if not check_network_connection("https://google.com", timeout=3):
#             title = "Connection Error"
#             message = ("An internet connection is not available, can not connect to the AFDB"
#                        " to download the corresponding PAE values.")
#             self.main_window.show_error_message(title, message)
#             return None
#
#         print(af_model)
#         af_file_dirpath = self.af_modeling_session_list[self.current_af_session][af_model]["af_session_dirpath"] + "/" + af_model + ".pdb"
#         print(af_file_dirpath)
#         uniprot_header = os.path.basename(os.path.normpath(self.af_modeling_session_list[self.current_af_session][af_model]["af_session_dirpath"])).rsplit('-', 2)[1]
#         print(uniprot_header)
#         model_element = self.af_modeling_session_list[self.current_af_session][af_model]["model_element"]
#         print(model_element)
#
#         # Link to AFDB
#         link_to_afdb_api = "https://alphafold.ebi.ac.uk/api/prediction/"
#
#         # Make full link with Uniprot-ID of interest
#         full_link = link_to_afdb_api + str(uniprot_header)
#         print(full_link)
#
#         # Open link and get data in JSON format
#         response = urllib.request.urlopen(full_link)
#         data = json.loads(response.read())
#
#         # Get url, from the JSON data, to download the document containing pae values
#         #TODO: SEE WHAT HAPPENS IF ENTRYid IS NOT FOUND
#         url_to_pae_doc = data[0]["paeDocUrl"]
#         print(url_to_pae_doc)
#
#         # Download PAE values corresponding to AlphaFold model in json format
#         response = requests.get(url_to_pae_doc)
#         pae_doc_dirpath = os.path.join(self.af_session_dirpath, "AF-" + uniprot_header + ".json")
#         pae_file = open(pae_doc_dirpath, "wb")
#         pae_file.write(response.content)
#         pae_file.close()
#
#         # To open the file containing PAE values in reading mode
#         with open(pae_doc_dirpath,'r') as file:
#             data = json.load(file)
#
#         # if the json file contains a list with one single element, that is a dictionary
#         if isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict):
#             # to obtain the dictionary 'pae_data' containing the values of PAE
#             pae_data = data[0]
#             # to get the values associated to the key 'predicted_aligned_error'
#             pae_values = pae_data.get("predicted_aligned_error")
#             #print(pae_values)
#
#         # if the json file directly contains a dictionary
#         else:
#             # to get the values associated to the key 'predicted_aligned_error'
#             pae_values = data["predicted_aligned_error"]
#
#
#         structure = Bio.PDB.PDBParser(QUIET=True).get_structure(af_model, af_file_dirpath)
#
#         res_numbers = []
#
#         for model in structure:
#
#             for chain in model:
#
#                 for residue in chain:
#
#                     if "CA" in residue:
#
#                         res_number = residue.get_id()[1]
#                         res_numbers.append(res_number)
#
#         # Remove residues not present in the specified range
#         pae_values_to_keep = remove_residues_not_in_range(pae_values, res_numbers)
#         #print(pae_values_to_keep)
#
#         model_element_list = [model_element]
#
#         ref_residues = []
#         ref_selectors = []
#
#         # Get the coordinates of the residues.
#         residues, coords_array, selectors = self.current_af_session.get_coords_array(model_element, "ca")
#
#         ref_residues = residues
#         #print(ref_residues)
#         ref_selectors = selectors
#         #print(ref_selectors)
#
#
#         # Initializes the distance map window.
#         cp = Contact_map_analysis_window_qt(self.main_window)
#         cp.initialize_map(pymod=self,
#                           data_array=pae_values_to_keep,
#                           pymod_elements=model_element_list,
#                           ref_residues=ref_residues,
#                           ref_selectors=ref_selectors,
#                           title="Predicted Aligned Error (PAE) for " + af_model,
#                           pixel_size=5,
#                           feature_type="pae",
#                           threshold=35,
#                           interaction_center="ca")
#         cp.show()
#         print("plotting pae")
