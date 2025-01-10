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

from pymol import cmd

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.modeling_protocols.ESMFold._gui_qt import ESM_search_window_qt
from pymod_lib.pymod_protocols.structural_databases_protocols import Associate_structure
from pymod_lib.pymod_protocols.similarity_searches_protocols.local_blast import LOC_BLAST_search

###################################################################################################
# ESM database search.                                                                              #
###################################################################################################

class ESM_search(PyMod_protocol):

    """
    Class to represent a database search on the AlphaFold DB (AFDB).
    """

    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)

        self.esm_search_protocol_initialization()

    def esm_search_protocol_initialization(self):
        self.esm_search_modeling_window_class_qt = ESM_search_window_qt

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

        selected_sequences = self.pymod.get_selected_sequences()

        # First check if at least one sequence is selected.
        if not len(selected_sequences) > 0:
            self.pymod.main_window.show_error_message("Selection Error", "Please select at least one target sequence to use AFDB search tool.")
            return None

        # First check if at least one sequence is selected.
        if len(selected_sequences) > 1:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only one target sequence to use AFDB search tool.")
            return None

        # Checks if all the selected sequences can be used to build a model.
        if False in [s.can_be_modeled() for s in selected_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only sequences that do not have a structure loaded in PyMOL.")
            return None

        self.build_esm_search_window()


    def build_esm_search_window(self):
        """
        Builds the ESMFold search tool window.
        """

        self.esm_search_window = self.esm_search_modeling_window_class_qt(parent=self.pymod.main_window,
                                                             protocol=self)
        self.esm_search_window.show()


    def launch_esm_search(self):

        # Get selected sequences in PyMod
        selected_sequence = (self.pymod.get_selected_sequences())[0]

        # Get header of the selected sequences
        header = selected_sequence.my_sequence
        print(header)

        self.connect_to_esm(header, "prova")

        self.esm_search_window.close()


    def connect_to_esm(self, sequence, name = None):

        ABS_PATH = os.path.abspath("./")
        """Predict protein structure with ESMFold
        Args:
            sequence (str): amino acid sequence
            name (str, optional): _description_. Defaults to None.
        """
        sequence = re.sub("[^A-Z:]", "", sequence.replace("/", ":").upper())
        sequence = re.sub(":+", ":", sequence)
        sequence = re.sub("^[:]+", "", sequence)
        sequence = re.sub("[:]+$", "", sequence)

        headers = {
            "Content-Type": "application/x-www-form-urlencoded",
        }

        response = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/", headers=headers, data=sequence
        )
        if not name:
            name = sequence[:3] + sequence[-3:]
        pdb_filename = os.path.join(ABS_PATH, name) + ".pdb"
        pdb_string = response.content.decode("utf-8")
        if pdb_string.startswith("HEADER"):
            with open(pdb_filename, "w") as out:
                out.write(pdb_string)
            print(f"Results saved to {pdb_filename}")
            plddt = self.cal_plddt(pdb_string)
            print("="*40)
            print("    pLDDT: "+"{:.2f}".format(plddt))
            print("="*40)
            cmd.load(pdb_filename)
        else:
            print(pdb_string)


    def cal_plddt(self, pdb_string: str):
        """read b-factors of ca
        Args:
            pdb_string (str): _description_
        """

        lines = pdb_string.split("\n")
        plddts = []
        for line in lines:
            if " CA " in line:
                plddt = float(line[60:66])
                plddts.append(plddt)
        if max(plddts) <= 1.0:
            plddts =[ plddt * 100 for plddt in plddts]
            print("Guessing the scale is [0,1], we scale it to [0, 100]")
        else:
            print("Guessing the scale is [0,100]")
        return sum(plddts) / len(plddts)
