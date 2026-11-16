# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for predicting the structure of a sequence with ESMFold in PyMod.

Two ways in, both predictions. Without a key the sequence is folded by the public
ESMAtlas service; with a key it is folded through the ESM API, which also does the
sequence design. Neither looks a structure up in a database of precomputed models:
the only thing fetched is the sequence itself, from UniProt, when the header of the
selection looks like an accession.
"""

#AF

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

from pymod_lib.pymod_seq.seq_conservation import get_conservation_symbol


from pymod_lib.pymod_protocols.alignment_protocols.clustalw import Clustalw_regular_alignment, Clustalw_profile_alignment
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol, AFDB_common
from pymod_lib.pymod_protocols.modeling_protocols.AlphaFold.__gui_qt import AFDB_search_window_qt
from pymod_lib.pymod_protocols.structural_databases_protocols import Associate_structure
from pymod_lib.pymod_protocols.similarity_searches_protocols.local_blast import LOC_BLAST_search
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import Generic_BLAST_search, Similarity_searches_results_window_qt
from pymod_lib.pymod_protocols.modeling_protocols.modeling_protocols_common import Fetch_models_file
from pymod_lib.pymod_protocols.structural_analysis_protocols.structural_divergence_plot import Structural_divergence_plot
from pymod_lib.pymod_protocols.structural_analysis_protocols.contact_map_analysis import Contact_map_analysis_window_qt

from pymod_lib.pymod_seq import seq_manipulation

from pymod_lib.pymod_plot.pymod_plot_qt import PyMod_plot_window_qt
from pymod_lib.pymod_residue import PyMod_residue

##### MODIFIED #####
from pymod_lib.pymod_os_specific import check_network_connection
##### END #####

##### MODIFIED #####
from pymol.Qt import QtWidgets, QtCore
##### END #####

#AF


import os
import sys
import shutil
import importlib
import pickle
import re
import urllib.request
import json
import requests
from Bio import SeqIO
import Bio.PDB
import numpy as np


from pymol import cmd

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.modeling_protocols.ESMFold._gui_qt import ESM_search_window_qt
from pymod_lib.pymod_protocols.structural_databases_protocols import Associate_structure
from pymod_lib.pymod_protocols.similarity_searches_protocols.local_blast import LOC_BLAST_search

from pymod_lib.pymod_os_specific import check_network_connection
from Bio.PDB import PDBParser
from pymod_lib.pymod_protocols.modeling_protocols.modeling_protocols_common import Fetch_models_file
import pymod_lib.pymod_structure as pmstr
import pymod_lib.pymod_vars as pmdt
from pymod_lib.pymod_plot.pymod_plot_qt import PyMod_plot_window_qt

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


###################################################################################################
# ESMFold structure prediction.                                                                   #
###################################################################################################

# How long to wait for the keyless ESMAtlas service to fold a sequence. It is the
# folding which takes the time: see the note at the request below.
ESMATLAS_TIMEOUT = 300


class ESM_search(PyMod_protocol):

    """
    Predicts the structure of the selected sequence with ESMFold.

    The name says search for historical reasons and the class is still called that so
    that nothing which refers to it has to change, but nothing is searched: the
    sequence is sent to ESMFold and folded.
    """

    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)

        self.esm_search_protocol_initialization()
        #############MODIFICATION on 02/07/2025##########
        self.selected_sequences = self.pymod.get_selected_sequences()
        ##############END####################
        #############MODIFICATION on 10/07/2025##########
        #to give a name to ESMFOld model generated
        self.model_filename= None
        ##############END####################
        #############MODIFICATION on 17/07/2025##########
        self.model_output = []
        #################END########################
        #############MODIFICATION on 23/07/2025##########
        self.models_count = 0
        #################END########################
        
    
    def esm_search_protocol_initialization(self):
        self.esm_search_modeling_window_class_qt = ESM_search_window_qt

    def launch_from_gui(self):
        """
        Called by the entry which predicts a structure without a key, under
        "Tools -> Structure Modeling -> ESMFold".
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
            self.pymod.main_window.show_error_message("Selection Error", "Please select at least one target sequence to use ESMFold search tool.")
            return None

        ####################MODIFIED on 18/02/2025######################
        # Check if any selected sequence is a nucleic acid (DNA/RNA).
        # If so, display an error message and prevent the ESMFold search tool from proceeding.
        if any(e.polymer_type in ["dna", "rna"] for e in selected_sequences):
            self.pymod.main_window.show_error_message(
                "Selection Error", 
                "Cannot use ESMFold search tool for nucleic acids (RNA/DNA)."
            )
            return None
        #####################END###############################


        # First check if only one sequence is selected. ESMFold can work only on single chain
        if len(selected_sequences) > 1:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only one target sequence.\nESMFold supports single-chain inputs only")
            return None
        
        ################MODIFICATION on 30/07/2025#################
        ##########WARNING and block when a sequence is greater than 400 AA

        for e in selected_sequences:
            #print(dir(e))
            seq = e.my_sequence  
            if len(seq) > 400:
                self.pymod.main_window.show_error_message(
                    "Selection Error",
                    f"Sequence too long: {len(seq)}\nMaximum sequence length supported by ESMFold is 400 aa."
                )
                return None
        ###################END####################


        # Checks if all the selected sequences can be used to build a model.
        if False in [s.can_be_modeled() for s in selected_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only sequences that do not have a structure loaded in PyMOL.")
            return None
        
        ###############MODIFICATION on 29/07/2025#################
        ########CHECK FOR THE INTERNET CONNECTION IN GENERAL appears as soon as you click on build model button######
        ##### MODIFIED #####
        #Keep it here?? and the server connection error after???
        """
        if not check_network_connection("https://google.com", timeout=3):
            title = "Connection Error"
            message = ("An internet connection is not available.")
            self.esm_search_window.close()
            self.pymod.main_window.show_error_message(title, message)
            

            return None
        """ 
        ##### END #####
        ##################END#################

        self.build_esm_search_window()


    def build_esm_search_window(self):
        """
        Builds the ESMFold search tool window.
        """

        self.esm_search_window = self.esm_search_modeling_window_class_qt(parent=self.pymod.main_window,
                                                             protocol=self)
        self.esm_search_window.show()

        
    def launch_esm_search(self):
        """Main method to launch ESMFold search and modeling."""
        
        ###################MODIFIED on 02/07/2025################################
        # Note: The sequence needs to be one single chain, for multiple chains use AlphaFold instead
        ###################MODIFIED on 011/07/2025################################
        # Setup directory structure for model output
        
        # Base directory path for storing models
        models_dir = os.path.join(self.pymod.current_project_dirpath, self.pymod.models_dirname)
        
        # Check for UniProt ID in sequence header
        self.check_existing_uniprot_id()
        
        # Get UniProt ID if available
        uniprot_id = self.uniprot_header if self.uniprot_header else None
        
        # Generate session index (incremental counter)
        session_index = self.pymod.performed_esm_modeling_count + 1
        
        # Create session-specific directory name
        if uniprot_id:
            model_subdir_name = f"ESM-{uniprot_id}-session"
        else:
            model_subdir_name = f"ESM-session{session_index}"
        
        # Full output directory path
        self.esm_output_dir_path = os.path.join(models_dir, model_subdir_name)
        
        # Create directory if it doesn't exist
        if not os.path.isdir(self.esm_output_dir_path):
            os.mkdir(self.esm_output_dir_path)
        #################END##################
        
        ######OLD CODE (commented out)
        #output_dir=self.esm_output_dir_path
        ################END################

        #################MODIFICATION on 23/07/2025################
        # Get the selected sequence (only first one since ESMFold is single-chain)
        selected_sequence = (self.pymod.get_selected_sequences())[0]

        # Get sequence string from selected element
        header = selected_sequence.my_sequence

        # Check if sequence comes from a file
        file_path = selected_sequence.structure_file_path

        if file_path is not None:
            # If file exists, use its name (without extension)
            filename = os.path.basename(file_path)
            filename_no_ext = os.path.splitext(filename)[0]
        else:
            filename = None

            # Generate incremental model name → ESM_model_1, ESM_model_2, ...
            # Initialize counter if it doesn't exist
            if not hasattr(self.pymod, "performed_esm_modeling_count"):
                self.pymod.performed_esm_modeling_count = 0

            # Increment counter and create name
            self.pymod.performed_esm_modeling_count += 1
            filename_no_ext = f"ESM_model_{self.pymod.performed_esm_modeling_count}"

        # Debug output
        print(f"file_path: {file_path}")
        print(f"filename + ext: {filename}")
        print(f"filename (no ext): {filename_no_ext}")
        #################END################     

        ######################################MODIFIED on 02/07/2025################################
        # Connect to ESMFold API with sequence and filename
        self.connect_to_esm(header, filename_no_ext)
        
        # Close search window after submission
        self.esm_search_window.close()
        
        # Track current session
        self.current_esm_session = self
        self.pymod.esm_modeling_session_list.append(self)
        
        # Debug output
        print(f"ESM sessions: {len(self.pymod.esm_modeling_session_list)}")
        print(f"[LAUNCH] self.pymod id: {id(self.pymod)}")
        print(f"[LAUNCH] ESM session list: {self.pymod.esm_modeling_session_list}")
        
        # Update GUI display
        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True, update_elements=True)
        
    ###################MODIFIED on 02/07/2025################################
    ##########Function got from AFDB init################
    def check_existing_uniprot_id(self):

        # Get header of the selected sequences
        for seq in self.selected_sequences:
            header = self.pymod.get_header(seq)

        # Search for Uniprot-ID pattern matching.
        # If NO valid Uniprot-ID is found in the sequence header, self.uniprot_header = None
        uniprot_id_search = re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", header, re.IGNORECASE)


        ##### MODIFIED #####
        
        if not check_network_connection("https://google.com", timeout=3):
            title = "Connection Error"
            message = ("An internet connection is not available.")
            self.esm_search_window.close()
            self.pymod.main_window.show_error_message(title, message)

            return None
    
        ##### END #####

        if not uniprot_id_search:
            self.uniprot_header = None
        else:
            self.uniprot_header = uniprot_id_search.group(0)
            path = os.path.join(self.pymod.temp_directory_dirpath, self.uniprot_header + ".fasta")
            fastaurl = "https://rest.uniprot.org/uniprotkb/" + self.uniprot_header + ".fasta"
            response = requests.get(fastaurl)
            file = open(path, "wb")
            file.write(response.content)
            file.close()
            with open(path) as fasta_handle:
                fasta_records = list(SeqIO.parse(fasta_handle, "fasta"))
            for record in fasta_records:
                # Then builds a PyMod_element object and add it to the 'pymod_elements_list'.
                self.new_model_element = self.pymod.build_pymod_element_from_seqrecord(record)
    ###################END##########################
            print(f"uniprot ID: {self.uniprot_header}")
            print(f"new_model: {self.new_model_element}")
        

    def obtain_pdb_string(self, sequence):
        """
        The coordinates for this sequence as PDB text, or None when they could not
        be obtained, in which case the reason has already been shown to the user.

        Here the sequence is folded by the keyless ESMAtlas service. A protocol
        which obtains a structure some other way overrides this one method and
        inherits everything which happens to the coordinates afterwards.
        """

        ##############MODIFICATION on 29/07/2025##############
        #ESMFold server connection error
        # Attempt to send the amino acid sequence to the ESMFold API for structure prediction.
        # Handles various possible errors including timeouts, connection failures,
        # server-side errors, and unexpected exceptions by showing an appropriate
        # error message to the user and safely closing the related search window.

        try:
            response = requests.post(
                # To simulate a DNS failure, you could use a fake URL like:
                # "https://api.esmatlas.fake/foldSequence/v1/pdb/",
                
                "https://api.esmatlas.com/foldSequence/v1/pdb/",
                
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                
                # To simulate other server errors, you might send invalid data like:
                # data = "INVALID_SEQUENCE",
                
                data=sequence,
                
                # To simulate a timeout error, you could set a very low timeout like:
                # timeout=0.001,
                
                # Folding is what takes the time here, not the connection. Measured
                # against the service: a 21 residue peptide it had never seen took
                # 20 seconds, and a 39 residue one was still running when the server
                # itself gave up at 30. Three seconds meant almost every sequence
                # which was not already cached on the server was reported as a
                # connection error, so the keyless prediction looked broken.
                timeout=ESMATLAS_TIMEOUT
            )
            
            # Raises an HTTPError if the response status code is 4xx or 5xx
            response.raise_for_status()
            
            # Decode the response content from bytes to a UTF-8 string
            pdb_string = response.content.decode("utf-8")

        except requests.exceptions.Timeout as e:
            # Handles cases where the request times out (exceeds the timeout limit)
            title = "Timeout Error"
            message = ("The ESMFold server did not answer within %s seconds.\n\n"
                       "The keyless service folds the sequence while PyMod waits, so"
                       " a long sequence can take a while and the server sometimes"
                       " gives up on its own. Trying again often works, since it"
                       " keeps what it has already folded. For sequences which"
                       " repeatedly fail, the prediction with an API key does not go"
                       " through this service at all."
                       "\n\n%s" % (ESMATLAS_TIMEOUT, e))
            
            # Close the ESMFold search window because the process cannot continue
            self.esm_search_window.close()
            
            # Show an error message to the user indicating a timeout happened
            self.pymod.main_window.show_error_message(title, message)
            
            # Exit the function early due to the error
            return None

        except requests.exceptions.ConnectionError as e:
            # Handles connection-related errors, such as DNS failures or no internet connection
            title = "Connection Error"
            message = "Failed to connect to ESMFold server.\nPlease check your internet connection or try again later.\n\n" + str(e)
            
            self.esm_search_window.close()
            self.pymod.main_window.show_error_message(title, message)
            return None

        except requests.exceptions.HTTPError as e:
            # The server answered, and the answer was an error. Which error decides
            # whether there is anything the user can do about it, so they are not all
            # reported in the same way: a 5xx is the service failing, and no amount of
            # correcting the input will help.
            status = getattr(e.response, "status_code", None)

            if status in (500, 502, 503, 504):
                title = "The ESMFold service is unavailable"
                message = ("The ESMAtlas server answered %s, which means the service"
                           " behind it is not folding sequences at the moment.\n\n"
                           "This is not a problem with your sequence or with PyMod:"
                           " the same answer comes back for a sequence of any length,"
                           " including ones the server has folded before. The service"
                           " is free and unattended, and it goes down from time to"
                           " time.\n\n"
                           "What you can do:\n"
                           "  - try again later, since these outages are usually"
                           " temporary;\n"
                           "  - predict the structure with an API key instead, through"
                           " 'Predict with ESMFold2 (API key)', which does not go"
                           " through this service;\n"
                           "  - or use Boltz-2, which predicts complexes as well;\n"
                           "  - or, if the sequence is in UniProt, fetch the model"
                           " which was computed already with AlphaFold ->"
                           " Fetch from database, which costs nothing.\n\n%s" % (status, e))

            elif status == 413 or (status == 400 and "length" in str(e).lower()):
                title = "The sequence is too long for the keyless service"
                message = ("The ESMAtlas server refused the sequence with %s.\n\n"
                           "The keyless service takes one protein chain of at most 400"
                           " residues. For a longer one, use the prediction with an API"
                           " key, or Boltz-2.\n\n%s" % (status, e))

            else:
                title = "ESMFold Server Error"
                message = ("The ESMFold server refused the request%s.\n\n%s"
                           % (" with %s" % status if status else "", e))

            self.esm_search_window.close()
            self.pymod.main_window.show_error_message(title, message)
            return None

        except requests.exceptions.RequestException as e:
            # Anything else which went wrong in the request itself.
            title = "ESMFold Server Error"
            message = "An error occurred while communicating with the ESMFold server:\n\n" + str(e)
            
            self.esm_search_window.close()
            self.pymod.main_window.show_error_message(title, message)
            return None

        except Exception as e:
            # Catches any other unexpected errors that are not related to the HTTP request
            title = "Unexpected Error"
            message = f"An unexpected error occurred during ESMFold execution:\n\n{str(e)}"
            
            self.esm_search_window.close()
            self.pymod.main_window.show_error_message(title, message)
            return None
        
        ############END#################

        return pdb_string


    def connect_to_esm(self, sequence, name=None):
        """
        Connects to ESMFold API (not database) for structure prediction.
        Processes sequence and naming conventions before submitting to ESMFold.
        """
        
        ########### MODIFICATION on 5/07/2025 ###############
        # 1. Build filename consistent with chain information
        chain_name = ""
        if name:
            match = re.search(r'chain_[A-Za-z0-9]+', name)
            if match:
                chain_name = match.group()  # Example: "chain_A"

        # 2. Clean sequence - remove invalid characters and standardize format
        sequence = re.sub("[^A-Z:]", "", sequence.replace("/", ":").upper())
        sequence = re.sub(":+", ":", sequence)  # Collapse multiple colons
        sequence = re.sub("^[:]+", "", sequence)  # Remove leading colons
        sequence = re.sub("[:]+$", "", sequence)  # Remove trailing colons

        # Generate default name from sequence if not provided
        if not name:
            name = sequence[:3] + sequence[-3:]  # First 3 + last 3 residues

        # Check for UniProt ID in sequence header
        self.check_existing_uniprot_id()

        # Extract chain identifier from name (looking for '_chain_X' pattern)
        chain_match = re.search(r"_chain_([A-Za-z0-9])", name)
        if chain_match:
            chain = chain_match.group(1)
        else:
            # Fallback to chain_name attribute or default to 'A'
            chain = getattr(self, "chain_name", None) or "A"

        # Use UniProt ID as base if available, otherwise use "model"
        header_base = self.uniprot_header if self.uniprot_header else "model"

        # Add chain suffix if not already present in header
        if not header_base.endswith(f"_chain_{chain}"):
            header = f"{header_base}_chain_{chain}"
        else:
            header = header_base

        # Get current session index for unique naming
        if hasattr(self.pymod, "esm_modeling_session_list"):
            session_index = len(self.pymod.esm_modeling_session_list) + 1
        else:
            session_index = 1

        # Create output filenames with session index prefix
        self.model_filename = f"{session_index}_ESM-{header}.pdb"
        self.model_name = f"{session_index}_ESM-{header}"
        self.full_output_path = os.path.join(self.esm_output_dir_path, self.model_filename)

        # 4. ESMFold API call (old version commented out)
        """OLD VERSION:
        # 4. ESMFold API call
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
        response = requests.post("https://api.esmatlas.com/foldSequence/v1/pdb/", 
                            headers=headers, data=sequence)
        pdb_string = response.content.decode("utf-8")
        """

        # The coordinates, from wherever this protocol takes them.
        pdb_string = self.obtain_pdb_string(sequence)
        if pdb_string is None:
            return None

    
        # 5. Validate and save the PDB file
        if not pdb_string.startswith("HEADER"):
            print("Error in ESMFold response:\n", pdb_string)
            return None

        with open(self.full_output_path, "w") as out:
            out.write(pdb_string)
        print(f"[✓] ESM model saved to {self.full_output_path}")

        # 6. Calculate pLDDT score
        plddt = self.cal_plddt(pdb_string)
        print("=" * 40)
        print("    pLDDT: " + "{:.2f}".format(plddt))
        print("=" * 40)


        #####################MODIFICATION on 17/07/2025#############################
        
        # 9. Add model info to self.model_output (similar to MODELLER output)
        try:
            # Create dictionary with current model info
            current_model_info = {
                'name': self.model_filename,
                'num': session_index,
                'DOPE score': None,     # Not calculated by ESMFold
                'pLDDT': round(plddt, 2),
            }

            # Initialize model_output list if it doesn't exist
            if not hasattr(self, "model_output") or self.model_output is None:
                self.model_output = []
                
            # Add model info if not already present
            if current_model_info not in self.model_output:
                self.model_output.append(current_model_info)
            print(f"[✓] model_output updated with: {current_model_info}")
        except Exception as e:
            print(f"[X] Error updating model_output: {e}")


        #####################MODIFICATION on 23/07/2025###############
        ###### Loading into PyMod and PyMOL #############
        model_output = self.model_output
        # Terminate session if no valid models were generated
        if not model_output:
            error_message = ("No models could be correctly predicted by ESMFold."
                            " Terminating the modeling session.")
            self.finish_modeling_session(successful=False, error_message=error_message)
            return None

        # Get color scheme for visualization
        colors = pmdt.pymol_light_colors_list

        # Process the generated model. ESMFold produces a single, single-chain
        # model, so we take the first (and only) entry of 'model_output'.
        model = model_output[0]

        # Color index for this model: the number of previously completed ESM
        # sessions (the current session has not been appended to the session
        # list yet, that happens back in 'launch_esm_search').
        model_file_number = len(self.pymod.esm_modeling_session_list)

        # The sequence that was folded, used below to align the parsed model
        # back to the original target sequence.
        try:
            original_target_aligned_sequence = self.selected_sequences[0].my_sequence.replace("-", "")
        except Exception:
            original_target_aligned_sequence = None

        model_file_name = model['name']  # e.g.: '1_ESM-P12235_chain_A.pdb'
        model_file_full_path = os.path.join(self.esm_output_dir_path, model_file_name)

        ################MODIFICATION on 24/07/2025##################
        # Clean model name by removing duplicate chain identifiers
        filename_no_ext = os.path.splitext(model_file_name)[0]

        # Remove any "_chain_X" suffix where X is a letter or number
        filename_no_ext = re.sub(r'_chain_[A-Za-z0-9]+$', '', filename_no_ext)

        pymod_model_name = filename_no_ext
        #######################END###########################

        # Determine if we're processing multiple chains
        self.multiple_chain_mode = len(self.pymod.esm_modeling_session_list) > 1

        # Set color scheme - single color for single chain, multiple colors for complexes
        if not self.multiple_chain_mode:
            self.single_chain_model_color = colors[model_file_number % len(colors)]
        else:
            self.single_chain_model_color = "white"
        self.list_of_model_chains_colors = colors
        # Parse the PDB model file
        self.parsed_model_file = pmstr.Parsed_model_pdb_file(
            self.pymod, model_file_full_path,
            output_directory=self.pymod.structures_dirpath,
            new_file_name=pymod_model_name,
            model_root_name=pymod_model_name  # keeps the exact name
        )

        parsed_model_file = self.parsed_model_file
        current_model_chains_elements = []

        # First, verify pymod_model_name is correctly defined
        pymod_elements = parsed_model_file.get_pymod_elements()

        print(f"[DEBUG] len(pymod_elements): {len(pymod_elements)}")
        # Process each chain in the model
        for chain_number, model_element in enumerate(pymod_elements):
            # Determine color for this chain
            color = (self.list_of_model_chains_colors[
                    (model_file_number + chain_number) % len(self.list_of_model_chains_colors)]
                    if self.multiple_chain_mode else self.single_chain_model_color)
            # Add element to PyMod with assigned color
            self.pymod.add_element_to_pymod(model_element, color=color)
            current_model_chains_elements.append(model_element)

            # Track sequence alignment if available. For a complete de novo
            # ESMFold model this aligns the parsed model back to the folded
            # target sequence; skip it if the target sequence is unavailable.
            if original_target_aligned_sequence and hasattr(model_element, 'my_sequence'):
                model_element.trackback_sequence(original_target_aligned_sequence)

        # Store the loaded model element on the session so the per-model plots
        # (e.g. the pLDDT plot in 'show_esm_plots') can access it later.
        self.new_model_element = (current_model_chains_elements[0]
                                  if current_model_chains_elements else None)

        return self.full_output_path

    def get_model_number(self):
        """
        Returns the current model count and logs debug information.
        
        Returns:
            int: The number of models generated in this session
        """
        return self.models_count


    def get_model_color(self, chain_number, multiple_chain_mode):
        """
        Determines the appropriate color for a model chain based on visualization mode.
        
        Args:
            chain_number (int): Index of the chain to color
            multiple_chain_mode (bool): Whether we're displaying multiple chains
            
        Returns:
            str: The color to use for this chain in visualization
        """
        # For single chain mode, use the predefined single color
        if not multiple_chain_mode:
            color = self.single_chain_model_color
        else:
            # For multi-chain mode, cycle through available colors
            color = self.list_of_model_chains_colors[
                chain_number % len(self.list_of_model_chains_colors)
            ]
            
        return color
    ###########################END#######################


    def cal_plddt(self, pdb_string: str) -> float:
        """
        Calculates the average pLDDT score from a PDB file string.
        
        Extracts B-factors (storing pLDDT scores) for all C-alpha atoms and returns
        the average score. Handles both [0,1] and [0,100] pLDDT scales automatically.

        Args:
            pdb_string (str): The PDB file content as a string

        Returns:
            float: The average pLDDT score (0-100 scale)

        Notes:
            - pLDDT (predicted Local Distance Difference Test) is a per-residue 
              confidence score where:
              * >90 = high confidence
              * 70-90 = good confidence  
              * 50-70 = low confidence
              * <50 = very low confidence
            - The B-factor column (columns 61-66) in PDB files stores pLDDT scores
        """
        lines = pdb_string.split("\n")
        plddts = []
        
        # Extract pLDDT scores from all C-alpha atoms
        for line in lines:
            if " CA " in line:  # C-alpha atom line
                plddt = float(line[60:66])  # Extract B-factor (pLDDT score)
                plddts.append(plddt)
        
        # Handle different pLDDT scales (0-1 or 0-100)
        if max(plddts) <= 1.0:
            plddts = [plddt * 100 for plddt in plddts]  # Convert 0-1 to 0-100 scale
            print("Detected pLDDT scores in [0,1] range - converted to [0,100] scale")
        else:
            print("Detected pLDDT scores in [0,100] range - using as-is")
        
        # Calculate and return average pLDDT
        return sum(plddts) / len(plddts)


def show_esm_plots(self, esm_session, type):
    """
    Visualizes ESMFold model quality metrics in interactive plots.
    
    Args:
        esm_session: The ESM modeling session containing the models
        type (str): Type of plot to generate (currently only 'plddt' supported)
    
    Currently supports:
    - pLDDT (per-residue confidence scores) visualization
    - Interactive residue selection in PyMOL
    """
    
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
        cp.build_plotting_area(
            messagebar_initial_text="Set the 'Interact on click' option to 'Yes' to click on the plot and to highlight residues in PyMOL.",
            update_messagebar=True,
            messagebar_text_on_update=messagebar_text_on_update_qt,
            on_click_action=_highlight_in_pymol_from_plddt_plot_qt,
            x_label_text="Positions",
            y_label_text="pLDDT"
        )

        
        for esm_session in self.esm_modeling_session_list:
            esm_model = esm_session.model_name
            esm_file_path = esm_session.full_output_path

            # All chain elements of the (possibly multi-chain) model, in chain
            # order. Forge fold sets 'new_model_elements'; legacy / keyless ESM
            # sessions only have a single 'new_model_element'.
            chain_elements = getattr(esm_session, "new_model_elements", None)
            if not chain_elements:
                single = getattr(esm_session, "new_model_element", None)
                chain_elements = [single] if single is not None else []
            chain_elements = [e for e in chain_elements if e is not None]
            if not chain_elements:
                continue

            # The PDB chains, in file order, paired one-to-one with the PyMod
            # chain elements (loaded in the same order).
            structure = Bio.PDB.PDBParser(QUIET=True).get_structure(esm_model, esm_file_path)
            pdb_chains = [chain for model in structure for chain in model]

            multichain = len(chain_elements) > 1

            for chain_index, model_element in enumerate(chain_elements):
                # Legacy / keyless ESM elements may not have a structure yet;
                # only then associate one (single chain "A"). Forge fold models
                # are already loaded with a correct structure, so we must NOT
                # overwrite it (doing so corrupts get_pymol_selector()).
                if not model_element.has_structure():
                    new_structure = pmstr.PyMod_structure(
                        file_name_root=model_element,
                        full_file_path=esm_file_path,
                        chain_file_path=esm_file_path,
                        chain_id="A",
                        original_structure_file_path=esm_file_path,
                        original_structure_id=ESM_search.get_model_number)
                    model_element.set_structure(new_structure)
                    new_structure.initial_chain_file_path = esm_file_path
                    new_structure.pymod_element = model_element

                if chain_index >= len(pdb_chains):
                    continue
                chain = pdb_chains[chain_index]
                objsel = model_element.get_pymol_selector()

                # Per-residue pLDDT (CA B-factor) for this chain.
                pLDDT_data = []
                res_numbers = []
                res_names = []
                pymol_selectors = []
                for residue in chain:
                    if "CA" not in residue:
                        continue
                    pLDDT_data.append(residue["CA"].get_bfactor())
                    res_number = residue.get_id()[1]
                    res_numbers.append(res_number)
                    res_names.append(residue.get_resname())
                    pymol_selectors.append("%s and resi %s" % (objsel, res_number))

                if not pLDDT_data:
                    continue

                residue_additional_data = [
                    {"residue_name": r_name,
                     "residue_pdb_position": r_position,
                     "pymol_selector": r_selector,
                     "export_label": "%s %s" % (r_name, r_position)}
                    for r_name, r_position, r_selector
                    in zip(res_names, res_numbers, pymol_selectors)]

                # Assign the pLDDT score to this chain element's residues.
                a_plddt = np.array(pLDDT_data)
                bins = np.linspace(a_plddt.min(), a_plddt.max(), num=10)
                inds = np.digitize(a_plddt, bins)
                plddt_items = [{"score": score, "interval": bin_id}
                               for score, bin_id in zip(a_plddt, inds)]
                for res, plddt_item in zip(model_element.get_polymer_residues(),
                                           plddt_items):
                    res.plddt_score = plddt_item
                model_element._has_plddt_scores = True

                # Only build widgets if the element does not already have them
                # (legacy / keyless case). Forge fold models are already
                # displayed, and rebuilding their widget group would orphan the
                # existing GUI row.
                if not getattr(model_element, "widget_group", None):
                    self.main_window.add_pymod_element_widgets(model_element)
                self.main_window.color_element_by_plddt(model_element)

                plot_label = (esm_model if not multichain
                              else "%s (chain %s)" % (esm_model, chain.id))
                cp.add_plot(res_numbers, pLDDT_data, label=plot_label,
                            additional_data=residue_additional_data)

        cp.show()
        print("plotting plddt")

    elif type == "pae":
        """
        Predicted Aligned Error (PAE) gives, for every pair of residues, the
        expected position error at residue x when the predicted and true
        structures are aligned on residue y. It is shown as an N x N heatmap.
        Unlike AlphaFold (which downloads PAE from the AFDB), the Forge fold
        service returns the PAE matrix directly, stored on the session as
        'pae_matrix'.
        """
        pae_values = getattr(esm_session, "pae_matrix", None)
        if not pae_values:
            self.main_window.show_error_message(
                "PAE not available",
                "This model does not have Predicted Aligned Error data.")
            return None

        # For a complex the PAE matrix spans every chain (in chain order), so
        # gather the reference residues/selectors across all chain elements.
        model_elements = getattr(esm_session, "new_model_elements", None)
        if not model_elements:
            single = getattr(esm_session, "new_model_element", None)
            model_elements = [single] if single is not None else []
        model_elements = [e for e in model_elements if e is not None]
        if not model_elements:
            self.main_window.show_error_message(
                "PAE not available",
                "The model element for this session is not available.")
            return None

        esm_model = getattr(esm_session, "model_name", "ESM model")

        # Reference residues/selectors for mapping plot clicks back to PyMOL,
        # concatenated over all chains in chain order. These are built directly
        # from each element's polymer residues (one entry per residue) so they
        # line up exactly with the PAE matrix, which is N x N over all residues
        # of the complex. Using 'get_coords_array' here is unsafe: it drops any
        # residue whose CA is not currently found in PyMOL, which would make
        # 'ref_residues' shorter than the matrix and break the plot's hover.
        ref_residues = []
        ref_selectors = []
        for model_element in model_elements:
            objsel = model_element.get_pymol_selector()
            for res in model_element.get_polymer_residues():
                ref_residues.append(res)
                ref_selectors.append("object %s and n. CA and i. %s"
                                     % (objsel, res.db_index))

        # As a final safety net, keep the matrix and the reference lists the
        # same length (the plot grid is sized from the matrix).
        n = min(len(pae_values), len(ref_residues))
        if n != len(pae_values):
            pae_values = [row[:n] for row in pae_values[:n]]
            ref_residues = ref_residues[:n]
            ref_selectors = ref_selectors[:n]

        cp = Contact_map_analysis_window_qt(self.main_window)
        cp.initialize_map(pymod=self,
                          data_array=pae_values,
                          pymod_elements=model_elements,
                          ref_residues=ref_residues,
                          ref_selectors=ref_selectors,
                          title="Predicted Aligned Error (PAE) for " + esm_model,
                          pixel_size=5,
                          feature_type="pae",
                          threshold=35,
                          interaction_center="ca")
        cp.show()
        print("plotting pae")


def _highlight_in_pymol_from_plddt_plot_qt(point_data):

    #Data for the option interact on click
    #ex. point data: {'residue_name': 'GLN', 'residue_pdb_position': 13, 'pymol_selector': '1_ESM-Q96DW6_chain_A and resi 13', 'export_label': 'GLN 13'}
    if point_data is not None:

        cmd.select("pymod_selection", point_data["pymol_selector"])
        cmd.center("pymod_selection")
