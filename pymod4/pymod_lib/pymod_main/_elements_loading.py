# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing a series of methods which allow to load PyMod elements in PyMod.
"""

import os
import re
import string

import shutil

from Bio import SeqIO

from pymol import cmd

from typing import Tuple, Dict, List, Set
from typing import Optional

from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt, askopenfile_qt
from pymod_lib.pymod_gui.specific_gui_components_qt import (Search_string_window_qt,
                                                            Add_feature_window_qt,
                                                            Edit_sequence_window_qt,
                                                            Rename_sequence_window_qt)

from pymod_lib import pymod_structure
from pymod_lib import pymod_vars
from pymod_lib.pymod_protocols import structural_databases_protocols
from pymod_lib.pymod_protocols.modeling_protocols.modeling_protocols_common import Fetch_models_file
from pymod_lib.pymod_seq.seq_manipulation import (check_correct_sequence,
                                                  get_invalid_typed_characters)
from pymod_lib.pymod_element_feature import Element_feature

from pymod_lib.pymod_exceptions import PyModInvalidFile


###################MODIFICATION On 22/05/2025-SARA###################################

#####WORKING ALSO ON MULTIPLE CHAIN HM#####################
def get_chains_IDs_pdb_file(pdb_file: str) -> List[str]:
    """
    Identify all protein chain IDs in a PDB file.
    
    Args:
        pdb_file: Path to the PDB file to analyze
        
    Returns:
        List of chain IDs (single characters) that contain standard amino acids
    """
    
    with open(pdb_file) as fin:
        lines = fin.readlines()

    # Set of standard 3-letter amino acid codes
    standard_aminoacids = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "SEC", "PYL"
    }

    # Step 1: Scan through all ATOM records to find protein chains
    protein_chains = []
    seen_residues = set()  # Track (chainID, resid) pairs to avoid duplicates
    
    for line in lines:
        if line.startswith("ATOM") and len(line) >= 26:
            resname = line[17:20].strip().upper()  # Extract residue name (columns 18-20)
            chainID = line[21]                     # Extract chain ID (column 22)
            resid = line[22:26].strip()           # Extract residue number (columns 23-26)
            
            # Check if this is a standard amino acid we haven't seen before
            if resname in standard_aminoacids and (chainID, resid) not in seen_residues:
                seen_residues.add((chainID, resid))
                if chainID not in protein_chains:  # Add chain ID if new
                    protein_chains.append(chainID)
                    
    return protein_chains

###################MODIFICATION On 26/06/2025-SARA###################################

def get_original_nucleic_chains(pdb_file: str) -> dict:
    """
    Identifies the original nucleic chains in the PDB file, both RNA and DNA,
    excluding mixed or non-nucleic chains.
    
    Args:
        pdb_file: Path to the PDB file
        
    Returns:
        Dictionary {chain_id: (num_residues, 'RNA' or 'DNA')} for pure nucleic chains only
    """
    with open(pdb_file) as fin:
        lines = fin.readlines()

    rna_resnames = {"A", "U", "G", "C", "RA", "RU", "RG", "RC"}
    dna_resnames = {"DA", "DT", "DG", "DC"}

    chain_res_types = {}    # chainID -> set of found types ("RNA", "DNA")
    chain_residues = {}     # chainID -> set of residues

    for line in lines:
        if len(line) < 27 or not line.startswith(("ATOM", "HETATM")):
            continue

        resname = line[17:20].strip().upper()
        chainID = line[21]
        resid = line[22:26].strip()

        if chainID not in chain_res_types:
            chain_res_types[chainID] = set()
            chain_residues[chainID] = set()

        if resname in rna_resnames:
            chain_res_types[chainID].add("RNA")
            chain_residues[chainID].add(resid)
        elif resname in dna_resnames:
            chain_res_types[chainID].add("DNA")
            chain_residues[chainID].add(resid)
        else:
            # Ignore non-nucleic residues
            pass

    result = {}
    # Consider only *pure* chains (only RNA or only DNA)
    protein_chains_ID = get_chains_IDs_pdb_file(pdb_file)

    for chainID, types in chain_res_types.items():
        if chainID not in protein_chains_ID:
            if len(types) == 1:  # only RNA or only DNA
                res_type = next(iter(types))  # "RNA" or "DNA"
                result[chainID] = (len(chain_residues[chainID]), res_type)

    return result

def merge_nucleic_acids_into_single_chain(pdb_file: str, chain_ID: str) -> None:
    """
    Merge all nucleic acids from different chains into a single specified protein chain in a PDB file.
    The result is written to a temporary file with a suffix based on the target chain.
    
    Args:
        pdb_file: Path to the input PDB file.
        chain_ID: The protein chain ID into which the nucleic acids will be merged.
    """

    output_file = f"{pdb_file}_temp_chain_{chain_ID}"

    # Read all lines from the input PDB file
    with open(pdb_file) as fin:
        lines = fin.readlines()

    last_resid = 0               # To track the last residue number in the target protein chain
    last_atom_serial = 0        # To track the last atom serial number in the target chain
    seen_residues = set()       # To avoid duplicate residue entries

    # First pass: Find the highest residue and atom serial numbers in the target chain
    for line in lines:
        if len(line) < 27:
            continue  # Skip lines too short to contain chain info

        if line.startswith(("ATOM", "HETATM", "TER")) and line[21] == chain_ID:
            try:
                serial = int(line[6:11])
                last_atom_serial = max(last_atom_serial, serial)
            except ValueError:
                continue

        if line.startswith("ATOM") and line[21] == chain_ID:
            try:
                resid = int(line[22:26])
                icode = line[26]
                key = (resid, icode)
                if key not in seen_residues:
                    seen_residues.add(key)
                    last_resid = max(last_resid, resid)
            except ValueError:
                continue

    # Define recognized nucleic acid residue names
    nucleic_resnames = {
    "A", "U", "G", "C", "T",        # singole lettere
    "DA", "DT", "DG", "DC",        # DNA
    "RA", "RU", "RG", "RC",         # RNA (a volte indicato così)
}


    new_hetatm_lines = []       # Will store new HETATM lines for merged nucleic acids
    resid_map = {}              # Map to assign new residue numbers to each unique nucleotide
    next_resid = last_resid + 1
    atom_serial = last_atom_serial + 1
    last_hetatm_line = None
    previous_chain = None

    # Second pass: Reassign nucleic acid atoms from other chains into the target chain
    for line in lines:
        if len(line) < 27:
            continue

        if line.startswith(("ATOM", "HETATM")) and line[21] != chain_ID:
            resname = line[17:20].strip().upper()
            if resname in nucleic_resnames:
                orig_chain = line[21]
                resid = line[22:26]
                icode = line[26]
                key = (orig_chain, resname, resid.strip(), icode)

                # Assign new residue number if not seen before
                if key not in resid_map:
                    resid_map[key] = next_resid
                    next_resid += 1
                new_resid = resid_map[key]

                # Insert TER line if the chain has changed
                if previous_chain and orig_chain != previous_chain:
                    ter_line = f"TER   {atom_serial:>5}      {resname} {chain_ID}{new_resid - 1:>4} \n"
                    new_hetatm_lines.append(ter_line)
                    atom_serial += 1

                # Create new HETATM line with updated serial, chain ID, and residue number
                newline = "HETATM" + line[6:]
                newline = f"HETATM{atom_serial:5d}" + newline[11:]
                newline = newline[:21] + chain_ID + newline[22:]
                newline = newline[:22] + f"{new_resid:>4}" + newline[26:]

                new_hetatm_lines.append(newline)
                last_hetatm_line = newline
                atom_serial += 1
                previous_chain = orig_chain

    # Final TER line after last HETATM block
    ter_line = ""
    if last_hetatm_line:
        resname = last_hetatm_line[17:20]
        chainID = last_hetatm_line[21]
        resSeq = last_hetatm_line[22:26]
        iCode = last_hetatm_line[26]
        ter_line = f"TER   {atom_serial:>5}      {resname} {chainID}{resSeq}{iCode}                                                      \n"
        atom_serial += 1

    # Write output file with merged content
    with open(output_file, "w") as fout:
        inserted = False
        for i, line in enumerate(lines):
            fout.write(line)

            # Insert HETATM block after the first TER line for the target chain
            if not inserted and line.startswith("TER"):
                previous_atom_is_chain = False
                for j in range(i - 1, -1, -1):
                    prev_line = lines[j]
                    if prev_line.startswith("ATOM") and prev_line[21] == chain_ID:
                        previous_atom_is_chain = True
                        break
                    elif prev_line.startswith("ATOM"):
                        break

                if previous_atom_is_chain:
                    for het_line in new_hetatm_lines:
                        fout.write(het_line)
                    if ter_line:
                        fout.write(ter_line)
                    inserted = True

################# MODIFICATION on 15/05/2026 - SARA #########################
def merge_nucleic_acids_into_protein_chains(pdb_file: str, output_file: Optional[str] = None, backup: bool = True) -> None:
    """
    Merges protein and nucleic acid chains together so MODELLER can read them.
    It combines the old working logic with new safety checks and clean prints.
    """
    if not os.path.isfile(pdb_file):
        raise FileNotFoundError(f"❌ File not found: '{pdb_file}'.")

    # Set default output file name if not provided
    output_file = output_file or pdb_file.replace(".pdb", "_merged.pdb")

    # Create a backup file if it already exists
    if backup and os.path.exists(output_file):
        shutil.copy(output_file, output_file + ".backup")
        print(f"📁 Backup of existing output created: {output_file}.backup")

    # --- START OF WORKING MERGE LOGIC ---
    with open(pdb_file) as fin:
        original_lines = fin.readlines()

    master_line = None
    for line in original_lines:
        if line.startswith("MASTER"):
            master_line = line
            break

    protein_chains_IDs = get_chains_IDs_pdb_file(pdb_file)
    
    # Get original nucleic acid info for the final summary print
    nucleic_chains_info = get_original_nucleic_chains(pdb_file)

    output_lines = []        
    used_lines_idx = set()   

    for chain_ID in protein_chains_IDs:
        chain_lines = []

        # Get all ATOM lines for the current protein chain
        for i, line in enumerate(original_lines):
            if len(line) > 21 and ((line.startswith("ATOM") and line[21] == chain_ID) or (line.startswith("TER") and line[21] == chain_ID)):
                chain_lines.append(line)
                used_lines_idx.add(i)

        # Create a temporary file to run the merge routine
        temp_pdb = f"{pdb_file}_temp_chain_{chain_ID}"
        merge_nucleic_acids_into_single_chain(pdb_file, chain_ID)

        with open(temp_pdb) as ftemp:
            modified_lines = ftemp.readlines()

        hetatm_start = False
        hetatm_block = []
        # Get the new HETATM lines created for this chain
        for line in modified_lines:
            if len(line) > 21 and line.startswith("HETATM") and line[21] == chain_ID:
                hetatm_start = True
                hetatm_block.append(line)
            elif hetatm_start and line.startswith("TER") and line[21] == chain_ID:
                hetatm_block.append(line)
                break

        # Save both protein lines and their merged nucleic lines
        output_lines.extend(chain_lines)
        output_lines.extend(hetatm_block)

        # Delete the temporary file
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)

    # Keep any other line that was not processed (like waters or other ligands)
    for i, line in enumerate(original_lines):
        if i not in used_lines_idx:
            if line.startswith(("ATOM", "HETATM", "TER")):
                output_lines.append(line)

    # Write the final PDB file for MODELLER
    with open(output_file, "w") as fout:
        fout.writelines(output_lines)
        if master_line:
            fout.write(master_line)
        fout.write("END\n")
    # --- END OF WORKING MERGE LOGIC ---

    # --- START OF SUMMARY PRINTS ---
    print(f"✅ Final merged file written to: {output_file}")
    print(f"✅ Protein chains processed: {protein_chains_IDs}")
    print(f"✅ Nucleic acid chains detected: {list(nucleic_chains_info.keys())}")
    print("📌 Merge successfully performed for MODELLER template preparation.")
    # --- END OF SUMMARY PRINTS ---

################ END ###################


class PyMod_elements_loading:

    ###############################################################################################
    # FILES MANAGMENT.                                                                            #
    ###############################################################################################

    #################################################################
    # Check correct files formats.                                  #
    #################################################################

    def is_sequence_file(self, file_path, file_format, show_error=True):
        """
        Try to open a sequence file using Biopython. Returns 'True' if the file is a valid file of
        the format specified in the 'file_format' argument.
        """
        valid_file = False
        file_handler = None
        try:
            file_handler = open(file_path,"r")
            r = list(SeqIO.parse(file_handler, file_format))
            if len(r) > 0:
                valid_file = True
        except:
            valid_file = False
        if file_handler != None:
            file_handler.close()
        return valid_file


    def is_valid_structure_file(self, file_name, format="pdb", show_error=True):
        """
        Validate structure files. Supports PDB, CIF, and other formats.
        """
        valid_structure = False
        
        try:
            file_handler = open(file_name, "r")
            content = file_handler.read()
            file_handler.close()
            
            if format.lower() == "pdb":
                # Check for PDB format - look for ATOM or HETATM records with coordinates
                for line in content.split('\n'):
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            if len(line) >= 54:
                                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                                valid_structure = True
                                break
                        except ValueError:
                            continue
                            
            elif format.lower() in ("cif", "mmcif"):
                # Check for CIF format - look for typical CIF headers and coordinate data
                if "_atom_site." in content or "data_" in content:
                    # Look for coordinate data in CIF format
                    lines = content.split('\n')
                    coordinate_section = False
                    for line in lines:
                        if "_atom_site." in line:
                            coordinate_section = True
                        elif coordinate_section and not line.startswith(('_', '#', 'loop_')):
                            # Try to parse a coordinate line
                            try:
                                parts = line.split()
                                if len(parts) >= 3:
                                    # Try to find coordinate columns (usually positions vary in CIF)
                                    coords_found = 0
                                    for part in parts:
                                        try:
                                            float(part)
                                            coords_found += 1
                                        except ValueError:
                                            continue
                                    if coords_found >= 3:
                                        valid_structure = True
                                        break
                            except:
                                continue
                                
            else:
                # For other formats, just check if file exists and is readable
                valid_structure = len(content.strip()) > 0
                
        except Exception as e:
            valid_structure = False
            
        if not valid_structure and show_error:
            title = "FileType Error"
            format_name = format.upper() if format else "structure"
            message = f"The selected file is not a valid {format_name} file."
            self.main_window.show_error_message(title, message)
            
        return valid_structure


    #################################################################
    # Load sequence files.                                          #
    #################################################################

    ###################MODIFICATION on 31/03/2025-SARA###################################
    """
    Distinguish different polymer types in the fasta files sequences 
    The identify_sequence_type function was added to classify sequences as DNA, RNA, or protein based on their composition. 
    The open_sequence_file method was updated to use this function, associating the correct polymer type with each sequence element.
    """
    def identify_sequence_type(self, sequence):
        """
        Identifies the type of sequence: DNA, RNA, or protein.
        """
        dna_bases = {'A', 'T', 'G', 'C'}
        rna_bases = {'A', 'U', 'G', 'C'}

        # Without the gaps and the residues which are not known: a sequence which was
        # already aligned when it was read was called a protein on the strength of its
        # gaps, whatever it was made of.
        sequence = "".join(c for c in sequence.upper() if c not in "-.X")
        if not sequence:
            return "protein"

        # Check if the sequence is DNA
        if all(base in dna_bases for base in sequence):
            return "dna"
        # Check if the sequence is RNA
        elif all(base in rna_bases for base in sequence):
            return "rna"
        # Otherwise, treat it as a protein
        else:
            return "protein"
    
    def open_sequence_file(self, file_full_path, file_format="fasta"):
        """
        Method for loading new sequences parsed from sequence files.
        It distinguishes between DNA, RNA, and protein sequences.
        """
        if not os.path.isfile(file_full_path):
            raise IOError("File does not exist: %s." % file_full_path)
        if not self.is_sequence_file(file_full_path, file_format):
            raise PyModInvalidFile("Cannot open an invalid '%s' file: %s." % (file_format, file_full_path))

        elements_to_return = []
        
        with open(file_full_path) as seq_handle:
            seq_records = list(SeqIO.parse(seq_handle, file_format))
        for record in seq_records:
            # Extract the sequence from the record
            sequence = str(record.seq)

            # Identify the sequence type
            sequence_type = self.identify_sequence_type(sequence)
            
            # Now associate polymer_type with the element
            element = self.build_pymod_element_from_seqrecord(record)
            element.polymer_type = sequence_type  # Associate polymer type with the element
            ###############MODIFICATION on 27/06/2025##############
            #######Return the file path################
            element.structure_file_path = file_full_path
            ##################END##############

            # Add the element to PyMod (this will now carry the polymer type metadata)
            e = self.add_element_to_pymod(element)
            elements_to_return.append(e)

            # Optionally print the polymer type for confirmation
            print(f"{sequence_type.upper()} sequence: {sequence}")

        return elements_to_return
    ###########END MODIFICATION#############

    def build_cluster_from_alignment_file(self, alignment_file, extension="fasta"):
        """
        Creates a cluster with all the sequences contained in an alignment file.
        """
        # Gets the sequences using Biopython.
        aligned_elements = []
        with open(alignment_file) as aln_handle:
            records = list(SeqIO.parse(aln_handle, extension))
        for record in records:
            new_child_element = self.build_pymod_element_from_seqrecord(record)
            self.add_element_to_pymod(new_child_element)
            aligned_elements.append(new_child_element)
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="alignment", child_elements=aligned_elements, algorithm="imported")
        return new_cluster

    

    #################################################################
    # Opening PDB files.                                            #
    #################################################################

   ##################### MODIFICATION ON 19/05/2026 - SARA ###########################
    def open_structure_file(self, pdb_file_full_path, file_format="pdb"):
        """
        Opens a PDB file, processes hybrid nucleic/protein entries, splits them into 
        individual chain files starting from atom serial 1, and strips out the 
        obsolete segment ID column to prevent parsing errors.
        Additionally, injects proper TER records between chain transitions in the main file.
        Force-syncs the finalized, clean content directly into __0_structure_full__.pdb.
        
        ULTRA-ACCURATE VERSION: Fixed pixel-perfect column alignment for HOH and residues.
        """

        # Check if the file format and path are valid before processing
        if not self.is_valid_structure_file(pdb_file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file." % file_format)

        # Extract directory path and filename components
        pdb_dir = os.path.dirname(pdb_file_full_path)
        pdb_filename = os.path.basename(pdb_file_full_path)

        # Sanitize filename by removing spaces, dashes, and invalid characters
        new_pdb_filename = pdb_filename.replace("-", "_").replace(" ", "_")
        new_pdb_filename = re.sub(r'[^\w\.]+', '_', new_pdb_filename)

        # Define the target path inside the local structures directory
        new_pdb_file_full_path = os.path.join(
            self.structures_dirpath,
            new_pdb_filename
        )

        # Copy the original file to the structures directory if it is not already there
        if os.path.abspath(pdb_file_full_path) != os.path.abspath(new_pdb_file_full_path):
            try:
                shutil.copy(pdb_file_full_path, new_pdb_file_full_path)
                pdb_file_full_path = new_pdb_file_full_path

            except OSError as e:
                print(f"Error copying file to structures directory: {e}")
                raise

        # An mmCIF becomes a PDB here, before anything reads it. What follows --
        # the search for nucleic acids, the splitting of the chains, the TER
        # records -- reads the text of the file by the columns of a PDB, and the
        # structure parser was being given the mmCIF as though it were one:
        # "invalid literal for int() with base 10: 'A'", the chain read where
        # the residue number should be.
        if file_format == "cif":
            converted_path = os.path.splitext(new_pdb_file_full_path)[0] + ".pdb"
            try:
                pymod_structure.convert_cif_to_pdb(new_pdb_file_full_path,
                                                   converted_path)
            except Exception as e:
                raise PyModInvalidFile(
                    "The mmCIF file could not be read: %s." % e)
            new_pdb_file_full_path = converted_path
            pdb_file_full_path = converted_path
            new_pdb_filename = os.path.basename(converted_path)

        # Helper internal function to scan for the presence of nucleic acids (DNA/RNA)
        def has_nucleic_acids(pdb_path):
            with open(pdb_path) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        resname = line[17:20].strip().upper()
                        if resname in {
                            "A", "U", "G", "C", "T",
                            "DA", "DT", "DG", "DC", "DU",
                            "RA", "RU", "RG", "RC"
                        }:
                            return True
            return False

        # If nucleic acids are detected, create a backup and merge them into the protein setup
        is_na_present = has_nucleic_acids(new_pdb_file_full_path)

        # Where the structure is kept aside untouched, before its nucleic acids are
        # merged into the protein chains. Named with splitext and not by stripping
        # ".pdb": rstrip takes a set of characters, not a suffix, so it ate any
        # trailing '.', 'p', 'd' or 'b' of the name itself -- "1abd.pdb" became "1a",
        # and two structures whose names differ only in those letters were kept aside
        # in the same file. "1aay.pdb" and the like are unaffected either way.
        original_copy_path = (
            os.path.splitext(new_pdb_file_full_path)[0] + ".original.pdb"
        )

        if is_na_present:
            print(">>> Nucleic acids detected. Preparing files in structures folder...")
            shutil.copy(new_pdb_file_full_path, original_copy_path)

            try:
                # Execute the external routine to merge nucleic chains into main chains
                merge_nucleic_acids_into_protein_chains(
                    pdb_file=original_copy_path,
                    output_file=new_pdb_file_full_path,
                    backup=False
                )
                print(
                    f"✅ PDB file modified for nucleotides directly in structures: "
                    f"{new_pdb_file_full_path}"
                )
                pdb_file_full_path = new_pdb_file_full_path

            except Exception as e:
                print(f"⚠️ Error while modifying the PDB file: {e}")

        # Determine which file version should be parsed initially by the application
        file_to_parse = new_pdb_file_full_path

        if is_na_present and os.path.exists(original_copy_path):
            file_to_parse = original_copy_path

        # Initialize the core parsed PDB object.
        # 'new_file_name' is what the elements are named after. It has to be the name
        # of the structure the user opened, not the name of the file being parsed:
        # for a structure with nucleic acids that file is the '.original.pdb' copy
        # kept aside above, and without this every chain would be called
        # "<name>.original_chain_A" instead of "<name>_chain_A".
        p = pymod_structure.Parsed_pdb_file(
            self,
            file_to_parse,
            output_directory=self.structures_dirpath,
            new_file_name=new_pdb_filename
        )

        # Store file paths and metadata inside the internal components
        dict_na = get_original_nucleic_chains(new_pdb_file_full_path)
        elements_to_return = []

        for element in p.get_pymod_elements():
            element.structure_file_path = new_pdb_file_full_path
            element.chain_len_na_info_dict = dict_na
            e = self.add_element_to_pymod(element)
            elements_to_return.append(e)

        # Handle naming conflicts if the root filename is already taken in the session
        original_pdb_filename = p.structure_file_name

        if original_pdb_filename in self.original_pdb_files_set:
            counter = 1
            while original_pdb_filename in self.original_pdb_files_set:
                original_pdb_filename = "%s_%s" % (
                    counter,
                    original_pdb_filename
                )
                counter += 1

            for e in elements_to_return:
                e.rename_file_name_root(original_pdb_filename)

        self.original_pdb_files_set.add(original_pdb_filename)

        # -------------------------------------------------------------------------------------
        # CHAIN FILES GENERATION
        # -------------------------------------------------------------------------------------
        # Extract the base name and gather all unique chain identifiers from the current file
        pdb_base_name = os.path.splitext(
            os.path.basename(new_pdb_file_full_path)
        )[0]

        if os.path.exists(new_pdb_file_full_path):
            with open(new_pdb_file_full_path, "r") as f_merged:
                merged_lines = f_merged.readlines()
        else:
            merged_lines = []

        all_chains = set()
        for line in merged_lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chain_char = line[21].strip()
                if chain_char:
                    all_chains.add(chain_char)

        # Process each chain individually to write its dedicated sub-file
        for target_chain in all_chains:
            custom_split_filename = (
                f"{pdb_base_name}.original_chain_{target_chain}.pdb"
            )
            custom_split_filepath = os.path.join(
                self.structures_dirpath,
                custom_split_filename
            )

            chain_specific_lines = []
            atom_serial = 1

            # Filter lines matching the target chain and reformat standard columns
            for line in merged_lines:
                if line.startswith("TER"):
                    continue

                if (
                    line.startswith(("ATOM", "HETATM"))
                    and len(line) > 21
                    and line[21] == target_chain
                ):
                    record_type = line[0:6]
                    atom_name = line[12:16]
                    alt_loc = line[16]
                    res_name = line[17:20]
                    chain_id = line[21]
                    res_seq = line[22:26]
                    i_code = line[26]
                    coords = line[27:54]
                    occupancy = line[54:60]
                    temp_factor = line[60:66]

                    elem_char = atom_name.strip()[0]
                    extra_fields = f"           {elem_char:<2}"

                    # Reconstruct the line ensuring exact column spacing
                    new_line = (
                        f"{record_type}{atom_serial:>5} "
                        f"{atom_name}{alt_loc}{res_name} "
                        f"{chain_id}{res_seq}{i_code}"
                        f"{coords}{occupancy}{temp_factor}"
                        f"{extra_fields:<14}\n"
                    )
                    chain_specific_lines.append(new_line)
                    atom_serial += 1

            # Append the required structural trailing records (TER, END) to the chain file
            if chain_specific_lines:
                last_serial = atom_serial
                last_line_parsed = chain_specific_lines[-1]
                last_res_name = last_line_parsed[17:20]
                last_res_seq = last_line_parsed[22:26]

                ter_line = "%-6s%5d      %3s %1s%4s\n" % (
                    "TER",
                    int(last_serial),
                    str(last_res_name).strip().rjust(3),
                    str(target_chain).strip(),
                    str(last_res_seq).strip().rjust(4)
                )

                chain_specific_lines.append(ter_line)
                chain_specific_lines.append("END\n")

                # Save the completed individual chain file
                with open(custom_split_filepath, "w") as f_out:
                    f_out.writelines(chain_specific_lines)

                print(
                    f"   -> Fixed columns and serials for chain "
                    f"{target_chain}: {custom_split_filepath}"
                )

                # Link the written file path to its corresponding internal element tracking variable
                for element in elements_to_return:
                    curr_id = None
                    # get_chain_id() is how a PyMod element reports its chain. Without
                    # it none of the attributes below ever matched, every element kept
                    # pointing at the whole structure, and a template built from it
                    # carried every chain of the complex instead of its own.
                    try:
                        curr_id = element.get_chain_id()
                    except Exception:
                        curr_id = None
                    if curr_id:
                        pass
                    elif hasattr(element, "chain") and element.chain:
                        curr_id = element.chain
                    elif hasattr(element, "chain_id") and element.chain_id:
                        curr_id = element.chain_id
                    elif (
                        hasattr(element, "chain_id_string")
                        and element.chain_id_string
                    ):
                        curr_id = element.chain_id_string

                    if curr_id and str(curr_id).strip() == target_chain:
                        element.structure_file_path = custom_split_filepath

        # -------------------------------------------------------------------------------------
        # MAIN HYBRID FILE POST-PROCESSING
        # -------------------------------------------------------------------------------------
        # Organize the layout of the consolidated file, separating standard and special lines
        if os.path.exists(new_pdb_file_full_path):
            with open(new_pdb_file_full_path, "r") as f_in:
                main_lines = f_in.readlines()

            # Separate standard macromolecular chains from special nucleic chains (e.g., X, Y)
            special_chains = ['X', 'Y']
            standard_chains = sorted(list(all_chains - set(special_chains)))

            # Catalogs to categorize macromolecule records and solvent records separately
            chain_macromol = {c: [] for c in all_chains}
            chain_water = {c: [] for c in all_chains}

            for line in main_lines:
                if line.startswith(("ATOM", "HETATM")):
                    chain_id = line[21]
                    res_name = line[17:20].strip().upper()
                    
                    # Separate water molecules (HOH) from standard polymeric residues
                    if res_name == "HOH":
                        if chain_id in chain_water:
                            chain_water[chain_id].append(line)
                    else:
                        if chain_id in chain_macromol:
                            chain_macromol[chain_id].append(line)

            formatted_main = []
            global_serial = 1

            # Internal helper to perform pixel-perfect column formatting for structural files
            def format_pdb_line(line, serial):
                record_type = line[0:6]
                atom_name = line[12:16]
                alt_loc = line[16]
                res_name = line[17:20]
                chain_id = line[21]
                res_seq = line[22:26]
                i_code = line[26]
                coords = line[27:54]
                occupancy = line[54:60]
                temp_factor = line[60:66]
                elem_char = atom_name.strip()[0]
                extra_fields = f"           {elem_char:<2}"
                
                return (
                    f"{record_type}{serial:>5} "
                    f"{atom_name}{alt_loc}{res_name} "
                    f"{chain_id}{res_seq}{i_code}"
                    f"{coords}{occupancy}{temp_factor}"
                    f"{extra_fields:<14}\n"
                )

            # Phase 1: Write macromolecules belonging to standard chains (Protein + Nucleic)
            for c in standard_chains:
                lines = chain_macromol[c]
                if not lines:
                    continue
                
                for line in lines:
                    formatted_main.append(format_pdb_line(line, global_serial))
                    global_serial += 1
                
                # Append a specific TER entry based on the final residue coordinates of the chain
                last_line = lines[-1]
                ter_line = "%-6s%5d      %3s %1s%4s\n" % (
                    "TER",
                    global_serial,
                    last_line[17:20].strip().rjust(3),
                    c,
                    last_line[22:26].strip().rjust(4)
                )
                formatted_main.append(ter_line)
                global_serial += 1
                
                # Add a blank TER line to separate standard macromolecular segments
                formatted_main.append("TER\n")

            # Phase 2: Write solvent blocks (water molecules) belonging to standard chains
            for c in standard_chains:
                w_lines = chain_water[c]
                if not w_lines:
                    continue
                
                for line in w_lines:
                    formatted_main.append(format_pdb_line(line, global_serial))
                    global_serial += 1
                
                # Add a trailing blank TER separator after each water block
                formatted_main.append("TER\n")

            # Phase 3: Write special terminal chains (X, Y) grouped with their respective waters
            for c in special_chains:
                # Write the core special macromolecule block
                lines = chain_macromol.get(c, [])
                if lines:
                    for line in lines:
                        formatted_main.append(format_pdb_line(line, global_serial))
                        global_serial += 1
                    
                    last_line = lines[-1]
                    ter_line = "%-6s%5d      %3s %1s%4s\n" % (
                        "TER",
                        global_serial,
                        last_line[17:20].strip().rjust(3),
                        c,
                        last_line[22:26].strip().rjust(4)
                    )
                    formatted_main.append(ter_line)
                    global_serial += 1
                    formatted_main.append("TER\n")

                # Write the water molecules directly connected to this special chain
                w_lines = chain_water.get(c, [])
                if w_lines:
                    for line in w_lines:
                        formatted_main.append(format_pdb_line(line, global_serial))
                        global_serial += 1
                    
                    formatted_main.append("TER\n")

            # Remove any trailing blank TER lines to clean up before the file closure record
            if formatted_main and formatted_main[-1] == "TER\n":
                formatted_main.pop()

            # Append the definitive structural END marker
            formatted_main.append("END\n")

            # Save the processed, newly structured content back into the main file
            with open(new_pdb_file_full_path, "w") as f_out:
                f_out.writelines(formatted_main)

            print(
                "   -> Main hybrid file finalized with custom layout: "
                f"{new_pdb_file_full_path}"
            )

        # -------------------------------------------------------------------------------------
        # FORCE-SYNC
        # -------------------------------------------------------------------------------------
        # Replicate the finalized, synchronized file structure across all active working directories
        possible_full_paths = [
            os.path.join(
                self.structures_dirpath,
                "__0_structure_full__.pdb"
            ),
            os.path.join(
                os.path.dirname(new_pdb_file_full_path),
                "__0_structure_full__.pdb"
            )
        ]

        if os.path.exists(new_pdb_file_full_path):
            try:
                with open(new_pdb_file_full_path, "r") as f_src:
                    clean_content = f_src.read()

                for full_path in possible_full_paths:
                    with open(full_path, "w") as f_dest:
                        f_dest.write(clean_content)
                    print(f"   ✅ Synchronized completely: {full_path}")

            except Exception as e:
                print(f"   ⚠️ Warning during __0_structure_full__.pdb force-sync: {e}")

        return elements_to_return
    ############################# END ##############################
    
    def color_struct(self):
        color_to_return = pymod_vars.pymod_regular_colors_list[self.color_index % len(pymod_vars.pymod_regular_colors_list)]
        self.color_index += 1
        return color_to_return


    #################################################################
    # Open files dialogs from PyMod.                                #
    #################################################################

    def choose_alignment_file(self):
        """
        Lets users choose an alignment file.
        """
        # Creates a PyQt widget that lets the user select multiple files.
        alignment_file_path = askopenfile_qt("Open an alignment file",
                                             name_filter="*fasta *fa *fas *aln *clu *sto *sth *phylip *phy *nexus *nex *gb *genbank *embl *swiss *pir",
                                             parent=self.get_qt_parent())
        if not alignment_file_path:
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(alignment_file_path)[1].replace(".","").lower()

        if extension in ("fasta", "fa", "fas"):
            extension = "fasta"
        elif extension in ("aln", "clu"):
            extension = "clustal"
        elif extension in ("sto", "sth"):
            extension = "stockholm"
        elif extension in ("phylip", "phy"):
            extension = "phylip"
        elif extension in ("nexus", "nex"):
            extension = "nexus"
        elif extension in ("gb", "genbank"):
            extension = "genbank"
        elif extension == "embl":
            extension = "embl"
        elif extension == "swiss":
            extension = "swiss"
        elif extension == "pir":
            extension = "pir"
        # Unknown format.
        else:
            title = "Format Error"
            message = "Unknown alignment file format: %s" % (extension)
            self.main_window.show_error_message(title, message)
            return (None, None)
        return alignment_file_path, extension


    def choose_structure_file(self):
        """
        Lets users choose a structure file.
        """
        # Creates a PyQt widget that lets the user select multiple files.
        open_filepath = askopenfile_qt("Open a structure file",
                                       name_filter="*pdb *ent *cif *mmcif *sdf *mol *mol2 *xyz",
                                       parent=self.get_qt_parent())
        if open_filepath == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(os.path.basename(open_filepath))[1].replace(".","").lower()
        return open_filepath, extension


    ###############################################################################################
    # EDIT SEQUENCE AND STRUCTURES.                                                               #
    ###############################################################################################

    def show_edit_sequence_window(self, pymod_element):
        """
        Edit a sequence.
        """
        self.edit_sequence_window = Edit_sequence_window_qt(self.main_window,
                                                            pymod_element=pymod_element,
                                                            title="Edit Sequence",
                                                            upper_frame_title="Edit your Sequence",
                                                            submit_command=self.edit_sequence_window_state)
        self.edit_sequence_window.show()
        # self.edit_sequence_window.resize(700, self.edit_sequence_window.sizeHint().height())


    def edit_sequence_window_state(self):
        """
        Accept the new sequence.
        """
        edited_sequence = self.edit_sequence_window.get_sequence()
        # When editing elements with a structure loaded in PyMOL, only indels can be added/removed.
        if self.edit_sequence_window.pymod_element.has_structure():
            if self.edit_sequence_window.pymod_element.my_sequence.replace("-", "") != edited_sequence.replace("-", ""):
                message = ("The amino acid sequence of an element with a 3D structure"
                           " loaded in PyMOL can not be edited. Only indels can be"
                           " added or removed.")
                self.main_window.show_error_message("Sequence Error", message)
                return None
        if not len(edited_sequence):
            self.main_window.show_error_message("Sequence Error", "Please submit a non empty string.")
            return None
        # The same alphabet "Add Raw Sequence" is typed in. What stood here allowed the
        # standard amino acids and nothing else, so the sequence of a DNA or RNA element
        # could not be edited at all: its uracils were refused as "not standard amino
        # acid characters".
        invalid_characters = get_invalid_typed_characters(edited_sequence)
        if invalid_characters:
            self.main_window.show_error_message(
                "Sequence Error",
                '%s %s not a residue. A sequence may hold the standard amino acid'
                ' letters, "X" for a residue which is not known, the bases of DNA and'
                ' RNA (A, T, C, G, U) and "-" for a gap.'
                % (", ".join('"%s"' % c for c in invalid_characters),
                   "is" if len(invalid_characters) == 1 else "are"))
            return None
        self.edit_sequence_window.pymod_element.set_sequence(edited_sequence, permissive=True)
        self.main_window.gridder(update_clusters=True, update_elements=True)
        self.edit_sequence_window.destroy()

    def show_rename_sequence_window(self, pymod_element):
        """
        Show dialog to rename a sequence element.
        """
        self.rename_sequence_window = Rename_sequence_window_qt(self.main_window,
                                                               pymod_element=pymod_element,
                                                               title="Rename Sequence",
                                                               upper_frame_title="Rename Sequence/Object",
                                                               submit_command=self.rename_sequence_window_state)
        self.rename_sequence_window.show()

    def rename_sequence_window_state(self):
        """
        Accept the new sequence name.
        """
        new_name = self.rename_sequence_window.get_new_name()
        
        # Validate the new name
        if not new_name:
            self.main_window.show_error_message("Name Error", "Please provide a non-empty name.")
            return None
            
        # Check for invalid characters (only letters, numbers, underscore allowed)
        import re
        if not re.match(r'^[a-zA-Z0-9_]+$', new_name):
            self.main_window.show_error_message("Name Error", "Invalid characters detected. Only letters, numbers and underscore (_) are allowed.")
            return None
            
        # Get the element
        element = self.rename_sequence_window.pymod_element
        old_header = element.my_header
        
        # Update the element header
        element.my_header = new_name
        
        # If element has a structure, also rename the PyMOL object and update ALL references
        if element.has_structure():
            old_pymol_selector = element.get_pymol_selector()
            
            # New PyMOL name is the same as new_name (already validated for allowed characters)
            new_pymol_name = new_name
            
            print(f"=== RENAME DEBUG START ===")
            print(f"Element: {element}")
            print(f"Old selector: '{old_pymol_selector}'")
            print(f"New name: '{new_pymol_name}'")
            
            try:
                # Rename in PyMOL first
                from pymol import cmd
                if old_pymol_selector in cmd.get_names("objects"):
                    cmd.set_name(old_pymol_selector, new_pymol_name)
                    print(f"✓ Renamed PyMOL object: '{old_pymol_selector}' → '{new_pymol_name}'")
                else:
                    print(f"⚠ PyMOL object '{old_pymol_selector}' not found in PyMOL")
                
                # Update ALL structure references in the element
                print(f"--- Updating structure references ---")
                
                # Set new pymol_selector attribute
                if not hasattr(element.structure, 'pymol_selector'):
                    element.structure.pymol_selector = None
                element.structure.pymol_selector = new_pymol_name
                print(f"✓ Set structure.pymol_selector = '{new_pymol_name}'")
                
                # Update chain_pdb_file_name
                if hasattr(element.structure, 'chain_pdb_file_name'):
                    old_file_name = element.structure.chain_pdb_file_name
                    new_file_name = old_file_name.replace(old_pymol_selector, new_pymol_name)
                    element.structure.chain_pdb_file_name = new_file_name
                    print(f"✓ Updated chain_pdb_file_name: '{old_file_name}' → '{new_file_name}'")
                
                # Update other selector references
                if hasattr(element.structure, 'chain_selector'):
                    element.structure.chain_selector = new_pymol_name
                    print(f"✓ Updated chain_selector: '{new_pymol_name}'")
                
                # Update element cache
                element._pymol_selector = new_pymol_name
                print(f"✓ Updated element._pymol_selector = '{new_pymol_name}'")
                
                # Update unique_index if needed
                if hasattr(element, 'unique_index') and old_pymol_selector in str(element.unique_index):
                    element.unique_index = element.unique_index.replace(old_pymol_selector, new_pymol_name)
                    print(f"✓ Updated unique_index")
                
                # Update initial_chain_file_path if it contains the old name
                if hasattr(element.structure, 'initial_chain_file_path'):
                    old_initial_path = element.structure.initial_chain_file_path
                    if old_pymol_selector in old_initial_path:
                        new_initial_path = old_initial_path.replace(old_pymol_selector, new_pymol_name)
                        element.structure.initial_chain_file_path = new_initial_path
                        print(f"✓ Updated initial_chain_file_path: '{old_initial_path}' → '{new_initial_path}'")
                
                # Test our get_pymol_selector method
                test_selector = element.get_pymol_selector()
                print(f"--- POST-RENAME TEST ---")
                print(f"element.get_pymol_selector() now returns: '{test_selector}'")
                
                if test_selector == new_pymol_name:
                    print(f"✅ SUCCESS: get_pymol_selector() returns correct new name")
                else:
                    print(f"❌ ERROR: get_pymol_selector() still returns old name!")
                    print(f"   Expected: '{new_pymol_name}'")
                    print(f"   Got: '{test_selector}'")
                
                print(f"=== RENAME DEBUG END ===")
                    
            except Exception as e:
                print(f"❌ Error during rename: {e}")
                import traceback
                traceback.print_exc()
        
        # Update PyMod interface
        self.main_window.gridder(update_elements=True)
        
        print(f"Renamed sequence: '{old_header}' → '{new_name}'")
        self.rename_sequence_window.destroy()


  
    def duplicate_sequence(self, element_to_duplicate):
        
        #Make a copy of a certain element.
        
        if element_to_duplicate.has_structure():
            p = pymod_structure.Parsed_pdb_file(self, element_to_duplicate.get_structure_file(basename_only=False),
                                                output_directory=self.structures_dirpath,
                                                new_file_name=pymod_vars.copied_chain_name % self.new_objects_index)
            self.new_objects_index += 1
            for element in p.get_pymod_elements():
                self.add_element_to_pymod(element, color=element_to_duplicate.my_color) # Add this to use the old color shceme of PyMod: color=self.color_struct()
        else:
            #duplicated_element = self.build_pymod_element(pmel.PyMod_sequence_element, element_to_duplicate.my_sequence, element_to_duplicate.my_header_root)
            duplicated_element = self.build_pymod_element_from_args(element_to_duplicate.my_header_root, element_to_duplicate.my_sequence)
            self.add_element_to_pymod(duplicated_element)
            return duplicated_element
    


    #################################################################
    # Add features.                                                 #
    #################################################################

    def show_add_feature_window(self, pymod_element, selected_residue):
        """
        Edit a feature to a residue (or series of residues).
        """

        self.add_feature_window = Add_feature_window_qt(self.main_window,
                                                        pymod_element=pymod_element,
                                                        selected_residue=selected_residue,
                                                        title="Add a Feature to %s" % pymod_element.compact_header,
                                                        upper_frame_title="Add a Feature to your Sequence",
                                                        submit_command=self.add_feature_window_state)
        self.add_feature_window.show()


    def add_feature_window_state(self):
        """
        Accept the new feature.
        """

        selected_element = self.add_feature_window.pymod_element

        #-----------------------------------
        # Get the parameters from the GUI. -
        #-----------------------------------

        # Gets the residue range.
        residue_range_str = self.add_feature_window.get_residue_range()

        selection_warning_title = "Selection Warning"
        examples_string = "'34' for a single residue, '12-20' for a residue range"
        selection_warning_message = "Please provide a valid string to select one ore more residues (use the format: %s)." % examples_string

        # Checks for an empty string.
        if not residue_range_str:
            self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
            return None

        # Checks for a single residue.
        sre_match = re.search("^\d+$", residue_range_str)
        if sre_match:
            try:
                residue_range = (int(residue_range_str), int(residue_range_str))
            except ValueError:
                self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
                return None

        # Checks for a residue range.
        else:
            sre_match = re.search("^(\d+\-\d+)$", residue_range_str)
            if not sre_match:
                self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
                return None

            try:
                min_res, max_res = residue_range_str.split("-")
                residue_range = (int(min_res), int(max_res))
            except (ValueError, IndexError):
                self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
                return None

            if residue_range[1] < residue_range[0]:
                self.main_window.show_warning_message(selection_warning_title,
                    "Invalid residue range. The index of the second residue in the range (%s) can not be smaller than the first one (%s)." % (residue_range[1], residue_range[0]))
                return None

        # Convert the 'db_index' values obtained from the GUI into 'seq_index' values.
        try:
            feature_start = selected_element.get_residue_by_db_index(residue_range[0]).seq_index
            feature_end = selected_element.get_residue_by_db_index(residue_range[1]).seq_index
         # The 'db_index' values do not correspond to any residue.
        except KeyError as e:
            self.main_window.show_warning_message(selection_warning_title,
                'The selected sequence does not have a residue with the following id: %s.' % e)
            return None


        # Gets the feature name.
        feature_name = self.add_feature_window.get_feature_name()
        sre_match = re.search("[^ a-zA-Z0-9_-]", feature_name)
        if not feature_name or sre_match:
            self.main_window.show_warning_message(selection_warning_title,
                'Please provide a valid "Feature Name" string (only alphanumeric characters and the "-", "_" and " " characters are allowed).')
            return None
            if len(feature_name) > 15:
                self.main_window.show_warning_message(selection_warning_title,
                    'Please provide a "Feature Name" string shorter than 15 characters.')
                return None


        # Gets the feature color.
        selected_rgb, selected_hex = self.add_feature_window.get_selected_colors()
        feature_color = self.add_new_color(selected_rgb, selected_hex)

        #------------------------------------------------
        # Actually adds the new feature to the element. -
        #------------------------------------------------

        # Adds the domain to the sequence.
        new_feature = Element_feature(id=None, name=feature_name,
                                      start=feature_start, end=feature_end,
                                      description=feature_name,
                                      feature_type='sequence',
                                      color=feature_color)
        selected_element.add_feature(new_feature)

        self.main_window.color_element_by_custom_scheme(selected_element, selected_feature=new_feature)

        if self.add_feature_window.get_select_in_pymol():
            selection_name = "%s_%s_%s" % (feature_name.replace(" ", "_"), selected_element.unique_index, selected_element.features_count)
            cmd.select(selection_name,
                       "object %s and resi %s-%s" % (selected_element.get_pymol_selector(), residue_range[0], residue_range[1]))
            selected_element.features_selectors_list.append(selection_name)

        self.add_feature_window.destroy()


    def delete_features_from_context_menu(self, pymod_element):
        for selection in pymod_element.features_selectors_list:
            try:
                cmd.delete(selection)
            except:
                pass
        pymod_element.clear_features()
        pymod_element.revert_original_color()
        self.main_window.color_element(pymod_element)


    #################################################################
    # Search subsequences.                                          #
    #################################################################

    def show_search_string_window(self, pymod_element):
        """
        Search for a sub-sequence in a protein sequence.
        """
        self.search_string_window = Search_string_window_qt(self.main_window,
                                                            pymod_elements=[pymod_element],
                                                            title="Search in %s" % pymod_element.compact_header,
                                                            upper_frame_title="Search a string in your sequence",
                                                            submit_command=self.search_string_window_state,
                                                            submit_button_text="Search")
        self.search_string_window.show()


    def search_string_window_state(self, event=None):

        selected_element = self.search_string_window.pymod_element

        #-----------------------------------
        # Get the parameters from the GUI. -
        #-----------------------------------

        # Get the string to search for in the selected sequence.
        search_string = self.search_string_window.get_search_string()
        _search_string = "" # " '" + search_string + "'"

        selection_warning_title = "Search Warning"
        # examples_string = "'34' for a single residue, '12-20' for a residue range"
        # selection_warning_message = "Please provide a valid string to select one ore more residues (use the format: %s)." % examples_string

        # Checks for an empty string.
        if not search_string:
            self.main_window.color_element(selected_element)
            self.search_string_window.show_results(self.search_string_window.default_message, state="empty")
            return None

        # Decides whether to use regex to search for a string.
        use_regex = self.search_string_window.get_regex_use()

        # Checks the search string.
        if not use_regex:
            # Search for invalid characters.
            search_sre = re.search("[^A-Za-z]", search_string)
            if search_sre:
                self.main_window.show_warning_message(selection_warning_title,
                    "Please provide a string containing only standard amino acid characters (example: 'ATGV').")
                return None

        else:
            # Test the regular expression.
            try:
                search_sre = re.search(search_string, "test_string")
            except re.error:
                self.main_window.show_warning_message(selection_warning_title,
                    "Please provide a string containing a valid regular expression.")
                return None

        # Gets the highlight color.
        highlight_color = self.search_string_window.get_highlight_color()


        #------------------------------------------------------------
        # Searches for the string in the selected protein sequence. -
        #------------------------------------------------------------

        selection_results_title = "Search Results"
        full_sequence = str(selected_element.my_sequence.replace("-", ""))

        # Actually executes the regular expression to look for sequences.
        finditer_results = list(re.finditer(search_string, full_sequence, flags=re.IGNORECASE))

        if not finditer_results:
            self.main_window.color_element(selected_element)
            results_message = "Pattern" + _search_string + " not found."
            self.search_string_window.show_results(results_message, state="not_found")
            return None

        # Builds a list of residue indices encompassing the matched strings.
        selected_residues = []
        for re_match in finditer_results:
            selected_residues.extend(list(range(re_match.start(), re_match.end())))
        selected_residues = set(selected_residues)

        # Stores the original colors and assign the temporary highlight color to the matched
        # substrings.
        original_color_scheme = selected_element.color_by
        for res in selected_element.get_polymer_residues():
            res.store_current_colors() # Stores the original colors.
            if res.seq_index in selected_residues:
                res.custom_color = highlight_color # Sets the highlight color.
            else:
                res.custom_color = res.get_default_color() # Use the original color.

        # Colors only the sequence and highlights the matching residues.
        self.main_window.color_element_by_custom_scheme(selected_element, use_features=False,
                                                        color_structure=False)

        # Restores the original colors. The next time the sequence is changed (for example
        # when inserting/removing gaps), the highlighted residues will be colored back with
        # their original color.
        for res in selected_element.get_polymer_residues():
            res.revert_original_colors()
        selected_element.color_by = original_color_scheme


        # Show a message with the summary of the results.
        def _get_results_string(n_results):
            if n_results == 1:
                return str(n_results) + " time."
            else:
                return str(n_results) + " times."

        results_message = ("Pattern" + _search_string + " found " + _get_results_string(len(finditer_results)))
        self.search_string_window.show_results(results_message, state="found")


        # Select in PyMOL.
        if self.search_string_window.get_select_in_pymol():
            residues_to_select_list = []
            for re_match in finditer_results:
                try:
                    pymol_start_res = selected_element.get_residue_by_index(re_match.start(), only_polymer=True).db_index
                    pymol_end_res = selected_element.get_residue_by_index(re_match.end(), only_polymer=True).db_index
                    residues_to_select_list.extend(list(range(pymol_start_res, pymol_end_res)))
                except IndexError:
                    pass
            if residues_to_select_list:
                selection_name = "pymod_search_string"
                selector_str = "object %s" % selected_element.get_pymol_selector()
                selector_str += " and resi " + self.main_window._join_residues_list(residues_to_select_list)
                cmd.select(selection_name, selector_str)


    #################################################################
    # Clusters.                                                     #
    #################################################################

    def update_cluster_sequences(self, cluster_element):
        """
        Updates the sequences of a cluster when some sequences are removed or added from the
        cluster.
        """
        children = cluster_element.get_children()
        if len(children) > 1:
            cluster_element.adjust_aligned_children_length()
            cluster_element.update_stars()
        else:
            if len(children) == 1:
                children[0].extract_to_upper_level()
            cluster_element.delete()


    def extract_selection_to_new_cluster(self):
        selected_sequences = self.get_selected_sequences()
        original_cluster_index = self.get_pymod_element_index_in_container(selected_sequences[0].mother) + 1
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="generic", child_elements=selected_sequences, algorithm="extracted")
        self.change_pymod_element_list_index(new_cluster, original_cluster_index)
        self.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)


    #################################################################
    # Transfer alignment files.                                     #
    #################################################################

    def transfer_alignment(self, alignment_element):
        """
        Changes the sequences of the elements contained in a PyMod cluster according to the
        information present in an externally supplied file (chosen by users through a file dialog)
        containing the same sequences aligned in a different way. Right now it supports transfer
        only for sequences having the exactly same sequences in PyMod and in the external alignment.
        """
        # Let users choose the external alignment file.
        openfilename, extension = self.choose_alignment_file()
        if None in (openfilename, extension):
            return False

        # Sequences in the aligment currently loaded into PyMod.
        aligned_elements = alignment_element.get_children()[:]

        # Sequences in the alignment files.
        external_records = list(SeqIO.parse(openfilename, extension))

        if len(external_records) < len(aligned_elements):
            title = "Transfer error"
            message = "'%s' has more sequences (%s) than the alignment in '%s' (%s) and the 'Transfer Alignment' function can't be used in this situation." % (alignment_element.my_header, len(aligned_elements), openfilename, len(external_records))
            self.main_window.show_error_message(title,message)
            return False

        correspondance_list = []
        # First try to find sequences that are identical (same sequence and same lenght) in both
        # alignments.
        for element in aligned_elements[:]:
            identity_matches = []
            for record in external_records:
                if str(element.my_sequence).replace("-","") == str(record.seq).replace("-",""):
                    match_dict = {"target-seq":element, "external-seq": record, "identity": True}
                    identity_matches.append(match_dict)
            if len(identity_matches) > 0:
                correspondance_list.append(identity_matches[0])
                aligned_elements.remove(identity_matches[0]["target-seq"])
                external_records.remove(identity_matches[0]["external-seq"])

        # Then try to find similar sequences among the two alignments. Right now this is not
        # implemented.
        # ...

        if not len(aligned_elements) == 0:
            title = "Transfer error"
            message = "Not every sequence in the target alignment has a corresponding sequence in the external alignment."
            self.main_window.show_error_message(title, message)
            return False

        # Finally transfer the sequences.
        for match in correspondance_list[:]:
            if match["identity"]:
                match["target-seq"].set_sequence(str(match["external-seq"].seq))
                correspondance_list.remove(match)

        self.main_window.gridder(update_clusters=True)


    def delete_cluster_dialog(self, cluster_element):
        title = "Delete Cluster?"
        message = "Are you sure you want to delete %s?" % (cluster_element.my_header)
        remove_cluster_choice = askyesno_qt(message=message, title=title, parent=self.get_qt_parent())
        if not remove_cluster_choice:
            return None
        title = "Delete Sequences?"
        message = "Would you like to delete all the sequences contained in the %s cluster? By selecting 'No', you will only extract them from the cluster." % (cluster_element.my_header)
        remove_children_choice = askyesno_qt(message=message, title=title, parent=self.get_qt_parent())

        # Delete both the cluster and its children.
        if remove_children_choice:
            cluster_element.delete()
        # Delete only the cluster element and extract the child sequences.
        else:
            children = cluster_element.get_children()
            for c in reversed(children[:]):
                c.extract_to_upper_level()
            cluster_element.delete()

        self.main_window.gridder(update_menus=True)


    #################################################################
    # Import PDB files.                                             #
    #################################################################

    def fetch_pdb_files(self, mode, target_selection):
        fp = structural_databases_protocols.Fetch_structure_file(self)
        fp.initialize_from_gui(mode, target_selection)
        fp.launch_from_gui()

    #################################################################
    # Import 3D-structures files of the models (both AlphaFold and ESMFold models) #
    #################################################################

    def fetch_models_files(self, mode, target_selection):
        fp = Fetch_models_file(self)
        fp.initialize_from_gui(mode, target_selection)

    # def associate_structure_from_popup_menu(self, target_element):
    #     """
    #     Launched when users press the 'Associate 3D Structure' from the leeft popup menu.
    #     """
    #     a = structural_databases_protocols.Associate_structure(self, target_element)
    #     a.launch_from_gui()
