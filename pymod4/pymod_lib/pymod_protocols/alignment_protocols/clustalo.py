# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.



####################MODIFICATION code Sara 08/06/2026############################
# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Clustal Omega alignment protocol module for PyMod.

================================================================================
DEBUG & FIX NOTES (Profile and Sequence-to-Profile Alignments):
================================================================================
- PROBLEM:
  Clustal Omega's profile alignment mode (--profile1/2) uses HHalign, which 
  fails if profiles are 100% identical or lack internal gaps. Injecting fake 
  gaps at the end breaks HHalign matrix logic. Additionally, PyMod sometimes 
  passes unaligned .fasta files instead of .aln files to these methods.

- SOLUTION:
  We check if the input profiles contain actual alignment gaps ('-'). If they 
  do not, it means the sequences inside are already aligned/identical. In this 
  case, instead of forcing a broken profile alignment, we dynamically merge the 
  files into a single fasta file and use Clustal Omega's regular alignment mode 
  (--infile). This completely avoids the strict HHalign validator, prevents 
  FATAL ERRORs, and stops downstream Segmentation Faults.
================================================================================
"""

import os

# Protocols.
from ._clustal_common import Clustal_regular_alignment, Clustal_profile_alignment

# GUI.
from ._base_alignment._gui import Regular_alignment_window_qt, Profile_alignment_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_entryfield_qt, PyMod_radioselect_qt


###################################################################################################
# GUI COMPONENTS                                                                                  #
###################################################################################################

class Clustalo_base_window_qt:
    """
    GUI window for Clustal Omega options.
    """
    def build_algorithm_options_widgets(self):
        # Option: Use full distance matrix
        self.use_full_dm_rds = PyMod_radioselect_qt(label_text="Use Full Distance Matrix",
                                                    buttons=('Yes', 'No'))
        self.use_full_dm_rds.setvalue("No")
        self.middle_formlayout.add_widget_to_align(self.use_full_dm_rds)

        # Option: Number of iterations
        self.clustalo_iterations_enf = PyMod_entryfield_qt(label_text="Combined Iterations",
                                                           value='0',
                                                           validate={'validator': 'integer',
                                                                     'min': 0, 'max': 5})
        self.middle_formlayout.add_widget_to_align(self.clustalo_iterations_enf)
        self.middle_formlayout.set_input_widgets_width("auto")

    def get_iterations_value(self):
        return self.clustalo_iterations_enf.getvalue(validate=True)

    def get_use_full_dm_value(self):
        return self.use_full_dm_rds.getvalue() == "Yes"


class Clustalo_regular_window_qt(Clustalo_base_window_qt, Regular_alignment_window_qt):
    pass


class Clustalo_profile_window_qt(Clustalo_base_window_qt, Profile_alignment_window_qt):
    pass


###################################################################################################
# ALIGNMENT PROTOCOLS                                                                             #
###################################################################################################

class Clustalomega_alignment:
    """
    Base class for Clustal Omega.
    """
    alignment_program = "clustalo"
    protocol_name = "clustalo"

    def additional_initialization(self):
        self.tool = self.pymod.clustalo

    def get_options_from_gui(self):
        self.params_from_gui = {}
        error_from_gui = False

        self.params_from_gui["use_full_dm"] = self.alignment_window.get_use_full_dm_value()

        try:
            self.params_from_gui["iterations"] = self.alignment_window.get_iterations_value()
        except Exception as e:
            error_from_gui = True
            error_message = str(e)

        if error_from_gui:
            self.pymod.main_window.show_error_message("Parameters Error", error_message)
            return False
        else:
            return True


class Clustalomega_regular_alignment(Clustalomega_alignment, Clustal_regular_alignment):
    """
    Standard multiple sequence alignment.
    """
    def get_alignment_window_class_qt(self):
        return Clustalo_regular_window_qt

    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_clustalo(sequences_to_align,
                          output_file_name=output_file_name,
                          extraoption="",
                          iterations=self.params_from_gui["iterations"],
                          use_full_dm=self.params_from_gui["use_full_dm"])

    def run_clustalo(self, sequences_to_align, output_file_name=None, extraoption="", iterations=0, use_full_dm=False):
        print("\n[CLUSTALO DEBUG] Starting regular alignment...")
        self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)

        input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
        output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
        guidetree_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".dnd")

        clustalo_exec = self.tool["exe_file_path"].get_value()
        command = [clustalo_exec,
                   "--infile", input_file_path,
                   "--outfile", output_file_path,
                   "--guidetree-out", guidetree_file_path,
                   "--force", "--outfmt", "clustal",
                   "--resno", "--output-order", "tree-order"]

        if any(e.polymer_type in ["dna", "rna"] for e in self.selected_elements):
            if any(e.polymer_type == "dna" for e in self.selected_elements):
                seqtype_flag = "--seqtype=dna"
            elif any(e.polymer_type == "rna" for e in self.selected_elements):
                seqtype_flag = "--seqtype=rna"
            command.append(seqtype_flag)

        if iterations != 0:
            command.extend(["--iter", str(iterations)])
        if use_full_dm:
            command.extend(["--full", "--full-iter"])
        
        print(f"[CLUSTALO DEBUG] Command: {' '.join(command)}")
        self.pymod.execute_subprocess(' '.join(command))


class Clustalomega_profile_alignment(Clustalomega_alignment, Clustal_profile_alignment):
    """
    Profile-to-Profile and Sequence-to-Profile alignment.
    """
    def get_alignment_window_class_qt(self):
        return Clustalo_profile_window_qt

    def _file_has_gaps(self, filepath):
        """ Returns True if the file contains alignment gaps ('-'). """
        if not os.path.exists(filepath):
            return False
        try:
            with open(filepath, 'r') as f:
                content = f.read()
            return '-' in content
        except:
            return False

    def _convert_to_clean_fasta_content(self, filepath):
        """ Reads an input file (.aln or .fasta) and returns a clean list of lines in FASTA format. """
        if not os.path.exists(filepath):
            return ""
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # If it's already a fasta file, return its content
            if any(line.startswith('>') for line in lines):
                return "".join(lines)
            
            # If it's a clustal file, convert blocks to simple fasta format strings
            sequences = {}
            for line in lines:
                if line.strip() and not line.startswith("CLUSTAL") and not line.strip().startswith(("*", ":", ".")):
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_id = parts[0]
                        seq_data = parts[1].replace('-', '') # Remove any gap for safety
                        if seq_id not in sequences:
                            sequences[seq_id] = ""
                        sequences[seq_id] += seq_data
            
            fasta_str = ""
            for seq_id, seq_data in sequences.items():
                fasta_str += f">{seq_id}\n{seq_data}\n"
            return fasta_str
        except:
            return ""

    def _build_fallback_regular_command(self, file1, file2, output_file_shortcut):
        """ Fuses two gapless or problematic inputs into one single standard ClustalOmega run. """
        print(f"[CLUSTALO FALLBACK] Converting inputs into a unified standard alignment run to bypass HHalign bug...")
        
        content1 = self._convert_to_clean_fasta_content(file1)
        content2 = self._convert_to_clean_fasta_content(file2)
        
        fallback_fasta_path = output_file_shortcut + "_merged.fasta"
        with open(fallback_fasta_path, 'w') as f:
            f.write(content1 + "\n" + content2)
            
        clustalo_path = self.tool["exe_file_path"].get_value()
        cline = [clustalo_path,
                 "--infile=" + fallback_fasta_path,
                 "--outfile=" + output_file_shortcut + ".aln",
                 "--outfmt=clustal", "--force",
                 "--resno", "--output-order", "tree-order"]
                 
        if self.params_from_gui.get("iterations", 0) != 0:
            cline.extend(["--iter", str(self.params_from_gui["iterations"])])
        if self.params_from_gui.get("use_full_dm", False):
            cline.extend(["--full", "--full-iter"])
            
        return " ".join(cline)

    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        # FALLBACK: If the profile has no real gaps, HHalign will crash. Fallback to standard multiple alignment mode.
        if not self._file_has_gaps(profile_file_shortcut):
            return self._build_fallback_regular_command(profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut)

        clustalo_path = self.tool["exe_file_path"].get_value()
        cline = [clustalo_path,
                 "--profile1=" + profile_file_shortcut,
                 "--outfile=" + output_file_shortcut + ".aln",
                 "--outfmt=clustal", "--force",
                 "--resno", "--output-order", "tree-order"]

        polymer_types = set()
        elements_source = getattr(self, 'elements_to_add', []) or getattr(self, 'selected_elements', [])

        for element in elements_source:
            if hasattr(element, 'get_children'):
                for child in element.get_children():
                    if hasattr(child, 'polymer_type'): polymer_types.add(child.polymer_type)
            elif hasattr(element, 'polymer_type'):
                polymer_types.add(element.polymer_type)
        if not polymer_types:
            polymer_types.add("protein")

        if "dna" in polymer_types or "rna" in polymer_types:
            if "dna" in polymer_types and "rna" not in polymer_types:
                seqtype_flag = "--seqtype=dna"
            elif "rna" in polymer_types and "dna" not in polymer_types:
                seqtype_flag = "--seqtype=rna"
            else:  
                self.pymod.main_window.show_error_message("Selection Error", "Cannot mix DNA and RNA.")
                return None
            cline.append(seqtype_flag)
        
        if hasattr(self, 'elements_to_add') and len(self.elements_to_add) > 1:
            cline.extend(['--infile', sequences_to_add_file_shortcut])
        else:
            cline.extend(['--profile2', sequences_to_add_file_shortcut])
    
        if self.params_from_gui.get("iterations", 0) != 0:
            cline.extend(["--iter", str(self.params_from_gui["iterations"])])
        if self.params_from_gui.get("use_full_dm", False):
            cline.extend(["--full", "--full-iter"])

        return ' '.join(cline)

    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        # FALLBACK: If either profile lacks real gaps, HHalign fails. Fallback to standard multiple alignment mode.
        if not self._file_has_gaps(profile1) or not self._file_has_gaps(profile2):
            return self._build_fallback_regular_command(profile1, profile2, output_file_shortcut)

        clustalo_path = self.tool["exe_file_path"].get_value()
        cline = [clustalo_path,
                 "--profile1=" + profile1,
                 "--profile2=" + profile2,
                 "--outfile=" + output_file_shortcut + ".aln",
                 "--outfmt=clustal", "--force"]

        polymer_types = set()
        if hasattr(self, 'selected_elements') and self.selected_elements:
            for element in self.selected_elements:
                if hasattr(element, 'get_children'):
                    for child in element.get_children():
                        if hasattr(child, 'polymer_type'): polymer_types.add(child.polymer_type)
                elif hasattr(element, 'polymer_type'):
                    polymer_types.add(element.polymer_type)
        else:
            polymer_types.add("protein")

        if "dna" in polymer_types or "rna" in polymer_types:
            if "dna" in polymer_types and "rna" not in polymer_types:
                seqtype_flag = "--seqtype=dna"
            elif "rna" in polymer_types and "dna" not in polymer_types:
                seqtype_flag = "--seqtype=rna"
            cline.append(seqtype_flag)

        if self.params_from_gui.get("iterations", 0) != 0:
            cline.extend(["--iter", str(self.params_from_gui["iterations"])])
        if self.params_from_gui.get("use_full_dm", False):
            cline.extend(["--full", "--full-iter"])

        return ' '.join(cline)
    
#########################END MODIFICATION ######################