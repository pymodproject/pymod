# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
ClustalW.
"""

import os

import subprocess

# Protocols.
from ._clustal_common import Clustal_regular_alignment, Clustal_profile_alignment
from pymod_lib.pymod_main._elements_interactions import PyMod_elements_interactions
from pymod_lib.pymod_element import PyMod_sequence_element, PyMod_cluster_element




# GUI.
from ._base_alignment._gui import Regular_alignment_window_qt, Profile_alignment_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_entryfield_qt, PyMod_radioselect_qt


class Clustalw_alignment:
    # General Clustal W alignments.
    protocol_name = "clustalw"
    
    def additional_initialization(self):
        self.tool = self.pymod.clustalw

    
class Clustalw_regular_alignment(Clustalw_alignment, Clustal_regular_alignment):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.polymer_type = "Protein"  # Default value

    #############################MODIFICATION on 01/04/2025###################################
    ###I moved this function from the class Clustalw_alignment to this class
    def get_options_from_gui(self):
        self.params_from_gui = {}
        error_message = None  

        #get the polymer_type
        polymer_type = self.get_polymer_type()

        if polymer_type is None:
            self.pymod.main_window.show_error_message("Parameters Error", "Polymer type could not be determined.")
            return None

        polymer_type = str(polymer_type)

        try:
            if polymer_type == "Protein":
                self.params_from_gui["selected_matrix"] = self.alignment_window.get_matrix_value()
            
            #Common parameters for all polymer types 
            self.params_from_gui["gapopen_value"] = float(self.alignment_window.get_gapopen_value())
            self.params_from_gui["gapextension_value"] = float(self.alignment_window.get_gapextension_value())
        
        except Exception as e:
            error_message= f"Invalid parameter: {str(e)}"
        
        if error_message:
            self.pymod.main_window.show_error_message("Parameters Error", error_message)
            return None
        else:
            return True
        
    def get_alignment_window_class_qt(self):
        
       return Clustalw_regular_window_qt 
        

    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        
        #############################MODIFICATION on 01/04/2025###################################
        #Run ClustalW with the defined polymer type automatically 
        polymer_type = self.get_polymer_type()
        
        if polymer_type=="Protein":
            self.run_clustalw(sequences_to_align,
                            output_file_name=output_file_name,
                            matrix=self.params_from_gui["selected_matrix"],
                            gapopen=self.params_from_gui["gapopen_value"],
                            gapext=self.params_from_gui["gapextension_value"],
                            polymer_type=polymer_type)
        else:
            self.run_clustalw(sequences_to_align,
                            output_file_name=output_file_name,
                            #matrix=self.params_from_gui["selected_matrix"], ###Not necessary for nucleic acids
                            gapopen=self.params_from_gui["gapopen_value"],
                            gapext=self.params_from_gui["gapextension_value"],
                            polymer_type=polymer_type)
            


    def run_clustalw(self, sequences_to_align, output_file_name, matrix="blosum", gapopen=10, gapext=0.2,polymer_type="None"):
        # This method allows to interact with the local ClustalW.

        # First build an input FASTA file containing the sequences to be aligned.
        self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)
        # Sets the full paths of input and output files.
        input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
        output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")

        #############################MODIFICATION on 01/04/2025###################################
        # Check if polymer_type is DNA or RNA and use the correct command.
        polymer_type = self.get_polymer_type()
    
        if polymer_type == "Protein":

            clustalw_path = self.pymod.clustalw["exe_file_path"].get_value()
            cline = [clustalw_path, f"-infile={input_file_path}", f"-outfile={output_file_path}", "-outorder=INPUT",
                     f"-matrix={matrix}", f"-gapopen={str(gapopen)}", f"-gapext={str(gapext)}"]
            
        else:
            # The penalties are the ones the window holds. They used to be overwritten
            # here with 15.0 and 6.66, the nucleotide defaults of ClustalW, which the
            # window already offers as the values it starts on: whatever was typed over
            # them was collected, carried down to here and dropped, and since the
            # figures were the same as the defaults nothing looked wrong unless the
            # user changed them.
            #
            # 'PROTEIN' and 'DNA' are the only two values -TYPE accepts ("-TYPE= :
            # PROTEIN or DNA" in its help), so RNA is sent as DNA. The IUB matrix which
            # goes with it scores the uracil as any other base.
            clustalw_path = self.pymod.clustalw["exe_file_path"].get_value()
            cline = [clustalw_path, f"-infile={input_file_path}", f"-outfile={output_file_path}",
                     "-outorder=INPUT", "-TYPE=DNA", "-DNAMATRIX=IUB", "-PWDNAMATRIX=IUB",
                     f"-PWGAPOPEN={gapopen}", f"-PWGAPEXT={gapext}",
                     f"-GAPOPEN={gapopen}", f"-GAPEXT={gapext}"]
        
        #############################MODIFICATION on 01/04/2025###################################
        # Show the command in the terminal for debugging
        print(f"Clustal W command: {' '.join(cline)}")
        self.pymod.execute_subprocess(' '.join(cline))
    

#######################################TO MODIFY in order to make profile alignmetn possible for nucleic acids###################
class Clustalw_profile_alignment(Clustalw_alignment, Clustal_profile_alignment):
    
    #############################MODIFICATION on 18/04/2025###################################
    ###added this function for Profile alignment too#####
    def get_options_from_gui(self):
        self.params_from_gui = {}
        error_message = None

        if not any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.selected_elements):
            polymer_type = "Protein"
        elif any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.selected_elements):
            if any(hasattr(e, "polymer_type") and e.polymer_type == "dna" for e in self.selected_elements):
                polymer_type = "DNA"
            elif any(hasattr(e, "polymer_type") and e.polymer_type == "rna" for e in self.selected_elements):
                polymer_type = "RNA"

        if polymer_type is None:
            self.pymod.main_window.show_error_message("Parameters Error", "Polymer type could not be determined.")
            return None

        polymer_type = str(polymer_type)

        try:
            if polymer_type == "Protein":
                self.params_from_gui["selected_matrix"] = self.alignment_window.get_matrix_value()

            self.params_from_gui["gapopen_value"] = float(self.alignment_window.get_gapopen_value())
            self.params_from_gui["gapextension_value"] = float(self.alignment_window.get_gapextension_value())

        except Exception as e:
            error_message= f"Invalid parameter: {str(e)}"
        
        if error_message:
            self.pymod.main_window.show_error_message("Parameters Error", error_message)
            return None
        else:
            return True

    def get_alignment_window_class_qt(self):
        return Clustalw_profile_window_qt

    #############################MODIFICATION on 18/04/2025###################################
    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        
        clustalw_path = self.tool["exe_file_path"].get_value()
        if not any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.selected_elements):
            polymer_type = "Protein"
        elif any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.selected_elements):
            if any(hasattr(e, "polymer_type") and e.polymer_type == "dna" for e in self.selected_elements):
                polymer_type = "DNA"
            elif any(hasattr(e, "polymer_type") and e.polymer_type == "rna" for e in self.selected_elements):
                polymer_type = "RNA"
        if polymer_type == "Protein":
            cline = [clustalw_path,
                     "-PROFILE1=" + profile_file_shortcut,
                     "-PROFILE2=" + sequences_to_add_file_shortcut,
                     "-SEQUENCES", "-OUTORDER=INPUT",
                     "-MATRIX=" + self.params_from_gui["selected_matrix"],
                     "-GAPOPEN=" + str(self.params_from_gui["gapopen_value"]), 
                     "-GAPEXT=" + str(self.params_from_gui["gapextension_value"]), 
                     "-OUTFILE=" + output_file_shortcut + ".aln"]
        else:
            # No matrix for a nucleic acid: the IUB one goes with -TYPE=DNA, which is
            # the only nucleotide type ClustalW accepts, RNA included. The penalties are
            # the ones the window holds, as on the protein branch above -- they used to
            # be overwritten here with the defaults of ClustalW.
            cline = [clustalw_path,
                     "-PROFILE1=" + profile_file_shortcut,
                     "-PROFILE2=" + sequences_to_add_file_shortcut,
                     "-SEQUENCES", "-OUTORDER=INPUT",
                     "-TYPE=DNA",
                     "-DNAMATRIX=IUB", "-PWDNAMATRIX=IUB",
                     "-GAPOPEN=" + str(self.params_from_gui["gapopen_value"]),
                     "-GAPEXT=" + str(self.params_from_gui["gapextension_value"]),
                     "-OUTFILE=" + output_file_shortcut + ".aln"]

        return ' '.join(cline)
    
    #############################MODIFICATION on 18/04/2025###################################
    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalw_path = self.tool["exe_file_path"].get_value()

        if not any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.selected_elements):
            polymer_type = "Protein"
        elif any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.selected_elements):
            if any(hasattr(e, "polymer_type") and e.polymer_type == "dna" for e in self.selected_elements):
                polymer_type = "DNA"
            elif any(hasattr(e, "polymer_type") and e.polymer_type == "rna" for e in self.selected_elements):
                polymer_type = "RNA"

        if polymer_type == "Protein":
            cline = [clustalw_path,
                     "-PROFILE1=" + profile1,
                     "-PROFILE2=" + profile2,
                     "-OUTORDER=INPUT",
                     "-MATRIX=" + self.params_from_gui["selected_matrix"],
                     "-GAPOPEN=" + str(self.params_from_gui["gapopen_value"]), 
                     "-GAPEXT=" + str(self.params_from_gui["gapextension_value"]), 
                     "-OUTFILE=" + output_file_shortcut + ".aln"]
        else:
            # As above: the penalties of the window, and DNA for every nucleic acid.
            cline = [clustalw_path,
                     "-PROFILE1=" + profile1,
                     "-PROFILE2=" + profile2,
                     "-OUTORDER=INPUT",
                     "-TYPE=DNA",
                     "-DNAMATRIX=IUB", "-PWDNAMATRIX=IUB",
                     "-GAPOPEN=" + str(self.params_from_gui["gapopen_value"]),
                     "-GAPEXT=" + str(self.params_from_gui["gapextension_value"]),
                     "-OUTFILE=" + output_file_shortcut + ".aln"]

        return ' '.join(cline)


###################################################################################################
# GUI.                                                                                            #
###################################################################################################


class Clustalw_base_window_qt:

    def build_algorithm_options_widgets(self):

        ########################MODIFICATION on 04/04/2025#######################################
        ##Access selected_elements through self.protocol 
        if not any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.protocol.selected_elements):
            polymer_type = "Protein"
        elif any(hasattr(e, "polymer_type") and e.polymer_type in ["dna", "rna"] for e in self.protocol.selected_elements):
            if any(hasattr(e, "polymer_type") and e.polymer_type == "dna" for e in self.protocol.selected_elements):
                polymer_type = "DNA"
            elif any(hasattr(e, "polymer_type") and e.polymer_type == "rna" for e in self.protocol.selected_elements):
                polymer_type = "RNA"

        try:
            
            #####################TOGLI PRIMA DEL COMMIT############################À
            #print("##################################")
            #print(self)
            #print(dir(self))
            #print(dir(self.protocol))
            #print(dir(self.protocol.selected_elements[0].polymer_type))

            if polymer_type=="Protein":
                
                self.clustal_matrices = ["Blosum", "Pam", "Gonnet", "Id"]
                self.clustal_matrices_dict = {"Blosum": "blosum", "Pam": "pam", "Gonnet": "gonnet", "Id": "id"}
                self.matrix_rds = PyMod_radioselect_qt(label_text="Scoring Matrix Selection",
                                                    buttons=self.clustal_matrices)
                self.matrix_rds.setvalue("Blosum")
                self.middle_formlayout.add_widget_to_align(self.matrix_rds)

                # Gap open entryfield.
                self.gapopen_enf = PyMod_entryfield_qt(label_text="Gap Opening Penalty",
                                                    value="10.0",
                                                    validate={'validator': 'real',
                                                                'min': 0, 'max': 1000})
                self.middle_formlayout.add_widget_to_align(self.gapopen_enf)

                # Gap extension entryfield.
                self.gapextension_enf = PyMod_entryfield_qt(label_text="Gap Extension Penalty",
                                                            value="0.2",
                                                            validate={'validator': 'real',
                                                                    'min': 0, 'max': 1000})
                self.middle_formlayout.add_widget_to_align(self.gapextension_enf)
    
            else:
                self.gapopen_enf = PyMod_entryfield_qt(label_text="Gap Opening Penalty",
                                            value="15.0",
                                            validate={'validator': 'real',
                                                        'min': 0, 'max': 1000})
                self.middle_formlayout.add_widget_to_align(self.gapopen_enf)

                # Gap extension entryfield.
                # 6.66, the nucleotide default of ClustalW itself, so that leaving the
                # box alone gives what the program would have done on its own.
                self.gapextension_enf = PyMod_entryfield_qt(label_text="Gap Extension Penalty",
                                                            value="6.66",
                                                            validate={'validator': 'real',
                                                                    'min': 0, 'max': 1000})
                self.middle_formlayout.add_widget_to_align(self.gapextension_enf)

        except Exception as e:
            print(f"Errore durante la costruzione dei widget: {e}")

        ################################END MODIFICATION####################################

    
    def get_matrix_value(self):
        return self.clustal_matrices_dict[self.matrix_rds.getvalue()]

    def get_gapopen_value(self):
        return self.gapopen_enf.getvalue(validate=True)

    def get_gapextension_value(self):
        return self.gapextension_enf.getvalue(validate=True)


class Clustalw_regular_window_qt(Clustalw_base_window_qt, Regular_alignment_window_qt):
    pass

class Clustalw_profile_window_qt(Clustalw_base_window_qt, Profile_alignment_window_qt):
    pass





"""ORIGINAL CODE 



class Clustalw_alignment:
    
    #General Clustal W alignments.
    

    protocol_name = "clustalw"

    def additional_initialization(self):
        self.tool = self.pymod.clustalw

    def get_options_from_gui(self):

        self.params_from_gui = {}
        error_from_gui = False

        try:
            self.params_from_gui["selected_matrix"] = self.alignment_window.get_matrix_value()
        except Exception as e:
            error_from_gui = True
            error_message = "Invalid Matrix."
        try:
            self.params_from_gui["gapopen_value"] = float(self.alignment_window.get_gapopen_value())
        except Exception as e:
            error_from_gui = True
            error_message = str(e) # "Invalid Gap Open Value."
        try:
            self.params_from_gui["gapextension_value"] = float(self.alignment_window.get_gapextension_value())
        except Exception as e:
            error_from_gui = True
            error_message = str(e) # "Invalid Gap Extension Value."

        if error_from_gui:
            self.pymod.main_window.show_error_message("Parameters Error", error_message)
            return False
        else:
            return True


class Clustalw_regular_alignment(Clustalw_alignment, Clustal_regular_alignment):

    def get_alignment_window_class_qt(self):
        return Clustalw_regular_window_qt


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_clustalw(sequences_to_align,
                          output_file_name=output_file_name,
                          matrix=self.params_from_gui["selected_matrix"],
                          gapopen=self.params_from_gui["gapopen_value"],
                          gapext=self.params_from_gui["gapextension_value"])


    def run_clustalw(self, sequences_to_align, output_file_name, matrix="blosum", gapopen=10, gapext=0.2):
        
        #This method allows to interact with the local ClustalW.
        

        # First build an input FASTA file containing the sequences to be aligned.
        self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)
        # Sets the full paths of input and output files.
        input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
        output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
        # Run an alignment with all the sequences using ClustalW command line, through Biopython.
        cline = ClustalwCommandline(self.pymod.clustalw["exe_file_path"].get_value(),
                infile=input_file_path, outfile=output_file_path, outorder="INPUT",
                matrix=matrix, gapopen=gapopen, gapext=gapext)
        # Show the command in the terminal for debugging
        print(f"Clustal W command: {cline}")
        self.pymod.execute_subprocess(str(cline))


class Clustalw_profile_alignment(Clustalw_alignment, Clustal_profile_alignment):

    def get_alignment_window_class_qt(self):
        return Clustalw_profile_window_qt


    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        clustalw_path = self.tool["exe_file_path"].get_value()
        cline='"'         +clustalw_path+'"'+ \
            ' -PROFILE1="'+profile_file_shortcut+'"'+ \
            ' -PROFILE2="'+sequences_to_add_file_shortcut+'" -SEQUENCES -OUTORDER=INPUT'+ \
            ' -MATRIX='   +self.params_from_gui["selected_matrix"] + \
            ' -GAPOPEN='  +str(self.params_from_gui["gapopen_value"]) + \
            ' -GAPEXT='   +str(self.params_from_gui["gapextension_value"]) + \
            ' -OUTFILE="' +output_file_shortcut+'.aln"'
        return cline


    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalw_path = self.tool["exe_file_path"].get_value()
        cline='"'          +clustalw_path+'"' \
            ' -PROFILE1="' +profile1+'"'+ \
            ' -PROFILE2="' +profile2+'" -OUTORDER=INPUT' \
            ' -MATRIX='    +self.params_from_gui["selected_matrix"]+ \
            ' -GAPOPEN='   +str(self.params_from_gui["gapopen_value"])+ \
            ' -GAPEXT='    +str(self.params_from_gui["gapextension_value"])+ \
            ' -OUTFILE="'  +output_file_shortcut+'.aln"'
        return cline


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Clustalw_base_window_qt:
    
    #Base class for ClustalW protocols.


    def build_algorithm_options_widgets(self):

        # Scoring matrix radioselect.
        self.clustal_matrices = ["Blosum", "Pam", "Gonnet", "Id"]
        self.clustal_matrices_dict = {"Blosum": "blosum", "Pam": "pam", "Gonnet": "gonnet", "Id": "id"}
        self.matrix_rds = PyMod_radioselect_qt(label_text="Scoring Matrix Selection",
                                               buttons=self.clustal_matrices)
        self.matrix_rds.setvalue("Blosum")
        self.middle_formlayout.add_widget_to_align(self.matrix_rds)

        # Gap open entryfield.
        self.gapopen_enf = PyMod_entryfield_qt(label_text="Gap Opening Penalty",
                                               value="10.0",
                                               validate={'validator': 'real',
                                                         'min': 0, 'max': 1000})
        self.middle_formlayout.add_widget_to_align(self.gapopen_enf)

        # Gap extension entryfield.
        self.gapextension_enf = PyMod_entryfield_qt(label_text="Gap Extension Penalty",
                                                    value="0.2",
                                                    validate={'validator': 'real',
                                                              'min': 0, 'max': 1000})
        self.middle_formlayout.add_widget_to_align(self.gapextension_enf)

        self.middle_formlayout.set_input_widgets_width("auto")


    def get_matrix_value(self):
        return self.clustal_matrices_dict[self.matrix_rds.getvalue()]

    def get_gapopen_value(self):
        return self.gapopen_enf.getvalue(validate=True)

    def get_gapextension_value(self):
        return self.gapextension_enf.getvalue(validate=True)


class Clustalw_regular_window_qt(Clustalw_base_window_qt, Regular_alignment_window_qt):
    pass

class Clustalw_profile_window_qt(Clustalw_base_window_qt, Profile_alignment_window_qt):
    pass


"""