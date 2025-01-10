# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os

from Bio import SeqIO

import pymod_lib.pymod_vars as pmdt
from pymod_lib.pymod_seq.seq_io import convert_sequence_file_format

from ._base_alignment._base_regular_alignment import Regular_sequence_alignment
from ._base_alignment._base_profile_alignment import Profile_alignment
from ._salign_common import SALIGN_alignment, SALIGN_regular_alignment

from ._base_alignment._gui import Regular_alignment_window_qt, Profile_alignment_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_radioselect_qt

try:
    import modeller
except:
    pass


class SALIGN_seq_alignment(SALIGN_alignment):
    """
    Mixin class for SALIGN sequence alignments (both regular and profile alignments).
    """

    alignment_program = "salign-seq"

    def additional_initialization(self):
        self.tool = self.pymod.modeller

    def get_options_from_gui(self):
        self.use_str_information = self.alignment_window.get_use_str_information_var()
        return True

    def check_structural_information(self):
        if self.use_str_information:
            for e in self.elements_to_align:
                if "X" in e.my_sequence:
                    print("- WARNING: could not use structural information in this alignment because of an unknown residue in '%s'." % (e.my_header))
                    self.use_str_information = False
                    break
                if e.polymer_type == "nucleic_acid":
                    print("- WARNING: could not use structural information in this alignment because '%s' is a nucleic acid element." % (e.my_header))
                    self.use_str_information = False
                    break

    def run_regular_alignment_program(self, sequences_to_align, output_file_name, use_parameters_from_gui=True, use_structural_information=False):
        self.check_structural_information()

        if use_parameters_from_gui:
            use_structural_information = self.use_str_information

        self.run_salign_malign(sequences_to_align, output_file_name, use_structural_information)

    def run_salign_malign(self, sequences_to_align, output_file_name, use_structural_information):
        """
        alignment.malign - align sequences
        alignment.align2d - sequence-structure alignment
        """

        shortcut_to_temp_files = os.path.join(self.pymod.alignments_dirpath, output_file_name)

        # The .pir file will be written in a different way if the user decides to use
        # structural information in the alignment.
        self.pymod.build_sequence_file(self.elements_to_align, output_file_name, file_format="pir",
            unique_indices_headers=True, use_structural_information=use_structural_information)

        modeller.log.minimal()
        env = modeller.environ()
        env.io.atom_files_directory = ['.', self.pymod.structures_dirpath]
        if self.use_hetatm:
            env.io.hetatm = True
        aln = modeller.alignment(env,
                                 file=shortcut_to_temp_files +".ali",
                                 alignment_format='PIR')
        if use_structural_information:
            env.libs.topology.read(file="$(LIB)/top_heav.lib")
            # # Structure sensitive variable gap penalty alignment:
            # aln.salign(auto_overhang=True,
            #     gap_penalties_1d=(-100, 0),
            #     gap_penalties_2d=(3.5,3.5,3.5,.2,4.,6.5,2.,0.,0.),
            #     gap_function=True, # structure-dependent gap penalty
            #     feature_weights=(1., 0., 0., 0., 0., 0.),
            #     similarity_flag=True,
            #     alignment_type='tree', #output='ALIGNMENT',
            #     dendrogram_file=shortcut_to_temp_files+".tree")
            aln.salign(rr_file='$(LIB)/as1.sim.mat',  # Substitution matrix used
                       # output='',
                       max_gap_length=20,
                       gap_function=True,              # If False then align2d not done
                       feature_weights=(1., 0., 0., 0., 0., 0.),
                       gap_penalties_1d=(-100, 0),
                       gap_penalties_2d=(3.5, 3.5, 3.5, 0.2, 4.0, 6.5, 2.0, 0.0, 0.0),
                       similarity_flag=True,   # Ensuring that the dynamic programming
                                               # matrix is not scaled to a difference matrix
                       alignment_type="tree",
                       dendrogram_file=shortcut_to_temp_files+".tree")


        else:
            aln.salign(auto_overhang=True, gap_penalties_1d=(-450, 0),
               alignment_type='tree', output='ALIGNMENT',
               dendrogram_file=shortcut_to_temp_files+".tree")
        aln.write(file=shortcut_to_temp_files +'.ali', alignment_format='PIR')

        # Convert output_file_name.ali to alignment_tmp.fasta.
        record = SeqIO.parse(shortcut_to_temp_files + ".ali", "pir")
        SeqIO.write(record, shortcut_to_temp_files + ".aln", "clustal")


    def salign_profile_profile_alignment(self, output_file_name="al_result", use_structural_information=False):

        profile1_name = self.profiles_to_join_file_list[0]+".ali"
        profile1_shortcut = os.path.join(self.pymod.alignments_dirpath, profile1_name)

        modeller.log.minimal()
        env = modeller.environ()
        env.io.atom_files_directory = ['.', self.pymod.structures_dirpath]
        if self.use_hetatm:
            env.io.hetatm = True
        env.libs.topology.read(file="$(LIB)/top_heav.lib")

        for profile2 in [os.path.join(self.pymod.alignments_dirpath,
            e+".ali") for e in self.profiles_to_join_file_list[1:]]:
            # cat profile2 to profile1 and return number of sequences
            # in the original profile1

            with open(profile1_shortcut,'r') as p_fh:
                ali_txt1 = p_fh.read()
            with open(profile2,'r') as p_fh:
                ali_txt2 = p_fh.read()
            align_block = len([e for e in ali_txt1.splitlines() if e.startswith('>')])


            with open(profile1_shortcut, 'w') as p_fh:
                p_fh.write(ali_txt1+ali_txt2)

            aln = modeller.alignment(env, file=profile1_shortcut, alignment_format="PIR")
            if use_structural_information:
                env.libs.topology.read(file='$(LIB)/top_heav.lib')
                aln.salign(rr_file='${LIB}/blosum62.sim.mat',
                    gap_penalties_1d=(-500, 0), output='',
                    align_block=align_block, #max_gap_length=20,
                    align_what='PROFILE', alignment_type="PAIRWISE",
                    comparison_type='PSSM',
                    gap_function=True,#structure-dependent gap penalty
                    feature_weights=(1., 0., 0., 0., 0., 0.),
                    gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),
                    similarity_flag=True,
                    substitution=True,smooth_prof_weight=10.0)
            else:
                aln.salign(rr_file='${LIB}/blosum62.sim.mat',
                gap_penalties_1d=(-500, 0), output='',
                align_block=align_block,   # no. of seqs. in first MSA
                align_what='PROFILE', alignment_type='PAIRWISE',
                comparison_type='PSSM',
                similarity_flag=True, substitution=True,
                smooth_prof_weight=10.0) # For mixing data with priors

            #write out aligned profiles (MSA)
            aln.write(file=profile1_shortcut, alignment_format="PIR")

        convert_sequence_file_format(profile1_shortcut, "pir", "clustal", output_filename=output_file_name)


    def update_aligned_sequences(self):
        self.update_aligned_sequences_inserting_modres()


###################################################################################################
# SALIGN sequence alignments.                                                                     #
###################################################################################################

class SALIGN_seq_regular_alignment(SALIGN_regular_alignment, SALIGN_seq_alignment, Regular_sequence_alignment):

    protocol_name = "salign-seq"

    def get_alignment_window_class_qt(self):
        return SALIGN_seq_regular_window_qt


###################################################################################################
# SALIGN sequence alignments.                                                                     #
###################################################################################################

class SALIGN_seq_profile_alignment(SALIGN_seq_alignment, Profile_alignment):

    protocol_name = "salign-seq"

    def get_alignment_window_class_qt(self):
        return SALIGN_seq_profile_window_qt


    def run_sequence_to_profile_alignment_program(self):

        # List of sequences of profile to be kept (target cluster)
        target_cluster_element = self.selected_clusters_list[self.target_cluster_index]
        alignment_to_keep_elements = target_cluster_element.get_children()

        # Used by generate_highest_identity_pairs_list
        self.selected_sequences_in_target_alignment = alignment_to_keep_elements

        # List of the selected sequences to be appended to target cluster.
        self.elements_to_add = [e for e in self.pymod.get_selected_sequences() if not e in alignment_to_keep_elements]

        #-----------------------------------------------------------------------------------------
        # Perform a first sequence alignment between all selected sequences and sequences in the -
        # target cluster.                                                                        -
        #-----------------------------------------------------------------------------------------
        initial_alignment_name = "all_temporary"
        self.elements_to_align = alignment_to_keep_elements + self.elements_to_add

        self.check_structural_information()

        # Perform sequence alignment even if sequence-structure alignment was requested, because the
        # former is signficantly faster.
        self.run_regular_alignment_program(self.elements_to_align, initial_alignment_name,
                                           use_parameters_from_gui=False, use_structural_information=False)

        #-----------------------------------------------------------------------------------------
        # For each sequence to be appended to the alignment, finds the most similiar sequence in -
        # the target cluster according to previous multiple sequence alignment.                  -
        #-----------------------------------------------------------------------------------------
        highest_identity_pairs_list = self.generate_highest_identity_pairs_list(initial_alignment_name)
        max_identity_list = [max(row) for row in highest_identity_pairs_list] # map(max, highest_identity_pairs_list)

        # Sort 'self.elements_to_add' according to 'max_identity_list'.
        max_identity_list, self.elements_to_add = zip(*sorted(zip(max_identity_list, self.elements_to_add),
                                                              reverse=True,
                                                              key=lambda t: t[0]))

        #-------------------------------------
        # Construct a PIR format input file. -
        #-------------------------------------
        self.profiles_to_join_file_list=[]
        profiles=[alignment_to_keep_elements]+[[e] for e in self.elements_to_add]

        for (i,children) in enumerate(profiles):
            file_name = "cluster_" + str(i)
            self.pymod.build_sequence_file(children, file_name,
                                           file_format="pir", remove_indels=False,
                                           use_structural_information=self.use_str_information,
                                           unique_indices_headers=True)
            self.profiles_to_join_file_list.append(file_name)

        #-----------------------------------------------------------------------------------
        # Sequentially apply profile-profile alignment to each element of elements_to_add. -
        #-----------------------------------------------------------------------------------
        profile_alignment_output = "al_result"
        self.salign_profile_profile_alignment(output_file_name=profile_alignment_output, use_structural_information=self.use_str_information)
        self.build_elements_to_align_dict(self.elements_to_align)
        self.protocol_output_file_name = profile_alignment_output


    def run_profile_to_profile_alignment_program(self):

        # Sequences in selected clusters will all be aligned. Sequences not in selected clusters,
        # will not be aligned.
        for cluster in self.selected_clusters_list:
            self.elements_to_align += list(cluster.get_children())

        self.check_structural_information()

        self.profiles_to_join_file_list=[] # two MSA files

        for (i,cluster) in enumerate(self.selected_clusters_list):
            file_name = "cluster_" + str(i) # Build FASTA with the MSAs.
            children = cluster.get_children()
            # Builds a series of alignment files for each selected cluster.
            # self.pymod.build_sequence_file(children, file_name, file_format="clustal", remove_indels = False, unique_indices_headers=True)
            self.pymod.build_sequence_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information=self.use_str_information, unique_indices_headers=True)
            self.profiles_to_join_file_list.append(file_name)

        profile_alignment_output = "al_result"
        output_file_shortcut=os.path.join(self.pymod.alignments_dirpath, profile_alignment_output)
        profile1=os.path.join(self.pymod.alignments_dirpath, self.profiles_to_join_file_list[0]+".aln")

        self.salign_profile_profile_alignment(profile_alignment_output, use_structural_information=self.use_str_information)

        self.build_elements_to_align_dict(self.elements_to_align)
        self.protocol_output_file_name = profile_alignment_output


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class SALIGN_seq_base_window_qt:

    def build_algorithm_options_widgets(self):
        if self.protocol.structures_are_selected:
            # Use structure information to guide sequence alignment.
            self.salign_seq_struct_rds = PyMod_radioselect_qt(label_text="Use structural information",
                                                        buttons=('Yes', 'No'))
            self.salign_seq_struct_rds.setvalue("No")
            self.middle_formlayout.add_widget_to_align(self.salign_seq_struct_rds)
            self.middle_formlayout.set_input_widgets_width("auto")

    def get_use_str_information_var(self):
        if self.protocol.structures_are_selected:
            return pmdt.yesno_dict[self.salign_seq_struct_rds.getvalue()]
        else:
            return False


class SALIGN_seq_regular_window_qt(SALIGN_seq_base_window_qt, Regular_alignment_window_qt):
    pass

class SALIGN_seq_profile_window_qt(SALIGN_seq_base_window_qt, Profile_alignment_window_qt):
    pass
