# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import shutil
import re

import pymod_lib.pymod_vars as pmdt
from ._base_alignment._base_regular_alignment import Regular_sequence_alignment
from ._base_alignment._base_profile_alignment import Profile_alignment


class Clustal_regular_alignment(Regular_sequence_alignment):

    def update_additional_information(self):
        """
        Sets the guide tree file path once the alignment has been performed.
        """

        if len(self.elements_to_align) > 2 and self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            # Builds a permanent copy of the original temporary .dnd file.
            temp_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, self.protocol_output_file_name+".dnd")
            new_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_%s_guide_tree.dnd" % (self.pymod.alignments_files_names, self.alignment_element.unique_index))
            shutil.copy(temp_dnd_file_path, new_dnd_file_path)

            # Edit the new .dnd file to insert the actual names of the sequences.
            dnd_file_handler = open(new_dnd_file_path, "r")
            dnd_file_lines = dnd_file_handler.readlines()
            dnd_file_handler.close()
            new_dnd_file_lines = []
            for line in dnd_file_lines:
                for m in re.findall(pmdt.unique_index_header_regex, line):
                    line = line.replace(m, self.elements_to_align_dict[m].my_header)
                new_dnd_file_lines.append(line)
            dnd_file_handler = open(new_dnd_file_path, "w")
            for line in new_dnd_file_lines:
                dnd_file_handler.write(line)
            dnd_file_handler.close()

            # # ClustalO produces a .dnd file without changing the ":" characters in the name of the
            # # PDB chains and this gives problems in displaying the names when using Phylo. So the
            # # ":" characters have to be changed in "_".
            # if self.protocol_name == "clustalo":
            #     old_dnd_file = open(new_dnd_file_path,"rU")
            #     new_dnd_file_content = ''
            #     for dnd_item in old_dnd_file.readlines():
            #         if re.search(r"_Chain\:?\:",dnd_item):
            #             Chain_pos=dnd_item.find("_Chain:")+6
            #             dnd_item=dnd_item[:Chain_pos]+'_'+dnd_item[Chain_pos+1:]
            #         new_dnd_file_content+=dnd_item
            #     old_dnd_file.close()
            #     new_dnd_file = open(new_dnd_file_path,"w")
            #     new_dnd_file.write(new_dnd_file_content)
            #     new_dnd_file.close()

            self.alignment_element.set_tree_file_path(new_dnd_file_path)


class Clustal_profile_alignment(Profile_alignment):

    def run_sequence_to_profile_alignment_program(self):
        """
        Align sequences to a target profile by clustalw/clustalo.
        """

        # List of sequences belonging to profile to be kept (target cluster).
        target_cluster_element = self.selected_clusters_list[self.target_cluster_index]
        target_profile_elements = target_cluster_element.get_children()

        # List of sequences to be appended to target cluster.
        self.elements_to_add = [e for e in self.pymod.get_selected_sequences() if not e in target_profile_elements]

        # create target cluster file
        profile_file_name = "cluster_0"
        profile_file_shortcut=os.path.join(self.pymod.alignments_dirpath, profile_file_name+".fasta")
        self.pymod.build_sequence_file(target_profile_elements, profile_file_name,
            file_format="fasta", remove_indels=False, unique_indices_headers=True)

        # create sequence file for sequences to be appended to target cluster
        sequences_to_add_file_name = "cluster_1"
        sequences_to_add_file_shortcut=os.path.join(self.pymod.alignments_dirpath, sequences_to_add_file_name+".fasta")
        self.pymod.build_sequence_file(self.elements_to_add, sequences_to_add_file_name,
            file_format="fasta", remove_indels=True, unique_indices_headers=True)

        # Output file name.
        sequence_to_profile_output = "al_result"
        output_file_shortcut = os.path.join(self.pymod.alignments_dirpath, sequence_to_profile_output)

        # Actually run the sequence to profile alignment.
        cline = self.prepare_sequence_to_profile_commandline(profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut)
        self.pymod.execute_subprocess(cline)

        # Converts the .aln output file into a .txt file, that will be used to update the sequences
        # loaded in PyMod.
        self.build_elements_to_align_dict(target_profile_elements+self.elements_to_add)
        self.protocol_output_file_name = sequence_to_profile_output


    def run_profile_to_profile_alignment_program(self):
        # Sequences in selected clusters will all be aligned. Sequences not in selected clusters,
        # will not be aligned.
        for cluster in self.selected_clusters_list:
            self.elements_to_align += list(cluster.get_children())

        self.profiles_to_join_file_list=[] # two MSA files

        for (i,cluster) in enumerate(self.selected_clusters_list):
            file_name = "cluster_" + str(i) # Build FASTA with the MSAs.
            children = cluster.get_children()

            # Builds a series of alignment files for each selected cluster.
            self.pymod.build_sequence_file(children, file_name, file_format="clustal", remove_indels = False, unique_indices_headers=True)
            self.profiles_to_join_file_list.append(file_name)

        profile_alignment_output = "al_result"

        output_file_shortcut=os.path.join(self.pymod.alignments_dirpath, profile_alignment_output)

        profile1=os.path.join(self.pymod.alignments_dirpath, self.profiles_to_join_file_list[0]+".aln")

        for profile2 in self.profiles_to_join_file_list[1:]:
            profile2=os.path.join(self.pymod.alignments_dirpath, profile2+".aln")
            cline = self.prepare_profile_to_profile_commandline(profile1, profile2, output_file_shortcut)
            self.pymod.execute_subprocess(cline)
            profile1=output_file_shortcut+'.aln'

        self.build_elements_to_align_dict(self.elements_to_align)
        self.protocol_output_file_name = profile_alignment_output
