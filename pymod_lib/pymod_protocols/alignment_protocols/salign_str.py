# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os

from Bio import SeqIO

from pymol import cmd

try:
    import modeller
except:
    pass

from ._salign_common import SALIGN_alignment, SALIGN_regular_alignment
from ._base_alignment._base_regular_alignment import Regular_structural_alignment
from ._base_alignment._gui import Regular_alignment_window_qt, Structural_alignment_base_window_qt


###################################################################################################
# SALIGN structural alignment.                                                                    #
###################################################################################################

class SALIGN_str_regular_alignment(SALIGN_regular_alignment, SALIGN_alignment, Regular_structural_alignment):

    alignment_program = "salign-str"
    protocol_name = "salign-str"

    def additional_initialization(self):
        self.tool = self.pymod.modeller

    def get_alignment_window_class_qt(self):
        return SALIGN_str_regular_window_qt

    def run_regular_alignment_program(self, sequences_to_align, output_file_name, use_parameters_from_gui=True, use_structural_information=False):
        if use_parameters_from_gui:
            pass
        self.run_salign_align3d(sequences_to_align, output_file_name)


    def run_salign_align3d(self, structures_to_align, output_file_name):
        """
        alignment.malign3d - align structures
        """

        # if len(structures_to_align)>2:
        #     self.build_salign_dendrogram_menu=True
        # else: # salign only output dendrogram_file when there are 3 sequences or more
        #     self.build_salign_dendrogram_menu=False

        shortcut_to_temp_files = os.path.join(self.pymod.current_project_dirpath,self.pymod.alignments_dirpath,output_file_name)
        struct_tup=list(range(0,len(structures_to_align)))
        for ii in range(0,len(structures_to_align)):
            struct_entry=structures_to_align[ii].get_pymol_selector()
            header = structures_to_align[ii].get_unique_index_header()
            chain_id=structures_to_align[ii].get_chain_id()
            struct_tup[ii]=(struct_entry,header,chain_id)

        # Change the working directory, so that the ouptut files will be created in the structures
        # directory.
        os.chdir(self.pymod.structures_dirpath)

        modeller.log.minimal()
        env = modeller.environ()
        aln = modeller.alignment(env)

        for (pdb_file_name, code, chain) in struct_tup:
            mdl = modeller.model(env, file=pdb_file_name,
                             model_segment=("FIRST:"+chain,"LAST:"+chain))
            aln.append_model(mdl, atom_files=pdb_file_name, align_codes=code)

        for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                            ((1., 0.5, 1., 1., 1., 0.), False, True),
                                            ((1., 1., 1., 1., 1., 0.), True, False)):
            aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                   rr_file="$(LIB)/as1.sim.mat", overhang=30,
                   gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                   gap_gap_score=0, gap_residue_score=0,
                   dendrogram_file= shortcut_to_temp_files + ".tree",
                   alignment_type="tree", feature_weights=weights,
                   improve_alignment=True, fit=True, write_fit=write_fit,
                   write_whole_pdb=whole,output="ALIGNMENT QUALITY")

        aln.write(file=shortcut_to_temp_files +".ali", alignment_format="PIR")

        aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
               gap_gap_score=0, gap_residue_score=0,
               dendrogram_file=shortcut_to_temp_files + '.tree',
               alignment_type='progressive', feature_weights=[0]*6,
               improve_alignment=False, fit=False, write_fit=True,
               write_whole_pdb=False,output='QUALITY')

        # Returns back to the project dir from the project/Structures directory.
        os.chdir(self.pymod.current_project_dirpath)

        # SALIGN does not superpose ligands. The generated "*_fit.pdb"
        # files are therefore ligandless. The following loop superposes
        # original structure to saligned structures, and replaces
        # "*_fit.pdb" files with the superposed liganded original structure.
        for pymod_element, (pdb_file_name_root, code, chain) in zip(structures_to_align, struct_tup):
            # Updates the name of the chains PDB files.
            fixed=os.path.join(self.pymod.structures_dirpath, pdb_file_name_root + "_fit.pdb")
            pymod_element.set_current_chain_file(os.path.join(self.pymod.current_project_dirpath, self.pymod.structures_dirpath, pdb_file_name_root + "_fit.pdb"))
            cmd.load(fixed,"salign_fixed_fit")
            if hasattr(cmd,"super"): # super is sequence-independent
                cmd.super(pdb_file_name_root,"salign_fixed_fit")
            else: # PyMOL 0.99 does not have cmd.super
                cmd.align(pdb_file_name_root,"salign_fixed_fit")
            cmd.set("retain_order", 1)
            cmd.save(fixed, pdb_file_name_root) # quick-and-dirty
            cmd.set("retain_order", 0)
            cmd.delete("salign_fixed_fit")

        # Convert the PIR format output file into a clustal format file.
        record = SeqIO.parse(shortcut_to_temp_files + '.ali', "pir")
        SeqIO.write(record, shortcut_to_temp_files + ".aln", "clustal")


    def update_aligned_sequences(self):
        self.update_aligned_sequences_inserting_modres()

    def update_additional_information(self):
        SALIGN_regular_alignment.update_additional_information(self)
        Regular_structural_alignment.update_additional_information(self)


class SALIGN_str_base_window_qt(Structural_alignment_base_window_qt):

    def build_algorithm_options_widgets(self):
        self.build_rmsd_option()
        self.middle_formlayout.set_input_widgets_width("auto")

class SALIGN_str_regular_window_qt(SALIGN_str_base_window_qt, Regular_alignment_window_qt):
    pass
