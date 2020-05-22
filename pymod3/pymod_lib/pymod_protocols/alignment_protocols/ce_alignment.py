# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for performing CE-alignments in PyMod.
"""

import os
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymol import cmd

from pymod_lib.pymod_seq.seq_star_alignment import save_cstar_alignment

from ._base_alignment._base_regular_alignment import Regular_structural_alignment
from ._base_alignment._gui import Regular_alignment_window_qt, Structural_alignment_base_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_combobox_qt


##################################################################################################
# CE alignment.                                                                                  #
##################################################################################################

class CEalign_alignment:

    alignment_program = "ce"

    def additional_initialization(self):
        self.tool = None

    def alignment_program_exists(self):
        """
        In PyMOL 2 CE-align is alway present.
        """
        return True

    def alignment_program_not_found(self):
        title = "CE-alignment Error"
        message = "CE-alignment is not available on your PyMod installation. If you want to use this function please see CE-alignment installation instructions on PyMod's User Guide."
        self.pymod.main_window.show_error_message(title, message)

    def update_aligned_sequences(self):
        self.update_aligned_sequences_inserting_modres()


class CEalign_regular_alignment(CEalign_alignment, Regular_structural_alignment):

    protocol_name = "ce"

    def get_alignment_window_class_qt(self):
        return CEalign_regular_window_qt


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):

        # Change the order of the 'sequences_to_align' so that the center star is the first element
        # in the list.
        _sequences_to_align = sequences_to_align[:]
        reference_id = self.alignment_window.get_reference_id()
        reference_element = _sequences_to_align.pop(reference_id)
        _sequences_to_align.insert(0, reference_element)

        # Run the alignment.
        self.run_ce_alignment(_sequences_to_align, output_file_name=output_file_name)


    def run_ce_alignment(self, structures_to_align, output_file_name, use_seq_info=False):
        """
        Used to launch Ce_align.
        """

        backup_list = structures_to_align[:]

        #-----------------------------------------------------------------------------
        # Align the first two structures and produces an ce_temp.txt alignment file. -
        #-----------------------------------------------------------------------------

        temp_ceali_prefix = "ce_temp"
        temp_ceali = temp_ceali_prefix + "_0"
        temp_ceali_list = [temp_ceali]

        current_elements_to_align = backup_list[0:2]

        self.ce_align(current_elements_to_align, output_file_name=temp_ceali)


        if len(backup_list) > 2:

            max_row = 0

            #-------------------------------------------------------------------
            # Align the rest of the structures to the first one progressively. -
            #-------------------------------------------------------------------

            mceali_count = 1
            for n in range(2, len(backup_list)):
                current_elements_to_align = [backup_list[0], backup_list[n]]
                temp_ceali_n = "%s_%s" % (temp_ceali_prefix, mceali_count)
                self.ce_align(current_elements_to_align, output_file_name=temp_ceali_n)

                temp_ceali_list.append(temp_ceali_n)
                mceali_count += 1

            #-------------------------------------------------------------------------
            # Join the alignments using a modification of the center star alignment. -
            #-------------------------------------------------------------------------

            # center_star_id = backup_list[0].get_unique_index_header()
            aligned_pairs = []
            seqs = []
            all_ids = []

            # Builds a list of pairwise alignments.
            for temp_ceali_i_id, temp_ceali_i in enumerate(temp_ceali_list):
                c0_recs = list(SeqIO.parse(os.path.join(self.pymod.alignments_dirpath, temp_ceali_i + ".aln"), "clustal"))
                aligned_pairs.append([str(c0_recs[0].seq), str(c0_recs[1].seq)])
                if temp_ceali_i_id == 0:
                    seqs.append(str(c0_recs[0].seq).replace("-", ""))
                    all_ids.append(c0_recs[0].id)
                seqs.append(str(c0_recs[1].seq).replace("-", ""))
                all_ids.append(c0_recs[1].id)

            for ei_id, ei in enumerate(backup_list[1:]):
                for ej in backup_list[ei_id+1:]:
                    if ei is ej:
                        continue
                    aligned_pairs.append(None)


            # Joins the alignments.
            save_cstar_alignment(seqs=seqs, all_ids=all_ids,
                                 pairwise_alis=aligned_pairs,
                                 output_filepath=os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln"))

            #-----------------------------------------------------------------------------------
            # Complete by cleaning up the temporary files and by creating a final output file. -
            #-----------------------------------------------------------------------------------

            for temp_ceali_n in temp_ceali_list:
                os.remove(os.path.join(self.pymod.alignments_dirpath, temp_ceali_n + ".aln"))


        else:
            shutil.move(os.path.join(self.pymod.alignments_dirpath, temp_ceali + ".aln"),
                        os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln"))


    def ce_align(self, elements_to_align, output_file_name=None, use_seq_info=False):
        """
        Actually performs the structural alignment.
        """

        #----------------------------------------------------
        # Run CE-alignment using the PyMOL built-in module. -
        #----------------------------------------------------

        retain_order = True
        if retain_order:
            cmd.set("retain_order", 1)

        sel1 = elements_to_align[0].get_pymol_selector()
        sel2 = elements_to_align[1].get_pymol_selector()

        # Sets temporary names.
        tsel1 = "t" + elements_to_align[0].get_unique_index_header()
        tsel2 = "t" + elements_to_align[1].get_unique_index_header()

        cmd.set_name(sel1, tsel1)
        cmd.set_name(sel2, tsel2)

        # Actually performs the alignment.
        a = cmd.cealign(target=tsel1, mobile=tsel2, object="pymod_temp_cealign")
        # cmd.center('%s and %s' % (tsel1, tsel2))
        # cmd.zoom('%s and %s' % (tsel1, tsel2))

        # Updates the names of the chains PDB files and saves these new files.
        saved_file1 = sel1 + "_aligned.pdb"
        saved_file2 = sel2 + "_aligned.pdb"
        # elements_to_align[0].structure.chain_pdb_file_name = saved_file1
        # elements_to_align[1].structure.chain_pdb_file_name = saved_file2
        cmd.save(os.path.join(self.pymod.structures_dirpath, saved_file1), tsel1)
        cmd.save(os.path.join(self.pymod.structures_dirpath, saved_file2), tsel2)

        # Finally saves the structural alignment between the sequences.
        cmd.save(os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln"), "pymod_temp_cealign")
        cmd.delete("pymod_temp_cealign")

        # Converts it in .aln format.
        recs = SeqIO.parse(os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln"), "clustal")
        new_recs = []
        for rec, pymod_element in zip(recs, (elements_to_align[0], elements_to_align[1])):
            new_rec_id = "_".join(rec.id[1:].split("_")[0:3])
            new_rec_seq = str(rec.seq).replace("?", "X") # Replaces modified residues.

            new_recs.append(SeqRecord(Seq(new_rec_seq), id=new_rec_id))

        SeqIO.write(new_recs, os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln"), "clustal")

        # Sets the names of the objects back to original ones.
        cmd.set_name(tsel1, sel1)
        cmd.set_name(tsel2, sel2)

        if retain_order:
            cmd.set("retain_order", 0)


###################################################################################################
# Classes for the GUI.                                                                            #
###################################################################################################

class CEalign_base_window_qt(Structural_alignment_base_window_qt):

    def build_algorithm_options_widgets(self):

        # Reference structure combobox.
        if len(self.protocol.selected_elements) > 2:
            structures_list = [element.my_header for element in self.protocol.selected_elements]
            self.reference_combobox = PyMod_combobox_qt(label_text="Reference Structure",
                                                        items=structures_list)
            self.reference_combobox.combobox.setCurrentIndex(0)
            self.middle_formlayout.add_widget_to_align(self.reference_combobox)

        # RMSD options.
        self.build_rmsd_option()

        self.middle_formlayout.set_input_widgets_width(130)


    def get_reference_id(self):
        try:
            return self.reference_combobox.get_index()
        except:
            print("- Warning: could not obtain the reference structure id.")
            return 0


class CEalign_regular_window_qt(CEalign_base_window_qt, Regular_alignment_window_qt):
    pass
