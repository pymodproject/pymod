# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.











##################Original CODE#############################

"""

# TODO:
#   - add the possibility to save trees in the phylip format.

import os

from pymod_lib import pymod_vars

from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_protocol_window_qt, PyMod_radioselect_qt

from ._evolutionary_analysis_base import Evolutionary_analysis_protocol


tree_building_alg_dict = {"Neighbor Joining": "nj", "UPGMA": "upgma"}

class Tree_building(Evolutionary_analysis_protocol):

    def launch_from_gui(self):
        
        #It will check if a software to build a tree is available on the user's machine.
        
        self.tree_building_software = None
        can_build_tree = False
        if self.pymod.clustalw.tool_file_exists():
            self.tree_building_software = "clustalw"
            can_build_tree = True
        elif self.pymod.muscle.tool_file_exists():
            self.tree_building_software = "muscle"
            can_build_tree = True
        if can_build_tree:
            self.build_tree_building_window()
        else:
            title = "Tree building Error"
            message = "In order to build a tree out of an alignment you need to install either ClustalW or MUSCLE and specify an existent executable path of one of these tools in the PyMod Options Window."
            self.pymod.main_window.show_error_message(title, message)


    def check_tree_constructor_module(self):
        try:
            import Bio.Phylo.TreeConstruction
            return True
        except ImportError:
            return False


    def build_tree_building_window(self):
        
        #Builds a window with options to build a tree out of an alignment.
        

        self.tree_building_window = Tree_building_window_qt(self.pymod.main_window,
            protocol=self,
            title="Options for Tree Building",
            upper_frame_title="Here you can modify options for Tree Building",
            submit_command=self.run_tree_building_software)
        self.tree_building_window.show()


    def run_tree_building_software(self):
        # Saves a temporary input alignment file.
        alignment_file_name = "alignment_tmp"
        alignment_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.fasta')
        self.pymod.save_alignment_fasta_file(alignment_file_name, self.input_cluster_element.get_children())

        # Get the parameters from the GUI.
        clustering_algorithm = self.get_clustering_algorithm()

        # Prepares to run the tree-building algorithm.
        commandline = ""
        output_file_path = None

        if self.tree_building_software == "clustalw":
            commandline =  '"%s"' % (self.pymod.clustalw["exe_file_path"].get_value())
            commandline += ' -TREE -INFILE="%s"' % (alignment_file_path)
            commandline += ' -OUTPUTTREE=phylip'
            if self.get_distance_correction_val():
                commandline += ' -KIMURA'
            if self.get_exclude_gaps_val():
                commandline += ' -TOSSGAPS'
            # if self.get_boostrap_val():
            #     commandline += ' -SEED='+str(random.randint(0,1000))
            #     commandline += ' -BOOTLABELS=node'
            if clustering_algorithm == "nj":
                commandline += ' -CLUSTERING=NJ'
            elif clustering_algorithm == "upgma":
                commandline += ' -CLUSTERING=UPGMA'
            output_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.ph')

        elif self.tree_building_software == "muscle":
            commandline =  '"%s"' % (self.pymod.muscle["exe_file_path"].get_value())
            commandline += ' -maketree -in %s' % (alignment_file_path)
            output_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.phy')
            commandline += ' -out %s' % (output_file_path)
            if clustering_algorithm == "nj":
                commandline += ' -cluster neighborjoining'
            elif clustering_algorithm == "upgma":
                pass

        # Actually runs the tree building algorithm.
        self.pymod.execute_subprocess(commandline)

        # Remove temporary files.
        new_tree_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_%s_align_tree.phy" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index))
        os.rename(output_file_path, new_tree_file_path)
        os.remove(alignment_file_path)
        self.tree_building_window.destroy()

        # Reads the output tree file with Phylo and displays its content using PyMod plotting
        # engine.
        self.pymod.show_tree(new_tree_file_path)


    def get_clustering_algorithm(self):
        return tree_building_alg_dict[self.tree_building_window.algorithm_rds.getvalue()]

    def get_boostrap_val(self):
        return pymod_vars.yesno_dict[self.tree_building_window.bootstrap_rds.getvalue()]

    def get_distance_correction_val(self):
        return pymod_vars.yesno_dict[self.tree_building_window.distance_correction_rds.getvalue()]

    def get_exclude_gaps_val(self):
        return pymod_vars.yesno_dict[self.tree_building_window.exclude_gaps_rds.getvalue()]


class Tree_building_window_qt(PyMod_protocol_window_qt):

    def build_protocol_middle_frame(self):

        # Add some options.
        self.algorithm_rds = PyMod_radioselect_qt(label_text="Clustering Algorithm",
                                                  buttons=list(sorted(tree_building_alg_dict.keys())))
        self.algorithm_rds.setvalue("Neighbor Joining")
        self.middle_formlayout.add_widget_to_align(self.algorithm_rds)

        if self.protocol.tree_building_software == "clustalw":
            # Kimura distance correction.
            self.distance_correction_rds = PyMod_radioselect_qt(label_text="Use Distance Correction",
                                                                buttons=('Yes', 'No'))
            self.distance_correction_rds.setvalue("No")
            self.middle_formlayout.add_widget_to_align(self.distance_correction_rds)

            # Toss gaps.
            self.exclude_gaps_rds = PyMod_radioselect_qt(label_text="Exclude Gaps",
                                                        buttons=('Yes', 'No'))
            self.exclude_gaps_rds.setvalue("No")
            self.middle_formlayout.add_widget_to_align(self.exclude_gaps_rds)

        self.middle_formlayout.set_input_widgets_width("auto", padding=10)

"""

############END Original CODE###############################



##############Add ETE code 21/05/2026-SARA###############

# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.



import os
import subprocess
import sys
import shutil
import re
import time
import glob
import colorsys
import random

from pymod_lib import pymod_vars
from pymod_lib.pymod_os_specific import pymod_which, get_exe_file_name
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_combobox_qt,
                                                          PyMod_hbox_option_qt,
                                                          askyesno_qt,
                                                          askopenfile_qt)

from ._evolutionary_analysis_base import (Evolutionary_analysis_protocol,
                                          is_a_nucleic_acid_alignment)
from ._tree_visualization import bubble_ly, ColumnMutation, layout


from Bio import Phylo
from Bio import AlignIO, SeqIO
from Bio.Phylo import TreeConstruction
from Bio.Phylo.Consensus import *
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


from ete3 import Tree, PhyloTree, AttrFace, TreeFace, TreeStyle, TextFace, NodeStyle, faces

# Import the Qt bindings through 'pymol.Qt', like the rest of PyMod, instead of
# importing PyQt5 directly. NOTE: the widget names below live in 'QtWidgets'.
# A previous version of this module imported some of them from 'PyQt5.QtGui',
# where they do not exist: that only worked because 'pymod_lib.pymod_plot.pyqtgraph'
# happens to be imported earlier and copies the 'QtWidgets' names into 'QtGui'
# for Qt4 compatibility.
from pymol.Qt import QtWidgets, QtCore

QObject = QtCore.QObject
QThread = QtCore.QThread
pyqtSignal = QtCore.pyqtSignal
pyqtSlot = QtCore.pyqtSlot

QWidget = QtWidgets.QWidget
QProgressBar = QtWidgets.QProgressBar
QPushButton = QtWidgets.QPushButton
QVBoxLayout = QtWidgets.QVBoxLayout
QDesktopWidget = QtWidgets.QDesktopWidget

tree_building_alg_dict = {"Neighbor Joining": "NJ", "UPGMA": "UPGMA", "Parsimony" : "Parsimony", "ML (RAxML)": "ML_raxml", "ML (IQ-TREE 2)": "ML_iqtree"}

calculator_matrix_dict = {"Identity" : False, "Blosum62" : False} #previous "Identity" : True , Common ancestor not working with identity clustering algorithm 

# Neither RAxML nor IQ-TREE2 will bootstrap fewer than this many sequences.
MIN_BOOTSTRAP_SEQUENCES = 4

# The substitution models RAxML is told to use. RAxML has to be told which kind
# of sequences it is reading: given an amino acid model, a nucleotide alignment
# runs to the end and comes back as a tree, because A, C, G and T are amino acid
# codes as well, but the rates behind that tree are those of protein evolution
# and the tree means nothing.
RAXML_PROTEIN_MODEL = "PROTGAMMAJTT"
RAXML_NUCLEOTIDE_MODEL = "GTRGAMMA"


def raxml_model_of(cluster_element):
    """
    The RAxML substitution model for the sequences of an alignment.
    """
    if is_a_nucleic_acid_alignment(cluster_element):
        return RAXML_NUCLEOTIDE_MODEL
    return RAXML_PROTEIN_MODEL

visualization_dict = {"Alignment" : "aln", "Common Ancestor" : "ancestors","Basic" : "bp", "Evolutionary events" : "evo",  "Mutation Events": "mut_all","Single Mutation Event": "mut_single", "Bootstrap Analysis" : "bootstrap"}

random_trees_dict_RAxML = {"100":100, "500":500, "1000": 1000}

random_trees_dict_IQTREE = {"1000":1000, "5000":5000, "10000": 10000}

boot_program_dict = {"IQ-TREE2":"iqtree2", "RAxML":"raxml"}



class Tree_building(Evolutionary_analysis_protocol):

    def launch_from_gui(self):
        """
        It will check if a software to build a tree is available on the user's machine.
        """
        ####################
        # Biopython module #
        ####################

        self.tree_building_software = None
        can_build_tree = False
        if self.check_tree_constructor_module():
            self.tree_building_software = "biopython"
            can_build_tree = True

        if can_build_tree:
            self.build_tree_building_window()
        else:
            title = "Tree building Error"
            message = "In order to build a tree out of an alignment you need to install Biopython."
            self.pymod.main_window.show_error_message(title, message)

    def check_tree_constructor_module(self):
        try:
            import Bio.Phylo.TreeConstruction
            return True
        except ImportError:
            return False

    def build_tree_building_window(self):
        """
        Builds a window with options to build a tree out of an alignment.
        """
        self.tree_building_window = Tree_building_window_qt(self.pymod.main_window,
            protocol=self,
            title="Options for Tree Building",
            upper_frame_title="Here you can modify options for Tree Building",
            submit_command=self.run_tree_building_software)
        self.tree_building_window.show()


    def run_tree_building_software(self):

        # Saves a temporary input alignment file.
        alignment_file_name = "alignment_tmp"
        self.alignment_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.fasta')
        self.pymod.save_alignment_fasta_file(alignment_file_name, self.input_cluster_element.get_children())
        self.pymod.tree_aln_id = self.input_cluster_element.cluster_id

        # Get the parameters from the GUI.
        self.clustering_algorithm = self.get_clustering_algorithm()

        # Convert id FASTA file in phylip format.
        # alignment file.
        self.new_aln_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.phylip')
        self.edit_FASTA_file(self.alignment_file_path, self.new_aln_file_path)
        with open(self.new_aln_file_path) as aln_fasta_handle:
            aln = AlignIO.read(aln_fasta_handle, 'fasta')

        # RAxML answers "TOO FEW SPECIES" and IQ-TREE2 "It makes no sense to
        # perform bootstrap with less than 4 sequences", both after being started
        # and both by way of a non-zero exit code. Asked here instead.
        if (self.get_visualization_mode() == "bootstrap"
                and len(aln) < MIN_BOOTSTRAP_SEQUENCES):
            self.pymod.main_window.show_error_message(
                "Bootstrap Analysis",
                "A bootstrap analysis needs at least %s sequences, and this"
                " alignment has %s. Build the tree with another view, or align"
                " more sequences first."
                % (MIN_BOOTSTRAP_SEQUENCES, len(aln)))
            return None

        # Convert FASTA to PHYLIP format useful for bootstrapping.
        self.new_aln_file_path_phy = os.path.join(self.pymod.alignments_dirpath, "alignment_tmp_boot.phylip")
        with open(self.new_aln_file_path) as fasta_handle:
            records = list(SeqIO.parse(fasta_handle, "fasta"))
        SeqIO.write(records, self.new_aln_file_path_phy, "phylip")

        # A distance read from an amino acid substitution matrix says nothing
        # about nucleotides: the matrix has no entry for a pair of bases, it
        # has one for the amino acids which happen to be written with the same
        # letters, so a tree comes out of it which means nothing.
        if (self.clustering_algorithm in ("NJ", "UPGMA")
                and self.set_calculator_value()[1] == "Blosum62"
                and is_a_nucleic_acid_alignment(self.input_cluster_element)):
            self.pymod.main_window.show_error_message(
                "Tree building Error",
                "The Blosum62 distance is not available for nucleic acid"
                " alignments: it is a matrix of amino acid substitutions."
                " Build the tree with a maximum likelihood method instead.")
            return None

        if self.clustering_algorithm == "Parsimony":
            # BioPython TreeConstruction module.
            scorer = TreeConstruction.ParsimonyScorer()
            searcher = TreeConstruction.NNITreeSearcher(scorer)
            constructor = TreeConstruction.ParsimonyTreeConstructor(searcher)
            pars_tree = constructor.build_tree(aln)

            # BioPhylo tree need to be converted in newick format.
            new_tree_file_name = "%s_%s_align_tree.phy" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
            self.new_tree_file_path = os.path.join(self.pymod.alignments_dirpath, new_tree_file_name)

            self.new_tree_file_name_output_file = "%s_%s_tree_output_PARSIMONY.nwk" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
            self.new_tree_file_path_nwk = os.path.join(self.pymod.alignments_dirpath, self.new_tree_file_name_output_file)

            Phylo.write(pars_tree, self.new_tree_file_path, "newick")
            self.remove_Inner_label(self.new_tree_file_path,self.new_tree_file_path_nwk)

            # Bootstrap analysis requirement.
            self.to_be_unrooted = True


        elif self.clustering_algorithm == "NJ":
            # BioPython DistanceTreeConstructor module.
            calculator_type = self.set_calculator_value()[1]
            calculator = self.set_calculator_value()[0]
            constructor = DistanceTreeConstructor(calculator, 'nj')
            nj_tree = constructor.build_tree(aln)

            # BioPhylo tree need to be converted in newick format.
            new_tree_file_name = "%s_%s_align_tree.phy" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
            self.new_tree_file_path = os.path.join(self.pymod.alignments_dirpath, new_tree_file_name)

            self.new_tree_file_name_output_file = "%s_%s_tree_output_NJ.nwk" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
            self.new_tree_file_path_nwk = os.path.join(self.pymod.alignments_dirpath, self.new_tree_file_name_output_file)

            Phylo.write(nj_tree, self.new_tree_file_path, "newick")
            self.remove_Inner_label(self.new_tree_file_path,self.new_tree_file_path_nwk)

            # Bootstrap analysis requirement.
            self.to_be_unrooted = False

        elif self.clustering_algorithm == "UPGMA":
            # BioPython DistanceTreeConstructor module.
            calculator_type = self.set_calculator_value()[1]
            calculator = self.set_calculator_value()[0]
            constructor = DistanceTreeConstructor(calculator, 'upgma')
            upgma_tree = constructor.build_tree(aln)

            new_tree_file_name = "%s_%s_align_tree.phy" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
            self.new_tree_file_path = os.path.join(self.pymod.alignments_dirpath, new_tree_file_name)

            # BioPhylo tree need to be converted in newick format.
            self.new_tree_file_name_output_file = "%s_%s_tree_output_UPGMA.nwk" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
            self.new_tree_file_path_nwk = os.path.join(self.pymod.alignments_dirpath, self.new_tree_file_name_output_file)

            Phylo.write(upgma_tree, self.new_tree_file_path, "newick")
            self.remove_Inner_label(self.new_tree_file_path, self.new_tree_file_path_nwk)

            # Bootstrap analysis requirement.
            self.to_be_unrooted = True

        elif self.clustering_algorithm == "ML_raxml":
            try:
                # RAxML does not work when run on sessions saved before its download.
                old_session = True
                self.tool = self.pymod.raxml

                # Output file name.
                self.new_tree_file_name_output_file = "%s_%s_tree_output_ML_RAxML.nwk" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
                self.new_tree_file_path_nwk = os.path.join(self.pymod.alignments_dirpath, self.new_tree_file_name_output_file)

                # RaXML executable file path.
                self.raxml_PTHREADS_AVX_exe_filepath = os.path.join(self.pymod.raxml["exe_dir_path"].get_value(), get_exe_file_name("raxmlHPC_PTHREADS_AVX"))
                self.iqtree2_exe_file_path = os.path.join(self.pymod.iqtree2["exe_file_path"].get_value())
                old_session = False

                # args for PopUpProgressB.
                self.clustering_algorithm = self.get_clustering_algorithm()
                self.visualization = self.get_visualization_mode()
                self.to_be_unrooted = False
                self.number_random_trees = None
                self.bootstrap_program = None
                try:
                    self.bootstrap_program = self.get_bootstrap_program()
                    if self.bootstrap_program == "raxml":
                        self.number_random_trees = self.get_random_trees_raxml()
                    else:
                        self.number_random_trees = self.get_random_trees_iqtree()
                except:
                    pass
                # If the tree has already been found (the command has already been run) there is no need to call the class with the method to run the program again.
                if not os.path.isfile(self.new_tree_file_path_nwk) or self.visualization == "bootstrap":
                    # This class calls a Worker class that starts the thread (RAxML Maximum-likelihood-Tree search) and the interactive progress bar.
                    self.popup_raxml = PopUpProgressB(pymod = self,
                                                show_tree = self.pymod.show_tree,
                                                main_window = self.pymod.main_window,
                                                alignments_files_names = self.pymod.alignments_files_names,
                                                input_cluster_element = self.input_cluster_element.unique_index,
                                                alignment_file_path_bootstrap = self.new_aln_file_path_phy,
                                                alignment_file_path_best_tree = self.alignment_file_path,
                                                alignment_best_tree = self.new_aln_file_path,
                                                alignments_dirpath =self.pymod.alignments_dirpath,
                                                tree_file_path = None,
                                                tree_id = self.new_tree_file_name_output_file,
                                                tree = None,
                                                num_random_trees = self.number_random_trees,
                                                clustering_algorithm = self.clustering_algorithm,
                                                bootstrap_program = self.bootstrap_program,
                                                raxml_exec = self.raxml_PTHREADS_AVX_exe_filepath,
                                                iqtree_exec = self.iqtree2_exe_file_path,
                                                tree_building_window=self.tree_building_window,
                                                visualization = self.visualization,
                                                new_aln_file_path = self.new_aln_file_path,
                                                to_be_unrooted = self.to_be_unrooted)
                    if self.visualization == "bootstrap":
                        # PopUpProgressB should wait until the user chooses a bootstrap_program.
                        if self.bootstrap_program == None:
                            self.pymod.main_window.show_warning_message("Bootstrap Analysis","Please select a Bootstrap Program")
                        else:self.popup_raxml.start_progress()

                    else:self.popup_raxml.start_progress()
                else:pass
            except:
                if old_session:
                    self.pymod.main_window.show_error_message("Program not found", "RAxML does not work when run on sessions saved before its download.")
                else:
                    self.pymod.main_window.show_error_message("Program not found", "RAxML was not found. In order to run it download the external tools via the 'HELP' --> 'PyMod components'.")


        elif self.clustering_algorithm == "ML_iqtree":
            try:
                # IQ-TREE2 does not work when run on sessions saved before its download.
                old_session = True
                self.tool = self.pymod.iqtree2

                # Output file name.
                self.new_tree_file_name_output_file = "%s_%s_tree_output_ML_IQTREE2.nwk" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index)
                self.new_tree_file_path_nwk = os.path.join(self.pymod.alignments_dirpath, self.new_tree_file_name_output_file)

                # Convert FASTA to PHYLIP format useful for bootstrapping.
                self.new_aln_file_path_phy = os.path.join(self.pymod.alignments_dirpath, "alignment_tmp_boot.phylip")
                with open(self.new_aln_file_path) as fasta_handle:
                    records = list(SeqIO.parse(fasta_handle, "fasta"))
                SeqIO.write(records, self.new_aln_file_path_phy, "phylip")

                #args for PopUpProgressB
                # RaXML executable file path.
                self.raxml_PTHREADS_AVX_exe_filepath = os.path.join(self.pymod.raxml["exe_dir_path"].get_value(), get_exe_file_name("raxmlHPC_PTHREADS_AVX"))
                # IQ-TREE 2 executable file path.
                self.iqtree2_exe_file_path = os.path.join(self.pymod.iqtree2["exe_file_path"].get_value())
                old_session = False
                self.clustering_algorithm = self.get_clustering_algorithm()
                self.visualization = self.get_visualization_mode()
                self.to_be_unrooted = False
                self.number_random_trees = None
                self.bootstrap_program = None
                try:
                    self.bootstrap_program = self.get_bootstrap_program()
                    if self.bootstrap_program == "raxml":
                        self.number_random_trees = self.get_random_trees_raxml()
                    else:
                        self.number_random_trees = self.get_random_trees_iqtree()
                except:
                    pass

                # If the tree has already been found (the command has already been run) there is no need to call the class with the method to run the program again.
                if not os.path.isfile(self.new_tree_file_path_nwk) or self.visualization == "bootstrap":
                    # This class calls a Worker class that starts the thread (IQTREE2 Maximum-likelihood-Tree search) and the interactive progress bar.
                    self.popup_iqtree = PopUpProgressB(pymod = self,
                                            show_tree = self.pymod.show_tree,
                                            main_window = self.pymod.main_window,
                                            alignments_files_names = self.pymod.alignments_files_names,
                                            input_cluster_element = self.input_cluster_element.unique_index,
                                            alignment_file_path_bootstrap = self.new_aln_file_path_phy,
                                            alignment_file_path_best_tree = self.alignment_file_path,
                                            alignment_best_tree = self.new_aln_file_path,
                                            alignments_dirpath =self.pymod.alignments_dirpath,
                                            tree_file_path = None,
                                            tree_id = self.new_tree_file_name_output_file,
                                            tree = None,
                                            num_random_trees = self.number_random_trees,
                                            clustering_algorithm = self.clustering_algorithm,
                                            bootstrap_program = self.bootstrap_program,
                                            raxml_exec = self.raxml_PTHREADS_AVX_exe_filepath,
                                            iqtree_exec = self.iqtree2_exe_file_path,
                                            tree_building_window=self.tree_building_window,
                                            visualization = self.visualization,
                                            new_aln_file_path = self.new_aln_file_path,
                                            to_be_unrooted = self.to_be_unrooted)
                    if self.visualization == "bootstrap":
                        # PopUpProgressB should wait until the user chooses a bootstrap_program.
                        if self.bootstrap_program == None:
                            answer = self.pymod.main_window.show_warning_message("Bootstrap Analysis","Please select a Bootstrap Program")

                        else:self.popup_iqtree.start_progress()

                    else:self.popup_iqtree.start_progress()

                else:pass
            except:
                if old_session:
                    self.pymod.main_window.show_error_message("Program not found", "IQ-TREE2 does not work when run on sessions saved before its download.")
                else:
                    self.pymod.main_window.show_error_message("Program not found", "IQ-TREE2 was not found. In order to run it download the external tools via the 'HELP' --> 'PyMod components'.")

        # A maximum likelihood run happens in a thread of its own, and shows its
        # tree when that thread finishes. Asking for it here as well means asking
        # for a file which the search has not written yet, which is why every such
        # run used to raise FileNotFoundError at this line -- swallowed, until the
        # failures started being reported. The exception is a tree which was built
        # before: the thread finds the file and shows nothing, so this does.
        if self.should_show_tree_now():
            try:
                self.show_plot_window(self.new_tree_file_path_nwk)
            except Exception as e:
                self.report_visualization_failure(e)

        #Remove temporary files.
        #os.remove(self.alignment_file_path)
        #self.tree_building_window.destroy()


        #######################################################################
        # Tree dictionary to be appended to the alignment clusters list
        try:
            tree = {"algorithm": [self.clustering_algorithm, calculator_type],
                    "filepath": self.new_tree_file_path_nwk}

        # Parsimony has not calculator option.
        except UnboundLocalError:
            tree = {"algorithm": self.clustering_algorithm ,
                    "filepath": self.new_tree_file_path_nwk}


        # Set "trees_list" attribute for old saved Pymod sessions.
        if not hasattr(self.input_cluster_element, "trees_list"):
            setattr(self.input_cluster_element, 'trees_list',[])

        self.input_cluster_element.trees_list.append(tree)
        ######################################################################

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                   update_menus=True, update_elements=True)

   
   ##############MODIFICATION on 31/07/2025################
    #Maybe put it outside any class or leave it here?


    ##############MODIFICATION SARA 11/04/2026###########################


    ########Check function ###############
    def is_tree_unrooted(self, tree):
        return tree.get_tree_root().is_leaf() is False and len(tree.children) > 2

    def handle_non_rooted_tree_error(self, tree):

        reply = QtWidgets.QMessageBox.question(
            self.pymod.get_qt_parent(),
            "Root Tree Automatically?",
            "The tree is currently unrooted.\n"
            "Do you want to root it using midpoint-like rooting?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.Yes
        )

        if reply != QtWidgets.QMessageBox.Yes:
            return None

        try:
            leaves = tree.get_leaves()

            if len(leaves) < 2:
                return tree

            # find most distant pair
            max_d = -1
            best_pair = (leaves[0], leaves[1])

            for i in range(len(leaves)):
                for j in range(i + 1, len(leaves)):
                    d = tree.get_distance(leaves[i], leaves[j])
                    if d > max_d:
                        max_d = d
                        best_pair = (leaves[i], leaves[j])

            # ALWAYS ROOT WITH A LEAF (robust)
            tree.set_outgroup(best_pair[0])

            return tree

        except Exception as e:
            self.pymod.main_window.show_error_message(
                "Rooting error",
                str(e)
            )
            return None



###################END####################


    """
    def handle_non_rooted_tree_error(self):
        
        #Prompts the user with a dialog asking whether to automatically root the tree 
        #using midpoint rooting if the tree is not rooted.

        #If the user agrees, attempts to root the tree and return its root.
        #If rooting fails, shows an error message.
        #If the user declines, returns None to indicate no rooting was performed.

        #TODO: Implement the actual midpoint rooting function (e.g., self.root_tree_midpoint()).
        

        # Show a message box asking the user whether to root the tree automatically
        reply = QtWidgets.QMessageBox.question(
            self.pymod.get_qt_parent(),  # Parent widget for the dialog
            "Root Tree Automatically?",  # Dialog title
            "The tree must be rooted to perform this analysis.\n"
            "Would you like to root it automatically using midpoint rooting?",  # Message text
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,  # Buttons: Yes and No
            QtWidgets.QMessageBox.Yes  # Default button: Yes
        )

        if reply == QtWidgets.QMessageBox.Yes:
            # User agreed to root the tree automatically
            # TODO: Call the function that performs midpoint rooting here, e.g.:
            # self.root_tree_midpoint()
            
            try:
                # Try to return the root of the now rooted tree
                return self.tree.root
            except TypeError:
                # If it still fails, show an error message to the user
                self.pymod.main_window.show_error_message(
                    "Error", "The Tree is still not rooted after midpoint rooting."
                )
                return None
        else:
            # User declined to root the tree; return None
            return None
    #################END######################

    """

    def show_plot_window(self,tree_file_path):
          self.visualization = self.get_visualization_mode()

          if self.clustering_algorithm == "ML_raxml" or self.clustering_algorithm == "ML_iqtree":
              tree_file_path = self.new_tree_file_path_nwk

          if self.visualization == "bp":
              self.pymod.show_tree(tree_file_path)

          elif self.visualization == "aln":
              # Visualization with alignment file.
              ts = TreeStyle()
              tree = PhyloTree(tree_file_path, format = 1)
              with open(self.alignment_file_path,'r') as f:
                  aln = f.read()
              tree.link_to_alignment(aln)

              tree.show(tree_style=ts)

          elif self.visualization == "mut_single":
              # Visualization with single mutation event and alignment file.
              tree = PhyloTree(tree_file_path, format = 1)
              with open(self.new_aln_file_path,'r') as f:
                  aln = f.read()
              tree_params = ColumnMutation()
              t, ts = tree_params.get_mutation_event(tree, aln, True)

              t.show(tree_style=ts)

          elif self.visualization == "mut_all":
              # Visualization with all mutation events and alignment file.
              tree = PhyloTree(tree_file_path, format = 1)
              with open(self.new_aln_file_path,'r') as f:
                  aln = f.read()
              tree_params = ColumnMutation()
              t, ts = tree_params.get_mutation_event(tree, aln, False)
              tree.show(tree_style = ts)
          
          ######ALIGNEMTN section##########
        
          elif self.visualization == "ancestors":
                
                ######OLD code 
                tree = PhyloTree(tree_file_path, format = 1)
                not_ancestor_list = []
                for n in tree.traverse():
                    if n.is_leaf():
                        sister_node = n.get_sisters()
                        for s in sister_node:
                            try:
                                ancestor = tree.get_common_ancestor(n.name,s.name)
                                nst = NodeStyle()
                                # Generate random colors.
                                h,s,l = random.random(), 0.5 + random.random()/2.0, 0.4 + random.random()/5.0
                                r= [int(256*i) for i in colorsys.hls_to_rgb(h,l,s)]
                                color = '#%02x%02x%02x' % (r[0],r[1],r[2])
                                nst["bgcolor"] = color
                                ancestor.set_style(nst)
                            except:
                                not_ancestor_list.append(n.name)

                ts = TreeStyle()
                ts.layout_fn = layout
                ts.show_leaf_name = False
                ts.mode = "c"
                ts.root_opening_factor = 10
                tree.show(tree_style = ts)

                """NEW
               # Visualization showing common ancestors.

                tree = PhyloTree(tree_file_path, format=1)

                # ROOT FIRST
                tree = Tree_building.handle_non_rooted_tree_error(self, tree)

                if tree is None:
                    return

                # ORIGINAL LOGIC
                not_ancestor_list = []
                for n in tree.traverse():
                    if n.is_leaf():
                        sister_node = n.get_sisters()
                        for s in sister_node:
                            try:
                                ancestor = tree.get_common_ancestor(n.name, s.name)
                                nst = NodeStyle()

                                # Generate random colors
                                h, s_val, l = random.random(), 0.5 + random.random()/2.0, 0.4 + random.random()/5.0
                                r = [int(256*i) for i in colorsys.hls_to_rgb(h, l, s_val)]
                                color = '#%02x%02x%02x' % (r[0], r[1], r[2])

                                nst["bgcolor"] = color
                                ancestor.set_style(nst)

                            except:
                                not_ancestor_list.append(n.name)

                ts = TreeStyle()
                ts.layout_fn = layout
                ts.show_leaf_name = False
                ts.mode = "c"
                ts.root_opening_factor = 10

                tree.show(tree_style=ts)

                """

                """OLD code CONNECTED TO evo 
                # Visualization showing common ancestors.
                try:
                    tree = PhyloTree(tree_file_path, format = 1)
                    not_ancestor_list = []
                    for n in tree.traverse():
                        if n.is_leaf():
                            sister_node = n.get_sisters()
                            for s in sister_node:
                                try:
                                    ancestor = tree.get_common_ancestor(n.name,s.name)
                                    nst = NodeStyle()
                                    # Generate random colors.
                                    h,s,l = random.random(), 0.5 + random.random()/2.0, 0.4 + random.random()/5.0
                                    r= [int(256*i) for i in colorsys.hls_to_rgb(h,l,s)]
                                    color = '#%02x%02x%02x' % (r[0],r[1],r[2])
                                    nst["bgcolor"] = color
                                    ancestor.set_style(nst)
                                except:
                                    not_ancestor_list.append(n.name)

                    ts = TreeStyle()
                    ts.layout_fn = layout
                    ts.show_leaf_name = False
                    ts.mode = "c"
                    ts.root_opening_factor = 10
                    tree.show(tree_style = ts)

                except TypeError:
                    title = "Evolutionary Events"
                    message = "The Tree is not rooted."
                    
                    # Show an initial error message indicating that the tree is not rooted
                    self.pymod.main_window.show_error_message(title, message)
                    
                    #### OLD BEHAVIOR
                    # Previously: silently passed or returned None, which could lead to misleading success messages
                    # pass
                    # return None

                    #### MODIFICATION on 31/07/2025 ###########
                    # Now we explicitly call a function to ask the user whether they want to root the tree automatically.
                    # If the user agrees, the tree will be midpoint-rooted and the root will be retrieved.
                    # If the user declines or if rooting fails, the function handles that internally.
                    ###########################################
                    
                    #Tree_building.handle_non_rooted_tree_error()  # Call the handler to prompt user interaction and attempt recovery

                    ##############MODIFICATION SARA 11/04/2026##################
                    tree = PhyloTree(self.input_cluster_element, format=1)

                    # FIX 1: define aln
                    aln = None

                    # FIX 2: try to load alignment
                    try:
                        with open(self.alignment_file_path, 'r') as f:
                            aln = f.read()
                    except FileNotFoundError:
                        pass  # or your dialog logic

                    # ROOTING
                    tree = Tree_building.handle_non_rooted_tree_error(self, tree)

                    if tree is None:
                        return

                    # ALIGNMENT CHECK
                    if aln is None:
                        return

                    # ANALYSIS
                    tree_params = ColumnMutation()
                    t, ts = tree_params.get_mutation_event(tree, aln, True)

                    t.show(tree_style=ts)

                    
                    ########END#################
                """
          elif self.visualization == "evo":
              
             
                tree = PhyloTree(tree_file_path, format=1)
                
                #############MODIFICATION SARA on 12/04/2026#################
                # ONLY ASK IF REALLY UNROOTED
                if self.is_tree_unrooted(tree):
                    tree = Tree_building.handle_non_rooted_tree_error(self, tree)

                    if tree is None:
                        return

                # DEBUG
                print("Tree loaded and rooted")

                # SAFE CALL (no assumption about output)
                try:
                    _ = tree.get_descendant_evol_events()
                except Exception as e:
                    print("EVO FAILED:", e)

                print("EVO COMPUTED")

                ts = TreeStyle()
                ts.show_leaf_name = True
                ts.mode = "c"

                print("SHOWING TREE")

                tree.show(tree_style=ts)

                ##################END ######################
                        
                """OLD MODIFCATON 
                # Visualization showing duplication and speciation events.
                tree = PhyloTree(tree_file_path, format = 1)
                try:
                    ts = TreeStyle()
                    # Get evolutionary events for each internal node (Duplication or Speciation).
                    events = tree.get_descendant_evol_events()
                    i = 0
                    for n in tree.traverse():
                        if not n.is_leaf():
                            ev = events[i]
                            if ev.etype == "S":
                                # Add node face.
                                n.add_face(TextFace("S",fgcolor = "blue", penwidth = 8), column = i, position = "branch-top")
                            elif ev.etype == "D":
                                # Add node face.
                                n.add_face(TextFace("D", fgcolor = "red", penwidth = 8), column = i, position = "branch-top")
                            i += 1
                    ts.mode = "c"
                    ts.show_leaf_name = True

                    tree.show(tree_style=ts)

                except TypeError:
                    title = "Evolutionary Events"
                    message = "The Tree is not rooted."
                    
                    # Show an initial error message indicating that the tree is not rooted
                    self.pymod.main_window.show_error_message(title, message)
                    
                    #### OLD BEHAVIOR
                    # Previously: silently passed or returned None, which could lead to misleading success messages
                    # pass
                    # return None

                    #### MODIFICATION on 31/07/2025 ###########
                    # Now we explicitly call a function to ask the user whether they want to root the tree automatically.
                    # If the user agrees, the tree will be midpoint-rooted and the root will be retrieved.
                    # If the user declines or if rooting fails, the function handles that internally.
                    ###########################################
                    
                    #Tree_building.handle_non_rooted_tree_error()  # Call the handler to prompt user interaction and attempt recovery
                    ###############END###############

                    ################MODIFICATION SARA 11/04/2026####################
                    tree = PhyloTree(self.input_cluster_element, format=1)

                    # FIX 1: define aln
                    aln = None

                    # FIX 2: try to load alignment
                    try:
                        with open(self.alignment_file_path, 'r') as f:
                            aln = f.read()
                    except FileNotFoundError:
                        pass  # or your dialog logic

                    # ROOTING
                    tree = Tree_building.handle_non_rooted_tree_error(self, tree)

                    if tree is None:
                        return

                    # ALIGNMENT CHECK
                    if aln is None:
                        return

                    # ANALYSIS
                    tree_params = ColumnMutation()
                    t, ts = tree_params.get_mutation_event(tree, aln, True)

                    t.show(tree_style=ts)


                    ###############END ##################
                    """
          elif self.visualization == "bootstrap":
              # Visualization showing bootstrap values.
              try:
                  if self.clustering_algorithm!="ML_raxml" and self.clustering_algorithm!="ML_iqtree":
                      # IQ-TREE2 & RAxML do not work when run on sessions saved before their download.
                      old_session = False
                      # Get bootstrap_program.
                      self.bootstrap_program = self.get_bootstrap_program()
                      old_session = True
                      # Get User input number of random trees.
                      if self.bootstrap_program == "raxml":
                          self.number_random_trees = self.get_random_trees_raxml()
                      else:
                          self.number_random_trees = self.get_random_trees_iqtree()

                      title = "Bootstrap Analysis"
                      message = ("The Bootstrap Analysis will be run via %s and it will produce %s random trees and the process may take some time depending on the size of your multiple alignment. Do you want to proceed?"
                                  % (self.bootstrap_program ,self.number_random_trees))
                      answer = askyesno_qt(title, message)
                      if answer:
                          tree = PhyloTree(tree_file_path, format = 1)

                          # UPGMA and Parsimony generate rooted trees which are not processable by RAxML (bootstrap analysis).
                          if self.to_be_unrooted:
                              tree.unroot() # converts a split into multifurcation
                              # Save the unrooted tree in a new file.
                              tree_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_%s_unrooted_tree.nwk" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index))
                              tree.write(format = 1, outfile = tree_file_path)
                          else:
                              pass

                          # RAxML executable file path.
                          self.raxml_PTHREADS_AVX_exe_filepath = os.path.join(self.pymod.raxml["exe_dir_path"].get_value(), get_exe_file_name("raxmlHPC_PTHREADS_AVX"))
                           # IQ-TREE 2 executable file path.
                          self.iqtree2_exe_file_path = os.path.join(self.pymod.iqtree2["exe_file_path"].get_value())
                          old_session = False

                          if os.path.isfile(self.raxml_PTHREADS_AVX_exe_filepath) and os.path.isfile(self.iqtree2_exe_file_path):
                              # This class calls a Worker class that starts the thread and the interactive progress bar.
                              self.popup = PopUpProgressB(pymod = self,
                                                          show_tree = None,
                                                          main_window = self.pymod.main_window,
                                                          alignments_files_names = self.pymod.alignments_files_names,
                                                          input_cluster_element = self.input_cluster_element.unique_index,
                                                          alignment_file_path_bootstrap = self.new_aln_file_path_phy,
                                                          alignment_file_path_best_tree = None,
                                                          alignment_best_tree = None,
                                                          alignments_dirpath =self.pymod.alignments_dirpath,
                                                          tree_file_path = tree_file_path,
                                                          tree_id = self.new_tree_file_name_output_file,
                                                          tree = tree,
                                                          num_random_trees = self.number_random_trees,
                                                          clustering_algorithm = None,
                                                          bootstrap_program = self.bootstrap_program,
                                                          raxml_exec = self.raxml_PTHREADS_AVX_exe_filepath,
                                                          iqtree_exec = self.iqtree2_exe_file_path,
                                                          tree_building_window = None,
                                                          visualization = self.visualization,
                                                          new_aln_file_path = None,
                                                          to_be_unrooted = None)


                              self.popup.start_progress()

                          else:
                              self.pymod.main_window.show_error_message("Program not found", "The bootstrap program was not found. In order to run it download the external tools via the 'HELP' --> 'PyMod components'.")
                              #pass
                              ######MODIFICATION on 31/07/2025###########
                              # Changed `pass` to `return None` to prevent a misleading message ("Tree has been imported successfully")
                              # from appearing even when the tree was not correctly processed.
                              return None
                              ########END#################
                      else:
                          #pass
                          ######MODIFICATION on 31/07/2025###########
                          # Changed `pass` to `return None` to prevent a misleading message ("Tree has been imported successfully")
                          # from appearing even when the tree was not correctly processed.
                          return None
                          ########END#################
                  # Bootstrap Analysis is called in the PopUpProgressB class.
                  else:pass
              except:
                  if old_session:
                      self.pymod.main_window.show_warning_message("Bootstrap Analysis","The bootstrap program does not work when run on sessions saved before its download.")
                  else: self.pymod.main_window.show_warning_message("Bootstrap Analysis","Please select a Bootstrap Program")


    def should_show_tree_now(self):
        """
        Whether the tree is to be shown from here, rather than by the thread of a
        maximum likelihood search when it finishes.

        Everything but RAxML and IQ-TREE2 builds its tree in this thread, and the
        file is on disk by now. Those two run in a thread and emit a signal when
        they are done -- unless the tree was built before and is still on disk, in
        which case the thread shows nothing and this is what shows it.
        """

        if self.clustering_algorithm not in ("ML_raxml", "ML_iqtree"):
            return True

        path = getattr(self, "new_tree_file_path_nwk", None)
        return bool(path) and os.path.isfile(path)


    def report_visualization_failure(self, error):
        """
        Says that the tree could not be shown, and why.

        The tree itself was built: what failed is the drawing of it. This used to
        be an 'except: pass', so a mode which could not open left the user with a
        finished run and nothing on screen and no message.
        """

        import traceback

        mode = self.get_visualization_mode()
        label = mode
        for name, code in visualization_dict.items():
            if code == mode:
                label = name
                break

        print("- The tree was built, but the '%s' view could not be shown:" % label)
        traceback.print_exc()

        self.pymod.main_window.show_error_message(
            "Tree Visualization Error",
            "The tree was built, but it could not be shown in the '%s' view:\n\n"
            "%s: %s\n\nThe tree itself is under the 'Trees' menu, and its 'Basic'"
            " view is drawn by PyMod rather than by the viewer of ETE."
            % (label, type(error).__name__, error))



    def edit_FASTA_file(self, input_file, output_file):
        '''
        Change the record id in a FASTA file using BioPython.
        '''
        with open(input_file) as original, open(output_file, 'w') as corrected:
            records = SeqIO.parse(original, 'fasta')
            # Duplicate names found errors.
            bad_chars = ['_', 'c', 'h', 'a', 'i', 'n'] #"_chain_"
            for record in records:
                if len(record.id) > 10:
                    if "chain" in record.id:
                        for i in bad_chars :
                            record.id = record.id.replace(i, '')
                    record.id = record.id[0:10]
                    record.description = record.id[0:10]
                SeqIO.write(record, corrected, 'fasta')

    def remove_Inner_label(self,input_file,output_file):
        '''
        Remove "Inner" node label from Phylo.write output.
        '''
        input = open(input_file, "r")
        output = open(output_file, "w")
        output.write(input.read().replace("Inner",""))
        input.close()
        output.close()

    def get_visualization_mode(self):
        return visualization_dict[self.tree_building_window.visualization_rds.getvalue()]

    def get_clustering_algorithm(self):
        return tree_building_alg_dict[self.tree_building_window.algorithm_rds.getvalue()]

    def get_calculator(self):
        return calculator_matrix_dict[self.tree_building_window.calculator_rds.getvalue()]

    def get_random_trees_raxml(self):
        return random_trees_dict_RAxML[self.tree_building_window.random_trees_rds_raxml.getvalue()]

    def get_random_trees_iqtree(self):
        return random_trees_dict_IQTREE[self.tree_building_window.random_trees_rds_iqtree.getvalue()]

    def get_bootstrap_program(self):
        return boot_program_dict[self.tree_building_window.boot_program_rds.getvalue()]

    def set_calculator_value(self):
        if self.get_calculator():
            calculator = DistanceCalculator('identity')
            calculator_type = "Identity"
        else:
            calculator = DistanceCalculator('blosum62')
            calculator_type = "Blosum62"
        return calculator, calculator_type


class PopUpProgressB(QWidget):

    is_pymod_window = True

    def __init__(self,pymod, show_tree, main_window, alignments_files_names,
                    input_cluster_element, alignment_file_path_bootstrap, alignment_file_path_best_tree, alignment_best_tree, alignments_dirpath, tree_file_path, tree_id, tree,
                     num_random_trees, clustering_algorithm, bootstrap_program, raxml_exec, iqtree_exec, tree_building_window, visualization, new_aln_file_path,to_be_unrooted):

        # Parent this window to the PyMod main window. A top level widget with no
        # parent would not inherit the stylesheet which PyMOL sets on its main
        # window, and this progressbar would look different from the ones of the
        # other PyMod protocols. 'QtCore.Qt.Window' keeps it a separate window
        # (a parented QWidget would otherwise be drawn inside its parent).
        super().__init__(main_window, QtCore.Qt.Window)
        self.pbar = QProgressBar(self)

        # When set to True, the progress bar shows the percentage of progress (only possible when running RAxML)
        percentage_bar = True
        if clustering_algorithm =="ML_raxml" or clustering_algorithm =="ML_iqtree":
            percentage_bar = False
            self.setWindowTitle('Searching Best ML tree')
        else: self.setWindowTitle('Generating Random Trees')

        # Setting maximum value of progress bar based on the number of random trees.
        if sys.platform != "win32" and bootstrap_program == "raxml" and  percentage_bar:
            self.pbar.setMaximum(num_random_trees)
        else:
            # Simple progressbar.
            self.pbar.setRange(0,0)

        self.layout = QVBoxLayout()
        self.layout.addWidget(self.pbar)
        self.setLayout(self.layout)

        # Centering QtWidgets ProgressBar.
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        self.resize(350, 100)

        # Stop Button
        self.button_stop = QPushButton('Stop', self)
        self.button_stop.clicked.connect(self.kill_thread)
        self.button_stop.move(240, 64)

        self.pymod = pymod
        self.pymod.show_tree = show_tree
        # NOTE: despite its name, 'pymod' here is the 'Tree_building' protocol
        # object, not the main PyMod object (the callers pass 'pymod=self').
        # A protocol only owns a 'pymod' attribute, so the assignment below
        # *creates* 'main_window' on it: it is read back a few lines down, when
        # building the Worker. Do not remove it.
        self.pymod.main_window = main_window
        self.pymod.alignments_files_names = alignments_files_names
        self.pymod.input_cluster_element.unique_index = input_cluster_element
        self.new_aln_file_path_phy = alignment_file_path_bootstrap
        self.alignment_file_path = alignment_file_path_best_tree
        self.new_aln_file_path = alignment_best_tree
        self.pymod.alignments_dirpath = alignments_dirpath
        self.tree_file_path = tree_file_path
        self.new_tree_file_name_output_file = tree_id
        self.tree = tree
        self.number_random_trees = num_random_trees
        self.clustering_algorithm = clustering_algorithm
        self.bootstrap_program = bootstrap_program
        self.raxml_PTHREADS_AVX_exe_filepath = raxml_exec
        self.iqtree2_exe_file_path = iqtree_exec
        self.tree_building_window = tree_building_window
        self.visualization = visualization
        self.new_aln_file_path = new_aln_file_path
        self.to_be_unrooted = to_be_unrooted

        # Class Thread and ProgressBar.
        self.thread = Worker(self.pymod, self.pymod.show_tree, self.pymod.main_window, self.pymod.alignments_files_names, self.pymod.input_cluster_element.unique_index,
                                    self.pymod.alignments_dirpath, self.new_aln_file_path_phy,self.alignment_file_path, self.new_aln_file_path, self.tree_file_path, self.new_tree_file_name_output_file,
                                     self.tree, self.number_random_trees, self.clustering_algorithm, self.bootstrap_program, self.raxml_PTHREADS_AVX_exe_filepath,
                                     self.iqtree2_exe_file_path, self.tree_building_window, self.visualization,self.new_aln_file_path, self.to_be_unrooted)

        self.thread.intReady.connect(self.on_count_changed)
        self.thread.finished.connect(self.thread.quit)
        self.thread.showTree.connect(self.show_tree)
        self.thread.closeProgressBar.connect(self.show_message)

    def start_progress(self):
        # Start Thread in the Worker class.
        self.show()
        self.thread.start()

    def kill_thread(self):
        '''
        Stop Thread and Subprocess via the 'STOP' button
        '''
        # Get child subprocess.
        self.subprocess = self.thread.get_subp()
        # Close ProgressBar.
        self.close()

        # The program which is building the tree is what has to be stopped: the
        # thread is only waiting for it, and ends as soon as it is gone. The
        # thread used to be killed first, with QThread.terminate(), which stops
        # a thread wherever it happens to be -- inside the interpreter, holding
        # the lock every other Python thread needs, PyMOL's own included. PyMOL
        # froze for good that way. Nothing is killed but the child now.
        try:
            if self.subprocess is not None:
                self.subprocess.kill()
        except Exception:
            pass
        try:
            self.thread.requestInterruption()
        except AttributeError:  # Qt versions without it
            pass
        self.thread.quit()
        self.thread.wait(500)
        # Remove incomplete files.
        try:
            # Remove temporary bootstrap files in pymod.alignments_dirpath.
            keywords = ["treefile","log","iqtree",".ckp.gz","splits.nex", "uniqueseq","RAxML"]
            for file in os.listdir(self.pymod.alignments_dirpath):
                for keyword in keywords:
                    if keyword in file:
                        os.remove(os.path.join(self.pymod.alignments_dirpath, file))
        except AttributeError:
            # Remove temporary IQTREE2 ML tree files in pymod.alignments_dirpath.
            file_name = ["bionj","iqtree","log","mldist","uniquese","boot","model.gz",".ckp.gz"]
            for file in os.listdir(self.pymod.alignments_dirpath):
                for keyword in file_name:
                    if keyword in file:
                        os.remove(os.path.join(self.pymod.alignments_dirpath, file))
        except :
            # Remove temporary RAXML ML tree search files in pymod.alignments_dirpath.
            for file in os.listdir(self.pymod.alignments_dirpath):
                    if file.startswith("RAxML"):
                        os.remove(os.path.join(self.pymod.alignments_dirpath, file))

    def on_count_changed(self, value):
        # Update ProgressBar.
        self.pbar.setValue(value)

    def show_tree(self):
        # Close ProgressBar
        self.close()
        if self.visualization == "bootstrap":
            if self.clustering_algorithm == "ML_raxml" or self.clustering_algorithm == "ML_iqtree":
                self.new_tree_file_path_nwk = self.thread.get_tree_file_path()
                self.tree = PhyloTree(self.new_tree_file_path_nwk, format = 1)

            # Ask which rendering (bubble or dendogram) should be used to represent bootstrap values.
            answer = askyesno_qt("Bootstrap Rendering", "Click 'Yes' if you want to display the bootstrap values as branches of a dendogram. By clicking 'No', the bootstrap values will be represented as spheres.")
            # Add bootstrap values/weights to internal nodes.
            i = 0
            self.bipartition_values = self.thread.get_bipartition_values()
            for n in self.tree.traverse():
                try:
                    if not n.is_leaf() and not n.is_root() and n.dist !=0 :
                        n.add_features(weight=int(self.bipartition_values[i]))
                        # Show branch weight.
                        support = faces.AttrFace("weight", fsize=9, ftype="Arial", fgcolor="darkred", formatter = "%0.3g")
                        n.add_face(support, i, position="branch-bottom")
                        i += 1
                except IndexError:
                    continue
            ts = TreeStyle()
            if not answer:
                ts.layout_fn = bubble_ly
                # Draw a tree in circular mode.
                ts.mode = "c"
            # add node names manually.
            ts.show_leaf_name = False
            # Show branch data.
            ts.show_branch_length = True
            ts.show_branch_support = False

            self.tree.show(tree_style=ts)

        elif self.clustering_algorithm == "ML_raxml" or self.clustering_algorithm == "ML_iqtree" and self.visualization != "bootstrap":
            self.new_tree_file_path_nwk = self.thread.get_tree_file_path()

            # Call the show_plot_window() method from Tree_building class.
            tree_ML = Show_tree_ML(self.pymod, self.pymod.show_tree, self.pymod.input_cluster_element, self.alignment_file_path, self.new_tree_file_path_nwk, self.new_aln_file_path,
                        self.tree_building_window,self.visualization, self.clustering_algorithm, self.to_be_unrooted)
            tree_ML.show_plot_window(self.new_tree_file_path_nwk)


    def show_message(self):
        # Alert message.
        title = "Bootstrap Analysis"
        message = "Bootstrapping analysis with %s and %s random trees has already been performed. You can visualize the final Tree via the 'Trees' menu" % (self.bootstrap_program, self.number_random_trees)
        self.pymod.main_window.show_error_message(title, message)
        # Close ProgressBar.
        self.close()


class Worker(QThread,QObject):
    finished = pyqtSignal()
    intReady = pyqtSignal(int)
    showTree = pyqtSignal()
    closeProgressBar = pyqtSignal()

    def __init__(self, pymod, show_tree, main_window, alignments_files_names, input_cluster_element, alignments_dirpath, alignment_file_path_bootstrap, alignment_file_path_best_tree,
                alignment_best_tree, tree_file_path, tree_id, tree,num_random_trees, clustering_algorithm, bootstrap_program, raxml_exec, iqtree_exec,
                tree_building_window, visualization, new_aln_file_path, to_be_unrooted):
        super().__init__()

        self.pymod = pymod
        self.pymod.show_tree = show_tree
        # As in 'PopUpProgressB', 'pymod' is the 'Tree_building' protocol object:
        # this assignment creates the attribute rather than overwriting PyMod's.
        self.pymod.main_window = main_window
        self.pymod.alignments_files_names = alignments_files_names
        self.pymod.input_cluster_element.unique_index = input_cluster_element
        self.pymod.alignments_dirpath = alignments_dirpath
        self.new_aln_file_path_phy = alignment_file_path_bootstrap
        self.alignment_file_path = alignment_file_path_best_tree
        self.new_aln_file_path = alignment_best_tree
        self.tree_file_path = tree_file_path
        self.new_tree_file_name_output_file = tree_id
        self.tree = tree
        self.number_random_trees = num_random_trees
        self.clustering_algorithm = clustering_algorithm
        self.bootstrap_program = bootstrap_program
        self.raxml_PTHREADS_AVX_exe_filepath = raxml_exec
        self.iqtree2_exe_file_path = iqtree_exec
        self.tree_building_window = tree_building_window
        self.visualization = visualization
        self.new_aln_file_path = new_aln_file_path
        self.to_be_unrooted = to_be_unrooted


    @pyqtSlot()
    def run (self):
        if self.bootstrap_program == "raxml" and self.clustering_algorithm!="ML_raxml" and self.clustering_algorithm!="ML_iqtree":
            '''
            Execute RAxML

            '''
            self.execute_raxml()


        elif self.bootstrap_program == "iqtree2" and self.clustering_algorithm!="ML_raxml" and self.clustering_algorithm!="ML_iqtree":
            '''
            Execute IQ-TREE 2
            '''
            self.execute_iqtree()

        else:
            if self.clustering_algorithm == "ML_raxml":
                tree_file_already_exist = False
                self.new_tree_file_path_nwk = os.path.join(self.pymod.alignments_dirpath, "%s_%s_tree_output_ML_RAxML.nwk" % (self.pymod.alignments_files_names, self.pymod.input_cluster_element.unique_index))
                # Check if RaXML maximum likelihood tree already exists (no need to run the ML-tree search again).
                if os.path.isfile(self.new_tree_file_path_nwk):
                    self.tree_file_path = self.new_tree_file_path_nwk
                    tree_file_already_exist = True

                else:
                    # RaXML to find best ML tree.
                    self.raxml_bestML_tree = "%s -m %s -s %s -p 54439 -w %s -# 20 -n ML20 -T auto" % (self.raxml_PTHREADS_AVX_exe_filepath, raxml_model_of(self.pymod.input_cluster_element), self.new_aln_file_path,self.pymod.alignments_dirpath)
                    self.execute_subprocess(str(self.raxml_bestML_tree), percentage = False)

                    # Rename RAxML output file name.
                    os.rename(os.path.join(self.pymod.alignments_dirpath,"RAxML_bestTree.ML20"), self.new_tree_file_path_nwk)
                    # Remove temporary files.
                    for file in os.listdir(self.pymod.alignments_dirpath):
                        if file.startswith("RAxML"):
                            os.remove(os.path.join(self.pymod.alignments_dirpath, file))

                    self.tree_file_path = self.new_tree_file_path_nwk

                    if self.bootstrap_program != "raxml" and self.bootstrap_program != "iqtree2" and tree_file_already_exist == False:
                        # 'Show Tree' signal.
                        self.showTree.emit()
                        # 'Stop Thread' signal.
                        self.finished.emit()


            elif self.clustering_algorithm == "ML_iqtree":
                tree_file_already_exist = False
                self.new_tree_file_path_nwk = os.path.join(self.pymod.alignments_dirpath, "%s_%s_tree_output_ML_IQTREE2.nwk" % (self.pymod.alignments_files_names, self.pymod.input_cluster_element.unique_index))
                # Check if IQTREE-2 maximum likelihood tree already exists (no need to run the ML-tree search again).
                if os.path.isfile(self.new_tree_file_path_nwk):
                    self.tree_file_path = self.new_tree_file_path_nwk
                    tree_file_already_exist = True

                else:
                    self.iqtree2_cline = "%s -s %s " % (self.iqtree2_exe_file_path, self.new_aln_file_path)
                    self.execute_subprocess(str(self.iqtree2_cline))

                    # Rename IQTREE2 output file name.
                    os.rename(os.path.join(self.pymod.alignments_dirpath,"%s.treefile" % self.new_aln_file_path), self.new_tree_file_path_nwk)
                    # Remove temporary IQTREE2 ML tree search files.
                    file_name = ["bionj","iqtree","log","mldist","uniquese","boot","model.gz",".ckp.gz"]
                    for file in os.listdir(self.pymod.alignments_dirpath):
                        for keyword in file_name:
                            if keyword in file:
                                os.remove(os.path.join(self.pymod.alignments_dirpath, file))

                    self.tree_file_path = self.new_tree_file_path_nwk

                    if self.bootstrap_program != "raxml" and self.bootstrap_program != "iqtree2" and tree_file_already_exist == False:
                        # 'Show Tree' signal.
                        self.showTree.emit()
                        # 'Stop Thread' signal.
                        self.finished.emit()

            if self.bootstrap_program == "raxml" and self.visualization == "bootstrap":
                '''
                Execute RAxML

                '''
                self.execute_raxml()



            elif self.bootstrap_program == "iqtree2" and self.visualization == "bootstrap":

                '''
                Execute IQ-TREE2

                '''
                tree = PhyloTree(self.tree_file_path, format = 1)

                # Convert FASTA to PHYLIP format.
                self.new_aln_file_path_phy = os.path.join(self.pymod.alignments_dirpath, "alignment_tmp_boot.phylip")
                with open(self.new_aln_file_path) as fasta_handle:
                    records = list(SeqIO.parse(fasta_handle, "fasta"))
                SeqIO.write(records, self.new_aln_file_path_phy, "phylip")

                self.execute_iqtree()


    def execute_raxml(self):
        # Check whether Bootstrap file already exists.
        random_trees_file_path = os.path.join(self.pymod.alignments_dirpath, "RAxML_bootstrap.%s_BOOT%s" %(self.new_tree_file_name_output_file,self.number_random_trees))
        if os.path.isfile(random_trees_file_path):
            self.closeProgressBar.emit()
        else:
            # Generates random trees from alignment.
            outputname = "%s_BOOT%s" %(self.new_tree_file_name_output_file,self.number_random_trees)
            self.raxml_cline_BOOT = "%s -f d -m %s -s %s -b 12463 -p 54439 -w %s -n %s -N %s" % (self.raxml_PTHREADS_AVX_exe_filepath, raxml_model_of(self.pymod.input_cluster_element), self.new_aln_file_path_phy, self.pymod.alignments_dirpath, outputname, self.number_random_trees)

            self.execute_subprocess(str(self.raxml_cline_BOOT))

            # Create file path for random trees.
            random_trees_file_path = os.path.join(self.pymod.alignments_dirpath, "RAxML_bootstrap.%s_BOOT%s" %(self.new_tree_file_name_output_file,self.number_random_trees))

            # Calculate bipartition values based on starting_tree.
            outputname =  "%s_consensus_BOOT%s" % (self.new_tree_file_name_output_file, self.number_random_trees)
            self.raxml_cline_BIPARTIONS = "%s -t %s -m %s -f b -w %s -n %s -z %s" % (self.raxml_PTHREADS_AVX_exe_filepath, self.tree_file_path, raxml_model_of(self.pymod.input_cluster_element), self.pymod.alignments_dirpath, outputname, random_trees_file_path)

            self.execute_subprocess(str(self.raxml_cline_BIPARTIONS), progressbar = False)

            # Rename RAxML consensus tree output in order to be found when visualizing the tree from the 'Trees' menu.
            consensus_tree_file_path = os.path.join(self.pymod.alignments_dirpath,"%s_consensusRAxML_BOOT%s" % (self.new_tree_file_name_output_file, self.number_random_trees))

            os.rename(os.path.join(self.pymod.alignments_dirpath,"RAxML_bipartitionsBranchLabels.%s_consensus_BOOT%s" % (self.new_tree_file_name_output_file,self.number_random_trees)),
                    consensus_tree_file_path)


            # Open consensus tree (generated by RAxML) and get confidence values (in square brackets).
            with open(consensus_tree_file_path,"r") as f:
                values = f.readlines()
                for v in values:
                    bipartition_pattern_values = re.findall(r"\[([0-9_]+)\]", v)
            self.bipartition_values = bipartition_pattern_values

            self.get_bipartition_values()

            # 'Show Tree' signal.
            self.showTree.emit()

            # 'Stop Thread' signal.
            self.finished.emit()

    def execute_iqtree(self):
        # Check whether Bootstrap file already exists.
        random_trees_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_consensusIQTREE2_BOOT%s" %(self.new_tree_file_name_output_file,self.number_random_trees))
        if os.path.isfile(random_trees_file_path):
            self.closeProgressBar.emit()
        else:
            
            ############MODIFICATION on 30/07/2025##########
            ## Explicitly specify the sequence type as amino acids with '-st AA'.
            # This is required when using an amino acid substitution model like JTT.
            # Without this flag, IQ-TREE defaults to DNA mode and will incorrectly interpret 'JTT' as a file name,
            # leading to the error: "ERROR: File not found JTT".
             
            """
            self.UFBOOT_cline = "%s -s %s -st AA -m JTT -t %s -B %s" % (
                self.iqtree2_exe_file_path,
                self.new_aln_file_path_phy,
                self.tree_file_path,
                self.number_random_trees
            )
            """
            ############END############
            
            #MODIFICATION SARA 29/03/2026
            #to handle both protein and NA

            self.UFBOOT_cline = "%s -s %s -m MFP -t %s -B %s" % (self.iqtree2_exe_file_path, self.new_aln_file_path_phy, self.tree_file_path, self.number_random_trees)

           #######END#####

            #OLD CODE not working
            #self.UFBOOT_cline = "%s -s %s -m JTT -t %s -B %s" % (self.iqtree2_exe_file_path, self.new_aln_file_path_phy, self.tree_file_path, self.number_random_trees)
            self.execute_subprocess(str(self.UFBOOT_cline))

            # Rename consensus tree file path.
            consensus_tree_file_path = os.path.join(self.pymod.alignments_dirpath,"%s_consensusIQTREE2_BOOT%s" % (self.new_tree_file_name_output_file, self.number_random_trees))
            os.rename(os.path.join(self.pymod.alignments_dirpath,"alignment_tmp_boot.phylip.contree"), consensus_tree_file_path)

            # Rename IQTREE2 info file path for this run.
            iqtree2_info_file_path = os.path.join(self.pymod.alignments_dirpath,"%s_IQTREE2_info_BOOT%s" % (self.new_tree_file_name_output_file, self.number_random_trees))
            os.rename(os.path.join(self.pymod.alignments_dirpath,"alignment_tmp_boot.phylip.iqtree"), iqtree2_info_file_path)

            # Remove temporary files.
            keywords = ["treefile","log",".ckp.gz","splits.nex", "uniqueseq"]
            for file in os.listdir(self.pymod.alignments_dirpath):
                for keyword in keywords:
                    if keyword in file:
                        os.remove(os.path.join(self.pymod.alignments_dirpath, file))

            # Open consensus tree (generated by IQTREE2) and get confidence values.
            with open(consensus_tree_file_path,"r") as f:
                values = f.readlines()
                for v in values:
                    self.bipartition_values = re.findall(r"\)([0-9_]+)\:",v)

            self.get_bipartition_values()

            # 'Show Tree' signal.
            self.showTree.emit()

            # 'Stop Thread' signal.
            self.finished.emit()

    def get_bipartition_values(self):
        return(self.bipartition_values)

    def get_tree_file_path (self):
        return(self.new_tree_file_path_nwk)

    def execute_subprocess(self,  commandline,
                           new_stdout=subprocess.PIPE, new_stderr=subprocess.PIPE,
                           new_shell=sys.platform,
                           verbose=True, progressbar = True, percentage = True):

        # "exec" causes cmd to inherit the shell process, instead of having the shell launch a child process.
        if sys.platform!="win32":
            self.subp = subprocess.Popen("exec " + commandline, stdout=new_stdout, stderr=new_stderr, shell=new_shell)
            if progressbar:
                # Update ProgressBar
                x = 1
                self.intReady.emit(x)
                while True:
                    # Print stdout process line by line.
                    lines = self.subp.stdout.readline()
                    print(lines)
                    if not lines:
                        break
                    if self.bootstrap_program == "raxml" and percentage:
                        # Progressbar signal outputs only when stdout writes bootstrap trees.
                        if len(lines) > 100 and len(lines) < 160:
                            x += 1
                            i = int(x)
                            time.sleep(0.5)
                            self.intReady.emit(i)
                        else: continue
                    else:
                        while True:
                            # Print stdout process line by line.
                            lines = self.subp.stdout.readline()
                            print(lines)
                            if not lines:
                                break
                            self.intReady.emit(0)
        else:
            self.subp = subprocess.Popen(commandline, stdout=new_stdout, stderr=new_stderr, shell=(sys.platform=="win32"))
            # Simple progessBar.
            if progressbar:
                while True:
                    self.intReady.emit(0)

        out_std, err_std = self.subp.communicate()
        returncode = self.subp.returncode

        if verbose:
            print("- Stdout:", out_std)

        if returncode != 0:
            if verbose:
                print("- Code:", returncode, ", Stderr:", err_std)

            def _readable(stream):
                if not stream:
                    return ""
                if isinstance(stream, bytes):
                    stream = stream.decode("utf-8", "replace")
                # The last lines are the ones which say what went wrong.
                lines = [line for line in stream.strip().split("\n") if line.strip()]
                return "\n".join(lines[-6:])

            detail = _readable(err_std) or _readable(out_std) or "no output"
            raise Exception("Subprocess returned non-zero return code: %s\n\n%s\n\n"
                            "The command was:\n%s" % (returncode, detail, commandline))
        else:
            os.system(commandline)


    def get_subp(self):
        '''
        Return subprocess
        '''
        return self.subp


class Show_tree_ML(Tree_building):
    def __init__(self,pymod, show_tree, input_cluster_element, alignment_file_path_best_tree, new_tree_file_path_nwk, new_aln_file_path, tree_building_window,visualization, clustering_algorithm, to_be_unrooted):
        super().__init__(self,pymod, input_cluster_element)

        self.pymod = pymod
        self.pymod.show_tree = show_tree
        self.pymod.input_cluster_element = input_cluster_element
        self.alignment_file_path = alignment_file_path_best_tree
        self.new_tree_file_path_nwk = new_tree_file_path_nwk
        self.new_aln_file_path = new_aln_file_path
        self.tree_building_window = tree_building_window
        self.visualization = visualization
        self.clustering_algorithm = clustering_algorithm
        self.to_be_unrooted = to_be_unrooted


class PyMod_form_item(QtWidgets.QWidget):
    pass

class Tree_building_window_qt(PyMod_protocol_window_qt, PyMod_form_item):

    def build_protocol_middle_frame(self):

        # Algorithm Tree-Building.
        self.algorithm_rds = PyMod_radioselect_qt(label_text="Clustering Algorithm",
                                                  buttons=list(sorted(tree_building_alg_dict.keys())))

        self.algorithm_rds.setvalue("Neighbor Joining")
        self.middle_formlayout.add_widget_to_align(self.algorithm_rds)

        self.calculator_rds = PyMod_radioselect_qt(label_text="Matrix Calculator",
                                                            buttons=("Identity","Blosum62"))
        self.calculator_rds.setvalue("Identity")

        self.algorithm_rds.get_button_at(0).clicked.connect(lambda a=None: PyMod_radioselect_qt.hide_widgets(self.calculator_rds))
        self.algorithm_rds.get_button_at(1).clicked.connect(lambda a=None: PyMod_radioselect_qt.hide_widgets(self.calculator_rds))
        self.algorithm_rds.get_button_at(2).clicked.connect(lambda a=None: PyMod_radioselect_qt.show_widgets(self.calculator_rds))
        self.algorithm_rds.get_button_at(3).clicked.connect(lambda a=None: PyMod_radioselect_qt.hide_widgets(self.calculator_rds))
        self.algorithm_rds.get_button_at(4).clicked.connect(lambda a=None: PyMod_radioselect_qt.show_widgets(self.calculator_rds))

        self.middle_formlayout.add_widget_to_align(self.calculator_rds)

        # Advanced options for visualization.
        self.advance_options_button = QtWidgets.QPushButton("Show Advanced Options")
        self.advance_options_button.clicked.connect(self.toggle_advanced_options_)
        self._advanced_options_label = QtWidgets.QLabel("")
        self.middle_formlayout.addRow(self.advance_options_button, self._advanced_options_label)

        # Rendering options.
        self.visualization_rds = PyMod_radioselect_qt(label_text="Visualization mode",
                                                  buttons=list(sorted(visualization_dict.keys())))
        self.visualization_rds.setvalue("Basic")

        self.middle_formlayout.add_widget_to_align(self.visualization_rds, advanced_option=True)
        self.middle_formlayout.set_input_widgets_width("auto", padding=10)


        # Bootstrap widgets are not been created yet.
        self.bootstrap_widget_enabled = False
        self.random_trees_rds_raxml_enabled = False
        self.random_trees_rds_iqtree_enabled = False


        # Show "Bootstrap Program " if "Bootstrap Analysis" is clicked.
        self.visualization_rds.get_button_at(2).clicked.connect(lambda a=None: self.add_widget(label_text = "Bootstrap program",
                                                                                                buttons=list(boot_program_dict.keys()), initial_value = None))
        self.visualization_rds.get_button_at(2).clicked.connect(lambda a=None: self.get_widget("bootstrap_program"))

        # 'Bootstrap program' visible only if Bootstrap Analysis is clicked.
        self.visualization_rds.get_button_at(0).clicked.connect(lambda a=None: self.on_click_hide_bootstrap_widget(0))
        self.visualization_rds.get_button_at(1).clicked.connect(lambda a=None: self.on_click_hide_bootstrap_widget(1))
        self.visualization_rds.get_button_at(3).clicked.connect(lambda a=None: self.on_click_hide_bootstrap_widget(3))
        self.visualization_rds.get_button_at(4).clicked.connect(lambda a=None: self.on_click_hide_bootstrap_widget(4))
        self.visualization_rds.get_button_at(5).clicked.connect(lambda a=None: self.on_click_hide_bootstrap_widget(5))
        self.visualization_rds.get_button_at(6).clicked.connect(lambda a=None: self.on_click_hide_bootstrap_widget(6))


    # Method for showing/hiding widgets for advanced options.
    def toggle_advanced_options_(self):
        if not self.showing_advanced_widgets:
            self.showing_advanced_widgets = True
            self.advance_options_button.setText("Hide Advanced Options")
            row_elements = []
            for row in self.middle_formlayout.widgets_to_align:
                if row.is_advanced_option:
                    row_elements.append(row)
                    row_elements[0].show_widgets()
        else:
            self.showing_advanced_widgets = False
            self.advance_options_button.setText("Show Advanced Options")
            for row in self.middle_formlayout.widgets_to_align:
                if row.is_advanced_option:
                    row.hide_widgets()

    ##########################
   # Dynamic widget creation #
   # for bootstrap analysis. #
   ###########################

    # Add Bootstrap buttons.
    def add_widget(self, label_text, buttons, initial_value):
        self.widget = PyMod_radioselect_qt(label_text=label_text, buttons=buttons)
        self.middle_formlayout.add_widget_to_align(self.widget, advanced_option=True)
        if not initial_value == None:
            self.widget.setvalue(initial_value)
        self.middle_formlayout.set_input_widgets_width("auto", padding=10)

    def get_widget(self, widget_id):
        if widget_id =="bootstrap_program":
            self.boot_program_rds =  self.widget
            self.bootstrap_widget_enabled = True
            self.show_child_widget_bootstrap()
        if widget_id =="random_trees_raxml":
            self.random_trees_rds_raxml = self.widget
            self.random_trees_rds_raxml_enabled = True
        if widget_id == "random_trees_iqtree":
            self.random_trees_rds_iqtree = self.widget
            self.random_trees_rds_iqtree_enabled = True

    # Show "number of random trees" buttons.
    def show_child_widget_bootstrap(self):
        if self.visualization_rds.get_button_at(2).isChecked():
            # Show IQ-TREE2 Random Trees button  if "UFBOOT" (IQ-TREE2) is clicked.
            self.boot_program_rds.get_button_at(0).clicked.connect(lambda a=None: self.add_widget(label_text="Number of random trees",
                                                                    buttons=list(random_trees_dict_IQTREE.keys()),initial_value="1000"))
            self.boot_program_rds.get_button_at(0).clicked.connect(lambda a=None: self.get_widget("random_trees_iqtree"))
            # Hide Random Trees RAxML.
            self.boot_program_rds.get_button_at(0).clicked.connect(lambda a=None: self.on_click_hide_random_trees_widget(0))

            # Show RaXML Random Trees button  if "RaXML" is clicked.
            self.boot_program_rds.get_button_at(1).clicked.connect(lambda a=None: self.add_widget(label_text="Number of random trees",
                                                                    buttons=list(random_trees_dict_RAxML.keys()),initial_value="100"))
            self.boot_program_rds.get_button_at(1).clicked.connect(lambda a=None: self.get_widget("random_trees_raxml"))
            # Hide Random Trees IQ-TREE2
            self.boot_program_rds.get_button_at(1).clicked.connect(lambda a=None: self.on_click_hide_random_trees_widget(1))

    def on_click_hide_random_trees_widget(self, button):
        # Hide Random Trees RAxML when IQ-TREE2 is selected.
        if button == 0:
            if self.random_trees_rds_raxml_enabled:
                PyMod_radioselect_qt.hide_widgets(self.random_trees_rds_raxml)
        # Hide Random Trees IQ-TREE2 when RAxML is selected.
        elif button == 1:
            if self.random_trees_rds_iqtree_enabled:
                PyMod_radioselect_qt.hide_widgets(self.random_trees_rds_iqtree)


    def on_click_hide_bootstrap_widget(self,i):
        # Hide Bootstrap widgets.
        if self.visualization_rds.get_button_at(i).isChecked():
            if self.bootstrap_widget_enabled:
                PyMod_radioselect_qt.hide_widgets(self.boot_program_rds)
            if self.random_trees_rds_raxml_enabled:
                PyMod_radioselect_qt.hide_widgets(self.random_trees_rds_raxml)
            if self.random_trees_rds_iqtree_enabled:
                PyMod_radioselect_qt.hide_widgets(self.random_trees_rds_iqtree)



class Tree_visualization(Evolutionary_analysis_protocol):

    def launch_from_gui(self):
        """
        It will check if a software to build a tree is available on the user's machine.
        """
        self.tree_visualization_software = None
        can_show_tree = False
        if self.check_tree_visualization_module():
            self.tree_visualization_software = "ete3"
            can_show_tree = True

        if can_show_tree:
            self.tree_visualization_window()
        else:
            title = "Tree Visualization Error"
            message = "To render the tree via the ETE-toolkit you need to install the package."
            self.pymod.main_window.show_error_message(title, message)

    def check_tree_visualization_module(self):
        try:
            import ete3
            return True
        except ImportError:
            return False

    def tree_visualization_window(self):

        """
        Builds a window with options to show a tree out of an alignment.
        """

        self.tree_visualization_window = Tree_visualization_qt(self.pymod.main_window,
            protocol=self,
            title="Options for Tree Visualization",
            upper_frame_title="Here you can modify options for Tree Visualization",
            submit_command=self.show_plot_window_from_list)
        self.tree_visualization_window.show()

    # The alignment a tree read from a file is drawn against. Three of the views need
    # one: the alignment beside the tree, and the two mutation events.
    #
    # They used to look for 'alignments/alignment_tmp.phylip' and ask for a file only
    # when it was not there. That path is the scratch file of the last tree built in
    # this session and has nothing to do with a tree read from a file, so whenever a
    # tree had been built earlier the imported one was drawn against an alignment of
    # other sequences, with nothing said. The file is in PHYLIP too, where ETE reads
    # what it is handed as FASTA.

    alignment_file_filter = "*.fasta *.fa *.fas *.aln *.phy *.phylip *.sth"

    def read_alignment_as_fasta(self, alignment_file_path):
        """
        The alignment in this file as FASTA text, which is what 'link_to_alignment'
        reads. The format is worked out by trying: a file chosen by hand can be in any
        of the ones PyMod itself writes. Returns None when it is not an alignment.
        """

        for file_format in ("fasta", "clustal", "phylip", "phylip-relaxed",
                            "stockholm", "nexus"):
            try:
                with open(alignment_file_path) as alignment_handle:
                    records = list(AlignIO.read(alignment_handle, file_format))
            except Exception:
                continue
            if records:
                return "".join(">%s\n%s\n" % (record.id, str(record.seq))
                               for record in records)
        return None

    def ask_for_the_alignment(self, tree):
        """
        Ask for the alignment this tree was built from, read it, and check that it
        belongs to the tree. Returns FASTA text, or None when there is nothing to draw
        with -- the user declined, or the file was not an alignment of these sequences.
        """

        if not askyesno_qt("Alignment Needed",
                           "This tree was read from a file, so PyMod does not know"
                           " which alignment it was built from.\n\nOpen the alignment"
                           " file?",
                           parent=self.pymod.get_qt_parent()):
            return None

        chosen_file = askopenfile_qt("Open the alignment of this tree",
                                     name_filter=self.alignment_file_filter,
                                     parent=self.pymod.get_qt_parent())
        if not chosen_file:
            return None
        alignment_file_path = getattr(chosen_file, "name", chosen_file)

        alignment = self.read_alignment_as_fasta(alignment_file_path)
        if alignment is None:
            self.pymod.main_window.show_error_message(
                "Alignment Error",
                "'%s' could not be read as an alignment."
                % os.path.basename(alignment_file_path))
            return None

        headers = set(line[1:].split()[0] for line in alignment.splitlines()
                      if line.startswith(">"))
        leaves = [leaf.name for leaf in tree.iter_leaves()]
        missing = [name for name in leaves if name not in headers]

        if len(missing) == len(leaves):
            self.pymod.main_window.show_error_message(
                "Alignment Error",
                "None of the %s leaves of this tree is in '%s', so it is not the"
                " alignment this tree was built from.\n\nThe tree calls its first"
                " leaves: %s."
                % (len(leaves), os.path.basename(alignment_file_path),
                   ", ".join(leaves[:3])))
            return None

        if missing:
            self.pymod.main_window.show_warning_message(
                "Alignment Warning",
                "%s of the %s leaves of this tree are not in that alignment, and will"
                " be drawn without a sequence: %s."
                % (len(missing), len(leaves), ", ".join(missing[:5])))

        return alignment


    def show_plot_window_from_list(self):

        visualization = self.get_visualization_mode()

        if visualization == "bp":
            self.pymod.show_tree(self.input_cluster_element)

        elif visualization == "aln":
            tree = PhyloTree(self.input_cluster_element, format=1)
            alignment = self.ask_for_the_alignment(tree)
            if alignment is None:
                return None
            tree.link_to_alignment(alignment)
            tree.show(tree_style=TreeStyle())

        elif visualization == "ancestors":

            ########OLD code with modification, erasing the try section###########

            tree = PhyloTree(self.input_cluster_element, format = 1)
            
            not_ancestor_list = []
            for n in tree.traverse():
                if n.is_leaf():
                    sister_node = n.get_sisters()
                    for s in sister_node:
                        try:
                            ancestor = tree.get_common_ancestor(n.name,s.name)
                            nst = NodeStyle()
                            # Generate random colors.
                            h,s,l = random.random(), 0.5 + random.random()/2.0, 0.4 + random.random()/5.0
                            r= [int(256*i) for i in colorsys.hls_to_rgb(h,l,s)]
                            color = '#%02x%02x%02x' % (r[0],r[1],r[2])
                            nst["bgcolor"] = color
                            ancestor.set_style(nst)
                        except:
                            not_ancestor_list.append(n.name)

            ts = TreeStyle()
            ts.layout_fn = layout
            ts.show_leaf_name = False
            ts.mode = "c"
            ts.root_opening_factor = 10
            tree.show(tree_style = ts)
            """OLD CODE connected with EVO 
            tree = PhyloTree(self.input_cluster_element, format = 1)
            try:
                not_ancestor_list = []
                for n in tree.traverse():
                    if n.is_leaf():
                        sister_node = n.get_sisters()
                        for s in sister_node:
                            try:
                                ancestor = tree.get_common_ancestor(n.name,s.name)
                                nst = NodeStyle()
                                # Generate random colors.
                                h,s,l = random.random(), 0.5 + random.random()/2.0, 0.4 + random.random()/5.0
                                r= [int(256*i) for i in colorsys.hls_to_rgb(h,l,s)]
                                color = '#%02x%02x%02x' % (r[0],r[1],r[2])
                                nst["bgcolor"] = color
                                ancestor.set_style(nst)
                            except:
                                not_ancestor_list.append(n.name)

                ts = TreeStyle()
                ts.layout_fn = layout
                ts.show_leaf_name = False
                ts.mode = "c"
                ts.root_opening_factor = 10
                tree.show(tree_style = ts)
            """

            """OLD CODE for EVO 
            except:
                #Here leave the old code simple return None or the new one?
               
                title = "Evolutionary Events"
                message = "The Tree is not rooted."
                
                # Show an initial error message indicating that the tree is not rooted
                self.pymod.main_window.show_error_message(title, message)
                
                #### OLD BEHAVIOR
                # Previously: silently passed or returned None, which could lead to misleading success messages
                # pass
                # return None

                #### MODIFICATION on 31/07/2025 ###########
                # Now we explicitly call a function to ask the user whether they want to root the tree automatically.
                # If the user agrees, the tree will be midpoint-rooted and the root will be retrieved.
                # If the user declines or if rooting fails, the function handles that internally.
                ###########################################
                
                #Tree_building.handle_non_rooted_tree_error()  # Call the handler to prompt user interaction and attempt recovery
                ##############END##############

                ##################MODIFICATION SARA ON 11/04/2026#######################
                tree = PhyloTree(self.input_cluster_element, format=1)

                # FIX 1: define aln
                aln = None

                # FIX 2: try to load alignment
                try:
                    with open(self.alignment_file_path, 'r') as f:
                        aln = f.read()
                except FileNotFoundError:
                    pass  # or your dialog logic

                # ROOTING
                tree = Tree_building.handle_non_rooted_tree_error(self, tree)

                if tree is None:
                    return

                # ALIGNMENT CHECK
                if aln is None:
                    return

                # ANALYSIS
                tree_params = ColumnMutation()
                t, ts = tree_params.get_mutation_event(tree, aln, True)

                t.show(tree_style=ts)


                ##################END################################
            """


        ###############àMODIFICATION SARA on 12/04/2026######################
        elif visualization == "evo":

                #####NEW#######12/04/2026


               # 1. LOAD TREE (NEWICK)
                tree = PhyloTree(self.input_cluster_element, format=1)

                # 2. ROOT TREE
                tree = Tree_building.handle_non_rooted_tree_error(self, tree)
                if tree is None:
                    return

                # 3. EVOLUTIONARY EVENTS (NO ALIGNMENT REQUIRED VERSION)
                try:
                    ts = TreeStyle()

                    events = tree.get_descendant_evol_events()

                    i = 0
                    for n in tree.traverse():
                        if not n.is_leaf():

                            if i >= len(events):
                                break

                            ev = events[i]

                            if ev.etype == "S":
                                n.add_face(TextFace("S", fgcolor="blue", penwidth=8),
                                        column=i, position="branch-top")

                            elif ev.etype == "D":
                                n.add_face(TextFace("D", fgcolor="red", penwidth=8),
                                        column=i, position="branch-top")

                            i += 1

                    ts.mode = "c"
                    ts.show_leaf_name = True

                    tree.show(tree_style=ts)

                except Exception as e:
                    self.pymod.main_window.show_error_message(
                        "Evolutionary Events Error",
                        str(e)
                    )
        ################END####################

        elif visualization == "mut_single":
            tree = PhyloTree(self.input_cluster_element, format=1)
            alignment = self.ask_for_the_alignment(tree)
            if alignment is None:
                return None
            shown_tree, tree_style = ColumnMutation().get_mutation_event(tree, alignment,
                                                                        True)
            shown_tree.show(tree_style=tree_style)

        elif visualization == "mut_all":
            tree = PhyloTree(self.input_cluster_element, format=1)
            alignment = self.ask_for_the_alignment(tree)
            if alignment is None:
                return None
            shown_tree, tree_style = ColumnMutation().get_mutation_event(tree, alignment,
                                                                        False)
            shown_tree.show(tree_style=tree_style)

        """OLD CODE REDUNDANT (NOT USE FOR NOW)
        
        if visualization == "evo":

            ##############NEW##############

             # 1. BUILD TREE
            tree = PhyloTree(self.input_cluster_element, format=1)

            # 2. ROOT TREE (IMPORTANT FIX)
            tree = Tree_building.handle_non_rooted_tree_error(self, tree)
            if tree is None:
                return

            # 3. LOAD ALIGNMENT
            try:
                with open(self.alignment_file_path, "r") as f:
                    aln = f.read()
            except Exception as e:
                self.pymod.main_window.show_error_message(
                    "Alignment Error",
                    str(e)
                )
                return

            # 4. EVOLUTIONARY EVENTS (SAFE)
            try:
                ts = TreeStyle()

                events = tree.get_descendant_evol_events()

                i = 0
                for n in tree.traverse():
                    if not n.is_leaf():

                        if i >= len(events):
                            break

                        ev = events[i]

                        if ev.etype == "S":
                            n.add_face(
                                TextFace("S", fgcolor="blue", penwidth=8),
                                column=i,
                                position="branch-top"
                            )

                        elif ev.etype == "D":
                            n.add_face(
                                TextFace("D", fgcolor="red", penwidth=8),
                                column=i,
                                position="branch-top"
                            )

                        i += 1

                ts.mode = "c"
                ts.show_leaf_name = True

                tree.show(tree_style=ts)

            except Exception as e:
                self.pymod.main_window.show_error_message(
                    "Evolutionary Events Error",
                    str(e)
                )
                return
            """
            
        """

        
            tree = PhyloTree(self.input_cluster_element, format = 1)
            try:
                ts = TreeStyle()
                # Get evolutionary events for each internal node.
                events = tree.get_descendant_evol_events()
                i = 0
                for n in tree.traverse():
                    if not n.is_leaf():
                        ev = events[i]
                        if ev.etype == "S":
                            # Add node face.
                            n.add_face(TextFace("S",fgcolor = "blue", penwidth = 8), column = i, position = "branch-top")
                        elif ev.etype == "D":
                            # Add node face.
                            n.add_face(TextFace("D", fgcolor = "red", penwidth = 8), column = i, position = "branch-top")
                        i += 1

                ts.mode = "c"
                ts.show_leaf_name = True
                tree.show(tree_style=ts)
            except TypeError:
                title = "Evolutionary Events"
                message = "The Tree is not rooted."
                
                # Show an initial error message indicating that the tree is not rooted
                self.pymod.main_window.show_error_message(title, message)
                
                #### OLD BEHAVIOR
                # Previously: silently passed or returned None, which could lead to misleading success messages
                # pass
                # return None

                #### MODIFICATION on 31/07/2025 ###########
                # Now we explicitly call a function to ask the user whether they want to root the tree automatically.
                # If the user agrees, the tree will be midpoint-rooted and the root will be retrieved.
                # If the user declines or if rooting fails, the function handles that internally.
                ###########################################
                
                #Tree_building.handle_non_rooted_tree_error()  # Call the handler to prompt user interaction and attempt recovery
                ##################END################


                ###################MODIFICATION SARA ON 11/04/2026##################
                tree = PhyloTree(self.input_cluster_element, format=1)

                # FIX 1: define aln
                aln = None

                # FIX 2: try to load alignment
                try:
                    with open(self.alignment_file_path, 'r') as f:
                        aln = f.read()
                except FileNotFoundError:
                    pass  # or your dialog logic

                # ROOTING
                tree = Tree_building.handle_non_rooted_tree_error(self, tree)

                if tree is None:
                    return

                # ALIGNMENT CHECK
                if aln is None:
                    return

                # ANALYSIS
                tree_params = ColumnMutation()
                t, ts = tree_params.get_mutation_event(tree, aln, True)

                t.show(tree_style=ts)


                #################END#############################
                """

                
        if visualization == "bootstrap":
            tree = PhyloTree(self.input_cluster_element, format = 1)

            # Find consensus_tree file/files.
            file_list = glob.glob(os.path.join(self.pymod.alignments_dirpath,"%s_consensus*" % self.input_cluster_element))

            # Get Alignment_id.
            aln_id = os.path.basename(self.input_cluster_element)

            if len(file_list) == 1:
                self.consensus_tree_file_path = file_list[0]
                self.show_bootstrap_tree(tree)
            # Bootstrap Analysis run with different number of random trees.
            elif len(file_list) > 1:
                self.consensus_tree_file_path = askopenfile_qt("Open an Bootstrap Tree file",
                                               initialdir = self.pymod.alignments_dirpath,
                                               name_filter=("%s_consensus*BOOT*" % aln_id),
                                               parent=self.pymod.get_qt_parent())
                self.show_bootstrap_tree(tree)

            elif len(file_list) == 0:
                title = "Bootstrap File Not Found"
                message = "To visualize the consensus tree, first run the Bootstrap Analysis via the 'Alignments'-->'Build Tree From Alignment'-->'Advanced Options' menu."
                self.pymod.main_window.show_error_message(title, message)
                pass

    def show_bootstrap_tree (self,tree):
        # UPGMA and Parsimony need to be unrooted.
        unrooted = True
        # Ask which rendering (bubble or dendogram) should be used to represent bootstrap values.
        answer = askyesno_qt("Bootstrap Rendering", "Click 'Yes' if you want to display the bootstrap values as branches of a dendogram. By clicking 'No', the bootstrap values will be represented as spheres.")
        while unrooted:
            try:
            # Open consensus tree (generated by RAxML) and get confidence values (in square brackets).
                with open(self.consensus_tree_file_path,"r") as f:
                    values = f.readlines()
                    if "RAxML" in str(self.consensus_tree_file_path):
                        for v in values:
                            # RAxML bipartition pattern.
                            bipartition_values = re.findall(r"\[([0-9_]+)\]", v)
                    else:
                        for v in values:
                            # IQ-TREE2 bipartition pattern.
                            bipartition_values = re.findall(r"\)([0-9_]+)\:",v)
                # Add bootstrap values/weights to internal nodes.
                i = 0
                for n in tree.traverse():
                    try:
                        if not n.is_leaf() and not n.is_root() and n.dist !=0 :
                            n.add_features(weight=int(bipartition_values[i]))
                            # Show branch weight.
                            support = faces.AttrFace("weight", fsize=9, ftype="Arial", fgcolor="darkred", formatter = "%0.3g")
                            n.add_face(support, i, position="branch-bottom")
                            i += 1
                    except IndexError:
                        continue
                unrooted = False
                ts = TreeStyle()
                if not answer:
                    ts.layout_fn = bubble_ly
                    ts.mode = "c"
                ts.show_leaf_name = False
                ts.show_branch_length = True
                ts.show_branch_support = False

                tree.show(tree_style=ts)

            except IndexError:
                tree.unroot() # converts a split into multifurcation
                # Save the unrooted tree in a new file.
                tree_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_unrooted_tree.nwk" % self.input_cluster_element)
                tree.write(format = 1, outfile = tree_file_path)
                tree = PhyloTree(tree_file_path, format = 1)
                continue


    def get_visualization_mode(self):
        return visualization_dict[self.tree_visualization_window.visualization_rds.getvalue()]


    def edit_FASTA_file(self, input_file, output_file):
        '''
        Change the record id in a FASTA file using BioPython.
        '''
        with open(input_file) as original, open(output_file, 'w') as corrected:
            records = SeqIO.parse(original, 'fasta')
            print(records)
            # Duplicate names found errors.
            bad_chars = ['_', 'c', 'h', 'a', 'i', 'n'] #"_chain_"
            for record in records:
                if len(record.id) > 10:
                    if "chain" in record.id:
                        for i in bad_chars :
                            record.id = record.id.replace(i, '')
                    record.id = record.id[0:10]
                    record.description = record.id[0:10]
                SeqIO.write(record, corrected, 'fasta')

class Tree_visualization_qt(PyMod_protocol_window_qt):

    def build_protocol_middle_frame(self):

        self.visualization_rds = PyMod_radioselect_qt(label_text="Visualization mode",
                                                  buttons=list(sorted(visualization_dict.keys())))
        self.visualization_rds.setvalue("Basic")
        self.middle_formlayout.add_widget_to_align(self.visualization_rds)

        self.middle_formlayout.set_input_widgets_width("auto", padding=10)


####################END ETE#####################################