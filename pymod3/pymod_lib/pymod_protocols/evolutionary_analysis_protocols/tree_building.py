# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

# TODO:
#   - add the possibility to save trees in the phylip format.

import os

from pymod_lib import pymod_vars

from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_protocol_window_qt, PyMod_radioselect_qt

from ._evolutionary_analysis_base import Evolutionary_analysis_protocol


tree_building_alg_dict = {"Neighbor Joining": "nj", "UPGMA": "upgma"}

class Tree_building(Evolutionary_analysis_protocol):

    def launch_from_gui(self):
        """
        It will check if a software to build a tree is available on the user's machine.
        """
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
