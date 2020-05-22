# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing a loop modeling protocol using MODELLER in PyMod.
"""

import os
import shutil

import pymod_lib.pymod_structure as pmstr

from .homology_modeling import MODELLER_homology_modeling, Modeling_cluster


###################################################################################################
# Classes for the loop modeling protocol.                                                         #
###################################################################################################

###############################################################################
# Main class.                                                                 #
###############################################################################

class MODELLER_loop_refinement(MODELLER_homology_modeling):
    """
    Class for performing loop refinement of a protein PDB structure.
    Implements the following protocol: https://salilab.org/modeller/manual/node36.html.
    """

    protocol_name = "modeller_loop_refinement"

    def modeling_protocol_initialization(self):
        self.modeling_mode = "loop_refinement"
        self.modeling_window_class_qt = Loop_refinement_window_qt


    def initialize_modeling_protocol_from_gui(self):

        #---------------------------------------------------------------------
        # Get the selected sequences to see if there is a correct selection. -
        #---------------------------------------------------------------------

        selected_sequences = self.pymod.get_selected_sequences()

        # First check if at least one sequence is selected.
        if not len(selected_sequences) > 0:
            self.pymod.main_window.show_error_message("Selection Error", "Please select at least one target sequence to use MODELLER.")
            return None

        # Checks if all the selected sequences can be used to build a model.
        if False in [s.has_structure() for s in selected_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only sequences that have a structure loaded in PyMOL.")
            return None

        # Checks that no nucleic acids have been selected.
        if True in [e.polymer_type == "nucleic_acid" for e in selected_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Can not perform loop refinement with structures with nucleic acids.")
            return None

        # Checks that all target chains come from the same PDB file.
        original_pdbs_set = set([t.get_structure_file(full_file=True) for t in selected_sequences])
        if len(original_pdbs_set) > 1:
            self.pymod.main_window.show_error_message("Selection Error", "In order to perform multiple chain loop modeling, please select only chains coming from the same original PDB file.")
            return None

        # This will contain a list of 'Modeling_cluster' objects.
        self.loop_modeling_targets = selected_sequences
        self.modeling_clusters_list = [Loop_modeling_cluster(self, t) for t in self.loop_modeling_targets]
        # Fixes the 'block_index' attribute so that each block index will correspond to the
        # numeric ID of the chain in the 3D model.
        for mc_idx, mc in enumerate(self.get_modeling_clusters_list(sorted_by_id=True)):
            mc.block_index = mc_idx

        # Define the modeling mode.
        self.multiple_chain_mode = len(self.loop_modeling_targets) > 1


        # Single chain homology modeling mode.
        if not self.multiple_chain_mode:
            self.build_modeling_window()

        # Multiple chains modeling requires additional controls.
        else:
            # This will be used to build the symmetry restraints options for loop modeling.
            self.initialize_multichain_modeling()
            self.build_modeling_window()


    def _get_multichain_options_from_gui(self):
        self.template_complex_name = "starting_complex.pdb"
        self.template_complex_modeller_name = self.template_complex_name[:-4]


    def check_modeling_cluster_parameters(self, modeling_cluster):
        """
        Checks the if there are any problems with the user-supplied parameters of a "modeling cluster"
        before starting the modeling process.
        """
        # Checks if there are some templates that have been selected.
        if modeling_cluster.loops_list == []:
            raise ValueError("You have to select at least one loop for target '%s' in order to perform refinement." % (modeling_cluster.target_name))


    def prepare_modeling_protocol_session_files(self):

        # Copies the target chains files in the modeling directory.
        for mc_i, mc in enumerate(self.modeling_clusters_list):
            """
            # Copy the templates structure files in the modeling directory.
            tar_file = mc.target.get_structure_file(basename_only=False)
            copied_tar_file = os.path.basename(tar_file)
            shutil.copy(tar_file, os.path.join(self.modeling_protocol.modeling_dirpath, copied_tar_file))
            mc.structure_file = copied_tar_file
            """

            tar_file = mc.target.get_structure_file(basename_only=False)
            copied_tar_file = os.path.basename(tar_file)
            shutil.copy(tar_file, os.path.join(self.modeling_dirpath, copied_tar_file))
            if not self.multiple_chain_mode and mc_i == 0:
                self.loop_refinement_starting_model = copied_tar_file


        # Prepares the starting complex file for multiple chain loop refinement.
        if self.multiple_chain_mode:
            list_of_starting_chains_files = []
            for t in self.get_targets_list(sorted_by_id=True):
                list_of_starting_chains_files.append(os.path.join(self.modeling_dirpath, t.get_structure_file()))
            pmstr.join_pdb_files(list_of_starting_chains_files, os.path.join(self.modeling_dirpath, self.template_complex_name))
            self.loop_refinement_starting_model = self.template_complex_name


###############################################################################
# Modeling clusters.                                                          #
###############################################################################

class Loop_modeling_cluster(Modeling_cluster):
    """
    Class to represent a single target chain for loop modeling. It is a simpler version of the
    'Modeling_cluster' class used for the homology modeling protocol.
    """

    def __init__(self, modeling_protocol, target):
        self.modeling_protocol = modeling_protocol
        self.target = target
        self.target_name = self.get_target_name()
        self.model_elements_list = []
        self.block_index = target.get_chain_numeric_id()
        self.model_chain_id = target.get_chain_id()
        self.restraints_options_dict = {}


    def set_options_from_gui(self):

        self.loops_list = []

        # Gets a list of elements like: ('GLN 2', 'LYS 6')
        # defining the starting and ending residues of  loop.
        _loops_list = self.user_loop_selector_frame.user_defined_loops

        for l in _loops_list:
            # For each loop, builds a tuple like: ("2", "6", "A")
            # where the third element if the chain of the target chain.
            res_1_pdb_id = l[0].split()[1]
            res_2_pdb_id = l[1].split()[1]
            self.loops_list.append((res_1_pdb_id, res_2_pdb_id, self.target.structure.chain_id))


    def use_water_in_cluster(self):
        return False

    def get_symmetry_chain_id(self):
        return self.target.get_chain_id()


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

from pymol.Qt import QtWidgets, QtCore

from pymod_lib.pymod_gui.shared_gui_components_qt import small_font_style
from .homology_modeling._gui_qt import (Modeling_window_qt, Optimization_level_frame_qt, Custom_obj_func_frame_qt,
                                        modeling_window_title_style, modeling_options_sections_style,
                                        modeling_options_subsections_style)


################################################################################
# Main window.                                                                 #
################################################################################

class Loop_refinement_window_qt(Modeling_window_qt):
    """
    Loop modeling window.
    """

    optimization_level_choices = ("Approximate", "Low", "Mid", "Default", "High", "Custom")

    def build_modeling_protocol_main_page(self):
        """
        Starts to insert content in the "Main" page.
        """

        # Builds a frame for each modeling_cluster.
        for (i, modeling_cluster) in enumerate(self.protocol.modeling_clusters_list):

            # Add a spacer to separate the sections for each modeling cluster.
            if i != 0 and self.protocol.multiple_chain_mode:
                spacer_frame = QtWidgets.QFrame()
                spacer_frame.setFrameShape(QtWidgets.QFrame.HLine)
                spacer_frame.setFrameShadow(QtWidgets.QFrame.Sunken)
                self.main_page_interior_layout.addRow(spacer_frame)

            show_symmetry_restraints_option = False

            # Widgets necessary to choose the templates for a single target sequence.
            modeling_option_label = QtWidgets.QLabel("Modeling options for target: %s" % (modeling_cluster.target_name))
            modeling_option_label.setStyleSheet(modeling_window_title_style)
            self.main_page_interior_layout.addRow(modeling_option_label)

            if self.protocol.multiple_chain_mode:

                # Use symmetry restraints option.
                if modeling_cluster.symmetry_restraints_id != None:

                    show_symmetry_restraints_option = True
                    symmetry_checkbox, symmetry_info = self.build_symmetry_restraints_option(modeling_cluster)


            if any((show_symmetry_restraints_option, )): # This might include other flags in future releases.
                additional_options_label = QtWidgets.QLabel("Restraints options")
                additional_options_label.setStyleSheet(modeling_options_sections_style)
                self.main_page_interior_layout.addRow(additional_options_label)
                if show_symmetry_restraints_option:
                    self.main_page_interior_layout.addRow(symmetry_checkbox, symmetry_info)


            # Build a series of widgets for loop ranges selection.
            modeling_option_label = QtWidgets.QLabel("Loop selection")
            modeling_option_label.setStyleSheet(modeling_options_sections_style)
            self.main_page_interior_layout.addRow(modeling_option_label)


            # Build a 'User_loop_selector_frame' object, which will contain information about
            # user-defined loops for this target.
            uls = User_loop_selector_frame_qt(parent=None, parent_window=self,
                                              modeling_cluster=modeling_cluster)
            modeling_cluster.user_loop_selector_frame = uls
            self.main_page_interior_layout.addRow(uls)


    def get_optimization_level_class(self):
        return Optimization_level_frame_loop_qt

    def get_custom_obj_func_frame_class(self):
        return Custom_obj_func_frame_loop_qt


################################################################################
# Selectors for loop ranges.                                                   #
################################################################################

class User_loop_selector_frame_qt(QtWidgets.QFrame):
    """
    Each modeling cluster will be used to build an object of this class. It will be used to let
    users define loop ranges in the model chains.
    """

    def __init__(self, parent, parent_window, modeling_cluster, *args, **configs):
        self.parent_window = parent_window
        self.protocol = self.parent_window.protocol
        self.modeling_cluster = modeling_cluster
        super(User_loop_selector_frame_qt, self).__init__(parent, *args, **configs)
        self.initUI()
        self.initialize_list()

        self.loop_selector_frame_layout.setAlignment(QtCore.Qt.AlignLeft)
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)


    def initUI(self):
        self.loop_selector_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.loop_selector_frame_layout)
        label_text = "Select two residues for a loop"
        self.modeling_cluster_loop_selection_label = QtWidgets.QLabel(label_text)
        self.modeling_cluster_loop_selection_label.setStyleSheet(modeling_options_subsections_style)
        self.loop_selector_frame_layout.addWidget(self.modeling_cluster_loop_selection_label, 0, 0, 1, 3)


    def initialize_list(self):
        """
        Build the initial row in the user-defined loops frame.
        """

        # The target chain sequence.
        self.target_seq = self.modeling_cluster.target.my_sequence

        # This is going to contain User_loop_combo objects.
        self.list_of_loop_combos = []

        # This list is going to contain info about loop ranges defined by the user through the
        # GUI.
        self.user_defined_loops = []

        # Creates the comboboxes.
        first = User_loop_combo_qt(self, self.modeling_cluster.target.get_polymer_residues())
        self.list_of_loop_combos.append(first)


    def add_new_user_loop(self):
        """
        This is called when the "Add" button to add a user-defined loop is pressed.
        """

        res_1 = self.list_of_loop_combos[-1].get_res1_combobox_val()
        res_2 = self.list_of_loop_combos[-1].get_res2_combobox_val()

        # Checks that both the comboboxes have been used to select a residue.
        if res_1 == "" or res_2 == "":
            txt = "You have to select two residues to define a loop."
            self.protocol.pymod.main_window.show_warning_message("Warning", txt)
            return None

        # Checks that the same residue has not been selected in both comboboxes.
        elif res_1 == res_2:
            txt = "You cannot select the same starting and ending residue to define a loop."
            self.protocol.pymod.main_window.show_warning_message("Warning", txt)
            return None

        res_1_pdb_id = int(res_1.split()[1])
        res_2_pdb_id = int(res_2.split()[1])

        if res_1_pdb_id >= res_2_pdb_id:
            txt = "The index of the starting residue of a loop (you selected %s) can not be higher than the ending residue index (you selected %s)." % (res_1_pdb_id, res_2_pdb_id)
            self.protocol.pymod.main_window.show_warning_message("Warning", txt)
            return None

        # If the two residues can be used to define a loop, then adds the new loop and updates
        # the frame with a new combobox row.
        new_ls_combo = User_loop_combo_qt(self, self.modeling_cluster.target.get_polymer_residues())

        # Activates the previous row and returns the name of the 2 selected residues.
        residues = self.list_of_loop_combos[-1].activate()

        # Finishes and adds the new row.
        self.list_of_loop_combos.append(new_ls_combo)

        # Adds the pair to the self.user_defined_loops, which is going to be
        # used in the launch_modelization() method.
        self.user_defined_loops.append(residues)


    def remove_user_loop(self, uls_to_remove):
        """
        This is called when the "Remove" button is pressed.
        """
        # Deactivate and get the right loop to remove.
        ls_to_remove = uls_to_remove.deactivate()
        # Finishes to adds the new row.
        self.list_of_loop_combos.remove(uls_to_remove)
        # Also removes the loop from the self.user_defined_loops.
        self.user_defined_loops.remove(ls_to_remove)


class User_loop_combo_qt:
    """
    Class for building in the Loop selection page in the modeling window a "row"
    with two comboboxes and a button to add or remove a new loops. This is very
    similar to the system used in the custom disulfide GUI.
    """

    # This is used in the constructor when a new combobox row is created.
    id_counter = 0

    def __init__(self, parent, polymer_residues):
        self.parent = parent

        # Selected have the "Add" button, unselected have the "Remove" button.
        self.selected = False

        User_loop_combo_qt.id_counter += 1
        self.row = User_loop_combo_qt.id_counter

        # The list of residues of the target sequence.
        self.polymer_residues = polymer_residues

        # The list of strings that is going to appear on the scrollable menus of the comboboxes.
        self.scrollable_res_list = []
        self.scrollable_res_dict = {}
        for res_i, res in enumerate(self.polymer_residues):
            res_label = "%s %s" % (res.three_letter_code, res.db_index)
            self.scrollable_res_list.append(res_label)
            self.scrollable_res_dict[res_label] = self.polymer_residues

        # Creates the first row with two comboboxes.
        self.create_combobox_row()


    combobox_padding = 30
    button_padding = 40

    def create_combobox_row(self):
        """
        Builds a row with two comboboxes an "Add" button for defining a loop.
        """

        # First residue combobox.
        self.res1_combobox = QtWidgets.QComboBox() # Select the starting residue:
        for item in self.scrollable_res_list:
            self.res1_combobox.addItem(item)
        self.res1_combobox.setEditable(False)
        self.res1_combobox.setStyleSheet(small_font_style)
        self.res1_combobox.setFixedWidth(self.res1_combobox.sizeHint().width() + self.combobox_padding)

        # Second residue combobox.
        self.res2_combobox = QtWidgets.QComboBox() # Select the ending residue:
        for item in self.scrollable_res_list:
            self.res2_combobox.addItem(item)
        self.res2_combobox.setEditable(False)
        self.res2_combobox.setStyleSheet(small_font_style)
        self.res2_combobox.setFixedWidth(self.res2_combobox.sizeHint().width() + self.combobox_padding)

        # "Add" button.
        self.new_loop_button = QtWidgets.QPushButton("Add")
        self.new_loop_button.setStyleSheet(small_font_style)
        self.new_loop_button.clicked.connect(self.press_add_button)
        self.new_loop_button.setFixedWidth(self.new_loop_button.sizeHint().width() + self.button_padding)

        self.parent.loop_selector_frame_layout.addWidget(self.res1_combobox, self.row, 0)
        self.parent.loop_selector_frame_layout.addWidget(self.res2_combobox, self.row, 1)
        self.parent.loop_selector_frame_layout.addWidget(self.new_loop_button, self.row, 2)

        User_loop_combo_qt.id_counter += 1


    def select_res(self, i):
        """
        This is launched by the combobox when some residue is selected.
        """
        pass


    def press_add_button(self):
        self.parent.add_new_user_loop()


    def activate(self):
        """
        This is going to return the information about which residue have been selected when the "Add"
        button is pressed.
        """
        self.selected = True

        # Removes the "Add" button.
        self.remove_widget(self.new_loop_button)
        self.new_loop_button = None

        # Removes both comboboxes, but before get their values.
        self.text1 = self.get_res1_combobox_val()
        self.remove_widget(self.res1_combobox)
        self.res1_combobox = None

        self.text2 = self.get_res2_combobox_val()
        self.remove_widget(self.res2_combobox)
        self.res2_combobox = None

        self.create_label_row()

        return self.text1, self.text2


    def remove_widget(self, widget):
        self.parent.loop_selector_frame_layout.removeWidget(widget)
        widget.deleteLater()


    def create_label_row(self):
        """
        Row with two labels and a "Remove" button.
        """

        # Create a first residue label that tells which residue has been selected.
        self.res1_label = QtWidgets.QLabel(self.text1)
        self.res1_label.setStyleSheet(small_font_style)
        self.res1_label.setAlignment(QtCore.Qt.AlignCenter)
        self.parent.loop_selector_frame_layout.addWidget(self.res1_label, self.row, 0)

        # Second residue label.
        self.res2_label = QtWidgets.QLabel(self.text2)
        self.res2_label.setStyleSheet(small_font_style)
        self.res2_label.setAlignment(QtCore.Qt.AlignCenter)
        self.parent.loop_selector_frame_layout.addWidget(self.res2_label, self.row, 1)

        # Adds the "Remove" button.
        self.remove_loop_button = QtWidgets.QPushButton("Remove")
        self.remove_loop_button.setStyleSheet(small_font_style)
        self.remove_loop_button.clicked.connect(self.press_remove_button)
        self.parent.loop_selector_frame_layout.addWidget(self.remove_loop_button, self.row, 2)


    def press_remove_button(self):
        self.parent.remove_user_loop(self)


    def deactivate(self):
        """
        This is going to return the information about which loop has been removed when "Remove"
        button is pressed.
        """
        self.selected = False

        self.remove_widget(self.res1_label)
        self.res1_label = None

        self.remove_widget(self.res2_label)
        self.res2_label = None

        self.remove_widget(self.remove_loop_button)
        self.remove_loop_button = None

        return self.text1, self.text2


    def get_res1_combobox_val(self):
        return self.res1_combobox.currentText()

    def get_res2_combobox_val(self):
        return self.res2_combobox.currentText()


################################################################################
# Options widgets.                                                             #
################################################################################

class Optimization_level_frame_loop_qt(Optimization_level_frame_qt):

    use_max_cg_iterations_enf = True
    use_vtfm_schedule_rds = False
    use_md_schedule_rds = True
    md_schedule_rds_default_val = "slow"
    use_repeat_optimization_enf = False
    use_max_obj_func_value_enf = False


class Custom_obj_func_frame_loop_qt(Custom_obj_func_frame_qt):

    use_loop_stat_pot_rds = True
    use_hddr_frame = False
    use_hdar_frame = True
    hdar_frame_text = "Dihedral angle restraints. Weight: "
    use_charmm22_frame = True
    use_soft_sphere_frame = True
    use_dope_frame = False
    use_stat_pot_frame = True
    use_lj_frame = True
    use_gbsa_frame = True

    initial_nb_cutoff_val = 7.0
