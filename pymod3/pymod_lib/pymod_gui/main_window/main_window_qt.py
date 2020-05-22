# Copyright 2020 by Giacomo Janson, Alessandro Paiardini. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
This module implements the main window of PyMod written in Qt. Other modules in this package provide
several main window widgets and mixins used to extend the class written in this file.
"""

import os
import sys

from pymol.Qt import QtWidgets, QtCore, QtGui

# Import the stylesheet.
from .aqua.qsshelper import QSSHelper

from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt, small_font_style
from ._main_menu_qt import PyMod_main_window_main_menu
from ._element_widgets_qt import PyMod_element_widgets_group_qt

import pymod_lib.pymod_vars as pmdt

from pymol import cmd


class PyMod_main_window_qt(QtWidgets.QMainWindow,
                           PyMod_main_window_main_menu):
    """
    A class for the PyQT PyMod main window.
    """

    is_pymod_window = True
    is_pymod_main_window = True

    def __init__(self, pymod, parent=None):
        super().__init__(parent)

        # Initial settings.
        self.pymod = pymod
        self.title = self.pymod.pymod_plugin_name + "." + self.pymod.pymod_revision
        self.statusBar_message = "Ready"
        self.left = 50
        self.top = 50
        self.width = 1100 # 2200, 1100
        self.height = 800 # 400
        self.bg_color = 'rgb(0, 0, 0)'
        self.text_color = "rgb(255, 255, 255)"

        # Fonts and styles.
        self.font_weight = "bold"
        self.font_size = 10 # 9, 10, 12
        self.font = "courier new"

        self.vertical_spacing = 1
        self.horizontal_spacing = 20

        # Creating a menu bar.
        self.make_main_menu() # create_menu_bar()

        # Creating central widget and setting it as central.
        self.central_widget = Centralwid(self)
        self.setCentralWidget(self.central_widget)
        self.update_qt_font()

        # Creating status bar.
        self.statusBar().showMessage(self.statusBar_message, 3000)
        self.statusBar().setSizeGripEnabled(1)

        # Initialize User Interface.
        self.initUI()

        # Set the layout.
        self.central_widget.id_form_layout.setFormAlignment(QtCore.Qt.AlignLeft)
        self.central_widget.id_form_layout.setVerticalSpacing(self.vertical_spacing)
        self.central_widget.sequence_ID_groupbox.setLayout(self.central_widget.id_form_layout)

        self.central_widget.seq_form_layout.setHorizontalSpacing(self.horizontal_spacing)
        self.central_widget.seq_form_layout.setFormAlignment(QtCore.Qt.AlignLeft)
        self.central_widget.seq_form_layout.setLabelAlignment(QtCore.Qt.AlignLeft)
        self.central_widget.seq_form_layout.setVerticalSpacing(self.vertical_spacing)
        self.central_widget.sequence_SEQ_groupbox.setLayout(self.central_widget.seq_form_layout)

        # Initializes a clipboard.
        self.clipboard = QtWidgets.QApplication.clipboard()


        #----------------
        # Key bindings. -
        #----------------

        # Select elements.
        self.select_all_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Ctrl+A"), self)
        self.select_all_shortcut.activated.connect(self.select_all_binding)

        self.show_all_structures_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Ctrl+S"), self)
        self.show_all_structures_shortcut.activated.connect(self.show_all_structures_binding)
        self.hide_all_structures_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Ctrl+H"), self)
        self.hide_all_structures_shortcut.activated.connect(self.hide_all_structures_binding)

        self.deselect_all_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Esc"), self)
        self.deselect_all_shortcut.activated.connect(self.deselect_all_binding)

        # Move elements.
        self.press_up_key_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("up"), self)
        self.press_up_key_shortcut.activated.connect(self.press_up_key_binding)

        self.press_down_key_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("down"), self)
        self.press_down_key_shortcut.activated.connect(self.press_down_key_binding)


    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        # Loads and sets the Qt stylesheet.
        module_path = sys.modules[__name__].__file__
        self.qss = QSSHelper.open_qss(os.path.join(os.path.dirname(module_path), 'aqua', 'aqua.qss'))
        self.setStyleSheet(self.qss)
        self.show()


    def closeEvent(self, evnt):
        if evnt.spontaneous():
            self.confirm_close(evnt)
        else:
            pass

    def confirm_close(self, evnt=None):
        """
        Asks confirmation when the main window is closed by the user.
        """
        answer = askyesno_qt("Exit PyMod?", "Are you really sure you want to exit PyMod?",
                             parent=self.pymod.get_qt_parent())
        if answer:
            self.pymod.inactivate_session()
            self.close() # Closes the PyMod main window.
        else:
            if evnt is not None:
                evnt.ignore()


    ###############################################################################################
    # Gridding system.                                                                            #
    ###############################################################################################

    def gridder(self, clear_selection=False, update_clusters=True, update_menus=True, update_elements=True):
        """
        Grids the PyMod elements (of both sequences and clusters) widgets in PyMod main window.
        When new elements are added to PyMod using the 'add_pymod_element_widgets' method of the
        'PyMod' class, their PyQT widgets are initializated, but they are not displayed in PyMod
        main window. This method actually displayes the widgets.
        """

        #---------------------------------------
        # Update clusters elements appearance. -
        #---------------------------------------

        if update_clusters:
            # First updates the cluster sequences and removes clusters with one or zero children.
            for cluster in self.pymod.get_cluster_elements():
                self.pymod.update_cluster_sequences(cluster)


        #----------------------------
        # Assigns the grid indices. -
        #----------------------------

        self.global_grid_row_index = 0
        self.global_grid_column_index = 0
        for pymod_element in self.pymod.root_element.get_children():
            self._set_descendants_grid_indices(pymod_element)


        #--------------------------------------------
        # Grids the widgets with their new indices. -
        #--------------------------------------------

        for pymod_element in self.pymod.root_element.get_children():
            self._grid_descendants(pymod_element, update_elements=update_elements)


        #---------------------------------------------
        # Updates other components of the PyMod GUI. -
        #---------------------------------------------

        if clear_selection:
            self.pymod.deselect_all_sequences()

        if update_menus:
            self.build_alignment_submenu()
            self.build_models_submenu()

        if self.pymod.DEVELOP:
            print("- Gridded. Row count id: %s. Row count seq: %s." % (self.central_widget.id_form_layout.rowCount(),
                                                                       self.central_widget.seq_form_layout.rowCount()))


    #################################################################
    # Set elements grid indices.                                    #
    #################################################################

    def _set_descendants_grid_indices(self, pymod_element):
        if pymod_element.is_mother():
            self.global_grid_column_index += 1
            self._set_element_grid_indices(pymod_element)
            for child_element in pymod_element.get_children():
                self._set_descendants_grid_indices(child_element)
            self.global_grid_column_index -= 1
        else:
            self._set_element_grid_indices(pymod_element)


    def _set_element_grid_indices(self, pymod_element):
        pymod_element.widget_group.old_grid_row_index = pymod_element.widget_group.grid_row_index

        if pymod_element.widget_group.show:
            pymod_element.widget_group.grid_row_index = self.global_grid_row_index
            self.global_grid_row_index += 1
        else:
            pymod_element.widget_group.grid_row_index = None


    #################################################################
    # Grid widgets in the PyMod main window.                        #
    #################################################################

    # Grid.
    def _grid_descendants(self, pymod_element, update_elements=False):
        # If the element is visible, grid it and update it (if necessary).
        if pymod_element.widget_group.show:
            # Shows/updates the widgets in PyMod main window.
            pymod_element.widget_group.grid_widgets(update_element_text=update_elements)
            # Proceed on to show the children of the element.
            if pymod_element.is_mother():
                for child_element in pymod_element.get_children():
                    self._grid_descendants(child_element, update_elements=update_elements)

        # If the element is hidden, only update it.
        else:
            if update_elements:
                self._update_widgets_recursively(pymod_element)

    # Update hidden widgets.
    def _update_widgets_recursively(self, pymod_element):
        self.update_widgets(pymod_element)
        if pymod_element.is_mother():
            for child_element in pymod_element.get_children():
                self._update_widgets_recursively(child_element)

    def update_widgets(self, pymod_element):
        pymod_element_widgets_group = pymod_element.widget_group
        pymod_element_widgets_group.sequence_text.update_text()
        pymod_element_widgets_group.header_entry.update_title()


    # Hide.
    def delete_element_widgets(self, pymod_element):
        """
        Remove the widgets of a PyMod element which has to be deleted.
        """
        pymod_element.widget_group.hide_widgets(save_status=True) # self.hide_element_widgets(pymod_element)


    ###############################################################################################
    # Handle PyMod data.                                                                          #
    ###############################################################################################

    def add_pymod_element_widgets(self, pymod_element):
        """
        Add widgets a 'pymod_element'. They will be displayed by default at the bottom of the main window.
        """
        pewp = PyMod_element_widgets_group_qt(main_window=self, pymod_element=pymod_element)
        pymod_element.widget_group = pewp


    ###############################################################################################
    # Key bindings.                                                                               #
    ###############################################################################################

    def select_all_binding(self):
        self.pymod.select_all_sequences()

    def deselect_all_binding(self):
        self.pymod.deselect_all_sequences()

    def show_all_structures_binding(self):
        self.pymod.show_all_structures_from_main_menu()

    def hide_all_structures_binding(self):
        self.pymod.hide_all_structures_from_main_menu()


    ####################################
    # Move elements up and down by one #
    # position in PyMod main window.   #
    ####################################

    def press_up_key_binding(self):
        self.move_elements_from_key_press("up")

    def press_down_key_binding(self):
        self.move_elements_from_key_press("down")


    def move_elements_from_key_press(self, direction):
        """
        Move 'up' or 'down' by a single position the selected elements in PyMod main window.
        This is not an efficient implementation, but it works. It should be rewritten to make it
        faster.
        """

        # Gets the elements to move.
        elements_to_move = self._get_elements_to_move()

        # Allow to move elements on the bottom of the list.
        if direction == "down":
            elements_to_move.reverse()

        # Temporarily adds 'None' elements to the list, so that multiple elements at the top or
        # bottom of container lists can be moved correctly.
        containers_set = set([e.mother for e in elements_to_move if not e.mother.selected]) # in elements_to_move
        for container in containers_set:
            container.list_of_children.append(None)
            container.list_of_children.insert(0, None)

        # Actually move the elements in their container lists.
        for element in elements_to_move:
            if not element.mother.selected:
                self.move_single_element(direction, element, element.mother.get_children())

        # Remove the 'None' values added before.
        for container in containers_set:
            container.list_of_children = [e for e in container.list_of_children if e != None]

        # Shows the the elements in the new order.
        if elements_to_move != []:
            self.gridder(update_clusters=False, update_menus=False, update_elements=False)


    def _get_elements_to_move(self):
        elements_to_move = []
        for e in self.pymod.get_selected_elements():
            # If the element is lead of a collapsed cluster, in order to move it in PyMod main
            # window, its mother will have to be moved.
            if not e in elements_to_move: # Prevents children from being included two times.
                elements_to_move.append(e)

        return elements_to_move # list(set(elements_to_move))


    def move_single_element(self, direction, element, container_list):
        """
        Move 'up' or 'down' by a single position a single element in a list.
        """
        change_index = 0
        old_index = container_list.index(element)
        if direction == "up":
            change_index -= 1
        elif direction == "down":
            # if old_index == len(container_list) - 1:
            #     return None
            change_index += 1
        self.pymod.change_pymod_element_list_index(element, old_index + change_index)


    def print_selected(self):
        ael = self.pymod.get_pymod_elements_list()
        sel = [i for i in ael if i.selected]
        print("\n# Selected elements: %s/%s." % (len(sel), len(ael)))
        for e in sel:
            print ("- %s. sel=%s; has_struct=%s." % (repr(e), e.selected, e.has_structure()))


    #################################################################
    # Check status of clusters.                                     #
    #################################################################

    def is_lead_of_collapsed_cluster(self, pymod_element):
        return pymod_element.is_child() and pymod_element.is_lead() and pymod_element.mother.widget_group._collapsed_cluster

    def is_collapsed_cluster(self, pymod_cluster):
        return pymod_cluster.is_cluster() and pymod_cluster.widget_group._collapsed_cluster


    ###############################################################################################
    # Color sequences in the main window.                                                         #
    ###############################################################################################

    #################################################################
    # Color the sequences and structures.                           #
    #################################################################

    def color_selection(self, mode, target_selection, color_scheme, regular_color=None):
        """
        Used to color a single sequence (and its structure) when "mode" is set to "single", to color
        mulitple sequences when "mode" is et to "multiple" or to color the list of the currently
        selected elements in the GUI if the mode is set to "selection".
        """

        # Builds a list of elements to be colored.
        elements_to_color = []
        if mode == "single":
            elements_to_color.append(target_selection)
        elif mode == "multiple":
            elements_to_color.extend(target_selection)
        elif mode == "selection":
            elements_to_color.extend(self.pymod.get_selected_sequences())
        elif mode == "all":
            elements_to_color.extend(self.pymod.get_all_sequences())

        # Actually proceeds to color the elements.
        selection_color_dict = {}
        for seq in elements_to_color:
            if color_scheme == "regular":
                self.color_element_by_regular_color(seq, regular_color)

            elif color_scheme == "polarity":
                color_dict = self.color_element_by_polarity(seq)

            elif color_scheme == "residue_type":
                color_dict = self.color_element_by_residue_type(seq)

            elif color_scheme == "secondary-observed":
                color_dict = self.color_element_by_obs_sec_str(seq)

            elif color_scheme == "secondary-predicted":
                color_dict = self.color_element_by_pred_sec_str(seq)

            # Colors elements with 3D structure according to the observed II str, elements with
            # predicted II str according to the prediction, and leaves the other elements unaltered.
            elif color_scheme == "secondary-auto":
                if seq.has_structure():
                    color_dict = self.color_element_by_obs_sec_str(seq)
                elif seq.has_predicted_secondary_structure():
                    color_dict = self.color_element_by_pred_sec_str(seq)

            elif color_scheme == "campo-scores":
                color_dict = self.color_element_by_campo_scores(seq)

            elif color_scheme == "entropy-scores":
                color_dict = self.color_element_by_entropy_scores(seq)

            elif color_scheme == "pair-conservation":
                color_dict = self.color_element_by_pc_scores(seq)

            elif color_scheme == "dope":
                color_dict = self.color_element_by_dope(seq)

            elif color_scheme == "domains":
                color_dict = self.color_element_by_domains(seq)

            elif color_scheme == "custom":
                color_dict = self.color_element_by_custom_scheme(seq)

            else:
                raise KeyError("Unknown color scheme: %s" % color_scheme)


    ##########################
    # Assigns color schemes. #
    ##########################

    def color_element_by_regular_color(self, element, color=None, color_structure=True):
        """
        Colors sequence by "regular" colors, that is, colors uniformly the sequence with some color.
        """
        element.color_by = "regular"
        if color != None:
            element.my_color = color
        return self.color_element(element, color_structure=color_structure)

    def color_element_by_polarity(self, element, color_structure=True):
        element.color_by = "polarity"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_residue_type(self, element, color_structure=True):
        element.color_by = "residue_type"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_obs_sec_str(self, element, color_structure=True):
        """
        Color elements by their observed secondary structure.
        """
        element.color_by = "secondary-observed"
        # If PyMOL has not been already used to assign sec str to this sequence.
        if not element.has_assigned_secondary_structure():
            self.pymod.assign_secondary_structure(element)
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_pred_sec_str(self, element, color_structure=True):
        """
        Colors according by secondary structure predicted by PSIPRED.
        """
        element.color_by = "secondary-predicted"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_campo_scores(self, element, color_structure=True):
        """
        Color by CAMPO scores.
        """
        element.color_by = "campo-scores"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_scr_scores(self, element, color_structure=True):
        """
        Color by CAMPO scores.
        """
        element.color_by = "scr-scores"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_entropy_scores(self, element, color_structure=True):
        """
        Color by CAMPO scores.
        """
        element.color_by = "entropy-scores"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_pc_scores(self, element, color_structure=True):
        """
        Color by pair conservation scores.
        """
        element.color_by = "pair-conservation"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_dope(self, element, color_structure=True):
        """
        Color by DOPE scores.
        """
        element.color_by = "dope"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_domains(self, element, color_structure=True):
        element.color_by = "domains"
        return self.color_element(element, color_structure=color_structure)


    def color_element_by_custom_scheme(self, element, color_structure=True, use_features=True, selected_feature=None):
        """
        Color by a custom color scheme (each residue is colored differently).
        """
        element.color_by = "custom"

        # Color by features.
        if use_features:
            # Show and color all the features of the element.
            if selected_feature is None:
                for feature in element.features["sequence"]:
                    feature.show = True
                for res_idx, res in enumerate(element.get_polymer_residues()):
                    res_features = res.features["sequence"]
                    if res_features:
                        # Sort the features by starting index, so that they can be colored properly.
                        for res_feature in sorted(res_features, key=lambda f: f.start):
                            res.custom_color = res_feature.color
                    else:
                        res.custom_color = res.get_default_color()
            # Show and color only the feature provided in the 'selected_feature' argument.
            else:
                selected_feature.show = True
                for res_idx, res in enumerate(element.get_polymer_residues()):
                    res_features = res.features["sequence"]
                    if res_features and selected_feature in res_features:
                        res.custom_color = selected_feature.color
                    else:
                        res.custom_color = res.get_default_color()
        else:
            pass

        return self.color_element(element, color_structure=color_structure)


    #################################################
    # Actually colors the sequences and structures. #
    #################################################

    def color_element(self, element, color_structure=True):
        """
        Colors the sequence entry when it is displayed by the 'gridder()' method or when the user
        changes the color scheme of a sequence. This can color also the PDB file of the element (if
        the element has one).
        """

        if element.color_by != "custom":
            # Hide the features of the element.
            for feature in element.features["sequence"]:
                feature.show = False
            # Stores the color scheme, so that the sequence can be colored with it again when the
            # "custom" coloring scheme has been applied.
            element.store_current_color()

        # Assignes a color to each residue.
        get_residue_color_seq = self.assign_residue_coloring_method(element, "sequence")
        get_residue_color_str = self.assign_residue_coloring_method(element, "structure")
        for r in element.get_polymer_residues():
            r.color_seq = get_residue_color_seq(r)
            r.color_str = get_residue_color_str(r)

        self.color_sequence_text(element)

        if color_structure and element.has_structure():
            self.color_structure(element)


    #######################
    # Coloring sequences. #
    #######################

    def color_sequence_text(self, element):
        """
        Colors the sequence entry in PyMod main window.
        """

        element.widget_group.sequence_text.update_text()


    ########################
    # Coloring structures. #
    ########################

    color_all_atoms_schemes = ("campo-scores", "scr-scores", "entropy-scores",
                               "pair-conservation")

    def color_structure(self, element, single_residue=False):
        """
        Colors the PDB structure of an element loaded into PyMOL.
        """
        chain_sel = element.get_pymol_selector()

        # Colors the structure according to some particular color scheme.
        if element.color_by != "regular":

            residues_to_color_dict = {}

            # Color residues one-by-one. Usually it it very slow in PyMOL.
            if single_residue:
                for res in element.get_polymer_residues():
                    color = res.color_str
                    residue_sel = res.get_pymol_selector()
                    cmd.color(color, residue_sel)
            # Color residues with the same color at the same time. Runs faster than the previous
            # method.
            else:
                for res in element.get_polymer_residues():
                    color = res.color_str # Gets the right color for the current residue.
                    if color in residues_to_color_dict:
                        residues_to_color_dict[color].append(res.db_index)
                    else:
                        residues_to_color_dict[color] = [res.db_index]

                for color in residues_to_color_dict:
                    cmd.color(color, chain_sel + " and resi " + self._join_residues_list(residues_to_color_dict[color]))

        # Colors all the residues of a structure with the same color.
        else:
            cmd.color(self.get_regular_structure_color(element.my_color), chain_sel)

        if not element.color_by in self.color_all_atoms_schemes:
            cmd.util.cnc(chain_sel) # Colors by atom.


    #--------------------------------
    # Compress PyMOL color strings. -
    #--------------------------------

    compress_color_strings = True

    def _join_residues_list(self, residues_ids):
        joined_list = "+".join([str(i) for i in residues_ids])
        if not self.compress_color_strings:
            return joined_list
        else:
            return self._compress_pymol_color_string(joined_list)

    def _compress_pymol_color_string(self, color_string):
        compressed_color_list = []
        for i, n in enumerate(color_string.split("+")):
            if i == 0:
                compressed_color_list.append([int(n)])
            else:
                if compressed_color_list[-1][-1] == int(n)-1:
                    compressed_color_list[-1].append(int(n))
                else:
                    compressed_color_list.append([int(n)])
        return "+".join([self._get_color_unit_string(e) for e in compressed_color_list])

    def _get_color_unit_string(self, indices_list):
        if len(indices_list) == 1:
            return str(indices_list[0])
        elif len(indices_list) == 2:
            return "%s-%s" % (indices_list[0], indices_list[-1])
        elif len(indices_list) > 2:
            return "%s-%s" % (indices_list[0], indices_list[-1])
        else:
            raise ValueError("Wrong 'indices_list' length: %s" % len(indices_list))


    ####################################
    # Get the colors for the residues. #
    ####################################

    def assign_residue_coloring_method(self, element, color_target):

        if element.color_by == "polarity":
            if color_target == "sequence":
                return self.get_polarity_sequence_color
            elif color_target == "structure":
                return self.get_polarity_structure_color
        if element.color_by == "residue_type":
            if color_target == "sequence":
                return self.get_residue_type_sequence_color
            elif color_target == "structure":
                return self.get_residue_type_structure_color
        elif element.color_by == "secondary-observed":
            if color_target == "sequence":
                return self.get_observed_sec_str_sequence_color
            elif color_target == "structure":
                return self.get_observed_sec_str_structure_color
        elif element.color_by == "secondary-predicted":
            if color_target == "sequence":
                return self.get_predicted_sec_str_sequence_color
            elif color_target == "structure":
                return self.get_predicted_sec_str_structure_color
        elif element.color_by == "campo-scores":
            if color_target == "sequence":
                return self.get_campo_sequence_color
            elif color_target == "structure":
                return self.get_campo_structure_color
        elif element.color_by == "scr-scores":
            if color_target == "sequence":
                return self.get_scr_sequence_color
            elif color_target == "structure":
                return self.get_scr_structure_color
        elif element.color_by == "entropy-scores":
            if color_target == "sequence":
                return self.get_entropy_sequence_color
            elif color_target == "structure":
                return self.get_entropy_structure_color
        elif element.color_by == "pair-conservation":
            if color_target == "sequence":
                return self.get_pc_sequence_color
            elif color_target == "structure":
                return self.get_pc_structure_color
        elif element.color_by == "dope":
            if color_target == "sequence":
                return self.get_dope_sequence_color
            elif color_target == "structure":
                return self.get_dope_structure_color
        elif element.color_by == 'domains':
            if color_target == "sequence":
                return self.get_domain_sequence_color
            elif color_target == "structure":
                return self.get_domain_structure_color
        elif element.color_by == "custom":
            if color_target == "sequence":
                return self.get_custom_sequence_color
            elif color_target == "structure":
                return self.get_custom_structure_color
        elif element.color_by == "regular":
            if color_target == "sequence":
                return self._get_regular_sequence_color
            elif color_target == "structure":
                return self._get_regular_structure_color
        else:
            raise KeyError("Unknown coloring method: %s." % element.color_by)


    # Regular colors.
    def _get_regular_sequence_color(self, res):
        return self.get_regular_sequence_color(res.pymod_element.my_color)

    def _get_regular_structure_color(self, res):
        return self.get_regular_structure_color(res.pymod_element.my_color)

    def get_regular_sequence_color(self, color):
        return self.pymod.all_colors_dict_tkinter[color]

    def get_regular_structure_color(self, color):
        return color


    # Residue polarity colors.
    def get_polarity_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter.get(self.form_residue_polarity_color_name(residue), "#ffffff")

    def get_polarity_structure_color(self, residue):
        return self.form_residue_polarity_color_name(residue)

    def form_residue_polarity_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_polarity_color_name, residue.one_letter_code)


    # Residue type colors.
    def get_residue_type_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter.get(self.form_residue_type_color_name(residue), "#ffffff")

    def get_residue_type_structure_color(self, residue):
        return self.form_residue_type_color_name(residue)

    def form_residue_type_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_residue_type_color_name, residue.one_letter_code)


    # Observed secondary structure colors.
    def get_observed_sec_str_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_observed_sec_str_color_name(residue)]

    def get_observed_sec_str_structure_color(self, residue):
        return self.form_observed_sec_str_color_name(residue)

    def form_observed_sec_str_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_obs_sec_str_name, residue.secondary_structure)


    # Predicted secondary structure colors.
    def get_predicted_sec_str_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_predicted_sec_str_color_name(residue)]

    def get_predicted_sec_str_structure_color(self, residue):
        return self.form_predicted_sec_str_color_name(residue)

    def form_predicted_sec_str_color_name(self, residue):
        return "%s_%s_%s" % (pmdt.pymol_psipred_color_name, residue.psipred_result["confidence"], residue.psipred_result["sec-str-element"])


    # CAMPO colors.
    def get_campo_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_campo_color_name(residue)]

    def get_campo_structure_color(self, residue):
        return self.form_campo_color_name(residue)

    def form_campo_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_campo_color_name, residue.campo_score["interval"])


    # SCR colors.
    def get_scr_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_scr_color_name(residue)]

    def get_scr_structure_color(self, residue):
        return self.form_scr_color_name(residue)

    def form_scr_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_scr_color_name, residue.scr_score["interval"])


    # Entropy scores.
    def get_entropy_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_entropy_color_name(residue)]

    def get_entropy_structure_color(self, residue):
        return self.form_entropy_color_name(residue)

    def form_entropy_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_entropy_color_name, residue.entropy_score["interval"])


    # Pair conservation colors.
    def get_pc_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_pc_color_name(residue)]

    def get_pc_structure_color(self, residue):
        return self.form_pc_color_name(residue)

    def form_pc_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_pc_color_name, residue.pc_score)


    # DOPE colors.
    def get_dope_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_dope_color_name(residue)]

    def get_dope_structure_color(self, residue):
        return self.form_dope_color_name(residue)

    def form_dope_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_dope_color_name, residue.dope_score["interval"])


    # Custom colors.
    def get_custom_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_custom_color_name(residue)]

    def get_custom_structure_color(self, residue):
        return self.form_custom_color_name(residue)

    def form_custom_color_name(self, residue):
        return residue.custom_color


    # Domains colors.
    def get_domain_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_domain_color_name(residue)]

    def get_domain_structure_color(self, residue):
        return self.form_domain_color_name(residue)

    def form_domain_color_name(self, residue):
        res_domains = residue.features["domains"]
        if len(res_domains) == 0:
            return 'grey70'
        elif len(res_domains) == 1:
            return res_domains[0].domain_color[0]
        else: # More than one domain per residue.
            return 'teal'


    #################################################################
    # Sequences font.                                               #
    #################################################################

    def update_font(self, new_font_type=None, new_font_size=None):
        """
        Updates the font of the sequences displayed in the main window. Called from the
        'Font Selection' menu items from the main window.
        """
        if new_font_type:
            self.font = new_font_type
        if new_font_size:
            self.font_size = int(new_font_size)

        id_frame_stylesheet = "QLabel {font: %spt %s; font-weight: %s; color: white}" % (self.font_size, self.font, self.font_weight)
        self.central_widget.sequence_ID_groupbox.setStyleSheet(id_frame_stylesheet)
        seq_frame_stylesheet = "QLabel {font: %spt %s; font-weight: %s; color: white}" % (self.font_size, self.font, self.font_weight)
        self.central_widget.sequence_SEQ_groupbox.setStyleSheet(seq_frame_stylesheet)
        self.update_qt_font()


    def update_qt_font(self):
        """
        Updates the font of the PyMod Qt main window.
        """
        self.qfont = QtGui.QFont(self.font, self.font_size)
        self.fm = QtGui.QFontMetrics(self.qfont)


    ###############################################################################################
    # Messages.                                                                                   #
    ###############################################################################################

    def show_info_message(self, title_to_show, message_to_show, parent="auto"):
        _parent = self.pymod.get_qt_parent() if parent == "auto" else parent
        QtWidgets.QMessageBox.information(_parent, title_to_show, message_to_show)

    def show_warning_message(self, title_to_show, message_to_show, parent="auto"):
        _parent = self.pymod.get_qt_parent() if parent == "auto" else parent
        QtWidgets.QMessageBox.warning(_parent, title_to_show, message_to_show)

    def show_error_message(self, title_to_show, message_to_show, parent="auto"):
        _parent = self.pymod.get_qt_parent() if parent == "auto" else parent
        QtWidgets.QMessageBox.critical(_parent, title_to_show, message_to_show)


class Centralwid(QtWidgets.QWidget):
    """
    Central PyQt window with 2 (left and right) frames and a kind of status bar.
    """

    def __init__(self, main_window):
        super().__init__()
        self.style = "background-color: rgb(0, 0, 0); color: rgb(255, 255, 255); font-weight: bold"
        self.main_window = main_window
        self.initUI()


    def initUI(self):

        #----------------------------
        # Left frame (for headers). -
        #----------------------------

        self.sequence_ID_groupbox = QtWidgets.QGroupBox('SEQUENCE ID')
        # self.sequence_ID_groupbox.setStyleSheet("QLabel {font: 14pt COURIER NEW font-weight: bold} ")
        id_frame_stylesheet = "QLabel {font: %spt %s; font-weight: %s; color: white}" % (self.main_window.font_size, self.main_window.font, self.main_window.font_weight)
        self.sequence_ID_groupbox.setStyleSheet(id_frame_stylesheet)

        self.id_form_layout = QtWidgets.QFormLayout()
        #self.left_scroll
        self.left_scroll = QtWidgets.QScrollArea()
        self.left_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.left_scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.left_scroll.resize(200, 400)
        self.left_scroll.setWidget(self.sequence_ID_groupbox) # sequence_ID_groupbox dentro left_scroll area
        self.left_scroll.setWidgetResizable(True)
        #self.left_frame
        self.left_frame = QtWidgets.QFrame(self)
        self.left_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.left_frame.resize(200, 400)
        self.left_frame.setStyleSheet(self.style)
        self.left_frame.setFrameShadow(QtWidgets.QFrame.Sunken)
        #self.left_frame_layout
        self.left_frame_layout = QtWidgets.QVBoxLayout(self)
        self.left_frame_layout.addWidget(self.left_scroll)
        self.left_frame.setLayout(self.left_frame_layout) # left_frame_layout dentro left_frame


        #-------------------------------
        # Right frame (for sequences). -
        #-------------------------------

        # This groupbox
        self.sequence_SEQ_groupbox = QtWidgets.QGroupBox('SEQUENCES')
        seq_frame_stylesheet = "QLabel {font: %spt %s; font-weight: %s; color: white}" % (self.main_window.font_size, self.main_window.font, self.main_window.font_weight)
        self.sequence_SEQ_groupbox.setStyleSheet(seq_frame_stylesheet)
        self.seq_form_layout = QtWidgets.QFormLayout()

        #self.right_scroll
        self.right_scroll = QtWidgets.QScrollArea()
        self.right_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.right_scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.right_scroll.resize(900, 400)
        self.right_scroll.setWidget(self.sequence_SEQ_groupbox) # sequence_ID_groupbox dentro left_scroll area
        self.right_scroll.setWidgetResizable(True)
        #self.right_frame
        self.right_frame = QtWidgets.QFrame(self)
        self.right_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.right_frame.resize(900, 400)
        self.right_frame.setStyleSheet(self.style)
        self.right_frame.setFrameShadow(QtWidgets.QFrame.Sunken)
        #self.right_frame_layout
        self.right_frame_layout = QtWidgets.QVBoxLayout(self)
        self.right_frame_layout.addWidget(self.right_scroll)
        self.right_frame.setLayout(self.right_frame_layout) # left_frame_layout dentro left_frame

        #connect the two Vertical Bars to move them togheter
        self.left_scroll.verticalScrollBar().valueChanged.connect(self.right_scroll.verticalScrollBar().setValue)
        self.right_scroll.verticalScrollBar().valueChanged.connect(self.left_scroll.verticalScrollBar().setValue)


        #----------------------------------
        # Bottom part of the main window. -
        #----------------------------------

        self.splitter1 = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.splitter1.addWidget(self.left_frame)
        self.splitter1.addWidget(self.right_frame)

        # creating sequence and position labels
        self.label_sequence = QtWidgets.QLabel(self)
        self.label_sequence.setText('Sequence:')
        self.label_sequence.setStyleSheet(small_font_style)
        self.textbox_sequence = QtWidgets.QLineEdit(self)
        self.textbox_sequence.setStyleSheet(self.style + "; " + small_font_style)
        self.textbox_sequence.setReadOnly(True)

        self.label_position = QtWidgets.QLabel(self)
        self.label_position.setText('Position:')
        self.label_position.setStyleSheet(small_font_style)
        self.textbox_position = QtWidgets.QLineEdit(self)
        self.textbox_position.setReadOnly(True)
        self.textbox_position.setStyleSheet(self.style + "; " + small_font_style)
        self.textbox_position.setMinimumWidth(675) # Width of the residues message bar width.

        # creating an horizontal layout with sequence and position labels
        self.text_layout = QtWidgets.QHBoxLayout()
        self.text_layout.addWidget(self.label_sequence)
        self.text_layout.addWidget(self.textbox_sequence)
        self.text_layout.addWidget(self.label_position)
        self.text_layout.addWidget(self.textbox_position)

        # creating a layout with sequence window and labels
        self.grid = QtWidgets.QVBoxLayout()
        self.grid.addWidget(self.splitter1)
        self.grid.addLayout(self.text_layout)
        self.setLayout(self.grid)
