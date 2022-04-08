# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing the widgets showed in the left pane of PyMod main window.
"""

from pymol.Qt import QtWidgets, QtCore, QtGui

from pymod_lib.pymod_gui.shared_gui_components_qt import add_qt_menu_command, open_color_dialog
from pymod_lib.pymod_seq.seq_manipulation import remove_gap_only_columns
import pymod_lib.pymod_vars as pmdt

from pymol import cmd


###################################################################################################
# Widget for the header label of an element.                                                      #
###################################################################################################

class MyQLabel_header(QtWidgets.QLabel):
    """
    A custom QLabel for headers of PyMod elements.
    """

    def __init__(self, parent_group):
        QtWidgets.QLabel.__init__(self, parent_group.pymod_element.my_header)
        self.parent_group = parent_group

        self.resize_to_content()

        # Set the style.
        self.setStyleSheet("color: red")

        self.set_default_cursor()


    def get_main_window(self):
        return self.parent_group.main_window

    # def get_pymod_element(self):
    #     return self.parent_group.pymod_element

    def set_default_cursor(self):
        # self.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))


    def resize_to_content(self):
        self.pixelsWide = self.get_main_window().fm.width(self.text())
        self.pixelsHigh = self.get_main_window().fm.height()


    def update_title(self):
        self.setText(self.parent_group.pymod_element.my_header)


    ###############################################################################################
    # Bindings for mouse events.                                                                  #
    ###############################################################################################

    def mousePressEvent(self, event):
        """
        Left click on an element header.
        """

        # self.setCursor(QtGui.QCursor(QtCore.Qt.OpenHandCursor))
        if event.buttons() == QtCore.Qt.LeftButton:
            self.parent_group.toggle_element()
        elif event.buttons() == QtCore.Qt.MiddleButton:
            self.click_structure_with_middle_button()


    def mouseMoveEvent(self, event):
        """
        Select/Unselect sequences hovering on their headers while the mouse left button is pressed.
        """

        # Only activates when the left button is being pressed.
        if event.buttons() == QtCore.Qt.LeftButton:
            highlighted_widget = QtWidgets.QApplication.widgetAt(self.mapToGlobal(event.pos()))

            # Checks if the widget hovered with the mouse belongs to the 'Header_entry' class and is not
            # the entry originally clicked.
            if isinstance(highlighted_widget, MyQLabel_header) and highlighted_widget != self:
                starting_element_state = self.parent_group.pymod_element.selected
                if starting_element_state and not highlighted_widget.parent_group.pymod_element.selected:
                    highlighted_widget.parent_group.toggle_element()
                elif not starting_element_state and highlighted_widget.parent_group.pymod_element.selected:
                    highlighted_widget.parent_group.toggle_element()


    def enterEvent(self, event):
        """
        Show information on the message bar of the main window.
        """
        self.parent_group.main_window.central_widget.textbox_sequence.setText(str(self.parent_group.pymod_element.my_header))

    def leaveEvent(self, event):
        """
        Called when the mouse leaves a header label.
        """
        self.parent_group.main_window.central_widget.textbox_sequence.setText("")

    def contextMenuEvent(self, event):
        """
        Shows the context menu when right clickings on the header of an element.
        """
        self.build_header_popup_menu()
        self.context_menu.exec_(self.mapToGlobal(event.pos()))


    ###############################################################################################
    # Context menu for the header entry.                                                          #
    ###############################################################################################

    def build_header_popup_menu(self):
        """
        Builds the popup menu that appears when users left-clicks with on the sequence header in the
        main window left pan.
        """

        self.context_menu = Header_context_menu(self)

        # Places the single sequence submenus in a separate submenu in order to distinguish it
        # from the selection submenu.
        if len(self.parent_group.pymod.get_selected_sequences()) > 1:
            self.single_sequence_context_submenu = QtWidgets.QMenu(self.parent_group.pymod_element.my_header, self.context_menu)
            self.context_menu.addMenu(self.single_sequence_context_submenu)
            self.single_element_target_submenu = self.single_sequence_context_submenu
        # Places the single sequence submenus in the main context menu.
        else:
            add_qt_menu_command(self.context_menu, self.parent_group.pymod_element.my_header, None)
            self.single_element_target_submenu = self.context_menu
        self.context_menu.addSeparator()

        # Builds a popup menu for sequence elements.
        if not self.parent_group.pymod_element.is_cluster():
            self.build_single_sequence_header_popup_menu()
        # For cluster elements (alignments or blast-search elements).
        else:
            self.build_cluster_popup_menu(self.single_element_target_submenu, mode="cluster", extra_spacer=True)

        self.update_left_popup_menu()


    def update_left_popup_menu(self):
        """
        Activates the "Selection" item when at least two elements are selected.
        In order to make this work the "Selection" item always has to be in the last position in all
        kind of menus.
        """
        if len(self.parent_group.pymod.get_selected_sequences()) > 1:
            # self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=NORMAL)
            self.build_selection_menu()
        elif not self.parent_group.pymod_element.is_cluster():
            # self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=DISABLED)
            pass


    def build_single_sequence_header_popup_menu(self):

        # Build the "Sequence" menu.
        self.build_sequence_menu()

        self.single_element_target_submenu.addSeparator()

        # Build the "Color" menu.
        self.build_color_menu()
        self.single_element_target_submenu.addSeparator()

        # Build the "Structure" menu.
        self.build_structure_menu()
        self.single_element_target_submenu.addSeparator()

        # Build the "Domains" menu.
        if self.parent_group.pymod_element.has_domains(only_original=True):
            self.build_domains_menu()
            self.single_element_target_submenu.addSeparator()

        # Build the "Features" menu.
        if self.parent_group.pymod_element.has_features():
            self.builds_features_menu()
            self.single_element_target_submenu.addSeparator()

        # Build the "Cluster Options" menu.
        if self.parent_group.pymod_element.is_child():
            self.build_cluster_options_menu()
            self.single_element_target_submenu.addSeparator()


    def build_sequence_menu(self):
        """
        Submenu with options for manipulating a sequence loaded in PyMod.
        """

        self.sequence_context_submenu = QtWidgets.QMenu('Sequence', self.single_element_target_submenu)
        self.single_element_target_submenu.addMenu(self.sequence_context_submenu)

        add_qt_menu_command(self.sequence_context_submenu, "Save Sequence to File", self.save_sequence_from_left_pane)
        add_qt_menu_command(self.sequence_context_submenu, "Copy Sequence to Clipboard", self.copy_sequence_to_clipboard)
        self.sequence_context_submenu.addSeparator()

        add_qt_menu_command(self.sequence_context_submenu, "Edit Sequence", self.edit_sequence_from_context)
        add_qt_menu_command(self.sequence_context_submenu, "Search Sub-sequence", self.search_string_from_context)
        self.sequence_context_submenu.addSeparator()

        add_qt_menu_command(self.sequence_context_submenu, "Duplicate Sequence", self.duplicate_sequence_from_the_left_pane)
        add_qt_menu_command(self.sequence_context_submenu, "Delete Sequence", self.delete_sequence_from_the_left_pane)


    def build_color_menu(self):
        """
        Color submenu containing all the option to color for a single sequence.
        """

        self.color_context_submenu = QtWidgets.QMenu('Color', self.single_element_target_submenu)
        self.single_element_target_submenu.addMenu(self.color_context_submenu)

        # A submenu to choose a single color used to color all the residues of a sequence.
        self.regular_colors_context_submenu = QtWidgets.QMenu('Color whole Sequence by', self.color_context_submenu)
        self.color_context_submenu.addMenu(self.regular_colors_context_submenu)
        self.build_regular_colors_submenu(self.regular_colors_context_submenu, "single")

        # Colors each kind of residue in a sequence in a different way.
        self.residues_colors_context_submenu = QtWidgets.QMenu('By Residue Properties', self.color_context_submenu)
        self.color_context_submenu.addMenu(self.residues_colors_context_submenu)
        add_qt_menu_command(self.residues_colors_context_submenu, "Residue Type", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "residue_type"))
        add_qt_menu_command(self.residues_colors_context_submenu, "Polarity", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "polarity"))


        # Secondary structure colors.
        if self.parent_group.can_be_colored_by_secondary_structure():
            # self.color_context_submenu.addSeparator()
            self.sec_str_color_submenu = QtWidgets.QMenu('By Secondary Structure', self.color_context_submenu)
            self.color_context_submenu.addMenu(self.sec_str_color_submenu)
            if self.parent_group.pymod_element.has_structure():
                add_qt_menu_command(self.sec_str_color_submenu, "Observed", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "secondary-observed"))
            if self.parent_group.pymod_element.has_predicted_secondary_structure():
                add_qt_menu_command(self.sec_str_color_submenu, "Predicted by PSIPRED", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "secondary-predicted"))


        # Conservation colors.
        if self.parent_group.pymod_element.has_campo_scores() or self.parent_group.pymod_element.has_entropy_scores():
            # self.color_context_submenu.addSeparator()
            self.conservation_colors_menu = QtWidgets.QMenu('By Conservation', self.color_context_submenu)
            self.color_context_submenu.addMenu(self.conservation_colors_menu)
            if self.parent_group.pymod_element.has_campo_scores():
                add_qt_menu_command(self.conservation_colors_menu, "CAMPO scores", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "campo-scores"))
            if self.parent_group.pymod_element.has_entropy_scores():
                add_qt_menu_command(self.conservation_colors_menu, "Entropy scores", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "entropy-scores"))


        # Energy colors.
        if self.parent_group.pymod_element.has_dope_scores():
            # self.color_context_submenu.addSeparator()
            self.energy_colors_menu = QtWidgets.QMenu('By Energy', self.color_context_submenu)
            self.color_context_submenu.addMenu(self.energy_colors_menu)
            add_qt_menu_command(self.energy_colors_menu, "DOPE scores", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "dope"))


        # Color by domain.
        if self.parent_group.pymod_element.has_domains():
            add_qt_menu_command(self.color_context_submenu, "By Domain", lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "domains"))


    def build_regular_colors_submenu(self, target_menu, color_mode, elements_to_color=None):
        if color_mode == "single":
            elements_to_color = self.parent_group.pymod_element

        # Regular color scheme of the element.
        add_qt_menu_command(target_menu, "Default Color", lambda: self.parent_group.main_window.color_selection(color_mode, elements_to_color, "regular"))

        # Build PyMOL color palette menu.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMOL Colors", pmdt.pymol_regular_colors_list)
        '''
        # Build PyMOL light color pallette.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMOL Light Colors", pmdt.pymol_light_colors_list)
        '''
        # Build PyMod color palette.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMod Colors", pmdt.pymod_regular_colors_list)

        target_menu.addSeparator()

        # Custom color selection.
        add_qt_menu_command(target_menu, "Pick Color", lambda a=None, m=color_mode, e=elements_to_color: self.open_color_dialog_from_context(m, e))


    def build_color_palette_submenu(self, color_mode, elements_to_color, target_submenu, title, list_of_colors):
        new_palette_submenu = QtWidgets.QMenu(title, target_submenu)
        target_submenu.addMenu(new_palette_submenu)

        for color in list_of_colors:
            add_qt_menu_command(parent=new_palette_submenu, label=color,
                                command=lambda a=None, c=color: self.parent_group.main_window.color_selection(color_mode, elements_to_color, "regular", c),
                                fg_color=self.parent_group.main_window.get_regular_sequence_color(color),
                                bg_color="#ABABAB")


    def open_color_dialog_from_context(self, color_mode, elements_to_color):
        """
        Gets from a 'QColorDialog' widget an HEX color string and colors a PyMod selection with it.
        """

        selected_color = open_color_dialog(color_format="all")
        if selected_color is not None:
            color_rgb, color_hex = selected_color
            new_color_name = self.parent_group.pymod.add_new_color(color_rgb, color_hex)
            self.parent_group.main_window.color_selection(color_mode, elements_to_color, "regular", new_color_name)


    def build_structure_menu(self):
        """
        Submenu for elements that have a structure loaded in PyMOL.
        """

        self.structure_context_submenu = QtWidgets.QMenu('Structure', self.single_element_target_submenu)

        if self.parent_group.pymod_element.has_structure():
            self.single_element_target_submenu.addMenu(self.structure_context_submenu)
            # add_qt_menu_command(self.structure_context_submenu, "PDB Chain Information", self.show_structure_info)
            add_qt_menu_command(self.structure_context_submenu, "Save PDB Chain to File", self.save_structure_from_left_pane)
            self.structure_context_submenu.addSeparator()
            add_qt_menu_command(self.structure_context_submenu, "Center Chain in PyMOL", self.center_chain_in_pymol_from_header_entry)
            # TODO: a switch would be nice.
            add_qt_menu_command(self.structure_context_submenu, "Show Chain in PyMOL", self.show_chain_in_pymol_from_header_entry)
            add_qt_menu_command(self.structure_context_submenu, "Hide Chain in PyMOL", self.hide_chain_in_pymol_from_header_entry)
            self.structure_context_submenu.addSeparator()
            add_qt_menu_command(self.structure_context_submenu, "Show Chain as Hedgehog", self.show_hedgehog_in_pymol_from_header_entry)
            add_qt_menu_command(self.structure_context_submenu, "Show Heteroatoms", self.show_het_in_pymol_from_header_entry)
            add_qt_menu_command(self.structure_context_submenu, "Hide Heteroatoms", self.hide_het_in_pymol_from_header_entry)
            self.single_element_target_submenu.addMenu(self.structure_context_submenu)
        else:
            if self.parent_group.pymod_element.pdb_is_fetchable():
                add_qt_menu_command(self.structure_context_submenu, "Fetch PDB File", lambda: self.parent_group.pymod.fetch_pdb_files("single", self.parent_group.pymod_element))
                # self.structure_context_submenu.addSeparator()
                # add_qt_menu_command(self.structure_context_submenu, "Associate 3D Structure", lambda: self.parent_group.pymod.associate_structure_from_popup_menu(self.parent_group.pymod_element))
                self.single_element_target_submenu.addMenu(self.structure_context_submenu)


    def build_domains_menu(self):
        self.domains_context_submenu = QtWidgets.QMenu('Domains', self.single_element_target_submenu)
        add_qt_menu_command(self.domains_context_submenu, "Split into Domains", lambda a=None, pe=self.parent_group.pymod_element: self.parent_group.pymod.launch_domain_splitting(pe))
        if self.parent_group.pymod_element.derived_domains_list:
            add_qt_menu_command(self.domains_context_submenu, "Fuse Domains Alignments", lambda a=None, pe=self.parent_group.pymod_element: self.parent_group.pymod.launch_domain_fuse(pe))
        self.single_element_target_submenu.addMenu(self.domains_context_submenu)


    def builds_features_menu(self):
        self.features_context_submenu = QtWidgets.QMenu('Features', self.single_element_target_submenu)
        add_qt_menu_command(self.features_context_submenu, "Show Features",
            lambda: self.parent_group.main_window.color_selection("single", self.parent_group.pymod_element, "custom"))
        add_qt_menu_command(self.features_context_submenu, "Delete Features",
            lambda: self.parent_group.pymod.delete_features_from_context_menu(self.parent_group.pymod_element))
        self.single_element_target_submenu.addMenu(self.features_context_submenu)


    def build_cluster_options_menu(self):
        """
        Submenu with options to manage a sequence within its cluster.
        """

        self.cluster_context_submenu = QtWidgets.QMenu('Cluster Options', self.single_element_target_submenu)
        self.single_element_target_submenu.addMenu(self.cluster_context_submenu)

        add_qt_menu_command(self.cluster_context_submenu, "Extract Sequence from Cluster", self.extract_from_cluster)
        # self.cluster_context_submenu.addSeparator()
        # if not self.parent_group.pymod_element.is_lead():
        #     add_qt_menu_command(self.cluster_context_submenu, "Make Cluster Lead", self.make_lead_from_left_menu)
        # else:
        #     add_qt_menu_command(self.cluster_context_submenu, "Remove Cluster Lead", self.remove_lead_from_left_menu)


    #######################################
    # Multiple elements (selection) menu. #
    #######################################

    def build_selection_menu(self):
        """
        Submenu with options for managing a selection.
        """

        self.selection_context_submenu = QtWidgets.QMenu('Selection', self.context_menu)
        self.context_menu.addMenu(self.selection_context_submenu)

        # Build the "Sequence" menu.
        self.build_selection_sequence_menu()
        self.selection_context_submenu.addSeparator()

        # Build the "Color" menu.
        self.build_selection_color_menu()

        # Build the "Structure" menu.
        if self.parent_group.pymod.all_sequences_have_structure() or self.parent_group.pymod.all_sequences_have_fetchable_pdbs():
            self.selection_context_submenu.addSeparator()
            self.build_selection_structure_menu()

        # Build the "Cluster" menu.
        if self.parent_group.pymod.all_sequences_are_children():
            self.selection_context_submenu.addSeparator()
            self.build_selection_cluster_menu()


    def build_selection_sequence_menu(self):

        self.selection_sequence_context_submenu = QtWidgets.QMenu('Sequences', self.selection_context_submenu)
        self.selection_context_submenu.addMenu(self.selection_sequence_context_submenu)

        add_qt_menu_command(self.selection_sequence_context_submenu, "Save Selection to File", self.save_selection_from_left_pane)
        add_qt_menu_command(self.selection_sequence_context_submenu, "Copy Selection to Clipboard", self.copy_selection)
        self.selection_sequence_context_submenu.addSeparator()
        add_qt_menu_command(self.selection_sequence_context_submenu, "Duplicate Selection", self.duplicate_selection)
        add_qt_menu_command(self.selection_sequence_context_submenu, "Delete Selection", self.delete_many_sequences)


    def build_selection_color_menu(self):

        self.selection_color_context_submenu = self.build_multiple_color_menu(mode="selection")
        self.selection_context_submenu.addMenu(self.selection_color_context_submenu)


    def build_multiple_color_menu(self, mode, cluster_target_menu=None):
        """
        Used to build the color menu of both Selection and cluster elements popup menus.
        """

        if mode == "selection":
            target_menu = self.context_menu
            color_selection_mode = "selection"
            color_selection_target = None
            sequences_list = self.parent_group.pymod.get_selected_sequences()
            # color_selection_target = sequences_list
            color_target_label = "Selection"

        elif mode == "cluster":
            target_menu = cluster_target_menu
            color_selection_mode = "multiple"
            if self.parent_group.pymod_element.is_cluster():
                color_selection_target = self.parent_group.pymod_element.get_children()
                sequences_list = self.parent_group.pymod_element.get_children()
            else:
                color_selection_target = self.parent_group.pymod_element.mother.get_children()
                sequences_list = self.parent_group.pymod_element.mother.get_children()
            color_target_label = "Cluster"


        # Builds the selection color menu.
        multiple_color_menu = QtWidgets.QMenu('Color', target_menu)

        # A submenu to choose a single color used to color all the residues of a sequence.
        multiple_regular_colors_menu = QtWidgets.QMenu("Color whole %s by" % (color_target_label), multiple_color_menu) # Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.build_regular_colors_submenu(multiple_regular_colors_menu, color_selection_mode, elements_to_color=color_selection_target)
        multiple_color_menu.addMenu(multiple_regular_colors_menu)
        # multiple_color_menu.addSeparator()

        # Colors each kind of residue in a sequence in a different way.
        multiple_residues_colors_menu = QtWidgets.QMenu('By Residue Properties', target_menu)
        add_qt_menu_command(multiple_residues_colors_menu, "Residue Type", command=lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "residue_type"))
        add_qt_menu_command(multiple_residues_colors_menu, "Polarity", command=lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "polarity"))
        multiple_color_menu.addMenu(multiple_residues_colors_menu)


        # Secondary structure colors.
        n_selected_seqs = len(sequences_list)
        n_structures = len([e for e in sequences_list if e.has_structure()])
        n_seq_with_predicted_sec_str = len([e for e in sequences_list if e.has_predicted_secondary_structure()])

        if n_structures + n_seq_with_predicted_sec_str == n_selected_seqs:
            # multiple_color_menu.addSeparator()
            multiple_sec_str_color_menu = QtWidgets.QMenu('By Secondary Structure', multiple_color_menu)
            multiple_color_menu.addMenu(multiple_sec_str_color_menu)

            # Available when all the selected sequences have a 3D structure.
            if n_structures == n_selected_seqs:
                add_qt_menu_command(multiple_sec_str_color_menu, "Observed", lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "secondary-observed"))

            # Available only if all the sequences have a predicted secondary structure.
            if n_seq_with_predicted_sec_str == n_selected_seqs:
                add_qt_menu_command(multiple_sec_str_color_menu, "Predicted by PSIPRED", lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "secondary-predicted"))

            # Available if there is at least one element with a 3D structure or a secondary
            # structure prediction.
            if not n_structures == n_selected_seqs:
                add_qt_menu_command(multiple_sec_str_color_menu, "Auto (Observed + Predicted)", lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "secondary-auto"))


        # Conservation colors.
        sel_has_campo_scores = all([e.has_campo_scores() for e in sequences_list])
        sel_has_entropy_scores = all([e.has_entropy_scores() for e in sequences_list])
        if sel_has_campo_scores or sel_has_entropy_scores:
            multiple_cons_colors_menu = QtWidgets.QMenu('By Conservation', multiple_color_menu)
            multiple_color_menu.addMenu(multiple_cons_colors_menu)
            if sel_has_campo_scores:
                add_qt_menu_command(multiple_cons_colors_menu, "CAMPO scores", lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "campo-scores"))
            if sel_has_entropy_scores:
                add_qt_menu_command(multiple_cons_colors_menu, "Entropy scores", lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "entropy-scores"))


        # Energy colors.
        if all([e.has_dope_scores() for e in sequences_list]):
            # multiple_color_menu.addSeparator()
            multiple_energy_colors_menu = QtWidgets.QMenu('By Energy', multiple_color_menu)
            multiple_color_menu.addMenu(multiple_energy_colors_menu)
            add_qt_menu_command(multiple_energy_colors_menu, "DOPE scores", lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "dope"))


        # Color by domain.
        if all([e.has_domains() for e in sequences_list]):
            add_qt_menu_command(multiple_color_menu, "By Domain", lambda: self.parent_group.main_window.color_selection(color_selection_mode, color_selection_target, "domains"))

        return multiple_color_menu


    def build_selection_structure_menu(self):
        self.selection_structure_context_submenu = QtWidgets.QMenu('Structures', self.selection_context_submenu)
        self.selection_context_submenu.addMenu(self.selection_structure_context_submenu)

        if self.parent_group.pymod.all_sequences_have_structure():
            add_qt_menu_command(self.selection_structure_context_submenu, "Show chains in PyMOL", self.show_selected_chains_in_pymol_from_popup_menu)
            add_qt_menu_command(self.selection_structure_context_submenu, "Hide chains in PyMOL", self.hide_selected_chains_in_pymol_from_popup_menu)
            self.selection_structure_context_submenu.addSeparator()
            add_qt_menu_command(self.selection_structure_context_submenu, "Show Chains as Hedgehog", self.show_hedgehog_chains_in_pymol_from_popup_menu)
            add_qt_menu_command(self.selection_structure_context_submenu, "Show Heteroatoms", self.show_het_chains_in_pymol_from_popup_menu)
            add_qt_menu_command(self.selection_structure_context_submenu, "Hide Heteroatoms", self.hide_het_chains_in_pymol_from_popup_menu)

        elif self.parent_group.pymod.all_sequences_have_fetchable_pdbs():
            add_qt_menu_command(self.selection_structure_context_submenu, "Fetch PDB Files", lambda: self.parent_group.pymod.fetch_pdb_files("selection", None))


    def build_selection_cluster_menu(self):

        self.selection_cluster_context_submenu = QtWidgets.QMenu('Cluster Options', self.selection_context_submenu)
        self.selection_context_submenu.addMenu(self.selection_cluster_context_submenu)

        add_qt_menu_command(self.selection_cluster_context_submenu, "Extract Sequences from their Clusters", self.extract_selection_from_cluster)

        selected_sequences = self.parent_group.pymod.get_selected_sequences()
        mothers_set = set([s.mother for s in selected_sequences])
        if len(mothers_set) == 1:
            if len(selected_sequences) < len(self.parent_group.pymod_element.mother.get_children()):
                add_qt_menu_command(self.selection_cluster_context_submenu, "Extract Sequences to New Cluster", self.extract_selection_to_new_cluster_from_left_menu)


    ##########################
    # Cluster elements menu. #
    ##########################

    def build_cluster_popup_menu(self, target_menu, mode="cluster", extra_spacer=False):
        self.build_cluster_edit_menu(target_menu)


    def build_cluster_edit_menu(self, target_menu):

        self.cluster_edit_context_submenu = QtWidgets.QMenu('Edit Cluster', target_menu)
        target_menu.addMenu(self.cluster_edit_context_submenu)

        add_qt_menu_command(self.cluster_edit_context_submenu, "Save Alignment To File", self.save_alignment_from_left_pan)
        self.cluster_edit_context_submenu.addSeparator()
        add_qt_menu_command(self.cluster_edit_context_submenu, "Delete Gap Only Columns", self.delete_gap_only_columns_from_left_pane)
        add_qt_menu_command(self.cluster_edit_context_submenu, "Transfer Alignment", self.transfer_alignment_from_left_pane)
        self.cluster_edit_context_submenu.addSeparator()
        add_qt_menu_command(self.cluster_edit_context_submenu, "Delete Cluster", self.delete_alignment_from_left_pane)


    ###############################################################################################
    # Sequence manipulation events.                                                               #
    ###############################################################################################

    #------------
    # Clusters. -
    #------------

    def extract_from_cluster(self):
        """
        Extracts an element from an alignment.
        """
        self.parent_group.pymod_element.extract_to_upper_level()
        self.parent_group.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)

    def extract_selection_from_cluster(self):
        selected_sequences = self.parent_group.pymod.get_selected_sequences()
        # Using 'reversed' keeps them in their original order once extracted.
        for e in reversed(selected_sequences):
            e.extract_to_upper_level()
        self.parent_group.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)

    def extract_selection_to_new_cluster_from_left_menu(self):
        # 'gridder' is called in this method.
        self.parent_group.pymod.extract_selection_to_new_cluster()

    def make_lead_from_left_menu(self):
        self.parent_group.pymod_element.set_as_lead()

    def remove_lead_from_left_menu(self):
        self.parent_group.pymod_element.remove_all_lead_statuses()


    #------------------
    # Save sequences. -
    #------------------

    def save_sequence_from_left_pane(self):
        """
        Save option in the popup menu, it saves a single sequence.
        """
        self.parent_group.pymod.sequence_save_dialog(self.parent_group.pymod_element)

    def save_selection_from_left_pane(self):
        self.parent_group.pymod.save_selection_dialog()


    #---------------------
    # Copy to clipboard. -
    #---------------------

    def copy_sequence_to_clipboard(self):
        """
        Copy option in the popup menu, copies a single sequence.
        """
        cb = self.parent_group.main_window.clipboard
        cb.clear(mode=cb.Clipboard)
        cb.setText(self.parent_group.pymod_element.my_sequence, mode=cb.Clipboard)


    def copy_selection(self):
        cb = self.parent_group.main_window.clipboard
        cb.clear(mode=cb.Clipboard)
        text_to_copy = ""
        for element in self.parent_group.pymod.get_selected_sequences():
            # TODO: Adapt it for WINDOWS.
            text_to_copy += element.my_sequence + "\n"
        cb.setText(text_to_copy, mode=cb.Clipboard)


    #---------------------------------
    # Edit sequences and structures. -
    #---------------------------------

    def edit_sequence_from_context(self):
        self.parent_group.pymod.show_edit_sequence_window(self.parent_group.pymod_element)

    def search_string_from_context(self):
        self.parent_group.pymod.show_search_string_window(self.parent_group.pymod_element)


    #------------------------------------------------
    # Build new sequences and delete old sequences. -
    #------------------------------------------------

    def duplicate_sequence_from_the_left_pane(self):
        """
        Duplicates a single sequence.
        """
        self.parent_group.pymod.duplicate_sequence(self.parent_group.pymod_element)
        self.parent_group.main_window.gridder()

    def duplicate_selection(self):
        for e in self.parent_group.pymod.get_selected_sequences():
            self.parent_group.pymod.duplicate_sequence(e)
        self.parent_group.main_window.gridder()


    def delete_sequence_from_the_left_pane(self):
        """
        Delete option in the popup menu.
        """
        self.parent_group.pymod_element.delete()
        self.parent_group.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)

    def delete_many_sequences(self):
        # Delete the selected sequences.
        for element in self.parent_group.pymod.get_selected_sequences():
            element.delete()
        # Empty cluster elements will be deleted in the 'gridder' method.
        self.parent_group.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)


    #------------------------------
    # Save and delete alignments. -
    #------------------------------

    def _get_cluster_from_popup_menu(self, pymod_element):
        # If the element is a cluster, return it.
        if pymod_element.is_cluster():
            return pymod_element
        # If it's the lead of a collapsed cluster, return its mother (a cluster element).
        else:
            return pymod_element.mother

    def save_alignment_from_left_pan(self):
        self.parent_group.pymod.alignment_save_dialog(self._get_cluster_from_popup_menu(self.parent_group.pymod_element))

    def delete_alignment_from_left_pane(self):
        self.parent_group.pymod.delete_cluster_dialog(self._get_cluster_from_popup_menu(self.parent_group.pymod_element))

    def delete_gap_only_columns_from_left_pane(self):
        remove_gap_only_columns(self._get_cluster_from_popup_menu(self.parent_group.pymod_element))
        self.parent_group.pymod.main_window.gridder(update_clusters=True, update_elements=True)

    def transfer_alignment_from_left_pane(self):
        self.parent_group.pymod.transfer_alignment(self._get_cluster_from_popup_menu(self.parent_group.pymod_element))


    #----------------------------
    # Save PDB chains to files. -
    #----------------------------

    def save_structure_from_left_pane(self):
        self.parent_group.pymod.save_pdb_chain_to_file_dialog(self.parent_group.pymod_element)


    ###############################################################################################
    # Interact with the chain in PyMOL.                                                           #
    ###############################################################################################

    def center_chain_in_pymol_from_header_entry(self):
        self.parent_group.pymod.center_chain_in_pymol(self.parent_group.pymod_element)

    def hide_chain_in_pymol_from_header_entry(self):
        self.parent_group.pymod.hide_chain_in_pymol(self.parent_group.pymod_element)

    def show_chain_in_pymol_from_header_entry(self):
        self.parent_group.pymod.show_chain_in_pymol(self.parent_group.pymod_element)

    def show_hedgehog_in_pymol_from_header_entry(self):
        self.parent_group.pymod.show_hedgehog_in_pymol(self.parent_group.pymod_element)

    def show_het_in_pymol_from_header_entry(self):
        self.parent_group.pymod.show_het_in_pymol(self.parent_group.pymod_element)

    def hide_het_in_pymol_from_header_entry(self):
        self.parent_group.pymod.hide_het_in_pymol(self.parent_group.pymod_element)

    def show_selected_chains_in_pymol_from_popup_menu(self):
        for e in self.parent_group.pymod.get_selected_sequences():
            self.parent_group.pymod.show_chain_in_pymol(e)

    def hide_selected_chains_in_pymol_from_popup_menu(self):
        for e in self.parent_group.pymod.get_selected_sequences():
            self.parent_group.pymod.hide_chain_in_pymol(e)

    def show_hedgehog_chains_in_pymol_from_popup_menu(self):
        for e in self.parent_group.pymod.get_selected_sequences():
            self.parent_group.pymod.show_hedgehog_in_pymol(e)

    def show_het_chains_in_pymol_from_popup_menu(self):
        for e in self.parent_group.pymod.get_selected_sequences():
            self.parent_group.pymod.show_het_in_pymol(e)

    def hide_het_chains_in_pymol_from_popup_menu(self):
        for e in self.parent_group.pymod.get_selected_sequences():
            self.parent_group.pymod.hide_het_in_pymol(e)

    def click_structure_with_middle_button(self, event=None):
        if self.parent_group.pymod_element.has_structure():
            # Shows the structure and centers if the sequence is selected in Pymod.
            if self.parent_group.pymod_element.selected:
                self.show_chain_in_pymol_from_header_entry()
                self.center_chain_in_pymol_from_header_entry()
            # If the sequence is not selected in Pymod, hide it in PyMOL.
            else:
                self.hide_chain_in_pymol_from_header_entry()


class Header_context_menu(QtWidgets.QMenu):
    """
    Class for the context menu appearing when the user clicks with the right button on the header
    of an element.
    """

    def __init__(self, parent):
        QtWidgets.QMenu.__init__(self, parent)
        self.parent = parent

        # Sets the context menu style.
        bg_color = "gray" # gray, #ABABAB
        context_menu_style = "color: white; font: %spt %s; font-weight: bold; padding: 0px; background: %s" % (self.parent.get_main_window().font_size, self.parent.get_main_window().font, bg_color)
        self.setStyleSheet(context_menu_style)
