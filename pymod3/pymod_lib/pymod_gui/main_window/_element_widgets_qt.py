# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing the "widget group" class for a PyMod element. A "widget group" is a class which
allows to control the "header entry" (a PyQT widget for the header of a PyMod element appearning in
the left pane of the PyMod main window) and "sequence entry" (a PyQT widget for the sequence of a
PyMod element appearning in the right pane of the PyMod main window).
The use of a "widget group" class, allows to simultaneously control the objects for the "header entry"
and "sequence entry".
"""

from pymol.Qt import QtCore

from ._header_entry_qt import MyQLabel_header
from ._sequence_text_qt import MyQLabel_sequence, Cluster_button


class PyMod_element_widgets_group_qt():
    """
    Class for coordinating the widgets belonging to a PyMod element. Represent a group of widgets
    belonging to a PyMod element.
    """

    unselected_color_dict = {True: "white", False: "red"}

    def __init__(self, main_window, pymod_element):

        self.pymod_element = pymod_element
        self.main_window = main_window
        self.pymod = self.main_window.pymod

        self.old_grid_row_index = None
        self.grid_row_index = None

        self.show = True
        self._collapsed_cluster = False


        #----------------------------
        # Builds the header widget. -
        #----------------------------

        self.header_entry = MyQLabel_header(parent_group=self)
        # By default widget are set as not visible. Visibility will be added in the 'gridder' method
        # of the PyMod main window.
        self.header_entry.setVisible(False)

        self.main_window.central_widget.id_form_layout.addRow(self.header_entry)


        #-----------------------------------
        # Builds the sequence text widget. -
        #-----------------------------------

        self.sequence_text = MyQLabel_sequence(parent_group=self)
        self.sequence_text.setVisible(False)
        self._show_sequence_text = True # Actually not used.


        #---------------------------------
        # Button and signs for clusters. -
        #---------------------------------

        self._cluster_button_state = True

        self.cluster_button = Cluster_button(parent_group=self)
        self.cluster_button.setVisible(False)

        self.main_window.central_widget.seq_form_layout.addRow(self.cluster_button, self.sequence_text)


    ###############################################################################################
    # Display widgets.                                                                            #
    ###############################################################################################

    #-----------------------
    # Control all widgets. -
    #-----------------------

    def grid_widgets(self, update_element_text=False):
        """
        Called only by the 'gridder' method.
        """

        # Shows the left pane widgets.
        self.grid_header(update_element_header=update_element_text)

        # Updates and shows the right pane widgets: adds the sequence of the element and its cluster
        # button/sign if present.
        self.grid_sequence(update_element_text=update_element_text)


    def hide_widgets(self, save_status=True):
        """
        Called by several methods.
        """
        if save_status:
            self.show = False
        self.hide_header()
        self.hide_sequence(save_status=save_status)


    #----------
    # Header. -
    #----------

    def grid_header(self, update_element_header=False):

        if update_element_header:
            self.header_entry.update_title()

        if self.old_grid_row_index != self.grid_row_index:

            i1 = self.main_window.central_widget.id_form_layout.itemAt(self.grid_row_index, 2)

            if i1 != None:
                if not i1.widget() is self.header_entry:
                    self.main_window.central_widget.id_form_layout.removeItem(i1)
                    i1.widget().setVisible(False)
                    self.main_window.central_widget.id_form_layout.removeWidget(self.header_entry)

            self.main_window.central_widget.id_form_layout.setWidget(self.grid_row_index, 2, self.header_entry)

        if not self.header_entry.isVisible():
            self.header_entry.setVisible(True)


    def hide_header(self):
        if self.header_entry.isVisible():
            self.header_entry.setVisible(False)


    #---------------------------------
    # Cluster button and child sign. -
    #---------------------------------

    def change_cluster_button_on_expand(self):
        self.cluster_button.setText('-')
        self._cluster_button_state = True
        if self.pymod_element.is_cluster():
            self._collapsed_cluster = False

    def change_cluster_button_on_collapse(self):
        self.cluster_button.setText('+')
        self._cluster_button_state = False
        if self.pymod_element.is_cluster():
            self._collapsed_cluster = True


    #-----------------
    # Sequence text. -
    #-----------------

    def grid_sequence(self, update_element_text=False, save_status=False):

        if save_status:
            self._show_sequence_text = True

        if update_element_text:
            self.sequence_text.update_text()

        # Change the widget position only if its 'grid_row_index' has changed.
        if self.old_grid_row_index != self.grid_row_index:

            # Get the sequence widget present in the row in which the current widget has to be placed.
            i1 = self.main_window.central_widget.seq_form_layout.itemAt(self.grid_row_index, 1)

            # If there is a widget in the target positio, it will be removed to make place for the
            # current widget.
            if i1 != None:
                # if not i1.widget() is self.sequence_text:

                    # Remove the target sequence widget.
                    self.main_window.central_widget.seq_form_layout.removeItem(i1)
                    i1.widget().setVisible(False)
                    # Remove the current sequence widget.
                    self.main_window.central_widget.seq_form_layout.removeWidget(self.sequence_text)

                    # Cluster buttons.
                    i2 = self.main_window.central_widget.seq_form_layout.itemAt(self.grid_row_index, 0)
                    # Remove the target cluster button widget.
                    self.main_window.central_widget.seq_form_layout.removeItem(i2)
                    i2.widget().setVisible(False)
                    # Remove the current cluster button widget.
                    self.main_window.central_widget.seq_form_layout.removeWidget(self.cluster_button)

            # Sets the current widgets.
            self.main_window.central_widget.seq_form_layout.setWidget(self.grid_row_index, 1, self.sequence_text)
            self.main_window.central_widget.seq_form_layout.setWidget(self.grid_row_index, 0, self.cluster_button)

        # Show the current widgets.
        if not self.sequence_text.isVisible():
            if not self.pymod_element.is_cluster():
                self.sequence_text.setVisible(True)
            else:
                self.sequence_text.setVisible(not self._collapsed_cluster)
        else:
            if self.pymod_element.is_cluster() and self._collapsed_cluster:
                self.sequence_text.setVisible(False)

        if not self.cluster_button.isVisible():
            self.cluster_button.set_visibility()


    def hide_sequence(self, save_status=True):

        if save_status:
            self._show_sequence_text = False

        if self.sequence_text.isVisible():
            self.sequence_text.setVisible(False)
        if self.cluster_button.isVisible():
            self.cluster_button.setVisible(False)


    ###############################################################################################
    # Selection of elements in the PyMod main window.                                             #
    ###############################################################################################

    def toggle_element(self):
        """
        Toggles elements selection state.
        """
        if self.pymod_element.selected: # Inactivate.
            self.deselect_element()
        else: # Activate.
            self.select_element()


    def deselect_element(self, deselect_all=False):
        """
        Deselects an element. The 'deselect_all' should be set to 'True' only when deselecting all
        elements from PyMod main menu.
        """
        if not deselect_all:
            self._deselect_recursively()
            if self.pymod_element.is_child():
                self._deselect_ancestry_recursively(is_in_cluster=True)
                self._color_headers_on_toggle()
        else:
            self._turn_selection_off()

    def select_element(self, select_all=False):
        """
        Selects an element.
        """
        if not select_all:
            self._select_recursively()
            if self.pymod_element.is_child():
                self._select_ancestry_recursively(is_in_cluster=True)
                self._color_headers_on_toggle()
        else:
            self._turn_selection_on()

    def select_collapsed_cluster_descendants(self):
        for descendant in self.pymod_element.get_descendants() + [self.pymod_element]:
            descendant.widget_group.select_element(select_all=True)

    def _deselect_recursively(self, is_in_cluster=False):
        """
        Deselect an element and all its children recursively.
        """
        self._turn_selection_off(is_in_cluster)
        if self.pymod_element.is_mother():
            for c in self.pymod_element.get_children():
                c.widget_group._deselect_recursively(is_in_cluster)

    def _select_recursively(self, is_in_cluster=False):
        """
        Select an element and all its children recursively.
        """
        self._turn_selection_on(is_in_cluster)
        if self.pymod_element.is_mother():
            for c in self.pymod_element.get_children():
                c.widget_group._select_recursively(is_in_cluster)


    def _deselect_ancestry_recursively(self, is_in_cluster=False):
        """
        Deselects the ancestry an element (that is, its mother and its mother's mother, and so on
        recursively).
        """
        if not self.pymod_element.is_child():
            return None
        mother = self.pymod_element.mother
        # Modify the mother and the siblings according to what happens to the children.
        if mother.selected:
            mother.widget_group._turn_selection_off(is_in_cluster=True)
        if mother.is_child():
            mother.widget_group._deselect_ancestry_recursively(is_in_cluster=True)

    def _select_ancestry_recursively(self, is_in_cluster=False):
        """
        Selects the ancestry of an element recursively.
        """
        # If the mother is not selected and if by selecting this child, all the children
        # are selected, also selects the mother.
        if not self.pymod_element.is_child():
            return None
        child = self.pymod_element
        mother = self.pymod_element.mother
        if not mother.selected:
            # If it is the last inactivated children in the cluster, by selecting it, all the
            # elements in the cluster will be selected and the mother will also be selected.
            siblings = child.get_siblings()
            if not False in [c.selected for c in siblings]:
                mother.widget_group._turn_selection_on()
        if mother.is_child():
            mother.widget_group._select_ancestry_recursively(is_in_cluster=False)


    def _set_header_color(self, color):
        self.header_entry.setStyleSheet("color: %s" % color)


    def _turn_selection_on(self, is_in_cluster=False):
        """
        Selects an element.
        """
        self.pymod_element.selected = True
        self._set_header_color("green")


    def _turn_selection_off(self, is_in_cluster=False):
        """
        Deselects an element.
        """
        self.pymod_element.selected = False
        self._set_header_color("red")


    def _color_headers_on_toggle(self):
          """
          Adjust the color of unselected headers in a cluster.
          """
          is_in_cluster = False
          ancestor = self.pymod_element.get_ancestor()
          if ancestor and not ancestor.selected:
              descendants = ancestor.get_descendants()
              for d in descendants:
                  if d.selected:
                      is_in_cluster = True
                      break

              # Descendants.
              for d in descendants:
                  if not d.selected:
                      d.widget_group._set_header_color(self.unselected_color_dict[is_in_cluster])

              # Ancestor.
              ancestor.widget_group._set_header_color(self.unselected_color_dict[is_in_cluster])


    ###############################################################################################
    # Control clusters.                                                                           #
    ###############################################################################################

    #-------------------------
    # Cluster button events. -
    #-------------------------

    def cluster_button_click(self, event):
        """
        Creates the mouse event for clicking cluster buttons. It is used to toggle the children of
        the cluster.
        """
        if self._cluster_button_state:
            self.collapse_cluster()
        else:
            self.expand_cluster()

    def _toggle_cluster_click(self, cluster_lead_action, cluster_element_action):

        pymod_cluster = self.pymod_element
        cluster_lead = pymod_cluster.get_lead()
        # If the cluster has a cluster lead.
        if cluster_lead:
            cluster_lead_action(cluster_lead)
        cluster_element_action(pymod_cluster)


    #--------------------------------
    # Expand and collapse clusters. -
    #--------------------------------

    # Expand.
    def expand_cluster(self):
        self._toggle_cluster_click(self._expand_cluster_lead, self._expand_cluster_element)


    def _expand_cluster_element(self, pymod_cluster):
        pymod_cluster_widgets_group = pymod_cluster.widget_group
        pymod_cluster_widgets_group.change_cluster_button_on_expand()
        # Shows the text of the collapsed cluster.
        pymod_cluster_widgets_group.grid_sequence(update_element_text=True, save_status=True)
        # Show the children of the collapsed cluster.
        for child in pymod_cluster.get_children():
            self._show_descendants(child)
        self.main_window.gridder(clear_selection=False, update_clusters=True, update_menus=False, update_elements=False)

    def _expand_cluster_lead(self, cluster_lead):
        pass


    def _show_descendants(self, pymod_element):
        pymod_element_widgets_group = pymod_element.widget_group
        if pymod_element.is_cluster():
            # If the element is not a collapsed cluster, then show it and all its children.
            if not pymod_element_widgets_group._collapsed_cluster:
                pymod_element_widgets_group.show = True
                for child in pymod_element.get_children():
                    self._show_descendants(child)
            # If the element is a collapsed cluster.
            else:
                pass
        else:
            pymod_element_widgets_group.show = True


    # Collapse.
    def collapse_cluster(self):
        self._toggle_cluster_click(self._collapse_cluster_lead, self._collapse_cluster_element)

    def _collapse_cluster_element(self, pymod_cluster):
        pymod_cluster_widgets_group = pymod_cluster.widget_group
        pymod_cluster_widgets_group.change_cluster_button_on_collapse()

        # Hide the cluster element sequence.
        pymod_cluster_widgets_group.hide_sequence()

        # Hide all the descendants widgets.
        for child in pymod_cluster.get_descendants():
            child.widget_group.hide_widgets()

        self.main_window.gridder(clear_selection=False, update_clusters=True, update_menus=False, update_elements=False)

    def _collapse_cluster_lead(self, cluster_lead):
        pass


    ###############################################################################################
    # Check coloring schemes.                                                                     #
    ###############################################################################################

    def can_be_colored_by_secondary_structure(self):
        """
        Returns 'True' if the element has an associated structure or has a secondary structure
        prediction.
        """
        return self.pymod_element.has_structure() or self.pymod_element.has_predicted_secondary_structure()
