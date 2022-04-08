# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing the widgets showing the sequences of PyMod elements in the right pane of PyMod
main window.
"""

import time # For development only.

from pymol.Qt import QtWidgets, QtCore, QtGui

from ._header_entry_qt import Header_context_menu
from pymod_lib.pymod_gui.shared_gui_components_qt import add_qt_menu_command, highlight_color
from pymod_lib.pymod_vars import psipred_element_dict

from pymol import cmd


###################################################################################################
# Widget for the sequence of an element.                                                          #
###################################################################################################

class MyQLabel_sequence(QtWidgets.QLabel):
    """
    A custom QLabel for sequences.
    """

    def __init__(self, parent_group):

        QtWidgets.QLabel.__init__(self, " ")

        self.setTextFormat(QtCore.Qt.RichText) # Set the text format as 'RichText' in order to display html.
        self.parent_group = parent_group

        self.setMouseTracking(True)
        self.resize_to_content()

        self.current_pos = 0
        self.anchor_residue_number = 0
        self.number_of_new_gaps = 0
        self.gap_insert_velocity = 0.01

        # Gaps (which are placed outside html tags) will be colored in white.
        self.setStyleSheet("color: #ffffff")

        self.drag_left_performed = None
        self.drag_right_performed = None

        self.set_default_cursor()


    def get_main_window(self):
        return self.parent_group.main_window

    def set_default_cursor(self):
        # self.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))

    def set_dragging_cursor(self):
        # self.setCursor(QtGui.QCursor(QtCore.Qt.ClosedHandCursor))
        self.setCursor(QtGui.QCursor(QtCore.Qt.OpenHandCursor))


    def build_html_seq(self):
        """
        This method takes the sequence of a PyMod element (stored in its 'my_sequence' attribute)
        and formats it as html text to be displayed in a QLabel widget of Qt. A sequence like:

            ARVLPI--CW

        where the 'ARVL' residues are colored in white (#ffffff) and where the 'PICW' residues are
        colored in red (#f43030) will be formatted as:

            <font color='#ffffff'>ARVL</font><font color='#f43030'>PI</font>--<font color='#f43030'>CW</font>

        each group of contigous residues with the same color will be included in the same 'font' tag.
        Note that groups of residues with the same color separated by gaps, will be placed in different
        tags. In this way, gaps ("-" characters) will always be placed outside the 'font' tags and will
        be colored in white.
        """

        t1 = time.time()

        # Get the residues of the sequence.
        if not self.parent_group.pymod_element.is_cluster():
            residues = self.parent_group.pymod_element.get_polymer_residues()
        else:
            residues = None

        # Define the function to color a residue.
        if self.parent_group.pymod_element.color_by == "regular":
            regular_color = self.parent_group.pymod.all_colors_dict_tkinter[self.parent_group.pymod_element.my_color]
            def _get_color(residue_count, residues):
                return regular_color
        else:
            def _get_color(residue_count, residues):
                return residues[residue_count].color_seq

        seq_text = ""
        previous_char = None
        previous_color = None
        residue_count = 0

        for p in self.parent_group.pymod_element.my_sequence:

            # Inserts residues.
            if p != "-":

                color = _get_color(residue_count, residues)
                residue_count += 1

                # Opens the first tag.
                if previous_color == None:
                    seq_text += "<font color='%s'>%s" % (color, p)
                    # seq_text += "<font color='" + color + "'>" + p

                # Handle the rest of the tags.
                else:
                    # Handle a residue.
                    if previous_char != "-":
                        # Continues to insert text in the current tag.
                        if color == previous_color:
                            seq_text += p
                        # Closes the current tag.
                        else:
                            seq_text += "</font><font color='%s'>%s" % (color, p)
                            # seq_text += "</font><font color='" + color + "'>" + p

                    # Opens a new tag when ending a gap.
                    else:
                        seq_text += "<font color='%s'>%s" % (color, p)
                        # seq_text += "<font color='" + color + "'>" + p

                previous_color = color

            # Inserts gaps.
            else:
                if previous_char == None or previous_char == "-":
                    seq_text += "-"
                else:
                    seq_text += "</font>-"

            previous_char = p

        # Closes the last tag.
        if p != "-":
            seq_text += "</font>"

        if self.parent_group.pymod.DEVELOP:
            print("- To build an html seq, it took:", time.time() - t1)

        return seq_text


    def update_text(self):
        self.setText(self.build_html_seq())


    def resize_to_content(self):
        self.pixelsWide = self.get_main_window().fm.boundingRect(self.text() + "  ").width();
        self.pixelsHigh = self.get_main_window().fm.height()


    def get_html_content(self, html_text):
        """
        Gets an html code, removes tags and return only the tag content.
        """
        return "".join([i for i in html_text if i.isupper() or i == "-"])


    ###############################################################################################
    # Mouse events to get sequence information.                                                   #
    ###############################################################################################

    def get_mousecurrent_position_on_string(self, event, curr_pos=None):
        """
        When hovering with the mouse over a sequence, it returns the alignment position highlighted by the mouse.
        """

        curr_pos = event.pos().x() # - 3

        # Extracts the sequence from the formatted text.
        seq_only_text = self.get_html_content(self.text())

        # Length in pixels of the sequence text (cosidering only sequence characters and not the
        # html tags used for coloring).
        lun_pixels = self.get_main_window().fm.width(seq_only_text)
        # Length of the sequence itself (including gaps).
        lun_seq = len(seq_only_text)

        # Position (index + 1) of the hovered character in the sequence.
        if lun_pixels != 0:
            seq_index = int(curr_pos/float(lun_pixels)*lun_seq) + 1
        else:
            return "start"

        if seq_index > lun_seq:
            return "end"
        else:
            return seq_index


    def get_highlighted_residue(self, alignment_pos_id):
        """
        Gets the highlighted position in the aligned sequence.
        """
        return self.parent_group.pymod_element.get_residue_by_index(alignment_pos_id, aligned_sequence_index=True)


    def leaveEvent(self, event):
        """
        Called when the mouse leaves a sequence label.
        """
        self.parent_group.main_window.central_widget.textbox_sequence.setText("")
        self.parent_group.main_window.central_widget.textbox_position.setText("")


    #---------------------------------------
    # Interact with the residues in PyMOL. -
    #---------------------------------------

    def contextMenuEvent(self, event):
        """
        Shows the context menu when right clickings on the sequence of an element.
        """

        if self.parent_group.pymod_element.is_cluster():
            return None

        alignment_id = self.get_mousecurrent_position_on_string(event)
        if alignment_id in ("start", "end"):
            return None
        alignment_id = alignment_id - 1

        self.context_menu = Header_context_menu(self)
        add_qt_menu_command(self.context_menu, self.parent_group.pymod_element.my_header, None)
        self.context_menu.addSeparator()
        if self.parent_group.pymod_element.has_structure():
            add_qt_menu_command(self.context_menu, "Select Residue in PyMOL", command=lambda a=None, aid=alignment_id: self.select_residue_in_pymol_from_sequence_text(aid))
            add_qt_menu_command(self.context_menu, "Center Residue in PyMOL", command=lambda a=None, aid=alignment_id: self.center_residue_in_pymol_from_sequence_text(aid))
            self.context_menu.addSeparator()
        add_qt_menu_command(self.context_menu, "Add Residue Feature", command=lambda a=None, aid=alignment_id: self.add_feature_from_sequence_text(aid))
        action = self.context_menu.exec_(self.mapToGlobal(event.pos()))

    def click_residue_with_middle_button(self, event):
        alignment_id = self.get_mousecurrent_position_on_string(event)
        if alignment_id in ("start", "end"):
            return None
        alignment_id = alignment_id - 1
        if self.parent_group.pymod_element.has_structure():
            if self.parent_group.pymod_element.my_sequence[alignment_id] != "-":
                self.select_residue_in_pymol_from_sequence_text(alignment_id)
                self.center_residue_in_pymol_from_sequence_text(alignment_id)


    def select_residue_in_pymol_from_sequence_text(self, alignment_id):
        res = self.get_highlighted_residue(alignment_id)
        cmd.select("pymod_selection", res.get_pymol_selector())
        # cmd.indicate("pymod_selection")

    def center_residue_in_pymol_from_sequence_text(self, alignment_id, event=None):
        res = self.get_highlighted_residue(alignment_id)
        cmd.center(res.get_pymol_selector())

    def add_feature_from_sequence_text(self, alignment_id, event=None):
        res = self.get_highlighted_residue(alignment_id)
        self.parent_group.pymod.show_add_feature_window(self.parent_group.pymod_element, res)


    ###############################################################################################
    # Mouse events for sequence manipulation.                                                     #
    ###############################################################################################

    def mousePressEvent(self, event):
        """
        Click on an element sequence.
        """

        # Manipulate sequences.
        if event.buttons() == QtCore.Qt.LeftButton:

            self.set_dragging_cursor()
            res_pos = self.get_mousecurrent_position_on_string(event)

            # Refreshes the number of new gaps.
            self.number_of_new_gaps = 0

            # Sets the position of the 'anchor' residue (the residue where the original click happened
            # before the user started to drag).
            if res_pos in ("start", "end"):
                self.anchor_residue_number = 0
                self.initialize_drag = False
                return None

            # Allow to drag sequences only when clicking on residues.
            self.anchor_residue_number = res_pos
            if self.parent_group.pymod_element.my_sequence[res_pos-1] in ("-", ".", ":", "*") or self.parent_group.pymod_element.is_cluster():
                self.initialize_drag = False
            else:
                self.initialize_drag = True

        # Quickly center the residue in PyMOL.
        elif event.buttons() == QtCore.Qt.MiddleButton:
            self.click_residue_with_middle_button(event)


    domain_str_attr = "name" # "description"

    def mouseMoveEvent(self, event):

        #--------------------------------------------------------------------------------------
        # Print only sequence name and position in mw.central_widget if no button is pressed. -
        #--------------------------------------------------------------------------------------

        if event.buttons() == QtCore.Qt.NoButton:

            # Get the highlighted residue object.
            is_residue = False
            alignment_pos = self.get_mousecurrent_position_on_string(event)
            if alignment_pos in ("start", "end"): # type(alignment_pos) is int:
                residue_information = alignment_pos
            else:
                alignment_pos_text = "Alignment Position: %s" % alignment_pos
                if self.parent_group.pymod_element.is_cluster():
                    residue_information = alignment_pos_text
                else:
                    if self.parent_group.pymod_element.my_sequence[alignment_pos-1] != "-":
                        pymod_res = self.get_highlighted_residue(alignment_pos-1)
                        residue_information = "%s %s - %s" % (pymod_res.three_letter_code, pymod_res.db_index, alignment_pos_text)
                        is_residue = True
                    else:
                        residue_information = "Gap - %s" % alignment_pos_text

            # Get additional information (depending on the color scheme of the element) to show in
            # the message bar.
            if is_residue:
                # Show the CAMPO score.
                if self.parent_group.pymod_element.color_by == "campo-scores":
                    residue_information += " - CAMPO score: %s" % (pymod_res.campo_score["campo-score"])
                # Show the entropy score.
                elif self.parent_group.pymod_element.color_by == "entropy-scores":
                    residue_information += " - Entropy score: %s" % (pymod_res.entropy_score["entropy-score"])
                # Show the prediction confidence.
                elif self.parent_group.pymod_element.color_by == "secondary-predicted":
                    prediction = pymod_res.psipred_result
                    pred_text = "%s %s" % (prediction["confidence"], psipred_element_dict[prediction["sec-str-element"]])
                    residue_information += " - PSIPRED confidence: %s" % (pred_text)
                # Show the DOPE score for the residue.
                elif self.parent_group.pymod_element.color_by == "dope":
                    score = pymod_res.dope_score["score"]
                    residue_information += " - DOPE score: %s" % (score)
                # Show the name of the domains.
                elif self.parent_group.pymod_element.color_by == "domains":
                    res_domains = pymod_res.features["domains"]
                    if len(res_domains) > 1: # More than one domain is assigned to the residue.
                        domain_list_text = " - ".join(["(%s) %s" % (d_i+1, getattr(d, self.domain_str_attr)) for d_i, d in enumerate(res_domains)])
                        residue_information += " - Domains: %s" % (domain_list_text)
                    elif len(res_domains) == 1: # Only one domain is assigned to the residue.
                        residue_information += " - Domain: %s" % (getattr(res_domains[0], self.domain_str_attr))
                # Show the name of the sequence features.
                elif self.parent_group.pymod_element.color_by == "custom":
                    res_features = [f for f in pymod_res.features["sequence"] if f.show]
                    if len(res_features) > 1: # More than one domain is assigned to the residue.
                        features_list_text = " - ".join(["(%s) %s" % (d_i+1, getattr(d, self.domain_str_attr)) for d_i, d in enumerate(res_features)])
                        residue_information += " - Features: %s" % (features_list_text)
                    elif len(res_features) == 1: # Only one domain is assigned to the residue.
                        residue_information += " - Feature: %s" % (getattr(res_features[0], self.domain_str_attr))

            # Show information on the message bars of the main window.
            self.parent_group.main_window.central_widget.textbox_sequence.setText(str(self.parent_group.pymod_element.my_header))
            self.parent_group.main_window.central_widget.textbox_position.setText(str(residue_information))


        #----------------------------------
        # We start dragging the sequence. -
        #----------------------------------

        elif event.buttons() == QtCore.Qt.LeftButton and self.initialize_drag:

            self.my_pos = event.pos().x()
            residue_index_of_mouse_over = self.get_mousecurrent_position_on_string(event, self.my_pos)

            if residue_index_of_mouse_over != 'start':

                if self.anchor_residue_number > 0:

                    # Insert a gap behind any residue but the last.
                    if residue_index_of_mouse_over != 'end':

                        # Perform no action when dragging to the left of the first residue.
                        if residue_index_of_mouse_over < 1:
                            return None

                        # Dragging the sequence on the right introduces a gap.
                        if residue_index_of_mouse_over >= self.anchor_residue_number + self.gap_insert_velocity:

                            self.anchor_residue_number = residue_index_of_mouse_over
                            self.insert_gap_in_text()
                            self.number_of_new_gaps += 1
                            self.drag_right_performed = True


                        # Dragging the sequence on the left remove previous gap, if present.
                        if residue_index_of_mouse_over < self.anchor_residue_number:

                            if (self.parent_group.pymod_element.my_sequence[residue_index_of_mouse_over-1] == "-" or
                                self.parent_group.pymod_element.my_sequence[residue_index_of_mouse_over] == "-"):

                                if self.parent_group.pymod_element.my_sequence[residue_index_of_mouse_over-1] == "-":
                                    self.delete_gap_in_text(residue_index_of_mouse_over, gap_offset=-1)
                                else:
                                    self.delete_gap_in_text(residue_index_of_mouse_over, gap_offset=0)

                            self.number_of_new_gaps += self.number_of_new_gaps - 1
                            self.anchor_residue_number = residue_index_of_mouse_over
                            self.drag_left_performed = True

                    # Insert a gap behind the residue.
                    else:
                        self.anchor_residue_number += 1
                        self.insert_gap_in_text()
                        self.number_of_new_gaps += 1
                        self.drag_right_performed = True

                    # If the sequence is a child, the alignment symbols of the mother have to be
                    # adjusted.
                    if (self.drag_left_performed or self.drag_right_performed) and self.parent_group.pymod_element.is_child():
                        # Updates the mother's sequence.
                        mother = self.parent_group.pymod_element.mother
                        mother.update_stars(adjust_elements=False)
                        # Update the mother's sequence widget.
                        mother.widget_group.sequence_text.update_text()


    def insert_gap(self, string, index=0):
        return string[:index] + '-' + string[index:]

    def delete_gap(self, string, index=0):
        return string[:index]+ string[index+1:]


    def insert_gap_in_text(self):

        #-----------------------------------
        # Inserts the gap in the sequence. -
        #-----------------------------------

        old_seq = self.parent_group.pymod_element.my_sequence
        new_seq = self.insert_gap(old_seq, self.anchor_residue_number - 2)
        self.parent_group.pymod_element.my_sequence = new_seq

        #---------------------------------------------------
        # Inserts the gap in the html text of the element. -
        #---------------------------------------------------

        old_html_seq = self.text()
        new_html_seq = self.build_html_seq()
        self.setText(new_html_seq)


    def delete_gap_in_text(self, residue_index_of_mouse_over, gap_offset=0):

        old_seq = self.parent_group.pymod_element.my_sequence
        new_seq = self.delete_gap(old_seq, residue_index_of_mouse_over + gap_offset)
        self.parent_group.pymod_element.my_sequence = new_seq

        old_html_seq = self.text()
        new_html_seq = self.build_html_seq()

        self.setText(new_html_seq)


    def mouseReleaseEvent(self, event):

        self.set_default_cursor()

        # If the element is in a cluster, modifies the sequence text of other elements of the
        # cluster.
        if self.parent_group.pymod_element.is_child() and (self.drag_left_performed or self.drag_right_performed):

            #######################################################################################
            # NOTE: an optimal way to do this would be to rstrip all the sequences, then to ljust #
            # them to the lenght of the "longest" one. However Tkinter is too slow to do this, it #
            # takes too much time to update all the sequences in big clusters at the same time,   #
            # therefore as long as Tkinter is used, the following code has to be applied. This    #
            # code prevents every sequence of a cluster from being updated every time an indel is #
            # added, and it tries to update only the necessary sequences.                         #
            #######################################################################################

            elements_to_update = [self.parent_group.pymod_element] + self.parent_group.pymod_element.get_siblings()

            max_len = max([len(e.my_sequence.rstrip("-")) for e in elements_to_update])
            for element in elements_to_update + [self.parent_group.pymod_element.mother]:
                element.widget_group.sequence_text.adjust_to_length(adjusted_len=max_len)

        self.drag_left_performed = None
        self.drag_right_performed = None


    def adjust_to_length(self, adjusted_len):
        if len(self.parent_group.pymod_element.my_sequence) > adjusted_len: # len(self.parent_group.pymod_element.my_sequence):
            self.rstrip_entry(adjusted_len)
        elif len(self.parent_group.pymod_element.my_sequence) < adjusted_len: # len(self.parent_group.pymod_element.my_sequence):
            self.ljust_entry(adjusted_len)

    def rstrip_entry(self, maxlength, update=True):
        self.parent_group.pymod_element.my_sequence = self.parent_group.pymod_element.my_sequence.rstrip("-").ljust(maxlength, "-")
        self.update_text()

    def ljust_entry(self, maxlength, update=True):
        self.parent_group.pymod_element.my_sequence = self.parent_group.pymod_element.my_sequence.ljust(maxlength, "-")
        self.update_text()


###################################################################################################
# Widget for the cluster button of an element.                                                    #
###################################################################################################

class Cluster_button(QtWidgets.QLabel):
    """
    A custom class for cluster buttons in the main window of PyMod.
    """

    fg_color = "white"
    bd_color = "gray"

    def __init__(self, parent_group):

        QtWidgets.QLabel.__init__(self, "-")
        self.parent_group = parent_group

        self.setMouseTracking(True)
        self.set_default_cursor()


    def set_default_cursor(self):
        # self.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))

    def mousePressEvent(self, event):
        if self.parent_group.pymod_element.is_cluster():
            if event.buttons() == QtCore.Qt.LeftButton:
                self.parent_group.cluster_button_click(None)


    def set_visibility(self):
        if self.parent_group.pymod_element.is_cluster():
            self.set_cluster_button_style()
            self.setVisible(True)
        elif self.parent_group.pymod_element.is_child():
            self.set_child_sign_style()
            self.setVisible(True)


    def set_cluster_button_style(self):
        self.setStyleSheet("color: %s; border: 0px solid %s; background-color: %s; padding: 0px 1px" % (self.fg_color, self.bd_color, highlight_color))
        if self.parent_group._cluster_button_state:
            self.setText("-")
        else:
            self.setText("+")
        self.setToolTip('Press to Expand/Shrink')


    def set_child_sign_style(self):
        """
        Shows an additional entry inside the right-frame for child elements.
        """
        self.setStyleSheet("color: %s; border-left: 1px solid %s; padding: 0px 0px" % (self.fg_color, self.bd_color))

        if self.parent_group.pymod_element.is_blast_query():
            child_sign = "q" # "|q"
        elif self.parent_group.pymod_element.is_lead():
            child_sign = "l" # "|l"
        elif self.parent_group.pymod_element.is_bridge():
            child_sign = "b" # "|b"
        else:
            child_sign = "_" # "|_"
        self.setText(child_sign)
        self.setToolTip(None)
