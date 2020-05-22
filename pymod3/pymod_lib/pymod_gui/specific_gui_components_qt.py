# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Classes for PyQt widgets used in specific parts of the PyMod GUI.
"""

import os

from pymol.Qt import QtWidgets, QtCore

from pymod_lib.pymod_seq.seq_manipulation import clean_white_spaces_from_input
from pymod_lib.pymod_gui.shared_gui_components_qt import open_color_dialog

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_tool_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_entryfield_button_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_plaintextedit_qt)
from pymod_lib.pymod_gui.shared_gui_components_qt import (options_title_style, small_font_style,
                                                          active_entry_style, inactive_entry_style,
                                                          inactive_bg_color, success_bg_color, failure_bg_color,
                                                          askyesno_qt, askopenfile_qt, askdirectory_qt)


#####################################################################
# PyMod options window.                                             #
#####################################################################

class PyMod_options_window_qt(PyMod_tool_window_qt):
    """
    Window showing a series of options for PyMod.
    """

    middle_layout_type = "qgrid"

    def __init__(self, parent, pymod, *args, **configs):
        self.pymod = pymod
        PyMod_tool_window_qt.__init__(self, parent, *args, **configs)

    def add_middle_frame_widgets(self):

        self.tools_params_dict = {}
        self.row_counter = 0

        # This list will be populated inside "build_tool_options_frame()".
        for single_tool in self.pymod.pymod_tools:
            # If the tool list of parameter widgets has some alignable widgets, adds them to the
            # option window list.
            self.display_options(single_tool)

        # self.middle_formlayout.set_input_widgets_width(200)

        return None


    def display_options(self, single_tool):
        """
        Displays at list of option in the the target_frame contained in a target widget.
        Used in the PyMod options window.
        """

        # Check that at least one parameter
        if not any([p.show_widget for p in single_tool.parameters]):
            return None

        self.tools_params_dict[single_tool.name] = {}

        # Grids a label with the name of the tool.
        tool_full_name_label = QtWidgets.QLabel(single_tool.full_name)
        tool_full_name_label.setStyleSheet(options_title_style)

        # print(dir(tool_full_name_label))
        # print(tool_full_name_label.font().size())

        self.middle_formlayout.addWidget(tool_full_name_label, self.row_counter, 0)
        self.row_counter += 1

        # Actually grids the parmater widgets.
        for parameter in single_tool.parameters:
            if not parameter.show_widget:
                continue
            # If the display options return a widget, adds it ot the Tool list of parameter widgets.
            w = self.display_paramenter_options(parameter)
            self.row_counter += 1


    def display_paramenter_options(self, parameter):
        """
        Used to display a series of widgets to choose the parameter value in PyMod
        options window. This is conditioned on the type of input widget of the
        'parameter'.
        """

        # Label with the name of the parameter.
        param_full_name_label = QtWidgets.QLabel(parameter.full_name)
        self.middle_formlayout.addWidget(param_full_name_label, self.row_counter, 0)

        # Do not display any input widget.
        if parameter.widget_type is None:
            pass

        # Display any a text entry input.
        elif parameter.widget_type == "path_entryfield":

            # Entry for the path.
            path_entryfield = QtWidgets.QLineEdit(str(parameter.get_starting_value()))

            if parameter.editable:
                path_entryfield.setStyleSheet(active_entry_style + "; " + small_font_style)
            else:
                path_entryfield.setStyleSheet(inactive_entry_style + "; " + small_font_style)
                path_entryfield.setEnabled(False)

            self.middle_formlayout.addWidget(path_entryfield, self.row_counter, 1)

            self.tools_params_dict[parameter.parent_tool.name][parameter.name] = path_entryfield


            # Button for browsing the path.
            if parameter.editable:
                path_browse_button = QtWidgets.QPushButton("Browse")
                path_browse_button.setStyleSheet(small_font_style)

                if parameter.path_type in ("file", "directory"):
                    path_browse_button.clicked.connect(lambda a=None, p=parameter: self.choose_path(p))
                else:
                    raise KeyError("Unknown 'path_type': %s" % str(parameter.path_type))

                self.middle_formlayout.addWidget(path_browse_button, self.row_counter, 2)


            # Button to automatically identifying the path.
            if parameter.auto_find:
                auto_find_button = QtWidgets.QPushButton("Auto Find")
                auto_find_button.setStyleSheet(small_font_style)
                self.middle_formlayout.addWidget(auto_find_button, self.row_counter, 3)
                auto_find_button.clicked.connect(lambda a=None, e=path_entryfield: parameter.auto_find_command(e))


        # Show the status of a parameter.
        elif parameter.widget_type == "show_status":

            text_to_show, status = parameter.get_status()
            status_entryfield = QtWidgets.QLineEdit(text_to_show)
            if status:
                status_entryfield.setStyleSheet("background-color: %s" % success_bg_color)
            else:
                status_entryfield.setStyleSheet("background-color: %s" % failure_bg_color)
            status_entryfield.setEnabled(False)
            status_entryfield.setFixedWidth(200)

            self.middle_formlayout.addWidget(status_entryfield, self.row_counter, 1)


        else:
            raise KeyError("Unkown 'widget_type': %s" % parameter.widget_type)


    def choose_path(self, parameter):
        """
        Called when users press the 'Browse' button in order to choose a path on their system.
        """

        new_path = None
        entry = self.tools_params_dict[parameter.parent_tool.name][parameter.name]
        current_path = entry.text()
        askpath_title = "Search for %s %s" % (parameter.parent_tool.full_name, parameter.full_name)

        # Let users choose a new path.
        if parameter.path_type == "file":
            new_path = askopenfile_qt(askpath_title, parent=self.pymod.get_qt_parent(),
                                  initialdir=os.path.dirname(current_path),
                                  initialfile=os.path.basename(current_path))

        elif parameter.path_type == "directory":
            new_path = askdirectory_qt(askpath_title, parent=self.pymod.get_qt_parent(),
                                       initialdir=current_path)

        # Updates the text in the Entry with the new path name.
        if new_path:
            entry.clear()
            entry.setText(new_path)


    def get_value_from_gui(self, parameter_obj):
        """
        Gets the value from the input widgets in the PyMod options window. This
        should be conditioned on the type class of the 'Tool_parameter' in order
        to be able to retrieve different types of input from the GUI.
        """
        return self.tools_params_dict[parameter_obj.parent_tool.name][parameter_obj.name].text()


#####################################################################
# Window for new sequences.                                         #
#####################################################################

class Raw_sequence_window_qt(PyMod_tool_window_qt):
    """
    Window with two text entries to add the sequence elements to PyMod. The two
    entries are one for the name and one for the sequence of the element.
    """

    build_name_entry = True
    build_sequence_entry = True

    def add_middle_frame_widgets(self):

        entry_style = "background-color: white; color: black; font-family: courier"

        # Creates an Entry for the name of the new sequence.
        if self.build_name_entry:
            self.seq_name_input = PyMod_entryfield_qt(label_text="Name:", value="", style=entry_style)
            self.middle_formlayout.add_widget_to_align(self.seq_name_input)

        # Creates an Entry widget for the sequence.
        if self.build_sequence_entry:
            self.seq_sequence_input = PyMod_plaintextedit_qt(label_text="Sequence:", value="", style=entry_style)
            self.middle_formlayout.add_widget_to_align(self.seq_sequence_input)


    def get_sequence(self):
        return clean_white_spaces_from_input(self.seq_sequence_input.getvalue()).upper()

    def get_sequence_name(self):
        return self.seq_name_input.getvalue()


class Edit_sequence_window_qt(Raw_sequence_window_qt):
    """
    Window editing the sequence of an element already loaded in PyMod. Does not
    allow to edit its name.
    """

    build_name_entry = False
    build_sequence_entry = True

    def __init__(self, parent, pymod_element, *args, **configs):
        Raw_sequence_window_qt.__init__(self, parent, *args, **configs)
        self.pymod_element = pymod_element
        self.seq_sequence_input.setvalue(self.pymod_element.my_sequence)


class Import_from_pymol_window_qt(PyMod_tool_window_qt):
    """
    Window showing a series of checkboxes to import PyMOL objects in PyMod.
    """

    def __init__(self, parent, selections_list, *args, **configs):

        self.selections_list = selections_list
        PyMod_tool_window_qt.__init__(self, parent, *args, **configs)

        # Builds a combobox for each PyMOL object to import.
        self.sele_checkbox_list = []
        for sele in selections_list:
            checkbox = QtWidgets.QCheckBox(sele)
            self.sele_checkbox_list.append(checkbox)
            self.middle_formlayout.addRow(checkbox)

    def get_objects_to_import(self):
        sele_list = []
        for sele, checkbox in zip(self.selections_list, self.sele_checkbox_list):
            if checkbox.isChecked():
                sele_list.append(sele)
        return sele_list


###############################################
# Window for selecting the "PyMod directory". #
###############################################

class Dir_selection_dialog_mixin:
    """
    Mixin class to be incorporated in all the directory selection dialogs.
    """

    def keyPressEvent(self, event):
        """
        By overriding this method, the dialog will not close when pressing the "esc" key.
        """
        if event.key() == QtCore.Qt.Key_Escape:
            pass
        else:
            QtWidgets.QDialog.keyPressEvent(self, event)

    def closeEvent(self, evnt):
        if evnt.spontaneous():
            title = "Exit PyMod?"
            message = "Are you really sure you want to exit PyMod?"
            answer = askyesno_qt(title, message, parent=self.pymod.get_qt_parent())
            if answer:
                self.close() # Closes the dialog.
                self.main_window.close() # Close the main window if the user exits the dialog.
            else:
                evnt.ignore()


class PyMod_dir_selection_dialog(Dir_selection_dialog_mixin, QtWidgets.QDialog):
    """
    Dialog to select the PyMod Directory. This is shown when launching PyMod for
    the first time.
    """

    is_pymod_window = True

    def __init__(self, app, pymod, confirm_close=True):

        QtWidgets.QDialog.__init__(self, parent=app)

        self.main_window = app
        self.pymod = pymod
        self.confirm_close = confirm_close

        self.initUI()


    def initUI(self):

        self.setWindowTitle('PyMod Directory Selection')

        self.vertical_layout = QtWidgets.QVBoxLayout()

        # Main label.
        self.label = QtWidgets.QLabel("Select a folder inside which to build the 'PyMod Directory'", self)
        self.vertical_layout.addWidget(self.label)


        # Entry and "Browse" button.
        self.horizontal_layout = QtWidgets.QHBoxLayout()

        self.main_entry = QtWidgets.QLineEdit(self.pymod.home_directory, self)
        self.main_entry.setStyleSheet("background-color: white; color: black")
        self.horizontal_layout.addWidget(self.main_entry)

        self.browse_button = QtWidgets.QPushButton("BROWSE", self)
        self.browse_button.clicked.connect(self.pymod_directory_browse_state)
        self.horizontal_layout.addWidget(self.browse_button)

        self.vertical_layout.addLayout(self.horizontal_layout)


        # "Submit" button.
        self.submit_button = QtWidgets.QPushButton("SUBMIT", self)
        self.submit_button.setFixedWidth(self.submit_button.sizeHint().width())
        self.submit_button.clicked.connect(self.on_submit_button_press)
        self.vertical_layout.addWidget(self.submit_button)
        self.vertical_layout.setAlignment(self.submit_button, QtCore.Qt.AlignCenter)


        # Set the layouts.
        self.setLayout(self.vertical_layout)


    def on_submit_button_press(self):
        self.pymod.pymod_directory_selection_state()


    def pymod_directory_browse_state(self):
        """
        Let users choose a new path.
        """
        new_path = askdirectory_qt(title="Select a folder in which to build the 'PyMod Directory'",
                                   initialdir=str(self.main_entry.text()),
                                   parent=self.pymod.get_qt_parent())
        if new_path: # Updates the text in the Entry with the new path.
            self.main_entry.setText(new_path)


##############################################
# Window for starting a new "PyMod session". #
##############################################

# TODO.


##############################################
# Window for adding a feature to a sequence. #
##############################################

class Add_feature_window_qt(PyMod_tool_window_qt):

    def __init__(self, parent, pymod_element, selected_residue, *args, **configs):

        self.pymod_element = pymod_element
        self.selected_residue = selected_residue
        PyMod_tool_window_qt.__init__(self, parent, *args, **configs)


    def add_middle_frame_widgets(self):

        # Entryfield for selecting the residues range.
        self.residue_range_enf = PyMod_entryfield_qt(label_text="Residue(s)",
                                                     value=str(self.selected_residue.db_index))
        self.middle_formlayout.add_widget_to_align(self.residue_range_enf)


        # Entryfield for the feature name.
        self.feature_name_enf = PyMod_entryfield_qt(label_text="Feature Name",
                                                    value="new feature")
        self.middle_formlayout.add_widget_to_align(self.feature_name_enf)


        # Widgets for choosing a color for the feature.
        self.selected_rgb = (1.0, 0.0, 0.0)
        self.selected_hex = '#ff0000'
        self.feature_color_enf = PyMod_entryfield_button_qt(label_text="Feature Color",
                                                            readonly=True,
                                                            button_text="Pick",
                                                            button_command=self.pick_color_dialog)
        self.middle_formlayout.add_widget_to_align(self.feature_color_enf)
        self.feature_color_enf.entry.setStyleSheet("background-color: %s" % self.selected_hex)


        # Select in PyMOL.
        if self.pymod_element.has_structure():
            self.select_in_pymol_rds = PyMod_radioselect_qt(label_text="Select in PyMOL",
                                                      buttons=("Yes", "No"))
            self.select_in_pymol_rds.setvalue("No")
            self.middle_formlayout.add_widget_to_align(self.select_in_pymol_rds)

        self.middle_formlayout.set_input_widgets_width(150)


    def get_residue_range(self):
        return self.residue_range_enf.getvalue()

    def get_feature_name(self):
        return self.feature_name_enf.getvalue()

    def pick_color_dialog(self):
        selected_color = open_color_dialog(color_format="all")
        if selected_color is not None:
            self.selected_rgb, self.selected_hex = selected_color
            # self.feature_color_enf.setvalue(self.selected_hex)
            self.feature_color_enf.entry.setStyleSheet("background-color: %s" % self.selected_hex)

    def get_selected_colors(self):
        return self.selected_rgb, self.selected_hex

    def get_select_in_pymol(self):
        return _get_select_in_pymol(self)


def _get_select_in_pymol(window):
    if window.pymod_element.has_structure():
        choice_value = window.select_in_pymol_rds.getvalue()
        if choice_value == "Yes":
            return True
        elif choice_value == "No":
            return False
        else:
            raise KeyError(choice_value)
    else:
        return False


################################################################
# Window for searching a string in a sequence loaded in PyMod. #
################################################################

class Search_string_window_qt(PyMod_tool_window_qt):

    inactive_results_color = inactive_bg_color
    found_results_color = success_bg_color
    not_found_results_color = failure_bg_color
    default_message = "Type a pattern and press Enter..."

    def __init__(self, parent, pymod_elements, *args, **configs):

        self.pymod_element = pymod_elements[0]
        PyMod_tool_window_qt.__init__(self, parent, *args, **configs)


    def add_middle_frame_widgets(self):

        #-------------------------------------
        # Form layout for input and results. -
        #-------------------------------------

        entries_width = 340

        # Entryfield for inserting a subsequence.
        self.search_string_enf = PyMod_entryfield_qt(label_text="Search For", value="",
                                                     enter_command=self.submit_command)
        self.search_string_enf.set_input_widget_width(entries_width)
        self.middle_formlayout.add_widget_to_align(self.search_string_enf, align=False)


        # Entryfield for showing the results.
        self.results_enf = PyMod_entryfield_qt(label_text="",
                                               value=self.default_message,
                                               readonly=True)
        self.set_results_entry_bg(self.inactive_results_color)
        self.results_enf.set_input_widget_width(entries_width)
        self.middle_formlayout.add_widget_to_align(self.results_enf, align=False)

        #---------------------------
        # Form layout for options. -
        #---------------------------

        # Use regular expressions.
        self.use_regex_rds = PyMod_radioselect_qt(label_text="Use Regex",
                                                  buttons=("Yes", "No"))
        self.use_regex_rds.setvalue("No")
        self.middle_formlayout.add_widget_to_align(self.use_regex_rds)

        # Highlight color selection.
        color_buttons = ("yellow", "red", "green", "cyan") # "violet"
        self.highlight_color_rds = PyMod_radioselect_qt(label_text='Highlight Color',
                                                        buttons=color_buttons)
        self.highlight_color_rds.setvalue('yellow')
        self.middle_formlayout.add_widget_to_align(self.highlight_color_rds)


        if self.pymod_element.has_structure():
            self.select_in_pymol_rds = PyMod_radioselect_qt(label_text='Select in PyMOL',
                                                            buttons=("Yes", "No"))
            self.select_in_pymol_rds.setvalue('No')
            self.middle_formlayout.add_widget_to_align(self.select_in_pymol_rds)

        self.middle_formlayout.set_input_widgets_width(width="auto")

    def get_search_string(self):
        return self.search_string_enf.getvalue()

    def get_regex_use(self):
        use_regex_val = self.use_regex_rds.getvalue()
        if use_regex_val == "Yes":
            return True
        elif use_regex_val == "No":
            return False
        else:
            raise KeyError(use_regex_val)

    def get_highlight_color(self):
        return self.highlight_color_rds.getvalue()

    def show_results(self, message, state):
        self.results_enf.setvalue(message)

        if state == "found":
            self.set_results_entry_bg(self.found_results_color)
        elif state == "not_found":
            self.set_results_entry_bg(self.not_found_results_color)
        elif state == "empty":
            self.set_results_entry_bg(self.inactive_results_color)
        else:
            raise KeyError(state)

    def set_results_entry_bg(self, color):
        self.results_enf.input.setStyleSheet("background-color: %s; color: black" % color)

    def get_select_in_pymol(self):
        return _get_select_in_pymol(self)
