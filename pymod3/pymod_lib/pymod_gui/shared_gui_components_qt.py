# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Common functions and variables used by the Qt GUI of PyMod.
"""

import os

from pymol.Qt import QtWidgets, QtCore, QtGui

from pymod_lib.pymod_vars import convert_hex_to_rgb


###############################################################################
# Menus.                                                                      #
###############################################################################

def add_qt_menu_command(parent, label, command=None, fg_color=None, bg_color=None):
    """
    Adds to a 'parent' QMenu widget a new menu item.
    """

    # Color the label of the menu item.
    if fg_color != None:
        action = QtWidgets.QWidgetAction(parent)
        label_widget = QtWidgets.QLabel(label)
        s = """QLabel {
            background-color: %s;
            color: %s;
            padding: 3px;
        }

        QLabel:hover {
            background-color: #466e82;
            color: %s;
        }""" % (bg_color, fg_color, fg_color)

        label_widget.setStyleSheet(s)
        action.setDefaultWidget(label_widget)

    # Don't color, use a regular 'QAction' object.
    else:
        action = QtWidgets.QAction(label, parent)

    # Associates a command.
    if command is not None:
        action.triggered.connect(command)
    parent.addAction(action)

    return action


###############################################################################
# Dialogs.                                                                    #
###############################################################################

def askyesno_qt(title, message, parent=None, buttons_text=None):
    """
    Wrapper to a Yes/no dialog in PyQt. If 'buttons_text' is 'None', the default
    "Yes" and "No" buttons will be used. If 'buttons_text' is a list with two
    strings, the first string will be the text of the "Yes" button and the second
    one will be the text of the "No" button.
    """

    # Use yes and no buttons.
    if buttons_text is None:
        answer = QtWidgets.QMessageBox.question(parent, title, message,
                                                QtWidgets.QMessageBox.Yes,
                                                QtWidgets.QMessageBox.No)
        return answer == QtWidgets.QMessageBox.Yes

    # Set custom text on the buttons.
    else:
        dialog = QtWidgets.QMessageBox(parent)
        dialog.setWindowTitle(title)
        dialog.setText(message)
        yesbutton = dialog.addButton(buttons_text[0], QtWidgets.QMessageBox.YesRole)
        nobutton = dialog.addButton(buttons_text[1], QtWidgets.QMessageBox.NoRole)
        answer = dialog.exec_()
        return dialog.clickedButton() is yesbutton


def askopenfile_qt(title, parent=None, initialdir="", initialfile=None, name_filter=""):
    """
    Wrapper to a show a "pick a file to open" dialog in PyQt.
    """
    askfile_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""

    askfile_dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)

    filepath = askfile_dialog.getOpenFileName(parent, title, _initialdir, name_filter)
    if isinstance(filepath, (tuple, list)):
        filepath = filepath[0]
    else:
        filepath = str(filepath)

    return filepath


def askopenfiles_qt(title, parent=None, initialdir="", name_filter=""):
    """
    Wrapper to a show a "pick multiple files to open" dialog in PyQt.
    """
    askfile_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""

    askfile_dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)

    _filepaths = askfile_dialog.getOpenFileNames(parent, title, _initialdir, name_filter)
    if isinstance(_filepaths, (tuple, list)):
        filepaths = _filepaths[0]
    else:
        filepaths = _filepaths

    return filepaths


def askdirectory_qt(title, parent=None, initialdir=""):
    """
    Wrapper to a show a "pick a directory" dialog in PyQt.
    """
    askdirectory_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""
    flags = QtWidgets.QFileDialog.ShowDirsOnly
    dirpath = str(askdirectory_dialog.getExistingDirectory(parent, title, _initialdir, flags))
    return dirpath


def asksaveasfile_qt(title, parent=None, initialdir="", name_filter="", check_existent=True):
    """
    Wrapper to a show a "pick a file to save" dialog in PyQt.
    """
    askfile_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""

    _filepath = askfile_dialog.getSaveFileName(parent, title, _initialdir, name_filter)
    if isinstance(_filepath, (tuple, list)):
        filepath = _filepath[0]
        sel_filter = _filepath[1]
    else:
        filepath = str(_filepath)
        sel_filter = ""

    if not filepath:
        return filepath

    if name_filter and check_existent:
        if sel_filter:
            extension = sel_filter[1:]
        else:
            return filepath

        # PyQt has already checked for the existance of the file.
        if filepath.endswith(extension):
            return filepath
        # PyQt may not have seen the file with the extension.
        else:
            if os.path.isfile(filepath + extension):
                title = "Save As"
                message = "File '%s' already exists. Do you want to replace it?" % (filepath + extension)
                choice = askyesno_qt(title, message, parent=parent)
                if choice:
                    return filepath + extension
                else:
                    return ""
            else:
                return filepath + extension
    else:
        return filepath


def open_color_dialog(color_format="all"):
    """
    Opens a Qt color dialog and return a string encoding the selected color.
    """

    if not color_format in ("rgb", "hex", "all"):
        raise KeyError("Unknown 'color_format': %s" % color_format)

    color = QtWidgets.QColorDialog.getColor()

    if color.isValid():
        color_hex = color.name()
        color_rgb = convert_hex_to_rgb(color_hex)
        if color_format == "rgb":
            return color_rgb
        elif color_format == "hex":
            return color_hex
        elif color_format == "all":
            return color_rgb, color_hex

    return None


###############################################################################
# Qt windows used in PyMod.                                                   #
###############################################################################

class PyMod_tool_window_qt(QtWidgets.QMainWindow):
    """
    Class for various types of windows in PyMod.
    """

    middle_layout_type = "qform"
    is_pymod_window = True

    def __init__(self, parent,
                 title="New PyMod Window",
                 upper_frame_title="New PyMod Window Sub-title",
                 submit_command=None, submit_button_text="Submit",
                 with_scroll=True,
                 # geometry=None
                 ):

        super(PyMod_tool_window_qt, self).__init__(parent)

        #------------------------
        # Configure the window. -
        #------------------------

        # Command executed when pressing on the main button of the window.
        self.submit_command = submit_command

        # Configure the window.
        self.setWindowTitle(title)
        # if geometry is not None:
        #     self.setGeometry(*geometry)

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        #---------------
        # Upper frame. -
        #---------------

        self.upper_frame_title = QtWidgets.QLabel(upper_frame_title)
        self.main_vbox.addWidget(self.upper_frame_title)


        #----------------
        # Middle frame. -
        #----------------

        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # The Vertical Box that contains other widgets to be displayed in the window.
        self.middle_vbox = QtWidgets.QVBoxLayout()
        self.middle_widget.setLayout(self.middle_vbox)

        # Scroll area which contains the widgets, set as the centralWidget.
        self.middle_scroll = QtWidgets.QScrollArea()

        # Scroll area properties.
        # self.middle_scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        # self.middle_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        # QFormLayout in the middle frame.
        if self.middle_layout_type == "qform":
            self.middle_formlayout = PyMod_QFormLayout()
            self.middle_vbox.addLayout(self.middle_formlayout)
        elif self.middle_layout_type == "qgrid":
            self.middle_formlayout = QtWidgets.QGridLayout()
            self.middle_vbox.addLayout(self.middle_formlayout)
        else:
            raise KeyError("Unknown 'middle_layout_type': %s" % middle_layout_type)

        self.add_middle_frame_widgets()

        self.main_vbox.addWidget(self.middle_scroll)


        #----------------
        # Bottom frame. -
        #----------------

        self.submit_command = submit_command
        if self.submit_command is not None:
            self.main_button = QtWidgets.QPushButton(submit_button_text)
            self.main_button.clicked.connect(lambda a=None: self.submit_command())
            self.main_vbox.addWidget(self.main_button)
            self.main_button.setFixedWidth(self.main_button.sizeHint().width())


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


    def add_middle_frame_widgets(self):
        """
        To be overriden in children classes. Add widgets to the 'middle_vbox' by using:
            self.middle_vbox.addWidget(widget)
        """
        pass


class PyMod_protocol_window_qt(PyMod_tool_window_qt):
    """
    Class for showing a window with options for several PyMod protocols.
    """

    def __init__(self, parent, protocol, *args, **configs):

        self.protocol = protocol
        PyMod_tool_window_qt.__init__(self, parent=parent, *args, **configs)
        self.showing_advanced_widgets = False
        self.additional_initialization()
        self.build_protocol_middle_frame()


    # Methods to be overriden in child classes.
    def build_protocol_middle_frame(self):
        pass

    def additional_initialization(self):
        pass

    # Methods for showing widgets for advanced options.
    def show_advanced_button(self):
        self.advance_options_button = QtWidgets.QPushButton("Show Advanced Options")
        self.advance_options_button.clicked.connect(self.toggle_advanced_options)
        self._advanced_options_label = QtWidgets.QLabel("")
        self.middle_formlayout.addRow(self.advance_options_button, self._advanced_options_label)


    def toggle_advanced_options(self):
        if not self.showing_advanced_widgets:
            self.showing_advanced_widgets = True
            self.advance_options_button.setText("Hide Advanced Options")
            for row in self.middle_formlayout.widgets_to_align:
                if row.is_advanced_option:
                    row.show_widgets()
        else:
            self.showing_advanced_widgets = False
            self.advance_options_button.setText("Show Advanced Options")
            for row in self.middle_formlayout.widgets_to_align:
                if row.is_advanced_option:
                    row.hide_widgets()


    def check_general_input(self):
        """
        Raises an exception if the input is not valid.
        """

        for row in self.middle_formlayout.widgets_to_align:
            if row.to_be_validated:
                if row.is_advanced_option and not self.showing_advanced_widgets:
                    continue
                if hasattr(row, "validate_input"):
                    row.validate_input()
        return None


    def show(self):
        PyMod_tool_window_qt.show(self)
        if hasattr(self, "advance_options_button"):
            self.advance_options_button.setFixedWidth(self.advance_options_button.sizeHint().width())
            for row in self.middle_formlayout.widgets_to_align:
                if row.is_advanced_option:
                    row.hide_widgets()


###############################################################################
# Qt widgets used in PyMod.                                                   #
###############################################################################

default_width_hint = 10


class PyMod_QFormLayout(QtWidgets.QFormLayout):
    """
    A custom 'QFormLayout' for many of PyMod input widgets.
    """

    def __init__(self, vertical_spacing=1, *args, **kwargs):
        QtWidgets.QFormLayout.__init__(self, *args, **kwargs)
        # self.setVerticalSpacing(vertical_spacing)
        self.widgets_to_align = []

    def add_widget_to_align(self, widget, advanced_option=False, validate=False, align=True):
        if align:
            self.widgets_to_align.append(widget)
        widget.is_advanced_option = advanced_option
        widget.to_be_validated = validate
        self.addRow(widget.label, widget.input)

    def set_input_widgets_width(self, width, min_width=60, max_width=200, padding=30):
        if width != "auto":
            for widget in self.widgets_to_align:
                widget.set_input_widget_width(width)
        else:
            widths = [widget.get_width_hint() for widget in self.widgets_to_align]
            _max_width = max(widths)
            if _max_width > max_width:
                _max_width = max_width
            if _max_width < min_width:
                _max_width = min_width
            for widget in self.widgets_to_align:
                widget.set_input_widget_width(_max_width+padding)


class PyMod_form_item(QtWidgets.QWidget):

    pass


class PyMod_entry_qt(QtWidgets.QLineEdit):

    def __init__(self, *args, **kwargs):
        super(PyMod_entry_qt, self).__init__(*args, **kwargs)
        self.pmw_validator = {}


    def set_pmw_validator(self, validator):
        self.pmw_validator = validator
        if not self.pmw_validator:
            return None
        if self.pmw_validator["validator"] == "integer":
            self.setValidator(QtGui.QIntValidator(self.pmw_validator["min"], self.pmw_validator["max"]))
        elif self.pmw_validator["validator"] == "real":
            self.setValidator(QtGui.QDoubleValidator(self.pmw_validator["min"], self.pmw_validator["max"], 9))
        else:
            raise KeyError("Unknown 'validator': %s" % self.pmw_validator["validator"])


    def getvalue(self, validate=False, option_name="Option"):

        # Just return the value from the GUI.
        if not validate:
            return self.text()

        # Returns the value from the GUI, but first validate it. If it can not be
        # validated, raises an exception.
        else:
            # A validator was provided.
            if self.pmw_validator:
                # Integers and floats.
                if self.pmw_validator["validator"] in ("integer", "real"):
                    if self.pmw_validator["validator"] == "integer":
                        try:
                            val = int(self.text())
                        except ValueError as e:
                            raise ValueError("Invalid input for '%s': could not convert to integer." % option_name)
                    else:
                        try:
                            val = float(self.text())
                        except ValueError as e:
                            raise ValueError("Invalid input for '%s': could not convert to float." % option_name)

                    if not self.pmw_validator["min"] <= val <= self.pmw_validator["max"]:
                        message = ("Invalid input for '%s'. The value must be in the following"
                                   " range: %s to %s." % (option_name,
                                                          self.pmw_validator["min"],
                                                          self.pmw_validator["max"]))
                        raise ValueError(message)
                    return val

                else:
                    raise KeyError("Unknown 'validator': %s" % self.pmw_validator["validator"])
            else:
                return self.text()

    def setvalue(self, value):
        self.setText(value)


class PyMod_entryfield_qt(PyMod_form_item):
    """
    Class for a entryfield widget in PyQt. Designed to be used in 'QFormLayout' GUIs.
    """

    def __init__(self, label_text="Input", value="",
                 readonly=False, style=None,
                 enter_command=None,
                 validate={}):
        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Entry.
        self.entry = PyMod_entry_qt(value)
        self.enter_command = enter_command
        if self.enter_command is not None:
            self.entry.returnPressed.connect(self.enter_command)
        if readonly:
            self.entry.setReadOnly(True)

        if style is not None:
            self.entry.setStyleSheet(style)
        else:
            self.entry.setStyleSheet(active_entry_style)

        self.validate = validate
        if self.validate:
            self.entry.set_pmw_validator(self.validate)

        self.input = self.entry


    def setvalue(self, value):
        self.entry.setvalue(value)

    def getvalue(self, validate=False):
        return self.entry.getvalue(validate=validate, option_name=self.label.text())

    def set_input_widget_width(self, width):
        self.entry.setFixedWidth(width)

    def get_width_hint(self):
        return default_width_hint # self.entry.sizeHint().width()

    def show_widgets(self):
        self.label.show()
        self.entry.show()

    def hide_widgets(self):
        self.label.hide()
        self.entry.hide()

    def validate_input(self):
        return self.getvalue(validate=True)


class PyMod_plaintextedit_qt(PyMod_form_item):
    """
    Class for a plain text edit widget in PyQt. Designed to be used in 'QFormLayout'
    GUIs.
    """

    def __init__(self, label_text="Input", value="", style=None):

        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Entry.
        self.entry = QtWidgets.QPlainTextEdit(value)
        if style is not None:
            self.entry.setStyleSheet(style)
        else:
            self.entry.setStyleSheet(active_entry_style)
        self.entry.setWordWrapMode(QtGui.QTextOption.WrapAnywhere)
        expanding_size_policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                                      QtWidgets.QSizePolicy.Expanding)
        expanding_size_policy.setVerticalStretch(1)
        self.entry.setSizePolicy(expanding_size_policy)

        self.input = self.entry


    def setvalue(self, value):
        self.entry.setPlainText(value)

    def getvalue(self, validate=False):
        return self.entry.toPlainText()

    def set_input_widget_width(self, width):
        self.entry.setFixedWidth(width)

    def get_width_hint(self):
        return default_width_hint # self.entry.sizeHint().width()

    def show_widgets(self):
        self.label.show()
        self.entry.show()

    def hide_widgets(self):
        self.label.hide()
        self.entry.hide()

    def validate_input(self):
        return self.getvalue(validate=True)


class PyMod_entryfield_button_qt(PyMod_entryfield_qt):
    """
    Class for a entryfield with a button widget in PyQt. Designed to be used in
    'QFormLayout' GUIs.
    """

    def __init__(self, button_text="Submit", button_command=None, *args, **kwargs):
        PyMod_entryfield_qt.__init__(self, *args, **kwargs)

        self.input = QtWidgets.QHBoxLayout()
        self.input.addWidget(self.entry)

        # Adds a button.
        self.button = QtWidgets.QPushButton(button_text)
        self.input.addWidget(self.button)
        self.button.setFixedWidth(self.button.sizeHint().width())
        self.button_command = button_command
        if self.button_command is not None:
            self.button.clicked.connect(self.button_command)

    def show_widgets(self):
        PyMod_entryfield_qt.show_widgets(self)
        self.button.show()

    def hide_widgets(self):
        PyMod_entryfield_qt.hide_widgets(self)
        self.button.hide()


class PyMod_radioselect_qt(PyMod_form_item):
    """
    Class for a radioselect widget in PyQt. Designed to be used in 'QFormLayout' GUIs.
    """

    def __init__(self, label_text="Input", buttons=[]):
        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Buttons.
        self.input = QtWidgets.QVBoxLayout()

        if not buttons:
            raise ValueError("Please provide a list of button names")
        if len(buttons) != len(set(buttons)):
            raise ValueError("Please provide a non redundant list of buttons")

        self.button_group = QtWidgets.QButtonGroup()
        self.buttons_names = []
        self.buttons_dict = {}

        for button_name in buttons:
            button = QtWidgets.QPushButton(button_name)
            button.setCheckable(True)
            self.input.addWidget(button)
            self.buttons_names.append(button_name)
            self.buttons_dict[button_name] = button
            self.button_group.addButton(button)

    def get_buttons(self):
        return self.button_group.buttons()

    def get_button_at(self, index):
        buttons = [b for b in self.get_buttons()]
        return buttons[index]

    def setvalue(self, value):
        self.buttons_dict[value].setChecked(True)

    def getvalue(self):
        checked_button = self.button_group.checkedButton()
        if checked_button is None:
            return None
        else:
            return checked_button.text()

    def set_input_widget_width(self, width):
        for button in self.button_group.buttons():
            button.setFixedWidth(width)

    def get_width_hint(self):
        return max([button.sizeHint().width() for button in self.button_group.buttons()])

    def show_widgets(self):
        self.label.show()
        for button in self.button_group.buttons():
            button.show()

    def hide_widgets(self):
        self.label.hide()
        for button in self.button_group.buttons():
            button.hide()


class PyMod_combobox_qt(PyMod_form_item):
    """
    Class for a combobox widget in PyQt. Designed to be used in 'QFormLayout' GUIs.
    """

    def __init__(self, label_text="Input", items=[]):
        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Combobox.
        self.combobox = QtWidgets.QComboBox()

        if not items:
            raise ValueError("Please provide a list of items for the combobox")
        self.items = items

        for item in self.items:
            self.combobox.addItem(item)
        self.combobox.setEditable(False)
        self.input = self.combobox

    def get(self):
        return self.combobox.currentText()

    def get_index(self):
        return self.combobox.currentIndex()

    def set_input_widget_width(self, width):
        self.combobox.setFixedWidth(width)

    def get_width_hint(self):
        return self.combobox.sizeHint().width()

    def show_widgets(self):
        self.label.show()
        self.combobox.show()

    def hide_widgets(self):
        self.label.hide()
        self.combobox.hide()


class PyMod_entrylabel_qt(PyMod_form_item):
    """
    Class for a label widget in PyQt. Designed to be used in 'QFormLayout' GUIs.
    """

    def __init__(self, label_text="Input", value=""):
        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Second Entry.
        self.right_label = QtWidgets.QLabel(value)
        self.right_label.setWordWrap(True)
        self.input = self.right_label

    def set_input_widget_width(self, width):
        self.right_label.setFixedWidth(width)

    def get_width_hint(self):
        return self.right_label.sizeHint().width()

    def show_widgets(self):
        self.label.show()
        self.right_label.show()

    def hide_widgets(self):
        self.label.hide()
        self.right_label.hide()


class PyMod_hbox_option_qt(PyMod_form_item):
    """
    Class for a combobox widget in PyQt. Designed to be used in 'QFormLayout' GUIs.
    """

    def __init__(self, label_text="Input"):
        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Hbox.
        self.hbox = QtWidgets.QHBoxLayout()
        self.input = self.hbox


    def set_auto_input_widget_width(self):
        for idx in range(0, self.hbox.count()):
            widget = self.hbox.itemAt(idx).widget()
            if hasattr(widget, "setFixedWidth"):
                widget.setFixedWidth(widget.sizeHint().width())

    def set_input_widget_width(self, width):
        pass

    def get_width_hint(self):
        return default_width_hint

    def show_widgets(self):
        self.label.show()
        for idx in range(0, self.hbox.count()):
            widget = self.hbox.itemAt(idx).widget()
            if hasattr(widget, "show"):
                widget.show()

    def hide_widgets(self):
        self.label.hide()
        for idx in range(0, self.hbox.count()):
            widget = self.hbox.itemAt(idx).widget()
            if hasattr(widget, "hide"):
                widget.hide()


class PyMod_spinbox_entry_qt(PyMod_form_item):
    """
    Class for a entryfield widget in PyQt. Designed to be used in 'QFormLayout' GUIs.
    """

    def __init__(self, label_text="Input", value=1,
                 spinbox_min=1, spinbox_max=100):
        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Spinbox.
        self.spinbox = QtWidgets.QSpinBox()
        self.spinbox_min = spinbox_min
        self.spinbox_max = spinbox_max
        self.spinbox.setRange(self.spinbox_min, self.spinbox_max)
        self.spinbox.setStyleSheet(active_entry_style)

        self.input = self.spinbox


    def setvalue(self, value):
        self.spinbox.setValue(value)

    def getvalue(self, validate=False):
        # Just return the value from the GUI.
        if not validate:
            return self.spinbox.value()

        # Returns the value from the GUI, but first validate it. If it can not be
        # validated, raises an exception.
        else:
            option_name = self.label.text()
            try:
                val = int(self.spinbox.value())
            except ValueError as e:
                raise ValueError("Invalid input for '%s': could not convert to integer." % option_name)

            if not self.spinbox_min <= val <= self.spinbox_max:
                message = ("Invalid input for '%s'. The value must be in the following"
                           " range: %s to %s." % (option_name,
                                                  self.spinbox_min,
                                                  self.spinbox_max))
                raise ValueError(message)
            return val

    def set_input_widget_width(self, width):
        self.spinbox.setFixedWidth(width)

    def get_width_hint(self):
        return self.spinbox.sizeHint().width()

    def show_widgets(self):
        self.label.show()
        self.spinbox.show()

    def hide_widgets(self):
        self.label.hide()
        self.spinbox.hide()

    def validate_input(self):
        return self.getvalue(validate=True)


class PyMod_scalebar_qt(PyMod_form_item):
    """
    Class for a scalerbar widget in PyQt. Designed to be used in 'QFormLayout' GUIs.
    """

    def __init__(self, label_text="New scalebar",
                 slider_value=1,
                 slider_from=1, slider_to=10,
                 slider_resoution=1,
                 # slider_digits=3,
                 slider_tickinterval=1,
                 slider_use_float=False, slider_use_float_val=100.0,
                 slider_binding=None,
                 slider_width=None):

        PyMod_form_item.__init__(self)

        # Label.
        self.label = QtWidgets.QLabel(label_text)

        # Layout for the input widget and its label.
        self.input = QtWidgets.QHBoxLayout()

        # Adds a slider.
        self.slider = QtWidgets.QSlider()
        self.slider.setOrientation(QtCore.Qt.Horizontal)
        self.slider.setTickPosition(QtWidgets.QSlider.TicksBelow)

        self.slider_use_float = slider_use_float
        self.slider_use_float_val = slider_use_float_val
        self.slider.slider_resoution = slider_resoution
        self.slider.setMinimum(round(self._get_slider_val(slider_from, internal=True)))
        self.slider.setMaximum(round(self._get_slider_val(slider_to, internal=True)))
        self.slider.setValue(round(self._get_slider_val(slider_value, internal=True)))
        self.slider.setTickInterval(round(self._get_slider_val(slider_tickinterval, internal=True)))
        self.slider.setSingleStep(round(self._get_slider_val(slider_resoution, internal=True)))
        self.slider.setPageStep(round(self._get_slider_val(slider_tickinterval, internal=True)))

        self.slider.valueChanged.connect(self._on_slider_change)
        self.slider.sliderPressed.connect(self._on_slider_pressed)
        self.slider.sliderReleased.connect(self._on_slider_release)
        self.on_drag = False

        self.input.addWidget(self.slider)

        # Add a label on the right of the slider.
        self.slider_label = QtWidgets.QLabel(str(slider_to))
        self.input.addWidget(self.slider_label)

        if slider_width:
            self.slider.setFixedWidth(slider_width)
            self.slider_label.setFixedWidth(self.slider_label.sizeHint().width())
        self.slider_label.setText(str(slider_value))

        self.slider_binding = slider_binding


    def show_widgets(self):
        self.label.show()
        self.slider.show()
        self.slider_label.show()

    def hide_widgets(self):
        self.label.hide()
        self.slider.hide()
        self.slider_label.hide()


    def _on_slider_change(self):
        val = self._get_slider_val(self.slider.value(), internal=False)
        self.slider_label.setText(str(val))
        if not self.on_drag:
            self.call_slider_biding()

    def _on_slider_pressed(self):
        self.on_drag = True

    def _on_slider_release(self):
        self.call_slider_biding()
        self.on_drag = False

    def call_slider_biding(self):
        if self.slider_binding is not None:
            self.slider_binding()


    def _get_slider_val(self, val, internal=True):
        if not self.slider_use_float:
            return val
        else:
            if internal:
                return val*self.slider_use_float_val
            else:
                return val/self.slider_use_float_val

    def getvalue(self):
        return self._get_slider_val(self.slider.value(), internal=False)

    def get(self):
        return self.getvalue()


    def set_input_widget_width(self, width):
        pass

    def set_auto_input_widget_width(self, width):
        pass

    def get_width_hint(self):
        return self.slider.sizeHint().width() + self.slider_label.sizeHint().width()


###############################################################################
# CSS Styles used in PyMod.                                                   #
###############################################################################

default_pt_size = QtGui.QFont().pointSize()


active_entry_style = "background-color: white; color: #333333"
inactive_entry_style = "background-color: #ccc; color: #7c7c7c"

inactive_bg_color = "#e3e3e3"
success_bg_color = "#98fb98"
failure_bg_color = "#f08080"

highlight_color = "#5ac8ff"

options_title_style = "font-size: %spt; color: %s" % (default_pt_size+1, highlight_color)
small_font_style = "font-size: %spt" % (default_pt_size-1)
large_font_style = "font-size: %spt" % (default_pt_size+1)
