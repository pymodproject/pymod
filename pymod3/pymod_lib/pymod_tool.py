# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Classes used to control PyMod tools (the external programs that PyMod uses).
"""

import os
from . import pymod_os_specific as pmos

from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt
from pymol.Qt import QtWidgets


#####################################################################
# Tools classes.                                                    #
#####################################################################

class Tool:
    """
    A base class to represent a Tool, an external program used by PyMod.
    """

    def __init__(self, tool_name, tool_full_name):
        # Id of the tool. Example: 'clustalo'.
        self.name = tool_name
        # Full name of the tool. Example: 'Clustal Omega'.
        self.full_name = tool_full_name
        # A list of 'Tool_parameter' class objects.
        self.parameters = []


    def show_message_method(self, title, message):
        """
        Method used to display messages raised when using the tool.
        """
        self.pymod.main_window.show_error_message(title, message)


    def initialize_parameters(self, parameter_list):
        """
        Used to populate the 'self.parameter' attribute with a list of object of the Tool_parameter
        class and its child classes.
        """
        self.parameters = [] # Reinitializes.
        for p in parameter_list:
            p.set_parent_tool(self)
            self.parameters.append(p)


    def __getitem__(self, parameter_name):
        for parameter in self.parameters:
            if parameter.name == parameter_name:
                return parameter
        raise KeyError("'%s' tool does not have a parameter named: '%s'" % (self.name, parameter_name))


class Executable_tool(Tool):

    #--------------------
    # Executable files. -
    #--------------------

    def tool_file_exists(self, parameter_name="exe_file_path"):
        """
        Check if the executable is defined in the PyMod options and if the file exists.
        """
        if self[parameter_name].get_value():
            return os.path.isfile(self[parameter_name].get_value())
        else:
            return False


    def tool_file_not_found(self, parameter_name="exe_file_path"):
        # If the exe file path is defined in the PyMod configurations but the file does not exists.
        if self[parameter_name].get_value():
            title = "%s is Missing" % (self.full_name)
            message = "The %s path specified in the PyMod Options Window ('%s') does not exists on your system. If you want to use this function please specify an existent %s executable path in the Options of PyMod." % (self.full_name, self[parameter_name].get_value(), self.full_name)
        # If the exe file path is not defined in the PyMod configurations.
        else:
            title = "%s Path is Missing" % (self.full_name)
            message = "%s executable path is missing from the PyMod Options. If you want to use this function please install %s and specify its path in the PyMod Options Window." % (self.full_name, self.full_name)
        self.show_message_method(title, message)


    #-------------------------------------
    # Directories with executable files. -
    #-------------------------------------

    def tool_dir_exists(self, parameter_name="exe_dir_path"):
        """
        Check if the executable is defined in the PyMod options and if the file exists.
        """
        if self[parameter_name].get_value():
            return os.path.isdir(self[parameter_name].get_value())
        else:
            return False


    def tool_dir_not_found(self, parameter_name="exe_dir_path", message=None):
        if self[parameter_name].get_value():
            title = "%s Path not Found" % (self.full_name)
            _message = "The %s executables directory specified in the PyMod Options Window ('%s') does not exists on your system. If you want to use this function please specify an existent %s executables directory in the Options of PyMod." % (self.full_name, self[parameter_name].get_value(), self.full_name)
        else:
            title = "%s Path is Missing" % (self.full_name)
            _message = "%s executable path is missing from the PyMod Options. If you want to use this function please install %s and specify its path in the PyMod Options Window." % (self.full_name, self.full_name)
        message = _message if message == None else message
        self.show_message_method(title, message)


#####################################################################
# Tool parameters classes.                                          #
#####################################################################

class Tool_parameter:
    """
    A base class to represent a parameter of a Tool. Its child classes will be used throughout the
    plugin to represent different types of options of PyMod tools.
    """

    widget_type = None

    def __init__(self, paramater_name, parameter_full_name,
                 default_value=None,
                 editable=True,
                 auto_find=False,
                 show_widget=True):
        # Full name of the parameter. This will be displayed in the PyMod options window.
        # Example: "Executable File".
        self.full_name = parameter_full_name
        # Id of the parameter. Example: "exe_full_path"
        self.name = paramater_name
        # The value of the parameter. Example: "/usr/bin/clustalo"
        self.value = None

        # Execution mode. It can be "local", "remote", "both", "mixed".
        # self.execution_mode = parameter_execution_mode

        # This will be set through the 'set_parent_tool()' method.
        self.parent_tool = None

        # If this argument is specified in the constructor, its value will be the one to be returned
        # by the 'get_starting_value()' method of this class. This is used in order to set the value
        # of the parameter in the PyMod options window when PyMod is started for the first time.
        self.default_value = default_value

        # If 'True' allows to edit the paramater values from the PyMod options window.
        self.editable = editable

        # If 'True' allows to automatically guess the value of the parameter.
        self.auto_find = auto_find

        # If 'False' the widgets for the parameter will not be shown (used for "hidden" options).
        self.show_widget = show_widget


    def set_parent_tool(self, parent_tool):
        self.parent_tool = parent_tool


    def __repr__(self):
        return self.value


    def get_value(self):
        return self.value


    def set_value(self, new_value):
        self.value = new_value


    def get_starting_value(self):
        """
        Used to set the value of options in the PyMod options window when it is built.
        """
        # If PyMod is launched for the first time try to get a default value.
        if self.value == None:
            # If the parameter default value was defined in the object constructor this method will
            # return this value stored in 'self.default_value'.
            # Alternatively, the 'guess_first_session_value()' method is called to guess the starting
            # value.
            if not self.default_value == None:
                starting_text = self.default_value
            else:
                starting_text = self.guess_first_session_value()
        # Else get the value obtained from the PyMod main configuration file.
        else:
            starting_text = self.value
        return starting_text


    ##############################################
    # Methods to be overridden in child classes. #
    ##############################################

    def guess_first_session_value(self):
        """
        This method will try to guess a parameter value when PyMod is started for the first time.
        This will be overridden in the child classes of the 'Tool_parameter' class in order to be
        able to retrieve different information.
        """
        return ""


    def auto_find_command(self):
        """
        Override in child classes.
        """
        pass


    def can_be_updated_from_gui(self):
        if not self.show_widget: # Does not interface with the GUI.
            return False
        if not self.editable: # The value is displayed in the GUI, but is not editable.
            return False
        if self.widget_type == "show_status": # Only show the status.
            return False
        return True


class Tool_path(Tool_parameter):
    """
    A class to represent a path on the user's machine (both a file or a directory).
    """

    path_type = None
    widget_type = "path_entryfield"


    def auto_find_command(self, etry):
        """
        Automatically find the tool value and ask whether to set it in the GUI.
        """

        auto_find_results = self.guess_first_session_value()
        tool_param_text = "%s '%s'" % (self.parent_tool.full_name, self.full_name)
        if auto_find_results:
            message = ("%s was found in one of the system's directories (at: %s)."
                       " Would you like to use this file and set its value in the PyMod"
                       " Options window?" % (tool_param_text, auto_find_results))
            answer = askyesno_qt("File Automatically Found", message,
                                 parent=self.parent_tool.pymod.get_qt_parent())
            if answer:
                etry.setText(str(auto_find_results))
        else:
            message = ("%s was not found in any of the system's directories. Please"
                       " manually specify its path in the PyMod Options window through"
                       " the 'Browse' button." % (tool_param_text))
            self.parent_tool.pymod.main_window.show_info_message("File Not Found", message)


class Tool_file(Tool_path):
    """
    A class to represent a file path on the user machine.
    """
    path_type = "file"


class Tool_exec_file(Tool_file):
    """
    A class to represent a executable file path on the user machine.
    """
    def guess_first_session_value(self):
        """
        Automatically finds exe_full_path on the user's machine. Most tools have the same
        executable file name.
        """
        value = pmos.pymod_which(self.parent_tool.name)
        # 'pymod_which' returns 'None' if it doesn't find the executable.
        if value == None:
            value = ""
        return value


class Tool_directory(Tool_path):
    """
    A class to represent a directory path on the user machine.
    """
    path_type = "directory"


class Tool_exec_directory(Tool_directory):
    """
    A class to represent a directory path on the user machine containing executable files.
    """
    def guess_first_session_value(self):
        """
        Similar to 'Tool_exec_file' 'guess_first_session_value' method, but this one returns the name
        of the directory of an executable file localized by 'pymod_which'.
        """
        value = pmos.pymod_which(self.parent_tool.name)
        if value == None:
            value = ""
        else:
            value = os.path.dirname(value)
        return value


###################################################################################################
# Program specific classes.                                                                       #
###################################################################################################

#####################################################################
# MODELLER tool classes.                                            #
#####################################################################

# Tool.
class Modeller_tool(Executable_tool):

    def check_exception(self):
        """
        Checks if MODELLER can be launched, either internally or externally. Returns a 'None' value
        if MODELLER can be imported, otherwise returns an error message containing the string
        describing the type of exception raised.
        """
        modeller_exception = pmos.check_importable_modeller(get_exception=True)
        if modeller_exception is not None:
            if issubclass(modeller_exception.__class__, ImportError): # isinstance(modeller_exception, ImportError):
                error_message = ("MODELLER is missing. Please installed MODELLER if you want to use"
                                 " this function.")
            else:
                error_message = ("MODELLER is present on your system, but it could not be imported"
                                 " in PyMod for the following reason:\n\n%s" % modeller_exception)
        else:
            error_message = None
        return error_message


    def run_internally(self):
        """
        Checks if users want to run MODELLER internally. In PyMod 3, it is always set to 'True'.
        """
        return self["use_importable_modeller"].get_value()


# Parameters.
class Use_importable_modeller(Tool_parameter):

    widget_type = "show_status"

    def get_status(self):
        # Text of the label to inform the user if MODELLER libs can be imported in PyMOL.
        text_to_show = "Status: "
        if pmos.check_importable_modeller():
            text_to_show += "Available"
            status = True
        else:
            text_to_show += "Not available"
            status = False
        return text_to_show, status


    def guess_first_session_value(self):
        return True
