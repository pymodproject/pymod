import os
import sys
import pymod_os_specific as pmos
import pymod_gui as pmgi
from Tkinter import *

###################################################################################################
# Classes used to control PyMod tools (the external programs that PyMod uses).                    #                                                    #
###################################################################################################

#####################################################################
# Tools classes.                                                    #
#####################################################################

class Tool:
    """
    A base class to represent a Tool, an external program used by PyMod.
    """

    def __init__(self, tool_name, tool_full_name, show_message_method=None):
        # Id of the tool. Example: 'clustalo'.
        self.name = tool_name
        # Full name of the tool. Example: 'Clustal Omega'.
        self.full_name = tool_full_name
        # This can be "local", "external", "both", "mixed".
        # self.execution_mode = tool_execution_mode
        # A list of 'Tool_parameter' class objects.
        self.parameters = []
        # Method used to display messages raised when using the tool.
        self.show_message_method = show_message_method


    def initialize_parameters(self, parameter_list):
        """
        Used to populate the 'self.parameter' attribute with a list of object of the Tool_parameter
        class and its child classes.
        """
        self.parameters = []
        for p in parameter_list:
            p.set_parent_tool(self)
            self.parameters.append(p)


    def __getitem__(self, parameter_name):
        for parameter in self.parameters:
            if parameter.name == parameter_name:
                return parameter
        raise Exception("%s does not have a parameter named: '%s'" % (self.name, parameter_name))


    def display_options(self, target_frame):
        """
        Displays at list of option in the the target_frame contained in a target widget.
        next_row = target_frame.grid_size()[1]
        """
        # A list of Tkinter widgets for choosing the tools parameters values.
        self.parameters_widgets_list = []

        tool_full_name_label = Label(target_frame, text=self.full_name, **pmgi.label_style_1)
        tool_full_name_label.pack(**pmgi.pack_options_1)
        for i, parameter in enumerate(self.parameters):
            w = parameter.display_options(target_frame)
            # If the display options return a widget, adds it ot the Tool list of parameter widgets.
            if w != None:
                self.parameters_widgets_list.append(w)


class Executable_tool(Tool):

    def get_exe_file_path(self):
        """
        Shortcut method that allows to get the executable file path by a simple call.
        """
        try:
            return self["exe_file_path"].get_value()
        except:
            return None


    def exe_exists(self):
        """
        Another shortcut method.
        """
        path_exists = False
        try:
            if self["exe_file_path"].path_exists():
                path_exists = True
        except:
            pass
        print "path_exists:", path_exists
        return path_exists


    def exe_not_found(self):
        title = "%s is Missing" % (self.full_name)
        message = "If you want to use this function please install %s and specify \nits path in the Options of PyMod." % (self.full_name)
        # message = "If you want to use this function \nyou need to install MUSCLE and specify the '.exe' file in the Options of PyMod."
        # message = "PyMod has encountered a problem with Modeller.\nPlease check your Modeller license key."
        # message = "If you want to use this function you need to install\n the standalone version of BLAST and specify \nthe 'psiblast.exe' file in the Options of PyMod"
        self.show_message_method(title, message)


    def tool_problem(self):
        title = "%s Error" % (self.full_name)
        message = "PyMod has encountered a problem with %s." % (self.full_name)
        self.show_message_method(title, message)


#####################################################################
# Tool parameters classes.                                          #
#####################################################################

class Tool_parameter:
    """
    A base class to represent a parameter of a Tool. Its child classes will be used throughout the
    plugin to represent different types of options of PyMod tools.
    """
    def __init__(self, paramater_name, parameter_full_name, parameter_default_value = None, run_after_selection=None):
        # Full name of the parameter. This will be displayed in the PyMod options window.
        # Example: "Executable File".
        self.full_name = parameter_full_name
        # Id of the parameter. Example: "exe_full_path"
        self.name = paramater_name
        # The value of the parameter. Example: "/usr/bin/clustalo"
        self.value = None

        # Method to call when the user changes the input. Usually it will remain 'None'.
        self.run_after_selection = run_after_selection

        # Execution mode. It can be "local", "remote", "both", "mixed".
        # self.execution_mode = parameter_execution_mode

        # This will be set through the 'set_parent_tool()' method.
        self.parent_tool = None

        # If this argument is specified in the constructor, its value will be the one to be returned
        # by the 'get_starting_value()' method of this class. This is used in order to set the value
        # of the parameter in the PyMod options window when PyMod is started for the first time.
        self.default_value = parameter_default_value


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
            # NOTE: def get_first_session_value(self):
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

    def display_options(self, target_frame):
        """
        Used to display a series of widgets to choose the parameter value in PyMod options window.
        This is going to be overriden in the child classes of this class, according to the kind of
        parameter.
        """
        parameter_full_name_label = Label(target_frame, text=self.full_name, **pmgi.label_style_2)
        parameter_full_name_label.pack(**pmgi.pack_options_1)
        # If the widget has to be aligned with other widgets in PyMod options window, then it has
        # to be returned by this method.
        return None


    def guess_first_session_value(self):
        """
        This method will try to guess a parameter value when PyMod is started for the first time.
        This will be overridden in the child classes of the 'Tool_parameter' class in order to be
        able to retrieve different information.
        """
        return ""


    def get_value_from_gui(self):
        """
        Gets the value from the input widgets in the PyMod options window. This will be overridden
        in the child classes of the 'Tool_parameter' class in order to be able to retrieve different
        types of input from the GUI.
        """
        return ""


class Tool_path(Tool_parameter):
    """
    A class to represent a path on the user's machine (both a file or a directory).
    """

    path_type = None

    def display_options(self, target_frame):
        """
        Displays a custom entryfield with a 'Browse' button to choose a path on the user's machine.
        """
        current_path_type = "file"
        if self.__class__.path_type in ("file", "directory"):
            current_path_type = self.__class__.path_type

        starting_text = str(self.get_starting_value())

        self.pef = pmgi.PyMod_path_entryfield(target_frame,
            label_text = self.full_name,
            label_style = pmgi.label_style_2,
            value = starting_text,
            path_type = current_path_type, run_after_selection = self.run_after_selection,
            askpath_title = "Search for %s %s" % (self.parent_tool.full_name, self.full_name) )
        self.pef.pack(**pmgi.pack_options_1)

        # Returns the widget so that it can be aligned by the 'align_widgets()' method.
        return self.pef


    def get_value_from_gui(self):
        return self.pef.getvalue()


    def path_exists(self):
        if os.path.exists(self.value):
            return True
        else:
            return False


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
        # Automatically finds exe_full_path on the user's machine. Most tools have always the same
        # executable file name.
        value = pmos.pymod_which(self.parent_tool.name)
        # pymod_which() can return None if it doesn't find the executable.
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
        Similar to 'Tool_exec_file' 'guess_first_session_value()' method, but this one returns the name
        of the directory of an executable file localized by 'pymod_which()'.
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

#########
# Tool. #
#########
class Modeller_tool(Executable_tool):
    importable_modeller = False

    def can_be_launched(self):
        """
        Checks if MODELLER can be launched, either internally or externally.
        """
        can_be_launched = False
        if self.run_internally():
            if self.importable_modeller:
                can_be_launched = True
        else:
            if self.exe_exists():
                can_be_launched = True
        return can_be_launched

    def run_internally(self):
        """
        Checks if users want to run MODELLER internally.
        """
        use_importable = self["use_importable_modeller"].get_value()
        if use_importable == "True":
            return True
        else:
            return False


###############
# Parameters. #
###############
class Use_importable_modeller(Tool_parameter):

    def __init__(self, paramater_name, parameter_full_name, parameter_default_value = None):
        Tool_parameter.__init__(self, paramater_name, parameter_full_name, parameter_default_value = parameter_default_value)

    def display_options(self, target_frame):
        # Text of the label to inform the user if MODELLER libs can be imported in PyMOL.
        text_to_show = "Current status: "
        if self.parent_tool.importable_modeller:
            text_to_show += "available"
        else:
            text_to_show += "unavailable"
        self.rds = pmgi.Use_importable_modeller_radioselect (
            target_frame, importable_modeller=self.parent_tool.importable_modeller, initial_value=self.value,
            label_text = self.full_name, label_style = pmgi.label_style_2, run_after_selection=self.update_status,
            value = text_to_show)
        self.rds.pack(**pmgi.pack_options_1)
        return self.rds

    def update_status(self):
        if self.rds.getvalue() == "True":
            self.parent_tool["exe_file_path"].pef.hide_path_selector()
        else:
            path_to_show = self.parent_tool["exe_file_path"].get_value()
            self.parent_tool["exe_file_path"].pef.show_path_selector(path_to_show)

    def guess_first_session_value(self):
        if self.parent_tool.importable_modeller:
            return "True"
        else:
            return "False"

    def get_value_from_gui(self):
        return self.rds.getvalue()


class Modeller_exec_file(Tool_exec_file):

    def display_options(self, target_frame):
        starting_text = str(self.get_starting_value())
        self.pef = pmgi.Modeller_exec_entryfield(target_frame,
            label_text = self.full_name,
            label_style = pmgi.label_style_2,
            value = starting_text,
            path_type = "file", run_after_selection = self.run_after_selection,
            askpath_title = "Search for %s %s" % (self.parent_tool.full_name, self.full_name))
        self.pef.pack(**pmgi.pack_options_1)

        if self.parent_tool["use_importable_modeller"].get_value() == "True":
            self.pef.hide_path_selector()
        else:
            self.pef.show_path_selector(starting_text)

        return self.pef

    def get_value_from_gui(self):
        return self.pef.getvalue()
