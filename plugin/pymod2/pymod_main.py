###########################################################################
# PyMod 2: PyMOL Front-end to MODELLER and various other bioinformatics tools.
# Copyright (C) 2016 Giacomo Janson, Chengxin Zhang, Alessandro Paiardini
# Copyright (C) 2011-2012 Emanuele Bramucci & Alessandro Paiardini,
#                         Francesco Bossa, Stefano Pascarella
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
###########################################################################

# Tkinter.
from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import Pmw

# Python standard library.
import warnings
import sys
import urllib
import urllib2
import gzip
import os
import shutil
import subprocess
import webbrowser
import re
import zipfile
import pickle

# PyMOL modules.
import pymol
from pymol import cmd, stored, selector

# NumPy.
import numpy

# Biopython.
from Bio import SeqIO
from Bio import AlignIO
global has_phylo
try:
    from Bio import Phylo
    if hasattr(Phylo, "draw"):
        has_phylo = True
    else:
        has_phylo = False
except:
    has_phylo = False
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
import Bio.PDB
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import PPBuilder
try:
    from Bio.Align.Applications import ClustalOmegaCommandline
except:
    from pymod_lib.pymod_sup import ClustalOmegaCommandline

# Matplotlib.
global matplotlib_found
try:
    import matplotlib
    matplotlib_found = True
except:
    matplotlib_found = False

# PyMod modules.
from pymod_lib import pymod_campo as campo # Python implementation of the CAMPO algorithm.
from pymod_lib import pymod_os_specific as pmos # Different OS compatibility-related code.
from pymod_lib import pymod_sequence_manipulation as pmsm # General biological sequence manipulation.
from pymod_lib import pymod_gui as pmgi # Part of the graphical user interface of PyMod.
from pymod_lib import pymod_vars as pmdt # PyMod data used throughout the plugin.
from pymod_lib import pymod_tool as pm_tool # Classes to represent tools used within PyMod.
from pymod_lib import pymod_plot as pplt # Basic plots for building DOPE profiles and showing distance trees.
from pymod_lib import pymod_sup as pmsp # Supplementary code for PyMod.
from pymod_lib import pymod_updater as pmup # Updates PyMod fetching the latest stable version via network.

# CE-alignment.
global ce_alignment_mode
try:
    # Try to import the ccealign module.
    from pymod_lib.ccealign import ccealign
    ce_alignment_mode = "plugin"
except:
    if pmos.check_pymol_builtin_cealign():
        ce_alignment_mode = "pymol"
    else:
       ce_alignment_mode = None

# Function that launches PyMod from the plugin menu of PyMOL.
def pymod_launcher(app, pymod_plugin_name, pymod_version, pymod_revision):
    global pymod
    pymod = PyMod(app, pymod_plugin_name, pymod_version, pymod_revision)
    pymod.start_new_project()


###################################################################################################
# Class that creates the plugin structure.                                                        #
###################################################################################################

class PyMod:
    """
    Class to represent the PyMod plugin.
    """

    ###############################################################################################
    # STARTUP OF THE PLUGIN AND STRUCTURE OF THE MAIN WINDOW.                                     #
    ###############################################################################################

    def __init__(self, app, pymod_plugin_name, pymod_version, pymod_revision):

        self.pymod_plugin_name = pymod_plugin_name
        self.pymod_version = pymod_version
        self.pymod_revision = pymod_revision
        # Builds the plugin main window.
        self.build_pymod_main_window(app)

        # This is the list where are going to be stored all the sequences displayed in the main
        # window represented as objects of the "PyMod_element" class.
        self.pymod_elements_list = []
        # An index that increases by one every time an element is added to the above list by using
        # the .add_element_to_pymod() method.
        self.unique_index = 0

        # -----
        # Prepare PyMod files and folders that will be created in the project main directory.
        # -----

        # PyMod directory. The main folder where all PyMod files (with the exception of the
        # configuration file) will be stored.
        self.pymod_directory_name = "pymod" # "pymod"
        self.pymod_temp_directory_name = "pymod_temp_directory"
        self.projects_directory_name = "projects"
        self.default_project_name = "new_pymod_project"
        self.external_tools_directory_name = "external_tools"
        self.data_directory_name = "data"
        self.blast_databases_directory_name = "blast_databases"
        self.blast_databases_directory_shortcut = os.path.join(self.data_directory_name, self.blast_databases_directory_name)
        self.temp_directory_name = "temp"

        # Structures.
        self.structures_directory = "structures"
        # A list of all PDB files loaded in PyMod by the user.
        self.pdb_list = []

        # Alignments.
        self.alignments_directory = "alignments"
        self.alignments_files_names = "alignment_" # or "alignment_temp"
        # Number of the alignments performed by the user.
        self.alignment_count = 0
        self.new_clusters_counter = 0

        # Images directory
        self.images_directory = "images"
        self.logo_image_counter = 0

        # Models.
        self.models_directory = "models"
        self.models_subdirectory = "model"
        # Attributes that will keep track of how many models the user builds in a PyMod session.
        self.performed_modeling_count = 0
        # This will keep track of how many multiple chains models the user builds in a PyMod session.
        self.multiple_chain_models_count = 0
        # This will contain ojects of the 'Modeling_session' class in order to build the 'Models'
        # submenu on the plugin's main menu.
        self.modeling_session_list = []
        # The maximum number of models that Modeler can produce at the same time.
        self.max_models_per_session = 1000
        self.multiple_chains_models_name = "MyMultiModel"

        # PSIPRED.
        self.psipred_directory = "psipred"

        # BLAST.
        self.similarity_searches_directory = "similarity_searches"
        self.temp_database_directory_name = "db_temp"
        self.blast_cluster_counter = 0

        # Gets the home directory of the user.
        self.home_directory = pmos.get_home_dir()
        self.current_project_name = None
        self.current_project_directory_full_path = None

        # Creates the preferences file in an hidden directory in the home directory.
        self.cfg_directory_path = os.path.join(self.home_directory,".pymod")
        self.cfg_file_name = "preferences.pkl"
        self.cfg_file_path = os.path.join(self.cfg_directory_path, self.cfg_file_name)

        # -----
        # Prepare PyMod tools.
        # -----
        self.pymod_tools = []

        # In order for the plugin to work, the name of the attribute containing a 'Tool' object must
        # must be the same of the first argument provided in the 'Tool' object constructor.

        # PyMod itself.
        self.pymod_plugin = pm_tool.Tool("pymod", self.pymod_plugin_name)
        self.pymod_plugin.initialize_parameters([pm_tool.Tool_directory("pymod_dir_path", "PyMod Directory", parameter_default_value = pmos.get_home_dir())])
        self.pymod_tools.append(self.pymod_plugin)

        # ClustalW.
        self.clustalw = pm_tool.Executable_tool("clustalw", "Clustal W", "local")
        self.clustalw.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.clustalw)

        # Clustal Omega.
        self.clustalo = pm_tool.Executable_tool("clustalo", "Clustal Omega")
        self.clustalo.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.clustalo)

        # MUSCLE.
        self.muscle = pm_tool.Executable_tool("muscle", "MUSCLE")
        self.muscle.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.muscle)

        # BLAST+ suite. Used to run PSI-BLAST and store BLAST sequence databases retrieved from
        # ftp://ftp.ncbi.nlm.nih.gov/blast/db/ .
        self.blast_plus = pm_tool.Executable_tool("blast_plus", "BLAST+ suite")
        self.blast_plus.initialize_parameters([pm_tool.Tool_exec_directory("exe_dir_path", "Executable Directory"),
                                               # A default directory where the database folders available for the
                                               # PSI-BLAST database selection are going to be located.
                                               pm_tool.Tool_directory("database_dir_path", "Database Directory")])
        self.pymod_tools.append(self.blast_plus)

        # PSIPRED.
        self.psipred = pm_tool.Executable_tool("psipred", "PSIPRED", "local")
        self.psipred.initialize_parameters([pm_tool.Tool_directory("exe_dir_path", "Executable Directory"),
                                            pm_tool.Tool_directory("data_dir_path", "Data Files Directory"),
                                            pm_tool.Tool_directory("database_dir_path", "BLAST Db Directory")])
        self.pymod_tools.append(self.psipred)

        # KSDSSP.
        self.ksdssp = pm_tool.Executable_tool("ksdssp", "KSDSSP")
        self.ksdssp.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.ksdssp)

        # MODELLER.
        self.modeller = pm_tool.Modeller_tool("modeller", "MODELLER")
        # Attempts to import MODELLER. If MODELLER can't be imported, its usage will be external to
        # the Python interpreter of PyMOL.
        self.import_modeller()
        # Then initializes the tool.
        self.modeller.initialize_parameters([pm_tool.Use_importable_modeller("use_importable_modeller", "Internal MODELLER"),
                                             pm_tool.Modeller_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.modeller)

        # Finish to initialize PyMod tools.
        for tool in self.pymod_tools:
            tool.show_message_method = self.show_error_message

        # -----
        # Prepares colors for PyMod and PyMOL.
        # -----
        # This is an index that will bw used to color each structure loaded into PyMod with a
        # different color taken from the list above. It will be used in "color_struct()".
        self.color_index = 0

        # Generates PSIPRED predictions colors for PyMOL.
        for c in pmdt.psipred_color_dict.keys():
            # Generates names like: 'pymod_psipred_8_H' (the name of the color with which residues
            # predicted in an helix with confidence score of 8 will be colored).
            color_name = "%s_%s_%s" % (pmdt.pymol_psipred_color_name, c[0], c[1])
            cmd.set_color(color_name, pmdt.psipred_color_dict[c])

        # Prepares CAMPO colors in PyMOL. Generates something like: 'pymod_campo_7'
        for c in pmdt.campo_color_dictionary.keys():
            color_name = "%s_%s" % (pmdt.pymol_campo_color_name, c)
            cmd.set_color(color_name, pmdt.campo_color_dictionary[c])

        for c in pmdt.dope_color_dict.keys():
            color_name = "%s_%s" % (pmdt.pymol_dope_color_name, c)
            cmd.set_color(color_name, pmdt.dope_color_dict[c])

        for c in pmdt.polarity_color_dictionary.keys():
            cmd.set_color(pmdt.pymol_polarity_color_name + c, pmdt.polarity_color_dictionary[c])

        # Builds the menu of the main window.
        self.make_main_menu()


    def start_new_project(self):
        """
        Starts up a new job.
        """
        self.initialize_project_data()
        # Cheks for PyMod configuration file.
        self.configuration_file_error = False

        # If it is not found, then treat this session as the first one and asks the user to input
        # the 'PyMod Directory' path before beginning the first PyMod job.
        if not os.path.isfile(self.cfg_file_path):
            self.show_first_time_usage_message()
            self.show_pymod_directory_selection_window()

        # The configuration file is found.
        else:
            try:
                # Check if there is 'pymod_temp_directory' left by the PyMod installer script.
                if not self.check_installer_script_temp_directory():
                    # Get values options for each PyMod tool and start a new PyMod job.
                    self.initialize_session()
                # If there is a 'pymod_temp_directory' (the installer script has been used before
                # this last PyMod session).
                else:
                    # The installer script was run before configuring PyMod for the first time (it
                    # left an empty configuratio file).
                    if self.check_empty_configuration_file():
                        self.show_first_time_usage_message()
                        self.show_pymod_directory_selection_window()
                    # The installer script was run after the first PyMod session (in order to
                    # install some missing tools).
                    else:
                        self.initialize_session()

            except Exception, e:
                self.show_configuration_file_error(e, "read")
                title = 'Configuration file repair'
                message = "Would you like to delete PyMod configuration file and build a new functional copy of it?"
                repair_choice = tkMessageBox.askyesno(title, message)
                self.configuration_file_error = True
                if repair_choice:
                    self.show_pymod_directory_selection_window()
                else:
                    self.main_window.destroy()


    def initialize_project_data(self):
        self.unique_index = 0
        self.pdb_list = []
        self.alignment_count = 0
        self.new_clusters_counter = 0
        self.logo_image_counter = 0
        self.performed_modeling_count = 0
        self.multiple_chain_models_count = 0
        self.modeling_session_list = []
        self.blast_cluster_counter = 0
        self.current_project_name = None
        self.current_project_directory_full_path = None
        self.color_index = 0
        self.pymod_elements_list = []


    def show_first_time_usage_message(self):
        title = "PyMod first session"
        message = "This is the first time you run PyMod. Please specify a folder inside which to build the 'PyMod Directory'. "
        message += "All PyMod files (such as its external tools executables, sequence databases and its project files) will be stored in this 'PyMod Directory' on your system."
        tkMessageBox.showinfo(title, message, parent=self.main_window)


    def get_parameters_from_configuration_file(self):
        """
        Updates the values of the PyMod Tools parameters according to the information in the main
        configuration file.
        """
        cfgfile = open(self.cfg_file_path, 'r')
        # Reads the pickle configuration file, where PyMod options are stored in a dictionary.
        pymod_config_data = pickle.load(cfgfile)
        for tool_name in pymod_config_data.keys():
            tool_object = self.get_tool_by_name(tool_name)
            for parameter_name in pymod_config_data[tool_name].keys():
                tool_object[parameter_name].value = pymod_config_data[tool_name][parameter_name]
        cfgfile.close()


    def check_pymod_directory(self):
        if os.path.isdir(self.pymod_plugin["pymod_dir_path"].get_value()):
            return True
        else:
            raise Exception("The project directory specified in PyMod configuration file ('%s') is missing. Please specify a new one." % (self.pymod_plugin["pymod_dir_path"].get_value()))


    def check_installer_script_temp_directory(self):
        """
        Checks a temporary directory built by the PyMod installer script is present. If, when
        starting up PyMod, the plugins detects this directory generated by the installer scritp, it
        will treat the session as the first one. After the first session a regular configuration
        file will be built.
        """
        if os.path.isdir(os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name)):
            return True
        else:
            return False


    def check_empty_configuration_file(self):
        """
        Checks if the PyMod configuration file is empty. The PyMod installer script generates an
        empty configuration file when it is run before ever running a first PyMod session.
        """
        if os.path.exists(self.cfg_file_path) and os.stat(self.cfg_file_path).st_size == 0:
            return True
        else:
            return False


    def get_tool_by_name(self, tool_name):
        for tool in self.pymod_tools:
            if tool.name == tool_name:
                return tool
        raise Exception("PyMod does not have a tool named '%s'" % (tool_name))


    def show_pymod_directory_selection_window(self):
        """
        Allows to select the 'PyMod Directory' on PyMod first run.
        """
        self.pymod_dir_window = pmgi.PyMod_base_window(self.main_window, "PyMod Directory")
        self.pymod_dir_window.geometry('-100+75')
        self.pymod_dir_window.config()
        self.pymod_dir_window.protocol("WM_DELETE_WINDOW", lambda: self.confirm_close(parent=self.pymod_dir_window))

        self.pymod_dir_window_label=Label(self.pymod_dir_window.main_frame, text= "Select a folder inside which to build the 'PyMod Directory'", **pmgi.label_style_0) # Select a folder (where you would like) to create the 'PyMod Directory'
        self.pymod_dir_window_label.pack(fill="x", pady=10, padx=10)

        self.pymod_dir_window_entry_frame = Frame(self.pymod_dir_window.main_frame, bg="black")
        self.pymod_dir_window_entry_frame.pack()

        home_directory = pmos.get_home_dir()
        self.pymod_dir_window_main_entry=Entry(self.pymod_dir_window_entry_frame, bg='white', width=30)
        self.pymod_dir_window_main_entry.insert(0, home_directory)
        self.pymod_dir_window_main_entry.pack(side="left")

        self.pymod_dir_window_browse_button=Button(self.pymod_dir_window_entry_frame, text="BROWSE",
        command=self.pymod_directory_browse_state, **pmgi.button_style_2)
        self.pymod_dir_window_browse_button.pack(side="left", pady=0, padx=5)

        self.pymod_dir_window_button_frame = Frame(self.pymod_dir_window.main_frame, bg="black")
        self.pymod_dir_window_button_frame.pack()

        self.pymod_dir_window_submit_button=Button(self.pymod_dir_window_button_frame, text="SUBMIT",
            command=self.pymod_directory_selection_state, **pmgi.button_style_1)
        self.pymod_dir_window_submit_button.pack(side="left", pady=10, padx=5)


    def pymod_directory_browse_state(self):
        current_path = self.pymod_dir_window_main_entry.get()
        # Lets users choose a new path.
        new_path = askdirectory(title = "Select a folder in which to build the 'PyMod Directory'",
                                initialdir=current_path, mustexist = True, parent = self.pymod_dir_window)
        # Updates the text in the Entry with the new path name.
        if new_path:
            self.pymod_dir_window_main_entry.delete(0, END)
            self.pymod_dir_window_main_entry.insert(0, new_path)


    def pymod_directory_selection_state(self):
        """
        This is called when the SUBMIT button on the "PyMod project" window is pressed.
        """
        try:
            # Check if the parent folder of the new PyMod directory exists.
            new_pymod_directory_parent = self.pymod_dir_window_main_entry.get()
            if not os.path.isdir(new_pymod_directory_parent):
                title = 'PyMod directory Error'
                message = "The path where you would like to create your 'PyMod Directory' does not exist on your system. Please select an existing path."
                self.show_error_message(title, message, parent_window=self.pymod_dir_window)
                return False

            # Check if a PyMod directory already exists in the parent folder.
            new_pymod_directory_path = os.path.join(new_pymod_directory_parent, self.pymod_directory_name)
            if os.path.exists(new_pymod_directory_path):
                title = 'PyMod directory Warning'
                message = "A folder named '%s' already exists on your system. Would you like to overwrite it and all its contents to build a new 'PyMod Directory'?" % (new_pymod_directory_path)
                overwrite = tkMessageBox.askyesno(title, message, parent=self.pymod_dir_window)
                # If the directory (or a file with the same name) already exists, try to remove it.
                if overwrite:
                    pmos.pymod_rm(new_pymod_directory_path)
                else:
                    return False

            # Check if the configuration file directory exist. If it does not exist, build it.
            if not os.path.exists(self.cfg_directory_path):
                os.mkdir(self.cfg_directory_path)

            # Builds the new PyMod directory and with its projects folder.
            os.mkdir(new_pymod_directory_path)
            os.mkdir(os.path.join(new_pymod_directory_path, self.projects_directory_name))

            # Chekc ff the Installer Bundle install_all.py script was run. If so, a
            # 'pymod_temp_directory' with external tools and data files has been built.
            pymod_first_run = False
            if self.check_installer_script_temp_directory() and self.check_empty_configuration_file():
                # This will move the content of the 'pymod_temp_directory' to the new PyMod
                # directory.
                temp_installation_path = os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name)
                for installation_directory in os.listdir(temp_installation_path):
                    source_directory = os.path.join(temp_installation_path,installation_directory)
                    target_directory = os.path.join(new_pymod_directory_path, installation_directory)
                    shutil.move(source_directory, target_directory)
                pymod_first_run = True
            # Builds empty external tools and data directories.
            else:
                os.mkdir(os.path.join(new_pymod_directory_path, self.external_tools_directory_name))
                os.mkdir(os.path.join(new_pymod_directory_path, self.data_directory_name))

            # Updates the configuration file.
            cfgfile = open(self.cfg_file_path, 'w')
            pymod_config_data = {}
            for tool in self.pymod_tools:
                new_tool_parameters = {}
                for parameter in tool.parameters:
                    # NOTE: new_tool_parameters.update({parameter.name: parameter.get_first_session_value()})
                    new_tool_parameters.update({parameter.name: parameter.get_starting_value()})
                pymod_config_data.update({tool.name: new_tool_parameters})
            # Define the paths of the 'PyMod Directory' and the 'BLAST Database Directory'.
            pymod_config_data["pymod"]["pymod_dir_path"] = new_pymod_directory_path
            pymod_config_data["blast_plus"]["database_dir_path"] = os.path.join(new_pymod_directory_path, self.data_directory_name, self.blast_databases_directory_name)

            # If an empty configuration file was built by the PyMod installer script, update it.
            if self.check_installer_script_temp_directory():
                for tool in pymod_config_data.keys():
                    for parameter_name in pymod_config_data[tool].keys():
                        if pymod_config_data[tool][parameter_name] == "":
                            if tool in ("clustalw", "muscle", "clustalo", "ksdssp"):
                                exe_file_name = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "bin", tool)
                                if pymod_first_run:
                                    exe_file_name = pmos.get_exe_file_name(exe_file_name)
                                pymod_config_data[tool][parameter_name] = exe_file_name
                            elif tool == "blast_plus":
                                if parameter_name == "exe_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "bin")
                            elif tool == "psipred":
                                if parameter_name == "exe_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "bin")
                                elif parameter_name == "data_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "data")
                                elif parameter_name == "database_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.data_directory_name, self.blast_databases_directory_name, "swissprot")
                # Finally remove pymod temp directory in the configuration directory.
                shutil.rmtree(os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name))

            pickle.dump(pymod_config_data, cfgfile)
            cfgfile.close()

        except Exception,e:
            title = "PyMod Directory Error"
            message = "Unable to write the PyMod configuration directory '%s' because of the following error: %s." % (self.cfg_directory_path, e)
            self.show_error_message(title, message, parent_window=self.pymod_dir_window)
            return False

        # Begin a new PyMod job.
        self.pymod_dir_window.destroy()
        self.get_parameters_from_configuration_file()
        # tkMessageBox.showinfo("PyMod first run", "You are about to begin your first PyMod session.", parent=self.main_window)
        self.new_job_state()


    def show_configuration_file_error(self, error, mode):
        action = None
        if mode == "read":
            action = "reading"
        elif mode == "write":
            action = "writing"
        title = "Configuration file error"
        message = "There was an error while %s the PyMod configuration file (located at '%s'). The error is: '%s'." % (action, self.cfg_file_path, error)
        self.show_error_message(title, message)


    def build_pymod_main_window(self, app):
        """
        Builds the structure of the PyMod main window.
        """
        self.in_a_job = False
        self.parent = app.root
        self.main_window = Toplevel(self.parent)
        self.main_window.title(self.pymod_plugin_name)
        self.main_window.resizable(1,1)
        self.main_window.geometry('800x320')

        # Asks confirmation when the main window is closed by the user.
        self.main_window.protocol("WM_DELETE_WINDOW", self.confirm_close)
        # Hides PyMod main window, it will be displayed once the user begins a new project by
        # inserting the project name in the 'new project' window.
        self.main_window.withdraw()

        # Creates a scrolled frame in the main window.
        scroll_frame = Pmw.ScrolledFrame(self.main_window, borderframe = 0, usehullsize = 1,
                                    horizflex = 'elastic', vertflex = 'elastic', hull_borderwidth = 0 )
        scroll_frame.configure(frame_background = 'black')
        scroll_frame.pack(fill = 'both', expand = 1)
        frame_main = scroll_frame.interior()
        frame_main.config()

        # Creates a paned widget in the scrolled frame 'frame_main'
        self.panes = Pmw.PanedWidget(frame_main, orient = 'horizontal', hull_borderwidth = 0)

        # Adds the left pane (where the name of the sequences are) and the right pane (where the
        # sequences are displayed)
        self.panes.add('left', size = 0.2)
        self.panes.add('right', size = 0.8)
        self.panes.pack(fill = 'both', expand = 1)

        # This method is defined later
        self.create_main_window_panes()

        # Creates the bottom frame that display the name of the sequence
        self.sequence_name_bar = Pmw.MessageBar(self.main_window,
            entry_width = 10,
            entry_relief='groove',
            entry_bg = 'black',
            labelpos = 'w',
            label_text = 'Sequence:',
            label_fg = 'white',
            label_background='black')
        self.sequence_name_bar.pack(side=LEFT, fill = 'x', expand = 1)

        # Creates the bottom frame that display the number and the name of the residue
        self.residue_bar = Pmw.MessageBar(self.main_window,
                entry_width = 50, # This could be modified.
                entry_relief='groove',
                labelpos = 'w',
                label_text = 'Position:',
                label_fg = 'white',
                label_background='black')
        self.residue_bar.pack(side=RIGHT)

        # Variables needed to make Pmw dialogs work on Ubuntu 14.04+.
        self.pmw_dialog_wait = True
        self.pmw_dialog_val = None


    def initialize_session(self):
        self.get_parameters_from_configuration_file()
        self.check_pymod_directory()
        self.new_job_state()


    def show_new_job_window(self):
        """
        Builds a window that let users choose the name of the new projects direcotory at the
        beginning of a PyMod session.
        """
        pass
        # self.in_a_job = True
        # self.new_dir_window = pmgi.PyMod_base_window(self.main_window, "New PyMod Project")
        # self.new_dir_window.geometry('-100+75')
        # self.new_dir_window.resizable(0,0)
        # self.new_dir_window.config()
        # self.new_dir_window.protocol("WM_DELETE_WINDOW", lambda: self.confirm_close(parent=self.new_dir_window))
        #
        # self.new_dir_window_label=Label(self.new_dir_window.main_frame, text= "Enter the name of your new PyMod project", **pmgi.label_style_0) # Create a new directory inside your PyMod projects folder
        # self.new_dir_window_label.pack(fill="x", pady=10, padx=10)
        #
        # self.new_dir_window_main_entry=Entry(self.new_dir_window.main_frame, bg='white', width=18)
        # self.new_dir_window_main_entry.insert(0, "new_pymod_project")
        # self.new_dir_window_main_entry.pack()
        #
        # self.new_dir_window_main_submit=Button(self.new_dir_window.main_frame, text="SUBMIT",
        #     command=self.new_job_state, **pmgi.button_style_1)
        # self.new_dir_window_main_submit.pack(pady=10)


    def new_job_state(self):
        """
        This is called when the SUBMIT button on the "New Job" window is pressed.
        """

        new_project_dir_name = self.default_project_name
        try:
            #------------------------------------
            # Writes the new project directory. -
            #------------------------------------
            pymod_dir_path = self.pymod_plugin["pymod_dir_path"].get_value()
            pymod_projects_dir_path = os.path.join(pymod_dir_path, self.projects_directory_name)
            new_project_dir_path = os.path.join(pymod_projects_dir_path, new_project_dir_name)
            # If the projects directory is not present, built it.
            pmos.pymod_mkdir(pymod_projects_dir_path)
            # If a directory with the same name of the new project directory exists, delete it.
            if os.path.isdir(new_project_dir_path):
                self.clean_project_directory(new_project_dir_path)
            else:
                pmos.pymod_mkdir(new_project_dir_path)

            #--------------------------------------------------------------------------------
            # Initializes a session and shows the main window, which was previously hidden. -
            #--------------------------------------------------------------------------------
            self.current_pymod_directory = pymod_dir_path
            self.current_project_name = new_project_dir_name
            self.current_project_directory_full_path = os.path.join(self.current_pymod_directory, self.projects_directory_name, self.current_project_name)
            os.chdir(self.current_project_directory_full_path)
            self.create_project_subdirectories()
            self.main_window.deiconify()
            self.launch_default()
            self.gridder()

        except Exception, e:
            message = "Unable to write directory '%s' because of the following error: %s." % (new_project_dir_name, e)
            self.show_error_message("Initialization error", message)
            return False


    def clean_project_directory(self, new_project_dir_path):
        """
        Remove all the files in the previously used directory.
        """
        for the_file in os.listdir(new_project_dir_path):
            file_path = os.path.join(new_project_dir_path, the_file)
            if os.path.isfile(file_path):
                os.unlink(file_path)
        self.remove_project_subdirectories(new_project_dir_path)


    def launch_default(self):
        """
        For development only. A the 'open_sequence_file', 'open_pdb_file' and
        'build_cluster_from_alignment_file' methods to import sequences when PyMod starts.
        """
        pass


    def create_main_window_panes(self):
        """
        This method allows to create the panes containing the names and sequences to display.
        """
        # Creates a scrolled frame inside the RIGHT pane of the paned frame
        self.rightpan = Pmw.ScrolledFrame(self.panes.pane('right'),
            hull_bg='black', frame_bg='black', usehullsize = 0, borderframe = 0,
            hscrollmode='static', hull_borderwidth = 0, clipper_bg='black')
        self.rightpan.pack(fill = 'both', expand = 1)

        # Creates a scrolled frame inside the LEFT pane of the paned frame
        self.leftpan = Pmw.ScrolledFrame(self.panes.pane('left'),
            hull_bg='black', frame_bg = 'black', hull_borderwidth = 0, usehullsize = 0,
            borderframe = 0, vscrollmode=NONE, hscrollmode='static', clipper_bg='black' )
        self.leftpan.pack(fill = 'both', expand = 1)

        # Allows to scroll both RIGHT and LEFT scrolled frame using only one ScrollBar.
        def vertview(*args):
            self.rightpan.yview(*args)
            self.leftpan.yview(*args)

        self.rightpan.configure(vertscrollbar_command = vertview)


    def make_main_menu(self):
        """
        This method is called at the beginning of the constructor in order to build the main menu of
        the main window.
        """
        self.menubar = Menu(self.main_window)

        # ---
        # "File" menu.
        # ---
        self.filemenu = Menu(self.menubar, tearoff = 0)
        self.sequence_menu = Menu(self.filemenu, tearoff = 0)
        self.filemenu.add_cascade(label = "Sequences", menu = self.sequence_menu)
        self.sequence_menu.add_command(label = "Open from File", command = self.open_file_from_the_main_menu)
        self.sequence_menu.add_command(label = "Add Raw Sequence", command = self.raw_seq_input)
        self.sequence_menu.add_command(label = "Import PyMOL Objects", command = self.import_selections)
        self.sequence_menu.add_separator()
        self.sequence_menu.add_command(label = "Save All", command = self.save_all_files_from_main_menu)
        self.filemenu.add_separator()

        # Submenu to open alignments.
        self.alignment_files_menu = Menu(self.filemenu, tearoff = 0)
        self.filemenu.add_cascade(label = "Alignment", menu = self.alignment_files_menu)
        self.alignment_files_menu.add_command(label = "Open from File", command = self.open_alignment_from_main_menu)
        self.filemenu.add_separator()

        # Projects submenu.
        self.projects_menu = Menu(self.filemenu, tearoff = 0)
        self.filemenu.add_cascade(label = "Projects", menu = self.projects_menu)
        self.projects_menu.add_command(label = "New Project", command = self.begin_new_project_from_main_menu)
        self.projects_menu.add_command(label = "Save Project as", command = self.save_new_project_from_main_menu)
        self.projects_menu.add_command(label = "Open Project", command = self.open_new_project_from_main_menu)
        self.filemenu.add_separator()

        self.filemenu.add_command(label = "Exit", command = self.confirm_close)
        self.menubar.add_cascade(label = "File", menu = self.filemenu)

        # ---
        # "Tools" menu.
        # ---
        self.tools_menu = Menu(self.menubar, tearoff = 0)

        # Sequence alignment tools.
        self.sequence_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Sequence Alignment", menu = self.sequence_alignment_menu)
        self.sequence_alignment_menu.add_command(label = "ClustalW",
            command = lambda program="clustalw": self.launch_regular_alignment_from_the_main_menu(program))
        self.sequence_alignment_menu.add_command(label = "Clustal Omega",
            command = lambda program="clustalo": self.launch_regular_alignment_from_the_main_menu(program))
        self.sequence_alignment_menu.add_command(label = "MUSCLE",
            command = lambda program="muscle": self.launch_regular_alignment_from_the_main_menu(program))
        self.sequence_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)",
            command = lambda program="salign-seq": self.launch_regular_alignment_from_the_main_menu(program))

        # Profile alignment tools.
        self.profile_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Profile Alignment", menu = self.profile_alignment_menu)
        self.profile_alignment_menu.add_command(label = "ClustalW",
            command = lambda program="clustalw": self.launch_profile_alignment_from_the_main_menu(program))
        self.profile_alignment_menu.add_command(label = "Clustal Omega",
            command = lambda program="clustalo": self.launch_profile_alignment_from_the_main_menu(program))
        self.profile_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)",
            command = lambda program="salign-seq": self.launch_profile_alignment_from_the_main_menu(program))

        # Structural alignment tools.
        self.structural_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Structural Alignment", menu = self.structural_alignment_menu)
        self.structural_alignment_menu.add_command(label = "Superpose", command = self.superpose)
        self.structural_alignment_menu.add_command(label = "CE Alignment",
            command = lambda program="ce": self.launch_regular_alignment_from_the_main_menu(program))
        self.structural_alignment_menu.add_command(label = "SALIGN (Structure Alignment)",
            command = lambda program="salign-str": self.launch_regular_alignment_from_the_main_menu(program))

        # Database search for homologous sequences.
        self.database_search_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Database Search", menu = self.database_search_menu)
        self.database_search_menu.add_command(label = "BLAST", command = self.launch_ncbiblast)
        self.database_search_menu.add_command(label = "PSI-BLAST", command = self.launch_psiblast)

        # Structural analysis.
        self.structural_analysis_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Structural Analysis", menu = self.structural_analysis_menu)
        self.structural_analysis_menu.add_command(label = "Ramachandran plot", command = self.ramachandran_plot)
        self.structural_analysis_menu.add_command(label = "Assess with DOPE", command = self.dope_from_main_menu)
        self.structural_analysis_menu.add_command(label = "PSIPRED", command = self.launch_psipred_from_main_menu)

        # Homology modeling (MODELLER).
        self.homology_modeling_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Homology Modeling", menu = self.homology_modeling_menu)
        self.homology_modeling_menu.add_command(label = "MODELLER", command = self.launch_modeller_from_main_menu)

        # Options.
        self.tools_menu.add_separator()
        self.tools_menu.add_command(label = "Options", command = self.show_pymod_options_window)

        # Adds the "Tools" menu to the main menu
        self.menubar.add_cascade(label = "Tools", menu = self.tools_menu)

        # ---
        # "Alignments" menu.
        # ---
        self.alignments_menu = Menu(self.menubar, tearoff = 0)
        # When the plugin is started there are no alignments.
        self.build_alignment_submenu()
        # Adds the "Alignments" menu to the main menu
        self.menubar.add_cascade(label = "Alignments", menu = self.alignments_menu)
        self.define_alignment_menu_structure()

        # ---
        # "Models" menu.
        # ---
        self.models_menu = Menu(self.menubar, tearoff = 0)
        # When the plugin is started there are no models.
        self.build_models_submenu()
        # Adds the "Alignments" menu to the main menu
        self.menubar.add_cascade(label = "Models", menu = self.models_menu)
        # self.define_models_menu_structure()

        # ---
        # "Selection" menu.
        # ---
        self.main_selection_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Selection", menu = self.main_selection_menu)
        # When the plugin is started there are no models.
        self.main_selection_menu.add_command(label = "Select All", command=self.select_all_from_main_menu)
        self.main_selection_menu.add_command(label = "Deselect All", command=self.deselect_all_from_main_menu)
        # Structures selection submenu.
        self.selection_structures_menu = Menu(self.main_selection_menu,tearoff=0)
        self.selection_structures_menu.add_command(label="Show All in PyMOL",command=self.show_all_structures_from_main_menu)
        self.selection_structures_menu.add_command(label="Hide All in PyMOL",command=self.hide_all_structures_from_main_menu)
        self.selection_structures_menu.add_separator()
        self.selection_structures_menu.add_command(label="Select All",command=self.select_all_structures_from_main_menu)
        self.selection_structures_menu.add_command(label="Deselect All",command=self.deselect_all_structures_from_main_menu)
        self.main_selection_menu.add_cascade(menu=self.selection_structures_menu, label="Structures")
        # Clusters selection submenu.
        self.selection_clusters_menu = Menu(self.main_selection_menu,tearoff=0)
        self.selection_clusters_menu.add_command(label="Expand All",command=self.expand_all_clusters_from_main_menu)
        self.selection_clusters_menu.add_command(label="Collapse All",command=self.collapse_all_clusters_from_main_menu)
        self.main_selection_menu.add_cascade(menu=self.selection_clusters_menu, label="Clusters")

        # ---
        # "Display" menu.
        # ---
        self.display_menu = Menu(self.menubar, tearoff = 0)

        # Color menu.
        self.main_color_menu = Menu(self.display_menu, tearoff = 0)
        self.main_color_menu.add_command(label = "By Regular Color Scheme", command=lambda: self.color_selection("all", None, "regular"))
        # Residues.
        self.main_residues_colors_menu = Menu(self.main_color_menu,tearoff=0)
        self.main_residues_colors_menu.add_command(label="Polarity",command=lambda: self.color_selection("all", None, "residue"))
        self.main_color_menu.add_cascade(menu=self.main_residues_colors_menu, label="By residue properties")
        # Secondary structure.
        self.main_color_menu.add_command(label="Secondary Structure",command=lambda: self.color_selection("all", None, "secondary-auto"))
        self.display_menu.add_cascade(menu=self.main_color_menu, label="Color all Sequences")

        # Font size menu.
        self.menu_sequence_font_size = StringVar()
        self.default_font_size = 12 # "14"
        self.menu_sequence_font_size.set(self.default_font_size)
        self.font_menu = Menu(self.display_menu, tearoff = 0)
        self.font_menu.add_radiobutton(label="6",value="6",variable=self.menu_sequence_font_size, command=self.gridder)
        self.font_menu.add_radiobutton(label="8",value="8",variable=self.menu_sequence_font_size, command=self.gridder)
        self.font_menu.add_radiobutton(label="10",value="10",variable=self.menu_sequence_font_size, command=self.gridder)
        self.font_menu.add_radiobutton(label="12",value="12",variable=self.menu_sequence_font_size, command=self.gridder)
        self.font_menu.add_radiobutton(label="14",value="14",variable=self.menu_sequence_font_size, command=self.gridder)
        self.font_menu.add_radiobutton(label="16",value="16",variable=self.menu_sequence_font_size, command=self.gridder)
        self.font_menu.add_radiobutton(label="18",value="18",variable=self.menu_sequence_font_size, command=self.gridder)
        self.display_menu.add_cascade(label = "Font size", menu = self.font_menu)

        # Adds the "Display" menu to the main menu.
        self.menubar.add_cascade(label = "Display", menu = self.display_menu)

        # ---
        # "Help" menu.
        # ---
        self.help_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Help", menu = self.help_menu)
        self.help_menu.add_command(label = "Online Documentation", command = self.open_online_documentation)
        self.help_menu.add_command(label = "About", command = self.show_about_dialog)
        self.help_menu.add_separator()
        self.help_menu.add_command(label = "Check for PyMod Updates", command = self.launch_pymod_update)

        self.main_window.config(menu = self.menubar)


    def define_alignment_menu_structure(self):
        """
        Define a series of tuples that are going to be used in order to build Alignments menus with
        different options according to the alignment algorithm used to build the alignment.
        """
        self.can_show_rmsd_matrix = ("ce","salign-str")
        self.can_show_guide_tree = ("clustalw","clustalo")
        self.can_show_dendrogram = ("salign-seq","salign-str")
        self.can_use_scr_find = ("ce","salign-str")


    def show_about_dialog(self):
        Pmw.aboutversion(self.pymod_version + "." + self.pymod_revision)
        Pmw.aboutcopyright('Copyright (C): 2016 Giacomo Janson, Chengxin Zhang,\nAlessandro Paiardini')
        Pmw.aboutcontact(
            'For information on PyMod %s visit:\n' % (self.pymod_version) +
            '  http://schubert.bio.uniroma1.it/pymod/documentation.html\n\n' +
            'Or send us an email at:\n' +
            '  giacomo.janson@uniroma1.it'
        )
        self.about = Pmw.AboutDialog(self.main_window, applicationname = self.pymod_plugin_name)
        self.about.show()


    def open_online_documentation(self):
        webbrowser.open("http://schubert.bio.uniroma1.it/pymod/documentation.html")


    def launch_pymod_update(self):
        # Gets the latest release number from network.
        try:
            update_found = pmup.check_for_updates(self.pymod_version, self.pymod_revision)
        except Exception, e:
            self.show_error_message("Connection Error", "Can not obtain the latest PyMod version number beacause of the following error: '%s'" % e)
            return False

        if not update_found:
            self.show_warning_message("Update Canceled", "Your PyMod version (%s.%s) is already up to date." % (self.pymod_version, self.pymod_revision))
            return False

        # Ask for update confirmation.
        title = "Update PyMod?"
        message = "Would you like to update your current PyMod version (%s.%s) to the latest stable one available online (%s)? You will need to restart PyMOL in order to use the new version." % (self.pymod_version, self.pymod_revision, update_found)
        answer = tkMessageBox.askyesno(title, message, parent=self.main_window)
        if not answer:
            return False

        # Fetches the latest stable version files of PyMod.
        try:
            plugin_zipfile_temp_name = pmup.fetch_plugin_zipfile()
        except Exception, e:
            self.show_error_message("Connection Error", "Can not fetch the latest PyMod files beacause of the following error: '%s'" % e)
            return False

        if not plugin_zipfile_temp_name:
            return False

        # Installs the new PyMod version.
        pymod_plugin_dir = os.path.dirname(os.path.dirname(__file__))
        update_results = pmup.update_pymod(plugin_zipfile_temp_name, pymod_plugin_dir)
        if update_results[0]:
            self.show_info_message("Update Successful", "PyMod has been updated to version %s. Please restart PyMOL in order to use the updated PyMod version." % (update_found))
        else:
            self.show_error_message("Update Failed", update_results[1])


    def show_popup_message(self, popup_type="warning", title_to_show="ALLERT", message_to_show="THIS IS AN ALLERT MESSAGE", parent_window=None, refresh=True):
        """
        Displays error or warning messages and refreshes the sequence window.
        """
        # show_error_message
        # show_warning_message
        if parent_window == None:
            parent_window = self.main_window

        if popup_type == "error":
            tkMessageBox.showerror(title_to_show, message_to_show, parent=parent_window)
        elif popup_type == "info":
            tkMessageBox.showinfo(title_to_show, message_to_show, parent=parent_window)
        elif popup_type == "warning":
            tkMessageBox.showwarning(title_to_show, message_to_show, parent=parent_window)

        if refresh:
            self.gridder()


    def show_info_message(self, title_to_show,message_to_show,parent_window=None,refresh=True):
        self.show_popup_message("info",title_to_show,message_to_show,parent_window,refresh)


    def show_warning_message(self, title_to_show,message_to_show,parent_window=None,refresh=True):
        self.show_popup_message("warning",title_to_show,message_to_show,parent_window,refresh)


    def show_error_message(self, title_to_show,message_to_show,parent_window=None,refresh=True):
        self.show_popup_message("error",title_to_show,message_to_show,parent_window,refresh)


    def general_error(self,e=''):
        title = "Unknown Error"
        message = "PyMod has experienced an unknown error:\n"+str(e)
        self.show_error_message(title,message)


    def confirm_close(self, parent=None):
        """
        Asks confirmation when the main window is closed by the user.
        """
        parent_window = None
        if parent == None:
            parent_window = self.main_window
        else:
            parent_window = parent
        answer = tkMessageBox.askyesno(message="Are you really sure you want to exit PyMod?", title="Exit PyMod?", parent=parent_window)
        if answer:
            self.main_window.destroy()


    def check_general_input(self, pymod_tool_window):
        """
        Checks if valid input has been supplied by users in PyMod tools windows.
        """
        for widget in pymod_tool_window.get_widgets_to_validate():
            if widget.getvalue() in ("","-"):
                title = "Input Error"
                message = "Please fill in the '%s' option with valid input." % (widget.component("label").cget("text"))
                self.show_error_message(title, message, parent_window=pymod_tool_window, refresh=False)
                return False
        return True


    def execute_subprocess(self, commandline, new_stdout = subprocess.PIPE, new_stderr = subprocess.PIPE, new_shell = (sys.platform!="win32"), print_stdinfo = False, executing_modeller=False):
        if print_stdinfo:
            print "Executing the following command:", commandline
        if not executing_modeller:
            try:
                subp = subprocess.Popen(commandline, stdout= new_stdout, stderr= new_stderr, shell= new_shell)
                out_std, err_std = subp.communicate()
                returncode = subp.returncode
                if returncode != 0:
                    raise Exception("Subprocess returned non-zero return code...")
                if print_stdinfo:
                    print "Stdout:", out_std
            except Exception, e:
                if print_stdinfo:
                    print "Exception:", e
                    print "Stderr:", err_std
                raise Exception("An error occurred while running the child process.")
        # Official PyMOL builds on Mac OS will crash if executing MODELLER through using the
        # 'subprocess' module. For this reason, the 'os' module will be used instead.
        else:
            os.system(commandline)

    def work_in_progress(self):
        raise Exception("Work in progress...")


    ###############################################################################################
    # PROGRAMS PATH AND PROJECTS MANAGMENT.                                                       #
    ###############################################################################################

    #################################################################
    # Import modules.                                               #
    #################################################################

    def import_modeller(self):

        # First check if systemwide MODELLER can be imported in PyMOL.
        if not pmos.check_importable_modeller():
            # If MODELLER can't be immediately imported, try to find the modlib directory import it.
            modeller_path = None
            if hasattr(self, "modeller"):
                modeller_path = self.modeller.get_exe_file_path()
            modeller_lib_path = pmos.find_modlib_path(modeller_path)
            if modeller_lib_path:
                sys.path.append(modeller_lib_path)

        # After having searched for 'modlib', try to actually import MODELLER.
        global importable_modeller
        try:
            global modeller
            global complete_pdb
            import modeller
            import modeller.automodel
            from modeller.scripts import complete_pdb
            importable_modeller = True
        except Exception, e:
            print e
            importable_modeller = False

        # Updates MODELLER tool status.
        if hasattr(self, "modeller"):
            self.modeller.importable_modeller = importable_modeller


    #################################################################
    # Methods used to build the widgets diplayed in the PyMod       #
    # Options window.                                               #
    #################################################################

    def show_pymod_options_window(self):
        """
        Builds a window that allows to modify some PyMod options.
        """
        # Builds the options window.
        self.pymod_options_window = pmgi.PyMod_tool_window(self.main_window,
            title = "PyMod Options",
            upper_frame_title = "Here you can modify options for PyMod",
            submit_command = lambda: self.set_pymod_options_state(),
            with_frame=True)
        self.pymod_options_window.resizable(1,1)
        self.pymod_options_window.geometry('580x600') # '500x600'

        # This list will be populated inside "build_tool_options_frame()".
        for single_tool in self.pymod_tools:
            single_tool.display_options(self.pymod_options_window.midframe)
            # If the tool list of parameter widgets has some alignable widgets, adds them to the
            # option window list.
            if single_tool.parameters_widgets_list != []:
                for w in single_tool.parameters_widgets_list:
                    self.pymod_options_window.add_widget_to_align(w)
        self.pymod_options_window.align_widgets(25)


    def set_pymod_options_state(self):
        """
        This function is called when the SUBMIT button is pressed in the PyMod options window.
        """
        old_projects_dir = self.pymod_plugin["pymod_dir_path"].get_value()
        new_projects_dir = self.pymod_plugin["pymod_dir_path"].get_value_from_gui()
        if not os.path.isdir(new_projects_dir):
            title = "Configuration Error"
            message = "The PyMod Projects Directory you specified ('%s') does not exist on your system. Please choose an existing directory." % (new_projects_dir)
            self.show_error_message(title, message, parent_window=self.pymod_options_window)
            return False

        # Saves the changes to PyMod configuration file.
        cfgfile = open(self.cfg_file_path, 'w')
        pymod_config_data = {}
        for tool in self.pymod_tools:
            new_tool_parameters = {}
            for parameter in tool.parameters:
                new_tool_parameters.update({parameter.name: parameter.get_value_from_gui()})
            new_tool_dict = {tool.name: new_tool_parameters}
            pymod_config_data.update(new_tool_dict)
        pickle.dump(pymod_config_data, cfgfile)
        cfgfile.close()

        # Then updates the values of the parameters of the tools contained in "self.pymod_tools"
        # so that they can be used in the current PyMod session.
        try:
            # Prevents the user from changing the project directory during a session.
            self.get_parameters_from_configuration_file()
            if old_projects_dir != new_projects_dir:
                title = "Configuration Updated"
                message = "You changed PyMod projects directory, the new directory will be used the next time you launch PyMod."
                self.show_warning_message(title, message, parent_window=self.pymod_options_window)
            self.pymod_options_window.destroy()

        except Exception,e:
            self.show_configuration_file_error(e, "read")
            self.main_window.destroy()


    #################################################################
    # Builds subdirectories in the current project directory.       #
    #################################################################

    def create_subdirectory(self, subdir):
        try:
            os.mkdir(subdir)
        except:
            pass

    def create_alignments_directory(self):
        self.create_subdirectory(self.alignments_directory)

    def create_images_directory(self):
        self.create_subdirectory(self.images_directory)

    def create_models_directory(self):
        self.create_subdirectory(self.models_directory)

    def create_model_subdirectory(self, model_subdirectory):
        self.create_subdirectory(model_subdirectory)

    def create_structures_directory(self):
        self.create_subdirectory(self.structures_directory)

    def create_psipred_directory(self):
        self.create_subdirectory(self.psipred_directory)

    def create_similarity_searches_directory(self):
        self.create_subdirectory(self.similarity_searches_directory)

    def create_temp_directory(self):
        self.create_subdirectory(self.temp_directory_name)


    def create_project_subdirectories(self):
        self.create_alignments_directory()
        self.create_images_directory()
        self.create_models_directory()
        self.create_structures_directory()
        self.create_psipred_directory()
        self.create_temp_directory()
        self.create_similarity_searches_directory()


    def remove_project_subdirectories(self, new_dir_name):
        """
        Removes the previously used subdirectories and all their content when users decide to
        overwrite an existing project's directory.
        """
        dirs_to_remove = (self.structures_directory, self.models_directory, self.alignments_directory, self.psipred_directory, self.similarity_searches_directory, self.images_directory, self.temp_directory_name)
        for single_dir in dirs_to_remove:
            dir_to_remove_path = os.path.join(new_dir_name, single_dir)
            if os.path.isdir(dir_to_remove_path):
                shutil.rmtree(dir_to_remove_path)


    #################################################################
    # Projects managment.                                           #
    #################################################################

    def begin_new_project_from_main_menu(self):
        answer = tkMessageBox.askyesno(message="Are you really sure you want to begin a new PyMod project? If you do not save your current project, its data will be permanently lost.", title="Begin New Project?", parent=self.main_window)
        if not answer:
            return None
        self.start_new_project()


    ##################
    # Save projects. #
    ##################

    def save_new_project_from_main_menu(self):
        save_project_full_path = asksaveasfilename(defaultextension="", filetypes=pmdt.pymod_session_extension, parent=self.main_window)
        if save_project_full_path == "":
            return None
        self.save_pymod_project(save_project_full_path)


    def save_pymod_project(self, project_arc_full_path):
        try:
            # Temporary directory which will become zipped.
            project_temp_dir_name = "pymod_project_temp_dir"
            project_temp_dir_path = os.path.join(self.current_pymod_directory, project_temp_dir_name)
            # Define the file path of the PyMod project archive file.
            project_arc_dir_path = os.path.dirname(project_arc_full_path)
            project_arc_name = "%s.pmse" % os.path.splitext(os.path.basename(project_arc_full_path))[0]
            project_name = self.default_project_name
        except:
            title = "Save Project Error"
            message = "Could not save the project."
            self.show_error_message(title, message)
            return None

        try:
            # Builds a temporary directory in which to store project files.
            pmos.pymod_mkdir(project_temp_dir_path, overwrite_dirs=True)

            # Saves a pickle file with the information about the PyMod project. This will remove
            # Tkinter objects stored in the 'PyMod' object, because they can't be pickled.
            project_pickled_file = open(os.path.join(project_temp_dir_path, "%s.pkl" % project_name), "w")
            PyMod_Pickler(project_pickled_file).dump(self)
            project_pickled_file.close()

            # Saves a PyMOL session.
            cmd.save(os.path.join(project_temp_dir_path, "%s.pse" % project_name))

            # Copies the current project files in the directory.
            src = self.current_project_directory_full_path
            dst = os.path.join(project_temp_dir_path, os.path.basename(self.current_project_directory_full_path))
            shutil.copytree(src, dst)

            # Builds a .zip file of the directory.
            src = project_temp_dir_path
            zpf = os.path.join(project_arc_dir_path, project_arc_name)
            pmos.zip_directory(src, zpf, use_dirname_root=False)
            shutil.rmtree(project_temp_dir_path)

        except:
            if os.path.isdir(project_temp_dir_path):
                shutil.rmtree(project_temp_dir_path)
            title = "Save Project Error"
            message = "Could not save the project file to path: %s" % (os.path.join(project_arc_dir_path, project_arc_name))
            self.show_error_message(title, message)


    ##################
    # Load projects. #
    ##################

    def open_new_project_from_main_menu(self):
        project_archive_file_path = askopenfilename(filetypes=pmdt.pymod_session_extension, multiple=False, parent=self.main_window)
        if project_archive_file_path == "":
            return None
        if not os.path.isfile(project_archive_file_path):
            return None
        self.open_pymod_project(project_archive_file_path)


    def open_pymod_project(self, project_archive_file_path):
        # If some errors happen here, continue the PyMod session.
        try:
            # Builds a temporary directory in which to store project files.
            project_name = self.default_project_name
            project_temp_dir_name = "pymod_project_temp_dir"
            project_temp_dir_path = os.path.join(self.current_pymod_directory, project_temp_dir_name)
            pmos.pymod_mkdir(project_temp_dir_path, overwrite_dirs=True)

            # Check if the file is a valid zip file.
            if not zipfile.is_zipfile(project_archive_file_path):
                raise Exception("The file is not a zip file.")
            zfh = open(project_archive_file_path, 'rb')
            zipfile_obj = zipfile.ZipFile(zfh)
            # Check if the file is a valid PyMod session file.
            files_to_check = ["%s.pkl" % project_name, "%s.pse" % project_name]
            if not set(files_to_check) < set(zipfile_obj.namelist()):
                zfh.close()
                raise Exception("The file is not a valid PyMod session file.")
            # Extract the file to a temporary directory.
            pmos.zipfile_extract_all(zipfile_obj, project_temp_dir_path)
            zfh.close()

        except Exception, e:
            self.load_project_failure(project_archive_file_path, e, project_temp_dir_path, continue_session=True)
            return None

        # If some errors happens here, close PyMod.
        try:
            # Unpickle the data of the saved PyMod project.
            project_pickle_file = open(os.path.join(project_temp_dir_path, "%s.pkl" % project_name), 'r')
            project_dict = PyMod_Unpickler(project_pickle_file).load().__dict__
            self.__dict__.update(project_dict)
            project_pickle_file.close()

            # Loads a PyMOL session.
            cmd.reinitialize()
            cmd.load(os.path.join(project_temp_dir_path, "%s.pse" % project_name))

            # Copies the current project files in the directory.
            os.chdir("..")
            shutil.rmtree(self.current_project_directory_full_path)
            src = os.path.join(project_temp_dir_path, os.path.basename(self.current_project_directory_full_path))
            dst = self.current_project_directory_full_path
            shutil.copytree(src, dst)
            os.chdir(self.current_project_directory_full_path)
            # Rebuilds missing directories in the project folder.
            dirs_to_write = (self.structures_directory, self.models_directory, self.alignments_directory, self.psipred_directory, self.similarity_searches_directory, self.images_directory, self.temp_directory_name)
            for dir_to_write in dirs_to_write:
                if not os.path.isdir(dir_to_write):
                    os.mkdir(dir_to_write)

            # Builds a .zip file of the directory.
            shutil.rmtree(project_temp_dir_path)

            # Updates PyMod main window.
            self.gridder()

        except Exception, e:
            self.load_project_failure(project_archive_file_path, e, project_temp_dir_path, continue_session=False)


    def load_project_failure(self, project_archive_file_path, error, project_temp_dir_path, continue_session=True):
        if os.path.isdir(project_temp_dir_path):
            shutil.rmtree(project_temp_dir_path)
        title = "Open Project Error"
        message = "Could not open the project file '%s': because of the following error: %s." % (project_archive_file_path, error)
        if not continue_session:
            message += " PyMod is now shutting down."
        self.show_error_message(title, message)
        if not continue_session:
            self.main_window.destroy()

    # Set the attributes to store when an object of this class is pickled.
    attributes_to_store = ["unique_index",
                           "alignment_count",
                           "new_clusters_counter",
                           "logo_image_counter",
                           "performed_modeling_count",
                           "multiple_chain_models_count",
                           "modeling_session_list",
                           "blast_cluster_counter",
                           "color_index",
                           "pdb_list",
                           "pymod_elements_list"]
    lists_to_pickle = {}
    attributes_to_pickle = []

    def __getstate__(self):
        """
        Used to build a dict with only the 'attributes_to_store' when pickling.
        """
        return pymod_pickle(self)


    ###############################################################################################
    # METHODS TO MANIPULATE THE ELEMENTS: POD (PyMod object model).                               #
    ###############################################################################################

    def get_mother_count(self):
        """
        Returns the total number of "mother" elements (both childess and with children) loaded in
        the pymod_elements_list.
        """
        count = 0
        for e in self.pymod_elements_list:
            if e.is_mother:
                count += 1
        return count

    def set_as_mother(self,element):
        element.is_mother = True
        element.show_children = True
        element.is_child = False
        element.child_index = 0

    def set_as_children(self,element):
        element.is_mother = False
        element.show_children = None
        element.is_child = True

    def add_element_to_pymod(self,element,level,mother_index=None,child_index=None,grid=False, color=None):
        """
        Used to add elements to the pymod_elements_list. Once an element is added to pymod_elements_list by
        this method, it will be displayed in the PyMod main window.
            element: the element object
            level: "mother" or "child"
            mother_index: if not provided it will be the last mother of the element list.
            child_index: if not provided it will be the last element of the cluster.
        """
        # Add a new mother.
        if level == "mother":
            self.set_as_mother(element)
            # Check that a child index is not be provided.
            if child_index != None:
                raise Exception("When adding a mother a child_index can't be provided.")
            # Append the new mother to the bottom of the mother list.
            if mother_index == None:
                element.mother_index = self.get_mother_count()
            # Insert the new mother somewhere in the list.
            # MAYBE REMOVE THIS AND USE ONLY move_mother TO PLACE IT INSIDE THE LIST LATER.
            else:
                # Checks that the index is correct. It works because the element has not been
                # appended to self.pymod_elements_list yet.
                self.check_mother_index(mother_index,adding_element_to_list=True)
                # This works just like the "insert_mother()" method.
                for e in self.pymod_elements_list:
                    if e.mother_index >= mother_index:
                        e.mother_index += 1
                element.mother_index = mother_index

        # Add a new chidren.
        elif level == "child":
            if mother_index == None:
                raise Exception("When adding a child element the index of its mother must be specified.")

            self.set_as_children(element)

            mother = self.get_mother_by_index(mother_index)
            last_child_index = len(self.get_children(mother)) + 1 # Get the last index of the mother.

            if child_index == None:
                element.mother_index = mother.mother_index
                element.child_index = last_child_index
            else:
                element.mother_index = mother.mother_index
                element.child_index = child_index

        # Assigns default colors for structures.
        if element.has_structure() and not element.is_model and color == None:
            color = self.color_struct()
        # Assigns white as color if the 'color' argument is not specified.
        else:
            if color == None:
                color = "white"
        element.my_color = color

        self.pymod_elements_list.append(element) # Adds the element to the pymod_elements_list.
        element.unique_index = self.unique_index
        self.unique_index += 1

        if grid:
            self.gridder()


    def check_mother_index(self, mother_index,adding_element_to_list=False):
        correct_index = False
        # print "result:",mother_index,self.get_mother_count()
        if mother_index > self.get_mother_count()-1 or mother_index < -1:
             if adding_element_to_list==False:
                 raise Exception("The mother index is out of limits.")
             elif adding_element_to_list==True:
                 raise Exception("The specified mother index must be lower than the current mother count. "+
                                  "If you want to give to this element the highest mother index possible "+
                                  "do not provide it with a mother_index when calling add_element_to_list.")
        else:
            correct_index = True
        return correct_index

    # Also make an insert_cluster(target_mother_index,cluster_mother_index) method.
    # This should be used only inside move_mother.
    def insert_mother(self,mother,new_mother_index):
        # Checks if the new_mother_index is ok.
        # self.check_mother_index(new_mother_index)
        for e in self.pymod_elements_list:
            if e != mother:
                if e.mother_index >= new_mother_index:
                    e.mother_index += 1
        mother.mother_index = new_mother_index

    # This should be used only inside move_mother and add_to_mother.
    def pop_mother(self,mother):
        for e in self.pymod_elements_list:
            if e != mother:
                if e.mother_index > mother.mother_index:
                    e.mother_index -= 1

    def move_mother(self,mother,new_mother_index):
        # If the mother has some children also move them.

        # Checks if the new_mother_index is ok.
        self.check_mother_index(new_mother_index)
        self.pop_mother(mother)
        self.insert_mother(mother,new_mother_index)


    def check_child_index(self, child,child_index,adding_element_to_list=False):
        correct_index = False
        mother = self.get_mother(child)
        if child_index > len(self.get_children(mother)) or child_index < -1 or child_index == 0:
                 raise Exception("The child index is out of limits.")
        else:
            correct_index = True
        return correct_index

    # This should be used only inside move_child.
    def insert_child(self,child,new_child_index):
        # Checks if the new_mother_index is ok.
        # self.check_child_index(new_child_index)
        for e in self.pymod_elements_list:
            if e.mother_index == child.mother_index and e.is_child and e != child:
                if e.child_index >= new_child_index:
                    e.child_index += 1
        child.child_index = new_child_index

    # This should be used only inside move_child.
    def pop_child(self,child):
        for e in self.pymod_elements_list:
            if e.mother_index == child.mother_index and e.is_child and e != child:
                if e.child_index > child.child_index:
                    e.child_index -= 1

    def move_child(self,child,new_child_index):
        # Checks if the new_mother_index is ok.
        self.check_child_index(child, new_child_index)
        self.pop_child(child)
        self.insert_child(child, new_child_index)


    # Returns the mother with that index.
    def get_mother_by_index(self, mother_index):
        self.check_mother_index(mother_index)
        e = None
        for element in self.pymod_elements_list:
            if element.is_mother and element.mother_index == mother_index:
                e = element
                break
        if e == None:
            raise Exception("Mother with mother_index " + str(mother_index) + "not found")
        else:
            return e

    def get_mother(self,child):
        mother = None
        for e in self.pymod_elements_list:
            if e.is_mother and e.mother_index == child.mother_index:
                mother = e
                break
        return mother

    # Returns the list of children of some mother.
    # Or call it just get_children(). Use as an argument the mother object itself.
    def get_children(self,mother): # get_mothers_children
       children_list = []
       for element in self.pymod_elements_list:
           if element.is_child and element.mother_index == mother.mother_index:
               children_list.append(element)
       return children_list


    def get_siblings(self,child):
        siblings_list = []
        for element in self.pymod_elements_list:
           if element.is_child and element.mother_index == child.mother_index and not element == child:
               siblings_list.append(element)
        return siblings_list


    def get_cluster_lead(self, cluster):
        lead = None
        for child in self.get_children(cluster):
            if child.is_lead:
                lead = child
                break
        return lead

    def cancel_memory(self,element):
        element.is_blast_query = False
        element.is_bridge = False
        element.is_lead = False

    # Appends an element to some cluster/mother.
    # Arguments (elements): mother, children
    def add_to_mother(self, mother, element, child_index=None):

        # Checks that the element to add is actually in the self.pymod_elements_list.
        # ...
        if mother.mother_index == element.mother_index:
            raise Exception("You can't add an element to some mother with the same index.")

        last_child_index = None
        # For mothers.
        if element.is_mother:
            self.pop_mother(element)
            self.set_as_children(element)
        elif element.is_child:
            self.pop_child(element)
            self.cancel_memory(element)
            # element.child_index = len(self.get_children(mother)) + 1
            pass

        # Get the last index of the mother.
        last_child_index = len(self.get_children(mother)) + 1
        if child_index != None:
            # Check if the child index is correct...
            # ...
            element.child_index = child_index
        else:
            element.mother_index = mother.mother_index
            element.child_index = last_child_index


    def extract_child(self,child):
        # This should check that the element is a child.
        self.pop_child(child)
        self.cancel_memory(child)
        self.set_as_mother(child)
        self.insert_mother(child, child.mother_index + 1)


    def extract_selection_to_new_cluster(self):
        selected_sequences = self.get_selected_sequences()
        self.new_clusters_counter += 1
        # A new 'PyMod_element' object to represent the new cluster.
        new_cluster_element = PyMod_element(
                "...", "New cluster %s" % (self.new_clusters_counter),
                element_type = "alignment", adjust_header=False,
                alignment_object = Alignment("new-cluster", self.new_clusters_counter))
        # Giving it the same mother index the original cluster, will place it above it in PyMod's
        # main window.
        self.add_element_to_pymod(
            new_cluster_element, level="mother",
            mother_index=selected_sequences[0].mother_index)
        # Adds the selected sequences to the new cluster.
        for seq in selected_sequences:
            self.add_to_mother(new_cluster_element, seq)
        self.set_initial_ali_seq_number(new_cluster_element)
        self.gridder()


    def mark_as_lead(self, new_lead):
        for child in self.get_siblings(new_lead):
            if child.is_lead:
                child.is_lead = False
        new_lead.is_lead = True


    def mark_as_query(self, query):
        query.is_blast_query = True
        self.mark_as_lead(query)


    def delete_element(self,element):
        if element.has_structure():
            self.delete_pdb_file(element)

        if element.is_mother:
            self.pop_mother(element)
            self.pymod_elements_list.remove(element)
        elif element.is_child:
            self.pop_child(element)
            self.pymod_elements_list.remove(element)


    def delete_pdb_file(self,element):
        # If the sequence has a PDB file loaded inside PyMOL, then delete it.
        try:
            cmd.delete(element.build_chain_selector_for_pymol())
        except:
            pass


    # Delete a whole cluster.
    def delete_whole_cluster(self,mother):
        # First delete the children.
        for c in self.get_children(mother):
            self.delete_element(c)
        # Then delete the mother.
        self.delete_element(mother)


    # Deletes alignments elements and extracts their children.
    def delete_alignment(self,alignment_element):
        children = reversed(sorted(self.get_children(alignment_element), key=lambda child:child.child_index))
        for child in children:
            self.extract_child(child)
        self.delete_element(alignment_element)


    def delete_clusters_with_one_child(self):
        """
        Finds cluster that have only one child and deletes them.
        """
        for e in self.pymod_elements_list:
            if e.is_mother and e.is_cluster_element():
                if len(self.get_children(e)) < 2:
                    self.delete_alignment(e)


    def duplicate_sequence(self, element_to_duplicate):
        if element_to_duplicate.has_structure():
            new_file_shortcut = os.path.join(self.structures_directory, element_to_duplicate.structure.original_pdb_file_name)
            target_chain_id = element_to_duplicate.structure.pdb_chain_id
            # Builds a 'Parsed_pdb_file' object.
            pdb_file = Parsed_pdb_file(os.path.abspath(new_file_shortcut))
            # Start parsing the PDB file.
            pdb_file.parse_pdb_file()
            # Builds 'Pymod_elements' objects for each chain present in the PDB file and adds the PDB
            # file to the record of PDB files loaded in PyMod.
            pdb_file.build_structure_objects(add_to_pymod_pdb_list = False)
            # Builds an element and load a structure only for target elements.
            for chain_id in pdb_file.get_chains_ids():
                if chain_id == target_chain_id:
                    new_element = pdb_file.get_chain_pymod_element(chain_id)
                    self.add_element_to_pymod(new_element, "mother")
                    self.load_element_in_pymol(new_element)
        else:
            c = PyMod_element(
                str(element_to_duplicate.my_sequence).replace("-",""), element_to_duplicate.my_header_fix,
                full_original_header= "Copy of " + element_to_duplicate.full_original_header,
                element_type="sequence")
            self.add_element_to_pymod(c,"mother")


    # Parses the clusterseq list and returns the element that has the unique_index of the argument.
    def get_element_by_unique_index(self, index):
        element = None
        for e in self.pymod_elements_list:
            if e.unique_index == index:
                element = e
                break
        return element


    def all_sequences_are_children(self, selection=None):
        """
        Returns True if all the elements selected by the user are children. A list of PyMod elements
        is not specified in the 'selection' argument, the target selection will be the list of
        sequences currently selected in the GUI.
        """
        if selection == None:
            selection = self.get_selected_sequences()
        if False in [e.is_child for e in selection]:
            return False
        else:
            return True

    def all_selected_elements_have_fetchable_pdbs(self, selection=None):
        """
        Returns True if all the elements selected by the user can be used to download a PDB file.
        """
        if selection == None:
            selection = self.get_selected_sequences()
        if False in [e.pdb_is_fetchable() for e in selection]:
            return False
        else:
            return True

    def all_sequences_have_structure(self, selection=None):
        """
        Returns True if all the elements selected by the user have structure loaded into PyMOL.
        """
        if selection == None:
            selection = self.get_selected_sequences()
        if False in [e.has_structure() for e in selection]:
            return False
        else:
            return True


    def gridder(self,clear_selection = True,rebuild_submenus=True):
        """
        Displays the PyMod elements currently loaded in the plugin main window.
        """
        # This will automatically delete all the cluster elements (alignments and blast-searches)
        # that have only one child.
        self.delete_clusters_with_one_child()
        # Updates the sequences of all the clusters.
        for c in self.get_cluster_elements():
            self.update_cluster_sequences(c)

        # Sorts the pymod_elements_list so that its elements will be displayed according to their
        # mother and child indices.
        self.pymod_elements_list.sort(key=lambda e: (e.mother_index,e.child_index))
        for element in self.pymod_elements_list:
            element.is_shown = False

        display_on_pymod_window = True

        if display_on_pymod_window:
            #| Destroy the RIGHT and LEFT panes (containing the displayed sequences)
            self.rightpan.destroy()
            self.leftpan.destroy()
            #| Creates 2 new RIGTH and LEFT panes
            self.create_main_window_panes()
            # Use the "create_entry()" method on the elemets to display them on the window to
            # display sequences.
            grid_index = 0
            current_font_size = self.menu_sequence_font_size.get()
            for (k,element) in enumerate(self.pymod_elements_list):
                # For mothers.
                if element.is_mother:
                    element.create_entry(grid_position = grid_index, font_size=current_font_size)
                    grid_index += 1
                # For children. Hide children in collapsed clusters.
                elif element.is_child and self.get_mother(element).show_children:
                    element.create_entry(grid_position = grid_index, font_size=current_font_size)
                    grid_index += 1

        if rebuild_submenus:
            self.build_alignment_submenu()
            self.build_models_submenu()

        if clear_selection:
            self.deselect_all_elements()


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
            elements_to_color.extend(self.get_selected_sequences())
        elif mode == "all":
            elements_to_color.extend(self.get_all_sequences())

        # Actually proceeds to color the elements.
        for seq in elements_to_color:
            if color_scheme == "regular":
                seq.color_element_by_regular_color(regular_color)
            elif color_scheme == "residue":
                seq.color_element_by_residue()
            elif color_scheme == "secondary-observed":
                seq.color_element_by_obs_sec_str()
            elif color_scheme == "secondary-predicted":
                seq.color_element_by_pred_sec_str()
            # Colors elements with 3D structure according to the observed II str, elements with
            # predicted II str according to the prediction, and leaves the other elements unaltered.
            elif color_scheme == "secondary-auto":
                if seq.has_structure():
                    seq.color_element_by_obs_sec_str()
                elif seq.has_predicted_secondary_structure():
                    seq.color_element_by_pred_sec_str()
            elif color_scheme == "campo-scores":
                seq.color_element_by_campo_scores()
            elif color_scheme == "dope":
                seq.color_element_by_dope()


    def deselect_all_elements(self):
        for e in self.pymod_elements_list:
            e.selected = False

    def get_all_sequences(self):
        all_elements = [e for e in self.pymod_elements_list if not e.is_cluster_element()]
        return all_elements

    def get_selected_elements(self, include_hidden_children=False):
        """
        Returns a list of all the elements (both sequences and clusters) selected by the user.
        """
        selected_elements = []
        if not include_hidden_children:
            selected_elements = [e for e in self.pymod_elements_list if e.selected]
        else:
            for e in self.pymod_elements_list:
                if e.selected:
                    selected_elements.append(e)
                    if e.is_lead_of_collapsed_cluster():
                        mother = self.get_mother(e)
                        selected_elements.append(mother)
                        for s in self.get_siblings(e):
                            selected_elements.append(s)
            for e in selected_elements:
                e.selected = True

        return selected_elements

    def get_selected_sequences(self):
        """
        Returns a list of all the sequences selected by the user.
        """
        selected_sequences = [e for e in self.pymod_elements_list if e.selected and not e.is_cluster_element()]
        return selected_sequences


    def get_cluster_elements(self,cluster_type = "all"):
        """
        Returns only those elements in pymod_elements_list with cluster_type = "alignment" or
        "blast-search".
        """
        cluster_elements = []
        for element in self.pymod_elements_list:
            if element.is_cluster_element():
                if cluster_type == "all":
                    cluster_elements.append(element)
                elif cluster_type == "alignment" and element.element_type == "alignment":
                    cluster_elements.append(element)
                elif cluster_type == "blast-search" and element.element_type == "blast-search":
                    cluster_elements.append(element)
        return cluster_elements


    ###############################################################################################
    # FILES MANAGMENT.                                                                            #
    ###############################################################################################

    #################################################################
    # Opening sequence files.                                       #
    #################################################################

    def open_file_from_the_main_menu(self):
        """
        This method is called when new sequences are loaded from the main menu.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        file_names = askopenfilename(filetypes=pmdt.supported_file_types, multiple=True, parent=self.main_window)
        for single_file_name in pmos.get_askopenfilename_tuple(file_names):
            extension = os.path.splitext(os.path.basename(single_file_name))[1].replace(".","")
            if extension.lower() == "fasta":
                if self.is_sequence_file(single_file_name, "fasta"):
                    self.open_sequence_file(single_file_name, "fasta")
            elif extension.lower() == "gp":
                if self.is_sequence_file(single_file_name, "genbank"):
                    self.open_sequence_file(single_file_name, "genbank")
            elif extension.lower() == "pdb":
                if self.is_pdb(single_file_name):
                    self.open_pdb_file(single_file_name)
            elif extension.lower() == "ent":
                if self.is_pdb(single_file_name):
                    self.open_pdb_file(single_file_name)


    def is_sequence_file(self, file_name, file_format, show_error=True):
        """
        Try to open a sequence file using Biopython. Returns 'True' if the file is a valid fasta file.
        """
        valid_file = False
        file_handler = None
        try:
            file_handler = open(file_name,"r")
            r = list(SeqIO.parse(file_handler, file_format))
            if len(r) > 0:
                valid_file = True
            else:
                valid_file = False
        except:
            valid_file = False
        if file_handler != None:
            file_handler.close()
        if not valid_file:
            title = "File Type Error"
            message = "The selected File is not a valid %s." % (pmdt.supported_sequence_file_types[file_format])
            if show_error:
                self.show_error_message(title,message)
        return valid_file


    def is_pdb(self,file_name, show_error=True):
        valid_pdb = False
        file_handler = open(file_name, "r")
        for line in file_handler.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    valid_pdb = True
                    break
                except:
                    pass
        file_handler.close()
        if not valid_pdb and show_error:
            title = "FileType Error"
            message = "The selected File is not a valid PDB."
            self.show_error_message(title,message)
        return valid_pdb


    def raw_seq_input(self):
        """
        Launched when the user wants to add a new sequence by directly typing it into a Text entry.
        """
        def show_menu(e):
            w = e.widget
            the_menu.entryconfigure("Paste",
            command=lambda: w.event_generate("<<Paste>>"))
            the_menu.tk.call("tk_popup", the_menu, e.x_root, e.y_root)

        # This is called when the SUBMIT button packed below is pressed.
        def submit():
            def special_match(strg, search=re.compile(r'[^A-Z-]').search):
                return not bool(search(strg))
            def name_match(strg, search2=re.compile(r'[^a-zA-Z0-9_]').search):
                return not bool(search2(strg))
            sequence = textarea.get(1.0, "end").replace('\n','').replace('\r','').replace(' ','').replace('\t','').upper()
            if special_match(sequence) and len(sequence):
                if len(seq_name.get()) and name_match(seq_name.get()):
                    c = PyMod_element(sequence, seq_name.get(),
                        element_type="sequence")
                    self.add_element_to_pymod(c,"mother")
                    self.raw_seq_window.destroy()
                    self.gridder()
                else:
                    title = 'Name Error'
                    message = 'Please Check The Sequence Name:\n  Only Letters, Numbers and "_" Allowed'
                    self.show_error_message(title,message,parent_window=self.raw_seq_window,refresh=False)
            else:
                title = 'Sequence Error'
                message = 'Please Check Your Sequence:\n Only A-Z and "-" Allowed'
                self.show_error_message(title,message,parent_window=self.raw_seq_window,refresh=False)

        self.raw_seq_window = pmgi.PyMod_tool_window(self.main_window,
            title = "Add Raw Sequence",
            upper_frame_title = "Type or Paste your Sequence",
            submit_command = submit)

        L1 = Label(self.raw_seq_window.midframe,font = "comic 12", text="Name:", bg="black", fg= "red")
        L1.grid( row=0, column=0, sticky="e", pady=5, padx=5)

        # Creates an Entry for the name of the new sequence.
        seq_name=Entry(self.raw_seq_window.midframe, bd=0, disabledforeground = 'red', disabledbackground = 'black',
                    selectbackground = 'black', selectforeground = 'white', width=60, font = "%s 12" % pmgi.fixed_width_font)
        seq_name.grid( row=0, column=1,columnspan=2, sticky="nwe", pady=5, )
        seq_name.focus_set()
        seq_name.bind("<Button-3><ButtonRelease-3>", show_menu)

        L2 = Label(self.raw_seq_window.midframe, text="Sequence: ", bg="black", fg= "red", font = "comic 12")
        L2.grid( row=1, column=0, sticky="ne", ipadx=0, padx=5)

        scrollbar = Scrollbar(self.raw_seq_window.midframe)
        scrollbar.grid(row=1, column=2, sticky="ns")

        # Creates an Entry widget for the sequence.
        textarea=Text(self.raw_seq_window.midframe, yscrollcommand=scrollbar.set,
                      font = "%s 12" % pmgi.fixed_width_font, height=10,
                      bd=0, foreground = 'black',
                      background = 'white', selectbackground='black',
                      selectforeground='white', width = 60)
        textarea.config(state=NORMAL)
        textarea.tag_config("normal", foreground="black")
        textarea.grid( row=1, column=1, sticky="nw", padx=0)
        textarea.bind("<Button-3><ButtonRelease-3>", show_menu)
        scrollbar.config(command=textarea.yview)

        the_menu = Menu(seq_name, tearoff=0)
        the_menu.add_command(label="Paste")


    def open_sequence_file(self, file_full_path, file_format="fasta"):
        """
        Method for opening primary sequence files (FASTA, and also others should be included).
        It only needs the path of the file to open.
        """
        fn = open(file_full_path, "rU")
        # Parses a fasta file through Biopython. This will automatically crop headers that have " "
        # (space) characters.
        for record in SeqIO.parse(fn, file_format):
            # Then builds a PyMod_element object and add it to the pymod_elements_list.
            c = self.build_pymod_element_from_seqrecord(record)
            self.add_element_to_pymod(c,"mother")
        fn.close()
        self.gridder()


    def build_pymod_element_from_seqrecord(self, seqrecord):
        """
        Gets Biopython a 'SeqRecord' class object and returns a 'PyMod_element' object corresponding
        to the it.
        """
        new_element = PyMod_element(str(seqrecord.seq), seqrecord.id,
                full_original_header=seqrecord.description,
                element_type="sequence")
        return new_element


    def build_pymod_element_from_hsp(self, hsp):
        """
        Gets a hsp dictionary containing a Biopython 'HSP' class object and returns a
        'PyMod_element' object corresponding to the subject in the HSP.
        """
        # Gives them the query mother_index, to make them its children.
        record_header = self.correct_name(hsp["title"])
        cs = PyMod_element(
            str(hsp["hsp"].sbjct),
            record_header,
            full_original_header=hsp["title"],
            element_type="sequence")
        return cs

    #################################################################
    # Opening alignment files.                                      #
    #################################################################

    def choose_alignment_file(self):
        """
        Lets users choose an alignment file.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        openfilename = askopenfilename(filetypes=pmdt.alignment_file_formats, multiple=False,parent=self.main_window)
        if openfilename == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(os.path.basename(openfilename))[1].replace(".","") # BUG.
        if extension == "fasta":
            pass
        elif extension == "aln":
           extension = "clustal"
        # Unknown format.
        else:
            title = "Format Error"
            message = "Unkwnown alignment file format: %s" % (extension)
            self.show_error_message(title,message)
            return (None, None)
        return openfilename, extension


    def open_alignment_from_main_menu(self):
        """
        Lets users import in Pymod an alignment stored in an external file.
        """
        openfilename, extension = self.choose_alignment_file()
        if not None in (openfilename, extension):
            self.build_cluster_from_alignment_file(openfilename, extension)


    def build_cluster_from_alignment_file(self,alignment_file, extension="fasta"):
        """
        Creates a cluster with all the sequences contained in an alignment file.
        """
        aligned_elements = []
        fh = open(alignment_file, "rU")
        records = SeqIO.parse(fh, extension)
        for record in records:
            e = self.build_pymod_element_from_seqrecord(record)
            aligned_elements.append(e)
        fh.close()
        # Creates an alignment element.
        self.alignment_count += 1
        alignment_name = self.set_alignment_element_name("imported",self.alignment_count)
        imported_alignment_element = PyMod_element("...", alignment_name,
            element_type = "alignment",
            alignment_object = Alignment("imported",self.alignment_count), adjust_header=False)
        self.add_element_to_pymod(imported_alignment_element, "mother")
        # Adds the sequences to the new alignment cluster.
        for element in aligned_elements:
            self.add_element_to_pymod(element, "child", mother_index=imported_alignment_element.mother_index)
        # Computes the stars of the new alignment element.
        self.update_stars(imported_alignment_element)
        # Sets the initial number of sequences in the alignment.
        self.set_initial_ali_seq_number(imported_alignment_element)
        self.gridder()


    def transfer_alignment(self,alignment_element):
        """
        Changes the sequences of the elements contained in a PyMod cluster according to the
        information presente in an externally supplied file (chosen by users through a file diaolog)
        containing the same sequences aligned in a different way. Right now it supports transfer
        only for sequences having the exactly same sequences in PyMod and in the external alignment.
        """
        # Let users choose the external alignment file.
        openfilename, extension = self.choose_alignment_file()
        if None in (openfilename, extension):
            return False

        # Sequences in the aligment currently loaded into PyMod.
        aligned_elements = self.get_children(alignment_element)

        # Sequences in the alignment files.
        fh = open(openfilename, "rU")
        external_records = list(SeqIO.parse(fh, extension))
        fh.close()

        if len(external_records) < len(aligned_elements):
            title = "Transfer error"
            message = "'%s' has more sequences (%s) than the alignment in '%s' (%s) and the 'Transfer Alignment' function can't be used in this situation." % (alignment_element.my_header, len(aligned_elements), openfilename, len(external_records))
            self.show_error_message(title,message)
            return False

        correspondance_list = []
        # First try to find sequences that are identical (same sequence and same lenght) in both
        # alignments.
        for element in aligned_elements[:]:
            identity_matches = []
            for record in external_records:
                if str(element.my_sequence).replace("-","") == str(record.seq).replace("-",""):
                    match_dict = {"target-seq":element, "external-seq": record, "identity": True}
                    identity_matches.append(match_dict)
            if len(identity_matches) > 0:
                correspondance_list.append(identity_matches[0])
                aligned_elements.remove(identity_matches[0]["target-seq"])
                external_records.remove(identity_matches[0]["external-seq"])

        # Then try to find similar sequences among the two alignments. Right now this is not
        # implemented.
        # ...

        if not len(aligned_elements) == 0:
            title = "Transfer error"
            message = "Not every sequence in the target alignment has a corresponding sequence in the external alignment."
            self.show_error_message(title,message)
            return False

        # Finally transfer the sequences.
        for match in correspondance_list[:]:
            if match["identity"]:
                match["target-seq"].my_sequence = str(match["external-seq"].seq)
                correspondance_list.remove(match)

        self.gridder()


    #################################################################
    # Opening PDB files.                                            #
    #################################################################

    def choose_structure_file(self):
        """
        Lets users choose a strcture file.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        openfilename = askopenfilename(filetypes=pmdt.structure_file_types, multiple=False,parent=self.main_window)
        if openfilename == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(os.path.basename(openfilename))[1].replace(".","")
        return openfilename, extension


    def open_pdb_file(self,pdb_file_full_path, grid=True):
        """
        Opens a PDB file (specified in 'pdb_file_full_path'), reads its content and imports in PyMod
        the sequences of the polypeptide chains and loads in PyMOL their 3D structures.
        """
        # Builds a 'Parsed_pdb_file' object.
        pdb_file = Parsed_pdb_file(pdb_file_full_path)
        # Copies the original PDB file to the Structures directory in the current project folder.
        pdb_file.copy_to_structures_directory()
        # Start parsing the PDB file.
        pdb_file.parse_pdb_file()
        # Builds 'Pymod_elements' objects for each chain present in the PDB file and adds the PDB
        # file to the record of PDB files loaded in PyMod.
        pdb_file.build_structure_objects(add_to_pymod_pdb_list = True)
        # Actually adds as mothers the PyMod elements to the 'pymod_elements_list'.
        for chain_element in pdb_file.get_chains_pymod_elements():
            new_element = chain_element
            self.add_element_to_pymod(new_element, "mother")
            self.load_element_in_pymol(new_element)
        # Displays the sequences of the PDB chains when all of them are loaded into Pymod.
        if grid:
            self.gridder()


    def load_element_in_pymol(self, element, mode = None):
        """
        Loads the PDB structure of the chain into PyMol.
        """
        chain_root_name = element.build_chain_selector_for_pymol()
        file_name_to_load = os.path.join(pymod.structures_directory, chain_root_name+".pdb")
        cmd.load(file_name_to_load)
        cmd.select("last_prot", chain_root_name)
        cmd.hide("everything", "last_prot")
        cmd.show("cartoon", "last_prot" ) # Show the new chain as a cartoon.
        if mode == "model":
            cmd.color("white", "last_prot")
        else:
            cmd.color(element.my_color, "last_prot")
        cmd.util.cnc("last_prot") # Colors by atom.
        cmd.center('last_prot')
        cmd.zoom('last_prot')
        cmd.delete("last_prot")


    class Select_chain_and_first_model(Select):
        """
        This is needed to write new PDB files with only a single chain of the first model of the
        biopython object parsed by Bio.PDB.PDBIO().
        """
        def __init__(self,chain_id):
            self.chain_id = chain_id

        def accept_chain(self, chain):
            if chain.get_id() == self.chain_id:
                return True
            else:
                return False

        def accept_model(self,model):
            if model.id==0:
                return True
            else:
                return False


    def fetch_pdb_files(self, mode, target_selection):
        """
        Function for downloading a PDB file from the sequences retrived from BLAST.
        """
        # Builds a list of structures to be fetched.
        self.structures_to_fetch = []
        if mode == "single":
            self.structures_to_fetch.append(target_selection)
        elif mode == "selection":
            self.structures_to_fetch.extend(self.get_selected_sequences())

        # Let the user choose the way in which to retrieve the structures.
        import_all_text = 'Import all chains'
        import_single_text = 'Import only the hit sequences fragments'
        self.import_mode_choices = {import_single_text: "single-chain", import_all_text: "multiple-chains"}
        self.fetch_pdb_dialog = Pmw.MessageDialog(self.main_window,
            title = 'Import Options',
            message_text = (
            "Please select the 3D structure import mode:\n\n"+
            "- Import in PyMod the structure of every chain of the PDB files.\n\n"+
            "- Import in PyMod only the structure of the hit sequences fragments identified by (PSI-)BLAST."
            ),
            buttons = (import_all_text, import_single_text) )
        self.fetch_pdb_dialog.component("message").configure(justify="left")
        self.fetch_pdb_dialog.configure(command=self.fetch_pdb_files_state)


    def fetch_pdb_files_state(self, dialog_choice):
        self.fetch_pdb_dialog.withdraw()
        # Interrupt the process if users close the dialog window.
        if not dialog_choice:
            return None
        import_mode = self.import_mode_choices[dialog_choice]

        # Begins to actually fetch the PDB files.
        for element in self.structures_to_fetch:
            element_header = element.my_header
            if element_header.split("|")[2] == "pdb":
                pdb_code = element_header.split("|")[3]
                if element_header.split("|")[4] != "":
                    pdb_chain = element_header.split("|")[4][0]
                else:
                    pdb_chain = None
                    import_mode = "multiple-chains"
            elif element_header.split("|")[4] == "pdb":
                pdb_code=element_header.split("|")[5]
                if element_header.split("|")[6][0] != "":
                    pdb_chain = element_header.split("|")[6][0]
                else:
                    pdb_chain = None
                    mport_mode = "multiple-chains"

            zipped_file = None

            # Retrieve the PDB file from the internet.
            try:
                zipped_file = urllib.urlretrieve('http://www.rcsb.org/pdb/files/'+ pdb_code + '.pdb.gz')[0]
            except:
                title = "Connection Error"
                message = "Can not access to the PDB database.\nPlease check your Internet access."
                self.show_error_message(title,message)
                return False

            open_zipped_file = gzip.open(zipped_file) # Uncompress the file while reading
            new_name = pdb_code + '.pdb' # Form the pdb output name
            pdb_file_shortcut = os.path.join(self.structures_directory, new_name)
            saved_file = open(pdb_file_shortcut, 'w')
            saved_file.write(open_zipped_file.read()) # Write pdb file
            open_zipped_file.close()
            saved_file.close()

            # Builds a 'Parsed_pdb_file' object.
            pdb_file = Parsed_pdb_file(os.path.abspath(pdb_file_shortcut))
            # Start parsing the PDB file.
            pdb_file.parse_pdb_file()

            # Load in PyMod only the chain corresponding to the hit sequence and adjust its legth to
            # the region identified by BLAST.
            if import_mode == "single-chain":
                if not self.associate_structure(pdb_file, pdb_chain, element):
                    self.show_associate_structure_error()

            # Load each chain found in the PDB file where the 3D structure of the hit sequence is
            # present. This is actually like opening a new PDB file with the 'open_pdb_file()'
            # method, except that in this case, the chains not corresponging to the hit sequence
            # are colored in gray.
            elif import_mode == "multiple-chains":
                # Builds 'Pymod_elements' objects for each chain present in the PDB file.
                pdb_file.build_structure_objects(add_to_pymod_pdb_list = True)
                if pdb_chain:
                    # Actually adds as mothers the PyMod elements to the 'pymod_elements_list'.
                    for chain_id in pdb_file.get_chains_ids():
                        new_element = pdb_file.get_chain_pymod_element(chain_id)
                        if chain_id != pdb_chain:
                            self.add_element_to_pymod(new_element, "mother", color="gray")
                        else:
                            # Deletes the original hit sequence retrieved by BLAST and replaces it with
                            # a new element with an associated structure loaded in PyMOL.
                            self.delete_element(element)
                            self.add_element_to_pymod(new_element, "mother")
                        self.load_element_in_pymol(new_element)
                else:
                    for chain_id in pdb_file.get_chains_ids():
                        new_element = pdb_file.get_chain_pymod_element(chain_id)
                        self.add_element_to_pymod(new_element, "mother")
                        self.load_element_in_pymol(new_element)
                    self.delete_element(element)

        self.gridder()


    def associate_structure_from_popup_menu(self, target_element):
        """
        Launched when users press the 'Associate 3D Structure' from the leeft popup menu.
        """
        # This will be set to 'True' once the users select a valid PDB file and press the 'SUBMIT'
        # button.
        self.select_associate_chain = False
        # This will contain the 'Parsed_pdb_file' object of the structure to associate.
        self.associate_pdb_file = None
        # The 'PyMod_element' object of the target seequence (the sequence to be associated with a
        # structure).
        self.associate_target_element = target_element

        # Builds a new window.
        self.associate_structure_window = pmgi.PyMod_tool_window(self.main_window,
            title = "Associate Structure",
            upper_frame_title = "Associate 3D Structure Options",
            submit_command = self.associate_structure_state)

        # An entryfield to select the structure file.
        self.structure_file_enf = pmgi.PyMod_path_entryfield(self.associate_structure_window.midframe,
            label_text = "Select Structure File",
            label_style = pmgi.label_style_1,
            path_type = "file",
            file_types = pmdt.structure_file_types,
            askpath_title = "Select Structure File")
        self.structure_file_enf.pack(**pmgi.pack_options_1)
        self.associate_structure_window.add_widget_to_align(self.structure_file_enf)
        self.associate_structure_window.add_widget_to_validate(self.structure_file_enf)

        self.associate_structure_window.align_widgets(15)


    def associate_structure_state(self):
        # Checks if a correct structure file has been provided as input.
        if not self.select_associate_chain:
            if not self.check_general_input(self.associate_structure_window):
                return False
            pdb_file_path = self.structure_file_enf.getvalue()

            if not self.is_pdb(pdb_file_path, show_error=False):
                title = "File Type Error"
                message = "Please select a valid PDB file."
                self.show_error_message(title,message, parent_window=self.associate_structure_window)
                return False
            # Removes the entryfield to select the structure file.
            self.structure_file_enf.pack_forget()

            # Parses the structure file.
            self.associate_pdb_file = Parsed_pdb_file(pdb_file_path)
            self.associate_pdb_file.copy_to_structures_directory()
            self.associate_pdb_file.parse_pdb_file()
            # Gets its chains.
            available_chains = self.associate_pdb_file.get_chains_ids()

            # Displays a combobox to select the chain id of corresponind to the structure to be
            # associated with the target sequence.
            self.chain_selection_cbx = pmgi.PyMod_combobox(self.associate_structure_window.midframe,
                label_text = 'Select Chain to Associate',
                label_style = pmgi.label_style_1,
                scrolledlist_items=available_chains)
            self.chain_selection_cbx.pack(**pmgi.pack_options_1)
            self.chain_selection_cbx.selectitem(0)
            self.associate_structure_window.add_widget_to_align(self.chain_selection_cbx)
            self.associate_structure_window.align_widgets(15)

            self.select_associate_chain = True

        # If a valid structure file has been provided, this will try to associate the structure of
        # the chain specified in the combobox to the target element.
        elif self.select_associate_chain:
            if not self.associate_structure(self.associate_pdb_file, self.chain_selection_cbx.get(), self.associate_target_element):
                self.show_associate_structure_error(parent_window = self.associate_structure_window)
                return False
            self.associate_structure_window.destroy()
            self.gridder()


    def show_associate_structure_error(self, parent_window = None):
        title = "Associate Structure Failure"
        message = "The amminoacid sequences of the target chain and the chain in the PDB structure do not match."
        self.show_error_message(title, message, parent_window = parent_window)


    def associate_structure(self, parsed_pdb_file, chain_id, pymod_element):
        """
        Gets a 'Parsed_pdb_file' object and a 'PyMod_element' object as arguments, and associates
        the structure with chain id specified in 'chain_id' to the PyMod element.
        """
        # Builds 'Pymod_elements' objects for each chain present in the PDB file.
        parsed_pdb_file.build_structure_objects(add_to_pymod_pdb_list = False)
        # Crops the structure.
        sequences_match = parsed_pdb_file.crop_structure_chain(chain_id, adjust_to_sequence = pymod_element.my_sequence)
        if not sequences_match:
            return False
        # Build a 'PyMod_element' object representing the cropped chain and transfer its
        # data to the element if the hit sequence.
        cropped_element = parsed_pdb_file.get_chain_pymod_element(chain_id)
        pymod_element.update_element(new_sequence=cropped_element.my_sequence, new_header=cropped_element.my_header, new_structure=cropped_element.structure)
        pymod_element.my_color = self.color_struct()
        self.load_element_in_pymol(pymod_element)
        parsed_pdb_file.add_to_pdb_list()
        return True


    def import_selections(self):
        """
        Method for importing PyMOL Selections into PyMod. It saves PyMOL objects selected by users
        to file, and loads it into PyMOL using 'open_pdb_file()'.
        """
        # Find all structures already loaded into PyMod: items in struct_list are excluded from
        # importable PyMOL object list.
        struct_list=[]
        for member in self.pymod_elements_list:
            if member.has_structure():
                struct_list.append(member.build_chain_selector_for_pymol())

        scrolledlist_items=[] # importable PyMOL objects
        for obj in cmd.get_names("objects"):
            if not obj in struct_list and cmd.get_type(obj) == "object:molecule":
                scrolledlist_items.append(str(obj))

        if not len(scrolledlist_items):
            if struct_list:
                self.show_error_message("No Importabled Object", "All PyMOL objects are already imported into PyMod.")
            else:
                self.show_error_message("No Importabled Object", "No PyMOL object to import.")
            return

        # Builds a new window.
        self.import_from_pymol_window = pmgi.PyMod_tool_window(self.main_window,
            title = "Import from PyMOL",
            upper_frame_title = "Load PyMOL Objects into PyMod",
            submit_command = self.import_selected_pymol_object)

        # Builds a combobox for each PyMOL object to import.
        self.combobox_frame = Frame(self.import_from_pymol_window.midframe, background='black')
        self.combobox_frame.pack(side = TOP, fill = BOTH, anchor="center", ipadx = 5, ipady = 5, pady=5)
        self.sele_var=dict() # whether a PyMOL object is selected
        self.sele_checkbutton=dict() # checkbuttons for object selection
        row=0 # vetical location of checkbuttons
        for sele in scrolledlist_items:
            self.sele_var[sele]=IntVar()
            self.sele_checkbutton[sele]=Checkbutton(self.combobox_frame,
                text=sele, variable=self.sele_var[sele],
                background="black", foreground="white",
                selectcolor="red", highlightbackground="black")
            self.sele_checkbutton[sele].grid(row=row,column=0,sticky='w')
            row+=1


    def import_selected_pymol_object(self):
        selected_num=0
        for sele in self.sele_var:
            if self.sele_var[sele].get():
                selected_num+=1
                filename=sele+".pdb"
                pdb_file_shortcut = os.path.join(self.temp_directory_name, filename)
                cmd.save(pdb_file_shortcut,sele)
                cmd.delete(sele)
                self.open_pdb_file(os.path.abspath(pdb_file_shortcut))
        if not selected_num:
            tkMessageBox.showerror( "Selection Error",
                "Please select at least one object to import.")
        else:
            self.import_from_pymol_window.destroy()


    def pymol_save(self, filepath, selection_name):
        """
        Calls the 'cmd.save' function of PyMOL, but retains the order of the atoms in the original
        PDB file.
        """
        try:
            old_retain_order = cmd.get("retain_order")
            cmd.set("retain_order", 1)
        except:
            pass
        cmd.save(filepath, selection_name)
        try:
            cmd.set("retain_order", old_retain_order)
        except:
            pass


    def show_pdb_info(self):
        self.work_in_progress()

    #################################################################
    # Saving files.                                                 #
    #################################################################

    def save_all_files_from_main_menu(self):
        """
        Saves all files in a single FASTA file.
        """
        if len(self.pymod_elements_list) != 0:
            self.save_selection(mode="all")
        else:
            self.show_error_message("Selection Error","There aren't any sequences currently loaded in PyMod.")


    def sequence_save(self, element):
        """
        Save a single sequence to a file.
        """
        remove_indels_choice = False
        if "-" in element.my_sequence:
            remove_indels_choice = tkMessageBox.askyesno(message="Would you like to remove indels from the sequence when saving it to a file?", title="Save File", parent=self.main_window)

        filepath=asksaveasfilename(filetypes=[("fasta","*.fasta")],parent=self.main_window)

        if not filepath == "":
            dirpath = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            self.build_sequences_file([element], filename, file_format="fasta", remove_indels=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def save_selection(self, mode="selection"):
        """
        Save selection in a single file.
        """
        # Builds the selection.
        selection = None
        if mode == "selection":
            selection = self.get_selected_sequences()
        elif mode == "all":
            selection = self.get_all_sequences()

        # Ask users if they want to include indels in the sequences to save.
        remove_indels_choice = False
        for e in selection:
            if "-" in e.my_sequence:
                remove_indels_choice = tkMessageBox.askyesno(message="Would you like to remove indels from the sequences when saving them to a file?", title="Save Selection", parent=self.main_window)
                break

        # Ask users to chose a directory where to save the files.
        filepath=asksaveasfilename(filetypes=[("fasta","*.fasta")],parent=self.main_window)
        if not filepath == "":
            dirpath = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            self.build_sequences_file(selection, filename, file_format="fasta", remove_indels=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def alignment_save(self, alignment_element):
        """
        Lets the user choose the path to which an alignment file is going to be saved, and saves
        an alignment file there.
        """
        save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = pmdt.alignment_file_formats, parent=pymod.main_window)
        alignment_file_name, extension = os.path.splitext(os.path.basename(save_file_full_path))
        extension = extension.replace(".","")

        if save_file_full_path != "":
            # The get all the aligned elements.
            aligned_elements = self.get_children(alignment_element)

            # Saves a file with all the sequences in the project "Alignments" directory.
            if extension == "fasta":
                self.save_alignment_fasta_file(alignment_file_name, aligned_elements)
            elif extension == "aln":
                self.build_sequences_file(aligned_elements, alignment_file_name, file_format="clustal", remove_indels=False)
            else:
                title = "Format Error"
                message = "Unknown alignment file format: %s" % (extension)
                self.show_error_message(title, message)
                return

            # Moves the saved file to the path chosen by the user.
            try:
                old_path = os.path.join(self.alignments_directory, alignment_file_name + "." + extension)
                os.rename(old_path, save_file_full_path)
            except:
                title = "File Error"
                message = "Could not save the alignment file to path: %s" % (save_file_full_path)
                self.show_error_message(title, message)


    def save_alignment_fasta_file(self, file_name, aligned_elements, first_element=None):
        """
        Saves in the Alignments directory a .fasta alignment file containing the sequences of the
        "aligned_elements".
        """
        self.build_sequences_file(aligned_elements, file_name, file_format="fasta", remove_indels=False,first_element=first_element)


    ###############################################################################################
    # HEADER AND SEQUENCES MANIPULATION.                                                          #
    ###############################################################################################

    def correct_name(self, name):
        """
        This function allows to rename a sequence with the same name (needed by ClustalW to work
        properly). Checks if there are other elements in the pymod_elements_list that have the same
        name. If there are, then append to the name of the sequence a string to diversifity it as a
        copy.
        """
        def get_new_name(name, n=1, name_root=None):
            if name_root == None:
                name_root = name
            if name in [e.my_header for e in self.pymod_elements_list]:
                new_name = str(n)+"_"+name_root
                return get_new_name(new_name, n+1, name_root)
            else:
                return name
        new_correct_name = get_new_name(name)
        return new_correct_name


    def build_header_string(self, unformatted_header):
        formatted_header = unformatted_header[0:150].replace(" ","_")
        formatted_header = formatted_header.replace("/","_")
        formatted_header = formatted_header.replace(":","_")
        return formatted_header


    def correct_sequence(self, sequence):
        sequence = sequence.replace("Z","X") # ambiguity, E or Q
        sequence = sequence.replace("B","X") # ambiguity, D or N
        sequence = sequence.replace("J","X") # ambiguity, I or L
        sequence = sequence.replace("O","X") # pyrrolysine
        sequence = sequence.replace("U","X") # selenocysteine
        sequence = sequence.replace(".","X") # selenocysteine
        return sequence


    def adjust_aligned_elements_length(self,elements,remove_right_indels=True):
        # First remove indels at the end of the sequences.
        if remove_right_indels:
            for e in elements:
                e.my_sequence = str(e.my_sequence).rstrip("-")
        # Then pad each sequence with the right number of indels to make them of the same length as
        # the longest sequence.
        max_length = max([len(e.my_sequence) for e in elements])
        for e in elements:
            e.my_sequence = str(e.my_sequence).ljust(max_length,"-")


    def update_cluster_sequences(self, cluster_element):
        """
        Updates the sequences of a cluster when some sequences are removed or added from the
        cluster.
        """
        children = self.get_children(cluster_element)
        self.adjust_aligned_elements_length(children)
        self.update_stars(cluster_element)


    def update_stars(self, cluster_element):
        stars = self.compute_stars(self.get_children(cluster_element))
        cluster_element.my_sequence = stars


    def one2three(self, letter):
        """
        Returns a three letter code for a residue corresponding to a one letter symbol.
        """
        if pmdt.prot_one_to_three_code.has_key(letter):
            return pmdt.prot_one_to_three_code[letter]
        else:
            return "???"

    def three2one(self, res, force_standard_parent=False):
        """
        Returns a one letter symbol corresponding to a three letter code for a residue. If
        'force_standard_parent' is set to 'True', if the three letter code of a modified residue
        is supplied, the method will attempt to return the code of the corresponding unmodified
        residue.
        """
        # Try to set modified residues of the protein with the letter belonging to the unmodified
        # residue
        if not force_standard_parent:
            # If it is any kind of known modified amminoacid set it as its original non
            # modified amminoacid
            if pmdt.code_standard.has_key(res):
                return pmdt.code_standard[res]
            else:
                return "X"
        else:
            pass


    def get_polymer_type(self, sequence):
        polymer_type = "protein"
        nucleotides = [nt for nt in pmdt.nucleic_acids_dictionary.keys()]
        list_of_three_letter_codes = [r.three_letter_code for r in sequence]
        for res in list_of_three_letter_codes:
            if res in nucleotides:
                polymer_type = "nucleic-acid"
                break
        return polymer_type


    def color_struct(self):
        if self.color_index > len(pmdt.regular_colours) - 1:
            self.color_index=0
        color_index_to_return = self.color_index
        self.color_index += 1
        return pmdt.regular_colours[color_index_to_return]


    ###############################################################################################
    # SELECTION MENU COMMANDS.                                                                    #
    ###############################################################################################
    def select_all_from_main_menu(self):
        for element in self.get_all_sequences():
            if not element.selected:
                element.toggle_element()

    def deselect_all_from_main_menu(self):
        for element in self.get_all_sequences():
            if element.selected:
                element.toggle_element()


    def show_all_structures_from_main_menu(self):
        for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
            element.show_chain_in_pymol()

    def hide_all_structures_from_main_menu(self):
        for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
            element.hide_chain_in_pymol()


    def select_all_structures_from_main_menu(self):
        for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
            if not element.selected:
                element.toggle_element()

    def deselect_all_structures_from_main_menu(self):
        for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
            if element.selected:
                element.toggle_element()


    def expand_all_clusters_from_main_menu(self):
        for element in self.get_cluster_elements():
            element.expand_cluster()

    def collapse_all_clusters_from_main_menu(self):
        for element in self.get_cluster_elements():
            element.collapse_cluster()


    ###############################################################################################
    # ALIGNMENT MENU AND ITS BEHAVIOUR.                                                           #
    ###############################################################################################

    def build_alignment_submenu(self):
        """
        Build an "Alignment N" voice in the "Alignments" submenu when alignment N is performed.
        """
        # Delete the old alignment submenu.
        self.alignments_menu.delete(0,500)

        # Then rebuilds it with the new alignments.
        alignment_list = self.get_cluster_elements()

        if alignment_list != []:
            for element in alignment_list:
                uid = element.unique_index
                alignment = element.alignment

                # Alignment menu for each cluster loaded in PyMod.
                alignment_submenu = Menu(self.alignments_menu, tearoff = 0)
                # Save to a file dialog.
                alignment_submenu.add_command(label = "Save to File",
                        command = lambda ui=uid: self.save_alignment_to_file_from_ali_menu(ui))
                alignment_submenu.add_separator()

                # Matrices submenu.
                matrices_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Matrices", menu = matrices_submenu)
                matrices_submenu.add_command(label = "Identity matrix",
                    command = lambda ui=uid: self.display_identity_matrix(ui))
                if alignment.algorithm in self.can_show_rmsd_matrix and alignment.rmsd_list != None:
                    matrices_submenu.add_command(label = "RMSD matrix",
                        command = lambda ui=uid: self.display_rmsd_matrix(ui))

                # Trees.
                if alignment.initial_number_of_sequence > 2:
                    trees_submenu = Menu(alignment_submenu, tearoff = 0)
                    alignment_submenu.add_cascade(label = "Trees", menu = trees_submenu)
                    if alignment.algorithm in self.can_show_guide_tree:
                        trees_submenu.add_command(label = "Show Guide Tree",
                            command = lambda ui=uid: self.show_guide_tree_from_alignments_menu(ui))
                    if alignment.algorithm in self.can_show_dendrogram and 0:
                        trees_submenu.add_command(label = "Show Dendrogram",
                            command = lambda ui=uid: self.show_dendrogram_from_alignments_menu(ui))
                    if len(self.get_children(element)) >= 2:
                        trees_submenu.add_command(label = "Build Tree from Alignment",
                            command = lambda ui=uid: self.build_tree_from_alignments_menu(ui))

                # Evolutionary conservation.
                evolutionary_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Evolutionary Conservation", menu = evolutionary_submenu)
                evolutionary_submenu.add_command(label = "CAMPO",
                    command = lambda ui=uid: self.build_campo_window(ui))
                if alignment.algorithm in self.can_use_scr_find and 0:
                    evolutionary_submenu.add_command(label = "SCR_FIND",
                        command = lambda ui=uid: self.build_scr_find_window(ui))

                # Render alignment.
                render_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Render Alignment", menu = render_submenu)
                render_submenu.add_command(label = "Generate Logo through WebLogo 3",
                    command = lambda ui=uid: self.build_logo_options_window(ui))
                render_submenu.add_command(label = "Launch ESPript in Web Browser",
                    command = lambda ui=uid: self.espript(ui))

                # Adds the alignment submenu to the PyMod main menu.
                label_text = element.my_header
                self.alignments_menu.add_cascade(label = label_text, menu = alignment_submenu)

        else:
            self.alignments_menu.add_command(label = "There aren't any alignments")


    def save_alignment_to_file_from_ali_menu(self,alignment_unique_id):
        self.alignment_save(self.get_element_by_unique_index(alignment_unique_id))

    #################################################################
    # CAMPO.                                                        #
    #################################################################

    def build_campo_window(self, alignment_unique_id):
        """
        Builds a window with opotions for the CAMPO algorithm.
        """
        self.input_alignment_element = self.get_element_by_unique_index(alignment_unique_id)

        current_pack_options = pmgi.pack_options_1
        current_label_options = pmgi.label_style_1

        # Builds the window.
        self.campo_window = pmgi.PyMod_tool_window(self.main_window,
            title = "CAMPO algorithm options",
            upper_frame_title = "Here you can modify options for CAMPO",
            submit_command = self.campo_state)

        # Scoring matrix combobox.
        self.campo_matrices = ["Blosum90","Blosum80","Blosum62","Blosum50","Blosum45","PAM30","PAM120","PAM250" ]
        self.campo_matrices_dict = {"Blosum62": "blosum62", "Blosum90": "blosum90","Blosum80":"blosum80",
                                    "Blosum50": "blosum50", "Blosum45":"blosum45",
                                    "PAM30": "pam30", "PAM120": "pam120", "PAM250": "pam250"}
        self.matrix_cbx = pmgi.PyMod_combobox(self.campo_window.midframe, label_text = 'Scoring Matrix Selection',label_style = current_label_options, scrolledlist_items=self.campo_matrices)
        self.matrix_cbx.pack(**current_pack_options)
        self.matrix_cbx.selectitem(2)
        self.campo_window.add_widget_to_align(self.matrix_cbx)

        # Gap open entryfield.
        self.campo_gap_penalty_enf = pmgi.PyMod_entryfield(
            self.campo_window.midframe,
            label_text = "Gap Score",
            label_style = current_label_options,
            value = '-1',
            validate = {'validator' : 'integer',
                        'min' : -1000, 'max' : 0})
        self.campo_gap_penalty_enf.pack(**current_pack_options)
        self.campo_window.add_widget_to_align(self.campo_gap_penalty_enf)

        # Gap extension entryfield.
        self.campo_gap_to_gap_score_enf = pmgi.PyMod_entryfield(
            self.campo_window.midframe,
            label_text = "Gap to Gap Score",
            label_style = current_label_options,
            value = '0',
            validate = {'validator' : 'integer',
                        'min' : -1000, 'max' : 0})
        self.campo_gap_to_gap_score_enf.pack(**current_pack_options)
        self.campo_window.add_widget_to_align(self.campo_gap_to_gap_score_enf)

        # Toss gaps.
        self.campo_exclude_gaps_rds = pmgi.PyMod_radioselect(self.campo_window.midframe, label_text = 'Toss gaps')
        for text in ('Yes', 'No'):
            self.campo_exclude_gaps_rds.add(text)
        self.campo_exclude_gaps_rds.setvalue('Yes')
        self.campo_exclude_gaps_rds.pack(**current_pack_options)
        self.campo_window.add_widget_to_align(self.campo_exclude_gaps_rds)

        self.campo_window.align_widgets(10)


    def campo_state(self):
        """
        Called when the "SUBMIT" button is pressed on the CAMPO window. Contains the code to compute
        CAMPO scores using the 'CAMPO' class.
        """
        # Saves a .fasta file for the alignment.
        aligned_sequences = self.get_children(self.input_alignment_element)
        self.save_alignment_fasta_file("temp", aligned_sequences)
        input_file_shortcut = os.path.join(self.alignments_directory,"temp.fasta")

        # Computes CAMPO scores by using the campo module.
        cbc = campo.CAMPO(input_file_shortcut,
                          mutational_matrix = self.campo_matrices_dict[self.matrix_cbx.get()],
                          gap_score = int(self.campo_gap_penalty_enf.getvalue()),
                          gap_gap_score = int(self.campo_gap_to_gap_score_enf.getvalue()),
                          toss_gaps = pmdt.yesno_dict[self.campo_exclude_gaps_rds.getvalue()])
        cbc.compute_id_matrix()
        cbc.run_CAMPO()

        # Gets the list of CAMPO score. There are as many values as positions in the alignment.
        campo_list = cbc.get_campo_items_list()

        # Assigns CAMPO scores to each one of the aligned sequences.
        for seq in aligned_sequences:
            seq.campo_scores = []
            for (r,v) in zip(seq.my_sequence,campo_list):
                if r != "-":
                    seq.campo_scores.append(v)
            seq.color_element_by_campo_scores()

        # Removes the temporary alignment file.
        os.remove(input_file_shortcut)
        self.campo_window.destroy()


    def build_scr_find_window(self, alignment_unique_id):
        pass


    #################################################################
    # Build and display sequence identity and RMSD matrices of      #
    # alignments.                                                   #
    #################################################################
    def display_identity_matrix(self,alignment_unique_id):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        # Get the cluster element.
        alignment_element = self.get_element_by_unique_index(alignment_unique_id)
        # Then get all its children (the aligned elements).
        aligned_elements = self.get_children(alignment_element)
        n = len(aligned_elements)

        # identity_matrix = [[None]*n]*n # [] # Builds an empty (nxn) "matrix".
        identity_matrix = []
        for a in range(n):
            identity_matrix.append([None]*n)

        # Computes the identities (or anything else) and builds the matrix.
        for i in range(len(aligned_elements)):
            for j in range(len(aligned_elements)):
                if j >= i:
                    sid = pmsm.compute_sequence_identity(aligned_elements[i].my_sequence,aligned_elements[j].my_sequence)
                    # This will fill "half" of the matrix.
                    identity_matrix[i][j] = sid
                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    identity_matrix[j][i] = sid

        identity_matrix = numpy.array(identity_matrix)

        # Build the list of sequences names.
        # Adjust it for PDB names...
        sequences_names = []
        for e in aligned_elements:
            sequences_names.append(e.get_compact_header())

        title = 'Identity matrix for ' + alignment_element.my_header
        self.show_table(sequences_names, sequences_names, identity_matrix, title)


    def display_rmsd_matrix(self,alignment_unique_id):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        # Get the cluster element.
        alignment_element = self.get_element_by_unique_index(alignment_unique_id)
        # Then get all its children (the aligned elements).
        aligned_elements = self.get_children(alignment_element)

        rmsd_list = alignment_element.alignment.rmsd_list
        rmsd_matrix_to_display = []
        n = len(aligned_elements)
        for a in range(n):
            rmsd_matrix_to_display.append([None]*n)

        for i,ei in enumerate(aligned_elements):
            for j,ej in enumerate(aligned_elements):
                if j >= i:
                    # This will fill "half" of the matrix.
                    rmsd = rmsd_list[(ei.unique_index,ej.unique_index)]
                    rmsd_matrix_to_display[i][j] = rmsd
                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    rmsd = rmsd_list[(ej.unique_index,ei.unique_index)]
                    rmsd_matrix_to_display[j][i] = rmsd

        # Build the list of sequences names.
        sequences_names = []
        for e in aligned_elements:
            sequences_names.append(e.get_compact_header())

        title = 'RMSD matrix for ' + alignment_element.my_header
        self.show_table(sequences_names, sequences_names, rmsd_matrix_to_display, title)


    def show_table(self, column_headers=None, row_headers=None, data_array=[], title = "New Table", columns_title = None, rows_title = None, number_of_tabs=2, width=800, height=450, rowheader_width=20):
        """
        Displayes in a new window a table with data from the bidimensional 'data_array' numpy array.
        """
        mfont = "monospace 10"
        number_of_newlines = 2

        # Builds a new window in which the table will be displayed.
        new_window = Toplevel(self.main_window)
        new_window.title(title)

        # Create a ScrolledText with headers.
        self.matrix_widget = Pmw.ScrolledText(
                new_window, borderframe = 1,
                usehullsize = True, hull_width = width, hull_height = height,
                columnheader = True, rowheader = True,
                text_padx = 20, text_pady = 20, Header_padx = 20,
                text_wrap='none', text_font = mfont, Header_font = mfont, # Header_foreground = 'blue',
                rowheader_width = rowheader_width, rowheader_pady = 20, rowheader_padx = 20 )

        # Create the row headers.
        for row in row_headers:
            row = str(row)
            self.matrix_widget.component('rowheader').insert('end', row+"\n"*number_of_newlines)

        # Create the column headers
        header_line = ''
        for column in column_headers:
            column_text = str(column) + "\t"*number_of_tabs
            header_line = header_line + column_text
        self.matrix_widget.component('columnheader').insert('0.0', header_line)

        # Enters the data.
        for i,row_items in enumerate(data_array):
            data_line = ""
            for item in row_items:
                column_text = str(item) + "\t"*number_of_tabs
                data_line += column_text
            if i != len(data_array):
                data_line += "\n"*number_of_newlines
            self.matrix_widget.insert('end', data_line)

        # Prevent users' modifying text and headers and packs it.
        self.matrix_widget.configure(text_state = 'disabled', Header_state = 'disabled')
        self.matrix_widget.pack(padx = 5, pady = 5, fill = 'both', expand = 1)


    #################################################################
    # Show guide trees and build trees out of alignments.           #
    #################################################################

    def show_guide_tree_from_alignments_menu(self,alignment_unique_id):
        """
        Shows the guide tree that was constructed in order to perform a multiple alignment.
        """
        # Gets the path of the .dnd file of the alignment.
        alignment_element = self.get_element_by_unique_index(alignment_unique_id)
        dnd_file_path = alignment_element.alignment.get_dnd_file_path()
        self.show_tree(dnd_file_path)


    def show_tree(self, tree_file_path):
        # Reads a tree file using Phylo.
        tree = Phylo.read(tree_file_path, "newick")
        tree.ladderize() # Flip branches so deeper clades are displayed at top
        # Displayes its content using PyMod plotting engine.
        pplt.draw_tree(tree, self.main_window)


    def show_dendrogram_from_alignments_menu(self,alignment_unique_id):
        """
        Shows dendrograms built by SALIGN.
        """
        # Gets the path of the .dnd file of the alignment.
        alignment_element = self.get_element_by_unique_index(alignment_unique_id)
        tree_file_path = alignment_element.alignment.get_dnd_file_path()
        pmsp.draw_salign_dendrogram(tree_file_path)


    ###################################
    # Tree building.                  #
    ###################################

    def build_tree_from_alignments_menu(self, alignment_unique_id):
        """
        Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        menu. It will check if a software to build a tree is available on the user's machine.
        """

        self.input_alignment_element = self.get_element_by_unique_index(alignment_unique_id)
        self.tree_building_software = None

        can_build_tree = False

        if self.clustalw.exe_exists():
            self.tree_building_software = "clustalw"
            can_build_tree = True
        elif self.muscle.exe_exists():
            self.tree_building_software = "muscle"
            can_build_tree = True

        if can_build_tree:
            self.build_tree_building_window()

        else:
            title = "Tree building Error"
            message = "In order to build a tree out of an alignment you need to install either ClustalW or MUSCLE."
            self.show_error_message(title, message)


    def check_tree_constructor_module(self):
        try:
            import Bio.Phylo.TreeConstruction
            return True
        except:
            return False


    def build_tree_building_window(self):
        """
        Builds a window with options to build a tree out of an alignment.
        """
        current_pack_options = pmgi.pack_options_1

        # Builds the window.
        self.tree_building_window = pmgi.PyMod_tool_window(self.main_window,
            title="Options for Tree Building",
            upper_frame_title="Here you can modify options for Tree Building",
            submit_command=self.run_tree_building_software)

        # Add some options.
        self.algorithm_rds = pmgi.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Clustering Algorithm')
        for alg_name in (sorted(pmdt.tree_building_alg_dict.keys())):
            self.algorithm_rds.add(alg_name)
        self.algorithm_rds.setvalue("Neighbor Joining")
        self.algorithm_rds.pack(**current_pack_options)
        self.tree_building_window.add_widget_to_align(self.algorithm_rds)

        if self.tree_building_software == "clustalw":
            # Kimura distance correction.
            self.distance_correction_rds = pmgi.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Use Distance Correction')
            for text in ('Yes', 'No'):
                self.distance_correction_rds.add(text)
            self.distance_correction_rds.setvalue('No')
            self.distance_correction_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.distance_correction_rds)

            # Toss gaps.
            self.exclude_gaps_rds = pmgi.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Exclude Gaps')
            for text in ('Yes', 'No'):
                self.exclude_gaps_rds.add(text)
            self.exclude_gaps_rds.setvalue('No')
            self.exclude_gaps_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.exclude_gaps_rds)

        self.tree_building_window.align_widgets(13)


    def run_tree_building_software(self):
        # Saves a temporary input alignment file.
        alignment_file_name = "alignment_tmp"
        alignment_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.fasta')
        self.save_alignment_fasta_file(alignment_file_name, self.get_children(self.input_alignment_element))
        clustering_algorithm = self.get_clustering_algorithm()

        commandline = ""
        output_file_path = None

        if self.tree_building_software == "clustalw":
            commandline =  '"%s"' % (self.clustalw.get_exe_file_path())
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
            output_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.ph')

        elif self.tree_building_software == "muscle":
            commandline =  '"%s"' % (self.muscle.get_exe_file_path())
            commandline += ' -maketree -in %s' % (alignment_file_path)
            output_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.phy')
            commandline += ' -out %s' % (output_file_path)
            if clustering_algorithm == "nj":
                commandline += ' -cluster neighborjoining'
            elif clustering_algorithm == "upgma":
                pass

        self.execute_subprocess(commandline)

        # Remove temporary files.
        ali_id = self.input_alignment_element.alignment.id
        new_tree_file_path = os.path.join(self.alignments_directory, self.alignments_files_names+str(ali_id)+"_align_tree" + ".phy")
        os.rename(output_file_path, new_tree_file_path)
        os.remove(alignment_file_path)

        self.tree_building_window.destroy()

        # Reads the output tree file with Phylo and displays its content using PyMod plotting
        # engine.
        self.show_tree(new_tree_file_path)


    def get_clustering_algorithm(self):
        return pmdt.tree_building_alg_dict[self.algorithm_rds.getvalue()]

    def get_boostrap_val(self):
        return pmdt.yesno_dict[self.bootstrap_rds.getvalue()]

    def get_distance_correction_val(self):
        return pmdt.yesno_dict[self.distance_correction_rds.getvalue()]

    def get_exclude_gaps_val(self):
        return pmdt.yesno_dict[self.exclude_gaps_rds.getvalue()]

    #################################################################
    # Methods for accessing the WebLogo web service.                #
    #################################################################

    def build_logo_options_window(self, alignment_unique_id):
        """
        Launched from the 'Alignments' menu on PyMod main menu. Displayes a window with a series of
        widgets through which users can define WebLogo parameters.
        """
        self.logo_window = pmgi.PyMod_tool_window(
            self.main_window,
            title = "WebLogo 3 web-application Options",
            upper_frame_title = "Here you can modify options for WebLogo 3",
            submit_command = self.logo_state,
            with_frame=True)

        #Units list.
        units_list=['Bits', 'Probability']
        #Units combobox.
        self.unit_combobox = pmgi.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Unit Selection',
            scrolledlist_items=units_list)
        self.unit_combobox.pack(**pmgi.pack_options_1)
        self.unit_combobox.selectitem(0)
        self.logo_window.add_widget_to_align(self.unit_combobox)

        #Color scheme list.
        colorscheme_list=['Auto', '(AA) Charge', '(AA) Chemistry', '(AA default) Hydrophobicity', '(NA) Classic', '(NA default) Base pairing']
        colorscheme_list.sort()
        #Color combobox.
        self.color_combobox = pmgi.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Color Scheme Selection',
            scrolledlist_items=colorscheme_list)
        self.color_combobox.pack(**pmgi.pack_options_1)
        self.color_combobox.selectitem(5)
        self.logo_window.add_widget_to_align(self.color_combobox)

        self.logo_al_element = self.get_element_by_unique_index(alignment_unique_id)
        self.AL_LENGTH = len(self.logo_al_element.my_sequence)

        #Sub-frame created to display entries for Logo Range option
        self.range_subframe = Frame(self.logo_window.midframe, background='black')
        self.range_subframe.pack(**pmgi.pack_options_1)
        #Logo Range Label
        self.logo_range_label=Label(self.range_subframe, text= "Logo Range", **pmgi.label_style_1 )
        self.logo_range_label.grid(row=0, column=0, sticky = "w", padx = (0,100))
        #Entry: Logo Start Position
        self.logo_start=Spinbox(self.range_subframe, from_=1, to=self.AL_LENGTH, width=5)
        self.logo_start.grid(row=0, column=1, sticky = "e")
        #Separator dash
        self.logo_range_dash=Label(self.range_subframe, font = "comic 10", height = 1,
                         text= " - ", background='black', fg='white')
        self.logo_range_dash.grid(row=0, column=2, sticky = "e")
        #Entry: Logo End Position
        self.logo_end=Spinbox(self.range_subframe, to=self.AL_LENGTH, width=5)
        self.logo_end.grid(row=0, column=3, sticky = "e")
        self.logo_end.insert(0, self.AL_LENGTH)
        self.logo_end.config(from_=2)

        # ADVANCED OPTIONS.
        self.logo_window.show_advanced_button()

        #Logo Format
        format_list=['PDF', 'PNG image']
        #Logo format combobox.
        self.format_combobox = pmgi.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Logo Format',
            scrolledlist_items=format_list)
        self.format_combobox.selectitem(0)
        self.logo_window.add_widget_to_align(self.format_combobox)
        self.logo_window.add_advanced_widget(self.format_combobox)

        #LOGO title entry.
        self.logo_title_enf = pmgi.PyMod_entryfield(self.logo_window.midframe,
            label_text = 'Logo Title',
            value = "")
        self.logo_window.add_widget_to_align(self.logo_title_enf)
        self.logo_window.add_advanced_widget(self.logo_title_enf)
        self.logo_window.add_widget_to_validate(self.logo_title_enf)

        #Stacks per line entry (default:80).
        self.logo_stacks_enf = pmgi.PyMod_entryfield(self.logo_window.midframe,
            label_text = 'Stacks per line',
            value = 80,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.logo_window.add_widget_to_align(self.logo_stacks_enf)
        self.logo_window.add_advanced_widget(self.logo_stacks_enf)
        self.logo_window.add_widget_to_validate(self.logo_stacks_enf)

        #Option: Scale stacks width.
        self.scale_width_rds = pmgi.PyMod_radioselect(self.logo_window.midframe, label_text = 'Scale stacks width')
        for text in ('Yes', 'No'):
            self.scale_width_rds.add(text)
        self.scale_width_rds.setvalue('No')
        self.logo_window.add_widget_to_align(self.scale_width_rds)
        self.logo_window.add_advanced_widget(self.scale_width_rds)

        #Option: Show error bars.
        self.show_error_rds = pmgi.PyMod_radioselect(self.logo_window.midframe, label_text = 'Show error bars')
        for text in ('Yes', 'No'):
            self.show_error_rds.add(text)
        self.show_error_rds.setvalue('No')
        self.logo_window.add_widget_to_align(self.show_error_rds)
        self.logo_window.add_advanced_widget(self.show_error_rds)

        self.logo_window.align_widgets(13)


    def check_logo_correct_parameters(self):
        '''
        Checks if the values that were insert in the LOGO window are correct.
        '''
        correct_input = True #This variable defines the status
        try:
            #checks if entries are integer numbers
            start=int(self.logo_start.get())
            end=int(self.logo_end.get())
            # Filters the BLAST record according to the advanced options.
            if self.logo_window.showing_advanced_widgets:
                stacks_pl = int(self.logo_stacks_enf.getvalue())
            #check on the logic of choosing extremities
            if start >= end:
                correct_input=False
                errortitle = "Input Error"
                errormessage = "Start value cannot be greater than the end value.\nPlease correct."
                self.show_error_message(errortitle, errormessage)
            elif start > self.AL_LENGTH or end > self.AL_LENGTH or start<0 or end<0:
                correct_input=False
                errortitle = "Input Error"
                errormessage = "Values cannot be greater than the sequence length and both must be greater then 0.\nPlease correct."
                self.show_error_message(errortitle, errormessage)
        except:
            correct_input=False
            errortitle = "Input Error"
            errormessage = "Non valid numeric input.\nPlease correct."
            self.show_error_message(errortitle, errormessage)
        return correct_input


    def logo_state(self):
        """
        This method is called when the 'Submit' button on the LOGO window is pressed. It runs a
        check on the entries, if they are correct it calls the getLogo() function
        """
        if not self.check_logo_correct_parameters():
            return False
        self.getLogo()


    def getLogo(self):
        '''
        Generates a LOGO of the alignment, by using WebLogo 3 site.
        Requires active Internet connection.
        '''
        #Units dictionary
        UNITS = {'Bits':'bits', 'Probability':'probability'}
        #Color scheme dictionary
        COLOR_SCHEME = {
            'Auto':'color_auto',
            '(NA default) Base pairing':'color_base_pairing',
            '(NA) Classic':'color_classic',
            '(AA default) Hydrophobicity':'color_hydrophobicity',
            '(AA) Chemistry':'color_chemistry',
            '(AA) Charge':'color_charge'
            }
        #Format dictionary
        FORMATS =  {'PNG image' : 'png_print',    'PDF' : 'pdf'}
        #switch format-extension
        extensions =  {'png_print': 'png',    'pdf' : 'pdf'}
        logo_yesno = {"Yes": "true", "No": "false"}

        #Options defined in the window
        LOGO_UNIT            = UNITS[self.unit_combobox.get()]
        LOGO_COLOR           = COLOR_SCHEME[self.color_combobox.get()]
        LOGO_RANGE_START     = self.logo_start.get()
        LOGO_RANGE_END       = self.logo_end.get()
        #Options defined in advanced options sub-window, not always visible. Here they are initialised.
        LOGO_FORMAT          = 'pdf'
        LOGO_TITLE           = ''
        LOGO_STACKS_PER_LINE = '80'
        LOGO_SCALE_STACKS    = 'false'
        LOGO_SHOW_ERRORBARS  = 'false'

        if self.logo_window.showing_advanced_widgets:
            LOGO_FORMAT          = FORMATS[self.format_combobox.get()]
            LOGO_TITLE           = self.logo_title_enf.getvalue()
            LOGO_STACKS_PER_LINE = self.logo_stacks_enf.getvalue()
            LOGO_SCALE_STACKS    = logo_yesno[self.scale_width_rds.getvalue()]
            LOGO_SHOW_ERRORBARS  = logo_yesno[self.show_error_rds.getvalue()]
        self.logo_window.destroy()

        print 'Running GetLogo...'

        #weblogo3 URL
        weblogourl = 'http://weblogo.threeplusone.com/create.cgi'

        #Sets fields and arguments collecting values from the LOGO options window
        values = {'unit_name': LOGO_UNIT, 'color_scheme': LOGO_COLOR,
                  'logo_start': LOGO_RANGE_START, 'logo_end'  : LOGO_RANGE_END,
                  'format': LOGO_FORMAT, 'logo_title': LOGO_TITLE,
                  'stacks_per_line': LOGO_STACKS_PER_LINE,
                  'show_xaxis': 'true', 'show_yaxis': 'true',
                  'show_ends': 'true', 'show_fineprint': 'true', }
        values_update_scale = {'scale_width': LOGO_SCALE_STACKS}
        values_update_errorbars = {'show_errorbars': LOGO_SHOW_ERRORBARS}

        if LOGO_SCALE_STACKS != 'false':
            values.update(values_update_scale)
        if LOGO_SHOW_ERRORBARS != 'false':
            values.update(values_update_errorbars)

        # Builds an url with the multiple alingment and WebLogo parameters and sends a request to
        # the WebLogo server.
        upload_response = self.upload_alignment(self.logo_al_element, weblogourl, 'sequences_file', other_values=values)

        #Check if valid response is given
        if upload_response:
            #Writes output content in a file with extension given by LOGO_FORMAT
            logofile = os.path.join(self.images_directory,'logo_' + str(self.logo_image_counter) + '.' + extensions[LOGO_FORMAT])
            lf = open(logofile, 'wb')
            print 'Creating file...'
            lf.write(upload_response)
            lf.close()
            self.logo_image_counter += 1
            pmos.open_document_with_default_viewer(logofile)
            print 'Done!'
        else:
            print 'No response. Aborted.'
            title = "Error"
            message = "No valid response from server"
            self.show_error_message(title,message)


    #################################################################
    # Methods for accessing the ESPript web service.                #
    #################################################################

    def espript(self, alignment_unique_id):
        '''
        Opens in the default browser the ESPript page, with the current alignment pre-loaded.
        Requires active Internet connection. It needs also the Schubert server to be reachable.
        '''
        # Prepares the target alignment element.
        self.espript_alignment_element = self.get_element_by_unique_index(alignment_unique_id)
        # A list of the header names of those aligned sequences with an associated 3D structure.
        self.espript_structures_list = ["None"]
        self.espript_structures_dict = {"None": None}
        for structure_element in filter(lambda e: e.has_structure(), self.get_children(self.espript_alignment_element)):
            self.espript_structures_list.append(structure_element.my_header)
            # Populates 'espript_structures_dict' so that the structures PDB file names can be
            # accessed by using as keys their header names.
            self.espript_structures_dict.update({structure_element.my_header: structure_element})
        if len(self.espript_structures_list) == 1:
            self.espript_state()
        else:
            self.show_espript_window()


    def show_espript_window(self):
        """
        Displayes a window with a combobox to let users select a strucure file of which the
        secondary structure information will be included in ESPript output.
        """
        self.espript_sec_str_window = pmgi.PyMod_tool_window(
            self.main_window,
            title = "ESPript Options",
            upper_frame_title = "Here you can modify options for ESPript",
            submit_command = self.espript_state )
        #Units combobox.
        self.espript_sec_str_combobox = pmgi.PyMod_combobox(self.espript_sec_str_window.midframe,
            label_text = 'Show Secondary Structure of',
            scrolledlist_items=self.espript_structures_list)
        self.espript_sec_str_combobox.pack(**pmgi.pack_options_1)
        self.espript_sec_str_combobox.selectitem(0)
        self.espript_sec_str_window.add_widget_to_align(self.espript_sec_str_combobox)
        self.espript_sec_str_window.align_widgets(15)


    def espript_state(self):
        """
        Uploads a sequence alignment file in fasta format on schubert (and optionally a structure
        file in the pdb format) and then opens a new tab on users' web browser with the ESPript page
        with the fasta (and the pdb) uploaded files a input.
        """
        schubert_url = 'http://schubert.bio.uniroma1.it/uploader/php_upload.php'
        schubert_folder_url = 'http://schubert.bio.uniroma1.it/uploader/uploads/'
        espript_basic_url = 'http://espript.ibcp.fr/ESPript/cgi-bin/ESPript.cgi?FRAMES=YES&amp;alnfile0='

        selected_structure_element = None
        if len(self.espript_structures_list) > 1:
            selected_structure_element = self.espript_structures_dict[self.espript_sec_str_combobox.get()]

        if selected_structure_element != None:
            upload_response = self.upload_alignment(self.espript_alignment_element, schubert_url, 'sequences_file', structure_element = selected_structure_element)
        else:
            upload_response = self.upload_alignment(self.espript_alignment_element, schubert_url, 'sequences_file')

        print 'Attempting to upload...'

        if len(self.espript_structures_list) > 1:
            self.espript_sec_str_window.destroy()

        #Checks if the upload is successful
        print upload_response
        if upload_response.startswith('TRUE'):
            if selected_structure_element == None:
                uploaded_alignment_file = upload_response[6:]
            else:
                uploaded_alignment_file, uploaded_structure_file= upload_response[6:].split(",")
            espript_url = espript_basic_url+schubert_folder_url+uploaded_alignment_file   #creates the URL
            if selected_structure_element != None:
                espript_url += ";struct1file0=%s%s" % (schubert_folder_url,uploaded_structure_file)
                espript_url += ";struct1chain0=%s" % (selected_structure_element.structure.pdb_chain_id)
            webbrowser.open(espript_url)    #opens the URL
            print 'Done'
        else:
            title = "Error"
            message = "Error while uploading the file. Please try again later or check your Internet connection."
            self.show_error_message(title,message)


    #################################################################
    # Common methods for interacting with web services.             #
    #################################################################

    def upload_alignment(self, alignment_element, url, form_upload_file_name, structure_element = None, other_values={}):
        '''
        This function creates a POST request to the URL 'url'. The 'form_upload_file_name' argument is the
        name of the form field that encodes the file to be uploaded. For instance: if in the upload form
        the field of the file is called "sequence_file", the form_upload_file_name argument has to be set to
        'sequence_file'. It's equivalent to the 'name' variable of the UNIX command curl:
            curl --form name=@content
        The function saves the current alignment and sends it to the server. It may also send other data,
        encoded in 'other_values' dictionary (a dictionary containing the parameters normally sent by compiling
        a form in the HTML page). This argument is optional and by default is an empty dictionary.
        Returns the response given by the server as a string.
        '''
        response_content = ''

        #Saves alignment in FASTA format
        alignment_file_name='alignment_tmp'
        self.save_alignment_fasta_file(alignment_file_name, self.get_children(alignment_element), first_element=structure_element)
        alignment_file_path=os.path.join(self.alignments_directory, alignment_file_name + '.fasta')

        #Copy file content to a string
        al_file = open(alignment_file_path)
        alignment_string = al_file.read()
        al_file.close()
        os.remove(alignment_file_path)
        print alignment_string

        values={form_upload_file_name: alignment_string}

        # Adds other values to the url.
        if other_values:
            values.update(other_values)
        # Uploads also a structure file.
        if structure_element != None:
            # values.update(other_values)
            structure_file = open(os.path.join(self.structures_directory, structure_element.structure.chain_pdb_file_name))
            structure_file_string = structure_file.read()
            dbref_line = "DBREF %s" % (structure_element.my_header).ljust(80, " ")
            structure_file_string = dbref_line + "\n" + structure_file_string
            structure_file.close()
            values.update({"structure_file": structure_file_string})

        user_agent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_6_8)'
        headers = { 'User-Agent' : user_agent }

        try:
            #Creates a request
            data = urllib.urlencode(values)
            req = urllib2.Request(url, data, headers=headers)
            #Gets server response and reads it
            response = urllib2.urlopen(req)
            response_content = response.read()
        except:
            response_content = ''
            title = "Connection Error"
            message = "Can not access the server.\nPlease check your Internet access."
            self.show_error_message(title,message)

        return response_content


    ###############################################################################################
    # MODELS MENU AND ITS BEHAVIOUR.                                                              #
    ###############################################################################################
    def build_models_submenu(self):
        """
        Build an "Modeling Session n" voice in the "Models" submenu once some models have been
        built.
        """
        self.models_menu.delete(0,500)

        if self.modeling_session_list != []:
            for modeling_session in self.modeling_session_list:
                modeling_session_submenu = Menu(self.models_menu, tearoff = 0)
                modeling_session_submenu.add_command(label = "DOPE Profile",
                    command = lambda ms=modeling_session: self.show_session_profile(ms))
                modeling_session_submenu.add_command(label = "Assessment Table",
                    command = lambda ms=modeling_session: self.show_assessment_table(ms))
                modeling_session_submenu.add_separator()
                # Adds the alignment submenu to the PyMod main menu.
                label_text = "Modeling Session %s" % (modeling_session.session_id)
                for full_model in modeling_session.full_models:
                    full_model_submenu = Menu(modeling_session_submenu, tearoff = 0)
                    full_model_submenu.add_command(label = "Save to File",
                        command = lambda fm=full_model: self.save_full_model_to_file(fm))
                    full_model_submenu.add_separator()
                    full_model_submenu.add_command(label = "DOPE Profile",
                        command = lambda fm=full_model: self.show_full_model_profile(fm))
                    full_model_submenu.add_command(label = "Assessment Values",
                        command = lambda fm=full_model: self.show_full_model_assessment_values(fm))
                    modeling_session_submenu.add_cascade(label = full_model.model_name, menu = full_model_submenu)
                self.models_menu.add_cascade(label = label_text, menu = modeling_session_submenu)
        else:
            self.models_menu.add_command(label = "There aren't any models")


    def show_session_profile(self, modeling_session):
        """
        Shows a DOPE profile of a modeling session.
        """
        self.show_dope_plot(modeling_session.session_profile)


    def show_assessment_table(self, modeling_session):
        self.show_table(**modeling_session.assessment_table_data)


    def save_full_model_to_file(self, full_model):
        save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = [("PDB","*.pdb")], parent=pymod.main_window)
        if save_file_full_path != "":
            # Moves the saved file to the path chosen by the user.
            try:
                old_path = full_model.original_file_path
                os.rename(old_path, save_file_full_path)
            except:
                title = "File Error"
                message = "Could not save the alignment file to path: %s" % (save_file_full_path)
                self.show_error_message(title, message)


    def show_full_model_profile(self, full_model):
        self.show_dope_plot(full_model.model_profile)


    def show_full_model_assessment_values(self, full_model):
        objfv = full_model.assessment_data[0]
        dopes = full_model.assessment_data[1]
        title = "Assessement Information"
        message = "Assessment Information for %s\n\nObjective Function Value: %s\n\nDOPE Score: %s" % (full_model.model_name, objfv, dopes)
        tkMessageBox.showinfo(title, message, parent=self.main_window)


    ###############################################################################################
    # SIMILARITY SEARCHES.                                                                        #
    ###############################################################################################

    #################################################################
    # Common methods for BLAST and PSI-BLAST.                       #
    #################################################################

    def check_blast_search_selection(self):
        """
        Checks that only one sequence is selected as a query for BLAST and stores the query PyMod
        element in 'self.blast_query_element'.
        """
        correct_selection = False
        selected_sequences = self.get_selected_sequences()
        if len(selected_sequences) == 1:
            # Gets the selected sequence. The main index will be used later to build the cluster.
            self.blast_query_element = selected_sequences[0]
            correct_selection = True
            # Let users decide how to import new sequences when the query is a child element (that
            # is, it is already present in a cluster).
            if self.blast_query_element.is_child:
                new_cluster_text = 'Build a new cluster'
                old_cluster_text = 'Expand old cluster'
                self.blast_search_choices = {new_cluster_text: "extract", old_cluster_text: "expand"}
                self.blast_dialog = Pmw.MessageDialog(self.main_window,
                    title = 'Import new sequences options',
                    message_text = (
                    "Please select how to import in PyMod the new sequences identified in the search:\n\n"+
                    "- Extract the query sequence from its cluster and build a new cluster with\n"+
                    "  the new hit sequences.\n\n"+
                    "- Expand the already existing cluster by appending to it the new hit sequences." ),
                    buttons = (new_cluster_text, old_cluster_text))
                self.blast_dialog.component("message").configure(justify="left")
                self.blast_dialog.configure(command=self.blast_dialog_state)
            else:
                if self.blast_version == "blast":
                    self.blast_state()
                elif self.blast_version == "psi-blast":
                    self.psiblast_state()
        else:
            title = "Selection Error"
            message = "Please select one sequence to perform a PSI-BLAST search"
            self.show_error_message(title, message)
            correct_selection = False

        return correct_selection


    def build_blast_window(self):
        """
        Builds a window containing the widget necessary to define the options for BLAST and
        PSI-BLAST searches.
        """
        current_pack_options = pmgi.pack_options_1
        current_label_options = pmgi.label_style_1

        self.blast_window = pmgi.PyMod_tool_window(self.main_window,
            title = "%s Search Options" % (pmdt.blast_algorithms_dictionary[self.blast_version]),
            upper_frame_title = "Here you can modify search options for %s" % (pmdt.blast_algorithms_dictionary[self.blast_version]),
            submit_command = self.blast_window_state,
            with_frame=True)
        self.blast_window.geometry("550x600")

        # ---
        # Simple options.
        # ---

        # Makes the user chose the folder where the BLAST database files are stored locally.
        if self.blast_version == "psi-blast":
            # A list containing information about the databases present in PyMod BLAST database
            # folder.
            self.list_of_databases_directories = self.build_blast_db_list()
            self.psiblast_database_rds = pmgi.PyMod_radioselect(self.blast_window.midframe, label_text = 'Database Selection')
            # Add the buttons to choose the database.
            # self.psiblast_database_rds.add("select...")
            for db in self.list_of_databases_directories:
                self.psiblast_database_rds.add(db["prefix"])
            # Adds a 'Browse' button in order to let users specify a custom database on their
            # system.
            self.interior = self.psiblast_database_rds.component('frame')
            self.choose_path_label = Label(self.interior, text="None", **pmgi.label_style_2)
            self.choose_path_label.grid(column=3,row=0, padx=(0,0))
            self.psiblast_database_rds.button(0).configure(command=self.choose_psiblast_db_dir)
            # Packs the PSI-BLAST database selection widget.
            self.psiblast_database_rds.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.psiblast_database_rds)

            self.psiblast_iterations_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
                label_text = "PSI-BLAST Iterations",
                label_style = current_label_options,
                value = 3,
                validate = {'validator' : 'integer', 'min' : 1, 'max' : 10} )
            self.psiblast_iterations_enf.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.psiblast_iterations_enf)
            self.blast_window.add_widget_to_validate(self.psiblast_iterations_enf)

        elif self.blast_version == "blast":
            self.ncbiblast_database_rds = pmgi.PyMod_radioselect(self.blast_window.midframe, label_text = 'Database Selection')
            for text, val in pmdt.ncbi_databases:
                self.ncbiblast_database_rds.add(text)
            self.ncbiblast_database_rds.setvalue('Pdb')
            self.ncbiblast_database_rds.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.ncbiblast_database_rds)

        # E-value selection.
        self.e_value_threshold_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "E-value Threshold",
            label_style = current_label_options,
            value = 10.0,
            validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
        self.e_value_threshold_enf.pack(**current_pack_options)
        self.blast_window.add_widget_to_align(self.e_value_threshold_enf)
        self.blast_window.add_widget_to_validate(self.e_value_threshold_enf)

        # Max hit number selection.
        self.max_hits_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "Max Number of Hits",
            label_style = current_label_options,
            value = 100,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000} )
        self.max_hits_enf.pack(**current_pack_options)
        self.blast_window.add_widget_to_align(self.max_hits_enf)
        self.blast_window.add_widget_to_validate(self.max_hits_enf)

        # ---
        # Advanced options.
        # ---
        self.blast_window.show_advanced_button()

        # Minimum id% on with query.
        self.min_id_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "Min ID% Threshold",
            label_style = current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.blast_window.add_widget_to_align(self.min_id_enf)
        self.blast_window.add_advanced_widget(self.min_id_enf)
        self.blast_window.add_widget_to_validate(self.min_id_enf)

        # Minimum coverage on the query.
        self.min_coverage_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "Min Coverage% Threshold",
            label_style = current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.blast_window.add_widget_to_align(self.min_coverage_enf)
        self.blast_window.add_advanced_widget(self.min_coverage_enf)
        self.blast_window.add_widget_to_validate(self.min_coverage_enf)

        # Organisms.
        # To be done.

        # PSI-BLAST minimum inclusion E-value.
        self.psiblast_min_inclusion_eval_default = 0.005
        if self.blast_version == "psi-blast":
            self.psiblast_eval_threshold_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
                label_text = "PSI-BLAST E-value Threshold",
                label_style = current_label_options,
                value = self.psiblast_min_inclusion_eval_default,
                validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
            self.blast_window.add_widget_to_align(self.psiblast_eval_threshold_enf)
            self.blast_window.add_advanced_widget(self.psiblast_eval_threshold_enf)
            self.blast_window.add_widget_to_validate(self.psiblast_eval_threshold_enf)

        # Use current cluster for PSI-BLAST PSSM.
        # if self.blast_query_element.is_child:
        #     self.use_current_pssm_rds = pmgi.PyMod_radioselect(self.blast_window.midframe, label_text = 'Use current cluster as PSSM')
        #     for text in ('Yes', 'No'):
        #         self.use_current_pssm_rds.add(text)
        #     self.use_current_pssm_rds.setvalue('No')
        #     # self.use_current_pssm_rds.pack(side = 'top', padx = 10, pady = 10, anchor="w")
        #     self.blast_window.add_widget_to_align(self.use_current_pssm_rds)
        #     self.blast_window.add_advanced_widget(self.use_current_pssm_rds)

        self.blast_window.align_widgets(input_widget_width=10)


    def blast_window_state(self):
        """
        This function is called when the 'SUBMIT' button in some BLAST main window is pressed.
        """
        # Do not proceed if users have not provided a correct set of input parameters through
        # the GUI.
        if not self.check_blast_input_parameters():
            return False

        # Performs a similarity search with the appropriate program.
        blast_status = None
        if self.blast_version == "blast":
            blast_status = self.run_ncbiblast()
        elif self.blast_version == "psi-blast":
            blast_status = self.run_psiblast()

        # Displays the window with results.
        if blast_status and self.check_blast_results():
            self.blast_window.destroy()
            self.show_blast_output_window()
        else:
            self.blast_window.destroy()

        # Removes temp files at the end of the whole process, both if hits were found or not.
        self.remove_blast_temp_files()

    def check_blast_input_parameters(self):
        """
        Checks if users have provide a set of valid input parameters in order to perform a search.
        """
        # Check if a valid database for PSI-BLAST was provided.
        if self.blast_version == "psi-blast":
            db_full_path = self.get_psiblast_database_from_gui()
            if db_full_path == None:
                title = "Input Error"
                message = "Please choose a valid database."
                self.show_error_message(title, message, parent_window=self.blast_window, refresh=None)
                return False
            if not pmos.verify_valid_blast_dbdir(db_full_path):
                title = "Input Error"
                message = "The database '%s' directory does not contain a valid set of database files." % (db_full_path)
                self.show_error_message(title, message, parent_window=self.blast_window, refresh=None)
                return False
        # Check all the other input fields.
        if not self.check_general_input(self.blast_window):
            return False
        # Returns 'True' only if all input parameters are valid.
        return True


    def check_blast_results(self):
        """
        Checks if at least one hsp was identified in the search and stores the results in
        'self.blast_record'.
        """
        # An attribute where is going to be stored a Biopython "Blast" class object.
        self.blast_record = None
        result_handle = open(os.path.join(self.similarity_searches_directory,self.xml_blast_output_file_name),"r")
        if self.blast_version == "blast":
            self.blast_record = NCBIXML.read(result_handle)
        elif self.blast_version == "psi-blast":
            self.blast_record = self.get_psi_blast_record(result_handle)
        result_handle.close()

        # Filters the BLAST record according to the advanced options.
        if self.blast_window.showing_advanced_widgets:
            for a in self.blast_record.alignments[:]:
                for hsp in a.hsps[:]:
                    # Gets the id% and the query span of the hsp.
                    hspd = self.get_hsp_info(hsp)
                    if hspd["id"]*100 < int(self.min_id) or hspd["query_span"]*100 < int(self.min_coverage):
                        a.hsps.remove(hsp)
                if len(a.hsps) == 0:
                    self.blast_record.alignments.remove(a)

        # Exit the whole process if no hits were found.
        if len(self.blast_record.alignments) == 0:
            self.show_warning_message("PSI-BLAST Message", "No hits weew found for PSI-BLAST for %s." % (self.blast_query_element.my_header))
            return False

        # Returns 'True' if some hits were found.
        return True


    def get_hsp_info(self, hsp, full_query_sequence = None):
        """
        Gets a Biopython HSP object and computes additional information on it and returns it as a
        dictionary.
        """
        # Gets the id% of the hsp.
        matches = float(len(hsp.query) - hsp.gaps)
        idp = float(hsp.identities)/matches

        # Gets the query span.
        qt = len(str(self.blast_query_element.my_sequence).replace("-",""))
        qs = hsp.query_start
        qe = len(str(hsp.query).replace("-","")) + qs
        query_span = float(qe - qs)/float(qt)

        # Gets the subject span.
        hs = hsp.sbjct_start
        he = len(str(hsp.sbjct).replace("-","")) + hs

        additional_infor_dict = {"id": idp, "query_span":query_span, "matches":matches, "query_end":qe, "sbjct_end":he}

        return additional_infor_dict


    def get_blast_output_basename(self):
        basename = (pmos.clean_file_name(self.blast_query_element.get_compact_header()) + "_" +
                    pmdt.blast_algorithms_dictionary[self.blast_version] + "_" +
                    "search_%s" % (self.blast_cluster_counter + 1) )
        return basename


    def remove_blast_temp_files(self):
        output_filename = self.get_blast_output_basename() + ".xml"
        try:
            os.rename(os.path.join(self.similarity_searches_directory, self.xml_blast_output_file_name),
                      os.path.join(self.similarity_searches_directory, output_filename))
            files_to_remove = filter(lambda f: not os.path.splitext(f)[-1] == ".xml", os.listdir(self.similarity_searches_directory))
            map(lambda f: os.remove(os.path.join(self.similarity_searches_directory,f)), files_to_remove)
        except:
            pass

    def show_blast_output_window(self):
        """
        Displays the window with results from BLAST in a new window.
        """
        # BLAST results window.
        self.blast_output_window=Toplevel(self.main_window)
        version_full_name = pmdt.blast_algorithms_dictionary[self.blast_version]
        self.blast_output_window.title("<< %s Output >>" % (version_full_name))
        # Freezes the parent.
        try:
            self.blast_output_window.grab_set()
        except:
            pass
        self.blast_output_window.resizable(1,1)
        self.blast_output_window.geometry('920x520') # '800x320', "920x520"

        # Main frame of the window.
        self.blast_results_main = Frame(self.blast_output_window, background='black')
        self.blast_results_main.pack(expand = YES, fill = BOTH)

        # An upper frame.
        self.blast_results_up_frame = Frame(self.blast_results_main, borderwidth=5, background='black', relief='groove', pady=15)
        self.blast_results_up_frame.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3)
        title_text = "%s Output for: %s\nPlease Select the Sequences to Import" % (version_full_name, self.blast_query_element.get_compact_header())
        self.blast_message = Label(self.blast_results_up_frame, font = "comic 12", height = 1,
                                   text= title_text, background='black', fg='white', pady = 2)
        self.blast_message.pack(ipady=10)

        # A middle scrolled frame.
        self.blast_middleframe = Pmw.ScrolledFrame(
            self.blast_results_main, hull_bg='black', frame_bg='black',
            usehullsize = 0, borderframe = 0, hscrollmode='dynamic',
            vscrollmode='dynamic', hull_borderwidth = 0, clipper_bg='black',)
        self.blast_middleframe.pack(side = TOP, fill = 'both', expand = 1)

        # A frame with some options to choose the hits to import in PyMod.
        self.blast_controls_frame = Frame(self.blast_middleframe.interior(),background='black')
        self.blast_controls_frame.pack(anchor="w")
        self.blast_select_all_button = Button(self.blast_controls_frame,text="Select All", command=self.blast_select_all, **pmgi.button_style_1)
        self.blast_select_all_button.pack(side="left", padx=(30,10),pady=(10,5))
        self.blast_select_none_button = Button(self.blast_controls_frame, text="Select None", command=self.blast_select_none, **pmgi.button_style_1)
        self.blast_select_none_button.pack(side="left", padx=10,pady=(10,5))
        self.blast_select_n_button = Button(self.blast_controls_frame, text="Select Top:", command=self.blast_select_n, **pmgi.button_style_1)
        self.blast_select_n_button.pack(side="left", padx=10,pady=(10,5))
        self.blast_select_n_enf = Pmw.EntryField(self.blast_controls_frame,
                                                 labelpos = None, value = '10',
                                                 validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000})
        self.blast_select_n_enf.component("entry").configure(width = 5)
        self.blast_select_n_enf.pack(side="left", padx=0, pady=(10,5))

        # A frame were the widgets to choose the hits to import are going to be displayed.
        self.blast_ouput_frame = Frame(self.blast_middleframe.interior(), background='black')
        self.blast_ouput_frame.pack(expand=True,fill="both")

        # A bottom frame with a 'SUBMIT' button.
        self.blast_submitframe = Frame(self.blast_results_main, background='black', height=20)
        self.blast_submitframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center", ipadx = 5, ipady = 5)
        self.blast_submit_button=Button(self.blast_submitframe, text="SUBMIT",
            command=self.import_sequences_from_blast, **pmgi.button_style_1)
        self.blast_submit_button.pack(side = BOTTOM, fill=BOTH, anchor=CENTER, pady=10)

        self.display_blast_hits()


    def blast_select_all(self):
        for chk in self.blast_sbjct_checkbuttons_list:
            chk.select()


    def blast_select_none(self):
        for chk in self.blast_sbjct_checkbuttons_list:
            chk.deselect()


    def blast_select_n(self):
        select_top = int(self.blast_select_n_enf.getvalue())
        if select_top != "":
            self.blast_select_none()
            select_top = int(self.blast_select_n_enf.getvalue())
            count = 0
            for chk in self.blast_sbjct_checkbuttons_list:
                chk.select()
                count += 1
                if count == select_top:
                    break


    def display_blast_hits(self, iteration = 1):
        """
        This used inside blast_output_selection to display for each hit some information and a
        checkbutton to select it for importing it inside Pymod.
        """
        header_options = {'background':'black', 'fg':'red', 'height':1, 'pady':10, 'font': 12}
        self.blast_seq_label=Label(self.blast_ouput_frame, text= "Name",**header_options)
        self.blast_seq_label.grid(row=0, column=0, sticky="n")
        self.blast_e_val_label=Label(self.blast_ouput_frame, text= "E-Value", **header_options)
        self.blast_e_val_label.grid(row=0, column=1, sticky="n")
        self.blast_iden_label=Label(self.blast_ouput_frame, text= "Identity", **header_options)
        self.blast_iden_label.grid(row=0, column=2, sticky="n")
        self.query_span_label=Label(self.blast_ouput_frame, text= "Query span", **header_options)
        self.query_span_label.grid(row=0, column=3, sticky="n")
        self.subject_span_label=Label(self.blast_ouput_frame, text= "Subject span", **header_options)
        self.subject_span_label.grid(row=0, column=4, sticky="n")

        # Displays in the results window the hits found in the .xml file (with results from BLAST)
        # that was parsed by Biopython.
        self.blast_output_row = 1
        # This is going to contain the list of values of each checkbutton.
        self.blast_states = []
        # List containing the checkbutton widgets.
        self.blast_sbjct_checkbuttons_list = []

        row_options = {'background':'black', 'fg':'white', 'height':1,'highlightbackground':'black'}

        for alignment in self.blast_record.alignments:
            for hsp in alignment.hsps:
                # Hit info.
                blast_var = IntVar()
                subject_name = alignment.title[:100] + "..." # title[:150]
                chk = Checkbutton(self.blast_ouput_frame,text=subject_name,
                    variable=blast_var, background='black', foreground = "white",selectcolor = "red",
                    height=1, padx=10, highlightbackground='black')
                chk.grid(row=self.blast_output_row, column=0, sticky = "nw", padx=10)
                self.blast_sbjct_checkbuttons_list.append(chk)

                # E-value info.
                evalue=Label(self.blast_ouput_frame, text= "%.2e" % (hsp.expect), **row_options)
                evalue.grid(row=self.blast_output_row, column=1, sticky = "nw", padx=10)

                # HSP identity info.
                hspd = self.get_hsp_info(hsp)
                id_text = str("%s/%s" % (hsp.identities,int(hspd["matches"]))) + str(" (%.1f" % (hspd["id"]*100) + "%)")
                identities=Label(self.blast_ouput_frame, text= id_text, **row_options)
                identities.grid(row=self.blast_output_row, column=2, sticky = "n", padx=10)

                # Query span info.
                span_info_text = "%s-%s (%.1f" % (hsp.query_start, hspd["query_end"], hspd["query_span"]*100) + "%)"
                span_info=Label(self.blast_ouput_frame, text = span_info_text, **row_options)
                span_info.grid(row=self.blast_output_row, column=3, sticky = "n", padx=10)

                # Subject span info.
                hspan_info_text = "%s-%s" % (hsp.sbjct_start, hspd["sbjct_end"])
                hspan_info=Label(self.blast_ouput_frame, text = hspan_info_text, **row_options)
                hspan_info.grid(row=self.blast_output_row, column=4, sticky = "n", padx=10)

                self.blast_output_row += 1
                self.blast_states.append(blast_var)


    def import_sequences_from_blast(self):
        """
        Called when the 'SUBMIT' button is pressed on some BLAST results window.
        """
        # For each hsp takes the state of its tkinter checkbutton.
        self.my_blast_map = map((lambda var: var.get()), self.blast_states)

        # If the user selected some hsp.
        if len(self.my_blast_map) > 0:
            self.build_hits_to_import_list()
            # This will actually import the sequences inside Pymod.
            self.build_blast_cluster()

        self.blast_output_window.destroy()


    def build_hits_to_import_list(self):
        """
        Builds a list containing those hits that were selected by the user in the BLAST results
        window.
        """
        # This will be used to build PyMod elements out of the subjects of the HSP identified by
        # BLAST.
        self.hsp_imported_from_blast = []
        self.total_hsp_counter = 0 # Counts the total number of hsp.
        self.total_fetched_hsp_counter = 0 # Counts the total number of fetched hsp.
        self.total_hit_counter = 0 # Counts the total number of hits.
        self.fetched_hit_counter = 0 # Counts the number of fetched hits.

        for alignment in self.blast_record.alignments:
            hsp_counter = 0 # Counts the number of hsp for a certain hit.
            fetched_hsp_counter = 0 # Counts the number of fetched hsp for a certain hit.
            fetch_hit = False
            for hsp in alignment.hsps:
                fetch_hsp = False
                if self.my_blast_map[self.total_hsp_counter] == 1:
                    fetch_hsp = True
                    fetch_hit = True
                if fetch_hsp:
                    # Appends the hits (subjects).
                    self.hsp_imported_from_blast.append({"hsp":hsp,"title":alignment.title})
                    hsp_counter += 1
                    fetched_hsp_counter += 1
                    self.total_hsp_counter+=1
                    self.total_fetched_hsp_counter += 1
                else:
                    self.total_hsp_counter+=1
            self.total_hit_counter += 1
            if fetch_hit:
                self.fetched_hit_counter += 1


    def get_list_of_aligned_sequences(self, aligned_elements):
        """
        Gets a list of 'PyMod_elements' objects and returns a list of their sequences.
        """
        return [e.my_sequence for e in aligned_elements]


    def build_blast_cluster(self):
        """
        Builds a cluster with the query sequence as a mother and retrieved hits as children.
        """

        self.blast_cluster_counter += 1
        # A new 'PyMod_element' object to represent the new BLAST cluster.
        blast_cluster_element = None
        # This will contain a 'Star_alignment'.
        ba = None

        if self.blast_query_element.is_child:
            if self.new_sequences_import_mode == "extract":
                self.blast_query_element.remove_indels()
                self.extract_child(self.blast_query_element)

            elif self.new_sequences_import_mode == "expand":
                # Updates the cluster element to a new "BLAST search" element.
                blast_cluster_name = "%s cluster %s (query: %s)" % (pmdt.blast_algorithms_dictionary[self.blast_version], self.blast_cluster_counter, self.blast_query_element.get_compact_header())
                blast_cluster_element = self.get_mother(self.blast_query_element)
                blast_cluster_element.my_header = blast_cluster_name
                blast_cluster_element.element_type = "blast-search"
                blast_cluster_element.alignment_object = Alignment("blast-pseudo-alignment", self.blast_cluster_counter)
                self.mark_as_query(self.blast_query_element)
                # Builds a star alignment.
                ba = pmsm.Star_alignment(self.blast_query_element.my_sequence)
                # Preapare the 'Star_alignment' object with sequence already present in the cluster.
                siblings = self.get_siblings(self.blast_query_element)
                list_of_aligned_sequences = self.get_list_of_aligned_sequences(siblings)
                ba.extend_with_aligned_sequences(list_of_aligned_sequences)
                # Adds new hit sequences to the cluster and generate the alignment.
                ba.build_blast_local_alignment_list([h["hsp"] for h in self.hsp_imported_from_blast])
                ba.generate_blast_pseudo_alignment()

        # If the query sequence is a mother, or was a child that was extracted from its cluster.
        if self.blast_query_element.is_mother:
            # Builds the "BLAST search" element. It has the same mother_index of the query.
            blast_cluster_name = "%s cluster %s (query: %s)" % (pmdt.blast_algorithms_dictionary[self.blast_version], self.blast_cluster_counter, self.blast_query_element.get_compact_header())
            blast_cluster_element = PyMod_element(
                "...", blast_cluster_name,
                element_type = "blast-search", adjust_header=False,
                alignment_object = Alignment("blast-pseudo-alignment", self.blast_cluster_counter))
            self.add_element_to_pymod(
                blast_cluster_element, level="mother",
                mother_index=self.blast_query_element.mother_index)
            # Adds the query to the new "BLAST cluster".
            self.add_to_mother(blast_cluster_element, self.blast_query_element)
            self.mark_as_query(self.blast_query_element)

            # Builds a star alignment.
            ba = pmsm.Star_alignment(self.blast_query_element.my_sequence)
            ba.build_blast_local_alignment_list([h["hsp"] for h in self.hsp_imported_from_blast])
            ba.generate_blast_pseudo_alignment()


        # The list of elements whose sequences will be updated according to the star alignment.
        elements_to_update = []
        # Begins with the query element.
        elements_to_update.append(self.blast_query_element)
        if self.blast_query_element.is_child:
            elements_to_update.extend(self.get_siblings(self.blast_query_element))
        # Then creates PyMod elements for all the imported hits and add them to the cluster.
        for h in self.hsp_imported_from_blast:
            # Gives them the query mother_index, to make them its children.
            cs = self.build_pymod_element_from_hsp(h)
            self.add_element_to_pymod(cs, level="child", mother_index=blast_cluster_element.mother_index)
            elements_to_update.append(cs)
        # Updates the sequences according to the BLAST pseudo alignment.
        ba.update_pymod_elements(elements_to_update)

        self.set_initial_ali_seq_number(blast_cluster_element)

        self.gridder()

    def blast_dialog_state(self, dialog_choice):
        self.blast_dialog.withdraw()
        if not dialog_choice:
            return None
        self.new_sequences_import_mode = self.blast_search_choices[dialog_choice]
        if self.blast_version == "blast":
            self.blast_state()
        elif self.blast_version == "psi-blast":
            self.psiblast_state()


    #################################################################
    # Regular BLAST.                                                #
    #################################################################

    def launch_ncbiblast(self):
        """
        Called when BLAST is launched from the main menu.
        """
        self.blast_version = "blast"
        self.check_blast_search_selection()


    def blast_state(self):
        self.build_blast_window()


    def run_ncbiblast(self):
        """
        This function allows to contact the NCBI BLAST server using Biopython.
        """
        self.xml_blast_output_file_name = "blast_out.xml"
        # Actually connects to the server.
        query_seq = str(self.blast_query_element.my_sequence)
        try:
            if self.blast_window.showing_advanced_widgets:
                self.min_id = self.min_id_enf.getvalue()
                self.min_coverage = self.min_coverage_enf.getvalue()
            result_handle = pmsp.qblast("blastp",
                self.get_ncbiblast_database(),
                query_seq,
                hitlist_size=self.max_hits_enf.getvalue(),
                expect=self.e_value_threshold_enf.getvalue())
        except:
            title = "Connection Error"
            message = 'Can not NCBI BLAST server.\nPlease check your Internet access.'
            self.show_error_message(title,message)
            return False
        blast_results = result_handle.read()
        # Saves an XML file that contains the results and that will be used to display them on
        # the results window.
        save_file = open(os.path.join(self.similarity_searches_directory, self.xml_blast_output_file_name), "w")
        save_file.write(blast_results)
        save_file.close()

        # In this way the results window can be opened.
        return True


    def get_ncbiblast_database(self):
        text = self.ncbiblast_database_rds.getvalue()
        for i in pmdt.ncbi_databases:
            if i[0] == text:
                return i[1]


    #################################################################
    # PSI-BLAST.                                                    #
    #################################################################

    def launch_psiblast(self):
        """
        Called when PSI-BLAST is called from the main menu.
        """
        self.blast_version = "psi-blast"
        self.check_blast_search_selection()


    def psiblast_state(self):
        if not self.blast_plus["exe_dir_path"].path_exists():
            self.blast_plus.exe_not_found()
            return False
        self.build_blast_window()


    def run_psiblast(self):
        """
        Launches a standalone version of PSI-BLAST installed locally when using the PSI-BLAST
        option in the plugin main menu.
        """
        # Builds a temporary file with the sequence of the query needed by psiblast.
        query_file_name = "query"
        self.build_sequences_file([self.blast_query_element], query_file_name, file_format="fasta", remove_indels=True, new_directory=self.similarity_searches_directory)

        # Sets some parameters in needed to run PSI-BLAST.
        ncbi_dir = self.blast_plus["exe_dir_path"].get_value()
        db_path = self.get_psiblast_database_from_gui()
        iterations = self.psiblast_iterations_enf.getvalue()
        evalue_cutoff = self.e_value_threshold_enf.getvalue()
        max_hits = self.max_hits_enf.getvalue()
        if self.blast_window.showing_advanced_widgets:
            evalue_inclusion_cutoff = self.psiblast_eval_threshold_enf.getvalue()
            self.min_id = self.min_id_enf.getvalue()
            self.min_coverage = self.min_coverage_enf.getvalue()
        else:
            evalue_inclusion_cutoff = self.psiblast_min_inclusion_eval_default

        # Buids PSI-BLAST command line parameters.
        self.xml_blast_output_file_name = "blast_out.xml"

        try:
            self.execute_psiblast(
                ncbi_dir = ncbi_dir,
                db_path = db_path,
                query = os.path.join(self.similarity_searches_directory, query_file_name+".fasta"),
                inclusion_ethresh = evalue_inclusion_cutoff,
                outfmt = 5,
                out = os.path.join(self.similarity_searches_directory, self.xml_blast_output_file_name),
                num_iterations = iterations,
                evalue = evalue_cutoff,
                max_target_seqs = max_hits)
        except:
            self.show_error_message("PSI-BLAST Error", "There was an error while running PSI-BLAST for %s." % (self.blast_query_element.my_header))
            return False
        # If everything went ok, return 'True', so that the results window can be opened.
        return True


    def execute_psiblast(self, ncbi_dir, db_path, query,
                               inclusion_ethresh=0.001, num_iterations=3,
                               evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        """
        Execute the locally installed PSI-BLAST. Used when running PSI-BLAST through the 'PSI-BLAST'
        command on the plugin main menu or when predicting secondary structures with PSIPRED.
        """
        # Gests the prefix of the database folder.
        moved_to_db_dir = False
        try:
            dp_prefix = pmos.get_blast_database_prefix(db_path)
            # Makes a temporary directory in the folder of the selected database.
            temp_output_dir_name = "__blast_temp__"
            os.mkdir(os.path.join(db_path, temp_output_dir_name))
            # Copies the .fasta file of the query in the temporary folder.
            query_file_name = os.path.split(query)[1]
            shutil.copy(query, os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves in the database directory.
            os.chdir(db_path)
            moved_to_db_dir = True
            # Sets the input and  output file names.
            temp_query_shortcut = os.path.join(temp_output_dir_name, query_file_name)
            temp_out_file_shortcut = None
            if out != None:
                temp_out_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out)[1])
            temp_out_pssm_file_shortcut = None
            if out_pssm != None:
                temp_out_pssm_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out_pssm)[1])

            # Builds the PSI-BLAST commandline.
            psiblast_command = self.build_psiblast_commandline(
                ncbi_dir = ncbi_dir,
                db_path = dp_prefix,
                query = temp_query_shortcut,
                inclusion_ethresh = inclusion_ethresh,
                num_iterations = num_iterations,
                evalue = evalue,
                max_target_seqs = max_target_seqs,
                num_alignments = num_alignments,
                out = temp_out_file_shortcut,
                outfmt = outfmt,
                out_pssm = temp_out_pssm_file_shortcut)

            # Execute PSI-BLAST.
            self.execute_subprocess(psiblast_command)

            # Goes back to the original directory.
            os.chdir(self.current_project_directory_full_path)
            # Removes the query temp file.
            os.remove(os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves the temporary files in the originally specified output directory.
            for output_file in os.listdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.move(os.path.join(db_path, temp_output_dir_name, output_file), os.path.split(query)[0])
            # Remove the temporary directory.
            shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
            return True

        except:
            # If something goes wrong while executing PSI-BLAST, go back to the project directory
            # and removes the temporary directory in the database folder, it it was built.
            if moved_to_db_dir:
                os.chdir(self.current_project_directory_full_path)
            if os.path.isdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
            raise Exception("There was some error while running PSI-BLAST with the input query: %s." % (query))


    def build_psiblast_commandline(self, ncbi_dir, db_path, query,
                                   inclusion_ethresh=0.001, num_iterations=3,
                                   evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        # blastdbcmd -db "\"Users\joeuser\My Documents\Downloads\mydb\"" -info
        # blastdbcmd -db ' "path with spaces/mydb" ' -info
        psiblast_path = pmos.build_commandline_path_string(os.path.join(ncbi_dir, pmos.get_exe_file_name("psiblast")))
        db_path = pmos.build_commandline_file_argument("db", db_path)
        query = pmos.build_commandline_file_argument("query", query)
        inclusion_ethresh = " -inclusion_ethresh %s" % (inclusion_ethresh)
        num_iterations = " -num_iterations %s" % (num_iterations)
        if evalue:
            evalue = " -evalue %s" % (evalue)
        else:
            evalue = ""
        if outfmt:
            outfmt = " -outfmt %s" % (outfmt) # 5 produces an .xml output file.
        else:
            outfmt = ""
        if out:
            out = pmos.build_commandline_file_argument("out", out)
        else:
            out = ""
        if max_target_seqs:
            max_target_seqs = " -max_target_seqs %s" % (max_target_seqs)
        else:
            max_target_seqs = ""
        if out_pssm:
            out_pssm = pmos.build_commandline_file_argument("out_pssm", out_pssm)
        else:
            out_pssm = ""
        if num_alignments:
            num_alignments = " -num_alignments %s" % (num_alignments)
        else:
            num_alignments = ""

        psiblast_command = (psiblast_path + db_path + query +
                            inclusion_ethresh + out + outfmt + out_pssm +
                            num_iterations + evalue + max_target_seqs +
                            num_alignments)

        return psiblast_command


    def build_blast_db_list(self):
        """
        Generates a list of dictionaries each containing information about the sequence databases
        present in the default BLAST database directory.
        """
        blast_db_dir = self.blast_plus["database_dir_path"].get_value()
        list_of_databases_directories = []
        # Information about the database that can be specified by users through the 'Browse'
        # button. This will contain something like:
        # {'prefix': 'swissprot', 'full-path': '/home/user/pymod/databases/swissprot'}
        list_of_databases_directories.append({"full-path":None, "prefix": "browse"})

        # If there are multiple directories containing dabase files with the same prefixes, this
        # will rename their prefixes so that the database radioselect will not have multiple buttons
        # with the same name.
        def get_new_prefix(prefix, list_of_databases_directories, n=1, prefix_root=None):
            if prefix_root == None:
                prefix_root = prefix
            if prefix in [dbd["prefix"] for dbd in list_of_databases_directories]:
                new_prefix = prefix_root + "-" + str(n)
                return get_new_prefix(new_prefix, list_of_databases_directories, n+1, prefix_root)
            else:
                return prefix
        if os.path.isdir(blast_db_dir):
            for path in os.listdir(blast_db_dir):
                full_path = os.path.join(blast_db_dir,path)
                if os.path.isdir(full_path):
                    if pmos.verify_valid_blast_dbdir(full_path):
                        prefix = pmos.get_blast_database_prefix(full_path)
                        prefix = get_new_prefix(prefix, list_of_databases_directories)
                        dbd = {"full-path": full_path, "prefix": prefix}
                        list_of_databases_directories.append(dbd)

        return list_of_databases_directories


    def choose_psiblast_db_dir(self):
        """
        Called when users want to manually choose a BLAST sequence database folder on their system.
        """
        current_path = self.blast_plus["database_dir_path"].get_value()
        new_path = None
        # Lets users choose a new path.
        new_path = askdirectory(title = "Search for a BLAST database directory", initialdir=current_path, mustexist = True, parent = self.blast_window)
        if new_path:
            if pmos.verify_valid_blast_dbdir(new_path):
                prefix = pmos.get_blast_database_prefix(new_path)
                # Updates the label with the new prefix name.
                self.choose_path_label.configure(text=prefix)
                self.list_of_databases_directories[0]["full-path"] = new_path
            else:
                self.choose_path_label.configure(text="None")
                self.list_of_databases_directories[0]["full-path"] = None
                title = "Selection Error"
                message = "The directory you specified does not seem to contain a valid set of sequence database files."
                self.show_error_message(title, message, parent_window = self.blast_window, refresh=False)
        # Selects the 'browse' button once users click on it.
        self.psiblast_database_rds.setvalue("browse")


    def get_psiblast_database_from_gui(self):
        button_name = self.psiblast_database_rds.getvalue()
        for dbd in self.list_of_databases_directories:
            if dbd["prefix"] == button_name:
                return dbd["full-path"]


    def get_psi_blast_record(self,result_handle):
        """
        Convert it to a list because when a using .parse(), Biopython returns a generator.
        """
        records = list(NCBIXML.parse(result_handle))
        return records[0]


    ###############################################################################################
    # ALIGNMENT BUILDING.                                                                         #
    ###############################################################################################

    #################################################################
    # Step 1/4 for performing an alignment from the main menu.      #
    # Methods to launch an alignment program and check if it can be #
    # used (for example, if it is installed on the user's machine). #
    #################################################################

    def launch_regular_alignment_from_the_main_menu(self, program):
        """
        Launched from the 'Sequence' or 'Structure Alignment' menus of thwe main window."
        """
        self.launch_alignment_program(program, "regular-alignment")

    def launch_profile_alignment_from_the_main_menu(self, program):
        """
        Launched from the 'Profile Alignment menu of the main window'.
        """
        self.launch_alignment_program(program, "profile-alignment")


    def launch_alignment_program(self, program, alignment_strategy):
        # This attribute will be used from now on in many other methods that PyMod needs to perform
        # an alignment.
        self.alignment_program = program

        # It can be either "regular-alignment" or "profile-alignment".
        self.alignment_strategy = alignment_strategy

        if self.alignment_program_exists(self.alignment_program):
            self.start_alignment()
        else:
            self.alignment_program_not_found(self.alignment_program)


    def alignment_program_exists(self, alignment_program):
        """
        Returns True if the program full path was specified in the PyMod Options window.
        """
        program_found = False
        if alignment_program == "clustalw":
            if self.clustalw.exe_exists():
                program_found = True
        elif alignment_program == "clustalo":
            if self.clustalo.exe_exists():
                program_found = True
        elif alignment_program == "muscle":
            if self.muscle.exe_exists():
                program_found = True
        elif alignment_program == "salign-seq":
            if self.modeller.can_be_launched():
                program_found = True
        elif alignment_program == "salign-str":
            if self.modeller.can_be_launched():
                program_found = True
        elif alignment_program == "ce":
            if self.ce_exists():
                program_found = True
        else:
            self.unrecognized_alignment_program(alignment_program)
        return program_found


    def ce_exists(self):
        if ce_alignment_mode in ("plugin", "pymol"):
            return True
        else:
            return False


    def alignment_program_not_found(self, alignment_program):
        """
        Displays an error message that tells the user that some program was not found.
        """
        if alignment_program == "clustalw":
            self.clustalw.exe_not_found()
        elif alignment_program == "clustalo":
            self.clustalo.exe_not_found()
        elif alignment_program == "muscle":
            self.muscle.exe_not_found()
        elif alignment_program in ("salign-seq", "salign-str"):
            self.modeller.exe_not_found()
        elif alignment_program == "ce":
            title = "CE-alignment Error"
            message = "CE-alignment is not available on your PyMod installation. If you want to use this function please see CE-alignment installation instructions on PyMod's User Guide."
            self.show_popup_message("error", title, message)
        else:
            self.unrecognized_alignment_program(alignment_program)

    def unrecognized_alignment_program(self,program):
        title = "Alignment error"
        message = "Unrecognized alignment program: %s..." % (program)
        self.show_popup_message("error", title, message)


    def define_alignment_program(self,alignment_program):
        """
        Used at the beginning of a lot of methods below, in order to use the value of
        self.alignment_program if the alignment_program argument of those methods is not specified
        when they are called.
        """
        if alignment_program == None:
            alignment_program = self.alignment_program
        return alignment_program

    #################################################################
    # Step 2/4 for performing an alignment from the main menu.      #
    # Methods to check if the user built a correct selection in     #
    # order to perform an alignment and the start the alignment.    #
    #################################################################

    def start_alignment(self):
        """
        This method will check if there is a correct selection in order to perform an
        alignment, and it will create a window with the alignment options if necessary.
        """

        # A list of all kind of elements (both sequences, alignment and blast-search) that were
        # selected by the user. This is going to be used in other methods too, later in the Pymod
        # alignment process.
        self.selected_elements = []
        # If among the selected sequences there are some leader sequences of some collapsed cluster,
        # ask users if they want to include their hidden siblings in the alignment.
        include_hidden_children_choice = None
        if True in [seq.is_lead_of_collapsed_cluster() for seq in self.get_selected_sequences()]:
            title = "Selection Message"
            message = "Would you like to include in the alignment the hidden sequences of the collapsed clusters?"
            include_hidden_children_choice = tkMessageBox.askyesno(title, message,parent=self.main_window)
        self.selected_elements = self.get_selected_elements(include_hidden_children=include_hidden_children_choice)

        # This will build a series of lists containing informations about which cluster was
        # selected by the user.
        self.build_cluster_lists()

        # List of alignment programs which use a window to let the user choose some of the algorithm
        # options.
        self.alignment_algorithms_with_options = ["clustalw", "clustalo", "ce", "salign-str"]
        # Alignment programs which do not let the user modify some of the algorithm options.
        self.alignment_algorithms_without_options = ["muscle"]
        # Try to assign salign-seq to one of the lists above. If the user is aligning some sequences
        # that have a structure loaded in PyMOL, salign-seq is going to be included in
        # "alignment_algorithms_with_options" beacause the "use structural information to guide the
        # alignment" will be displayed.
        if [e for e in self.get_selected_sequences() if e.has_structure()]:
            self.alignment_algorithms_with_options.append("salign-seq")
        else:
            self.alignment_algorithms_without_options.append("salign-seq")

        # ---
        # For regular alignments.
        # ---
        if self.alignment_strategy == "regular-alignment":
            # First check if the selection is correct.
            if self.check_alignment_selection():
                # Ask if the user wants to proceed with rebuild-entire-old-alignment or extract-siblings
                # if needed.
                if self.check_sequences_level():
                    # Programs that need a window to display their options.
                    if self.alignment_program in self.alignment_algorithms_with_options:
                        self.show_alignment_window()
                    elif self.alignment_program in self.alignment_algorithms_without_options:
                        if self.clusters_are_involved:
                            self.show_alignment_window()
                        else:
                            # Proceeds directly with the alignment whithout showing a window with
                            # alignment options.
                            self.alignment_state()
                    else:
                        self.unrecognized_alignment_program(self.alignment_program)
                else:
                    self.finish_alignment()
            else:
                self.selection_not_valid()

        # ---
        # For profile alignments.
        # ---
        elif self.alignment_strategy == "profile-alignment":
            if self.check_profile_alignment_selection():
                if self.check_sequences_level():
                    self.show_alignment_window()
                else:
                    self.finish_alignment()
            else:
                self.selection_not_valid()


    def build_cluster_lists(self):
        """
        This will build the self.involved_cluster_elements_list, which will contain the elements
        belonging to cluster that were either selected entirely or with at least one selected child.
        """
        # A set that will contain all the mother_indices of the involved clusters (clusters that
        # have at least one selected element).
        self.involved_clusters_mi_list = set()
        # A set that will contain all the mother_indices of the selected clusters (clusters that
        # have all of their sequences selected).
        self.selected_clusters_mi_list = set()
        # A set that will contain all the mother_indices of all the selected childless mothers.
        self.childless_mothers_mi_list = set()

        for e in self.get_selected_elements():
            if e.is_cluster_element() or e.is_child:
                self.involved_clusters_mi_list.add(e.mother_index)
                if e.is_cluster_element():
                    self.selected_clusters_mi_list.add(e.mother_index)
            else:
                self.childless_mothers_mi_list.add(e.mother_index)

        # These are going to be used later. Build PyMod elements lists out of the sets defined
        # above.
        self.involved_cluster_elements_list = []
        for mother_index in sorted(list(self.involved_clusters_mi_list)):
            self.involved_cluster_elements_list.append(self.get_mother_by_index(mother_index))

        self.selected_cluster_elements_list = []
        for mother_index in sorted(list(self.selected_clusters_mi_list)):
            self.selected_cluster_elements_list.append(self.get_mother_by_index(mother_index))


    def check_alignment_selection(self):
        """
        Checks if the elements selected by the user can be aligned in a "regular-alignment".
        """
        correct_selection = False

        if self.alignment_program in pmdt.sequence_alignment_tools:
            # Checks that there are at least two sequences.
            if len(self.selected_elements) > 1:
                correct_selection = True
        elif self.alignment_program in pmdt.structural_alignment_tools:
            # Checks that there are at least two selected elements.
            if len(self.selected_elements) > 1:
                # And that only sequences with structures are selected.
                if not False in [e.has_structure() for e in self.get_selected_sequences()]:
                    correct_selection = True
        else:
            self.unrecognized_alignment_program(self.alignment_program)

        return correct_selection


    def check_sequences_level(self):
        """
        This method is used to ask the user a confirmation before performing an alignment in certain
        situations (for example when building an alignment only with sequences belonging to the same
        cluster).
        """
        proceed_with_alignment = False
        self.clusters_are_involved = False

        # ---
        # For regular alignments.
        # ---
        if self.alignment_strategy == "regular-alignment":

            self.rebuild_single_alignment_choice = False
            self.extract_siblings_choice = False

            if len(self.involved_clusters_mi_list) == 1 and len(self.childless_mothers_mi_list) == 0:
                # If there is only one cluster selected with all its elements: the user might want to
                # rebuild an alignment with all its elements, ask confirmation.
                if self.involved_clusters_mi_list == self.selected_clusters_mi_list:
                    title = "Rebuild alignment?"
                    message = "Would you like to rebuild the alignment with all its sequences?"
                    proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.main_window)
                    self.rebuild_single_alignment_choice = proceed_with_alignment
                else:
                    title = "Extract children?"
                    message = "Would you like to extract the selected children and build a new alignment?"
                    proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.main_window)
                    self.extract_siblings_choice = proceed_with_alignment

            elif len(self.involved_clusters_mi_list) > 0:
                self.clusters_are_involved = True
                proceed_with_alignment = True

            elif len(self.involved_clusters_mi_list) == 0:
                proceed_with_alignment = True

        # ---
        # For profile alignments.
        # ---
        elif self.alignment_strategy == "profile-alignment":
            proceed_with_alignment = True
            self.clusters_are_involved = True

        return proceed_with_alignment


    def check_profile_alignment_selection(self):
        """
        Checks if the selected elements can be used to perform a profile alignment.
        """
        # This will be set to True if there is an adequate selection in order to perform an alignment
        # with at least one profile.
        correct_selection = False
        # This will be set to True if there is an adequate selection in order to align two profiles.
        self.can_perform_ptp_alignment = False

        # Checks if there is at least one cluster which is entirely selected.
        number_of_selected_clusters = len(self.selected_clusters_mi_list)
        number_of_involved_clusters = len(self.involved_clusters_mi_list)
        number_childless_mothers = len(self.childless_mothers_mi_list)

        if number_of_selected_clusters > 0:
            # If there is only one selected cluster.
            if number_of_involved_clusters == 1 and number_of_selected_clusters == 1:
                # Check if there is at least one selected sequence outside the selected cluster.
                if number_childless_mothers > 0:
                    correct_selection = True
            # Two selected clusters.
            elif number_of_involved_clusters == 2:
                # If there aren't any other selected sequences a profile to profile alignment can be
                # performed.
                if number_of_selected_clusters == 2 and number_childless_mothers == 0:
                    self.can_perform_ptp_alignment = True
                    correct_selection = True
                else:
                    correct_selection = True
            # Can a profile to profile alignment be performed?
            elif number_of_involved_clusters >= 3:
                correct_selection = True
        else:
            pass

        return correct_selection


    def check_only_one_selected_child_per_cluster(self,cluster_element):
        """
        Returns True if the cluster element has only one selected child. This is used in
        "check_alignment_joining_selection()" and other parts of the PyMod class (while checking
        the selection for homology modeling).
        """
        if len([child for child in self.get_children(cluster_element) if child.selected]) == 1:
            return True
        else:
            return False


    def check_alignment_joining_selection(self):
        """
        Used to check if there is a right selection in order to perform the Alignment Joiner
        algorithm to join two or more clusters.
        """

        correct_selection = False
        if len(self.involved_cluster_elements_list) > 1:
            # Check that there is only one selected children per cluster.
            too_many_children_per_cluster = False
            for cluster in self.involved_cluster_elements_list:
                if not self.check_only_one_selected_child_per_cluster(cluster):
                    too_many_children_per_cluster = True
                    break

            if too_many_children_per_cluster:
                correct_selection = False
            else:
                correct_selection = True
        else:
            correct_selection = False

        return correct_selection


    def selection_not_valid(self):
        """
        Called to inform the user that there is not a right selection in order to perform an
        alignment.
        """
        title, message = "", ""

        if self.alignment_strategy == "regular-alignment":
            if self.alignment_program in pmdt.sequence_alignment_tools:
                title = "Selection Error"
                message = "Please select two or more sequences for the alignment."
                self.show_error_message(title, message)
            elif self.alignment_program in pmdt.structural_alignment_tools:
                title = "Structures Selection Error"
                message = "Please Select Two Or More\nStructures."
                self.show_error_message(title, message)
            else:
                self.unrecognized_alignment_program(self.alignment_program)

        elif self.alignment_strategy == "profile-alignment":
            title = "Selection Error"
            message = "Please select at least one entire cluster and some other sequences in order to perform a profile alignment."
            self.show_error_message(title, message)


    #################################################################
    # Structure of the windows showed when performing an alignment. #
    #################################################################

    def show_alignment_window(self):
        """
        This method builds the structure of the alignment options window.
        """
        # Builds the window.
        self.alignment_window = pmgi.PyMod_tool_window(self.main_window,
            title = " %s Options " % (pmdt.alignment_programs_full_names_dictionary[self.alignment_program]),
            upper_frame_title = "Here you can modify options for %s" % (pmdt.alignment_programs_full_names_dictionary[self.alignment_program]),
            submit_command = self.alignment_state)
        # Put into the middle frame some options to change the alignment parameters.
        self.build_alignment_window_middle_frame()


    def build_alignment_window_middle_frame(self):
        """
        The middle frame of the window will contain:
            - a frame with widgets to choose the alignment mode.
            - a frame with widgets to change the alignment algorithm parameters.
        """
        # Options to choose the alignment mode.
        if self.clusters_are_involved:
            self.build_alignment_mode_frame()

        # Options to choose the parameters of the alignment algoirthm being used.
        if self.alignment_program in self.alignment_algorithms_with_options:
            self.alignment_options_frame = pmgi.PyMod_frame(self.alignment_window.midframe)
            self.alignment_options_frame.grid(row=1, column=0, sticky = W+E+N+S)
            if self.alignment_program == "clustalw":
                self.build_clustalw_options_frame()
            elif self.alignment_program == "clustalo":
                self.build_clustalo_options_frame()
            elif self.alignment_program == "salign-seq":
                self.build_salign_seq_options_frame()
            elif self.alignment_program == "salign-str":
                self.build_salign_str_options_frame()
            elif self.alignment_program == "ce":
                self.build_ce_align_options_frame()


    ###################################
    # Part for building the a frame   #
    # containing the options for      #
    # choosing the alignment mode.    #
    ###################################

    def build_alignment_mode_frame(self):
        """
        Builds a frame with some options to choose the alignment mode.
        """
        self.alignment_mode_frame = pmgi.PyMod_frame(self.alignment_window.midframe)
        self.alignment_mode_frame.grid(row=0, column=0, sticky = W+E+N+S,pady=(0,10))

        self.alignment_mode_row = 0
        self.alignment_mode_label = Label(self.alignment_mode_frame, font = "comic 12", height = 1,
                    text= "Alignment Mode", background='black', fg='red',
                    borderwidth = 1, padx = 8)
        self.alignment_mode_label.grid(row=self.alignment_mode_row, column=0, sticky = W)

        self.alignment_mode_radiobutton_var = StringVar()

        # Regular alignments.
        if self.alignment_strategy == "regular-alignment":
            self.build_regular_alignment_mode_frame()

        # Profile alignments.
        elif self.alignment_strategy == "profile-alignment":
            self.build_profile_alignment_mode_frame()


    def build_regular_alignment_mode_frame(self):
        # ---
        # Build a new alignment using the selected sequences.
        # ---
        self.alignment_mode_radiobutton_var.set("build-new-alignment")
        self.alignment_mode_row += 1
        new_alignment_rb_text = "Build a new alignment from scratch using the selected sequences."
        self.new_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=new_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="build-new-alignment", background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_build_new_alignment_radio)
        self.new_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))

        # ---
        # Alignment joiner.
        # ---
        # This can be performed only if there is one selected child per cluster.
        if len(self.involved_cluster_elements_list) > 1 and self.check_alignment_joining_selection():
            self.alignment_mode_row += 1
            alignment_joiner_rb_text = "Join the alignments using the selected sequences as bridges (see 'Alignment Joining')."
            self.join_alignments_radiobutton = Radiobutton(self.alignment_mode_frame, text=alignment_joiner_rb_text, variable=self.alignment_mode_radiobutton_var, value="alignment-joining",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_alignment_joiner_radio)
            self.join_alignments_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))

        # ---
        # "Keep previous alignment".
        # ---
        # Right now it can be used only when the user has selected only one cluster.
        # This alignment mode might be used also for multiple clusters, but right now this is
        # prevented in order to keep the alignment modes selection as simple and as intuitive as
        # possible. If the user wants to append to a cluster some sequences that are contained
        # in another cluster using this method, he/she should firtst extract them from their
        # their original cluster. In order to let the user use this option also for multiple
        # clusters, change the == 1 into >= 1 in the below condition.
        if len(self.involved_cluster_elements_list) == 1:
            keep_alignment_rb_text = None
            # Shows a different label for the checkbutton if there is one or more clusters involved.
            if len(self.involved_cluster_elements_list) > 1:
                keep_alignment_rb_text = "Keep only one alignment and align to its selected sequences the remaining ones"
            elif len(self.involved_cluster_elements_list) == 1:
                target_cluster_name = self.involved_cluster_elements_list[0].my_header
                keep_alignment_rb_text = "Keep '%s', and align to its selected sequences the remaining ones." % (target_cluster_name)

            self.alignment_mode_row += 1
            self.keep_previous_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=keep_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="keep-previous-alignment",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',justify=LEFT,anchor= NW, command=self.click_on_keep_previous_alignment_radio)
            self.keep_previous_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))

            # Only if there are multiple clusters involved it displays a combobox to select the
            # target alignment.
            if len(self.involved_cluster_elements_list) > 1:
                # Frame with the options to control the new alignment. It will be gridded in the
                # click_on_keep_previous_alignment_radio() method.
                self.keep_previous_alignment_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.involved_cluster_elements_list, label_text = "Alignment to keep:")


    def build_profile_alignment_mode_frame(self):
        # ---
        # Perform a profile to profile alignment.
        # ---
        if self.can_perform_ptp_alignment:
            self.alignment_mode_radiobutton_var.set("profile-to-profile")
            self.alignment_mode_row += 1
            profile_profile_rb_text = "Profile to profile: perform a profile to profile alignment."
            self.profile_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=profile_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="profile-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_profile_to_profile_radio)
            self.profile_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))

        else:
            self.alignment_mode_radiobutton_var.set("sequence-to-profile")

        # ---
        # Perform sequence to profile alignment.
        # ---
        sequence_profile_rb_text = None
        build_target_profile_frame = False
        # Shows a different label for the checkbutton if there is one or more clusters involved.
        if len(self.selected_cluster_elements_list) > 1:
            sequence_profile_rb_text = "Sequence to profile: align to a target profile the rest of the selected sequences."
            build_target_profile_frame = True
        elif len(self.selected_cluster_elements_list) == 1:
            profile_cluster_name = self.involved_cluster_elements_list[0].my_header
            sequence_profile_rb_text = "Sequence to profile: align the selected sequence to the target profile '%s'." % (profile_cluster_name)

        # Radiobutton.
        self.alignment_mode_row += 1
        self.sequence_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=sequence_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="sequence-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_sequence_to_profile_radio)
        self.sequence_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))

        # If there is more than one selected cluster, then build a frame to let the user choose
        # which is going to be the target profile.
        if build_target_profile_frame:
            # Frame with the options to choose which is going to be the target profile.
            self.target_profile_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.involved_cluster_elements_list, label_text = "Target profile:")
            # If the profile to profile option is available, the "target_profile_frame" will be
            # hidden until the user clicks on the "sequence_to_profile_radiobutton".
            if not self.can_perform_ptp_alignment:
                self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))


    def click_on_build_new_alignment_radio(self):
        if len(self.involved_cluster_elements_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()

    def click_on_alignment_joiner_radio(self):
        if len(self.involved_cluster_elements_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()

    def click_on_keep_previous_alignment_radio(self):
        if len(self.involved_cluster_elements_list) > 1:
            self.keep_previous_alignment_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()

    def click_on_profile_to_profile_radio(self):
        if hasattr(self,"target_profile_frame"):
            self.target_profile_frame.grid_remove()

    def click_on_sequence_to_profile_radio(self):
        if self.can_perform_ptp_alignment:
            self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))


    ###################################
    # Part for building frames with   #
    # algorithm-specific options.     #
    ###################################

    # -----
    # ClustalW.
    # -----
    def build_clustalw_options_frame(self):

        widgets_to_align = []

        # Scoring matrix radioselect.
        self.matrix_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Scoring Matrix Selection')
        self.clustal_matrices = ["Blosum", "Pam", "Gonnet", "Id"]
        self.clustal_matrices_dict = {"Blosum": "blosum", "Pam": "pam", "Gonnet": "gonnet", "Id": "id"}
        for matrix_name in (self.clustal_matrices):
            self.matrix_rds.add(matrix_name)
        self.matrix_rds.setvalue("Blosum")
        self.matrix_rds.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.matrix_rds)

        # Gap open entryfield.
        self.gapopen_enf = pmgi.PyMod_entryfield(
            self.alignment_options_frame,
            label_text = "Gap Opening Penalty",
            value = '10',
            validate = {'validator' : 'integer',
                        'min' : 0, 'max' : 1000})
        self.gapopen_enf.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.gapopen_enf)

        # Gap extension entryfield.
        self.gapextension_enf = pmgi.PyMod_entryfield(
            self.alignment_options_frame,
            label_text = "Gap Extension Penalty",
            value = '0.2',
            validate = {'validator' : 'real',
                        'min' : 0, 'max' : 1000})
        self.gapextension_enf.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.gapextension_enf)

        Pmw.alignlabels(widgets_to_align, sticky="nw")
        pmgi.align_input_widgets_components(widgets_to_align, 10)


    def get_clustalw_matrix_value(self):
        return self.clustal_matrices_dict[self.matrix_rds.getvalue()]

    def get_gapopen_value(self):
        return self.gapopen_enf.getvalue()

    def get_gapextension_value(self):
        return self.gapextension_enf.getvalue()

    # -----
    # Clustal Omega.
    # -----
    def build_clustalo_options_frame(self):
        self.extraoption=Label(self.alignment_options_frame, font = "comic 12",
                           height=1, text="Extra Command Line Option",
                           background='black', fg='red',
                           borderwidth = 1, padx = 8)
        self.extraoption.grid(row=10, column=0, sticky = "we", pady=20)

        self.extraoption_entry=Entry(self.alignment_options_frame,bg='white',width=10)
        self.extraoption_entry.insert(0, "--auto -v")
        self.extraoption_entry.grid(row=10,column=1,sticky="we",
                                    pady=20)

        self.extraoption_def=Label(self.alignment_options_frame, font = "comic 10",
                               height = 1,
                               text= "--outfmt clustal --force",
                               background='black', fg='white',
                               borderwidth = 1, padx = 8)
        self.extraoption_def.grid(row=10,column=2,sticky="we",pady=20)

    # -----
    # SALIGN sequence alignment.
    # -----
    def build_salign_seq_options_frame(self):
        # Use structure information to guide sequence alignment.
        self.salign_seq_struct_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Use structure information')
        for option in ("Yes","No"):
            self.salign_seq_struct_rds.add(option)
        self.salign_seq_struct_rds.setvalue("No")
        self.salign_seq_struct_rds.pack(side = 'top', anchor="w", pady = 10)
        self.salign_seq_struct_rds.set_input_widget_width(10)

    def get_salign_seq_str_alignment_var(self):
        # This method may be called in situations where there aren't sequences with structures
        # among the ones to be aligned and when "salign_seq_struct_alignment_var" is not
        # properly initialized (beacause "build_salign_seq_options_frame" was not called). The
        # condition below is needed to provide 0 value (do not use structural informations) in
        # this kind of situation.
        if "salign-seq" in self.alignment_algorithms_with_options:
            salign_seq_str_alignment_value = pmdt.yesno_dict[self.salign_seq_struct_rds.getvalue()]
        else:
            salign_seq_str_alignment_value = False
        return salign_seq_str_alignment_value

    def build_salign_str_options_frame(self):
        self.compute_rmsd_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Compute RMSD Matrix')
        for option in ("Yes","No"):
            self.compute_rmsd_rds.add(option)
        self.compute_rmsd_rds.setvalue("Yes")
        self.compute_rmsd_rds.pack(side = 'top', anchor="w", pady = 10)
        self.compute_rmsd_rds.set_input_widget_width(10)

    # -----
    # CE alignment.
    # -----
    def build_ce_align_options_frame(self):
        option_widgets_to_align = []

        self.ce_use_seqinfo_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Use Sequence Information')
        for option in ("Yes","No"):
            self.ce_use_seqinfo_rds.add(option)
        self.ce_use_seqinfo_rds.setvalue("No")
        self.ce_use_seqinfo_rds.pack(side = 'top', anchor="w", pady = 10)
        option_widgets_to_align.append(self.ce_use_seqinfo_rds)

        self.compute_rmsd_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Compute RMSD Matrix')
        for option in ("Yes","No"):
            self.compute_rmsd_rds.add(option)
        self.compute_rmsd_rds.setvalue("Yes")
        self.compute_rmsd_rds.pack(side = 'top', anchor="w", pady = 10)
        option_widgets_to_align.append(self.compute_rmsd_rds)

        pmgi.align_set_of_widgets(option_widgets_to_align, input_widget_width=10)


    def get_ce_align_use_seqinfo_value(self):
        return pmdt.yesno_dict[self.ce_use_seqinfo_rds.getvalue()]

    #################################################################
    # Step 3/4 for performing an alignment from the main menu.      #
    # Methods to launch an alignment and to update the sequences    #
    # loaded in PyMod once the alignment is complete.               #
    #################################################################

    def alignment_state(self):
        """
        This method is called either by the "start_alignment()" method or when the 'SUBMIT' button
        in some alignment window is pressed. It will first define the alignment mode according to
        the choices made by the user. Then, depending on the alignment strategy and the alignment
        mode, it will execute all the steps necessary to perform the alignment.
        """
        # Gets the parameters from the GUI in order to chose the kind of alignment to perform.
        self.alignment_mode = self.define_alignment_mode()

        # This list is going to be used inside
        #     - "create_alignment_element"
        #     - "update_aligned_sequences"
        #     - "align_and_keep_previous_alignment"
        #     - "perform_alignment_joining"
        #     - "salign_sequence_profile_alignment"
        #     - "profile_profile_alignment"
        # methods.
        self.elements_to_align = []

        if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            self.elements_to_align = self.get_selected_sequences()
            self.perform_alignment(self.elements_to_align)

        elif self.alignment_mode == "keep-previous-alignment":
            self.elements_to_align = [] # It will be populated inside self.align_and_keep_previous_alignment().
            self.align_and_keep_previous_alignment()

        elif self.alignment_mode == "alignment-joining":
            self.elements_to_align = self.get_selected_sequences()
            self.perform_alignment_joining()

        elif self.alignment_mode == "sequence-to-profile":
            self.elements_to_align = []
            self.perform_sequence_to_profile_alignment()

        elif self.alignment_mode == "profile-to-profile":
            self.elements_to_align = []
            self.profile_profile_alignment()

        else:
            title = "Alignment Error"
            message = "Unrecognized alignment mode: %s" % (self.alignment_mode)
            self.show_error_message(title,message)
            return

        self.create_alignment_element()
        self.update_aligned_sequences()
        self.finish_alignment()


    def define_alignment_mode(self):
        """
        Gets several parameters from the user interface in order to define the alignment mode.
        """
        alignment_mode = None

        # ---
        # Regular alignments.
        # ---
        if self.alignment_strategy == "regular-alignment":

            if self.rebuild_single_alignment_choice:
                alignment_mode = "rebuild-old-alignment"

            elif self.extract_siblings_choice:
                alignment_mode = "build-new-alignment"

            elif self.clusters_are_involved:
                # It can be either "rebuild-old-alignment" or "keep-previous-alignment".
                alignment_mode = self.alignment_mode_radiobutton_var.get()
                # Takes the index of the target cluster.
                self.target_cluster_index = None

                # Takes the index of the target cluster for the "keep-previous-alignment" mode.
                if alignment_mode == "keep-previous-alignment":
                    # If there is only one cluster involved its index its going to be 0.
                    if len(self.involved_cluster_elements_list) == 1:
                        self.target_cluster_index = 0 # Cluster index.
                    # Get the index of the cluster from the combobox.
                    elif len(self.involved_cluster_elements_list) > 1:
                        target_cluster_name = self.keep_previous_alignment_frame.get_selected_cluster()
                        self.target_cluster_index = self.keep_previous_alignment_frame.get_selected_cluster_index(target_cluster_name)

            else:
                alignment_mode = "build-new-alignment"

        # ---
        # Profile alignments.
        # ---
        elif self.alignment_strategy == "profile-alignment":
            # It can be either "sequence-to-profile" or "profile-to-profile".
            alignment_mode = self.alignment_mode_radiobutton_var.get()
            # Takes the index of the target cluster.
            self.target_cluster_index = None

            # Takes the index of the target cluster for the "keep-previous-alignment" mode.
            if alignment_mode == "sequence-to-profile":
                # If there is only one cluster involved its index its going to be 0.
                if len(self.selected_cluster_elements_list) == 1:
                    self.target_cluster_index = 0 # Cluster index.
                # Get the index of the cluster from the combobox.
                elif len(self.selected_cluster_elements_list) > 1:
                    target_cluster_name = self.target_profile_frame.get_selected_cluster()
                    self.target_cluster_index = self.target_profile_frame.get_selected_cluster_index(target_cluster_name)

        return alignment_mode


    ###################################
    # Builds a PyMod cluster element  #
    # that will contain as children   #
    # all the aligned sequences.      #
    ###################################

    def create_alignment_element(self):
        """
        A method to create a PyMod element for the alignment and to build a cluster to contain the
        aligned sequences.
        """
        # ---
        # Build a new alignment.
        # ---
        if self.alignment_mode == "build-new-alignment":

            # If this is set to 0, the new cluster containing the newly aligned sequences will be
            # placed at the top of the element list in the PyMod main window.
            # If it is set to None, the new cluster will be placed at the bottom of the list.
            lowest_mother_index = None

            # Place the new cluster in the list of elements in the main window at the same point of
            # the aligned sequence with the lowest mother_index.
            # Finds the mother index that the new alignment element is going to have.
            mothers_list = [e for e in self.elements_to_align if e.is_mother and not e.is_cluster_element()]
            children_list = [e for e in self.elements_to_align if e.is_child]
            if mothers_list != []:
                lowest_mother_index = min([e.mother_index for e in mothers_list]) # CHECK! # min(mothers_list,key=lambda el: el.mother_index).mother_index
            # If there are only children, build the new alignment where the child with the lower
            # index was.
            else:
                lowest_mother_index = min([e.mother_index for e in children_list])# CHECK! # min(children_list,key=lambda el: el.mother_index).mother_index

            # Actually creates the new PyMod alignment element.
            self.alignment_count += 1
            ali_name = self.set_alignment_element_name(pmdt.alignment_programs_full_names_dictionary[self.alignment_program],self.alignment_count)
            ali_object = self.build_alignment_object(self.alignment_program, self.alignment_count)
            self.alignment_element = PyMod_element("...", ali_name, element_type="alignment", alignment_object=ali_object, adjust_header=False)
            self.add_element_to_pymod(self.alignment_element, "mother", mother_index = lowest_mother_index)

            # Move all the elements in the new cluster.
            for element in sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
                self.add_to_mother(self.alignment_element,element)

        # ---
        # Rebuilds an old alignment.
        # ---
        elif self.alignment_mode == "rebuild-old-alignment":
            self.alignment_element = self.involved_cluster_elements_list[0]
            old_id = self.alignment_element.alignment.id
            self.alignment_element.alignment = self.build_alignment_object(self.alignment_program, old_id)
            if self.alignment_element.element_type == "alignment":
                self.alignment_element.my_header = self.set_alignment_element_name(pmdt.alignment_programs_full_names_dictionary[self.alignment_program], old_id)
            elif self.alignment_element.element_type == "blast-search":
                old_cluster_name = self.alignment_element.my_header
                self.alignment_element.my_header = self.updates_blast_search_element_name(old_cluster_name, pmdt.alignment_programs_full_names_dictionary[self.alignment_program])

        # ---
        # Expand an already existing cluster with new sequences.
        # ---
        elif self.alignment_mode in ("keep-previous-alignment", "sequence-to-profile"):

            # Gets the target cluster element.
            self.alignment_element = None

            if self.alignment_mode == "keep-previous-alignment":
                self.alignment_element = self.involved_cluster_elements_list[self.target_cluster_index]
            elif self.alignment_mode == "sequence-to-profile":
                self.alignment_element = self.selected_cluster_elements_list[self.target_cluster_index]

            # Appends new sequences to the target cluster.
            for element in self.elements_to_add:
                self.add_to_mother(self.alignment_element,element)

            # Updates the alignment element with new information about the new alignment.

            # Creates an alignment object with the same id of the alignment that was kept.
            ali_to_keep_id = self.alignment_element.alignment.id
            self.alignment_element.alignment = self.build_alignment_object("merged", ali_to_keep_id)

            alignment_description = None
            if self.alignment_mode == "keep-previous-alignment":
                alignment_description = "merged with %s" % (pmdt.alignment_programs_full_names_dictionary[self.alignment_program])
            elif self.alignment_mode == "sequence-to-profile":
                alignment_description = "built with sequence-to-profile with %s" % (pmdt.alignment_programs_full_names_dictionary[self.alignment_program])
            alignment_description = "merged"

            # Changes the name of the alignment element.
            if self.alignment_element.element_type == "alignment":
                self.alignment_element.my_header = self.set_alignment_element_name(alignment_description,ali_to_keep_id)
            elif self.alignment_element.element_type == "blast-search":
                # self.alignment_element.my_header = self.set_alignment_element_name(alignment_description,ali_to_keep_id)
                pass

        # ---
        # Join two or more existing clusters.
        # ---
        elif self.alignment_mode in ("alignment-joining", "profile-to-profile"):

            # Find the right mother index in order to build the new cluster where one of the
            # original ones was placed.
            lowest_mother_index = 0
            mothers_list = [e for e in self.selected_elements[:]+self.involved_cluster_elements_list[:] if e.is_mother]
            # Use min or max to build the new cluster respectively where the top o bottom original
            # cluster were.
            lowest_mother_index = min([e.mother_index for e in mothers_list]) # CHECK! # min(mothers_list,key=lambda el: el.mother_index).mother_index

            # Build a new "Alignment" class object.
            self.alignment_count += 1

            alignment_description = None
            if self.alignment_mode == "alignment-joining":
                alignment_description = "joined by using " + pmdt.alignment_programs_full_names_dictionary[self.alignment_program]
            elif self.alignment_mode == "profile-to-profile":
                alignment_description = "joined by using " + pmdt.alignment_programs_full_names_dictionary[self.alignment_program] + "profile to profile alignment"
            alignment_description = "joined by using " + pmdt.alignment_programs_full_names_dictionary[self.alignment_program]

            ali_name = "Joined " + self.set_alignment_element_name(alignment_description, self.alignment_count)
            ali_object = self.build_alignment_object(self.alignment_program+"-joined", self.alignment_count)

            # Builds the new "PyMod_element" object for the new alignment.
            self.alignment_element = PyMod_element("...", ali_name, element_type="alignment", alignment_object=ali_object, adjust_header=False)
            self.add_element_to_pymod(self.alignment_element, "mother", mother_index=lowest_mother_index)

            # Move all the sequences in the new cluster.
            new_elements = []
            bridges_list = []
            # First appends the mothers (if any) to the new cluster.
            for e in self.selected_elements:
                if e.is_mother and not e.is_cluster_element():
                    new_elements.append(e)
                    bridges_list.append(e)
            # Then appends the children.
            for cluster in self.involved_cluster_elements_list:
                for c in self.get_children(cluster):
                    new_elements.append(c)
                    if c.selected:
                        bridges_list.append(c)
            for element in sorted(new_elements,key=lambda el: (el.mother_index,el.child_index)):
                self.add_to_mother(self.alignment_element,element)

            # Marks the bridges so that they are displayed with a "b" in their cluster.
            if self.alignment_mode == "alignment-joining":
                for b in bridges_list:
                    b.is_bridge = True

        self.set_initial_ali_seq_number(self.alignment_element)

    def set_initial_ali_seq_number(self, cluster_element):
        number_of_seqs = len(self.get_children(cluster_element))
        cluster_element.alignment.initial_number_of_sequence = number_of_seqs


    def set_alignment_element_name(self, alignment_description, alignment_id="?"):
        """
        Builds the name of a new alignment element. This name will be displayed on PyMod main
        window.
        """
        alignment_name = "Alignment " + str(alignment_id) + " (%s)" % (alignment_description)
        return alignment_name


    def updates_blast_search_element_name(self, old_cluster_name, alignment_program, alignment_id="?"):
        new_name = old_cluster_name # old_cluster_name.rpartition("with")[0] + "with %s)" % (alignment_program)
        return new_name


    def build_alignment_object(self,alignment_program, alignment_id="?"):
        """
        This method will build an "Alignment" class object when a new alignment is performed.
        """
        # Builds the new "Alignment" class object.
        alignment_object = Alignment(alignment_program, alignment_id)

        # For certain alignment programs builds a .dnd or .tree file that can be used to display
        # trees from the "Alignment" menu of the PyMod main window.
        if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment") and len(self.elements_to_align) > 2:
            # Clustal programs produce guide trees in newick format.
            if self.alignment_program in ("clustalw", "clustalo"):

                temp_dnd_file_path = os.path.join(self.alignments_directory, self.current_alignment_file_name+".dnd")
                new_dnd_file_path = os.path.join(self.alignments_directory,self.alignments_files_names+str(alignment_id)+"_guide_tree"+".dnd")
                shutil.copy(temp_dnd_file_path, new_dnd_file_path)

                # ClustalO produces a .dnd file without changing the ":" characters in the name of the
                # PDB chains and this gives problems in displaying the names when using Phylo. So the
                # ":" characters have to be changed in "_".
                if self.alignment_program == "clustalo":
                    old_dnd_file = open(new_dnd_file_path,"rU")
                    new_dnd_file_content = ''
                    for dnd_item in old_dnd_file.readlines():
                        if re.search(r"_Chain\:?\:",dnd_item):
                            Chain_pos=dnd_item.find("_Chain:")+6
                            dnd_item=dnd_item[:Chain_pos]+'_'+dnd_item[Chain_pos+1:]
                        new_dnd_file_content+=dnd_item
                    old_dnd_file.close()
                    new_dnd_file = open(new_dnd_file_path,"w")
                    new_dnd_file.write(new_dnd_file_content)
                    new_dnd_file.close()
                alignment_object.set_dnd_file_path(new_dnd_file_path)

            # SALIGN algorithms produce dendrograms.
            elif self.alignment_program.startswith("salign"):
                temp_tree_file_path = os.path.join(self.alignments_directory, self.current_alignment_file_name+".tree")
                new_tree_file_path = os.path.join(self.alignments_directory,self.alignments_files_names+str(alignment_id)+"_dendrogram"+".tree")
                shutil.copy(temp_tree_file_path, new_tree_file_path)
                alignment_object.set_dnd_file_path(new_tree_file_path)

        return alignment_object


    def compute_rmsd_list(self, aligned_elements):
        rmsd_list =  {}
        for i, ei in enumerate(self.elements_to_align):
            for j, ej in enumerate(self.elements_to_align):
                if j > i:
                    rmsd = self.get_rmsd(ei,ej)
                    # This will fill "half" of the matrix.
                    rmsd_list.update({(ei.unique_index, ej.unique_index): rmsd})
                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    rmsd_list.update({(ej.unique_index, ei.unique_index): rmsd})
                if j == i:
                    rmsd_list.update({(ei.unique_index, ej.unique_index): 0.0})
        return rmsd_list


    def get_rmsd(self, element_1, element_2):
        """
        Takes two 'PyMod_elements' objects and computes a RMSD deviation between their structures
        loaded in PyMOL. The RMSD is computed between a list of residues pairs defined by the
        alignment currently existing in PyMod between the two sequences.
        """
        list_of_matching_ids_1 = []
        list_of_matching_ids_2 = []
        ali_id = 0
        for id_1, id_2 in zip(element_1.my_sequence, element_2.my_sequence):
            if id_1 != "-" and id_2 != "-":
                list_of_matching_ids_1.append(element_1.get_pdb_index(ali_id))
                list_of_matching_ids_2.append(element_2.get_pdb_index(ali_id))
            ali_id += 1

        objsel_1 = element_1.build_chain_selector_for_pymol()
        objsel_2 = element_2.build_chain_selector_for_pymol()
        list_of_distances = []

        for resid_1, resid_2 in zip(list_of_matching_ids_1, list_of_matching_ids_2):
            res1_arg = "object %s and n. CA and i. %s" % (objsel_1, resid_1)
            res2_arg = "object %s and n. CA and i. %s" % (objsel_2, resid_2)
            d = 0.0
            try:
                d = cmd.get_distance(res1_arg, res2_arg)
            except:
                print "# ERROR!"
                print res1_arg, res2_arg
            list_of_distances.append(d)

        # Remove outliers: sometimes CE-align aligns residues that, even if actually homologous,
        # are found distant from each other, such as residues in proteins' flexible N- or
        # C-terminus.
        """
        from scipy import stats
        n = len(list_of_distances)
        mean = numpy.mean(list_of_distances)
        std = numpy.std(list_of_distances)
        for d in list_of_distances[:]:
            tval = (d - mean)/std
            pval = stats.t.sf(numpy.abs(tval), n-1)*2
            remove = "keep"
            if pval*n <= 0.5:
                list_of_distances.remove(d)
                remove = "remove"
            print 't-val = %6.3f p-val = %6.4f, %s' % (tval, pval, remove)
        """
        for d in list_of_distances[:]:
            if d >= 6.5:
                list_of_distances.remove(d)
        rmsd = numpy.sqrt(numpy.sum(numpy.square(list_of_distances))/len(list_of_distances))

        return rmsd


    ###################################
    # Updates the sequences of the    #
    # aligned elements and then       #
    # display them in PyMod.          #
    ###################################

    def update_aligned_sequences(self, remove_temp_files=True):
        """
        Called when an alignment is performed. It updates the sequences with the indels obtained in the
        alignment. And also deletes the temporary files used to align the sequences.
        """

        if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):

            # Sequence alignment tools have an output file that can be easily used by the
            # "display_ordered_sequences()" moethod.
            if self.alignment_program in pmdt.sequence_alignment_tools:
                self.display_ordered_sequences()

            # Structural alignments tools might have output files which need the
            # "display_hybrid_al" method in order to be displayed in the PyMod main window.
            elif self.alignment_program in pmdt.structural_alignment_tools:
                if self.alignment_program == "ce":
                    if len(self.elements_to_align) == 2:
                        self.display_ordered_sequences()
                    elif len(self.elements_to_align) > 2:
                        self.display_hybrid_al(self.current_alignment_file_name)
                elif self.alignment_program == "salign-str":
                    self.display_ordered_sequences()
                # Add information to build a root mean square deviation matrix. These RMSD will be
                # computed only once, when the structural alignment first built.
                if pmdt.yesno_dict[self.compute_rmsd_rds.getvalue()]:
                    rmsd_list = self.compute_rmsd_list(self.elements_to_align)
                    self.alignment_element.alignment.set_rmsd_list(rmsd_list)

        elif self.alignment_mode in ("keep-previous-alignment", "sequence-to-profile"):
            self.display_hybrid_al()

        elif self.alignment_mode in ("alignment-joining", "profile-to-profile"):
            self.display_hybrid_al()

        if remove_temp_files:
            self.remove_alignment_temp_files()

        self.gridder()


    def display_ordered_sequences(self):
        """
        Some alignments programs will produce an output file with the aligned sequences placed
        in the same order of the "elements_to_align" list. This method takes the newly aligned
        sequence from these output files and updates the sequences in the "elements_to_align".
        """

        # Gets from an alignment file the sequences with their indels produced in the alignment.
        handle = open(os.path.join(self.alignments_directory, self.current_alignment_file_name+".aln"), "rU")
        records = list(SeqIO.parse(handle, "clustal"))
        handle.close()

        # Updates the Sequences.
        for a,element in enumerate(self.elements_to_align):
            element.my_sequence = self.correct_sequence(str(records[a].seq))


    def display_hybrid_al(self,alignment_file_name="al_result"):
        """
        Actually updates the sequences aligned by using the "alignments_joiner()" method.
        """
        try:
            output_file_path = os.path.join(self.alignments_directory, alignment_file_name +".txt")
            f=open(output_file_path, "r")
            for element in range(len(self.pymod_elements_list)):
                for line in f:
                    # if self.pymod_elements_list[element].my_header_fix[:12].replace(":", "_")==line.split()[0][:12].replace(":", "_"):
                    if self.pymod_elements_list[element].my_header[:12].replace(":", "_")==line.split()[0][:12].replace(":", "_"):
                        new_sequence = self.correct_sequence(line.split()[1])
                        self.pymod_elements_list[element].my_sequence = new_sequence
                f.seek(0)
            f.close()
        except Exception,e:
            self.general_error(e)


    def remove_alignment_temp_files(self):
        """
        Used to remove the temporary files produced while performing an alignment.
        """
        def check_file_to_keep(file_basename):
            file_name = os.path.splitext(file_basename)[0]
            # The only files that are not going to be deleted are guide tree or tree files generated
            # from an alignment. They will be kept in order to be accessed by users who wants to
            # inspect the trees. Their names are going to be built like the following:
            #     (self.alignments_files_name) + (alignment_id) + ("_guide_tree" or "_align_tree")
            # resulting in:
            #     alignment_n_guide_tree.dnd or alignment_n_guide_tree.dnd
            if file_name.startswith(self.alignments_files_names) and (file_name.endswith("guide_tree") or file_name.endswith("align_tree") or file_name.endswith("dendrogram")):
                return False
            else:
                return True

        files_to_remove = filter(lambda file_basename: check_file_to_keep(file_basename), os.listdir(self.alignments_directory))
        for file_basename in files_to_remove:
            file_path_to_remove = os.path.join(self.alignments_directory,file_basename)
            os.remove(file_path_to_remove)


    ###################################
    # Finish the alignment.           #
    ###################################

    def finish_alignment(self):
       try:
           self.alignment_window.destroy()
       except:
           pass


    #################################################################
    # Step 4/4 for performing an alignment from the main menu.      #
    # Methods to perform different "regular" alignment modes.       #
    #################################################################

    ###################################
    # Methods to perform a regular    #
    # alignment.                      #
    ###################################

    def perform_alignment(self, sequences_to_align, alignment_name=None, alignment_program=None, use_parameters_from_gui=True):
        """
        Perform a new sequence (or structural) alignment with the algorithm provided in the
        "alignment_program" argument. This method can be used in other parts of the plugin
        independently of the whole process initiated when performing an alignment using the commands
        in the 'Tools' menu in the PyMod main menu.
        """

        # Generate a name for the current alignment, which will be used to name its files. If the
        # "name" argument is set to "None", the "set_current_alignment_file_name" will automatically
        # generate a name for it.
        self.current_alignment_file_name = self.set_current_alignment_file_name(alignment_name)
        # The same goes for the alignment program.
        alignment_program = self.define_alignment_program(alignment_program)

        if alignment_program == "clustalw":
            if use_parameters_from_gui:
                self.run_clustalw(sequences_to_align,
                              alignment_name=self.current_alignment_file_name,
                              matrix=self.get_clustalw_matrix_value(),
                              gapopen=int(self.get_gapopen_value()),
                              gapext=float(self.get_gapextension_value()) )
            else:
                self.run_clustalw(sequences_to_align, alignment_name=self.current_alignment_file_name)

        elif alignment_program == "clustalo":
            if use_parameters_from_gui:
                self.run_clustalo(sequences_to_align, extraoption=self.extraoption_entry.get(),
                              alignment_name=self.current_alignment_file_name)
            else:
                self.run_clustalo(sequences_to_align, alignment_name=self.current_alignment_file_name)

        elif alignment_program == "muscle":
            self.run_muscle(sequences_to_align, alignment_name=self.current_alignment_file_name)

        elif alignment_program == "salign-seq":
            if use_parameters_from_gui:
                self.salign_malign(sequences_to_align,
                                   alignment_name=self.current_alignment_file_name,
                                   use_structural_information= self.get_salign_seq_str_alignment_var() )
            else:
                self.salign_malign(sequences_to_align,
                                   alignment_name=self.current_alignment_file_name,
                                   use_structural_information=False)

        elif alignment_program == "salign-str":
            self.salign_align3d(sequences_to_align, alignment_name=self.current_alignment_file_name)

        elif alignment_program == "ce":
            if use_parameters_from_gui:
                self.run_ce_alignment(sequences_to_align,
                                      ce_alignment_file_name=self.current_alignment_file_name,
                                      use_seq_info=self.get_ce_align_use_seqinfo_value())
            else:
                self.run_ce_alignment(sequences_to_align, ce_alignment_file_name=self.current_alignment_file_name)

        else:
            self.unrecognized_alignment_program(self.alignment_program)


    ###################################
    # Methods for the "keep previous  #
    # alignment" mode.                #
    ###################################

    def align_and_keep_previous_alignment(self):
        """
        Align all selected elements to some cluster. Briefly, what it does is the following:
            - perform a multiple alignment between all the selected sequences in the target cluster
              and the other selected sequences to add to the alignment, using the algorithm chosen
              by the user.
            - for each of the sequences to add, find the sequence in the cluster which is less
              distant from it (in sequence alignments in terms of sequence identity) and estabilish
              conserved pairs.
            - align individually each conserved pair, and the merge the alignments with the original
              alignment of the target cluster.
        This mode is useful when the user is manually building an alignment and wants to append
        some sequences to some cluster by aligning them to a specific sequence in the target
        alignment.
        """

        # Gets the target cluster element (the alignment that has to be kept).
        target_cluster_element = self.involved_cluster_elements_list[self.target_cluster_index]
        # List of the sequences elements that belong to the target cluster.
        alignment_to_keep_elements = self.get_children(target_cluster_element)
        # List of the selected sequence in the target cluster.
        self.selected_sequences_in_target_alignment = [e for e in alignment_to_keep_elements if e.selected]

        # Checks if the there are multiple selected sequence in the target cluster.
        multiple_selected_seq_in_target_alignment = False
        if len(self.selected_sequences_in_target_alignment) > 1:
            multiple_selected_seq_in_target_alignment = True

        # List of the selected sequences that have to be appended to the target cluster.
        self.elements_to_add = []
        for e in self.selected_elements:
            if not e.is_cluster_element() and not e in alignment_to_keep_elements:
                self.elements_to_add.append(e)

        # ---
        # Perform a first alignment between all the selected sequences (belonging to the target
        # cluster and external).
        # ---

        self.initial_alignment_name = "all_temporary"

        self.elements_to_align = self.selected_sequences_in_target_alignment[:]+self.elements_to_add[:]

        # For sequence alignment algorithms, perform the first multiple alignment with the same
        # algorithtm.
        if self.alignment_program in pmdt.sequence_alignment_tools:
            self.perform_alignment(
                self.elements_to_align,
                alignment_name=self.initial_alignment_name,
                alignment_program=None,
                use_parameters_from_gui=True)

        # For structural alignment algorithms, perform the first multiple alignment with a sequence
        # alignment algorithm with default parameters.
        elif self.alignment_program in pmdt.structural_alignment_tools:
            self.perform_alignment(
                self.elements_to_align,
                alignment_name=self.initial_alignment_name,
                alignment_program="clustalw",
                use_parameters_from_gui=False)

        # ---
        # Generate the highest identity pair list.
        # ---
        highest_identity_pairs_list = self.generate_highest_identity_pairs_list(self.initial_alignment_name)

        # List of filenames of the pairwise alignments of each sequence from "elements_to_add" to
        # most similiar sequence in the "selected_sequences_in_target_alignment".
        self.highest_identity_pairs_alignment_list=[]

        # Performs the alignments and stores the name of the output files names (they will be .aln
        # files) in the list above.
        self.highest_identity_pairs_alignment_list = self.align_highest_identity_pairs_list(highest_identity_pairs_list)

        # ---
        # Actually joins all the alignments.
        # ---

        # First builds the al_result.txt file with the target alignment, this is needed by
        # "alignments_joiner()" method used below.
        self.merged_alignment_output = "al_result" # align_output.txt
        self.build_sequences_file(alignment_to_keep_elements, self.merged_alignment_output,
                                  file_format="pymod", remove_indels=False)

        # Performs the alignments joining progressively.
        for comp in self.highest_identity_pairs_alignment_list:
            self.alignments_joiner(
                os.path.join(self.alignments_directory, self.merged_alignment_output + ".txt"),
                os.path.join(self.alignments_directory, comp + ".aln"))

        # The temporary files needed to peform this alignment will be deleted inside the
        # update_aligned_sequences() method.


    def generate_highest_identity_pairs_list(self, initial_alignment_name):
        """
        For each sequence to add to the alignment, finds the nearest selected sequence (in terms
        of sequence identity) of the target cluster according to the information of previous
        multiple alignment between all the sequences.
        """
        # Reads the output file of the alignment and stores  in a variable a list of its biopython
        # record objects.
        initial_alignment_file = open(os.path.join(self.alignments_directory, initial_alignment_name + ".aln"), "rU")
        initial_alignment_records = list(SeqIO.parse(initial_alignment_file, "clustal"))
        initial_alignment_file.close()

        # A list that is going to contain as many rows as the sequence to add to the alignment and
        # as many columns as the selected sequences in target alignment.
        pair_list=[]
        index = 0
        # Parses the records in the fasta file with the initial alignment just generated.
        for element in initial_alignment_records:
            for sequence in self.elements_to_add:
                # If the sequence in the list is the same of some element in the records.
                if (element.id[:12] == sequence.my_header[:12] or
                    element.id[:12] == sequence.my_header[:12].replace(':', '_')):
                    pair_list.append([])
                    # Parses the list for sequences fo the alignment to keep.
                    for structure in initial_alignment_records:
                        for struct in self.selected_sequences_in_target_alignment:
                            if (structure.id[:12] == struct.my_header.replace(':', '_')[:12] or
                                structure.id[:12] == struct.my_header[:12]): # For muscle.
                                identity = pmsm.compute_sequence_identity(element.seq,structure.seq)
                                pair_list[index].append(identity)
                    index += 1
        return pair_list


    def align_highest_identity_pairs_list(self,pair_list):
        alignment_list = []
        for seq_counter, compared in enumerate(pair_list):
            pair_to_align=[]
            for num in range(len(self.selected_sequences_in_target_alignment)):
                # For each sequence to add perform an aligment to the sequence to which it has the
                # highest identity according to the initial alignment.
                if compared[num]==max(compared):
                    aligned_pair_name = "temp_seq_" + str(seq_counter)
                    pair_to_align.append(self.selected_sequences_in_target_alignment[num])
                    pair_to_align.append(self.elements_to_add[seq_counter])
                    self.perform_alignment(
                            pair_to_align,
                            alignment_name=aligned_pair_name,
                            use_parameters_from_gui=True)
                    alignment_list.append(aligned_pair_name)
                    break
        return alignment_list


    def alignments_joiner(self, al1, al2, output_file_name="al_result"):
        """
        The algorithm that actually builds the joined alignment.
        The first file is an alignment file in "PyMod" format, the second is an alignment file in
        .aln (clustal) format.
        """

        # Take the sequences of the CE-aligned structures.
        struct=open(al1, "r")
        structs=[]
        for structure in struct.readlines(): # Maybe just take the sequences instead of the whole line of the file.
            structs.append([structure])
        struct.close()

        # Take the sequences of the non CE-aligned elements to be aligned.
        mot_and_sons_1 = open(al2, "rU")
        records1 = list(SeqIO.parse(mot_and_sons_1, "clustal"))
        mot_and_sons_1.close()

        # Finds the sequence that the .txt and .aln alignments have in common, the "bridge" sequence.
        for ax in range(len(structs)):
            for line in range(len(records1)):

                # Finds the bridge.
                if (structs[ax][0].split()[0].replace(':', '_') == records1[line].id or
                    structs[ax][0].split()[0] == records1[line].id or
                    structs[ax][0].split()[0].replace("_",":").replace(":","_",1) == records1[line].id):

                    # Builds a list with as many sub-list elements as the sequences aligned in the
                    # .txt alignment.
                    seq1=[]
                    for s1 in range(len(structs)):
                        seq1.append([])

                    # Builds a list with as many sub-list elements as the sequences aligned in the
                    # .ali alignment.
                    seq2=[]
                    for s2 in range(len(records1)):
                        seq2.append([])

                    index1=0
                    index2=0

                    index1_max=len(structs[ax][0].split()[1])
                    index2_max=len(records1[line].seq)

                    # ---
                    # This is basically an implementation of the part of the "center star" alignment
                    # method that adds new sequences to the center star by padding indels when
                    # needed. Here the "bridge" sequence is the "center star".
                    # ---

                    # This catches the exception thrown when one of the indices of the sequences goes out of range.
                    try:
                        # Start to parse the same bridge sequence in the two alignments.
                        for aa in range(10000):
                            # If the indices are referring to the same residue in the two "versions"
                            # of the bridge sequence.
                            if structs[ax][0].split()[1][index1] == records1[line].seq[index2]:
                                for son in range(len(structs)):
                                    seq1[son].append(structs[son][0].split()[1][index1])
                                for son2 in range(len(records1)):
                                    seq2[son2].append(records1[son2].seq[index2])
                                index1+=1
                                index2+=1

                            # If one of the sequences have an indel.
                            if structs[ax][0].split()[1][index1] == '-' and records1[line].seq[index2] != '-':
                                for son in range(len(structs)):
                                    seq1[son].append(structs[son][0].split()[1][index1])
                                for son2 in range(len(records1)):
                                        seq2[son2].append('-')
                                index1+=1

                            # If one of the sequences have an indel.
                            if structs[ax][0].split()[1][index1] != '-' and records1[line].seq[index2] == '-':
                                for son in range(len(structs)):
                                    seq1[son].append('-')
                                for son2 in range(len(records1)):
                                    seq2[son2].append(records1[son2].seq[index2])
                                index2+=1

                    except:
                        stopped_index1=index1
                        stopped_index2=index2
                        if index1>=index1_max:
                            for son in range(len(structs)):
                                for a in range(index2_max-index2):
                                    seq1[son].append('-')
                            for son2 in range(len(records1)):
                                new_index2=stopped_index2
                                for b in range(index2_max-stopped_index2):
                                    seq2[son2].append(records1[son2].seq[new_index2])
                                    new_index2+=1
                        if index2>=index2_max:
                            for son in range(len(records1)):
                                for a in range(index1_max-index1):
                                    seq2[son].append('-')
                            for son2 in range(len(structs)):
                                new_index1=stopped_index1
                                for b in range(index1_max-stopped_index1):
                                    seq1[son2].append(structs[son2][0].split()[1][new_index1])
                                    new_index1+=1

                    # Write the results to the al_result.txt file.
                    f=open(os.path.join(self.alignments_directory, output_file_name) + ".txt", "w")
                    for seq_file1 in range(0,ax+1):
                        print >>f, structs[seq_file1][0].split()[0], "".join(seq1[seq_file1])
                    for seq_file2 in range(len(records1)):
                        if seq_file2 != line:
                            print >>f, records1[seq_file2].id, "".join(seq2[seq_file2])
                    for seq_file1_again in range(ax+1, len(structs)):
                        print >>f, structs[seq_file1_again][0].split()[0], "".join(seq1[seq_file1_again])

                    f.close()
                    break # This stops the cicle when a bridge sequence has been found.
                else:
                    pass
                    # raise Exception("A bridge sequence was not found in the two aligments...")


    ###################################
    # Methods to perform the          #
    # "alignment joining" mode.       #
    ###################################

    def perform_alignment_joining(self):

        # Prepares alignment files containing the alignments which have to be joined.
        self.alignments_to_join_file_list=[]

        for (i,cluster) in enumerate(self.involved_cluster_elements_list):
            # Build the .fasta files with the alignments.
            file_name = "cluster_" + str(i)
            children = self.get_children(cluster)
            self.build_sequences_file(children, file_name, file_format="clustal", remove_indels=False)
            self.alignments_to_join_file_list.append(file_name)

        # Finds the bridges.
        # If the bridges are specified by the user.
        user_selected_bridges = True
        bridges_list =  []
        if user_selected_bridges:
            children_list = [e for e in self.elements_to_align if e.is_child]
            mothers_list = [e for e in self.elements_to_align if e.is_mother and not e.is_cluster_element()]
            bridges_list = children_list[:] + mothers_list[:]
        # If the bridges are to be found by Pymod. To be implemented.
        else:
            # Perform an initial alignment between all the selected sequences.
            pass

        # Performs an alignment between the "bridges".
        self.elements_to_align = bridges_list
        self.bridges_alignment_name = "bridges_alignment"
        self.perform_alignment(self.elements_to_align, self.bridges_alignment_name, use_parameters_from_gui=True)

        # Builds an al_result.txt file for this alignment.
        self.alignment_joining_output = "al_result"
        self.convert_alignment_format(self.bridges_alignment_name +".aln", self.alignment_joining_output)

        # Actually joins the alignments and produces a final .txt file with the result.
        for alignment_file_name in self.alignments_to_join_file_list:
            self.alignments_joiner(
                os.path.join(self.alignments_directory, self.alignment_joining_output + ".txt"),
                os.path.join(self.alignments_directory, alignment_file_name + ".aln"))

        # The temporary file will be deleted inside the update_aligned_sequences() method later.


    ###################################
    # Methods to perform sequence to  #
    # profile alignments.             #
    ###################################

    def perform_sequence_to_profile_alignment(self):
        """
        Method used to initialize sequence to profile alignments.
        """
        if self.alignment_program in ("clustalw","clustalo"):
            self.clustal_sequence_profile_alignment()
        elif self.alignment_program == "salign-seq":
            self.salign_sequence_profile_alignment()


    def clustal_sequence_profile_alignment(self):
        """
        Align sequences to a target profile by clustalw/clustalo.
        """

        # List of sequences belonging to profile to be kept (target cluster).
        target_cluster_element = self.selected_cluster_elements_list[self.target_cluster_index]
        target_profile_elements = self.get_children(target_cluster_element)

        # List of sequences to be appended to target cluster.
        self.elements_to_add = [e for e in self.get_selected_sequences() if not e in target_profile_elements]

        # create target cluster file
        profile_file_name = "cluster_0"
        profile_file_shortcut=os.path.join(self.alignments_directory, profile_file_name+".fasta")
        self.build_sequences_file(target_profile_elements, profile_file_name,
                                  file_format="fasta", remove_indels=False)

        # create sequence file for sequences to be appended to target cluster
        sequences_to_add_file_name = "cluster_1"
        sequences_to_add_file_shortcut=os.path.join(self.alignments_directory, sequences_to_add_file_name+".fasta")
        self.build_sequences_file(self.elements_to_add, sequences_to_add_file_name,
                                  file_format="fasta", remove_indels=True)

        # Output file name.
        sequence_to_profile_output = "al_result"
        output_file_shortcut = os.path.join(self.alignments_directory, sequence_to_profile_output)

        if self.alignment_program=="clustalw":
            clustalw_path = self.clustalw.get_exe_file_path()
            cline='"'         +clustalw_path+'"'+ \
                ' -PROFILE1="'+profile_file_shortcut+'"'+ \
                ' -PROFILE2="'+sequences_to_add_file_shortcut+'" -SEQUENCES -OUTORDER=INPUT'+ \
                ' -MATRIX='   +self.get_clustalw_matrix_value() + \
                ' -GAPOPEN='  +self.get_gapopen_value() + \
                ' -GAPEXT='   +self.get_gapextension_value() + \
                ' -OUTFILE="' +output_file_shortcut+'.aln"'

        elif self.alignment_program=="clustalo":
            clustalo_path = self.clustalo.get_exe_file_path()
            cline='"'           +clustalo_path+'"'+ \
                ' --profile1="' +profile_file_shortcut+'"'+ \
                ' --outfile="'  +output_file_shortcut+'.aln"'+ \
                ' --outfmt=clustal --force'+ \
                ' ' +self.extraoption_entry.get()
            if len(self.elements_to_add)>1:
                cline+=' --infile="'  +sequences_to_add_file_shortcut+'"'
            else:
                cline+=' --profile2="'+sequences_to_add_file_shortcut+'"'

        self.execute_subprocess(cline)

        # Converts the .aln output file into a .txt file, that will be used to update the sequences
        # loaded in PyMod.
        self.convert_alignment_format(sequence_to_profile_output +".aln", sequence_to_profile_output)


    def salign_sequence_profile_alignment(self):

        # List of sequences of profile to be kept (target cluster)
        target_cluster_element = self.selected_cluster_elements_list[self.target_cluster_index]
        alignment_to_keep_elements = self.get_children(target_cluster_element)

        # Used by generate_highest_identity_pairs_list
        self.selected_sequences_in_target_alignment = alignment_to_keep_elements

        # List of the selected sequences to be appended to target cluster.
        self.elements_to_add = [e for e in self.get_selected_sequences() if not e in alignment_to_keep_elements]

        # ---
        # Perform a first sequence alignment between all selected sequences
        # and sequences in target cluster.
        # ---
        initial_alignment_name = "all_temporary"
        self.elements_to_align = alignment_to_keep_elements + self.elements_to_add

        # Perform sequence alignment even if sequence-structure alignment was requested, because the
        # former is signficantly faster.
        self.salign_malign(self.elements_to_align, initial_alignment_name, use_structural_information=False)

        # ---
        # For each sequence to be appended to the alignment, finds the most
        # similiar sequence in the target cluster according to previous
        # multiple sequence alignment. Compute a similarity list containing
        # as many rows as the number of sequences to be added, and as many
        # columns as the number of sequences in target cluster
        # ---

        highest_identity_pairs_list=self.generate_highest_identity_pairs_list(initial_alignment_name)
        max_identity_list=map(max,highest_identity_pairs_list)
        # sort self.elements_to_add according to max_identity_list
        max_identity_list, self.elements_to_add = zip(*sorted(
            zip(max_identity_list,self.elements_to_add),reverse=True))

        # ---
        # Construct PIR format input file
        # ---
        self.alignments_to_join_file_list=[]
        profiles=[alignment_to_keep_elements]+[[e] for e in self.elements_to_add]

        use_str_info = self.get_salign_seq_str_alignment_var()

        for (i,children) in enumerate(profiles):
            file_name = "cluster_" + str(i)
            self.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information = use_str_info)
            self.alignments_to_join_file_list.append(file_name)

        # ---
        # Sequentially apply profile-profile alignment to each element
        # of elements_to_add
        # ---
        self.alignment_joining_output = "al_result"

        self.salign_profile_profile_alignment(output_file_name=self.alignment_joining_output, use_structural_information=use_str_info)


    ###################################
    # Profile to profile alignments.  #
    ###################################

    def profile_profile_alignment(self):

        # Sequences in selected clusters will all be aligned. Sequences not in selected clusters,
        # will not be aligned.
        for cluster in self.selected_cluster_elements_list:
            self.elements_to_align +=list(self.get_children(cluster))

        # This will be used later in the "create_alignment_element()" method.
        self.selected_elements=[e for e in self.selected_elements if e.is_cluster_element() or e in self.elements_to_align]

        self.alignments_to_join_file_list=[] # two MSA files

        use_str_info = False
        if self.alignment_program == "salign-seq":
            use_str_info = self.get_salign_seq_str_alignment_var()

        for (i,cluster) in enumerate(self.selected_cluster_elements_list):
            file_name = "cluster_" + str(i) # Build FASTA with the MSAs.
            children = self.get_children(cluster)

            # Builds a series of alignment files for each selected cluster.
            if self.alignment_program.startswith("clustal"):
                self.build_sequences_file(children, file_name, file_format="clustal", remove_indels = False)

            elif self.alignment_program.startswith("salign"):
                self.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information=use_str_info)

            self.alignments_to_join_file_list.append(file_name)

        self.alignment_joining_output = "al_result"

        if self.alignment_program=="clustalw":
            self.clustal_profile_profile_alignment(matrix=self.get_clustalw_matrix_value(),
                gapopen=int(self.get_gapopen_value()),
                gapext=float(self.get_gapextension_value()),
                output_file_name=self.alignment_joining_output)

        elif self.alignment_program=="clustalo":
            self.clustal_profile_profile_alignment(
                output_file_name=self.alignment_joining_output,
                extraoption=self.extraoption_entry.get())

        elif self.alignment_program.startswith("salign"):
            self.salign_profile_profile_alignment(
                output_file_name=self.alignment_joining_output,use_structural_information=use_str_info)


    def clustal_profile_profile_alignment(self, matrix="blosum", gapopen=10, gapext=0.2, output_file_name="al_result", extraoption=''):

        output_file_shortcut=os.path.join(self.alignments_directory,
            output_file_name)

        profile1=os.path.join(self.alignments_directory,
            self.alignments_to_join_file_list[0]+".aln")

        for profile2 in self.alignments_to_join_file_list[1:]:
            profile2=os.path.join(self.alignments_directory,
                profile2+".aln")

            if self.alignment_program=="clustalo":
                clustalo_path = self.clustalo.get_exe_file_path()
                cline='"'           +clustalo_path+'"' \
                    ' --profile1="' +profile1+'"'+ \
                    ' --profile2="' +profile2+'"'+ \
                    ' --outfile="'  +output_file_shortcut+'.aln"' \
                    ' --outfmt=clustal --force' \
                    ' ' +extraoption
            elif self.alignment_program=="clustalw":
                clustalw_path = self.clustalw.get_exe_file_path()
                cline='"'          +clustalw_path+'"' \
                    ' -PROFILE1="' +profile1+'"'+ \
                    ' -PROFILE2="' +profile2+'" -OUTORDER=INPUT' \
                    ' -MATRIX='    +matrix+ \
                    ' -GAPOPEN='   +str(gapopen)+ \
                    ' -GAPEXT='    +str(gapext)+ \
                    ' -OUTFILE="'  +output_file_shortcut+'.aln"'

            profile1=output_file_shortcut+'.aln'

            self.execute_subprocess(cline)

        self.convert_alignment_format(output_file_name +".aln", output_file_name)


    def salign_profile_profile_alignment(self,output_file_name="al_result",use_structural_information=False):

        profile1_name = self.alignments_to_join_file_list[0]+".ali"
        profile1_shortcut=os.path.join(self.alignments_directory,profile1_name)

        if self.modeller.run_internally():
            modeller.log.minimal()
            env = modeller.environ()
            env.io.atom_files_directory = ['.', self.structures_directory]
            env.io.hetatm = True
            env.libs.topology.read(file="$(LIB)/top_heav.lib")

            for profile2 in [os.path.join(self.alignments_directory,
                e+".ali") for e in self.alignments_to_join_file_list[1:]]:
                # cat profile2 to profile1 and return number of sequences
                # in the original profile1

                ali_txt1=open(profile1_shortcut,'rU').read()
                ali_txt2=open(profile2,'rU').read()
                align_block=len([e for e in ali_txt1.splitlines() \
                    if e.startswith('>')])
                open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)

                aln = modeller.alignment(env, file=profile1_shortcut, alignment_format="PIR")
                if use_structural_information:
                    env.libs.topology.read(file='$(LIB)/top_heav.lib')
                    aln.salign(rr_file='${LIB}/blosum62.sim.mat',
                        gap_penalties_1d=(-500, 0), output='',
                        align_block=align_block, #max_gap_length=20,
                        align_what='PROFILE', alignment_type="PAIRWISE",
                        comparison_type='PSSM',
                        gap_function=True,#structure-dependent gap penalty
                        feature_weights=(1., 0., 0., 0., 0., 0.),
                        gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),
                        similarity_flag=True,
                        substitution=True,smooth_prof_weight=10.0)
                else:
                    aln.salign(rr_file='${LIB}/blosum62.sim.mat',
                    gap_penalties_1d=(-500, 0), output='',
                    align_block=align_block,   # no. of seqs. in first MSA
                    align_what='PROFILE', alignment_type='PAIRWISE',
                    comparison_type='PSSM',
                    similarity_flag=True, substitution=True,
                    smooth_prof_weight=10.0) # For mixing data with priors

                #write out aligned profiles (MSA)
                aln.write(file=profile1_shortcut, alignment_format="PIR")
        else: # create salign_profile_profile.py for external modeller

            for profile2 in [os.path.join(self.alignments_directory,
                e+".ali") for e in self.alignments_to_join_file_list[1:]]:
                # cat profile2 to profile1 and return number of sequences
                # in the original profile1
                ali_txt1=open(profile1_shortcut,'rU').read()
                ali_txt2=open(profile2,'rU').read()
                align_block=len([e for e in ali_txt1.splitlines() if e.startswith('>')])
                open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)

                config=open("salign_profile_profile.py", "w")
                print >>config, "import modeller"
                print >>config, "modeller.log.verbose()"
                print >>config, "env = modeller.environ()"
                print >>config, "env.io.atom_files_directory = ['.', '"+self.structures_directory+"']"
                print >>config, "env.io.hetatm = True"
                print >>config, "aln = modeller.alignment(env, file='%s', alignment_format='PIR')"%(profile1_shortcut)
                if use_structural_information:
                    print >>config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
                    print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), gap_penalties_2d=(0.35,1.2,0.9,1.2,0.6,8.6,1.2,0.0,0.0), similarity_flag=True, substitution=True,smooth_prof_weight=10.0)"%(align_block)
                else:
                    print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', similarity_flag=True, substitution=True, smooth_prof_weight=10.0) "%(align_block)
                print >>config, "aln.write(file='%s', alignment_format='PIR')"%(profile1_shortcut)
                config.close()

                cline=self.modeller.get_exe_file_path()+" salign_profile_profile.py"
                self.execute_subprocess(cline)

            os.remove("salign_profile_profile.py")

        self.convert_alignment_format(profile1_name, output_file_name)


    #################################################################
    # Methods for launching specific alignment programs.            #
    #################################################################

    ###################################
    # ClustalW.                       #
    ###################################

    def run_clustalw(self, sequences_to_align, alignment_name=None, matrix="blosum", gapopen=10, gapext=0.2):
        """
        This method allows to interact with the local ClustalW.
        """

        if self.clustalw.exe_exists():

            # First build an input .fasta file containing the sequences to be aligned.
            self.build_sequences_file(sequences_to_align, alignment_name)
            # Sets the full paths of input and output files.
            input_file_path = os.path.join(self.alignments_directory, alignment_name + ".fasta")
            output_file_path = os.path.join(self.alignments_directory, alignment_name + ".aln")

            # Run an alignment with all the sequences using ClustalW command line, through Biopython.
            cline = ClustalwCommandline(self.clustalw.get_exe_file_path(),
                    infile=input_file_path, outfile=output_file_path, outorder="INPUT",
                    matrix=matrix, gapopen=gapopen, gapext=gapext)

            self.execute_subprocess(str(cline))

        else:
            self.alignment_program_not_found("clustalw")


    ###################################
    # MUSCLE.                         #
    ###################################

    def run_muscle(self, sequences_to_align, alignment_name=None):
        """
        This method allows to interact with the local MUSCLE.
        """

        if self.muscle.exe_exists():

            self.build_sequences_file(sequences_to_align, alignment_name)

            # infasta - input FASTA for muscle
            # outfasta_tree - output FASTA from muscle, in tree order
            # outaln - output ALN in input order
            infasta=os.path.join(self.alignments_directory, alignment_name + ".fasta")
            outfasta_tree=os.path.join(self.alignments_directory, alignment_name + "_tree.fasta")
            outaln=os.path.join(self.alignments_directory, alignment_name + ".aln")

            cline = MuscleCommandline( self.muscle.get_exe_file_path(),
                input= infasta, out = outfasta_tree)
            self.execute_subprocess(str(cline))

            ''' muscle does not respect sequence input order. e.g.
            (infile)                       (outfile_tree)
            >1tsr.pdb                      >1tsr.pdb
            SSSVPSQKTYQGS                  SSSVPSQKTYQGS---
            >4ibs.pdb      MUSCLE v3.8.31  >1TSR_Chain:A
            VPSQKTYQGSYGF  =============>  SSSVPSQKTYQGS---
            >1TSR_Chain:A                  >4ibs.pdb
            SSSVPSQKTYQGS                  ---VPSQKTYQGSYGF

            The following code rearranges sequence order in `outfile_tree`
            to match that of `infile`. '''
            try:
                record_outtree=[record for record in SeqIO.parse(outfasta_tree,"fasta")]
                record_outtree_id=[record.id for record in record_outtree]
                record_outtree_ungap=[]
                record_outaln=[]
                for i,record_in in enumerate(SeqIO.parse(infasta,"fasta")):
                    record_index=-1
                    try:
                        record_index=record_outtree_id.index(record_in.id)
                    except:
                        if not len(record_outtree_ungap):
                            record_outtree_ungap=[str(record.seq.ungap('-').ungap('.')) for record in record_outtree]
                        try:
                            record_index=record_outtree_ungap.index(str(record_in.seq.ungap('-').ungap('.')))
                        except:
                            pass
                    if record_index<0:
                        print "Warning! Cannot find record.id="+str(record.id)
                        record_index=i
                    record_outaln.append(record_outtree[record_index])
            except Exception,e:
                print "ERROR!"+str(e)
                record_outaln=SeqIO.parse(outfasta_tree,"fasta")

            SeqIO.write(record_outaln, outaln, "clustal")

        else:
            self.alignment_program_not_found("muscle")


    ###################################
    # ClustalO.                       #
    ###################################

    def run_clustalo(self, sequences_to_align, alignment_name=None, extraoption=""):

        if self.clustalo.exe_exists():
            self.build_sequences_file(sequences_to_align, alignment_name)

            input_file_path = os.path.join(self.alignments_directory, alignment_name + ".fasta")
            output_file_path = os.path.join(self.alignments_directory, alignment_name + ".aln")
            guidetree_file_path = os.path.join(self.alignments_directory, alignment_name + ".dnd")

            cline = ClustalOmegaCommandline(
                self.clustalo.get_exe_file_path(),
                infile= input_file_path,
                outfile= output_file_path,
                guidetree_out=guidetree_file_path,
                force=True, outfmt="clustal")

            # Run MSA with all sequences using CLustalO command line.
            cline = str(cline) + ' ' + extraoption
            self.execute_subprocess(cline)

        else:
            self.alignment_program_not_found("clustalo")


    ###################################
    # CE-alignment.                   #
    ###################################

    def run_ce_alignment(self, structures_to_align,ce_alignment_file_name=None, use_seq_info=False):
        """
        Used to launch Ce_align.
        """

        # If there are just two selected sequences, just call self.CE_align().
        if len(structures_to_align) == 2:
            current_elements_to_align = structures_to_align[:]
            # Just produce as output an .aln file that will be used by the
            # "display_ordered_sequences" method.
            self.CE_align(current_elements_to_align,output_format="aln",ce_output_file_name=ce_alignment_file_name, use_seq_info=use_seq_info)

        # Multiple structural alignment: Ce_align two sequences per round.
        else:
            backup_list= structures_to_align[:]

            # Align the first two structures and produces an ce_temp.txt alignment file.
            temp_ce_alignment = "ce_temp"
            current_elements_to_align = backup_list[0:2]
            self.CE_align(current_elements_to_align,ce_output_file_name=temp_ce_alignment,output_format="txt", use_seq_info=use_seq_info)

            # Align the rest of the structures to the first one progressively.
            for n in range(2,len(backup_list)):
                current_elements_to_align = [backup_list[0],backup_list[n]]
                self.CE_align(current_elements_to_align,ce_output_file_name=temp_ce_alignment,output_format="aln", use_seq_info=use_seq_info)
                txt_file_path = os.path.join(self.alignments_directory, temp_ce_alignment + ".txt")
                aln_file_path = os.path.join(self.alignments_directory, temp_ce_alignment + ".aln")
                self.alignments_joiner(txt_file_path, aln_file_path, output_file_name = temp_ce_alignment )

            # Complete by cleaning up the temporary files and by creating a final output file.
            os.remove(os.path.join(self.alignments_directory, temp_ce_alignment + ".aln"))

            # In this cases pymod will need a .txt format alignment file.
            if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
                # Creates the final alignment file. It will be deleted inside
                # update_aligned_sequences() method.
                shutil.copy(os.path.join(self.alignments_directory, temp_ce_alignment + ".txt"),
                            os.path.join(self.alignments_directory, ce_alignment_file_name + ".txt") )

            # In this other cases pymod will need an .aln alignment file, so an .aln file has to be
            # built from a .txt file.
            elif self.alignment_mode in ("alignment-joining", "keep-previous-alignment"):
                # Parses the .txt alignment file.
                fh = open(os.path.join(self.alignments_directory, temp_ce_alignment+".txt"),"r")
                sequences = []
                for line in fh.readlines():
                    sequences.append(line.split())
                fh.close()

                # And then builds an .aln file copy of the alignment.
                records = []
                for s in sequences:
                    rec = SeqRecord(Seq(s[1]), id=s[0])
                    records.append(rec)
                handle = open(os.path.join(self.alignments_directory, ce_alignment_file_name+".aln"), "w")
                SeqIO.write(records, handle, "clustal")
                handle.close()

            os.remove(os.path.join(self.alignments_directory, temp_ce_alignment + ".txt"))


    def CE_align(self, elements_to_align,use_seq_info=False,ce_output_file_name=None,output_format="txt", ce_mode=0):
        """
        Actually performs the structural alignment.
        ce_mode:
            - 0: use sequence information to drive the structural alignment.
            - 1: don't use sequence informtaion.
        ce_output_file_name: the name of the alignment.
        output_format:
            - "txt": it will produce a pymod format alignment file.
            - "aln": it will produce an .ali alignment file in clustal format.
        """
        # Run CE-alignment using the external module.
        if ce_alignment_mode == "plugin":
            ############################################################################
            #
            #  Copyright (c) 2007, Jason Vertrees.
            #  All rights reserved.
            #
            #  Redistribution and use in source and binary forms, with or without
            #  modification, are permitted provided that the following conditions are
            #  met:
            #
            #      * Redistributions of source code must retain the above copyright
            #      notice, this list of conditions and the following disclaimer.
            #
            #      * Redistributions in binary form must reproduce the above copyright
            #      notice, this list of conditions and the following disclaimer in
            #      the documentation and/or other materials provided with the
            #      distribution.
            #
            #  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
            #  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
            #  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
            #  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
            #  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
            #  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
            #  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
            #  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
            #  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
            #  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
            #  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
            #
            #############################################################################

            # Takes the header and the chain id of the selected chains.
            def prepare_data_for_ce_alignment(element,n):
                sel_header = element.my_header
                chain_id = element.my_header.split(':')[-1]
                sel_name = element.structure.chain_pdb_file_name_root
                sel = element.structure.chain_pdb_file_name_root # element.my_header.replace(":", "_")
                return sel_header, chain_id, sel_name, sel

            #########################################################################
            def simpAlign( mat1, mat2, name1, name2, mol1=None, mol2=None, align=0, L=0 ):
                # check for consistency
                assert(len(mat1) == len(mat2))

                # must alway center the two proteins to avoid
                # affine transformations.  Center the two proteins
                # to their selections.
                COM1 = numpy.sum(mat1,axis=0) / float(L)
                COM2 = numpy.sum(mat2,axis=0) / float(L)
                mat1 = mat1 - COM1
                mat2 = mat2 - COM2

                # Initial residual, see Kabsch.
                E0 = numpy.sum( numpy.sum(mat1 * mat1,axis=0),axis=0) + numpy.sum( numpy.sum(mat2 * mat2,axis=0),axis=0)

                #
                # This beautiful step provides the answer.  V and Wt are the orthonormal
                # bases that when multiplied by each other give us the rotation matrix, U.
                # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
                V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(mat2), mat1))

                # we already have our solution, in the results from SVD.
                # we just need to check for reflections and then produce
                # the rotation.  V and Wt are orthonormal, so their det's
                # are +/-1.
                reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
                if reflect == -1.0:
                    S[-1] = -S[-1]
                    V[:,-1] = -V[:,-1]

                RMSD = E0 - (2.0 * sum(S))
                RMSD = numpy.sqrt(abs(RMSD / L))

                if ( align == 0 ):
                    return RMSD;

                assert(mol1 != None)
                assert(mol2 != None)

                #U is simply V*Wt
                U = numpy.dot(V, Wt)

                # rotate and translate the molecule
                mat2 = numpy.dot((mol2 - COM2), U) + COM1
                stored.sel2 = mat2.tolist()

                # let PyMol know about the changes to the coordinates
                cmd.alter_state(1,name2,"(x,y,z)=stored.sel2.pop(0)")

                if False:
                    print "NumAligned=%d" % L
                    print "RMSD=%f" % RMSD

            def cealign( sel1, sel2, verbose=1 ):
                winSize = 8
                # FOR AVERAGING
                winSum = (winSize-1)*(winSize-2) / 2;
                # max gap size
                gapMax = 30

                # make the lists for holding coordinates
                # partial lists
                stored.sel1 = []
                stored.sel2 = []
                # full lists
                stored.mol1 = []
                stored.mol2 = []

                # now put the coordinates into a list
                # partials

                # -- REMOVE ALPHA CARBONS
                sel1 = sel1 + " and n. CA"
                sel2 = sel2 + " and n. CA"
                # -- REMOVE ALPHA CARBONS

                cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
                cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")

                # full molecule
                mol1 = cmd.identify(sel1,1)[0][0]
                mol2 = cmd.identify(sel2,1)[0][0]

                # put all atoms from MOL1 & MOL2 into stored.mol1
                cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
                cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")

                if ( len(stored.mol1) == 0 ):
                        print "ERROR: Your first selection was empty."
                        return
                if ( len(stored.mol2) == 0 ):
                        print "ERROR: Your second selection was empty."
                        return

                # call the C function
                alignString = ccealign( (stored.sel1, stored.sel2) )

                if ( len(alignString) == 1 ):
                    if ( len(alignString[0]) == 0 ):
                        print "\n\nERROR: There was a problem with CEAlign's C Module.  The return value was blank."
                        print "ERROR: This is obviously bad.  Please inform a CEAlign developer.\n\n"
                        return

                bestPathID = -1
                bestPathScore = 100000
                bestStr1 = ""
                bestStr2 = ""

                # for each of the 20 possible alignments returned
                # we check each one for the best CE-Score and keep
                # that one.  The return val of ccealign is a list
                # of lists of pairs.
                for curAlignment in alignString:
                    seqCount = len(curAlignment)
                    matA = None
                    matB = None

                    if ( seqCount == 0 ):
                            continue;

                    for AFP in curAlignment:
                            first, second = AFP
                            if ( matA == None and matB == None ):
                                matA = [ stored.sel1[first-1] ]
                                matB = [ stored.sel2[second-1] ]
                            else:
                                matA.append( stored.sel1[first-1] )
                                matB.append( stored.sel2[second-1] )

                    curScore = simpAlign( matA, matB, mol1, mol2, stored.mol1, stored.mol2, align=0, L=len(matA) )

                    #########################################################################
                    # if you want the best RMSD, not CE Score uncomment here down
                    #########################################################################
                    #if ( curScore < bestPathScore ):
                            #bestPathScore = curScore
                            #bestMatA = matA
                            #bestMatB = matB
                    #########################################################################
                    # if you want the best RMSD, not CE Score uncomment here up
                    #########################################################################

                    #########################################################################
                    # if you want a proven, "better" alignment use the CE-score instead
                    # Uncomment here down for CE-Score
                    #########################################################################
                    internalGaps = 0.0;
                    for g in range(0, seqCount-1):
                        if (not curAlignment[g][0] + 1 == curAlignment[g+1][0]):
                                internalGaps += curAlignment[g+1][0]
                        if ( not curAlignment[g][1] + 1 == curAlignment[g+1][1] ):
                                internalGaps += curAlignment[g+1][1]

                        aliLen = float( len(curAlignment))
                        numGap = internalGaps;
                        curScore = float((curScore/aliLen)*(1.0+(numGap/aliLen)));

                    if ( curScore < bestPathScore ):
                        bestPathScore = curScore
                        bestMatA = matA
                        bestMatB = matB
                    #########################################################################
                    # if you want a proven, "better" alignment use the CE-score instead
                    # Uncomment here UP for CE-Score
                    #########################################################################

                # align the best one string
                simpAlign(bestMatA, bestMatB, mol1, mol2, stored.mol1, stored.mol2, align=1, L=len(bestMatA))

            # Performs an alignment between the two sequences.
            def NWrun(s1,s2,pdb1,pdb2,bio_str1,bio_str2,selection1_header,selection2_header,sequence_information=False):
                sc = pmsp.initScore()
                # set up the box
                L = pmsp.setUp(s1,s2,sc, pdb1, pdb2, bio_str1, bio_str2)
                #aaNames,m = BLOSUM.loadMatrix(fn='blosum50.txt')
                aaNames,m = pmsp.loadMatrix()
                pmsp.doScoring(L,s1,s2,m,sc,sequence_information)
                seq1,seq2 = pmsp.trackback(L,s1,s2,m)

                # Builds an output file with the sequences aligned according to the results of the
                # structural alignment.
                pymod_elements=[PyMod_element(record_seq=seq1, record_header=selection1_header, adjust_header=False),
                                PyMod_element(record_seq=seq2, record_header=selection2_header, adjust_header=False)]
                if output_format == "txt":
                    self.build_sequences_file(elements=pymod_elements, sequences_file_name=ce_output_file_name,
                        file_format="pymod", remove_indels=False)
                elif output_format == "aln":
                    self.build_sequences_file(elements=pymod_elements, sequences_file_name=ce_output_file_name,
                        file_format="clustal", remove_indels=False)
            #########################################################################

            sel1_header, chain_id1, sel1_name, sel1 = prepare_data_for_ce_alignment(elements_to_align[0],1)
            sel2_header, chain_id2, sel2_name, sel2 = prepare_data_for_ce_alignment(elements_to_align[1],2)

            cealign (sel1, sel2)
            cmd.show('cartoon', sel1 + ' or ' + sel2)
            cmd.center('visible')
            cmd.zoom('visible')

            # Updates the names of the chains PDB files.
            saved_file1 = sel1_name + "_aligned.pdb"
            saved_file2 = sel2_name + "_aligned.pdb"

            elements_to_align[0].structure.chain_pdb_file_name = saved_file1
            elements_to_align[1].structure.chain_pdb_file_name = saved_file2

            # And saves these new files.
            aligned_pdb_file1 = os.path.join(self.structures_directory, saved_file1)
            aligned_pdb_file2 = os.path.join(self.structures_directory, saved_file2)
            cmd.save(aligned_pdb_file1, sel1)
            cmd.save(aligned_pdb_file2, sel2)

            # Finally retrieves the structural alignment between the sequences.
            structure1 = Bio.PDB.PDBParser().get_structure(sel1, aligned_pdb_file1)
            structure2 = Bio.PDB.PDBParser().get_structure(sel2, aligned_pdb_file2)

            s1 = Seq(str(elements_to_align[0].my_sequence).replace("-",""))
            s2 = Seq(str(elements_to_align[1].my_sequence).replace("-",""))

            working_dir = os.getcwd()

            # This will generate the alignment output file with the results of the structural
            # alignment.
            NWrun(s1, s2,
                  os.path.join(working_dir, aligned_pdb_file1), os.path.join(working_dir, aligned_pdb_file2),
                  structure1, structure2,
                  sel1_header, sel2_header,
                  sequence_information = use_seq_info)

        # Run CE-alignment using the PyMOL built-in module.
        elif ce_alignment_mode == "pymol":
            sel1 = elements_to_align[0].build_chain_selector_for_pymol()
            sel2 = elements_to_align[1].build_chain_selector_for_pymol()

            cmd.cealign(target=sel1, mobile=sel2, object="pymod_temp_cealign")

            cmd.center('visible')
            cmd.zoom('visible')

            # Updates the names of the chains PDB files.
            saved_file1 = elements_to_align[0].structure.chain_pdb_file_name_root + "_aligned.pdb"
            saved_file2 = elements_to_align[1].structure.chain_pdb_file_name_root + "_aligned.pdb"

            elements_to_align[0].structure.chain_pdb_file_name = saved_file1
            elements_to_align[1].structure.chain_pdb_file_name = saved_file2

            # And saves these new files.
            aligned_pdb_file1 = os.path.join(self.structures_directory, saved_file1)
            aligned_pdb_file2 = os.path.join(self.structures_directory, saved_file2)
            self.pymol_save(aligned_pdb_file1, sel1)
            self.pymol_save(aligned_pdb_file2, sel2)

            # Finally saves the structural alignment between the sequences.
            cmd.save(os.path.join(self.alignments_directory, ce_output_file_name+".aln"),"pymod_temp_cealign")
            cmd.delete("pymod_temp_cealign")

            # Converts it in .txt format.
            if output_format == "txt":
                self.convert_alignment_format(ce_output_file_name+".aln", ce_output_file_name)


    ###################################
    # Modeller-based alignment        #
    # algorithms.                     #
    ###################################

    def salign_malign(self, sequences_to_align, alignment_name=None, use_structural_information=False):
        """
        alignment.malign - align sequences
        alignment.align2d - sequence-structure alignment
        """

        if self.modeller.can_be_launched():

            shortcut_to_temp_files= os.path.join(self.alignments_directory,alignment_name)
            # The .pir file will be written in a different way if the user decides to use
            # structural information in the alignment.
            self.build_sequences_file(self.elements_to_align, alignment_name, file_format="pir", use_structural_information=use_structural_information)

            if self.modeller.run_internally():
                modeller.log.minimal()
                env = modeller.environ()
                env.io.atom_files_directory = ['.', self.structures_directory]
                env.io.hetatm = True
                aln = modeller.alignment(env,
                                         file=shortcut_to_temp_files +".ali",
                                         alignment_format='PIR')
                if use_structural_information:
                    env.libs.topology.read(file="$(LIB)/top_heav.lib")
                    # Structure sensitive variable gap penalty alignment:
                    aln.salign(auto_overhang=True,
                        gap_penalties_1d=(-100, 0),
                        gap_penalties_2d=(3.5,3.5,3.5,.2,4.,6.5,2.,0.,0.),
                        gap_function=True, # structure-dependent gap penalty
                        feature_weights=(1., 0., 0., 0., 0., 0.),
                        similarity_flag=True,
                        alignment_type='tree', #output='ALIGNMENT',
                        dendrogram_file=shortcut_to_temp_files+".tree")
                else:
                    aln.salign(auto_overhang=True, gap_penalties_1d=(-450, 0),
                       alignment_type='tree', output='ALIGNMENT',
                       dendrogram_file=shortcut_to_temp_files+".tree")
                aln.write(file=shortcut_to_temp_files +'.ali', alignment_format='PIR')

            else:
                # create salign_multiple_seq.py to enable external modeller execution
                config=open("salign_multiple_seq.py", "w")
                print >> config, "import modeller"
                print >> config, "modeller.log.verbose()"
                print >> config, "env = modeller.environ()"
                print >> config, "env.io.atom_files_directory = ['.', '"+self.structures_directory+"']"
                print >> config, "env.io.hetatm = True"
                print >> config, "aln = modeller.alignment(env,file='%s', alignment_format='PIR')" % (shortcut_to_temp_files + ".ali")
                if use_structural_information:
                    print >> config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
                    print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-100, 0), gap_penalties_2d=(3.5,3.5,3.5,0.2,4.0,6.5,2.0,0.0,0.0), gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), similarity_flag=True, alignment_type='tree', dendrogram_file='%s')" %(shortcut_to_temp_files+".tree")
                else:
                    print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50), dendrogram_file='%s', alignment_type='tree', output='ALIGNMENT')" %(shortcut_to_temp_files+".tree")
                print >> config, "aln.write(file='"+shortcut_to_temp_files+".ali', alignment_format='PIR')"
                print >> config, ""
                config.close()

                cline=self.modeller.get_exe_file_path()+" salign_multiple_seq.py"
                self.execute_subprocess(cline)
                os.remove("salign_multiple_seq.py") # remove this temporary file.

            # convert alignment_name.ali to alignment_tmp.fasta
            record=SeqIO.parse(open(shortcut_to_temp_files + ".ali"),"pir")
            SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")

        else:
            self.alignment_program_not_found("salign-seq")


    def salign_align3d(self, structures_to_align, alignment_name=None):
        """
        alignment.malign3d - align structures
        """

        if self.modeller.can_be_launched():

            # if sys.platform=="win32":
            #     sys.path.append(self.modeller.get_exe_file_path()+"\modlib")
            if len(structures_to_align)>2:
                self.build_salign_dendrogram_menu=True
            else: # salign only output dendrogram_file when there are 3 sequences or more
                self.build_salign_dendrogram_menu=False

            shortcut_to_temp_files = os.path.join(self.current_project_directory_full_path,self.alignments_directory,alignment_name)
            struct_tup=range(0,len(structures_to_align))
            for ii in range(0,len(structures_to_align)):
                struct_entry=structures_to_align[ii].structure.chain_pdb_file_name_root
                header = structures_to_align[ii].my_header
                chain_id=structures_to_align[ii].structure.pdb_chain_id
                struct_tup[ii]=(struct_entry,header,chain_id)

            # Change the working directory, so that the ouptut files will be created in the structures
            # directory.
            os.chdir(self.structures_directory)

            if self.modeller.run_internally():
                modeller.log.minimal()
                env = modeller.environ()
                aln = modeller.alignment(env)

                for (pdb_file_name, code, chain) in struct_tup:
                    mdl = modeller.model(env, file=pdb_file_name,
                                     model_segment=("FIRST:"+chain,"LAST:"+chain))
                    aln.append_model(mdl, atom_files=pdb_file_name, align_codes=code)

                for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1., 0.5, 1., 1., 1., 0.), False, True),
                                        ((1., 1., 1., 1., 1., 0.), True, False)):
                    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                           rr_file="$(LIB)/as1.sim.mat", overhang=30,
                           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                           gap_gap_score=0, gap_residue_score=0,
                           dendrogram_file= shortcut_to_temp_files + ".tree",
                           alignment_type="tree", feature_weights=weights,
                           improve_alignment=True, fit=True, write_fit=write_fit,
                           write_whole_pdb=whole,output="ALIGNMENT QUALITY")

                aln.write(file=shortcut_to_temp_files +".ali", alignment_format="PIR")

                aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
                       rr_file='$(LIB)/as1.sim.mat', overhang=30,
                       gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                       gap_gap_score=0, gap_residue_score=0,
                       dendrogram_file=shortcut_to_temp_files + '.tree',
                       alignment_type='progressive', feature_weights=[0]*6,
                       improve_alignment=False, fit=False, write_fit=True,
                       write_whole_pdb=False,output='QUALITY')

            else: # except:
                # create salign_multiple_struc.py for external modeller execution

                config=open("salign_multiple_struc.py", "w")
                print >> config, "import modeller"
                print >> config, "modeller.log.verbose()"
                print >> config, "env = modeller.environ()"
                print >> config, "aln = modeller.alignment(env)"
                for (pdb_file_name, code, chain) in struct_tup:
                    print >> config, "mdl = modeller.model(env, file='"+pdb_file_name+"', model_segment=('FIRST:"+chain+"','LAST:"+chain+"'))"
                    print >> config, "aln.append_model(mdl, atom_files='"+pdb_file_name+"', align_codes='"+code+"')"
                print >> config, "for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True), ((1., 0.5, 1., 1., 1., 0.), False, True), ((1., 1., 1., 1., 1., 0.), True, False)):"
                print >> config, "    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='tree', feature_weights=weights, improve_alignment=True, fit=True, write_fit=write_fit, write_whole_pdb=whole, output='ALIGNMENT QUALITY')" % (shortcut_to_temp_files)
                print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
                print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
                print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
                print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
                print >> config, ""
                config.close()

                cline=self.modeller.get_exe_file_path()+" salign_multiple_struc.py"
                self.execute_subprocess(cline)
                os.remove("salign_multiple_struc.py") # Remove this temp file.

            # Returns back to the project dir from the project/Structures directory.
            os.chdir(self.current_project_directory_full_path)

            # SALIGN does not superpose ligands. The generated "*_fit.pdb"
            # files are therefore ligandless. The following loop superposes
            # original structure to saligned structures, and replaces
            # "*_fit.pdb" files with the superposed liganded original structure.
            for (pdb_file_name_root, code, chain) in struct_tup:
                fixed= os.path.join(self.structures_directory,pdb_file_name_root + "_fit.pdb")
                cmd.load(fixed,"salign_fixed_fit")
                if hasattr(cmd,"super"): # super is sequence-independent
                    cmd.super(pdb_file_name_root,"salign_fixed_fit")
                else: # PyMOL 0.99 does not have cmd.super
                    cmd.align(pdb_file_name_root,"salign_fixed_fit")
                self.pymol_save(fixed,pdb_file_name_root) # quick-and-dirty
                cmd.delete("salign_fixed_fit")

            # Updates the name of the chains PDB files.
            for element in structures_to_align:
                element.structure.chain_pdb_file_name = element.structure.chain_pdb_file_name_root+"_fit.pdb"

            # Convert the PIR format output file into a clustal format file.
            record=SeqIO.parse(open(shortcut_to_temp_files + '.ali',"rU"),"pir")
            SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")

        else:
            self.alignment_program_not_found("salign-str")


    #################################################################
    # File managment during the alignment process.                  #
    #################################################################

    def set_current_alignment_file_name(self, alignment_file_name):
        """
        If the "alignment_file_name" argument is set to "None" (this happens when performing a new
        alignment from the PyMod main menu), the "set_current_alignment_file_name" will
        automatically generate a name for it, using the standard "self.alignments_files_names"
        value.
        """
        if alignment_file_name == None:
            # Alignment files ending with the unique_id of the alignment are going to be created.
            alignment_file_name = "temp_" + self.alignments_files_names + str(self.unique_index)
        else:
            pass

        return alignment_file_name


    def build_sequences_file(self, elements, sequences_file_name, new_directory=None, file_format="fasta", remove_indels=True, use_structural_information=False, same_length=True, first_element=None):
        """
        Builds a sequence file (the format is specified in the alignment_"format" argument) that will
        contain the sequences supplied in the "elements" which has to contain a list of
        "PyMod_element" class objects.
        """

        alignment_extension = pmdt.alignment_extensions_dictionary[file_format]

        # Full path of the file that will contain the alignment.
        a_fh = None
        if new_directory == None:
            alignment_file_path = os.path.join(self.alignments_directory, sequences_file_name + "." + alignment_extension)
        else:
            alignment_file_path = os.path.join(new_directory, sequences_file_name + "." + alignment_extension)
        a_fh = open(alignment_file_path, 'w')

        if same_length:
            self.adjust_aligned_elements_length(elements)

        if first_element != None:
            elements.remove(first_element)
            elements.insert(0, first_element)

        if file_format == "fasta":
            for element in elements:
                self.prepare_fasta_file(a_fh, element.my_header, element.my_sequence, remove_indels)

        elif file_format == "pir":
            for child in elements:
                header=child.my_header.replace(':','_')
                sequence = str(child.my_sequence)
                if remove_indels:
                    sequence = sequence.replace("-","")
                sequence += '*'
                structure=''
                if hasattr(child.structure,"chain_pdb_file_name_root") and use_structural_information:
                    structure=child.structure.chain_pdb_file_name_root
                    chain=child.structure.pdb_chain_id
                if not structure: # sequence
                    print >> a_fh, ">P1;"+ header
                    print >> a_fh, "sequence:"+header+":::::::0.00:0.00"
                else: # structure
                    print >> a_fh, ">P1;"+header+chain
                    print >> a_fh, "structure:"+structure+":.:"+chain+":.:"+chain+":::-1.00:-1.00"
                for ii in xrange(0,len(sequence),75):
                    print >> a_fh, sequence[ii:ii+75].replace("X",".")

        elif file_format == "clustal":
            if remove_indels:
                records = [SeqRecord(Seq(str(child.my_sequence).replace("-","")), id=child.my_header) for child in elements]
            else:
                records = [SeqRecord(Seq(str(child.my_sequence)), id=child.my_header) for child in elements]
            SeqIO.write(records, a_fh, "clustal")

        elif file_format == "pymod":
            for element in elements:
                self.prepare_pymod_alignment_file(a_fh, element.my_header, element.my_sequence, remove_indels)

        else:
            self.show_error_message("File Error", "Unrecognized sequence file format '%s'." % (file_format))

        a_fh.close()


    def prepare_fasta_file(self,output_file_handler, header, sequence, remove_indels=True):
        """
        This function allows to write a file in FASTA format containing the sequences provided
        in a for cicle in the "sequence" argument.
        """
        sequence = str(sequence)
        if remove_indels:
            sequence = sequence.replace("-","")
        #| Write an output in Fasta format to the output_file_handler given as argument
        print >> output_file_handler , ">"+header
        for i in xrange(0, len(sequence), 60):
            print >> output_file_handler , sequence[i:i+60]
        print >> output_file_handler , ""


    def prepare_pymod_alignment_file(self, output_file_handler, header, sequence, remove_indels=True):
        sequence = str(sequence)
        if remove_indels:
            sequence = sequence.replace("-","")
        print >> output_file_handler, header, sequence


    def convert_alignment_format(self, input_file_name, output_file_name):
        """
        Converts an alignment file specified in the "input_file_name" argument in an alignment file
        in the pymod (.txt) format.
        """
        input_handle = open(os.path.join(self.alignments_directory, input_file_name), "rU")
        output_handle = open(os.path.join(self.alignments_directory, output_file_name + ".txt"), "w")

        # Gets the name of the format from the input file extension.
        original_extension = os.path.splitext(input_file_name)[1].replace(".","")
        extension_id = pmdt.alignment_extensions_dictionary.values().index(original_extension)
        format_name = pmdt.alignment_extensions_dictionary.keys()[extension_id]

        alignment = AlignIO.read(input_handle, format_name)
        for element in alignment:
            print >> output_handle, element.id, element.seq
        input_handle.close()
        output_handle.close()


    ###############################################################################################
    # CONSERVATION ANALYSIS TOOLS.                                                                #
    ###############################################################################################

    #################################################################
    # Methods to compute conservation values in an alignment        #
    # column.                                                       #
    #################################################################

    def get_conservation_symbol(self,column):
        """
        Calculate the conservation symbol for an alignment based on the positively scoring groups that
        occur in the Gonnet Pam250 matrix (see: http://www.clustal.org/download/clustalx_help.html).
        Takes as an input an alignment column represented by a list.

        "*" indicates positions which have a single, fully conserved residue.

        ":" indicates that one of the following 'strong' groups is fully conserved:
        STA, NEQK, NHQK, NDEQ, QHRK, MILV, MILF, HY, FYW

        "." indicates that one of the following 'weaker' groups is fully conserved:
        CSA, ATV, SAG, STNK, STPA, SGND, SNDEQK, NDEQHK, NEQHRK, FVLIM, HFY
        """

        symbol = "-"
        # If there is a gap in that position of the alignment, then the symbol is automatically a
        # "-".
        if "-" in column:
            symbol = "-"
        else:
            # If it is really slow try to use frozenset.
            residues = set(column)
            # All residues in the column are identical: full conservation.
            if len(residues) == 1:
                symbol = "*"
            else:
                # Strong conservation.
                if residues.issubset("STA"):
                    symbol = ":"
                elif residues.issubset("NEQK"):
                    symbol = ":"
                elif residues.issubset("NHQK"):
                    symbol = ":"
                elif residues.issubset("NDEQ"):
                    symbol = ":"
                elif residues.issubset("QHRK"):
                    symbol = ":"
                elif residues.issubset("MILV"):
                    symbol = ":"
                elif residues.issubset("MILF"):
                    symbol = ":"
                elif residues.issubset("HY"):
                    symbol = ":"
                elif residues.issubset("FYW"):
                    symbol = ":"
                # Weak conservation.
                elif residues.issubset("CSA"):
                    symbol = "."
                elif residues.issubset("ATV"):
                    symbol = "."
                elif residues.issubset("SAG"):
                    symbol = "."
                elif residues.issubset("STNK"):
                    symbol = "."
                elif residues.issubset("STPA"):
                    symbol = "."
                elif residues.issubset("SGND"):
                    symbol = "."
                elif residues.issubset("SNDEQK"):
                    symbol = "."
                elif residues.issubset("NDEQHK"):
                    symbol = "."
                elif residues.issubset("NEQHRK"):
                    symbol = "."
                elif residues.issubset("FVLIM"):
                    symbol = "."
                elif residues.issubset("HFY"):
                    symbol = "."
        return symbol


    def compute_stars(self,elements):
        """
        Computes the "stars" of an alignment.
        """
        sequences = [str(e.my_sequence) for e in elements]
        minimum_length = min([len(s) for s in sequences])

        stars = "" # Computed by Pymod.

        for i in range(0,minimum_length):
            column = []
            symbol = None
            for s in sequences:
                column.append(s[i])
            stars += self.get_conservation_symbol(column)

        return stars


    ###############################################################################################
    # STRUCTURAL ANALYSIS TOOLS.                                                                  #
    ###############################################################################################

    #################################################################
    # Compute the DOPE (Discrete optimized protein energy) of a     #
    # polypeptidic chain using MODELLER.                            #
    #################################################################

    def dope_from_main_menu(self):
        """
        Called when users decide calculate DOPE of a structure loaded in PyMod.
        """
        # Checks if the DOPE profiles can be computed.
        selection = self.get_selected_sequences()
        if not self.modeller.can_be_launched():
            self.show_error_message("MODELLER Error", "MODELLER is missing. In order to compute DOPE scores of a structure, MODELLER has to be installed.")
            return False
        if len(selection) == 0:
            self.show_error_message("Selection Error", "Please select at least one structure to assess.")
            return False
        if not self.all_sequences_have_structure():
            self.show_error_message("Selection Error", "Please select only elements that have a 3D structure currently loaded in PyMOL.")
            return False
        if len(set([seq.mother_index for seq in selection])) != 1:
            self.show_error_message("Selection Error", "You can assess multiple structures DOPE only if they are aligned in the same cluster.")
            return False

        # Ask users if they would like to color the sequences according to their DOPE values.
        title = "Color Option"
        message = "Would you like to color the selected sequences by their DOPE values, once they have been calculated?"
        color_by_dope_choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)

        # Initializes MODELLER.
        if self.modeller.run_internally():
            env = modeller.environ()
            env.io.atom_files_directory = []
            env.io.atom_files_directory.append(".")
            env.io.hetatm = True
            env.io.water = True
            env.libs.topology.read(file='$(LIB)/top_heav.lib')
            env.libs.parameters.read(file='$(LIB)/par.lib')
        else:
            env = None

        # Actually computes the DOPE scores of the polypeptide chains in the user selection.
        for element in selection:
            self.compute_dope(element,env=env)

        # Assigns to each residue of the selected chains a correspoding color according to its DOPE.
        self.assign_dope_items(selection)

        # Color the elements.
        if color_by_dope_choice:
            for element in selection:
                element.color_element_by_dope()
        self.gridder()

        # Shows the DOPE profiles plot.
        dope_graph_mode = None
        if len(selection) == 1:
            dope_graph_mode = "single"
        elif len(selection) >= 2:
            dope_graph_mode = "multiple"

        # Prepares the data to show in the plot.
        dope_plot_data = self.prepare_dope_plot_data(selection, mode = dope_graph_mode)

        # Shows the plot.
        self.show_dope_plot(dope_plot_data)


    def compute_dope(self, element, env=None):
        # Prepares the input for MODELLER.
        e_file_name = element.structure.chain_pdb_file_name_root
        e_file_shortcut = os.path.join(pymod.structures_directory, e_file_name)
        e_profile_file_shortcut = os.path.join(pymod.structures_directory, e_file_name+".profile")
        # Computes the DOPE of the 3D structure of the chain of the 'element'.
        self.compute_dope_of_structure_file(e_file_shortcut, e_profile_file_shortcut, env=env)
        # Reads the output file produced by MODELLER with the DOPE score of the chain of the
        # 'element'.
        dope_scores = self.get_dope_profile(e_profile_file_shortcut)
        element.set_dope_scores(dope_scores)


    def compute_dope_of_structure_file(self, str_file_path, profile_file_path, env=None):
        """
        Uses MODELLER to compute the DOPE of a polypeptidic chain, and ouptuts the results in
        'profile_file_path'. When 'env' is set to 'None', MODELLER will be initialized. If
        MODELLER has already been initialized, the its 'env' varibale can be passed in this
        argument so that it is not initialized again.
        """
        if self.modeller.run_internally():
            if env == None:
                env = modeller.environ()
                env.io.atom_files_directory = []
                env.io.atom_files_directory.append(".")
                env.io.hetatm = True
                env.io.water = True
                env.libs.topology.read(file='$(LIB)/top_heav.lib')
                env.libs.parameters.read(file='$(LIB)/par.lib')
            modstr = complete_pdb(env, str_file_path)
            # Assess with DOPE.
            s = modeller.selection(modstr).only_std_residues() # only_het_residues, only_std_residues, only_water_residues
            # Gets the DOPE score.
            score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',file=profile_file_path, normalize_profile=True, smoothing_window=15)
        else:
            # Builds the MODELLER script file to be executed externally.
            dope_profile_script_file_name = "dope_profile-script.py"
            dope_profile_temp_out_name = "dope_profile_temp_out.txt"
            dope_profile_script_file = open(dope_profile_script_file_name,"w")
            print >> dope_profile_script_file, "import modeller"
            print >> dope_profile_script_file, "from modeller.scripts import complete_pdb"
            print >> dope_profile_script_file, "env = modeller.environ()"
            print >> dope_profile_script_file, "env.io.atom_files_directory = []"
            print >> dope_profile_script_file, "env.io.atom_files_directory.append('.')"
            print >> dope_profile_script_file, "env.io.hetatm = True"
            print >> dope_profile_script_file, "env.io.water = True"
            print >> dope_profile_script_file, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
            print >> dope_profile_script_file, "env.libs.parameters.read(file='$(LIB)/par.lib')"
            print >> dope_profile_script_file, "modstr = complete_pdb(env, '%s')" % str_file_path
            print >> dope_profile_script_file, "s = modeller.selection(modstr).only_std_residues()"
            print >> dope_profile_script_file, "score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',file='%s', normalize_profile=True, smoothing_window=15)" % profile_file_path
            print >> dope_profile_script_file, "\n# Needed to compute DOPE in PyMod when MODELLER is run externally from PyMOL."
            print >> dope_profile_script_file, "dope_profile_out_file = open('%s','w')" % dope_profile_temp_out_name
            print >> dope_profile_script_file, "dope_profile_out_file.write(str(score))"
            print >> dope_profile_script_file, "dope_profile_out_file.close()"
            dope_profile_script_file.close()
            # Executes the script.
            cline = self.modeller.get_exe_file_path() + " " + dope_profile_script_file_name
            self.execute_subprocess(cline)
            # Gets the score from the output generated by the script and cleans up temporary files.
            dope_profile_out_file = open(dope_profile_temp_out_name, "r")
            score = float(eval(dope_profile_out_file.readline()))
            dope_profile_out_file.close()
            os.remove(dope_profile_temp_out_name)
            os.remove(dope_profile_script_file_name)

        return score


    def get_dope_profile(self, profile_file_name, seq=None):
        """
        Read 'profile_file' into a Python list, and add gaps corresponding to the alignment
        sequence 'seq'.
        """
        profile_file = open(profile_file_name,"r")
        vals = []
        # Adds None values to gaps the Modeller way.
        for line in profile_file.readlines():
            res_three_letter_code = line[8:11]
            # Read all non-comment and non-blank lines from the file:
            if not line.startswith('#') and len(line) > 10:
                # Initially do not exclude also water molecules (named 'TIP3') and heteroresidues
                # from the graph.
                spl = line.split()
                vals.append(float(spl[-1]))

        profile_file.close()
        # Add a 'None' value at position '0', so that we effectively count from 1.
        # vals.insert(0, None)
        return vals


    def assign_dope_items(self, selection):
        # Builds a list of all DOPE values of the residues in the selection.
        ldope = []
        for chain_element in selection:
            ldope.extend(chain_element.dope_scores)
        # Takes the min and max values among all the selected residues.
        min_value = min(ldope)
        max_value = max(ldope)
        # An array with the equally sapced limits generated with the list above.
        bins = numpy.array(numpy.linspace(min_value, max_value, num=10))
        for chain_element in selection:
            # An array with all the DOPE values of a single chain in the selection.
            adope = numpy.array(chain_element.dope_scores)
            # An array with the id of the bins where those values reside.
            inds = numpy.digitize(adope, bins)
            # Returns a list like:
            # [(-0.052, 4), (-0.03, 3), (-0.04, 5), (-0.04, 6), (-0.041, 7), (-0.042, 8), (-0.043, 10), ...]
            # which contains for all standard residues of a polypeptidic chain a tuple. The
            # first value of the tuple is the DOPE score of that residues, the second is the id
            # (going from 1 to 10) of the bin where that value resides.
            chain_element.dope_items = []
            for dope_score, bin_id in zip(adope, inds):# zip(ldope, inds):
                chain_element.dope_items.append({"dope-score":dope_score, "interval": bin_id})


    def prepare_dope_plot_data(self, selection, start_from=0, mode="single"):
        """
        Takes a selection of 'PyMod_elemet' objects, takes their DOPE scores and returns the data in
        a dictionary which can be supplied as an argument to the 'show_dope_plot()' in order to
        display it in a plot.
        """
        dope_plot_data = []
        for element in selection:
            # Prepares a list with the PyMOL additional data for each residue of the chain, so that
            # when clicking on some point, the corresponding residue will be highlighted in PyMOL,
            # and the message bar of the plot will be updated.
            residues_names = [res.three_letter_code for res in element.structure.get_all_residues_list()]
            residues_pdb_positions = [res.pdb_position for res in element.structure.get_all_residues_list()]
            pymol_selectors = [element.build_residue_selector_for_pymol(res.pdb_position) for res in element.structure.get_all_residues_list()]
            residue_additional_data = []
            for r_name, r_position, r_selector in zip(residues_names, residues_pdb_positions, pymol_selectors):
                residue_additional_data.append({"residue_name": r_name,
                                     "residue_pdb_position": r_position,
                                     "pymol_selector": r_selector,
                                     "export_label": "%s %s"%(r_name, r_position)})
            element_dope_scores = element.dope_scores[:]
            # If the sequences in the selection are aligned, adjust the profiles by inserting 'None'
            # values for gap positions.
            if mode == "multiple":
                # Insert gaps into the profile corresponding to those in seq:
                # R: seq = str(seq).replace("X","-")
                ri = 0
                seq = element.my_sequence
                for i, res in enumerate(seq):
                    if res != "-":
                        # If the first residue is preceeded by some indels.
                        first_residue_with_preceeding_gaps = False
                        if ri == 0 and i != 0:
                            n_gaps = i
                            for gap in range(n_gaps):
                                element_dope_scores.insert(ri, None)
                                residue_additional_data.insert(ri, {"export_label": "Gap"})
                            ri += 1 + n_gaps
                            first_residue_with_preceeding_gaps = True
                        # Applies to the rest of residues in the sequence.
                        if not first_residue_with_preceeding_gaps:
                            n_gaps = pmsm.get_leading_gaps(seq, i)
                            for gap in range(n_gaps):
                                element_dope_scores.insert(ri+1, None)
                                residue_additional_data.insert(ri+1, {"export_label": "Gap"})
                            ri += 1 + n_gaps

            # For DOPE plots of multiple chains models and templates.
            for g in range(start_from):
                element_dope_scores.insert(0, None)
                residue_additional_data.insert(0, {None:None})

            # Prepares the data.
            dope_plot_data.append({"dope_scores": element_dope_scores,
                                   "additional_data": residue_additional_data,
                                   "label": element.get_compact_header()})# element.my_header[0:15]})

        return dope_plot_data


    def show_dope_plot(self, selection_dope_plot_data):
        x_label_text = None
        message_bar_text_on_update = None
        if len(selection_dope_plot_data) > 1:
            x_label_text = "Alignment position"
            message_bar_text_on_update = "Selected: %s %s of __plot_name__ (alignment position: __x__), DOPE value: __y__"
        else:
            x_label_text = "Residue position"
            message_bar_text_on_update = "Selected: %s %s of __plot_name__, DOPE value: __y__"
        cp = pplt.Custom_plot_window(self.main_window, title="DOPE Profile")
        cp.build_plotting_area(message_bar_initial_text = "Click on the plot to highlight corresponding residues in PyMOL.",
                               update_message_bar=True,
                               message_bar_text_on_update=message_bar_text_on_update,
                               message_bar_vars_on_update=("residue_name","residue_pdb_position"),
                               on_click_action=self.highlight_in_pymol_from_dope_plot,
                               x_label_text=x_label_text,
                               y_label_text="DOPE score")
        for chain_dope_data in selection_dope_plot_data:
            # Adds the new plot corresponding to the current chain.
            cp.add_plot(range(1, len(chain_dope_data["dope_scores"])+1), chain_dope_data["dope_scores"],
                        label=chain_dope_data["label"],
                        additional_data=chain_dope_data["additional_data"])
        cp.show()


    def highlight_in_pymol_from_dope_plot(self, point, plot):
        cmd.select("pymod_selection", point.additional_data["pymol_selector"])
        cmd.center("pymod_selection")


    #################################################################
    # Energy minimization using MODELLER.                           #
    #################################################################

    def energy_minimization(self, model_file_path, parameters_dict, env=None, use_hetatm=True, use_water=True, check_structure=True):
        model_file_directory = os.path.dirname(model_file_path)
        model_file_name = os.path.basename(model_file_path)
        opt_code = model_file_name[:-4]+"_optB"
        #----------------------------------------------------
        if self.modeller.run_internally():
            if env == None:
                env = modeller.environ()
                if use_hetatm:
                    env.io.hetatm = True
                    if use_water:
                        env.io.water = True
                env.libs.topology.read(file='$(LIB)/top_heav.lib')
                env.libs.parameters.read(file='$(LIB)/par.lib')

            # This will optimize stereochemistry of a given model, including non-bonded contacts.
            old_dynamic_coulomb = env.edat.dynamic_coulomb
            env.edat.dynamic_coulomb = True
            old_dynamic_lennard = env.edat.dynamic_lennard
            env.edat.dynamic_lennard = True
            old_contact_shell = env.edat.contact_shell
            env.edat.contact_shell = parameters_dict["non_bonded_cutoff"]
            mdl = complete_pdb(env, model_file_path)
            mdl.write(file=os.path.join(model_file_directory, opt_code+'.ini'))
            # Select all atoms:
            atmsel = modeller.selection(mdl)
            # Generate the restraints:
            if parameters_dict["restraints"]["bond"]:
                mdl.restraints.make(atmsel, restraint_type='bond', spline_on_site=False)
            if parameters_dict["restraints"]["angle"]:
                mdl.restraints.make(atmsel, restraint_type='angle', spline_on_site=False)
            if parameters_dict["restraints"]["dihedral"]:
                mdl.restraints.make(atmsel, restraint_type='dihedral', spline_on_site=False)
            if parameters_dict["restraints"]["improper"]:
                mdl.restraints.make(atmsel, restraint_type='improper', spline_on_site=False)
            if parameters_dict["restraints"]["coulomb"]:
                mdl.restraints.make(atmsel, restraint_type='coulomb', spline_on_site=False)
            if parameters_dict["restraints"]["lj"]:
                mdl.restraints.make(atmsel, restraint_type='lj', spline_on_site=False)
            mdl.restraints.write(file=os.path.join(model_file_directory, opt_code+'.rsr'))
            mpdf = atmsel.energy()

            class SteepestDescent(modeller.optimizers.state_optimizer):
                """
                Very simple steepest descent optimizer, in Python, as reported at:
                http://www.salilab.org/modeller/9v4/manual/node252.html
                """
                # Add options for our optimizer
                _ok_keys = modeller.optimizers.state_optimizer._ok_keys + ('min_atom_shift', 'min_e_diff', 'step_size', 'max_iterations')
                def __init__(self, step_size=0.0001, min_atom_shift=0.01, min_e_diff=1.0, max_iterations=None, **vars):
                    modeller.optimizers.state_optimizer.__init__(self, step_size=step_size,
                                             min_atom_shift=min_atom_shift,
                                             min_e_diff=min_e_diff,
                                             max_iterations=max_iterations, **vars)

                def optimize(self, atmsel, **vars):
                    # Do normal optimization startup
                    modeller.optimizers.state_optimizer.optimize(self, atmsel, **vars)
                    # Get all parameters
                    alpha = self.get_parameter('step_size')
                    minshift = self.get_parameter('min_atom_shift')
                    min_ediff = self.get_parameter('min_e_diff')
                    maxit = self.get_parameter('max_iterations')
                    # Main optimization loop
                    state = self.get_state()
                    (olde, dstate) = self.energy(state)
                    while True:
                        for i in range(len(state)):
                            state[i] -= alpha * dstate[i]
                        (newe, dstate) = self.energy(state)
                        if abs(newe - olde) < min_ediff:
                            print "Finished at step %d due to energy criterion" % self.step
                            break
                        elif self.shiftmax < minshift:
                            print "Finished at step %d due to shift criterion" % self.step
                            break
                        elif maxit is not None and self.step >= maxit:
                            print "Finished at step %d due to step criterion" % self.step
                            break
                        if newe < olde:
                            alpha *= 2
                        else:
                            alpha /= 2
                        olde = newe
                        self.next_step()
                    self.finish()

            # Open a file to get basic stats on each optimization.
            trcfil = file(os.path.join(model_file_directory, opt_code+'.D00000001'),'w')
            # Create optimizer objects and set defaults for all further optimizations.
            if parameters_dict["steepest_descent"]["use"]:
                sd = SteepestDescent(max_iterations=parameters_dict["steepest_descent"]["cycles"]) # Optimize with our custom optimizer.
                sd.optimize(atmsel, actions=modeller.optimizers.actions.trace(5))
            if parameters_dict["conjugate_gradients"]["use"]:
                cg = modeller.optimizers.conjugate_gradients(output='REPORT')
                # Run CG on the all-atom selection; write stats every 5 steps.
                cg.optimize(atmsel, max_iterations=parameters_dict["conjugate_gradients"]["cycles"], actions=modeller.optimizers.actions.trace(5, trcfil))
            if parameters_dict["quasi_newton"]["use"]:
                qn = modeller.optimizers.conjugate_gradients(output='REPORT')
                qn.optimize(atmsel, max_iterations=parameters_dict["quasi_newton"]["cycles"], actions=modeller.optimizers.actions.trace(5, trcfil))
            if parameters_dict["molecular_dynamics"]["use"]:
                md = modeller.optimizers.molecular_dynamics(output='REPORT')
                # Run MD; write out a PDB structure (called 'model_name.D9999xxxx.pdb')
                # every 10 steps during the run, and write stats every 10 steps.
                md.optimize(atmsel,
                    temperature=parameters_dict["molecular_dynamics"]["temperature"],
                    max_iterations=parameters_dict["molecular_dynamics"]["cycles"],
                    actions=modeller.optimizers.actions.trace(10, trcfil))
                    # actions=[modeller.optimizers.actions.write_structure(10, opt_code+'.D9999%04d.pdb'),
                    #          modeller.optimizers.actions.trace(10, trcfil)])

            mpdf = atmsel.energy()
            mdl.write(file=os.path.join(model_file_directory, opt_code+'.pdb'))

            env.edat.dynamic_lennard = old_dynamic_lennard
            env.edat.dynamic_coulomb = old_dynamic_coulomb
            env.edat.contact_shell = old_contact_shell
        #----------------------------------------------------

        #####################################################
        else:
            optimize_script_file_path = os.path.join(model_file_directory, "optimize.py")
            optimize_fh = open(optimize_script_file_path, "w")
            print >> optimize_fh, "import modeller, modeller.optimizers"
            print >> optimize_fh, "from modeller.scripts import complete_pdb"
            print >> optimize_fh, "env = modeller.environ()"
            if use_hetatm:
                print >> optimize_fh, "env.io.hetatm = True"
                if use_water:
                    print >> optimize_fh, "env.io.water = True"
            print >> optimize_fh, "env.edat.dynamic_sphere = True"
            print >> optimize_fh, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
            print >> optimize_fh, "env.libs.parameters.read(file='$(LIB)/par.lib')"
            print >> optimize_fh, 'code = "%s"' % opt_code
            print >> optimize_fh, 'env.edat.dynamic_coulomb = True'
            print >> optimize_fh, 'env.edat.dynamic_lennard = True'
            print >> optimize_fh, 'env.edat.contact_shell = %s' % parameters_dict["non_bonded_cutoff"]
            print >> optimize_fh, 'mdl = complete_pdb(env, "%s")' % os.path.join(model_file_directory, model_file_name)
            print >> optimize_fh, "mdl.write(file='%s')" % os.path.join(model_file_directory, opt_code+'.ini')
            print >> optimize_fh, "atmsel = modeller.selection(mdl)"
            # Generate the restraints:
            if parameters_dict["restraints"]["bond"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='bond', spline_on_site=False)"
            if parameters_dict["restraints"]["angle"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='angle', spline_on_site=False)"
            if parameters_dict["restraints"]["dihedral"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='dihedral', spline_on_site=False)"
            if parameters_dict["restraints"]["improper"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='improper', spline_on_site=False)"
            if parameters_dict["restraints"]["coulomb"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='coulomb', spline_on_site=False)"
            if parameters_dict["restraints"]["lj"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='lj', spline_on_site=False)"
            print >> optimize_fh, "mdl.restraints.write(file='%s')" % os.path.join(model_file_directory, opt_code+'.rsr')
            print >> optimize_fh, "mpdf = atmsel.energy()"
            print >> optimize_fh, 'class SteepestDescent(modeller.optimizers.state_optimizer):'
            print >> optimize_fh, '   """'
            print >> optimize_fh, '   Very simple steepest descent optimizer, in Python, as reported at:'
            print >> optimize_fh, '   http://www.salilab.org/modeller/9v4/manual/node252.html'
            print >> optimize_fh, '   """'
            print >> optimize_fh, "   _ok_keys = modeller.optimizers.state_optimizer._ok_keys + ('min_atom_shift', 'min_e_diff', 'step_size', 'max_iterations')"
            print >> optimize_fh, "   def __init__(self, step_size=0.0001, min_atom_shift=0.01, min_e_diff=1.0, max_iterations=None, **vars):"
            print >> optimize_fh, '       modeller.optimizers.state_optimizer.__init__(self, step_size=step_size,'
            print >> optimize_fh, '                                min_atom_shift=min_atom_shift,'
            print >> optimize_fh, '                                min_e_diff=min_e_diff,'
            print >> optimize_fh, '                                max_iterations=max_iterations, **vars)'
            print >> optimize_fh, "   def optimize(self, atmsel, **vars):"
            print >> optimize_fh, "        modeller.optimizers.state_optimizer.optimize(self, atmsel, **vars)"
            print >> optimize_fh, "        alpha = self.get_parameter('step_size')"
            print >> optimize_fh, "        minshift = self.get_parameter('min_atom_shift')"
            print >> optimize_fh, "        min_ediff = self.get_parameter('min_e_diff')"
            print >> optimize_fh, "        maxit = self.get_parameter('max_iterations')"
            print >> optimize_fh, "        state = self.get_state()"
            print >> optimize_fh, "        (olde, dstate) = self.energy(state)"
            print >> optimize_fh, "        while True:"
            print >> optimize_fh, "            for i in range(len(state)):"
            print >> optimize_fh, "                state[i] -= alpha * dstate[i]"
            print >> optimize_fh, "            (newe, dstate) = self.energy(state)"
            print >> optimize_fh, "            if abs(newe - olde) < min_ediff:"
            print >> optimize_fh, '                print "Finished at step" + str(self.step) + "due to energy criterion"'
            print >> optimize_fh, '                break'
            print >> optimize_fh, '            elif self.shiftmax < minshift:'
            print >> optimize_fh, '                print "Finished at step" + str(self.step) + "due to shift criterion"'
            print >> optimize_fh, '                break'
            print >> optimize_fh, '            elif maxit is not None and self.step >= maxit:'
            print >> optimize_fh, '                print "Finished at step" + str(self.step) + "due to step criterion"'
            print >> optimize_fh, '                break'
            print >> optimize_fh, '            if newe < olde:'
            print >> optimize_fh, '                alpha *= 2'
            print >> optimize_fh, '            else:'
            print >> optimize_fh, '                alpha /= 2'
            print >> optimize_fh, '            olde = newe'
            print >> optimize_fh, '            self.next_step()'
            print >> optimize_fh, '        self.finish()'
            print >> optimize_fh, "trcfil = file('%s','w')" % os.path.join(model_file_directory, opt_code+'.D00000001')
            if parameters_dict["steepest_descent"]["use"]:
                print >> optimize_fh, "sd = SteepestDescent(max_iterations=%s)" % parameters_dict["steepest_descent"]["cycles"]
                print >> optimize_fh, "sd.optimize(atmsel, actions=modeller.optimizers.actions.trace(5))"
            if parameters_dict["conjugate_gradients"]["use"]:
                print >> optimize_fh, "cg = modeller.optimizers.conjugate_gradients(output='REPORT')"
                print >> optimize_fh, "cg.optimize(atmsel, max_iterations=%s, actions=modeller.optimizers.actions.trace(5, trcfil))" % parameters_dict["conjugate_gradients"]["cycles"]
            if parameters_dict["quasi_newton"]["use"]:
                print >> optimize_fh, "qn = modeller.optimizers.conjugate_gradients(output='REPORT')"
                print >> optimize_fh, "qn.optimize(atmsel, max_iterations=%s, actions=modeller.optimizers.actions.trace(5, trcfil))" % parameters_dict["quasi_newton"]["cycles"]
            if parameters_dict["molecular_dynamics"]["use"]:
                print >> optimize_fh, "md = modeller.optimizers.molecular_dynamics(output='REPORT')"
                print >> optimize_fh, "md.optimize(atmsel, temperature=%s, max_iterations=%s, actions=modeller.optimizers.actions.trace(10, trcfil))" % (parameters_dict["molecular_dynamics"]["temperature"], parameters_dict["molecular_dynamics"]["cycles"])
            print >> optimize_fh, "mpdf = atmsel.energy()"
            print >> optimize_fh, "mdl.write(file='%s')" % os.path.join(model_file_directory, opt_code+'.pdb')
            optimize_fh.close()
            cline= "%s %s" % (self.modeller.get_exe_file_path(), optimize_script_file_path)
            self.execute_subprocess(cline, executing_modeller=True)
            os.remove(optimize_script_file_path)
        #####################################################

        # Checks if all the atomic coordinates of the refined structure are valid.
        if check_structure:
            optmized_structure_file_name = os.path.join(model_file_directory, opt_code+'.pdb')
            fh = open(optmized_structure_file_name, "rU")
            try:
                parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(optmized_structure_file_name, fh)
                fh.close()
            except:
                fh.close()
                return None

        return opt_code+'.pdb'


    #################################################################
    # Ramachandran plot.                                            #
    #################################################################

    def ramachandran_plot(self):
        """
        PROCHEK style Ramachandran Plot.
        """
        selected_sequences = self.get_selected_sequences()

        # If there is only one selected sequence and it has a structure loaded into PyMOL.
        if len(selected_sequences) == 1 and selected_sequences[0].has_structure():
            if not len(str(selected_sequences[0].my_sequence).replace('-','')):
                tkMessageBox.showerror("Selection Error",
                    "No residue for Ramachandran Plot generation")
                return

            PDB_file=[]
            title=''
            filename = os.path.join(self.structures_directory, selected_sequences[0].structure.chain_pdb_file_name_root + ".pdb")
            header = selected_sequences[0].my_header
            if os.path.isfile(filename):
                PDB_file.append(filename)
                if title:
                    title=title+", "+header
                else:
                    title=header

            if PDB_file:
                self.ramachandran_option(PDB_file,title)

        else:
            self.show_error_message("Selection Error","Please select one structure to display its Ramachandran Plot.")


    def ramachandran_option(self,PDB_file,title): # choose kind of aa to plot

        self.child=Toplevel(self.main_window)
        self.child.resizable(0,0)
        #self.child.geometry('400x500-10+40')
        self.child.title("<< Ramachandran Plot Options >>")
        self.child.config()
        try:
            self.child.grab_set()
        except:
            pass

        self.ch_main = Frame(self.child, background="black")
        self.ch_main.pack(expand = YES, fill = BOTH)

        self.upperframe = Frame(self.ch_main, borderwidth=5,
            background="black", relief="groove", pady=15)
        self.upperframe.pack(side = TOP, expand = NO, fill = X,
                              ipadx = 3, ipady = 3, pady=15)

        self.midframe = Frame(self.ch_main, background="black")
        self.midframe.pack(side=TOP,fill=BOTH,anchor="n",ipadx=5,ipady=5)

        self.lowerframe = Frame(self.ch_main, background="black")
        self.lowerframe.pack(side = BOTTOM, expand = NO, fill = Y,
            anchor="center", ipadx = 5, ipady = 5)

        self.mess=Label(self.upperframe, font = "comic 12", height = 1,
            text="Options for Ramachandran Plot",
            background="black", fg="white", pady = 2)
        self.mess.pack(fill="x")

        self.aa_sele_label=Label(self.midframe, font="comic 12", height=1,
            text= "Select Amino Acids", background="black", fg="red",
            borderwidth = 1, padx = 20)
        self.aa_sele_label.grid(row=0, column=0, sticky = W+E+N+S)

        def show_select_single_aa_frame():
            self.select_single_aa_frame.grid(row=2, column=1, sticky = "w")

        def hide_select_single_aa_frame():
            self.select_single_aa_frame.grid_remove()

        self.aa_sele_options_var = StringVar()
        self.aa_sele_options_var.set("all") # ['all','single']

        self.use_all_aa = Radiobutton(self.midframe,
            text="Use all amino acids", variable=self.aa_sele_options_var,
            value="all", background="black", foreground = "white",
            selectcolor = "red", highlightbackground="black",
            command=hide_select_single_aa_frame)
        self.use_all_aa.grid(row=0, column=1, sticky='w')

        self.select_single_aa = Radiobutton(self.midframe,
            text="Select amino acids", variable=self.aa_sele_options_var,
            value="single", background="black", foreground = "white",
            selectcolor = "red", highlightbackground="black",
            command=show_select_single_aa_frame)
        self.select_single_aa.grid(row=1, column=1, sticky='w')

        self.select_single_aa_frame=Frame(self.midframe,background="black")

        AA_one_letter_list='ACDEFGHIKLMNPQRSTVWY'
        self.aa_sele_var=dict()
        self.aa_checkbutton=dict()
        for i,aa in enumerate(AA_one_letter_list):
            self.aa_sele_var[aa]=IntVar()
            aa_freq=str(self.get_selected_sequences()[0].my_sequence
                ).count(aa)
            self.aa_checkbutton[aa]=Checkbutton(
                self.select_single_aa_frame,
                text=self.one2three(aa)+" ("+str(aa_freq)+")",
                variable=self.aa_sele_var[aa], background="black",
                foreground="white",selectcolor="red",
                highlightbackground="black")
            self.aa_checkbutton[aa].grid(row=i%10,column=i/10,sticky='w')
            # Only enable selection of aa present in primary sequence
            if not aa_freq:
                self.aa_checkbutton[aa].config(state=DISABLED)


        def state():
            AA_list=None
            title_append=''
            if self.aa_sele_options_var.get()=="single":
                AA_list=''
                for aa in AA_one_letter_list:
                    if self.aa_sele_var[aa].get():
                        AA_list+=aa
                if AA_list:
                    title_append=" (Amino Acid: "+AA_list+")"
                else:
                    tkMessageBox.showerror("Selection Error",
                        "No residue for Ramachandran Plot generation")
                    return
            pmsp.ramachandran(PDB_file,title+title_append,AA_list=AA_list)
            self.child.destroy()

        self.submit=Button(self.lowerframe, text="SUBMIT", command=state,
            relief="raised", borderwidth="3", bg="black", fg="white")
        self.submit.pack()


    #################################################################
    # Secondary structure assignment.                               #
    #################################################################

    def assign_secondary_structure(self, element):
        if element.has_structure():
            if hasattr(self, "ksdssp") and self.ksdssp.exe_exists():
                self.assign_with_ksdssp(element)
            else:
                self.assign_with_pymol_dss(element)


    def assign_with_ksdssp(self, element):
        # Runs ksdssp.
        dssptext=pmsp.runKSDSSP(os.path.join(self.structures_directory, element.structure.chain_pdb_file_name), ksdssp_exe=self.ksdssp.get_exe_file_path())
        # Parses ksdssp's output, that is, an series of pdb format 'HELIX' and 'SHEET' record lines.
        dsspout = dssptext.split("\n")
        helices = set() # A set to store the sequence numbers of the residues in helical conformation.
        sheets = set() # A set to store the sequence numbers of the residues in sheet conformation.
        for line in dsspout:
            if line.startswith("HELIX"):
                new_residues_set = set(range(int(line[21:25]), int(line[33:37])+1))
                helices.update(new_residues_set)
            elif line.startswith("SHEET"):
                new_residues_set = set(range(int(line[22:26]), int(line[33:37])+1))
                sheets.update(new_residues_set)
        # Assigns to the PyMod element the observed secondaey structure observed using ksdssp.
        element.pymol_dss_list = []
        for residue in element.structure.get_all_residues_list():
            if residue.pdb_position in helices:
                element.pymol_dss_list.append("H")
                rsel = element.build_residue_selector_for_pymol(residue.pdb_position)
                cmd.alter(rsel,"ss='H'") # Set the residue new conformation in PyMOL.
            elif residue.pdb_position in sheets:
                element.pymol_dss_list.append("S")
                rsel = element.build_residue_selector_for_pymol(residue.pdb_position)
                cmd.alter(rsel,"ss='S'") # Set the residue new conformation in PyMOL.
            else:
                element.pymol_dss_list.append("L")
                rsel = element.build_residue_selector_for_pymol(residue.pdb_position)
                cmd.alter(rsel,"ss='L'") # Set the residue new conformation in PyMOL.
        # Updated PyMOL.
        cmd.rebuild()


    def assign_with_pymol_dss(self, element):
        """
        Uses PyMOL's DSS algorithm to assign the secondary structure to a sequence according to atom
        coordinates of its PDB file.
        """
        selection = "object %s and n. CA" % element.build_chain_selector_for_pymol()
        stored.resi_set = set()
        stored.temp_sec_str = []
        stored.pymol_info = []
        stored.pymod_resi_set = set([res.pdb_position for res in element.structure.get_all_residues_list()])
        def include_sec_str_val(ca_tuple):
            if not ca_tuple[1] in stored.resi_set and ca_tuple[1] in stored.pymod_resi_set:
                stored.temp_sec_str.append(ca_tuple[0])
                stored.resi_set.add(ca_tuple[1])
                stored.pymol_info.append(ca_tuple)
        stored.include_val = include_sec_str_val
        cmd.iterate(selection, "stored.include_val((ss, resv))")
        # print stored.pymol_info
        # print [res.pdb_position for res in element.structure.get_all_residues_list()]
        element.pymol_dss_list = list(stored.temp_sec_str)
        if not (len(element.pymol_dss_list) == len(element.structure.get_all_residues_list())):
            pass


    #################################################################
    # Run PSIPRED.                                                  #
    #################################################################

    def launch_psipred_from_main_menu(self):
        """
        Called when users decide to predict the secondary structure of a sequence using PSIPRED.
        """
        selection = self.get_selected_sequences()
        if len(selection) == 0:
            self.show_error_message("PSIPRED Error", "Please select at least one sequence to be analyzed with PSIPRED.")
            return False
        if self.check_psipred_parameters():
            for sequence in selection:
                # Actually calls the method that launches PSIPRED.
                prediction_successful = self.run_psipred(sequence)
                if prediction_successful:
                    sequence.color_by = "secondary-predicted"
                    sequence.color_element(on_grid=False,color_pdb=True)


    def check_psipred_parameters(self): # predict_secondary_structure(self, elements=None):
        """
        Checks that the files needed to run PSIPRED exists on users' machines.
        """
        # First checks for PSIPRED installation.
        if not self.psipred["exe_dir_path"].path_exists():
            self.psipred.exe_not_found()
            return False

        # Then checks for PSIPRED datafiles.
        if not self.psipred["data_dir_path"].path_exists():
            title = "PSIPRED error"
            message = "PSIPRED 'data' directory not found! Please specify it in the PSIPRED options in the options window of PyMod."
            self.show_error_message(title,message)
            return False

        # Checks for PSI-BLAST on the user's system.
        if not self.blast_plus["exe_dir_path"].path_exists():
            self.blast_plus.exe_not_found()
            return False

        # And finally checks for a BLAST database.
        if not self.psipred["database_dir_path"].path_exists():
            self.show_error_message("PSIPRED error", "A directory containing a BLAST database was not found! Please specify it in the PSIPRED options in the options window of PyMod.")
            return False

        dbpath = self.psipred["database_dir_path"].get_value()
        if not pmos.verify_valid_blast_dbdir(dbpath):
            self.show_error_message("PSIPRED Error", "The database '%s' directory does not contain a valid set of database files." % (dbpath))
            return False

        return True


    def run_psipred(self, element):
        """
        Actually runs PSIPRED, collects its results and map them on the sequences in PyMod main
        window.
        """
        print_output = True
        sequence_header = element.my_header
        if print_output:
            print "Beginning PSIPRED prediction for:", sequence_header

        # The name of the BLAST database file.
        # If the database files are contained in a folder like this: /home/user/pymod/databases/swissprot/swissprot
        dbpath = self.psipred["database_dir_path"].get_value() # e.g.: /home/user/pymod/databases/swissprot
        dbprefix = pmos.get_blast_database_prefix(dbpath) # e.g.: swissprot
        if print_output:
            print "dbpath:", dbpath

        # Where the NCBI programs have been installed.
        ncbidir = self.blast_plus["exe_dir_path"].get_value()
        if print_output:
            print "ncbidir:", ncbidir

        # Where the PSIPRED V2 programs have been installed.
        execdir = self.psipred["exe_dir_path"].get_value()
        if print_output:
            print "execdir:", execdir

        # Where the PSIPRED V2 data files have been installed.
        datadir = self.psipred["data_dir_path"].get_value()
        if print_output:
            print "datadir",datadir

        # Write the temporary input fasta file, setting its basename.
        basename = "psipred_temp"
        if print_output:
            print "basename: ", basename
        self.build_sequences_file([element], basename, file_format="fasta", remove_indels=True, new_directory=self.psipred_directory)

        # ---
        # Execute PSI-BLAST.
        # ---
        if print_output:
            print "Running PSI-BLAST with sequence", basename ,"..."
        try:
            self.execute_psiblast(
                ncbi_dir = ncbidir,
                db_path = dbpath,
                query = os.path.join(self.psipred_directory, basename+".fasta"),
                inclusion_ethresh = 0.001,
                out_pssm = os.path.join(self.psipred_directory, basename+".chk"),
                out = os.path.join(self.psipred_directory, basename+".blast"),
                num_iterations = 3,
                num_alignments = 0)
            # psiblast_output = open("%s.blast" % os.path.join(self.psipred_directory, basename),"w")
            # self.execute_subprocess(psiblast_command, new_stdout=psiblast_output)
            # psiblast_output.close()

        except:
            if print_output:
                print "FATAL: Error whilst running psiblast - script terminated!"
            self.show_error_message("PSIPRED Error", "There was an error while running PSI-BLAST, so PSIPRED cannot perform a prediction for %s." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        # ---
        # Execute chkparse.
        # ---
        # !WORKING! The problem is here.
        if print_output:
            print "Predicting secondary structure..."
        # query = pmos.build_commandline_file_argument("query", query, "fasta")
        chkdir_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("chkparse"))) + " " +
                          pmos.build_commandline_path_string("%s.chk" % os.path.join(self.psipred_directory, basename)))
        try:
            chkdir_output = open("%s.mtx" % os.path.join(self.psipred_directory, basename),"w")
            self.execute_subprocess(chkdir_command, new_stdout=chkdir_output)
            chkdir_output.close()
        except:
            if print_output:
                print "FATAL: Error whilst running chkdir - script terminated!"
            self.show_error_message("PSIPRED Error", "No homologous sequences were found by PSI-BLAST for %s, so PSIPRED cannot perform a prediction for this sequence." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        # ---
        # Execute PSIPRED pass 1.
        # ---
        psipass1_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipred"))) + " " +
                            pmos.build_commandline_path_string("%s.mtx" % os.path.join(self.psipred_directory, basename)) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat")) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat2")) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat3")))
        try:
            psipass1_output = open("%s.ss" % os.path.join(self.psipred_directory, basename),"w")
            self.execute_subprocess(psipass1_command, new_stdout=psipass1_output)
            psipass1_output.close()
        except:
            if print_output:
                print "FATAL: Error whilst running psipred 1 - script terminated!"
            self.show_error_message("PSIPRED Error", "There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        # ---
        # Execute PSIPRED pass 2.
        # ---
        psipass2_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipass2"))) + " " +
                            "%s 1 1.0 1.0" % pmos.build_commandline_path_string(os.path.join(datadir,"weights_p2.dat")) + " " +
                            pmos.build_commandline_path_string("%s.ss2" % os.path.join(self.psipred_directory, basename)) + " " +
                            pmos.build_commandline_path_string("%s.ss" % os.path.join(self.psipred_directory, basename)))
        try:
            psipass2_output = open("%s.horiz" % os.path.join(self.psipred_directory, basename),"w")
            self.execute_subprocess(psipass2_command, new_stdout=psipass2_output)
            psipass2_output.close()
        except:
            if print_output:
                print "FATAL: Error whilst running psipass 2 - script terminated!"
            self.show_error_message("PSIPRED Error", "There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        # ---
        # Clean up PSIPRED files.
        # ---
        if print_output:
            print "Cleaning up ..."

        # Remove temporary files.
        self.remove_psipred_temp_files()

        # Renames the output files.
        output_files_name = pmos.clean_file_name(element.my_header)
        for ext in pmdt.psipred_output_extensions:
            os.rename(os.path.join(self.psipred_directory, basename+ext),
                      os.path.join(self.psipred_directory, output_files_name+ext))

        if print_output:
            print "Final output files:" + output_files_name + ".ss2 " + output_files_name + ".horiz"
        print "Finished."

        # ---
        # Parses the results from .horiz output file.
        # ---
        results_file = open(os.path.join(self.psipred_directory, output_files_name+".horiz"),"r")
        confs = "" # String for confidence scores of each residue.
        preds = "" # String for the secondary structure elements prediction of each residue.
        for l in results_file.readlines():
            if l.startswith("Conf:"):
                rl = l[6:66].replace(" ","").replace("\n","")
                confs += rl
            elif l.startswith("Pred:"):
                rl = l[6:66].replace(" ","").replace("\n","")
                preds += rl
        results_file.close()

        # Actually stores in the PyMod elements the results.
        element.psipred_elements_list = []
        for c, e in zip(confs, preds):
            element.psipred_elements_list.append({"confidence":int(c),"sec-str-element":e})

        return True


    def remove_psipred_temp_files(self):
        try:
            files_to_remove = filter(lambda f: not os.path.splitext(f)[1] in pmdt.psipred_output_extensions, os.listdir(self.psipred_directory))
            map(lambda f: os.remove(os.path.join(self.psipred_directory,f)) , files_to_remove)
        except:
            pass

    #################################################################
    # Superpose.                                                    #
    #################################################################

    def superpose(self):
        """
        Called from the main menu. This will superpose to a 'fixed' structure (the first one in the
        selection) one or more 'mobile' structures.
        """
        correct_selection = False
        structures_to_superpose = self.get_selected_sequences()
        if len(structures_to_superpose) >= 2:
            if not False in [e.has_structure() for e in structures_to_superpose]:
                correct_selection = True
        if correct_selection:
            for i in range(1, len(structures_to_superpose)):
                sel1 = structures_to_superpose[0].build_chain_selector_for_pymol()
                sel2 = structures_to_superpose[i].build_chain_selector_for_pymol()
                self.superpose_in_pymol(sel2, sel1)
        else:
            self.show_error_message("Selection Error","Please select at least two structures before superposing")


    def superpose_in_pymol(self, selector_1, selector_2, save_superposed_structure=True):
        """
        align mobile, target
        """
        if hasattr(cmd,"super"): # super is sequence-independent
            cmd.super(selector_1, selector_2)
            # cmd.cealign(target=selector_1, mobile=selector_2)
        else: # PyMOL 0.99 does not have cmd.super
            cmd.align(selector_1, selector_2)
        if save_superposed_structure:
            self.pymol_save(os.path.join(self.structures_directory, selector_1+".pdb"), selector_1)



    ################################################################################################
    # HOMOLOGY MODELING.                                                                           #
    ################################################################################################

    def launch_modeller_from_main_menu(self):
        """
        This method is called when the "MODELLER" option is clicked in the "Tools" menu.
        """
        # Try to find if Modeller is installed on the user's computer.
        if self.modeller.can_be_launched():

            # Get the selected sequences to see if there is a correct selection.
            selected_sequences = self.get_selected_sequences()

            # First check if at least one sequence is selected.
            if len(selected_sequences) > 0:

                # Checks if all the selected sequences can be used to build a model.
                if not False in [s.can_be_modeled() for s in selected_sequences]:

                    # Checks that all the selected sequences are currently aligned to some other
                    # sequence (aligned sequences should always have .is_child == True). Only
                    # sequences aligned to some template can be modeled.
                    if not False in [e.is_child for e in selected_sequences]:

                        # Using the 'build_cluster_list()' method build the
                        # 'self.involved_cluster_elements_list' just like when performing an alignment.
                        self.build_cluster_lists()

                        # This will contain a list of 'Modeling_cluster' objects. These objects will
                        # be used from now on to store all the informations on who to build models.
                        self.modeling_clusters_list = []

                        # Checks other conditions in each cluster (that is, checks if there is only
                        # one target sequence selected per cluster and if it is aligned to some
                        # suitable template).
                        correct_selection_in_clusters = False
                        for cluster_element in self.involved_cluster_elements_list:
                            if self.check_correct_modeling_selection_in_cluster(cluster_element):
                                # If the selection for the cluster is correct, it builds a
                                # 'Modeling_cluster' object and stores it in the
                                # 'self.modeling_clusters_list'.
                                mco = Modeling_cluster(cluster_element)
                                self.modeling_clusters_list.append(mco)
                                correct_selection_in_clusters = True
                            else:
                                correct_selection_in_clusters = False
                                # Stops if there is just one cluster with an incorrect selection.
                                break

                        # Proceeds only if every modeling clusters has a correct selection.
                        if correct_selection_in_clusters:
                            # If there is only one cluster involved.
                            if len(self.modeling_clusters_list) == 1:
                                self.build_modeling_window()

                            # Multiple chains modeling requires the identification of "template
                            # complexes" and additional controls.
                            else:
                                # This will build the 'self.template_complex_list'.
                                self.initialize_multichain_modeling()
                                # Proceeds only if there are is at least one suitable "template
                                # complex".
                                if len(self.template_complex_list) > 0:
                                    # Finally builds the modeling window.
                                    self.build_modeling_window()
                                else:
                                    title = "Selection Error"
                                    message = "There isn't any suitable 'Template Complexes' to perform multiple chain homology modeling."
                                    self.show_error_message(title,message)
                    else:
                        title = "Selection Error"
                        message = "Please select only target sequences that are currently aligned to some structure."
                        self.show_error_message(title,message)
                else:
                    title = "Selection Error"
                    message = "Please select only sequences that do not have a structure loaded in PyMOL."
                    self.show_error_message(title,message)
            else:
                title = "Selection Error"
                message = "Please select at least one target sequence to use Modeller."
                self.show_error_message(title,message)

        # If MODELLER is not istalled.
        else:
            self.modeller.exe_not_found()


    def check_correct_modeling_selection_in_cluster(self,cluster_element):
        """
        This will check if there is only one target sequence selected per cluster and if it is
        currently aligned to some suitable template.
        """
        # This will be set as True only if all the necessaries conditions are met.
        correct_selection = False

        # Checks that only one sequence per cluster is selected.
        if not self.check_only_one_selected_child_per_cluster(cluster_element):
            title = "Selection Error"
            message = "Please select only one target sequence in the following cluster: %s" % (cluster_element.my_header)
            self.show_error_message(title,message)
            return False

        # This will be needed to inform the user about which cluster has an incorrect selection.
        target_name = None
        templates_temp_list = []
        # Look if there is at least one suitable template aligned to the target sequence.
        for sequence in self.get_children(cluster_element):
            if not sequence.selected and sequence.has_structure() and sequence.is_model != True:
                templates_temp_list.append(sequence)
            # Takes the name of the target sequence.
            if sequence.selected:
                target_name = sequence.my_header

        # Checks if some templates have been found.
        if len(templates_temp_list) > 0:
            # Checks the presence of nucleic acids templates: currently they are not
            # supported by PyMOd.
            if "nucleic-acid" in [t.polymer_type for t in templates_temp_list]:
                title = "Selection Error"
                message = "Template '%s' is a nucleic acid chain. PyMod currently does not support nucleic acid templates, so the modelization cannot be performed." % (t.my_header)
                self.show_error_message(title,message)
                return False
            else:
                return True
        else:
            title = "Selection Error"
            message = "The target sequence %s in the following cluster is currently not aligned to any PDB structure." % (target_name)
            self.show_error_message(title,message)
            return False


    def initialize_multichain_modeling(self):
        """
        This method will prepare data needed to perform multichain modeling. It will:
            - identify suitable template complexes
            - check if there are target sequences with the same sequence, so that symmetry restraints
              can be applied to them when using Modeller.
        """
        # ---
        # Generates modeling clusters dictionaries.
        # ---
        # They will be needed to check if a suitable "template complex" can be used. A cluster with
        # the following templates:
        #     - 1HHO_Chain:A, 2DN2_Chain:A, 2DN2_Chain:B
        # will generate a dictionary with the following structure:
        #     - {"1HHO.pdb":1, "2DN2.pdb":2}
        # The keys are the original PDB files, and the values are the number of chains in the cluster
        # which belong to that PDB structure.
        for mc in self.modeling_clusters_list:
            for t in mc.structure_list:
                if t.structure.original_pdb_file_name in mc.dictionary.keys():
                    mc.dictionary[t.structure.original_pdb_file_name] += 1
                else:
                    mc.dictionary.update({t.structure.original_pdb_file_name:1})

        # ---
        # Checks if there are suitable "template complexes".
        # ---
        # A "teplate complex" is available only if in each selected cluster there is at least ONE
        # chain coming from the same original PDB file. For example, with these two cluster:
        #     - cluster 1: <1HHO_Chain:A>, 2DN2_Chain:A
        #     - cluster 2: <1HHO_Chain:B>, 3EOK_Chain:A
        # the "template complex" is 1HHO.
        codes_per_cluster = []
        for mc in self.modeling_clusters_list:
            codes_per_cluster.append(set([k for k in mc.dictionary.keys()]))
        self.template_complex_list = list(set.intersection(*codes_per_cluster))
        self.template_complex_list.sort() # Sorts the list alphabetically.

        # ---
        # Checks if there are some target chains with the same sequence, so that the user may apply
        # symmetry restraints to them.
        # ---
        # Builds a Symmetry_restraints_groups object that is going to be used to keep track of
        # modeling clusters that have a target sequence with the same sequence.
        self.symmetry_restraints_groups = Symmetry_restraints_groups_list()
        for mc in self.modeling_clusters_list:
            seq = str(mc.target.my_sequence).replace("-","")
            # Adds a "symmetry restraints group" for each group of target sequences that share the
            # exact same sequence.
            if not seq in [g.id for g in self.symmetry_restraints_groups.get_groups()]:
                self.symmetry_restraints_groups.add_group(seq)
            self.symmetry_restraints_groups.get_group_by_id(seq).add_cluster(mc)
        # Also assigns "symmetry ids" to each modeling cluster.
        for mc in self.modeling_clusters_list:
            seq = str(mc.target.my_sequence).replace("-","")
            if seq in [g.id for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2)]:
                mc.set_symmetry_id(seq)
            else:
                mc.set_symmetry_id(None)


    def build_modeling_window(self):
        """
        Builds the modeling window.
        """
        # Window with all the options for Modeller.
        self.modeling_window=Toplevel(self.main_window)
        self.modeling_window.resizable(1,1)
        self.modeling_window.title("<< MODELLER Options >>")
        self.modeling_window.config()
        try:
            self.modeling_window.grab_set()
        except:
            pass
        # Frame that occupies all the window. It is going to contain 3 frames: upper, middle and
        # lower.
        self.ch_main = Frame(self.modeling_window, background='black')
        self.ch_main.pack(expand = YES, fill = BOTH)

        # Builds the upper frame with the title.
        self.upperframe = Frame(self.ch_main, borderwidth=5, background='black', relief='groove', pady=15)
        self.upperframe.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3, pady=15)
        self.mess=Label(self.upperframe, text= "Here you can modify options for MODELLER", **pmgi.label_style_0)
        self.mess.pack(fill="x")

        # Builds the middle frame where there are going to be a Notebook and its tabs.
        self.midframe = Frame(self.ch_main, background='red')
        self.midframe.pack(side = TOP, fill = BOTH, expand =1)
        # Notebook in the middle frame
        self.notebook = Pmw.NoteBook(self.midframe, borderwidth = 2)
        self.notebook.pack(fill = BOTH, expand = 1)
        # Configures the hull (the Canvas that contains the whole Notebook). This will determine the
        # size of the Notebook Pages.
        self.notebook.component("hull").configure(background="black",width=680,height=500)
        # Builds the different Notebook pages.
        self.build_main_page()
        self.build_disulfides_page()
        self.build_options_page()

        # Builds the lower frame of the modeling window where the "SUBMIT" button is.
        self.lowerframe = Frame(self.ch_main, background='black',bd=2,relief=GROOVE)
        self.lowerframe.pack(side = BOTTOM,expand=1,fill=BOTH, anchor="center",ipadx=5,ipady=5)
        # This is the button on the modellization window that when is pressed calls
        # the state() function above.
        self.submit=Button(self.lowerframe, text="SUBMIT", command=self.perform_modelization, **pmgi.button_style_1)
        self.submit.pack(pady=10)


    def build_main_page(self):
        """
        Add and configure the 'Main' page and tab fo the modeling window.
        """
        self.templates_page = self.notebook.add('Main')
        self.notebook.tab('Main').focus_set()
        self.notebook.page(0).configure(bg="black")
        self.notebook.tab(0).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
        # Frame with a scrollbar for the templates.
        self.templates_frame = Pmw.ScrolledFrame(
            self.templates_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
            horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
            borderframe = 1, usehullsize = 1, frame_background = 'black')
        # Change the border of the frame, bd=2 looks bad.
        self.templates_frame.component("borderframe").configure(bd=1)
        self.templates_frame.pack(side= TOP, fill = 'both', expand = 1)
        # This is the actual Frame where the content of the tab is going to be packed.
        self.templates_frame_interior = self.templates_frame.interior()
        self.templates_frame_interior.configure(bd=0,pady=20)

        # ---
        # Starts to insert content in the "Main" page.
        # ---

        # If the user choose to build a multiple chain model, it displays an additional option to
        # let the user choose his/her "template complex".
        if len(self.modeling_clusters_list) > 1:
            # An additional frame for the "template complex" selection.
            self.template_complex_selection_frame = Frame(self.templates_frame_interior,borderwidth=0, background='black', relief='groove', pady=15)
            self.template_complex_selection_frame.pack(side="top", anchor="w")
            self.template_complex_selection_label=Label(self.template_complex_selection_frame, text= "Template Complex selection: ", **pmgi.modeling_window_title_style)
            self.template_complex_selection_label.grid(row=0, column=0,sticky = W+N)

            # The user can choose the "template complex" with some Radiobuttons.
            self.template_complex_var = StringVar()
            # Initialize by default with the first PDB in the list.
            self.template_complex_var.set(self.template_complex_list[0])

            # Display some information to explain what is a "template complex".
            information = (
            "Select the PDB file containing the complex on which you would like to base the building\n"+
            "of your multiple chain model. The relative orientation in space and the interfaces of\n"+
            "your model's chains will be based on the architecture of the Template Complex.")

            self.template_complex_message = Label(self.template_complex_selection_frame, text= information, **pmgi.modeling_window_explanation)
            self.template_complex_message.grid(row=1, column=0, sticky = "w")

            for (tc_i,tc) in enumerate(self.template_complex_list):
                tcb = Radiobutton(self.template_complex_selection_frame, text=tc, variable=self.template_complex_var, value=tc, **pmgi.modeling_window_rb_big)
                tcb.grid(row=tc_i+2, column=0, sticky = "w",padx=(20,0))

        # Builds a frame for each modeling_cluster.
        for (i, modeling_cluster) in enumerate(self.modeling_clusters_list):

            if len(self.modeling_clusters_list) > 1:
                spacer_frame = Frame(self.templates_frame_interior, background='black',height = 2,bd=1,relief=GROOVE)
                spacer_frame.pack(side="top", padx = 20, anchor="w", fill="x")

            # A frame that will contain all the widgets necessary to choose the templates for a
            # single target sequence.
            modeling_cluster_frame = Frame(self.templates_frame_interior, borderwidth=0, background='black', relief='groove', pady=5, padx=0)
            modeling_cluster_frame.pack(side="top", anchor="w", pady=(0,10))

            ######################################################
            # This frame should contain also other options like: #
            #     - loop refinement                              #
            #     - secondary structure assignment to the model  #
            #     - others...                                    #
            ######################################################
            modeling_option_label = Label(modeling_cluster_frame, text= "Modeling options for target: %s" % (modeling_cluster.target_name), **pmgi.modeling_window_title_style)
            modeling_option_label.pack(side="top", anchor="w")

            additional_options_label=Label(modeling_cluster_frame, text= "Restraints options", **pmgi.modeling_options_sections_style)
            additional_options_frame = Frame(modeling_cluster_frame, **pmgi.target_box_style)
            show_additional_options = False
            if len(self.modeling_clusters_list) > 1:
                # Use symmetry restraints option.
                if modeling_cluster.symmetry_id != None:
                    symmetry_frame = Frame(additional_options_frame,background='black',bd=0,relief=GROOVE)
                    symmetry_frame.pack(side=LEFT)

                    symmetry_label = Label(symmetry_frame, text= "Use simmetry restraints for this chain:", **pmgi.modeling_window_option_style)
                    symmetry_label.grid(row=0, column=0,sticky= N+W)
                    use_symmetry_var = IntVar()
                    symmetry_chk = Checkbutton(symmetry_frame, text="", variable=use_symmetry_var, **pmgi.modeling_window_checkbutton)
                    symmetry_chk.grid(row=0, column=1,sticky= N+W)

                    symmetry_information = "Show Info"
                    symmetry_info = Button(symmetry_frame, text=symmetry_information, command= lambda: self.show_symmetry_info(modeling_cluster), relief="raised",borderwidth=0, bg="black", highlightbackground='black', fg="white", pady = 0, anchor = "w")
                    symmetry_info.grid(row=0, column=2,sticky= N+W)
                    modeling_cluster.set_symmetry_var(use_symmetry_var)

                    show_additional_options = True

            if show_additional_options:
                additional_options_label.pack(side="top", anchor="w")
                additional_options_frame.pack(side="top", anchor="w", padx = (30,0),pady=(0,5))

            # Builds a frame for each structure aligned to the target sequence of the current
            # modeling cluster.
            template_label=Label(modeling_cluster_frame, text= "Template selection", **pmgi.modeling_options_sections_style)
            template_label.pack(side="top", anchor="w")
            modeling_cluster.structure_frame_list = []
            for (si,structure) in enumerate(modeling_cluster.structure_list):
                # This object is not a tkinter one, but contains as attributes many of them.
                structure_frame = pmgi.Structure_frame(self, structure,modeling_cluster.target,modeling_cluster_frame,si,i)
                # Builds a frame for each template structure.
                structure_frame.build_frame()
                # Append the current "structure_frame" to the list of the current modeling cluster.
                # Storing this object will also store the value of each Checkbox, Radiobutton and
                # entry found inside it.
                modeling_cluster.structure_frame_list.append(structure_frame)


    def switch_all_hetres_checkbutton_states(self,het_radio_button_state):
        """
        Launched when the user activates/inactivates the "Include HetAtoms" in the Options page in
        the modeling window.
        """
        for mc in self.modeling_clusters_list:
            mc.switch_hetres_checkbutton_states(het_radio_button_state)


    def show_symmetry_info(self, modeling_cluster):
        """
        Displays informations about which target sequence shares the same sequence of other targets.
        """
        mc_list = self.symmetry_restraints_groups.get_group_by_id(modeling_cluster.symmetry_id).list_of_clusters
        mc_list = filter(lambda x: not x is modeling_cluster ,mc_list)
        message1 = "The target '%s' shares the same sequence with these other targets:" % (modeling_cluster.target_name)
        seqs = reduce(lambda x,y: x+",\n"+y, [mc.target_name for mc in mc_list])
        message2 = "so you may apply symmetry restraints for them."
        tkMessageBox.showinfo("Symmetry restraints information", message1 + "\n\n" + seqs + "\n\n" + message2, parent=self.modeling_window)


    def build_disulfides_page(self):
        """
        Add the "Disulfides" page to the modeling window notebook.
        """
        # Disulfide page.
        self.disulfides_bridges_page = self.notebook.add('Disulfides')
        self.notebook.page(1).configure(bg="black")
        self.notebook.tab(1).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
        # Frame with a scrollbar for the disulfides options.
        self.disulfides_scrolled_frame = Pmw.ScrolledFrame(
            self.disulfides_bridges_page, vscrollmode = "static", hscrollmode = "dynamic",
            horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
            borderframe = 1, usehullsize = 1, frame_background = 'black')
        # Same as for the Frame for the templates above.
        self.disulfides_scrolled_frame.component("borderframe").configure(bd=1)
        self.disulfides_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
        self.disulfides_container = self.disulfides_scrolled_frame.interior()
        self.disulfides_container.configure(bd=0,pady=20)

        # Part for the "Disulfides" page.
        self.disulfides_frame = pmgi.Disulfides_frame(self, self.disulfides_container)
        # If at least one cluster has a target with at least two CYS residues, then build the
        # disulfide page with all its options.
        if self.check_targets_with_cys():
            self.disulfides_frame.build_template_dsb_frame()
            # User defined dsb. Each target is going to have a frame to define additional dsb.
            self.disulfides_frame.build_user_defined_dsb_frame()
            for (mci,mc) in enumerate(self.modeling_clusters_list):
                self.disulfides_frame.build_modeling_cluster_users_dsb_frame(mc,mci)
            self.disulfides_frame.build_auto_dsb_frame()
        else:
            self.disulfides_frame.build_no_dsb_frame()


    def check_targets_with_cys(self):
        """
        Check if there is at least one modeling cluster with a target sequence with at least two CYS
        residues.
        """
        if True in [mc.target_with_cys for mc in self.modeling_clusters_list]:
           return True
        else:
            return False


    def build_options_page(self):
        """
        Add the "Options" page on modeling window notebook.
        """
        # Options page.
        self.options_page = self.notebook.add('Options')
        self.notebook.page(2).configure(bg="black")
        self.notebook.tab(2).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray25")
        # Frame with a scrollbar for the options.
        self.options_scrolled_frame = Pmw.ScrolledFrame(
            self.options_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
            horizflex = 'expand', vertflex = 'expand',  hull_borderwidth = 0,
            borderframe = 1, usehullsize = 1, frame_background = 'black')
        # Same as for the Frame for the templates above.
        self.options_scrolled_frame.component("borderframe").configure(bd=1)
        self.options_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
        self.options_frame = self.options_scrolled_frame.interior()
        self.options_frame.configure(bd=0,pady=20)

        grid_widgets = True
        self.options_frame_grid_options = {"padx": 10, "pady": 10, "sticky": "w"}
        # Start to insert modeling options widgets.
        option_widgets_to_align = []
        # Option to chose the number of models that Modeller has to produce.
        self.max_models_enf = pmgi.PyMod_entryfield(
            self.options_frame, label_text = "Models to Build", value = 1,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : self.max_models_per_session})
        if grid_widgets:
            self.max_models_enf.grid(row=0, **self.options_frame_grid_options)
        else:
            self.max_models_enf.pack(**pmgi.pack_options_1)

        option_widgets_to_align.append(self.max_models_enf)

        # Option to choose if Modeller is going to include HETATMs.
        self.exclude_heteroatoms_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Exclude Heteroatoms')
        for choice in ("Yes", "No"):
            self.exclude_heteroatoms_rds.add(choice)
        self.exclude_heteroatoms_rds.setvalue("No")
        if grid_widgets:
            self.exclude_heteroatoms_rds.grid(row=1, **self.options_frame_grid_options)
        else:
            self.exclude_heteroatoms_rds.pack(**pmgi.pack_options_1)
        option_widgets_to_align.append(self.exclude_heteroatoms_rds)
        self.exclude_heteroatoms_rds.button(0).configure(command=lambda: self.switch_all_hetres_checkbutton_states(0)) # Yes, inactivate.
        self.exclude_heteroatoms_rds.button(1).configure(command=lambda: self.switch_all_hetres_checkbutton_states(1)) # No, activate.

        # Option to choose the level of optimization for Modeller.
        self.optimization_level_choices = ("Low", "Default", "Mid", "High")
        self.optimization_level_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Optimization Level')
        for choice in self.optimization_level_choices:
            self.optimization_level_rds.add(choice)
        self.optimization_level_rds.setvalue("Default")
        if grid_widgets:
            self.optimization_level_rds.grid(row=2, **self.options_frame_grid_options)
        else:
            self.optimization_level_rds.pack(**pmgi.pack_options_1)
        option_widgets_to_align.append(self.optimization_level_rds)

        # Option to choose the level of additional energy optimization.
        self.energy_minimization_choices = ("None", "Use")
        self.energy_minimization_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Additional Energy Minimization')
        for choice in self.energy_minimization_choices:
            self.energy_minimization_rds.add(choice)
        self.energy_minimization_rds.setvalue("None")
        if grid_widgets:
            self.energy_minimization_rds.grid(row=3, **self.options_frame_grid_options)
        else:
            self.energy_minimization_rds.pack(**pmgi.pack_options_1)
        self.energy_minimization_rds.button(0).configure(command=self.hide_energy_minimization_frame)
        self.energy_minimization_rds.button(1).configure(command=self.show_energy_minimization_frame)
        option_widgets_to_align.append(self.energy_minimization_rds)
        self.energy_minimization_frame = pmgi.Energy_minimization_frame(self.options_frame)
        # if grid_widgets:
        #     self.energy_minimization_frame.grid(row=4, **self.options_frame_grid_options)
        # else:
        #     self.energy_minimization_frame.pack()

        # Option to choose the way to color the models.
        self.color_models_choices = ("Default", "DOPE Score") # Delta DOPE e b-factor
        self.color_models_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Color Models by')
        for choice in self.color_models_choices:
            self.color_models_rds.add(choice)
        self.color_models_rds.setvalue("Default")
        if grid_widgets:
            self.color_models_rds.grid(row=5, **self.options_frame_grid_options)
        else:
            self.color_models_rds.pack(**pmgi.pack_options_1)
        option_widgets_to_align.append(self.color_models_rds)

        # Option to choose whether to super models to template.
        self.superpose_models_to_templates_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Superpose Models to Templates')
        for choice in ("Yes", "No"):
            self.superpose_models_to_templates_rds.add(choice)
        self.superpose_models_to_templates_rds.setvalue("Yes")
        if grid_widgets:
            self.superpose_models_to_templates_rds.grid(row=6, **self.options_frame_grid_options)
        else:
            self.superpose_models_to_templates_rds.pack(**pmgi.pack_options_1)
        option_widgets_to_align.append(self.superpose_models_to_templates_rds)

        pmgi.align_set_of_widgets(option_widgets_to_align)


    def show_energy_minimization_frame(self):
        self.energy_minimization_frame.grid(row=4, padx= (25,0), pady= 10, sticky= "w")
        self.energy_minimization_rds.setvalue("Use")
        self.options_scrolled_frame.reposition()

    def hide_energy_minimization_frame(self):
        self.energy_minimization_frame.grid_forget()
        self.energy_minimization_rds.setvalue("None")
        self.options_scrolled_frame.reposition()


    def perform_modelization(self):
        """
        This method is called when the 'SUBMIT' button in the modelization window is pressed. It
        contains the code to instruct Modeller on how to perform the modelization.
        """
        try:
            self._perform_modelization()
        except Exception, e:
            self.modeling_session_failure(e)


    def modeling_session_failure(self, error_message):
        try:
            title = "Modeling Session Error"
            message = "PyMod has encountered the following error while running MODELLER: %s" % error_message
            self.show_error_message(title, message)
            if os.path.isdir(self.model_subdir):
                shutil.rmtree(self.model_subdir)
        except:
            self.show_error_message("Modeling Session Error", "PyMod has encountered an unknown error in the modeling session: %s" % error_message)
        os.chdir(self.current_project_directory_full_path)


    def _perform_modelization(self):

        # -----
        # Takes input supplied by users though the GUI.
        # -----
        # First builds the list of templates (with all their options) the user choosed.
        self.build_templates_list()

        # Starts the modeling process only if the user has supplied correct parameters.
        self.exclude_hetatms = self.exclude_heteroatoms_rds.getvalue()
        self.optimization_level = self.optimization_level_rds.getvalue()
        self.additional_optimization_level = self.energy_minimization_rds.getvalue()
        if self.additional_optimization_level == "Use":
            try:
                nb_cutoff = float(self.energy_minimization_frame.non_bondend_cutoff_rds.getvalue())
            except:
                nb_cutoff = 4.0
            self.additional_optimization_dict = {
                "steepest_descent": {"use": self.energy_minimization_frame.steepest_descent_frame.use_var.get(), "cycles": int(self.energy_minimization_frame.steepest_descent_frame.iterations_enf.getvalue())},
                "conjugate_gradients": {"use": self.energy_minimization_frame.conjugate_gradients_frame.use_var.get(), "cycles": int(self.energy_minimization_frame.conjugate_gradients_frame.iterations_enf.getvalue())},
                "quasi_newton": {"use": self.energy_minimization_frame.quasi_newton_frame.use_var.get(), "cycles": int(self.energy_minimization_frame.quasi_newton_frame.iterations_enf.getvalue())},
                "molecular_dynamics": {"use": self.energy_minimization_frame.md_frame.use_var.get(), "cycles": int(self.energy_minimization_frame.md_frame.iterations_enf.getvalue()), "temperature": int(self.energy_minimization_frame.md_frame.temperature_enf.getvalue())},
                "restraints": {"bond": self.energy_minimization_frame.bonds_checkbutton.getvalue(),
                               "angle": self.energy_minimization_frame.angles_checkbutton.getvalue(),
                               "dihedral": self.energy_minimization_frame.dihedrals_checkbutton.getvalue(),
                               "improper": self.energy_minimization_frame.impropers_checkbutton.getvalue(),
                               "coulomb": self.energy_minimization_frame.lj_checkbutton.getvalue(),
                               "lj": self.energy_minimization_frame.coulomb_checkbutton.getvalue()},
               "non_bonded_cutoff": nb_cutoff}
        self.superpose_to_templates = self.superpose_models_to_templates_rds.getvalue()
        self.color_by_dope_choice = self.color_models_rds.getvalue()

        if not self.check_all_modelization_parameters():
            # "Please Fill all the Fields"
            title = "Input Error"
            message = self.modelization_parameters_error
            self.show_error_message(title, message, self.modeling_window, refresh=False)
            return False

        # The modeling window can be destroyed.
        self.modeling_window.destroy()

        # ---
        # Builds a list with the "knowns" for MODELLER and sets the name of the target sequences.
        # ---
        self.all_templates_namelist = []
        self.modeller_target_name = ""

        # If there is only one chain to model.
        if len(self.modeling_clusters_list) == 1:
            self.all_templates_namelist = self.modeling_clusters_list[0].templates_namelist
            self.modeller_target_name = self.modeling_clusters_list[0].target_name

        # For multiple chains modeling.
        elif len(self.modeling_clusters_list) > 1:
            for mc in self.modeling_clusters_list:
                for t_i,t in enumerate(mc.templates_list):
                    if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                        # Includes the "template complex" name only once.
                        if t.structure.original_pdb_file_name[:-4] not in self.all_templates_namelist:
                            self.all_templates_namelist.append(t.structure.original_pdb_file_name[:-4])
                    else:
                        self.all_templates_namelist.append(mc.templates_namelist[t_i])
            self.modeller_target_name = self.multiple_chains_models_name

        # -----
        # Prepares the directories where MODELLER's output will be generated.
        # -----

        # The absolute path of the models directory.
        models_dir = os.path.join(self.current_project_directory_full_path, self.models_directory)
        # Name of the model subdirectory where Modeller output files are going to be placed.
        model_subdir_name = "%s_%s_%s" % (self.models_subdirectory, self.performed_modeling_count, self.modeller_target_name)
        # The absolute path of the model subdirectory.
        model_subdir = os.path.join(models_dir, model_subdir_name)
        self.model_subdir = model_subdir
        self.create_model_subdirectory(model_subdir)

        # The current directory has to be changed beacause in Modeller the user can't change the
        # output directory, it has to be the current directory.
        os.chdir(model_subdir)

        # -----
        # Start setting options for MODELLER.
        # -----
        #------------------------------------------------------------
        if self.modeller.run_internally():
            modeller.log.verbose()
            env = modeller.environ()
            env.io.atom_files_directory = []
            env.io.atom_files_directory.append(os.path.join(self.current_project_directory_full_path, self.structures_directory))
            env.io.atom_files_directory.append(".")
        #------------------------------------------------------------

        #############################################################
        # If set to True, creates a file "my_model.py" that can be used by command line MODELLER to
        # perform the modellization. It is necessary when using MODELLER as an external coomand line
        # tool, that is, when using MODELLER on PyMOL version which can't import the systemwide
        # 'modeller' library.
        write_modeller_script = True
        if write_modeller_script:
            self.modeller_script = open("my_model.py", "w")
            print >> self.modeller_script, "import modeller"
            print >> self.modeller_script, "import modeller.automodel"
            print >> self.modeller_script, "\n"
            print >> self.modeller_script, "modeller.log.verbose()"
            print >> self.modeller_script, "env = modeller.environ()"
            if not self.modeller.run_internally():
                env = None
            print >> self.modeller_script, "env.io.atom_files_directory = []"
            print >> self.modeller_script, "env.io.atom_files_directory.append(str('%s/%s'))" % (self.current_project_directory_full_path, self.structures_directory)
            print >> self.modeller_script, "env.io.atom_files_directory.append('.')" + "\n"
        #############################################################

        # ---
        # Sets heteroatoms and water options.
        # ---
        use_hetatm = False
        use_water = False
        # If the user wants to include hetero-atoms and water molecules.
        if self.exclude_hetatms == "No":
            # Use hetatm by default in this case.
            #--------------------------------------------------------
            if self.modeller.run_internally():
                env.io.hetatm = True
            #--------------------------------------------------------
            use_hetatm = True
            # Use water only if the user choosed to include water molecules from some template.
            found_water = False
            for mc in self.modeling_clusters_list:
                for (i,t) in enumerate(mc.templates_list):
                    if t.structure.water_state == 1:
                        found_water = True
                        break
            if found_water:
                #----------------------------------------------------
                if self.modeller.run_internally():
                    env.io.water = True
                #----------------------------------------------------
                use_water = True
            else:
                #----------------------------------------------------
                if self.modeller.run_internally():
                    env.io.water = False
                #----------------------------------------------------
                use_water = False

            #########################################################
            if write_modeller_script:
                print >> self.modeller_script, "env.io.hetatm = True"
                if found_water:
                    print >> self.modeller_script, "env.io.water = True"
            #########################################################

        # If the user doesn't want to include hetero-atoms and water.
        else:
            #--------------------------------------------------------
            if self.modeller.run_internally():
                env.io.hetatm = False
                env.io.water = False
            #--------------------------------------------------------

        # ---
        # Creates a file with the alignment in the PIR format.
        # ---
        self.pir_align(os.path.join(model_subdir,"align-multiple.ali"), hetatm=use_hetatm, water=use_water)

        # ---
        # Defines a custom class to use some additional Modeller features.
        # ---
        # This class is going to be used to build the "a" object used to perform the actual
        # homology modelization. It is going to inherit everything from the automodel class
        # but is going to have dynamically redifined routines to make it possible to:
        #   - include user defined disulfide bridges in the model
        #   - exclude template disulfide bridges in the model
        #   - build multichain models with symmetries restraints
        #   - rename the chains in multichain models
        #------------------------------------------------------------
        if self.modeller.run_internally():
            class MyModel(modeller.automodel.automodel):
                pass
        #------------------------------------------------------------

        #############################################################
        if write_modeller_script:
            print >> self.modeller_script, "\n"+"class MyModel(modeller.automodel.automodel):"
        #############################################################

        # ---
        # If there are some targets with at least two CYS residues, it decides whether to use
        # template disulfides or to let the CYS residues in a "reduced" state.
        # ---
        if self.check_targets_with_cys():
            # Decide to use template disulfides or not.
            if True in [mc.has_structures_with_disulfides() for mc in self.modeling_clusters_list]:
                # If the user choosed to use templates disulfides bridges.
                if self.disulfides_frame.use_template_dsb_var.get():
                    # Modeller will automatically use the patch_ss_templates() method of the
                    # automodel class.
                    pass
                # Don't use template dsbs: leave the model CYS residues that in the template are
                # engaged in a disulfied bridge in a "reduced" state.
                else:
                    #------------------------------------------------
                    if self.modeller.run_internally():
                        # This will not create any dsbs in the model by not disulfide patching.
                        def default_patches(self, aln):
                            pass
                        # Dynamically assigns the method.
                        setattr(MyModel, 'default_patches', default_patches)
                    #------------------------------------------------

                    #################################################
                    # Or write it to the modeller script.
                    if write_modeller_script:
                        print >> self.modeller_script, "\n"
                        print >> self.modeller_script, "    def default_patches(self,aln):"
                        print >> self.modeller_script, "        pass"+"\n"
                    #################################################
        # ---
        # Part for multichain models and user defined disulfide bridges, which requires to
        # the special_patches() method override.
        # ---
        if self.check_targets_with_cys():
            self.all_user_defined_dsb = [sel.user_defined_disulfide_bridges for sel in self.disulfides_frame.user_dsb_selector_list]

        #------------------------------------------------------------
        if self.modeller.run_internally():
            def special_patches(self, aln):

                # When building a multichain model it uses the special patches method to rename the
                # chains and give the residues the right ids.
                if len(pymod.modeling_clusters_list) > 1:
                    # Rename the chains. When Modeller builds a multichain model with heteroatoms
                    # and/or water it places them in additional chains. The following code will
                    # rename these extra chains in the right way.
                    segments = [s for s in pymod.target_segment_list if s.use]
                    for chain, segment in zip(self.chains, segments):
                        # print "seq. " + chain.name + " : " + segment
                        chain.name = segment.chain_id

                    # Renumber the residues in the new chains starting from 1. When Modeller builds
                    # a multichain model it doesn't restart to count residues from 1 when changing
                    # chain. The following code renumbers the residues in the correct way.
                    count_dictionary = {}
                    for chain in self.chains:
                        if chain.name not in count_dictionary.keys():
                            count_dictionary.update({chain.name: 1})
                    for chain in self.chains:
                        for num, residue in enumerate(chain.residues):
                            residue.num = '%d' % (count_dictionary[chain.name])
                            count_dictionary[chain.name] += 1

                # Informs Modeller on how to build custom disulfide bridges.
                if True in [mc.target_with_cys for mc in pymod.modeling_clusters_list]:
                    # If the user wants to use some custom dsb.
                    if pymod.disulfides_frame.use_user_defined_dsb_var.get():
                        # Gets the list of user defined dsb for each modeling cluster (if the
                        # target of the modeling cluster doesn't have at least two cys residues
                        # it will have an [] empty list).
                        for (mci,mc) in enumerate(pymod.modeling_clusters_list):
                            # If some user-defined disulfide bridges have been created by the user then get that
                            # information from self.dsb_page.
                            # Populate the self.user_defined_dsb list.
                            for dsb in pymod.all_user_defined_dsb[mci]:
                                # For example CYS321.
                                cys1 = dsb[0][3:]
                                cys2 = dsb[1][3:]
                                # If a bridge has the same cys: <class '_modeller.ModellerError'>: unqang__247E> Internal error:
                                # Redefine the routine to include user defined dsb.
                                if len(pymod.modeling_clusters_list) > 1:
                                    chain = mc.get_template_complex_chain().structure.pdb_chain_id
                                    self.patch(residue_type="DISU", residues=(self.chains[chain].residues[cys1], self.chains[chain].residues[cys2]))
                                else:
                                    self.patch(residue_type="DISU", residues=(self.residues[cys1], self.residues[cys2]))

                    # If the user wants Modeller to build automatically the dsb.
                    if pymod.disulfides_frame.auto_dsb_var.get():
                        # Adds disulfides bridges for cys that are sufficently close.
                        self.patch_ss()

            # Dynamically assigns the method.
            setattr(MyModel, 'special_patches', special_patches)
        #------------------------------------------------------------

        #############################################################
        if write_modeller_script:
            print >> self.modeller_script, "    def special_patches(self, aln):"
            if len(self.modeling_clusters_list) > 1:
                segments = [s for s in self.target_segment_list if s.use]
                print >> self.modeller_script, "        # Rename the chains so that hetatms and water are assigned in the right way."
                print >> self.modeller_script, "        segments = " + repr([s.chain_id for s in segments])
                print >> self.modeller_script, "        for chain, segment in zip(self.chains, segments):"
                print >> self.modeller_script, "            chain.name = segment" + "\n"
                print >> self.modeller_script, "        # Renumber the residues in the new chains starting from 1."
                print >> self.modeller_script, "        count_dictionary = {}"
                print >> self.modeller_script, "        for chain in self.chains:"
                print >> self.modeller_script, "            if chain.name not in count_dictionary.keys():"
                print >> self.modeller_script, "                count_dictionary.update({chain.name: 1})"
                print >> self.modeller_script, "        for chain in self.chains:"
                print >> self.modeller_script, "            for num, residue in enumerate(chain.residues):"
                print >> self.modeller_script, "                residue.num = '%d' % (count_dictionary[chain.name])"
                print >> self.modeller_script, "                count_dictionary[chain.name] += 1" + "\n"

            if self.check_targets_with_cys():
                for (mci,mc) in enumerate(self.modeling_clusters_list):
                    for dsb in self.all_user_defined_dsb[mci]:
                        # For example CYS321.
                        cys1 = dsb[0][3:]
                        cys2 = dsb[1][3:]
                        if len(self.modeling_clusters_list) > 1:
                            chain = mc.get_template_complex_chain().structure.pdb_chain_id
                            print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.chains['%s'].residues['%s'], self.chains['%s'].residues['%s']))" % (chain,cys1,chain,cys2)
                        else:
                            print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.residues['%s'], self.residues['%s']))" % (cys1,cys2)
                if self.disulfides_frame.auto_dsb_var.get():
                    print >> self.modeller_script, "        self.patch_ss()"
        #############################################################

        # ---
        # Apply simmetry restraints to target chains that have the same sequence.
        # ---
        if len(pymod.modeling_clusters_list) > 1:
            groups_to_use = [g for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2) if g.use]
            if len(groups_to_use) > 0:
                # Define group of chains on which symmetry restraints have to be applied.
                list_of_groups = []
                for srg in groups_to_use:
                    list_of_chains = []
                    for mcl in srg.list_of_clusters:
                        if mcl.symmetry_var.get() == 1:
                            list_of_chains.append(mcl.model_chain_id)
                    list_of_groups.append(list_of_chains)

                list_of_symmetry_restraints = []
                for list_of_chains in list_of_groups:
                    s = []
                    for c in range(len(list_of_chains)):
                        i1 = list_of_chains[c]
                        i2 = None
                        if c < len(list_of_chains) - 1:
                            i2 = list_of_chains[c+1]
                        else:
                            pass
                        if i2!=None:
                            s.append([i1,i2])
                    list_of_symmetry_restraints.append(s)

                #----------------------------------------------------
                if self.modeller.run_internally():
                    def special_restraints(self, aln):
                        # Constrain chains to be identical (but only restrain
                        # the C-alpha atoms, to reduce the number of interatomic distances
                        # that need to be calculated):
                        for symmetry_restraints_group in list_of_symmetry_restraints:
                            for s in symmetry_restraints_group:
                                s1 = modeller.selection(self.chains[s[0]]).only_atom_types('CA')
                                s2 = modeller.selection(self.chains[s[1]]).only_atom_types('CA')
                                self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))
                    setattr(MyModel, 'special_restraints', special_restraints)

                    def user_after_single_model(self):
                        # Report on symmetry violations greater than 1A after building
                        # each model:
                        self.restraints.symmetry.report(1.0)
                    setattr(MyModel, 'user_after_single_model', user_after_single_model)
                #----------------------------------------------------

                #####################################################
                if write_modeller_script:
                    print >> self.modeller_script, "    def special_restraints(self, aln):"
                    for si,symmetry_restraints_group in enumerate(list_of_symmetry_restraints):
                         print >> self.modeller_script, "        # Symmetry restraints group n. %d." % (si+1)
                         for s in symmetry_restraints_group:
                             print >> self.modeller_script, "        s1 = modeller.selection(self.chains['" +s[0] + "']).only_atom_types('CA')"
                             print >> self.modeller_script, "        s2 = modeller.selection(self.chains['" +s[1] + "']).only_atom_types('CA')"
                             print >> self.modeller_script, "        self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))"
                    print >> self.modeller_script, "\n"+"    def user_after_single_model(self):"
                    print >> self.modeller_script, "        self.restraints.symmetry.report(1.0)"
                #####################################################

        #############################################################
        if write_modeller_script:
            print >> self.modeller_script, "\n"+"        pass"+"\n"
        #############################################################

        # ---
        # Creates the "a" object to perform the modelization.
        # ---
        #------------------------------------------------------------
        if self.modeller.run_internally():
            a = MyModel(env,
                        alnfile = os.path.join(model_subdir, "align-multiple.ali"),          # alignment filename
                        knowns = tuple([str(tmpn) for tmpn in self.all_templates_namelist]), # codes of the templates
                        sequence = self.modeller_target_name)                                # code of the target
                        #, assess_methods=(modeller.automodel.assess.DOPE))
        #------------------------------------------------------------

        #############################################################
        if write_modeller_script:
            print >> self.modeller_script, "a =  MyModel("
            print >> self.modeller_script, "    env,"
            print >> self.modeller_script, "    alnfile =  'align-multiple.ali',"
            print >> self.modeller_script, "    knowns = " + repr(tuple([str(tmpn) for tmpn in self.all_templates_namelist])) + ","
            print >> self.modeller_script, "    sequence = '%s')" % (self.modeller_target_name)
        #############################################################

        # ---
        # Sets other Modeller options and finally build the model.
        # ---
        if self.optimization_level == "Low":
            #--------------------------------------------------------
            if self.modeller.run_internally():
                # Low VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.very_fast
                # Low MD optimization:
                a.md_level = modeller.automodel.refine.very_fast
            #--------------------------------------------------------

            #########################################################
            if write_modeller_script:
                print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.very_fast"
                print >> self.modeller_script, "a.md_level = modeller.automodel.refine.very_fast"
            #########################################################

        elif self.optimization_level == "Default":
            # a.library_schedule = modeller.automodel.autosched.normal
            # a.max_var_iterations = 200
            # a.md_level = modeller.automodel.refine.very_fast
            # a.repeat_optimization = 2
            # a.max_molpdf = 1e7
            pass

        elif self.optimization_level == "Mid":
            #--------------------------------------------------------
            if self.modeller.run_internally():
                # Thorough VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.fast
                a.max_var_iterations = 300
                # Mid MD optimization:
                a.md_level = modeller.automodel.refine.fast
                # Repeat the whole cycle 2 times.
                a.repeat_optimization = 2
            #--------------------------------------------------------

            #########################################################
            if write_modeller_script:
                print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.fast"
                print >> self.modeller_script, "a.max_var_iterations = 300"
                print >> self.modeller_script, "a.md_level = modeller.automodel.refine.fast"
                print >> self.modeller_script, "a.repeat_optimization = 2"
            #########################################################

        elif self.optimization_level == "High":
            #--------------------------------------------------------
            if self.modeller.run_internally():
                # Very thorough VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.slow
                a.max_var_iterations = 300
                # Thorough MD optimization:
                a.md_level = modeller.automodel.refine.slow
                # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
                a.repeat_optimization = 2
                a.max_molpdf = 1e6
            #--------------------------------------------------------

            #########################################################
            if write_modeller_script:
                print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.slow"
                print >> self.modeller_script, "a.max_var_iterations = 300"
                print >> self.modeller_script, "a.md_level = modeller.automodel.refine.slow"
                print >> self.modeller_script, "a.repeat_optimization = 2"
                print >> self.modeller_script, "a.max_molpdf = 1e6"
            #########################################################

        # ---
        # Determines how many models to build.
        # ---

        #############################################################
        if write_modeller_script:
            print >> self.modeller_script, "a.starting_model= 1"
            print >> self.modeller_script, "a.ending_model = " + str(self.ending_model_number)
            print >> self.modeller_script, "a.make()"
            # Saves an output file that will be red by PyMod when MODELLER is executed externally.
            if not self.modeller.run_internally():
                print >> self.modeller_script, "\n###################################"
                print >> self.modeller_script, "# Needed to run MODELLER externally from PyMOL."
                print >> self.modeller_script, "modeller_outputs_file = open('modeller_saved_outputs.txt','w')"
                print >> self.modeller_script, "modeller_outputs_file.write('[')"
                print >> self.modeller_script, "for model in a.outputs:"
                print >> self.modeller_script, "    model_copy = model.copy()"
                print >> self.modeller_script, "    model_copy.pop('pdfterms')"
                print >> self.modeller_script, "    modeller_outputs_file.write('%s,' % (repr(model_copy)))"
                print >> self.modeller_script, "modeller_outputs_file.write(']')"
                print >> self.modeller_script, "modeller_outputs_file.close()"
            self.modeller_script.close()

        if not self.modeller.run_internally():
            cline=self.modeller.get_exe_file_path() + " my_model.py"
            self.execute_subprocess(cline)
            # Builds the 'a.outputs' when MODELLER was executed externally by reading an output file
            # that was generated in the MODELLER script that was executed externally from PyMOL.
            modeller_outputs_file = open("modeller_saved_outputs.txt","r")
            class Empty_automodel:
                outputs = None
            a = Empty_automodel()
            a.outputs = eval(modeller_outputs_file.readline())
            modeller_outputs_file.close()
        #############################################################

        #------------------------------------------------------------
        if self.modeller.run_internally():
            a.starting_model = 1 # index of the first model
            a.ending_model = int(self.ending_model_number) # index of the last model
            # This is the method that launches the modl building phase.
            a.make()
        #------------------------------------------------------------

        # Changes back the working directory to the project main directory.
        os.chdir(self.current_project_directory_full_path)

        # ---
        # Cycles through all models built by MODELLER to import them into PyMod and PyMOL.
        # ---
        # Perform additional energy minimization.
        if self.additional_optimization_level == "Use":
            list_of_optimized_structures_names = []
            for model in a.outputs:
                optimized_structure_name = self.energy_minimization(model_file_path=os.path.join(model_subdir, model['name']),
                                                         parameters_dict=self.additional_optimization_dict, env=env,
                                                         use_hetatm = self.exclude_hetatms == "No", use_water=found_water)
                list_of_optimized_structures_names.append(optimized_structure_name)
            # Checks if there were some problems in the additional energy minimization process.
            if not None in list_of_optimized_structures_names:
                for model, opt_str_name in zip(a.outputs, list_of_optimized_structures_names):
                    model["name"] = opt_str_name
            else:
                title = "Energy Minimization Error"
                message = "There was an error in the additional energy minimization performed by MODELLER, therefore the final models will not be optimized using the additional energy minimization protocol you selected."
                self.show_error_message(title, message, refresh=True)

        self.models_file_name_dictionary = {}
        model_file_number = 1
        for model in a.outputs:
            ###########################################################################
            # Builds Structure objects for each of the model's chains and loads their #
            # structures in PyMOL.                                                    #
            ###########################################################################
            # Gets the file name generated by MODELLER.
            model_pdb_file_name = model['name']
            model_file_shortcut = os.path.join(model_subdir, model_pdb_file_name)
            model_file_shortcut_in_str_dir = os.path.join(self.structures_directory, model_pdb_file_name)
            # Builds a new file name for the model.
            model_name = str(self.get_model_number()+1)+"_"+self.modeller_target_name
            self.models_file_name_dictionary.update({model['name'] : model_name})
            # Parses tge PDB fiel fo the model.
            model_pdb_file = Parsed_pdb_file(model_file_shortcut)
            model_pdb_file.copy_to_structures_directory()
            model_pdb_file.parse_pdb_file()
            model_pdb_file.build_structure_objects(add_to_pymod_pdb_list = False, new_pdb_file_name=model_name)
            # A list of 'PyMod_element' objects to which the modeling output is going to be
            # assigned. It will be populated with an element for each chain of the model.
            model_chain_elements = []
            for mc in self.modeling_clusters_list:
                # If this is the first model built for this target, then assigns the 'Structure'
                # object to the 'PyMod_element' object of the target sequence.
                if not mc.target.is_model:
                    structure_to_assign = None
                    if len(self.modeling_clusters_list) > 1:
                        structure_to_assign = model_pdb_file.get_chain_structure(mc.model_chain_id)
                    else:
                        structure_to_assign = model_pdb_file.chains_structure_objects[0]
                    mc.target.update_element(new_structure=structure_to_assign)
                    mc.target.is_model = True
                    self.load_element_in_pymol(mc.target)
                    model_chain_elements.append(mc.target)
                    mc.model_elements_list.append(mc.target)
                # If the target sequence has already a model, then insert other models in the
                # 'pymod_elements_list' as new indipendent elements.
                else:
                    element_to_assign = None
                    if len(self.modeling_clusters_list) > 1:
                        element_to_assign = model_pdb_file.get_chain_pymod_element(mc.model_chain_id)
                    else:
                        element_to_assign = model_pdb_file.chains_pymod_elements[0]
                    new_element = element_to_assign
                    self.add_element_to_pymod(new_element, "mother", color="white")
                    self.load_element_in_pymol(new_element)
                    model_chain_elements.append(new_element)
                    mc.model_elements_list.append(new_element)

            ###########################################
            # Superpose models to templates in PyMOL. #
            ###########################################
            if self.superpose_to_templates:
                tc_temp_pymol_name = "template_complex_temp"
                mc_temp_pymol_name = "model_complex_temp"
                # Just superpose the model's chain to the first template.
                if len(self.modeling_clusters_list) == 1:
                    # Superpose in PyMOL the model to its first template.
                    super_template = self.modeling_clusters_list[0].templates_list[0].build_chain_selector_for_pymol()
                    # Builds only a selector for the first and only chain models.
                    model_selector = model_chain_elements[0].build_chain_selector_for_pymol()
                    self.superpose_in_pymol(model_selector, super_template)
                # Superposing is more complex, and follows a different strategy.
                else:
                    # Loads the full template complex file.
                    if model_file_number == 1:
                        template_complex_shortcut = os.path.join(self.structures_directory, self.template_complex.pdb_file_name)
                        cmd.load(template_complex_shortcut, tc_temp_pymol_name)
                        # Superpose each separated chain of the template complex to the corresponding
                        # chains of the full template complex.
                        for mc in self.modeling_clusters_list:
                            chain_id = mc.model_chain_id
                            template_complex_chain = mc.get_template_complex_chain().build_chain_selector_for_pymol()
                            self.superpose_in_pymol(template_complex_chain, "%s and chain %s" % (tc_temp_pymol_name, chain_id), save_superposed_structure=True)
                    # Loads the full model complex file.
                    cmd.load(model_file_shortcut, mc_temp_pymol_name)
                    # Superpose the full model complex file on the template complex using PyMOL.
                    self.superpose_in_pymol(mc_temp_pymol_name,tc_temp_pymol_name, save_superposed_structure=False)
                    # Saves the new superposed file in the structures directory.
                    self.pymol_save(model_file_shortcut_in_str_dir,mc_temp_pymol_name)
                    # Superpose single model chains to the correspondig one of the full model
                    # complex.
                    for me in model_chain_elements:
                        chain_id = me.structure.pdb_chain_id
                        model_chain = me.build_chain_selector_for_pymol()
                        self.superpose_in_pymol(model_chain, "%s and chain %s" % (mc_temp_pymol_name, chain_id), save_superposed_structure=True)
                    # Cleans up.
                    cmd.delete(mc_temp_pymol_name)

            model_file_number += 1
            self.increase_model_number()

        # Finish to clean up.
        if self.superpose_to_templates and len(self.modeling_clusters_list) > 1:
            cmd.delete(tc_temp_pymol_name)

        #####################################
        # Quality assessment of the models. #
        #####################################

        # Starts to build the 'current_modeling_session' which will be used to build a new item on
        # the 'Models' submenu on the main window.
        current_modeling_session = Modeling_session(self.performed_modeling_count + 1)

        # Create DOPE profile for this modeling session (for all models built in this session).
        session_plot_data = []
        # Create DOPE profile for each separated model built in this session.
        for model in a.outputs:
            fmo = Full_model(os.path.join(model_subdir, model['name']))
            single_model_profile = []
            fmo.model_profile = single_model_profile
            current_modeling_session.full_models.append(fmo)

        # Computes the DOPE scores.
        alignment_lenght = 0
        assessed_structures_list = []
        # This for cycle is used to add extra 'None' values in multiple chains profiles. In this way
        # if, for example, there is model with chains 'A' and 'B', in the matplotlib plot the
        # profile of chain 'B' will be put just after the end of the profile of chain 'A'.
        for mc in sorted(self.modeling_clusters_list, key = lambda mc: mc.block_index):
            mc.adjust_model_elements_sequence()
            # Actually computes the DOPE profile of the templates.
            for template in mc.templates_list:
                self.compute_dope(template, env=env)
                assessed_structures_list.append(template) # template_dope_data
                template_dope_data = self.prepare_dope_plot_data([template], start_from=alignment_lenght, mode="multiple")
                # Stores the templates also profiles to each 'Full_model' object, so that the
                # profile of the templates can be inspected by accessing the 'Models' menu.
                for fmo in current_modeling_session.full_models:
                    fmo.model_profile.append(template_dope_data[0])
                session_plot_data.append(template_dope_data[0])
            # Computes the DOPE profile of the models.
            for model_element, fmo in zip(mc.model_elements_list, current_modeling_session.full_models):
                self.compute_dope(model_element, env=env)
                model_dope_data = self.prepare_dope_plot_data([model_element], start_from=alignment_lenght, mode="multiple")
                assessed_structures_list.append(model_element)
                # Stores the templates profiles to each 'Full_model' object, so that the profile of
                # the models can be accessed in the 'Models' menu.
                fmo.model_profile.append(model_dope_data[0])
                session_plot_data.append(model_dope_data[0])
            alignment_lenght += len(mc.target.my_sequence)

        # Gets the objective function and DOPE scores values for each full model (the model
        # comprising all the chains) built.
        assessment_data = []
        list_of_models_names = []
        column_headers = ["Objective Function Value", "DOPE score"]
        for model, fmo in zip(a.outputs, current_modeling_session.full_models):
            # Gets the Objective function values.
            model_pdb_file_name = model['name']
            list_of_models_names.append(model_pdb_file_name)
            model_file_shortcut = os.path.join(model_subdir, model_pdb_file_name)
            model_file = open(model_file_shortcut, "r")
            obj_funct_value = float(model_file.readlines()[1][39:].replace(" ",""))
            model_file.close()
            # Gets the DOPE values.
            model_profile_shortcut = os.path.join(model_subdir, model_pdb_file_name[:-4]+".profile")
            dope_score = self.compute_dope_of_structure_file(model_file_shortcut, model_profile_shortcut,env=env)
            obj_funct_value, dope_score = round(obj_funct_value, 3), round(dope_score, 3)
            assessment_data.append([obj_funct_value, dope_score])
            fmo.assessment_data = [obj_funct_value, dope_score]

        # Prepares data to show a table with objective function values and DOPE scores for each
        # model.
        assessment_table_args = {"column_headers": column_headers, "row_headers": list_of_models_names, "data_array": assessment_data, "title": "Assessment of Models", "number_of_tabs": 4, "width": 850, "height" :420, "rowheader_width": 25}
        current_modeling_session.assessment_table_data = assessment_table_args
        current_modeling_session.session_profile = session_plot_data
        self.modeling_session_list.append(current_modeling_session)
        self.performed_modeling_count += 1

        self.assign_dope_items(assessed_structures_list)

        # Colors the models and templates according to their DOPE values. This follows the same
        # method used in the 'dope_from_main_menu()' method.
        if self.color_by_dope_choice == "DOPE Score":
            for element in assessed_structures_list:
                element.color_element_by_dope()

        self.gridder()

        # Finally shows the table and the previously built DOPE profile comprising DOPE curves
        # of every model and templates.
        self.show_table(**assessment_table_args)
        self.show_dope_plot(session_plot_data)


    #################################################################
    # Prepares input for MODELLER.                                  #
    #################################################################
    def build_templates_list(self):
        """
        Add to the Modeling_clusters objects information about which templates to use according to
        the parameters supplied by users.
        """
        for (mc_index, modeling_cluster) in enumerate(self.modeling_clusters_list):
            # Gets the values of each template checkbutton (it will be 0 if the structure was
            # not selected or 1 if it was selected).
            pdb_chains_map = map(lambda var: var.get(), modeling_cluster.get_use_as_template_states())

            # Begins a for cycle that is going to get the structures to be used as templates.
            modeling_cluster.templates_list = []
            modeling_cluster.templates_namelist = []
            for (a, structure_frame) in enumerate(modeling_cluster.structure_frame_list):

                # Selects only structures that were selected by the user to be used as templates.
                if pdb_chains_map[a] == 1:

                    # For every selected structure takes the HETRES option.
                    single_structure_hetres_option_value = structure_frame.hetres_options_var.get()
                    # And the values of each HETRES checkbutton.
                    single_structure_hetres_map = map(lambda var_e: var_e.get(), structure_frame.structure_hetres_states)
                    # Do the same with the water checkbutton.
                    water_state =  structure_frame.water_state.get()

                    # Adds some information about the modeling options to the elements.
                    modeling_cluster.structure_list[a].structure.set_modeling_information(
                        seq_min = 1, # structure_frame.from_enf.getvalue(),
                        seq_max = 10000, # structure_frame.to_enf.getvalue(),
                        hetres_option = single_structure_hetres_option_value,
                        hetres_map = single_structure_hetres_map,
                        water_state = water_state)

                    # Populate each modeling_cluster "template_list" with the elements selected by
                    # the user from the "structure_list".
                    modeling_cluster.templates_list.append(modeling_cluster.structure_list[a])
                    modeling_cluster.templates_namelist.append(modeling_cluster.structure_list[a].structure.chain_pdb_file_name.replace(":","_")[:-4])

                    # IF THE ORIGINAL PDB FILES ARE TO BE USED:
                    #     - self.struct_list[a].structure.original_chain_pdb_file_name.replace(":","_")
                    # In the original Pymod it was:
                    #     - "1UBI_Chain_A" for non ce-aligned seqs
                    #     - "1UBI_Chain_A_aligned.pdb" for aligned seqs
                    # These codes must be the same in the .ali file and when assigning the "knowns".
                    # If it the names don't contain the .pdb extension, Modeller will still find the
                    # right files.


    def check_all_modelization_parameters(self):
        """
        This will be used before launching Modeller to check:
            - if the parameters of each modeling clusters are correct
            - when performing multichain modeling
                - if there is exactly 1 template complex chain selected in each cluster
                - if symmetry restraints buttons are selected properly
        """
        # Checks if the parameters of all the "modeling clusters" are correct.
        for mc in self.modeling_clusters_list:
            if not self.check_modeling_cluster_parameters(mc):
                return False

        # Checks if there are only correct sequences.
        for mc in self.modeling_clusters_list:
            if not pmsm.check_correct_sequence(mc.target.my_sequence):
                self.modelization_parameters_error = "Target sequence '%s' contains an invalid character in its sequence (%s) and MODELLER can't modelize it." % (mc.target.my_header, pmsm.get_invalid_characters_list(mc.target.my_sequence)[0])
                return False
            for t in mc.templates_list:
                if not pmsm.check_correct_sequence(t.my_sequence):
                    self.modelization_parameters_error = "Template '%s' contains an invalid character in its sequence (%s) and MODELLER can't use it as a template." % (t.my_header, pmsm.get_invalid_characters_list(t.my_sequence)[0])
                    return False

        # If each "modeling cluster" has correct parameters, when performing multiple chain modeling,
        # there are other conditions that must be satisfied.
        if len(self.modeling_clusters_list) > 1:

            # First finds the PDB_file object of the "template complex" selected by the user.
            self.template_complex = None
            for p in self.pdb_list:
                if p.pdb_file_name == self.template_complex_var.get():
                    self.template_complex = p

            # Then perform additional controls for each modeling cluster and also get the list of
            # the "target complex" chains selected by the user.
            self.template_complex_selected_chain_list = []
            for mc in self.modeling_clusters_list:

                # Gets the "template complex" chains selected in the current modeling cluster.
                template_complex_selected_chains_in_cluster = []
                for t in mc.templates_list:
                    if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                        template_complex_selected_chains_in_cluster.append(t.structure.pdb_chain_id)
                self.template_complex_selected_chain_list.extend(template_complex_selected_chains_in_cluster)

                # Check if the current cluster has a selected chain from the "target complex".
                if len(template_complex_selected_chains_in_cluster) == 0:
                    self.modelization_parameters_error = "Please select AT LEAST one chain from the 'Template Complex' (%s) as a template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
                    return False

                # Checks if in some cluster there is more than one selected template belonging to the
                # "template complex". This is needed for because ONLY one chain belonging to the
                # "template complex" can be selected by ther user in each cluster.
                if len(template_complex_selected_chains_in_cluster) > 1:
                    self.modelization_parameters_error = "Please select ONLY one chain from the 'Template Complex' (%s) as template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
                    return False

            # Finally checks if the symmetries checkbuttons are selected properly.
            if not self.check_symmetry_vars():
                return False

        # Checks if a correct value in the max models entry has been supplied.
        if not self.check_max_model_entry_input():
            return False

        if self.additional_optimization_level == "Use":
            if not self.check_energy_minimization_parameters():
                return False

        # Returns 'True' only if all parameters are correct.
        return True


    def check_max_model_entry_input(self):
        """
        Checks the "max_models_entry" input and gets its value.
        """
        self.ending_model_number = self.max_models_enf.getvalue()
        if self.ending_model_number == "":
            self.modelization_parameters_error = "Non valid input in the 'Models to calculate' entry!"
            return False
        else:
            return True


    def check_energy_minimization_parameters(self):
        if False in [frame.check_parameters() for frame in self.energy_minimization_frame.minimization_algorithms_frames]:
            self.modelization_parameters_error = "Invalid parameters in the Additional Energy Minimization Options!"
            return False
        if not 1 in [frame.use_var.get() for frame in self.energy_minimization_frame.minimization_algorithms_frames]:
            self.modelization_parameters_error = "Please select at least one Additional Energy Minimization algorithm!"
            return False
        if not 1 in [checkbutton.getvalue() for checkbutton in self.energy_minimization_frame.list_of_parameters_checkbuttons]:
            self.modelization_parameters_error = "Please select at least one feature to minimize!"
            return False
        if 1 in [self.energy_minimization_frame.lj_checkbutton.getvalue(), self.energy_minimization_frame.coulomb_checkbutton.getvalue()]:
            if self.energy_minimization_frame.non_bondend_cutoff_rds.getvalue() == "":
                self.modelization_parameters_error = "Please insert a non bonded cutoff value!"
                return False
        return True


    def check_modeling_cluster_parameters(self, modeling_cluster):
        """
        Checks the if there are any problems with the user-supplied parameters of a "modeling cluster"
        before starting the modeling process.
        """
        self.modelization_parameters_error = ""
        # Checks if there are some templates that have been selected.
        if modeling_cluster.templates_list == []:
            self.modelization_parameters_error = "You have to select at least one template for target '%s' in order to build a model!" % (modeling_cluster.target_name)
            return False
        if not self.check_templates_limits_input(modeling_cluster):
            return False
        return True


    def check_templates_limits_input(self,modeling_cluster):
        """
        Checks the sequence limits entries. It will only return 'True' if the input provided by the
        user is correct.
        """
        for template in modeling_cluster.templates_list:
            if template.structure.seq_min == "" or template.structure.seq_max == "":
                self.modelization_parameters_error = "Non valid input in the 'From - to' entries of template %s!" % (template.my_header)
                return False
            template.structure.seq_min = int(template.structure.seq_min)
            template.structure.seq_max = int(template.structure.seq_max)
            if template.structure.seq_max < template.structure.seq_min:
                self.modelization_parameters_error = "The upper sequence limit (%s) can't be greater than the lower one (%s) for template %s!" % (template.structure.seq_max, template.structure.seq_min, template.my_header)
                return False
            if template.structure.seq_max == template.structure.seq_min:
                self.modelization_parameters_error = "The upper and lower sequence limits of template %s can't be equal!" % (template.my_header)
                return False
        return True


    def check_symmetry_vars(self):
        correct_symmetry_vars = True
        for srg in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2):
            si = len([mc for mc in srg.list_of_clusters if mc.symmetry_var.get() == 1])
            if si == 1:
                correct_symmetry_vars = False
                self.modelization_parameters_error = "In order to impose symmetry restraints you need select the 'Apply symmetry restraints' option for at least two targets with the same sequence (you selected this option only for target '%s')." % (mc.target_name)
                break
            elif si > 1:
                srg.use = True
            else:
                srg.use = False
        return correct_symmetry_vars


    # ---
    # This function creates alignments in a PIR format: this is entirely rewrtitten from the
    # original PyMod version.
    # The custom sequence object should already contain all the info created inside this method.
    # NOTE! The really useful rjust() or ljust() string methods could be used here for making
    # the code much more compact!
    # ---
    def pir_align(self, alignment_file_name, hetatm=False, water=False):

        fn=open(alignment_file_name, "w")

        for (mc_i,mc) in enumerate(self.modeling_clusters_list):

            # Finds the longest sequence among the templates and the target.
            # NOTE: Remove this, there should not be a need for it.
            maximum_length = max([len(t.my_sequence) for t in mc.templates_list[:]+[mc.target]])
            maximum_water_molecules_number = max([t.structure.water_molecules_count for t in mc.templates_list])

            # This is going to contain all the template sequences, because they are going to be used in
            # order to insert "modified residues" in the target sequence.
            complete_template_list = []
            # This list is needed for the same purpose.
            all_modified_residues_positions = []

            # Finds the template selected to include water molecules (there can be at most one per
            # modeling cluster).
            if water:
                for (i,t) in enumerate(mc.templates_list):
                    if t.structure.water_state == 1:
                        mc.set_water_molecules_number(t.structure.water_molecules_count)
                        mc.use_water_in_cluster = True
                        if (len(self.modeling_clusters_list) > 1 and
                            t.structure.original_pdb_file_name == self.template_complex.pdb_file_name):
                            mc.use_template_complex_waters = True
                        else:
                            mc.use_template_complex_waters = False
                        break

            # ---
            # Starts by building the template sequences.
            # ---
            for (i,template) in enumerate(mc.templates_list):
                # Adjust hetres options according to the user's preferences.
                if template.structure.hetres_option != 2:
                    for (k,h) in enumerate(template.structure.hetres_map):
                        # If the user selected the "Use all heteroatomic residues" use all hetres for this
                        # template.
                        if template.structure.hetres_option == 1:
                            template.structure.hetres_map[k] = 1
                        # Don't use any hetres if the user selected "Do not use any heteroatomic residue".
                        elif template.structure.hetres_option == 3:
                            template.structure.hetres_map[k] = 0

                # The sequence to be printed on the alignment file.
                template_sequence = str(template.my_sequence)

                # Not really necessary.
                if len(template_sequence) < maximum_length:
                    template_sequence += "-"*(maximum_length-len(template_sequence))

                # ---
                # Part for the modified residues.
                # ---

                # Modified residues position in the alignment.
                modified_residues_positions = []

                # Converts the sequence in a list to change the correspodding residues in "."s or "-"s.
                template_sequence_list = list(template_sequence)
                # Mark modified residues as "."
                for (k,res) in enumerate(template.structure.pdb_chain_sequence):
                    if res.residue_type == "het" and res.hetres_type == "modified-residue":
                        het_position_in_alignment = pmsm.get_residue_id_in_aligned_sequence(template_sequence,res.id)
                        modified_residues_positions.append(het_position_in_alignment)
                        # If the user want to use HETRES ("env.io.hetatm = True"), includes "."
                        # characters for modified residues.
                        if hetatm == True:
                            template_sequence_list[het_position_in_alignment] = "."
                        # If the user doens't want to use HETRES ("env.io.hetatm = False"), marks
                        # modified residues as "-". Modeller wants modified residues to be marked as
                        # "-" when "env.io.hetatm = False".
                        else:
                            template_sequence_list[het_position_in_alignment] = "-"

                # Reconverts the list in a string.
                template_sequence = "".join(template_sequence_list)
                all_modified_residues_positions.append(modified_residues_positions)

                # ---
                # Part for the sequence limits residues.
                # ---
                # Adjust the template sequences according to the limits supplied by the user through the
                # "From - To" entries.
                template_sequence = list(template_sequence)
                for (k,position) in enumerate(template_sequence):
                    if k < template.structure.seq_min - 1 or k > template.structure.seq_max - 1:
                        # Converts residues out of the interval chosen by the user in "-" characters
                        # so that Modeller doesn't see them.
                        template_sequence[k] = "-"
                template_sequence = "".join(template_sequence)

                # ---
                # Part for ligand hetres and water molecules.
                # ---
                template_ligands_sequence = ""
                if hetatm == True:
                   for (j,t) in enumerate(mc.templates_list):
                        # These lists are going to contain the checkbox states of HETRES in order to
                        # include at the end of the alignment the right number and combination of "."
                        # and "-".
                        # This is going to contain the checkbox of "ligands" HETRES
                        ligands = []
                        # Right now t.hetres_map contains the HETRES checkbox states, for example
                        # something like: [0, 0, 1, 0].
                        for (k,h) in enumerate(t.structure.hetero_residues):
                            if (h.hetres_type == "ligand"):
                                ligands.append(t.structure.hetres_map[k])
                        for li,h in enumerate(ligands):
                            # If the template is the right one.
                            if j == i:
                                template_ligands_sequence += "."
                            # If it is not the right one adds as many gaps for each HETRES found
                            # in the other templates.
                            else:
                                template_ligands_sequence += "-"

                # Include water molecules in the alignment. Each water molecule count as a residue and
                # is indicated with a "w".
                template_water_sequence = ""
                if water:
                    if len(self.modeling_clusters_list) > 1 and template.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                            template_water_sequence += "w"*template.structure.water_molecules_count
                    else:
                        if mc.use_water_in_cluster:
                            if template.structure.water_state == 1:
                                template_water_sequence += "w"*template.structure.water_molecules_count
                            else:
                                template_water_sequence += "-"*mc.water_molecules_number

                complete_template_list.append(template_sequence)
                # Sets the sequence
                pir_template = PIR_alignment_sequence(template_sequence,template_ligands_sequence,template_water_sequence)
                template.set_pir_alignment_sequence(pir_template)

            # ---
            # Target.
            # ---
            target_sequence = mc.target.my_sequence
            # Adjust to the longest sequence.
            if len(target_sequence) < maximum_length:
                target_sequence += "-"*(maximum_length-len(target_sequence))
            target_ligands_sequence = ""
            target_water_sequence = ""
            # Adds HETRES.
            if hetatm == True:
                # Adds modified residues if they were selected.
                # Make it better, it creates strange residues if more than one template in the alignment
                # has in the same position a modified residue...
                for (i,template_sequence) in enumerate(complete_template_list):
                    # [[hetres],[state]] it has as many elements as hetres in that template.
                    single_template_hetres = []
                    # Gets from the checkbox states maps only the ??? belonging to modified residues.
                    for (k,h) in enumerate(mc.templates_list[i].structure.hetres_map):
                        if mc.templates_list[i].structure.hetero_residues[k].hetres_type == "modified-residue":
                            single_template_hetres.append([mc.templates_list[i].structure.hetero_residues[k] , mc.templates_list[i].structure.hetres_map[k]])
                    # Sets a "-" or "." character in the target sequence for each modified residue.
                    for (k,t) in enumerate(single_template_hetres):
                        # Converts the sequence in a list to change the correspodding residues.
                        target_sequence = list(target_sequence)
                        for mr in all_modified_residues_positions[i]:
                             # If the checkbox was selected.
                             if t[1] == 1:
                                 for (x,residue) in enumerate(target_sequence):
                                     if residue != "-" or residue != ".":
                                         target_sequence[mr] = "."
                        # Reconverts the list in a string.
                        target_sequence = "".join(target_sequence)

                # Adds ligands if they were selected by the user.
                for t in mc.templates_list:
                    for (i,h) in enumerate(t.structure.hetero_residues):
                        if h.hetres_type == "ligand":
                            if t.structure.hetres_map[i] == 1:
                                target_ligands_sequence +="."
                            else:
                                target_ligands_sequence +="-"

            # Includes water molecules if some template was selected to include them.
            if water and mc.use_water_in_cluster:
                target_water_sequence += "w"*mc.water_molecules_number

            pir_target = PIR_alignment_sequence(target_sequence,target_ligands_sequence,target_water_sequence)
            mc.target.set_pir_alignment_sequence(pir_target)

        # template_block = reduce(lambda x,y: x+"/"+y, complete_template_list)

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        # ---
        # Single chain models.
        # ---
        if len(self.modeling_clusters_list) == 1:

            # First write to the file the template blocks.
            mc = self.modeling_clusters_list[0]
            for (i,template) in enumerate(mc.templates_list):

                # Writes the first line of the template. Modeller does not like names with ":" character
                # but they have been removed when populating self.templates_namelist.
                template_code = mc.templates_namelist[i]
                template_chain = template.structure.pdb_chain_id
                print >> fn , ">P1;"+template_code
                print >> fn , "structure:%s:.:%s:.:%s::::" % (template_code,template_chain,template_chain)
                # Print one the alignment file 60 characters-long lines.
                template_sequence = self.get_pir_formatted_sequence(template.pir_alignment_sequence.get_single_chain(use_hetres = hetatm, use_water = mc.use_water_in_cluster))
                print >> fn, template_sequence

            # There is still a problem with multiple templates and modified residues.
            # Then writes the target block.
            print >> fn , ">P1;"+self.modeller_target_name# mc.target_name
            print >> fn , "sequence:"+self.modeller_target_name+":.:.:.:.::::"
            target_sequence = self.get_pir_formatted_sequence(mc.target.pir_alignment_sequence.get_single_chain(use_hetres = hetatm, use_water = mc.use_water_in_cluster))
            print >> fn, target_sequence
            mc.set_block_index(0)

        # ---
        # Mulitple chains models.
        # ---
        elif len(self.modeling_clusters_list) > 1:

            # First find the right order of the modeling clusters.
            modeling_cluster_new_index = 0
            for (c_i,chain) in enumerate(self.template_complex.chains_list):
                for mc in self.modeling_clusters_list:
                    for t in mc.templates_list:
                        if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                            if t.structure.pdb_chain_id == chain:
                                mc.set_block_index(c_i) # modeling_cluster_new_index)
                                mc.set_model_chain_id(t.structure.pdb_chain_id)
                                # template_complex_chains.append(chain)
                                modeling_cluster_new_index += 1

            # ---
            # Builds the block structure.
            # ---
            segment_structure_to_print = []
            for segment in self.template_complex.segment_structure:
                if segment[1] == "A":
                    segment_structure_to_print.append(segment)
                elif segment[1] == "H" and hetatm:
                    segment_structure_to_print.append(segment)
                elif segment[1] == "W" and water:
                    segment_structure_to_print.append(segment)

            ogb = Original_Block()
            # Template complex segments.
            for segment in segment_structure_to_print:
                chain = self.template_complex.get_chain_by_id(segment[0])
                modeling_cluster = None
                for mc in self.modeling_clusters_list:
                    if chain in mc.templates_list:
                        modeling_cluster = mc
                        break
                sg = Segment(segment,modeling_cluster)
                ogb.add_segment(sg)
            # Extra hetatms segments.
            for mc in sorted(self.modeling_clusters_list, key = lambda mc:mc.block_index):
                if mc.has_ligands() and not mc.template_complex_chain_has_ligands():
                    sg = Segment((None,"H"),mc)
                    ogb.add_segment(sg)
            # Extra water segments.
            for mc in sorted(self.modeling_clusters_list, key = lambda mc:mc.block_index):
                if mc.use_water_in_cluster and not mc.use_template_complex_waters:
                    sg = Segment((None,"W"),mc)
                    ogb.add_segment(sg)

            # Now generates the blocks.
            list_of_blocks = []
            # Template complex block.
            list_of_blocks.append(ogb.generate_template_complex_block(self.template_complex))

            # Other templates.
            for mc in self.modeling_clusters_list:
                for t_i, t in enumerate(mc.templates_list):
                    if t.structure.original_pdb_file_name != self.template_complex.pdb_file_name:
                        list_of_blocks.append(ogb.generate_additional_template_block(t))
            # Target.
            list_of_blocks.append(ogb.generate_target_block())
            self.target_segment_list = ogb.get_target_segment_list()

            # Prints the whole alignment file.
            for bl in list_of_blocks:
                print >> fn, bl.first_line
                print >> fn, bl.second_line
                print >> fn, self.get_pir_formatted_sequence(reduce(lambda s1,s2: s1+"/"+s2,bl.segment_list),multi=True)

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        fn.close()


    def get_pir_formatted_sequence(self,sequence,multi=False):
        formatted_sequence = ""
        for s in xrange(0,len(sequence),60):
            # For all the lines except the last one.
            if (len(sequence) - s) > 60:
                formatted_sequence += sequence[s:s+60] + "\n"
            # For the last line.
            else:
                if not multi:
                    formatted_sequence += sequence[s:]+"*"+"\n"
                else:
                    formatted_sequence += sequence[s:]+"/*"+"\n"
        return formatted_sequence


    def get_model_number(self):
        model_number = 0
        if len(self.modeling_clusters_list) > 1:
            model_number = self.multiple_chain_models_count
        else:
            model_number = self.modeling_clusters_list[0].target.models_count
        return model_number

    def increase_model_number(self):
        if len(self.modeling_clusters_list) > 1:
            self.multiple_chain_models_count += 1
        else:
            self.modeling_clusters_list[0].target.models_count += 1


#####################################################################
# Classes for the graphical user interface.                         #
#####################################################################

class Cluster_selection_frame(Frame):
    """
    Class used to build a frame containing the widgets necessary to select a cluster from a
    combobox. This is used in the alignment options window.
    """

    def __init__(self, parent_widget, involved_cluster_elements_list, label_text):
        # Frame with the options to control the new alignment.
        Frame.__init__(self,parent_widget, background='black', pady=5, padx=5, bd=0, relief='groove')
        # Builds the lists to be displayed in the comboboxes.
        self.involved_clusters_combobox_list = [e.my_header for e in involved_cluster_elements_list]
        # Label.
        self.target_alignment_label = Label(self, fg="white" , text= label_text, background='black', padx = 20)
        self.target_alignment_label.grid(row=0, column=0, sticky = "w")

        # Combobox.
        self.target_alignment_combobox = Pmw.ComboBox(
                        self,
                        labelmargin = None, labelpos = None,
                        scrolledlist_items = self.involved_clusters_combobox_list,
                        history = 0 )
        # Make the combobox entries not editable.
        self.target_alignment_combobox.component("entryfield").component("entry").configure(
                        state='readonly', readonlybackground= "white", width=30,
                        fg="black", bg="white")
        self.target_alignment_combobox.grid(row = 0,column = 1)
        self.target_alignment_combobox.selectitem(0) # Selects the first cluster of the list.

    def get_selected_cluster(self):
        """
        Gets the name of the cluster selected in the combobox.
        """
        return self.target_alignment_combobox.get()

    def get_selected_cluster_index(self, cluster_name):
        return self.involved_clusters_combobox_list.index(cluster_name)


class Parsed_pdb_file:
    """
    Class used to parse PDB files and to build from them 'PyMod_element' and 'Structure' objects
    that are going to be used in PyMod.
    """

    def __init__(self,original_file_full_path, new_pdb_file_name=None):
        """
        Just prepare attributes that will be filled when parsing the PDB file.
        """
        # Absolute path of the original PDB file.
        self.original_file_full_path = original_file_full_path

        # Name of PDB file that will be copied in the Structures directory from the original one.
        if new_pdb_file_name == None:
            # If a PDB is named "mypdbfile.pdb", this will be set to "mypdbfile".
            self.pdb_file_base_name = os.path.splitext(os.path.basename(self.original_file_full_path))[0]
            self.copied_pdb_file_name = pymod.build_header_string(self.pdb_file_base_name)+".pdb"
        else:
            # The base name takes the value of the 'new_pdb_file_name' argument.
            self.pdb_file_base_name = new_pdb_file_name
            self.copied_pdb_file_name = pymod.build_header_string(new_pdb_file_name)+".pdb"

        # ---
        # Data taken by parsing the PDB file without biopython.
        # ---
        # This is going to contain information taken from the HETRES lines in the PDB file.
        # It will contain dictionary items like: {"name":"GLC", "position":0,"chain":"A"}
        self.structure_hetres = []
        # This is going to contain info from the MODRES section of the PDB file. It will also contain
        # dictionary items.
        self.structure_modres = []

        # This is going to contain the chains upper limits, taken from the SEQRES part of the PDB
        # file. It will be needed for choosing if an HETRES is a modified residue or a ligand.
        # It will contain dictionary object in this form: {"chain":"A","upper_limit": 502}
        self.structure_chains_upper_limits = []

        # This will contain info about a structure disulfides, taken from the SSBOND section.
        # Its elements are going to be "Disulfide_bridge" objects.
        self.structure_disulfides = []

        # The PDB code taken from the HEADER line, if present.
        self.pdb_code = None

        # This is a list of elements like. It will be used to build the 'segment_structure', which
        # is going to be used to build .pir alignment files for multichain modeling.
        self.list_of_residue_ids = []
        # This is needed in order to represent the order in which specific kind of atoms appear in
        # the PDB file. It will contain a list of tuples with the following information:
        # ('chain_id', 'kind_of_atoms'). To make an example, the file with PDB code 1T5A (which
        # contains four polypeptidic chains: 'A', 'B', 'C' and 'D') will have the following
        # 'segment_structure':
        # [('A', 'A'), ('B', 'A'), ('C', 'A'), ('D', 'A'), # 'A' (main chain atoms) of the four chains come first.
        #  ('A', 'H'), ('B', 'H'), ('C', 'H'), ('D', 'H'), # 'H' (hetero-atoms) of the four chains come after.
        #  ('A', 'W'), ('B', 'W'), ('C', 'W'), ('D', 'W')] # 'W' (water atoms) of the chains come last.
        # This is actually the standard order in which different kind of atoms appear in PDB files,
        # but sometimes users might use custom PDB files with a non-standard order, and this
        # attribute is needed to build a .pir alignment file with different types of residues in the
        # same order of the PDB file.
        self.segment_structure = []

        # ---
        # Data taken by parsing the PDB file through biopython.
        # ---
        self.parsed_pdb_file_chains = []
        self.parsed_biopython_structure = None


    def copy_to_structures_directory(self, overwrite_in_str_dir=True):
        """
        Copies the original PDB file to the structure directory. This is needed to build multiple
        chain models, because the original PDB file retains information about the quaternary
        structure.
        """
        pdb_file_shortcut = os.path.join(pymod.structures_directory, self.copied_pdb_file_name)
        # Overwrite the file if it already exists in the project folder.
        if os.path.isfile(pdb_file_shortcut):
            os.remove(pdb_file_shortcut)
        shutil.copy(self.original_file_full_path, pdb_file_shortcut)


    def parse_pdb_file(self):
        self.parse_pdb_file_custom()
        self.parse_pdb_file_using_biopython()


    def parse_pdb_file_custom(self):
        """
        Parses the PDB file for some info that biopython doesn't know how to get.
        """
        fh = open(self.original_file_full_path,"r")
        pdb_file_content = fh.readlines()

        self.found_header_line = False
        self.found_het_lines = False
        self.found_modres_lines = False
        self.found_ssbond_lines = False
        self.found_seqres_lines = False

        for (i,line) in enumerate(pdb_file_content):

            # First line. Assigns the pdb_code.
            if i == 0:
                rline = line.replace(" ","").rstrip("\r\n")
                # Only 'HEADER' line defines PDB id.
                if line.startswith("HEADER"):
                    self.pdb_code = rline[-4:]

            # Finds information about hetero-residues.
            # Finds the HET lines.
            # HET    NZS  A   1      27
            if line[0:4] == "HET ":
                hetres = {
                   "name" : line[7:10],
                   "position" : int(line[13:17]) ,
                   "chain" : line[12]}
                self.structure_hetres.append(hetres)
                if not self.found_het_lines:
                    self.found_het_lines = True
            # Finds the MODRES.
            # 0           12  16 18   23
            # MODRES 1SPU PAQ A  466  TYR
            if line[0:6] == "MODRES":
                modres = {
                   "name" : line[12:15],
                   "position" : int(line[18:22]) ,
                   "chain" : line[16],
                   "original-residue": line[23:27]}
                self.structure_modres.append(modres)
                if not self.found_modres_lines:
                    self.found_modres_lines = True

            # Finds information about disulfides bridges.
            if line[0:6] == "SSBOND":
                # 0        9     15 18         29 32                                        74
                # SSBOND   1 CYS A  320    CYS A  404                          1555   1555  2.03
                # Part for the first residue involved in the bond.
                dsb = Disulfide_bridge(
                    cys1_pdb_number = int(line[17:21]),
                    cys1_chain = line[15],
                    cys2_pdb_number = int(line[31:35]),
                    cys2_chain = line[29],
                    distance = 0,
                    number_in_pdb = int(line[8:10]) )
                self.structure_disulfides.append(dsb)

            # Finds information about the sequences in the SEQRES lines.
            # Finds the sequence maximum limit. It tells which hetres is a modified residues and
            # which not.
            if line[0:6] == "SEQRES":
                chain = line[11]
                # This is needed to represent each chain only once in structure_chains_upper_limits:
                # SEQRES for a single chain often spans more than one line ine the PDB file.
                if not chain in [e["chain"] for e in self.structure_chains_upper_limits]:
                    self.structure_chains_upper_limits.append({"chain":chain,"upper_limit": int(line[12:17])})

            # Get the "residue id". This is needed to build multichain models.
            resid = None
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]
                if line.startswith("ATOM"):
                    residue_type = "A"
                    resid = (chain,residue_type)
                elif line.startswith("HETATM"):
                    if line[17:20] == "HOH":
                        residue_type = "W"
                    else:
                        residue_type = "H"
                    resid = (chain,residue_type)
            if resid != None:
                self.list_of_residue_ids.append(resid)

        fh.close()

        # Builds the segment structure, needed to build .pir alignment files for multichain
        # models.
        for i,res in enumerate(self.list_of_residue_ids):
            if not res in self.segment_structure:
                if res[1] != "H":
                    self.segment_structure.append(res)
                else:
                    modres = False
                    for j,r in enumerate(self.list_of_residue_ids[i:]):
                        if r[0] == res[0] and r[1] == "A":
                            modres = True
                            break
                    if not modres:
                        self.segment_structure.append(res)


    def parse_pdb_file_using_biopython(self):
        """
        Parses the PDB file using biopython to retrieve the sequence of each chain in the file.
        """
        warnings.simplefilter("ignore")

        # Creates a biopython pdb object and starts to take informations from it.
        fh = open(os.path.join(pymod.structures_directory, self.copied_pdb_file_name), "rU")
        self.parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(self.pdb_file_base_name, fh)

        # Starts to iterate through the models in the biopython object.
        for model in self.parsed_biopython_structure.get_list():

            for chain in model.get_list():
                # These dictionaries will be filled with data needed to build Clusterseq and
                # Structure objects for each chain of the PDB file.
                parsed_chain = {
                    "original_id": None, # Chain ID in the PDB file.
                    "id": None, # The ID assigned in PyMod.
                    "max_residue": 0,
                    # Current chain HETATM residues.
                    "hetero_residues": [],
                    # Number of water molecules in the current chain.
                    "water_counter": 0,
                    "disulfide_bridges": [],
                    # This string is going to contain the sequence that will appear on the PyMod
                    # main window.
                    "sequence": "",
                    # This will be used to build an object of the PDB_chain_sequence class and it will
                    # contain a list of 'PDB_residue' objects for each residue of the current chain.
                    "pdb_sequence":[]
                }

                # Assigns a blank "X" chain id for PDB structures that do not specify chains id.
                parsed_chain["original_id"] = chain.id
                if chain.id != " ":
                    parsed_chain["id"] = chain.id
                elif chain.id == " ":
                    chain.id = "X"
                    parsed_chain["id"] = "X"

                # Finds the "chain_max_residue" for the current chain according to the informantion
                # in SEQRES part of the PDB file parsed before.
                for c in self.structure_chains_upper_limits:
                    if c["chain"] == parsed_chain["original_id"]:
                        parsed_chain["max_residue"] = c["upper_limit"]

                # Finds disulfide bridges in the current chain.
                for dsb in self.structure_disulfides:
                    if dsb.cys1_chain == parsed_chain["original_id"] or dsb.cys2_chain == parsed_chain["original_id"]:
                        parsed_chain["disulfide_bridges"].append(dsb)

                # Starts to build the sequences by parsing through every residue of the chain.
                for (id_position,residue) in enumerate(chain):
                    # Gets the 3 letter name of the current residue.
                    resname = residue.get_resname()
                    # get_id() returns something like: ('H_SCN', 1101, ' ').
                    residue_id = residue.get_id()
                    # Hetfield. Example: 'H_SCN' for an HETRES, while ' ' for a normal residue.
                    hetfield = residue_id[0]
                    # Number of the residue according to the PDB file.
                    pdb_position = residue_id[1]

                    # For HETATM residues.
                    if hetfield[0] == "H":
                        # Finds if the residue is a ligand (for example an ion, a metabolite or
                        # any kind of small molecule) or a modified residue belonging to protein
                        # (for example a phosphoserine, a seleno-cysteine, or a marked residue).
                        # This difference is important because when building a model and reading
                        # an alignment file, MODELLER wants modified residues to be inserted within
                        # the rest of the primary sequence of the template (because the HETATM lines
                        # for modified residues of a chain are usually inserted somewhere between
                        # the ATOM lines of the chain), while the "ligands" must be at the end of
                        # the template sequence.

                        # Builds a PDB_residue object.
                        one_letter_symbol = "X"
                        new_hetero_residue = PDB_residue(one_letter_symbol, resname, id_position, pdb_position, residue_type="het", chain_id=chain.id)

                        # Check if the current HETRES is a modres according to the info in the MODRES
                        # fields in the PDB file.
                        is_modified_residue = False
                        if self.found_het_lines and self.found_modres_lines:
                            for m in self.structure_modres:
                                if m["chain"] == parsed_chain["original_id"] and int(m["position"]) == int(pdb_position):
                                    is_modified_residue = True
                        # If the PDB file did not had an header, then try a different approach.
                        else:
                            pass

                        if is_modified_residue:
                            new_hetero_residue.set_hetres_type("modified-residue")
                            # If the HETRES is a "modified-residue" that is part of the protein
                            # try to assign to it its original unmodified amminoacid letter.
                            parsed_chain["sequence"] += one_letter_symbol
                        else:
                            new_hetero_residue.set_hetres_type("ligand")

                        # Sets its full name: it still needs to be implemented.
                        new_hetero_residue.set_hetres_full_name()

                        parsed_chain["pdb_sequence"].append(new_hetero_residue)
                        parsed_chain["hetero_residues"].append(new_hetero_residue)

                    # For water molecules.
                    elif hetfield == "W":
                        one_letter_symbol = "X"
                        new_water_molecule = PDB_residue(one_letter_symbol,resname,id_position,pdb_position,residue_type="water")
                        parsed_chain["pdb_sequence"].append(new_water_molecule)
                        parsed_chain["water_counter"] += 1

                    # For standard amminoacidic residues. Adds them to the primary sequence.
                    else:
                        one_letter_symbol = pymod.three2one(resname)
                        new_residue = PDB_residue(one_letter_symbol,resname,id_position,pdb_position,residue_type="standard")
                        parsed_chain["sequence"] += one_letter_symbol

                        # For cysteines involved in disulfide bridges: assigns cys_seq_numbers (the
                        # id of the cys in the sequence stored by pymod, they are usually different
                        # from the cys_number, the residue number in the PDB file).
                        # Bu it would be better to manually check for disulfied bonds in the PDB.
                        if one_letter_symbol == "C":
                            for dsb in parsed_chain["disulfide_bridges"]:
                                if dsb.bridge_type == "intrachain":
                                    if dsb.cys1_pdb_number == pdb_position:
                                        new_residue.set_disulfide_bridge(dsb)
                                        dsb.set_cys_seq_number(1,id_position)
                                    elif dsb.cys2_pdb_number == pdb_position:
                                        new_residue.set_disulfide_bridge(dsb)
                                        dsb.set_cys_seq_number(2,id_position)
                        parsed_chain["pdb_sequence"].append(new_residue)

                self.parsed_pdb_file_chains.append(parsed_chain)

            # Stops after having parsed the first "model" in the biopython "structure".
            # This is needed to import only the first model of NMR PDB files, which have multiple
            # models.
            break

        fh.close()
        warnings.simplefilter("always")


    def get_chains_ids(self):
        return [parsed_chain["id"] for parsed_chain in self.parsed_pdb_file_chains]


    def build_structure_objects(self, add_to_pymod_pdb_list = False, new_pdb_file_name=None):
        """
        Builds PyMod_element and Structure objects for each chain in the parsed PDB file.
        """
        self.chains_structure_objects = []
        self.chains_pymod_elements = []

        for parsed_chain in self.parsed_pdb_file_chains:
            # Builds header of the chain: it will be used to identify the Chain in PyMod (it will be
            # the header displayed to the user inside PyMod main window).
            header_name = None
            # PDB file that actually have a code will be named inside PyMod by using the their code.
            if self.pdb_code != None:
                header_name = str(self.pdb_code)+"_Chain:"+str(parsed_chain["id"])
            # If the name of the PDB files of the chains is specified, use it to name the files.
            elif new_pdb_file_name != None:
                header_name = new_pdb_file_name.replace(":","_")+"_Chain_"+str(parsed_chain["id"])
            # PDB files that don't have a code will be named inside python according to their
            # filename.
            else:
                header_name = str(self.pdb_file_base_name)+"_Chain:"+str(parsed_chain["id"])

            # The header_name is then used to build those two variables.
            corrected_record_header = pymod.correct_name(header_name.replace(":","_"))
            # Set the name of PDB file that will contain the current chain.
            current_chain_pdb_file_name_root = corrected_record_header

            # Builds a Structure object that will be given to the PyMod_element object that is
            # going to be created below.
            structure_object = Structure(
                biopython_structure = None,
                pdb_chain_sequence = parsed_chain["pdb_sequence"],
                pdb_chain_id = parsed_chain["id"],
                original_pdb_file_name = self.copied_pdb_file_name,
                chain_pdb_file_name_root = current_chain_pdb_file_name_root,
                hetero_residues = parsed_chain["hetero_residues"],
                water_molecules_count = parsed_chain["water_counter"],
                disulfides = parsed_chain["disulfide_bridges"])
            self.chains_structure_objects.append(structure_object)

            # Saves a PDB file with only the current chain of the first model of the structure.
            warnings.simplefilter("ignore")
            io=Bio.PDB.PDBIO()
            io.set_structure(self.parsed_biopython_structure)
            chain_selection = pymod.Select_chain_and_first_model(parsed_chain["id"])
            new_chain_pdb_file_name = os.path.join(pymod.structures_directory, current_chain_pdb_file_name_root+".pdb")
            io.save(new_chain_pdb_file_name, chain_selection)
            warnings.simplefilter("always")

            # Actually creates the PyMod_element objects that are going to store all the data of
            # the structure chains.
            element_sequence = parsed_chain["sequence"]
            chain_polymer_type = pymod.get_polymer_type(parsed_chain["pdb_sequence"])
            c = PyMod_element(
                    element_sequence, header_name,
                    structure = structure_object,
                    element_type = "sequence",
                    polymer_type = chain_polymer_type)
            # Appends the object to the list returned by this method.
            self.chains_pymod_elements.append(c)

        if add_to_pymod_pdb_list:
            self.add_to_pdb_list()

    # ---
    # Returns objects of all chains in the PDB file.
    # ---
    def get_chains_pymod_elements(self):
        return self.chains_pymod_elements


    def get_chains_structures(self):
        return self.chains_structure_objects

    # ---
    # Returns single object of single chains.
    # ---
    def get_chain_pymod_element(self, chain_id):
        """
        After the 'self.build_structure_objects()' method has been called, this can be used
        to return pymod elements of the corresponding to the chains specified in the 'chain_id'
        argument.
        """
        for element in self.chains_pymod_elements:
            if element.structure.pdb_chain_id == chain_id:
                return element


    def get_chain_structure(self, chain_id):
        """
        Used to return 'Structure' objects of the corresponding to the chains specified in the
        'chain_id' argument.
        """
        for structure in self.chains_structure_objects:
            if structure.pdb_chain_id == chain_id:
                return structure


    def crop_structure_chain(self, chain_id, adjust_to_sequence = None, add_to_pymod_pdb_list=True):
        """
        Once 'build_structure_objects()' has been used, this will edit the 'PyMod_element' and
        'Structure' objects corresponding to the 'chain_id' according to the sequence provided in
        the 'adjust_to_sequence' argument.
        Usually this is used when fetching a PDB file corresponding to some hit from a BLAST search,
        because hits in HSPs may have a shorter sequence with respect to the full PDB chain.
        This method can be called to crop the full 3D chain according to the hit sequence in the
        HSP (provided in the 'adjust_to_sequence' argument).
        """
        # Get the 'Structure' object of the chain specified in the 'chain_id' argument.
        t_element = self.get_chain_pymod_element(chain_id)
        t_sequence = t_element.my_sequence
        # And get the gapless sequence to which to adjust the cropped structure.
        h_sequence = str(adjust_to_sequence).replace("-","")
        # Align the two sequences using dynamic programming.
        ali = pmsm.global_pairwise_alignment(h_sequence, t_sequence, toss_modres=True)

        # If the sequences do not match, interrupt the process.
        if ali["id"] < 99.9:
            return False

        # Gets information about matching and missing residues in the two aligned sequences.
        pc = 0 # Alignment position counter.
        hc = 0 # Target residue counter.
        tc = 0 # PDB structure residue counter.
        matching_positions = [] # list of matching positions.
        missing_positions = [] # list of missing residues in the pdb structure with respect to the target sequence.
        for hr, tr in zip(ali["seq1"], ali["seq2"]):
            if hr != "-" and tr != "-" and hr == tr:
                matching_positions.append({"pc":pc,"hc":hc,"tc":tc})
            if tr == "-" and hr != "-":
                missing_positions.append({"pc":pc,"hc":hc,"tc":tc})
            if hr != "-":
                hc += 1
            if tr != "-":
                tc += 1
            pc += 1

        # Gets the starting and ending positions (using the PDB numeration) that will be used to
        # crop the 3D structure.
        start_position = None
        end_position = None
        for (id_position, residue) in enumerate(t_element.structure.get_all_residues_list()):
            if id_position == matching_positions[0]["tc"]:
                start_position = residue.pdb_position
            if id_position == matching_positions[-1]["tc"]:
                end_position = residue.pdb_position

        # Use PyMOL to build the new cropped structure.
        structure_root_name = t_element.structure.chain_pdb_file_name_root
        structure_file_shortcut = os.path.join(pymod.structures_directory, structure_root_name)
        # First loads the full PDB structure of the chain in PyMOL.
        cmd.load(structure_file_shortcut+".pdb", "full")
        # Select amminoacidic residues ranging from the starting and ending positions which define
        # the fragment to "excise".
        # cmd.select("ppfragment", "resi %s-%s and object full and not hetatm" % (start_position, end_position))
        cmd.select("ppfragment", "resi %s-%s and object full and not hetatm" % (start_position, end_position))
        # Join the selections and save a file PDB file of the cropped fragment.
        pdb_basename = os.path.splitext(t_element.structure.original_pdb_file_name)[0]
        cropped_structure_file_shortcut = os.path.join(pymod.structures_directory, pdb_basename)
        pymod.pymol_save("%s_cropped.pdb" % (cropped_structure_file_shortcut), "ppfragment")
        # Clean up the selections.
        cmd.delete("full")
        cmd.delete("ppfragment")

        # Builds a 'Parsed_pdb_file' object for the PDB file of the structure just saved.
        cpdb_file = Parsed_pdb_file(os.path.abspath("%s_cropped.pdb" % (cropped_structure_file_shortcut)))
        cpdb_file.parse_pdb_file()
        cpdb_file.build_structure_objects(add_to_pymod_pdb_list = False)

        # Updates the 'self.chains_pymod_elements' and 'self.chains_structure_objects' with the
        # sequence and structural information of the excised fragment.
        for ei, olement in enumerate(self.chains_pymod_elements):
            if olement.structure.pdb_chain_id == chain_id:
                self.chains_pymod_elements[ei] = cpdb_file.get_chain_pymod_element(chain_id)
                # Updates the sequence of the fragment to keep it in frame with the original
                # sequence provided in 'adjust_to_sequence' by including the target sequence indels.
                list_of_missing_positions = [p["hc"] for p in missing_positions]
                new_sequence = []
                adc = 0
                for i, p in enumerate(list(adjust_to_sequence)):
                    if adc in list_of_missing_positions:
                        new_sequence.append("-")
                    else:
                        # Right now modified residues are not included in the cropped structures,
                        # this prevents them from being included in the chain sequence.
                        if p != "X":
                            new_sequence.append(p)
                        else:
                            new_sequence.append("-")
                    if p != "-":
                        adc += 1
                new_sequence = "".join(new_sequence)
                self.chains_pymod_elements[ei].my_sequence = new_sequence

        for si, ostructure in enumerate(self.chains_structure_objects):
            if ostructure.pdb_chain_id == chain_id:
                self.chains_structure_objects[si] = cpdb_file.get_chain_structure(chain_id)
                break

        if add_to_pymod_pdb_list:
            self.add_to_pdb_list()

        # The sequence of the structure and the target sequences match.
        return True


    def add_to_pdb_list(self):
        # All the chains of the structure.
        if not self.copied_pdb_file_name in [p.pdb_file_name for p in pymod.pdb_list]:
            original_pdb_chains = map(lambda c: c.id, self.parsed_biopython_structure.get_chains())
            new_pdb = PDB_file(self.copied_pdb_file_name, original_pdb_chains, self.segment_structure, self.chains_pymod_elements)
            pymod.pdb_list.append(new_pdb)


# -----
# Modeling clusters.
# -----
class Modeling_cluster:
    def __init__(self,cluster):

        # This is actually the target sequence that is selected by the user.
        self.target = [c for c in pymod.get_children(cluster) if c.selected][0]
        # This is used to load the model into PyMOL when Modeller has done its job.
        self.target_name = pmos.clean_file_name(self.target.compact_header)

        # self.model_color=target.my_color
        self.aligned_elements_list = pymod.get_siblings(self.target)

        # Another for cycle to look for templates aligned to the target sequence.
        self.structure_list = []
        for aligned_sequence in self.aligned_elements_list:
            if aligned_sequence.has_structure() and aligned_sequence.is_model != True:
                # Populates the struct_list with templates to be displayed in the
                # modeling window.
                self.structure_list.append(aligned_sequence) # struct_list

        # This will contain a list objects from the Structure_frame class.
        self.structure_frame_list = []

        # Important list that is going to contain informations about the templates.
        self.templates_list = []
        # This list will be used to inform Modeller about which are the "known" sequences.
        # It will contain the headers of the templates.
        self.templates_namelist = []

        self.water_molecules_count = 0
        self.use_water_in_cluster = False
        self.use_template_complex_waters = False

        self.disulfides_frame = None
        self.target_with_cys = None
        if self.target.my_sequence.count("C") >= 2:
            self.target_with_cys = True
        else:
            self.target_with_cys = False

        self.symmetry_id = None
        self.apply_symmetry_restraints = None

        self.symmetry_var = None

        self.dictionary = {}

        self.model_elements_list = []

    def get_use_as_template_states(self):
        use_as_template_var_list = []
        for sf in self.structure_frame_list:
            use_as_template_var_list.append(sf.use_as_template_var)
        return use_as_template_var_list

    def set_block_index(self,index):
        self.block_index = index

    def set_water_molecules_number(self,n):
        self.water_molecules_number = n

    def has_structures_with_disulfides(self):
        disulfides = None
        if True in [e.structure.has_disulfides() for e in self.structure_list]:
            disulfides = True
        else:
            disulfides = False
        return disulfides

    def set_symmetry_id(self,symmetry_id):
        self.symmetry_id = symmetry_id

    def set_model_chain_id(self,chain_index):
        self.model_chain_id = chain_index

    def set_symmetry_var(self,symmetry_var):
        self.symmetry_var = symmetry_var

    def get_template_complex_chain(self):
        """
        Returns the 'PyMod_element' object if the template complex chain of this modeling cluster.
        """
        c = None
        for t in self.templates_list:
            if t in pymod.template_complex.get_pymod_elements():
                c = t
                break
        return c

    def has_ligands(self):
        ligands = False
        for t in self.templates_list:
            ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
            if ligand_count > 0:
                ligands = True
                break
        return ligands

    def template_complex_chain_has_ligands(self):
        ligands = False
        for t in self.templates_list:
            if t in pymod.template_complex.get_pymod_elements():
                ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
                if ligand_count > 0:
                    ligands = True
                break
        return ligands

    def switch_hetres_checkbutton_states(self,het_radio_button_state):
        for sf in self.structure_frame_list:
            # Activate.
            if het_radio_button_state == 1:
                sf.hetres_radiobutton_state = 1
                sf.activate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.activate_het_checkbuttons()
            # Inactivate.
            if het_radio_button_state == 0:
                sf.hetres_radiobutton_state = 0
                sf.inactivate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.inactivate_het_checkbuttons()

    def adjust_model_elements_sequence(self, remove_gaps = False):
        self.backup_target_sequence = None
        if not remove_gaps:
            self.backup_target_sequence = self.target.my_sequence
        else:
            self.backup_target_sequence = str(self.target.my_sequence).replace("-","")
        for model_element in self.model_elements_list:
            if model_element.unique_index != self.target.unique_index:
                model_element.my_sequence = self.backup_target_sequence

# ---
# PIR alignment class.
# ---
class PIR_alignment_sequence:
    def __init__(self,main_segment,ligands_segment,water_segment):
        self.main_segment = main_segment
        self.ligands_segment = ligands_segment
        self.water_segment = water_segment

    def get_single_chain(self,use_hetres = True, use_water = True):
        if use_hetres and use_water:
            return self.main_segment+self.ligands_segment+self.water_segment
        elif use_hetres and not use_water:
            return self.main_segment+self.ligands_segment
        else:
            return self.main_segment


# NOTE: Maybe just make a dictionary for this.
class Modeling_block:
    def __init__(self,pdb):
        self.first_line = ""
        self.second_line = ""
        self.pdb = pdb
        self.segment_list = []

class Original_Block:
    def __init__(self):
        self.segment_list = []

    def add_segment(self,segment):
        self.segment_list.append(segment)

    def generate_template_complex_block(self,template_complex=None):
        template_complex_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_template_complex()
            template_complex_segment_list.append(seq)
        tcbl = Modeling_block(template_complex.pdb_file_name)
        tcbl.first_line = ">P1;" + template_complex.pdb_file_name[:-4]
        tcbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (template_complex.pdb_file_name[:-4], "FIRST", template_complex.chains_list[0], "END", template_complex.chains_list[-1])
        tcbl.segment_list = template_complex_segment_list
        return tcbl

    def generate_additional_template_block(self,template=None):
        template_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_template(template)
            template_segment_list.append(seq)
        tid = template.structure.chain_pdb_file_name.replace(":","_")[:-4]
        tbl = Modeling_block(tid)
        tbl.first_line = ">P1;" + tid
        first_delimiter = "FIRST" # "FIRST"
        second_delimiter = "." # "LAST"
        tbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (tid, first_delimiter, ".", second_delimiter, ".")
        tbl.segment_list = template_segment_list
        return tbl

    def generate_target_block(self,target=None):
        self.target_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_target()
            self.target_segment_list.append(seq)
        tgid = pymod.modeller_target_name
        tgb = Modeling_block(tgid)
        tgb.first_line = ">P1;" + tgid
        tgb.second_line = "sequence:%s:%s:%s:%s:%s::::" % (tgid, "FIRST", ".", "LAST", ".")
        tgb.segment_list = self.target_segment_list
        return tgb

    def get_target_segment_list(self):
        target_segment_list = []
        for sg in self.segment_list:
            seg = sg.get_target_segment()
            target_segment_list.append(seg)
        return target_segment_list

# It's reduntant with Modeling_segment.
class Segment:

    def __init__(self,segment_type=None,modeling_cluster=None):
        self.segment_type = segment_type
        self.modeling_cluster = modeling_cluster
        self.get_template_complex_chain()

    def get_template_complex_chain(self):
        self.template_complex_chain = None
        for template_complex_chain in pymod.template_complex.get_pymod_elements():
            if self.segment_type[0] == template_complex_chain.structure.pdb_chain_id:
                self.template_complex_chain = template_complex_chain
                break

    # ---
    # Template complex.
    # ---
    def get_template_complex(self):
        seq = None
        # Template complex segments.
        if self.segment_type[0] != None:
            # Main segments.
            if self.segment_type[1] == "A":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    seq = self.template_complex_chain.pir_alignment_sequence.main_segment
                else:
                    seq = str(self.template_complex_chain.my_sequence).replace("-","")
            # Ligands segments.
            elif self.segment_type[1] == "H":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    seq = self.template_complex_chain.pir_alignment_sequence.ligands_segment
                else:
                    ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
                    seq = "."*ligand_count
            # Water segments.
            elif self.segment_type[1] == "W":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    # NOTE: just try to use
                    # seq = "w"*self.template_complex_chain.structure.water_molecules_count
                    # also here.
                    seq = self.template_complex_chain.pir_alignment_sequence.water_segment
                else:
                    seq = "w"*self.template_complex_chain.structure.water_molecules_count
        # Additional templates segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                # NOTE: write a method to get the ligands segment lenght in a mc.
                seq = "-"*len(self.modeling_cluster.templates_list[0].pir_alignment_sequence.ligands_segment)
            elif self.segment_type[1] == "W":
                seq = "-"*self.modeling_cluster.water_molecules_number

        return seq


    def get_template_complex_main_segment_in_alignment(self):
        seq = None
        if self.segment_type[0] in pymod.template_complex_selected_chain_list:
            seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.main_segment)
        else:
            seq = "-"*len(str(self.template_complex_chain.my_sequence).replace("-",""))
        return seq

    def get_template_complex_ligand_segment_in_alignment(self):
        seq = None
        if self.segment_type[0] in pymod.template_complex_selected_chain_list:
            seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.ligands_segment)
        else:
            ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
            seq = "-"*ligand_count
        return seq

    def get_template_complex_water_segment_in_alignment(self):
        pass

    # ---
    # Other templates.
    # ---
    def get_template(self,template):
        seq = None
        # Template complex segments.
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.main_segment
                else:
                    seq = self.get_template_complex_main_segment_in_alignment()

            elif self.segment_type[1] == "H":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.ligands_segment
                else:
                    seq = self.get_template_complex_ligand_segment_in_alignment()

            elif self.segment_type[1] == "W":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.water_segment)
                else:
                    seq = "-"*self.template_complex_chain.structure.water_molecules_count

        # Additional templates segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                if template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.ligands_segment
                else:
                    seq = "-"*len(template.pir_alignment_sequence.ligands_segment)
            elif self.segment_type[1] == "W":
                if template in self.modeling_cluster.templates_list:
                    if template.structure.water_state == 1:
                        seq = "w"*self.modeling_cluster.water_molecules_number
                    else:
                        seq = "-"*self.modeling_cluster.water_molecules_number
                else:
                    seq = "-"*self.modeling_cluster.water_molecules_number
        return seq

    # ---
    # Target.
    # ---
    def get_target(self):
        seq = None
        chain = self.get_template_complex_chain()
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                if self.modeling_cluster != None:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.main_segment
                else:
                    seq = self.get_template_complex_main_segment_in_alignment()

            elif self.segment_type[1] == "H":
                if self.modeling_cluster != None:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
                else:
                    seq = self.get_template_complex_ligand_segment_in_alignment()
            elif self.segment_type[1] == "W":
                if self.modeling_cluster != None:
                    if self.modeling_cluster.use_template_complex_waters:
                        seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment
                    else:
                        seq = "-"*self.template_complex_chain.structure.water_molecules_count
                else:
                    seq = "-"*self.template_complex_chain.structure.water_molecules_count
        # Extra segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
            elif self.segment_type[1] == "W":
                if self.modeling_cluster.use_template_complex_waters:
                    seq = "-"*chain.structure.water_molecules_count
                else:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment

        return seq

    def get_target_segment(self):
        seg = Modeling_segment(None) #self.modeling_cluster.target.pir_alignment_sequence.main_segment)
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None:
                    seg.set_segment_state(True)
                else:
                    seg.set_segment_state(False)
            elif self.segment_type[1] == "H":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None:
                    seg.set_segment_state(True)
                else:
                    seg.set_segment_state(False)
            elif self.segment_type[1] == "W":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None and self.modeling_cluster.use_template_complex_waters:
                        seg.set_segment_state(True)
                else:
                        seg.set_segment_state(False)

        # Extra segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                seg.set_chain_id(self.segment_type[0])
                seg.set_segment_state(True)
            elif self.segment_type[1] == "W":
                if self.modeling_cluster.use_template_complex_waters:
                    seg.set_chain_id(self.segment_type[0])
                    seg.set_segment_state(False)
                else:
                    seg.set_chain_id(self.segment_type[0])
                    seg.set_segment_state(True)
        return seg


# NOTE: use this also for templates.
class Modeling_segment(str):
    def set_segment_type(self,segment_type):
        self.type = segment_type
    def set_segment_state(self,use):
        self.use = use
    def set_chain_id(self,chain_id):
        self.chain_id = chain_id


# -----
# This class will be used to store a list of Symmetry_restraint_group objects.
# -----
class Symmetry_restraints_groups_list:

    def __init__(self):
        self.list_of_groups = []

    # Returns a list of Symmetry_restraints_group objects that contain at least the number of
    # sequences specified in the argument.
    def get_groups(self, min_number_of_sequences=0):
        return [g for g in self.list_of_groups if len(g.list_of_clusters) >= min_number_of_sequences]

    def add_group(self,symmetry_id):
        srg = Symmetry_restraints_group(symmetry_id)
        self.list_of_groups.append(srg)

    def get_group_by_id(self,symmetry_id):
        for g in self.list_of_groups:
            if g.id == symmetry_id:
                return g

# -----
# When performing multichain modeling, this will be used to identify a "symmetry restraints group",
# a group of target sequences that share the exact same sequence. By keeping track of these groups,
# PyMod can let the user apply symmetry restraints to those chains when using Modeller.
# -----
class Symmetry_restraints_group:
    def __init__(self,symmetry_id):
        # The "id" is just the target sequence stripped of all indels.
        self.id = symmetry_id
        # This will contain a list of Modeling_cluster objects that contain a target sequence
        # with the same sequence as the "id".
        self.list_of_clusters = []
        # This will be set to True if the user decides to apply symmetry restraints to this group
        # of target sequences.
        self.use = False
    # Adds a Modeling_cluster object to the group list_of_clusters.
    def add_cluster(self, modeling_cluster):
        self.list_of_clusters.append(modeling_cluster)


###################################################################################################
# CLASSES TO PREPRESENT PROTEIN STRUCTURES AND THEIR PROPERTIES.                                  #
###################################################################################################

class Structure:
    """
    # Class for respresenting the structure of a macromolecule loaded in PyMod.
    """
    def __init__(self,biopython_structure, pdb_chain_sequence, pdb_chain_id, original_pdb_file_name, chain_pdb_file_name_root, hetero_residues, water_molecules_count, disulfides):

        self.biopython_structure = biopython_structure
        # This will contain a list of 'PDB_residue' object for each residue of the chain.
        self.pdb_chain_sequence = pdb_chain_sequence
        self.pdb_chain_id = pdb_chain_id
        self.original_pdb_file_name = original_pdb_file_name  # Will contain something like: 1UBI.pdb
        # This corresponds to the name of the PyMOL object if this chain. In this way this attrbiute
        # cna be used to build selectors to interact with the objecty loaded in PyMOL.
        self.chain_pdb_file_name_root = chain_pdb_file_name_root # Will contain something like: 1UBI_Chain_A
        self.original_chain_pdb_file_name = chain_pdb_file_name_root+".pdb"

        # This variable is going to be updated when the structure is subject to a superposition or a
        # structural alignment. Initially it will be the same as original_chain_pdb_file_name but
        # when the structure is superposed or subjected to a structural alignment it will be updated
        # to something like 1UBI_Chain_A_aligned.pdb.
        self.chain_pdb_file_name = self.original_chain_pdb_file_name
        # This will contain iformation about hetero-atomic residues.
        self.hetero_residues = hetero_residues
        # Number of water molecules in the chain.
        self.water_molecules_count = water_molecules_count
        # A list that contains all disulfides bridges objects of the chain.
        self.disulfides = disulfides

    def set_modeling_information(self, seq_min, seq_max, hetres_option, hetres_map, water_state):
        self.seq_min = seq_min
        self.seq_max = seq_max
        self.hetres_option = hetres_option
        self.hetres_map = hetres_map
        self.water_state = water_state

    def clear_modeling_information(self):
        self.seq_min = 0
        self.seq_max = 0
        self.hetres_option = None
        self.hetres_map = None
        self.water_state = None

    def has_disulfides(self):
        disulfides = None
        if self.disulfides != []:
            disulfides = True
        else:
            disulfides = False
        return disulfides

    def has_waters(self):
        water = None
        if self.water_molecules_count > 0:
            water = True
        else:
            water = False
        return water


    def get_standard_residues_list(self, mode="list"):
        if mode == "list":
            return [res for res in self.pdb_chain_sequence if res.residue_type == "standard"]
        elif mode == "string":
            return reduce(lambda x,y: str(x)+str(y), [res for res in self.pdb_chain_sequence])
        elif mode == "tuple":
            return [(res.three_letter_code, res.pdb_position) for res in self.pdb_chain_sequence if res.residue_type == "standard"]


    def get_all_residues_list(self):
        return [res for res in self.pdb_chain_sequence if res.residue_type == "standard" or (res.residue_type == "het" and res.hetres_type == "modified-residue")]


    def get_chain_sequence(self):
        return reduce(lambda x,y: str(x)+str(y), self.get_all_residues_list())


class PDB_residue:
    """
    A class to represent the residues of a PDB chain.
    """
    def __init__(self,symbol,three_letter_code, id_position,pdb_position,residue_type="standard",chain_id=None):

        # Single letter identifier. Example: "Y" for tyrosine.
        self.symbol = symbol
        # Three letter identifier. Example: "TYR" for tyrosine. This is usefule for identifying
        # hetero-atomic residues. For example for a N-acetylglucosamine residue its value is: "NAG".
        self.three_letter_code = three_letter_code
        # Id of the residue. Goes from 0 up to the last residue. Also indels have an id.
        self.id = id_position
        # Number of the residue in the PDB file.
        self.pdb_position = pdb_position
        # If the electron density map is interpretable.
        self.pdb_present = True # False is missing

        # Not really necessary.
        # Id of the residue's chain in the PDB file.
        self.chain_id = chain_id
        # Can either be 'standart', 'het', 'water'.
        self.residue_type = residue_type

        # ---
        # Attributes for hetres.
        # ---
        # The attribute 'type' can either be "modified-residue" (like a phosphorylated
        # serine or threonine, or any non standard or covalalently modified residue) or
        # "ligand" [any metal ion, molecule or other stuff that is made of HETATMs (except water)
        # which is not covalentely bound to the protein].
        self.hetres_type = None
        # This info can only be extracted from the orgininal PDB file.
        self.hetres_full_name = None # It needs a method

        # ---
        # For disulfides.
        # ---
        self.disulfide_bridge = None


    def set_hetres_type(self,hetres_type):
        self.hetres_type = hetres_type


    def set_hetres_full_name(self):
        self.full_name = "complete name"
        # They can have really long names like:
        # HET    GD9  A2058      35
        # HETNAM     GD9 2-(1H-INDAZOL-4-YL)-6-{[4-(METHYLSULFONYL)PIPERAZIN-1-
        # HETNAM   2 GD9  YL]METHYL}-4-MORPHOLIN-4-YL-THIENO[3,2-D]PYRIMIDINE


    def set_disulfide_bridge(self,dsb):
        self.disulfide_bridge = dsb


    def __repr__(self):
        return self.symbol



class Disulfide_bridge:
    """
    # Class for disulfide bridges.
    """
    def __init__(self, cys1_pdb_number=0,cys1_chain="",cys2_pdb_number=0,cys2_chain="",distance=0,number_in_pdb=0, disulfide_id = None):

        # The id of the first cysteine residue in the PDB file.
        self.cys1_pdb_number = cys1_pdb_number
        self.cys1_chain = cys1_chain
        self.cys2_pdb_number = cys2_pdb_number
        self.cys2_chain = cys2_chain
        self.distance = distance
        self.number_in_pdb = number_in_pdb

        # The following attribute can be either "intrachain" (the 2 cys belong to the same
        # polypeptide chain) or "interchain" (the 2 cys belong to two different polypeptide chains).
        self.bridge_type = ""
        if self.cys1_chain == self.cys2_chain:
            self.bridge_type = "intrachain"
        elif self.cys1_chain != self.cys2_chain:
            self.bridge_type = "interchain"

        self.disulfide_id = disulfide_id

    # This sets the number of the cysteine in the sequence created by Pymod. Most often it is
    # different from the number from the PDB file.
    def set_cys_seq_number(self,cys,number):
        if cys == 1:
            self.cys1_seq_number = number
        elif cys == 2:
            self.cys2_seq_number = number


class PDB_file:

    def __init__(self,pdb_file_name,chains_list,segment_structure,clusterseq_elements):
        self.pdb_file_name = pdb_file_name
        # The list of chains ids ["A", "B", ...] in the original PDB file. This is needed when
        # performing multiple chain modeling.
        self.chains_list = chains_list
        self.segment_structure = segment_structure
        # A list of clusterseq elements for each chain of the structure.
        self.clusterseq_elements = clusterseq_elements

    def get_chain_by_id(self,chain_id):
        """
        Returns the 'PyMod_element' object corresponding to the chain_id argument.
        """
        chain_element = None
        for element in self.get_pymod_elements():
            if element.structure.pdb_chain_id == chain_id:
                chain_element = element
                break
        return chain_element

    def get_pymod_elements(self):
        return self.clusterseq_elements

    #----------------------------------------------------------
    # Methods needed when saving and loading a PyMod project. -
    #----------------------------------------------------------

    def get_unique_index(self, element):
        return element.unique_index

    attributes_to_store = ["pdb_file_name",
                           "chains_list",
                           "segment_structure",
                           "chains_list",
                           "clusterseq_elements"]
    lists_to_pickle = {} # {"chains_list": ("get_unique_index", "unpickle")}
    attributes_to_pickle = []

    def __getstate__(self):
        return pymod_pickle(self)


class Modeling_session:
    """
    Class for containing information on modeling sessions.
    """
    def __init__(self, session_id):
        self.session_id = session_id
        self.assessment_table_data = None
        self.session_profile = None
        self.full_models = []


class Full_model:
    """
    Class for containing information on models built in a modeling session. Object of this class
    will be contained in the '.full_models' attribute of 'Modeling_session' class objects.
    """
    def __init__(self, original_file_path):
        self.original_file_path = original_file_path
        self.model_name = os.path.basename(self.original_file_path)[:-4]
        self.model_profile = None
        self.assessment_data = None


class Alignment:
    """
    Class for alignments.
    """
    def __init__(self, alignment_algorithm, alignment_id, initial_number_of_sequence=None):
        """
        alignment_algorithm: the algorithm used to perform the alignment
        alignment_id: an int value that identifies the alignmente object
        """
        self.algorithm = alignment_algorithm
        self.id = alignment_id
        self.initial_number_of_sequence = initial_number_of_sequence
        self.rmsd_list = None

    def set_dnd_file_path(self, dnd_file_path):
        self.dnd_file_path = dnd_file_path

    def get_dnd_file_path(self):
        return self.dnd_file_path

    def set_rmsd_list(self, rmsd_list):
        self.rmsd_list = rmsd_list


###################################################################################################
# PyMod_element class.                                                                            #
###################################################################################################

class PyMod_element: # ClusterSeq
    """
    A class that stores all the informations of a sequence. The objects of this class are the sequences
    (both from FASTA and PDB files) that appear on the left column or alignments elements.
    """
    def __init__(self, record_seq, record_header, full_original_header=None, element_type="sequence", color="white", structure=None, alignment_object=None, adjust_header=True, adjust_sequence=True, polymer_type="protein"):
        """
        record_seq: sequence to be displayed.
        record_header: from this it will be built the title to be displayed.
        fixed header: ...
        full original header: ...
        """

        # It is a unique id to identify each element. It is given by the "unique_index" of the
        # "pymod" object. This values will be assigned when the element is added to the list
        # of PyMod element objects by the '.add_element_to_pymod' method of the PyMod class.
        self.unique_index = None
        self.mother_index = None
        self.child_index = None

        # Mother elements are those elements outside clusters or cluster elements themselves
        # (alignments and similarity searches).
        self.is_mother = None
        # Children elements are those elements found inside clusters.
        self.is_child = None

        # This is needed to display the sequences in the right order. It will be defined inside the
        # "create_entry()" method.
        self.grid_position = None

        # Forms the header name.
        self.set_header_name(record_header, adjust_header)
        # The full original header.
        if full_original_header != None:
            self.full_original_header = full_original_header
        else:
            self.full_original_header = record_header

        # Sets the 'my_sequence' attribute. The primary sequence of an element. If the sequence is
        # changed or aligned by the user, it will be modified to include indels.
        self.set_sequence(record_seq, adjust_sequence)

        # A string that will be used to build .pir alignment files to be used as input by MODELLER.
        self.pir_alignment_string = ""

        self.annotations = {}

        # Its value is False when the sequence header is not selected, it becomes True when the user
        # selects the sequence by left-clicking on its header.
        self.selected = False

        # This will be set to 'True' when the lement is displayed. When an element gets hidden
        # (for example every time the 'gridder' method of the 'PyMod' class is called to refresh
        # PyMod's main window), thi will be set again to 'False'.
        self.is_shown = False

        # This are needed to build the mother buttons when the mothers have children.
        self.show_children = None
        self.button_state = '-'
        self.mother_button_color = "gray"

        # Defines the way the sequence has to be colored.
        #     - regular: the sequence will be colored with the color defined in self.my_color
        #     - residue:
        #     - secondary-observed:
        #     - secondary-predicted:
        #     - campo:
        #     - dope-score:
        #     - delta-dope-score:
        #     - modeller-b-factor:
        #     - scr:
        self.color_by = "regular"

        # Name of the color of the sequence when its "color_by" attribute is set to "regular".
        self.my_color = color

        # This will contain a list containing ids specifying the kind of secondary structure assigned
        # by PyMOL DSS algorithm to each residue of the sequence.
        self.pymol_dss_list = []
        self.psipred_elements_list = []
        self.campo_scores = []
        self.dope_scores = []
        self.dope_items = []

        # This variable defines if the object is a:
        #     - sequence from a sequence or structure file: "sequence"
        #     - an alingment : "alignment"
        #     - a similariry search cluster: "blast-cluster"
        self.element_type = element_type # element_type

        # True if the sequence has been used for a BLAST search as a query.
        self.is_blast_query = False

        # Lead sequences are children that will appear in top of their clusters and when the cluster
        # is collapsed will be the only one to be displayed.
        self.is_lead = False

        # If the sequence has been used as a bridge to perform an alignment joining.
        self.is_bridge = False

        # This attribute is going to contain an object of the Structure class. If the element does
        # not have a structure its value will be 0.
        self.structure = structure

        # Alignment elements and BLAST-search elements will have an Alignment class object stored in
        # this attribute.
        self.alignment = alignment_object

        # Modeling part.
        self.is_model = False # For models. It should be made when the modeling process has been done.
        # When the user builds a single chain model of this sequence this counter will increase by 1.
        # If the user builds more models at the same time, this will increase by 1 for each model
        # built.
        self.models_count = 0

        self.polymer_type = polymer_type


    def set_header_name(self, header, adjust_header=True):
        if adjust_header:
            self.my_header_fix = pymod.build_header_string(header)
            # Just the header. This will be displayed in PyMod main window.
            self.my_header = pymod.correct_name(self.my_header_fix)
        else:
            self.my_header_fix = header
            self.my_header = header
        # A compact header.
        self.compact_header = self.get_compact_header(self.my_header)


    def get_compact_header(self, header=None):
        """
        This is needed to build a compact version of the headers that will be used in various
        instances (for example when the identity matrix table for some alignment is displayed, or
        when assigning the name of Modeller ouputs files) in order to save space.
        """
        compact_header = ""
        if header == None:
            header = self.my_header

        def crop_header(h):
            reduced_length = 11 # Like the old Pymod.
            return h[0:reduced_length]

        # Uniprot (sp/tr) and gi entries.
        # From something like "sp|Q9H492|MLP3A_HUMAN" it will build "sp|Q92622".
        if header.startswith("sp|") or header.startswith("tr|") or header.startswith("gi|"):
            so = re.search( r'(\w{2})\|(\S{6,9})\|',header)
            if so:
                compact_header = so.group(1)+"|"+so.group(2) # +"|"
            else:
                compact_header = crop_header(header)
        # Sequences imported from PDB files using the open_pdb_file() method.
        elif header[-8:-1] == "_Chain_":
            if len(header) == 12: # If it is the name of sequence with a 4 letter PDB id.
                compact_header=header[0:4]+"_"+header[-1]
            else:
                compact_header=crop_header(header[0:-8])+"_"+header[-1]
        # Other kind of headers. Just crop the header.
        else:
            compact_header = crop_header(header)

        return compact_header


    def set_sequence(self, sequence, adjust_sequence=True):
        if adjust_sequence:
            self.my_sequence = pymod.correct_sequence(sequence)
        else:
            self.my_sequence = sequence


    def set_pir_alignment_sequence(self, pir_alignment_sequence):
        self.pir_alignment_sequence = pir_alignment_sequence


    def set_dope_scores(self, dope_scores):
        self.dope_scores = []
        for score, residue in zip(dope_scores, self.structure.get_all_residues_list()):
            self.dope_scores.append(score)


    def update_element(self, new_sequence=None, new_header=None, new_structure=None):
        if new_sequence != None:
            self.set_sequence(new_sequence)
            self.update_element_data()
        if new_header != None:
            self.set_header_name(new_header)
        if new_structure != None:
            self.structure = new_structure


    def update_element_data(self):
        self.pymol_dss_list = []
        self.psipred_elements_list = []
        self.campo_scores = []
        self.dope_scores = []
        self.dope_items = []
        self.color_by = "regular"


    def remove_indels(self):
        self.my_sequence = str(self.my_sequence).replace("-","")


    def is_cluster_element(self):
        if self.is_mother and self.element_type in ("alignment", "blast-search"):
            return True
        else:
            return False


    def has_structure(self):
        if self.structure != None:
            return True
        else:
            return False


    def can_be_modeled(self):
        """
        Check if the element can be used to build a model with Modeller. Only sequences not coming
        from a PDB file chain imported by the user can be used.
        """
        r = False
        # Exlude cluster elements (alignments and BLAST-searches).
        if not self.is_cluster_element():
           if self.has_structure():
               # If the element have as a structure a model, than it can be used to build new models.
               if self.is_model:
                   r = True
               else:
                   r = False
           # The element is a primary sequence imported by the user.
           else:
               r = True
        return r


    def has_models(self):
        pass


    def has_assigned_secondary_structure(self):
        if self.pymol_dss_list != []:
            return True
        else:
            return False

    def has_predicted_secondary_structure(self):
        if self.psipred_elements_list != []:
            return True
        else:
            return False


    def has_campo_scores(self):
        if self.campo_scores != []:
            return True
        else:
            return False


    def is_being_showed(self):
        """
        Return 'False' is the sequence is not currently being showed, that is, if it is hidden
        inside its collapsed cluster.
        """
        return False


    def is_lead_of_collapsed_cluster(self):
        if self.is_child and self.is_lead:
            if not pymod.get_mother(self).show_children:
                return True
            else:
                return False
        else:
            return False


    def create_entry(self, grid_position=0 ,font_size="14"):
        """
        This method allows to display the sequence and the header in the window. It is called in the
        "Grid()" method of the "PyMod" class.
        """

        self.grid_position = grid_position

        # For inserting indels in the sequence through the mouse.
        self.mousepos = 0
        self.seqpos   = -1

        # Sets the font type and size.
        self.sequence_font_type = pmgi.fixed_width_font
        self.sequence_font_size = font_size
        self.sequence_font = self.sequence_font_type + " " + self.sequence_font_size # The default one is "courier 14".
        self.bg_color = "black"

        # Display the sequence text widget.
        if not self.is_cluster_element():
            self.build_element_header_entry()
            self.build_element_sequence_entry()

        elif self.is_cluster_element():
            if self.show_children:
                self.build_element_header_entry()
                self.build_element_sequence_entry()

            elif not self.show_children:
                cluster_lead = pymod.get_cluster_lead(self)
                if cluster_lead != None:
                    cluster_lead.create_entry(grid_position = grid_position, font_size=font_size)
                else:
                    self.build_element_header_entry()
            self.build_cluster_button()


    def build_element_sequence_entry(self, grid_position=None):
        """
        Creates a tkinter Text inside the right pane to display the sequence of the element.
        """
        if grid_position == None:
            grid_position = self.grid_position

        self.sequence_entry = Text(pymod.rightpan.interior(),
            font = self.sequence_font,
            cursor = "hand2",
            wrap=NONE,
            height=1,
            borderwidth=0,
            highlightcolor=self.bg_color,
            highlightbackground=self.bg_color,
            foreground = self.my_color,
            background = self.bg_color,
            exportselection=0,
            selectbackground= self.bg_color,
            selectforeground=self.my_color,
            selectborderwidth=0,
            width = len(self.my_sequence)) # The length of the entry is equal to the length of the sequence.
        try:
            self.sequence_entry.configure(inactiveselectbackground=self.bg_color)
        except:
            pass

        # Enters some sequence in the Text widget built above and colors it according to the element
        # current color scheme.
        self.build_text_to_display()

        # Modifier that allows to display the symbol '|_' of a child sequence.
        if self.is_child:
            self.sonsign = StringVar()
            self.sonsign.set("|_")
            if self.is_blast_query:
                self.sonsign.set("|q")
            elif self.is_lead:
                self.sonsign.set("|l")
            elif self.is_bridge:
                self.sonsign.set("|b")

            # Creates a sequence entry inside the right-frame.
            trattino=Entry(pymod.rightpan.interior(), font = self.sequence_font, cursor = "hand2",
                           textvariable=self.sonsign, bd=0, state = DISABLED,
                           disabledforeground = 'white', disabledbackground = self.bg_color,
                           highlightbackground= self.bg_color, justify = LEFT, width = 2 )
            trattino.grid(column =0, row = grid_position, sticky='nw', padx=0, pady=0,ipadx=0,ipady=0)

        # Actually grids the sequence Text widget.
        self.sequence_entry.grid(column = 2, row = grid_position, sticky='nw', padx=5)

        # Builds the sequence popup menu and binds events to it.
        self.build_right_popup_menu()
        self.bind_events_to_sequence_entry()


    def build_cluster_button(self, grid_position = None):
        """
        Method that creates a button for displaying/hiding a cluster sequences. This is called inside
        pymod.gridder().
        """
        if grid_position == None:
            grid_position = self.grid_position

        plusminus=StringVar()
        plusminus.set(self.button_state)

        # It's not a Button widget, it's an Entry widget (more customizable).
        button=Entry(pymod.rightpan.interior(), font = self.sequence_font,
                     cursor = "hand2", textvariable=plusminus,
                     relief="ridge", bd=0,
                     state = DISABLED, disabledforeground = 'white',
                     disabledbackground = self.mother_button_color, highlightbackground='black',
                     justify = CENTER, width = 1 )

        button.grid(column = 0, row = grid_position, sticky='nw', padx=5, pady=0,ipadx=3,ipady=0)
        # Binding the mouse event
        button.bind("<Button-1>", self.cluster_button_click)


    def cluster_button_click(self,event):
        """
        Creates the mouse event for the clicking the Button. It is used to toggle the children of
        the cluster.
        """
        if self.show_children:
            self.collapse_cluster()
        elif not self.show_children:
            self.expand_cluster()

    def expand_cluster(self):
        self.button_state = '-'
        self.mother_button_color = "gray"
        self.show_children = True
        pymod.gridder()

    def collapse_cluster(self):
        self.button_state = '+'
        self.mother_button_color = "red"
        self.show_children = False
        pymod.gridder()


    def build_element_header_entry(self):
        """
        Creates a tkinter Entry inside the left pane to display the header of the sequence.
        """
        # This is used only here to set the textvarialble of the entry as the header of the sequence.
        self.header_entry_var = StringVar()
        self.header_entry_var.set(self.my_header)

        self.header_entry = Entry(pymod.leftpan.interior(),
            font = self.sequence_font,
            cursor = "hand2",
            textvariable= self.header_entry_var,
            bd=0,
            highlightcolor='black',
            highlightbackground= self.bg_color,
            state = DISABLED,
            disabledforeground = 'red',
            disabledbackground = self.bg_color,
            selectbackground = 'green',
            justify = LEFT,
            width = int(len(self.header_entry_var.get())) )

        # MAYBE THIS DOES NOT HAVE TO BE ASSIGNED EVERY TIME THE ENTRIES ARE DISPLAYED.
        # Left menu object building and binding of the mouse events to the entries.
        self.build_left_popup_menu()
        self.bind_events_to_header_entry()
        self.header_entry.grid(row = self.grid_position, sticky='nw')

        # Marks the element as being 'showed' in PyMod's main window.
        self.is_shown = True


    ###############################################################################################
    # METHODS USED WHEN INTERACTING WITH THE HEADER ENTRY (IN THE LEFT PANE).                     #
    ###############################################################################################

    #################################################################
    # Mouse events and their bindings for the header Entry.         #
    #################################################################

    def bind_events_to_header_entry(self):
        self.header_entry.bind("<Button-1>", self.on_header_left_click)
        self.header_entry.bind("<Motion>", self.display_protname)
        if self.has_structure():
            self.header_entry.bind("<Button-2>", self.click_structure_with_middle_button)
        self.header_entry.bind("<ButtonRelease-3>", self.on_header_right_click)

    # Select/Unselect a protein clicking on its name on the left pane.
    def on_header_left_click(self,event):
        self.toggle_element()

    # Allows to show the protein name in the bottom frame 'pymod.sequence_name_bar'
    def display_protname(self,event):
            protein_name = self.full_original_header # self.header_entry.get()
            pymod.sequence_name_bar.helpmessage(protein_name)

    def click_structure_with_middle_button(self,event=None):
            # Shows the structure and centers if the sequence is selected in Pymod.
            if self.selected:
                """
                active_in_pymol = True
                if active_in_pymol:
                    # Centers the structure.
                    self.center_chain_in_pymol()
                    else:
                """
                self.show_chain_in_pymol()
                self.center_chain_in_pymol()
            # If the sequence is not selected in Pymod, hide it in PyMOL.
            else:
                self.hide_chain_in_pymol()

    # A popup menu in the left frame to interact with the sequence
    def on_header_right_click(self,event):
        try:
            self.header_entry["disabledbackground"] = 'grey'
            self.popup_menu_left.tk_popup(event.x_root, event.y_root, 0)
        except:
            pass
        #popup_menu.grab_release()
        self.header_entry["disabledbackground"] = 'black'


    ###################################
    # Selection of elements.          #
    ###################################

    # Toggles elements.
    def toggle_element(self):
        if self.is_mother:
            self.toggle_mother_element()
        elif self.is_child:
            if self.is_lead_of_collapsed_cluster():
                self.toggle_lead_element()
            else:
                self.toggle_child_element()

    # Toggle a mother element.
    def toggle_mother_element(self):
        # Inactivate.
        if self.selected:
            self.deselect_element()
            # Deselects also all the children.
            if self.is_cluster_element():
                for c in pymod.get_children(self):
                        if c.selected:
                            c.deselect_element()
        # Activate.
        else:
            self.select_element()
            # Activate also all the children!
            if self.is_cluster_element():
                for c in pymod.get_children(self):
                        if not c.selected:
                            c.select_element()

    # Toggle a child element.
    def toggle_child_element(self):
        mother = pymod.get_mother(self)
        # Inactivate.
        if self.selected:
            # Modify the mother and the siblings according to what happens to the children.
            if not mother.selected:
                siblings = pymod.get_siblings(self)
                # If it is not the last activated children in the cluster.
                if True in [s.selected for s in siblings]:
                    mother.deselect_element(is_in_cluster=True)
                    self.deselect_element(is_in_cluster=True)
                # If it is the last children to be inactivated.
                else:
                    mother.deselect_element()
                    for s in siblings:
                        s.deselect_element()
                    self.deselect_element()
            else:
                self.deselect_element(is_in_cluster=True)
                mother.deselect_element(is_in_cluster=True)

        # Activate.
        else:
            self.select_element()
            # If the mother is not selected and if by selecting this child, all the children
            # are selected, also selects the mother.
            if not mother.selected:
                # If it is the last inactivated children in the cluster, then by selecting it all the
                # elements in the cluster are selected and the mother is also selected.
                siblings = pymod.get_siblings(self)
                if not False in [c.selected for c in siblings]:
                    mother.select_element()
                else:
                    # Used to make the mother "gray".
                    mother.deselect_element(is_in_cluster=True)
                    # Used to make the siblings "gray".
                    for s in siblings:
                        if not s.selected:
                            s.deselect_element(is_in_cluster=True)


    # Toggle a lead element when a cluster is collapsed.
    def toggle_lead_element(self):
        if self.selected:
            self.deselect_element()
        else:
            self.select_element()

    # Selects an element.
    def select_element(self,is_in_cluster=False):
        self.selected = True
        if self.is_shown:
            if not is_in_cluster:
                self.header_entry["disabledforeground"] = 'green'
            else:
                self.header_entry["disabledforeground"] = 'green'

    # Deselects an element.
    def deselect_element(self, is_in_cluster=False):
        self.selected = False
        if self.is_shown:
            if not is_in_cluster:
                self.header_entry["disabledforeground"] = 'red'
            else:
                self.header_entry["disabledforeground"] = 'ghost white'

    # The two following methods are used only when the user clicks on the mother of a collapsed
    # cluster. They will select/deselect its children, without changing their color.
    def select_hidden_child(self):
        self.selected = True

    def deselect_hidden_child(self):
        self.selected = False


    ###################################
    # Interaction with PyMOL.         #
    ###################################

    # -----
    # Chains.
    # -----
    def select_chain_in_pymol(self,selection_name="pymod_selection"):
        sel = self.build_chain_selector_for_pymol(None)
        cmd.select(selection, sel)

    def center_chain_in_pymol(self):
        sel = self.build_chain_selector_for_pymol(None)
        cmd.center(sel)

    def hide_chain_in_pymol(self):
        # Use enable or disable?
        sel = self.build_chain_selector_for_pymol()
        cmd.disable(sel)

    def show_chain_in_pymol(self):
        sel = self.build_chain_selector_for_pymol()
        cmd.enable(sel)

    def show_selected_chains_in_pymol(self):
        for e in pymod.get_selected_sequences():
            e.show_chain_in_pymol()

    def hide_selected_chains_in_pymol(self):
        for e in pymod.get_selected_sequences():
            e.hide_chain_in_pymol()

    # Gets a selector for the object chain in PyMOL.
    def build_chain_selector_for_pymol(self,chain_id=None,pdb_header=None):
        # It should be updated, because some structural alignment or superposition stuff alters the
        # name of the sequences.
        assert(self.has_structure())
        selector = self.structure.chain_pdb_file_name_root
        return selector

    # -----
    # Residues.
    # -----
    # The following two methods are called when the user interacts with the sequence popup menu.
    def select_residue_in_pymol(self):
        res_id = self.get_highlighted_residue_position(pdb_position=True)
        sel = self.build_residue_selector_for_pymol(res_id)
        cmd.select("pymod_selection", sel)
        # cmd.color("white","pymod_selection")

    def center_residue_in_pymol(self,event=None):
        res_id = self.get_highlighted_residue_position(pdb_position=True)
        sel = self.build_residue_selector_for_pymol(res_id)
        cmd.center(sel)

    # Gets the correspondig selector in PyMOL.
    def build_residue_selector_for_pymol(self,res_id,pdb_header=None):
        # Selectors that work:
        #     /1UBI_Chain_A//A/LEU`43/CA
        #     1UBI_Chain_A and resi 43
        selector = self.build_chain_selector_for_pymol() + " and resi " + str(res_id)
        return selector


    #################################################################
    # Structure and methods of the left pane popup menu.            #
    #################################################################

    ###################################
    # Menus building.                 #
    ###################################

    # -----
    # Single elements menus.
    # -----
    def build_left_popup_menu(self):
        """
        Builds the popup menu that appears when the user clicks with the left button on the
        sequence title in the left pan.
        """

        self.popup_menu_left = Menu(pymod.parent, tearoff=0, bg='white', activebackground='black', activeforeground='white', postcommand=self.update_left_popup_menu)

        # Builds a popup menu for sequence elements.
        if not self.is_cluster_element():
            self.build_sequence_left_popup_menu()
        # For cluster elements ("alignment" or "blast-search" elements).
        else:
            self.build_cluster_popup_menu(self.popup_menu_left, mode="cluster", extra_spacer=True)

        # Selection menu. It will be activated only when there is more than one seleceted
        # sequence and the user clicks on some element with the mouse left-button.
        self.selection_menu = Menu(self.popup_menu_left, tearoff=0, bg='white',
            activebackground='black', activeforeground='white')
        self.popup_menu_left.add_cascade(menu=self.selection_menu, label="Selection", state=DISABLED)


    def build_sequence_left_popup_menu(self):

        # If the sequence is a lead of a cluster build the "Cluster" menu, to manage the cluster.
        if self.is_lead_of_collapsed_cluster():
            self.cluster_lead_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.build_cluster_popup_menu(self.cluster_lead_menu, mode="lead")
            self.popup_menu_left.add_cascade(menu=self.cluster_lead_menu, label="Cluster")
            self.popup_menu_left.add_separator()

        # Build the "Sequence" menu.
        self.build_sequence_menu()
        self.popup_menu_left.add_separator()

        # Build the "Color" menu.
        self.build_color_menu()
        self.popup_menu_left.add_separator()

        # Build the "Structure" menu.
        self.build_structure_menu()
        self.popup_menu_left.add_separator()

        # Build the "Cluster Options" menu.
        if self.is_child:
            self.build_cluster_options_menu()
            self.popup_menu_left.add_separator()


    def build_cluster_popup_menu(self, target_menu, mode="cluster", extra_spacer=False):
        self.build_cluster_edit_menu(target_menu)
        target_menu.add_separator()
        self.build_cluster_color_menu(target_menu)
        if extra_spacer:
            target_menu.add_separator()


    def update_left_popup_menu(self):
        """
        Activates the "Selection" item when at least two elements are selected.
        In order to make this work the "Selection" item always has to be in the last position in all
        kind of menus.
        """
        selected_sequences = pymod.get_selected_sequences()
        if len(selected_sequences) > 1: # and not self.is_cluster_element():
            self.popup_menu_left.entryconfig(self.popup_menu_left.index(END), state=NORMAL)
            self.build_selection_menu()
        elif not self.is_cluster_element():
            self.popup_menu_left.entryconfig(self.popup_menu_left.index(END), state=DISABLED)


    def build_sequence_menu(self):
        """
        Submenu with options for manipulating a sequence loaded in PyMod.
        """
        self.sequence_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.sequence_menu.add_command(label="Save Sequence to File", command=self.save_sequence_from_left_pane)
        self.sequence_menu.add_command(label="Copy Sequence to Clipboard", command=self.copy_sequence_to_clipboard)
        self.sequence_menu.add_separator()
        edit_command = None
        if not self.has_structure():
            edit_command = self.edit_sequence
        else:
            edit_command = self.edit_structure
        if not self.has_structure():
            self.sequence_menu.add_command(label="Edit Sequence", command=edit_command)
            self.sequence_menu.add_separator()
        self.sequence_menu.add_command(label="Duplicate Sequence",command=self.duplicate_sequence_from_the_left_pane)
        self.sequence_menu.add_command(label="Delete Sequence", command=self.delete_sequence_from_the_left_pane)
        self.popup_menu_left.add_cascade(menu=self.sequence_menu, label="Sequence")


    def build_color_menu(self):
        """
        Color submenu containing all the option to color for a single sequence.
        """
        self.color_menu = Menu(self.popup_menu_left,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        self.regular_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        for color in pmdt.regular_colours:
            self.regular_colors_menu.add_command(label=color, command = lambda c=color: pymod.color_selection("single",self,"regular",c))
        self.color_menu.add_cascade(menu=self.regular_colors_menu, label="Color whole Sequence by")
        self.color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        self.residues_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.residues_colors_menu.add_command(label="Polarity",command=lambda: pymod.color_selection("single", self, "residue"))
        self.color_menu.add_cascade(menu=self.residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        if self.can_be_colored_by_secondary_structure():
            self.color_menu.add_separator()
            self.sec_str_color_menu = Menu(self.color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            if self.has_structure():
                self.sec_str_color_menu.add_command(label="Observed", command=lambda: pymod.color_selection("single", self, "secondary-observed"))
            if self.has_predicted_secondary_structure():
                self.sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: pymod.color_selection("single", self, "secondary-predicted"))
            self.color_menu.add_cascade(menu=self.sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if self.can_be_colored_by_conservation():
            self.color_menu.add_separator()
            self.conservation_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: pymod.color_selection("single", self, "campo-scores"))
            self.color_menu.add_cascade(menu=self.conservation_colors_menu, label="By Convservation")

        # Energy colors.
        if self.can_be_colored_by_energy():
            self.color_menu.add_separator()
            self.energy_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.energy_colors_menu.add_command(label="DOPE scores",command=lambda: pymod.color_selection("single", self, "dope"))
            self.color_menu.add_cascade(menu=self.energy_colors_menu, label="By Energy")

        self.popup_menu_left.add_cascade(menu=self.color_menu, label="Color")


    def can_be_colored_by_secondary_structure(self):
        """
        Returns True if the element has an associated structure or has a secondary structure
        prediction.
        """
        if self.has_structure() or self.has_predicted_secondary_structure():
            return True
        else:
            return False

    def can_be_colored_by_conservation(self):
        if self.has_campo_scores():
            return True
        else:
            return False

    def can_be_colored_by_energy(self):
        if self.dope_items != []:
            return True
        else:
            return False


    def build_structure_menu(self):
        """
        Submenu for elements that have a structure loaded in PyMOL.
        """
        self.structure_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.has_structure():
            self.structure_menu.add_command(label="Center Chain in PyMOL", command=self.center_chain_in_pymol)
            # A switch could be nice.
            self.structure_menu.add_command(label="Show Chain in PyMOL", command=self.show_chain_in_pymol)
            self.structure_menu.add_command(label="Hide Chain in PyMOL", command=self.hide_chain_in_pymol)
            self.structure_menu.add_separator()
            self.structure_menu.add_command(label="PDB Chain Information", command=pymod.show_pdb_info)
        else:
            if self.pdb_is_fetchable():
                self.structure_menu.add_command(label="Fetch PDB File", command = lambda: pymod.fetch_pdb_files("single", self))
                self.structure_menu.add_separator()
            self.structure_menu.add_command(label="Associate 3D Structure", command=lambda: pymod.associate_structure_from_popup_menu(self))
        self.popup_menu_left.add_cascade(menu=self.structure_menu, label="Structure")


    def build_cluster_options_menu(self):
        """
        Submenu with options to manage a sequence within its cluster.
        """
        self.cluster_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_menu.add_command(label="Extract Sequence from Cluster", command=self.extract_from_cluster)
        if not self.is_lead:
            self.cluster_menu.add_separator()
            self.cluster_menu.add_command(label="Make Cluster Lead", command=self.make_lead_from_left_menu)
        self.popup_menu_left.add_cascade(menu=self.cluster_menu, label="Cluster Options")


    # -----
    # Selection menu.
    # -----
    def build_selection_menu(self):
        """
        Submenu with optios for managing a selection.
        """
        # Refreshes the menu each time the user clicks with the left mouse button on some sequence.
        self.selection_menu.delete(0,500)

        # Build the "Sequence" menu.
        self.build_selection_sequence_menu()
        self.selection_menu.add_separator()

        # Build the "Color" menu.
        self.build_selection_color_menu()

        # Build the "Structure" menu.
        if pymod.all_sequences_have_structure() or pymod.all_selected_elements_have_fetchable_pdbs():
            self.selection_menu.add_separator()
            self.build_selection_structure_menu()

        # Build the "Cluster" menu.
        if pymod.all_sequences_are_children():
            self.selection_menu.add_separator()
            self.build_selection_cluster_menu()


    def build_selection_sequence_menu(self):
        self.selection_sequence_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_sequence_menu.add_command(label="Save Selection to File", command=self.save_selection_from_left_pane)
        self.selection_sequence_menu.add_command(label="Copy Selection to Clipboard", command=self.copy_selection)
        self.selection_sequence_menu.add_separator()
        self.selection_sequence_menu.add_command(label="Duplicate Selection",command=self.duplicate_selection)
        self.selection_sequence_menu.add_command(label="Delete Selection", command=self.delete_many_sequences)
        self.selection_menu.add_cascade(menu=self.selection_sequence_menu, label="Sequences")


    def build_selection_color_menu(self):
        self.selection_color_menu = self.build_multiple_color_menu(mode="selection")
        self.selection_menu.add_cascade(menu=self.selection_color_menu, label="Color")


    def build_multiple_color_menu(self, mode, cluster_target_menu=None):
        """
        Used to build the color menu of both Selection and cluster elements popup menus.
        """
        target_menu = None
        color_selection_mode = None
        color_selection_target = None
        sequences_list = None
        color_target_label = None

        if mode == "selection":
            target_menu = self.selection_menu
            color_selection_mode = "selection"
            color_selection_target = None
            sequences_list = pymod.get_selected_sequences()
            color_target_label = "Selection"
        elif mode == "cluster":
            target_menu = cluster_target_menu
            color_selection_mode = "multiple"
            color_selection_target = pymod.get_children(self.get_cluster())
            sequences_list = pymod.get_children(self.get_cluster())
            color_target_label = "Cluster"

        multiple_color_menu = Menu(target_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        multiple_regular_colors_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        for color in pmdt.regular_colours:
            multiple_regular_colors_menu.add_command(label=color, command = lambda c=color: pymod.color_selection(color_selection_mode, color_selection_target, "regular",c))
        multiple_color_menu.add_cascade(menu=multiple_regular_colors_menu, label="Color whole %s by" % (color_target_label))
        multiple_color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        multiple_residues_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        multiple_residues_colors_menu.add_command(label="Polarity",command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "residue"))
        multiple_color_menu.add_cascade(menu=multiple_residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        n_selected_seqs = len(sequences_list)
        n_structures = len([e for e in sequences_list if e.has_structure()])
        n_seq_with_predicted_sec_str = len([e for e in sequences_list if e.has_predicted_secondary_structure()])

        if n_structures > 0 or n_seq_with_predicted_sec_str > 0:
            multiple_color_menu.add_separator()
            multiple_sec_str_color_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            # Available when all the selected sequences have a 3D structure.
            if n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Observed", command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "secondary-observed"))
            # Available only if all the sequences have a predicted secondary structure.
            if n_seq_with_predicted_sec_str == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "secondary-predicted"))
            # Available if there is at least one element with a 3D structure or a secondary
            # structure prediction.
            if not n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Auto (Observer + Predicted)", command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "secondary-auto"))
            multiple_color_menu.add_cascade(menu=multiple_sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if not False in [e.can_be_colored_by_conservation() for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_conservation_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "campo-scores"))
            multiple_color_menu.add_cascade(menu=multiple_conservation_colors_menu, label="By Convservation")

        # Energy colors.
        if not False in [e.can_be_colored_by_energy() for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_energy_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_energy_colors_menu.add_command(label="DOPE scores",command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "dope"))
            multiple_color_menu.add_cascade(menu=multiple_energy_colors_menu, label="By Energy")

        return multiple_color_menu


    def build_selection_structure_menu(self):
        self.selection_structure_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if pymod.all_sequences_have_structure():
            self.selection_structure_menu.add_command(label="Show chains in PyMOL", command=self.show_selected_chains_in_pymol)
            self.selection_structure_menu.add_command(label="Hide chains in PyMOL", command=self.hide_selected_chains_in_pymol)
            self.selection_structure_menu.add_separator()
            self.selection_structure_menu.add_command(label="Remove 3D Structures")
        elif pymod.all_selected_elements_have_fetchable_pdbs():
            self.selection_structure_menu.add_command(label="Fetch PDB Files", command=lambda: pymod.fetch_pdb_files("selection", None))
        self.selection_menu.add_cascade(menu=self.selection_structure_menu, label="Structures")


    def build_selection_cluster_menu(self):
        self.selection_cluster_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_cluster_menu.add_command(label="Extract Sequences from their Clusters", command=self.extract_selection_from_cluster)
        selected_sequences = pymod.get_selected_sequences()
        mother_indices_set = set([e.mother_index for e in selected_sequences])
        if len(mother_indices_set) == 1:
            mother = pymod.get_mother_by_index(list(mother_indices_set)[0])
            children = pymod.get_children(mother)
            if len(selected_sequences) < len(children):
                self.selection_cluster_menu.add_command(label="Extract Sequences to New Cluster", command=self.extract_selection_to_new_cluster_from_left_menu)
        self.selection_menu.add_cascade(menu=self.selection_cluster_menu, label="Cluster Options")


    # -----
    # Menu for cluster elements (alignments and similarity searches clusters).
    # -----
    def build_cluster_edit_menu(self, target_menu):
        self.cluster_edit_menu = Menu(target_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_edit_menu.add_command(label="Save Alignment To File", command=self.save_alignment_from_the_left_pan)
        self.cluster_edit_menu.add_separator()
        self.cluster_edit_menu.add_command(label="Transfer Alignment", command=self.transfer_alignment_from_the_left_pane)
        self.cluster_edit_menu.add_separator()
        self.cluster_edit_menu.add_command(label="Delete Cluster", command=self.delete_alignment_from_the_left_pane)
        target_menu.add_cascade(menu=self.cluster_edit_menu, label="Edit Cluster")


    def build_cluster_color_menu(self, target_menu):
        self.cluster_color_menu = self.build_multiple_color_menu(mode="cluster", cluster_target_menu = target_menu)
        target_menu.add_cascade(menu=self.cluster_color_menu, label="Color Cluster")


    ###################################
    # Sequence manipulation.          #
    ###################################

    # Extracts an element from an alignment.
    def extract_from_cluster(self):
        pymod.extract_child(self)
        pymod.gridder()

    def extract_selection_from_cluster(self):
        for e in pymod.get_selected_sequences():
            pymod.extract_child(e)
        pymod.gridder()

    def extract_selection_to_new_cluster_from_left_menu(self):
        pymod.extract_selection_to_new_cluster()

    def make_lead_from_left_menu(self):
        pymod.mark_as_lead(self)
        pymod.gridder()

    def save_sequence_from_left_pane(self):
        """
        Save option in the popup menu, it saves a single sequence.
        """
        pymod.sequence_save(self)

    def save_selection_from_left_pane(self):
        pymod.save_selection()

    # Copy option in the popup menu, copies a single sequence.
    def copy_sequence_to_clipboard(self):
        pymod.parent.clipboard_clear()
        pymod.parent.clipboard_append(self.my_sequence)# self.entry.get("1.0", END))

    # Copy selection
    def copy_selection(self):
        pymod.parent.clipboard_clear()
        text_to_copy = ""
        for element in pymod.pymod_elements_list:
            if element.selected and not element.is_cluster_element():
                # Adapt it for WINDOWS.
                text_to_copy += element.my_sequence + "\n"
        pymod.parent.clipboard_append(text_to_copy)


    def edit_sequence(self):
        child=Toplevel(pymod.main_window)
        child.resizable(0,0)
        #  self.child.geometry('400x500-10+40')
        child.title("<< Edit Sequence >>")
        child.config()
        try:
            child.grab_set()
        except:
            pass
        ch_main = Frame(child, background='black')
        ch_main.pack(expand = YES, fill = BOTH)

        midframe = Frame(ch_main, background='black')
        midframe.pack(side = TOP, fill = BOTH, anchor="n",
                              ipadx = 5, ipady = 5)

        lowerframe = Frame(ch_main, background='black')
        lowerframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center",
                              ipadx = 5, ipady = 5)

        L1 = Label(midframe,font = "comic 12", text="", bg="black", fg= "red")
        L1.grid( row=0, column=0, sticky="e", pady=5, padx=5)

        scrollbar = Scrollbar(midframe)
        scrollbar.grid(row=1, column=2, sticky="ns")

        textarea=Text(midframe, yscrollcommand=scrollbar.set, font = "comic 12",
                      height=10, bd=0, foreground = 'black', background = 'white',
                      selectbackground='black', selectforeground='white', width = 60 )
        textarea.config(state=NORMAL)
        textarea.tag_config("normal", foreground="black")
        textarea.insert(END, self.my_sequence)
        textarea.grid( row=1, column=1, sticky="nw", padx=0)

        scrollbar.config(command=textarea.yview)

        def submit():
            edited_sequence = textarea.get(1.0, "end").replace('\n','').replace('\r','').replace(' ','').replace('\t','').upper()
            if edited_sequence == "":
                pymod.show_error_message("Sequence Error", "Please submit a non empty string.", parent_window=child)
                return None
            if not pmsm.check_correct_sequence(edited_sequence):
                pymod.show_error_message("Sequence Error", "Please provide a sequence with only standard amino acid characters.", parent_window=child)
                return None
            self.my_sequence = edited_sequence
            pymod.gridder()
            child.destroy()

        sub_button=Button(lowerframe, text="SUBMIT", command=submit, relief="raised", borderwidth="3", bg="black", fg="white")
        sub_button.pack()


    def edit_structure(self):
        pass


    # Duplicates a single sequence.
    def duplicate_sequence_from_the_left_pane(self):
        pymod.duplicate_sequence(self)
        pymod.gridder()

    # Duplicate selection
    def duplicate_selection(self):
        for e in pymod.get_selected_sequences():
            pymod.duplicate_sequence(e)
        pymod.gridder()


    # Delete option in the popup menu. When multiple sequences have to be deleted the parameter
    # multiple is se to True.
    def delete_sequence_from_the_left_pane(self):
        pymod.delete_element(self)
        pymod.gridder()

    # Deletes many sequences.
    def delete_many_sequences(self):
        # First delete clusters that were entirely selected.
        to_delete_list = [e for e in pymod.pymod_elements_list if e.selected and e.is_cluster_element()]
        for element in to_delete_list:
            pymod.delete_whole_cluster(element)
        to_delete_list = [e for e in pymod.pymod_elements_list if e.selected]
        # Then delete other selected elements.
        for element in to_delete_list:
            pymod.delete_element(element)
        pymod.gridder()


    ###################################
    # Color sequences.                #
    ###################################

    def build_text_to_display(self):
        """
        This method displayes the sequence of an element by inserting it the ".sequence_entry" Text
        widget. It is called by "create_entry" method when "gridder" moethod of the PyMod class is
        called.
        """
        self.sequence_entry.tag_add("normal", "2.0")
        self.sequence_entry.insert(END, self.my_sequence,"normal")
        self.color_element(on_grid=True,color_pdb=False)
        self.sequence_entry.config(state=DISABLED)


    def color_element(self, on_grid=False, color_pdb=True):
        """
        Colors the sequence entry when it is displayed by the pymod.gridder() method or when the user
        changes the color scheme of a sequence. This can color also the PDB file of the element (if the
        element has one). In PyMod the PDB file is not colored when pymod.gridder() is called, it is
        colored only when the user decides to change the sequence color scheme.
        """
        # The 'residues_to_color' variable permits to color all the residues of a polypeptidic chain
        # (standard + modified) if set to 'all' or only standard residues if set to 'standard'.
        residues_to_color = "all"
        if self.is_shown:
            self.color_sequence_entry(on_grid, residues_to_color=residues_to_color)
            # Needed to display new colors on the color menu.
            if not on_grid:
                self.build_left_popup_menu()
        if color_pdb and self.has_structure():
            self.color_pdb_structure(residues_to_color=residues_to_color)


    def color_sequence_entry(self, on_grid, residues_to_color="all"):
        """
        Colors the sequence entry in PyMod main window.
        """
        # Colors the sequence according to some particular color scheme.
        if self.color_by != "regular":
            sequence = self.my_sequence
            for (i,aa) in enumerate(sequence):
                # Changing the foreground to "white" is needed to color indels with white.
                self.sequence_entry.tag_config("normal", foreground="white", background="black")
                # Creates a series of tags, used to color each residue differently.
                if aa != "-":
                    if residues_to_color == "standard" and aa == "X":
                        pass
                    else:
                        tag_name = aa+str(i)
                        self.sequence_entry.tag_add(tag_name, "1.%d" % (i)) # "1.0"
                        color = self.get_residue_color(residue_ids = (i,aa), color_target="sequence", residues_to_color=residues_to_color)
                        self.sequence_entry.tag_config(tag_name, foreground=color)
        # Just color each residue of the sequence with the same color.
        else:
            # First find if there are some tags in the Text (that is, if the sequence was colored
            # according to some particular color scheme that colors each residue differently).
            tags_to_delete = [tag for tag in self.sequence_entry.tag_names() if (tag != "normal" and tag != "sel")]
            # In order to color all the sequence with the same color the previously created tags
            # have to be deleted, so prooced to delete them.
            if tags_to_delete != []:
                for tag in tags_to_delete:
                    self.sequence_entry.tag_delete(tag)
            self.sequence_entry.tag_config("normal", foreground=self.my_color, background="black")


    def color_pdb_structure(self, residues_to_color="all"):
        """
        Colors the PDB structure of an element loaded into PyMOL.
        """
        chain_sel = self.build_chain_selector_for_pymol()
        # Colors the structure according to some particular color scheme.
        if self.color_by != "regular":
            # Gets the standard amminoacidic residues of the PDB file: they should be the same
            # residues of the sequence displayed in PyMod main window.
            residues = None
            if residues_to_color == "all":
                residues = filter(lambda r: r.residue_type == "standard" or r.hetres_type == "modified-residue", self.structure.pdb_chain_sequence)
            elif residues_to_color == "standard":
                residues = filter(lambda r: r.residue_type == "standard", self.structure.pdb_chain_sequence)
                # When the standard residues will be colored, this will make the modified residues
                # mantain a white color.
                cmd.color("white",chain_sel)
            for (i,res) in enumerate(residues):
                # Gets the right color for the current residue.
                color = self.get_residue_color(residue_ids = (i,str(res)), color_target="structure", residues_to_color=residues_to_color)
                # Builds a PyMOL selector for the current residue.
                residue_sel = self.build_residue_selector_for_pymol(res.pdb_position)
                # Finally colors the residue in PyMOL.
                cmd.color(color, residue_sel)
        # Colors all the residues of a structure with the same color.
        else:
            cmd.color(self.my_color,chain_sel)
        # Colors by atom.
        if self.color_by == "regular":
            cmd.util.cnc(chain_sel)


    def get_residue_color(self, residue_ids, color_target, residues_to_color="all"):
        """
        Gets the color of a residue in the sequence according to the color scheme defined by the
        "color_by" attribute. "residue_ids" is a tuple containing in [0] the numerical id of the
        residue in the aligned sequence, and in [1] the character of the residue.
        """
        # Name of the color to return.
        color = None

        # Needed to prepare colors for Tkinter widgets.
        def convert_to_tkinter_rgb(rgb_tuple):
            rgb_tuple = [i*255 for i in rgb_tuple]
            return '#%02x%02x%02x' % tuple(rgb_tuple)

        # Index of the residue to be colored.
        residue_index = None
        if color_target == "structure":
            residue_index = residue_ids[0]
        elif color_target == "sequence":
            if residues_to_color == "all":
                residue_index = pmsm.get_residue_id_in_gapless_sequence(self.my_sequence, residue_ids[0])
            elif residues_to_color == "standard":
                residue_index = pmsm.get_residue_id_in_gapless_sequence(self.my_sequence, residue_ids[0])
                residue_index -= self.my_sequence[0:residue_ids[0]].count("X")

        # Colors the sequence by residues.
        if self.color_by == "residue":
            if residue_ids[1] != "-":
                if color_target == "sequence":
                    color = convert_to_tkinter_rgb(pmdt.polarity_color_dictionary.get(residue_ids[1]))
                elif color_target == "structure":
                    # Generate a name to recall the color stored in PyMOL.
                    color = pmdt.pymol_polarity_color_name + residue_ids[1]
            else:
                color = "white"

        # Colors according by secondary structure assigned by PyMOL.
        elif self.color_by == "secondary-observed":
            if residue_ids[1] != "-":
                try:
                    color = pmdt.sec_str_color_dict.get(self.pymol_dss_list[residue_index])
                except:
                    color = "white"
            else:
                color = "white"

        # Colors according by secondary structure predicted by PSIPRED.
        elif self.color_by == "secondary-predicted":
            if residue_ids[1] != "-":
                psipred_vals = self.psipred_elements_list[residue_index]
                if color_target == "sequence":
                    psipred_tuple = (psipred_vals["confidence"], psipred_vals["sec-str-element"])
                    rgb = pmdt.psipred_color_dict[psipred_tuple]
                    color = convert_to_tkinter_rgb(rgb)
                elif color_target == "structure":
                    color = "%s_%s_%s" % (pmdt.pymol_psipred_color_name, psipred_vals["confidence"], psipred_vals["sec-str-element"])
            else:
                color = "white"

        # Color by CAMPO scores.
        elif self.color_by == "campo-scores":
            if residue_ids[1] != "-":
                # This string will be compatible with Tkinter.
                if color_target == "sequence":
                    color = convert_to_tkinter_rgb(pmdt.campo_color_dictionary[self.campo_scores[residue_index]["interval"]])
                # The id of the color (for example "campo1") needed to color the residue in PyMOL.
                elif color_target == "structure":
                    color = "%s_%s" % (pmdt.pymol_campo_color_name, self.campo_scores[residue_index]["interval"])
            else:
                color = "white"

        # Color by DOPE values.
        elif self.color_by == "dope":
            if residue_ids[1] != "-":
                # This string will be compatible with Tkinter.
                dope_vals = self.dope_items[residue_index]
                if color_target == "sequence":
                    rgb = pmdt.dope_color_dict[dope_vals["interval"]]
                    color = convert_to_tkinter_rgb(rgb)
                elif color_target == "structure":
                    color = "%s_%s" % (pmdt.pymol_dope_color_name, dope_vals["interval"])
            else:
                color = "white"

        return color


    # -----
    # Assigns new colors to the sequences using the menus.
    # -----

    def color_element_by_regular_color(self,color=None):
        """
        Colors sequence by "regular" colors, that is, colors uniformly the sequence with some color.
        """
        self.color_by = "regular"
        if color != None:
            self.my_color = color
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_residue(self):
        """
        Colors by residue properties.
        """
        self.color_by = "residue"
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_obs_sec_str(self):
        """
        Color elements by their observed secondary structure.
        """
        self.color_by = "secondary-observed"
        # If PyMOL has not been already used to assign sec str to this sequence.
        if not self.has_assigned_secondary_structure():
            pymod.assign_secondary_structure(self)
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_pred_sec_str(self):
        """
        # Color element by its predicted secondary structure.
        """
        # If PSI-PRED has not been already used to predict sec str for this sequence.
        if not self.has_predicted_secondary_structure():
            if pymod.predict_secondary_structure(self):
                self.color_by = "secondary-predicted"
        else:
            self.color_by = "secondary-predicted"
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_campo_scores(self):
        """
        Color by CAMPO scores.
        """
        self.color_by = "campo-scores"
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_dope(self):
        """
        Color by DOPE scores.
        """
        self.color_by = "dope"
        self.color_element(on_grid=False,color_pdb=True)


    ###################################
    # PDB files.                      #
    ###################################

    def pdb_is_fetchable(self):
        try:
            if self.header_entry.get().split("|")[2]=="pdb" or self.header_entry.get().split("|")[4]=="pdb":
                return True
            else:
                return False
        except:
            return False


    ###################################
    # Cluster elements.               #
    ###################################

    def get_cluster(self):
        if self.is_mother and self.is_cluster_element():
            return self
        elif self.is_child:
            return pymod.get_mother(self)


    def save_alignment_from_the_left_pan(self):
        pymod.alignment_save(self.get_cluster())


    def delete_alignment_from_the_left_pane(self):
        title = "Delete Cluster?"
        message = "Are you sure you want to delete %s?" % (self.get_cluster().my_header)
        choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)
        if choice:
            pymod.delete_alignment(self.get_cluster())
        pymod.gridder()


    def transfer_alignment_from_the_left_pane(self):
        pymod.transfer_alignment(self.get_cluster())


    ###############################################################################################
    # METHODS USED WHEN INTERACTING WITH THE SEQUENCE ENTRY (IN THE RIGHT PANE).                  #
    ###############################################################################################

    #################################################################
    # Mouse events and their bindings for the sequence Entry.       #
    #################################################################

    def bind_events_to_sequence_entry(self):
        self.sequence_entry.bind("<Leave>", self.leave_entry)
        self.sequence_entry.bind("<Motion>", self.set_messagebar_info)
        self.sequence_entry.bind("<Button-1>", self.on_sequence_left_click)
        self.sequence_entry.bind("<ButtonRelease-1>", self.on_sequence_left_release)
        self.sequence_entry.bind("<Enter>", self.enter_entry)
        # Centers and selects the residue in PyMOL by clicking on it with the middle mouse button.
        if self.has_structure():
            self.sequence_entry.bind("<Button-2>", self.click_residue_with_middle_button)
        self.sequence_entry.bind("<ButtonRelease-3>", self.on_sequence_right_click)

    def leave_entry(self,event):
            self.sequence_entry.unbind("<B1-Motion>")

    def set_messagebar_info(self,event):
        """
        Allows to show the protein name and the position of the aa residues in the bottom frames
        'pymod.sequence_name_bar', 'pymod.residue_bar'.
        """
        # Residues position (id +1) in the sequence.
        sequence_position = self.get_highlighted_residue_position()
        current_residue = self.sequence_entry.get(CURRENT)
        is_residue = None
        if current_residue in pmdt.protein_residues_set:
            is_residue = True
        else:
            is_residue = False

        # Include the position in sequences (and the PDB position for structures).
        position_text = ""
        if is_residue:
            if self.has_structure():
                sequence_index = self.get_highlighted_residue_position()
                residue_identifier = self.structure.pdb_chain_sequence[sequence_index-1].three_letter_code
                pdb_index = self.get_highlighted_residue_position(pdb_position=True)
                position_text = residue_identifier + " " + str(pdb_index) # + " (" + str(sequence_position) + ")"
            else:
                residue_identifier = pymod.one2three(current_residue)
                position_text = residue_identifier + " " + str(sequence_position)

        # Also include the position for alignments.
        if self.is_child:
            if is_residue:
                position_text += " - "
            position_text += ("Alignment: " + str(self.get_highlighted_residue_position(res_alignment_id=True)) )

        # Get additional information (depending on the sequence current color scheme) to show in the
        # message bar.
        if is_residue:
            if self.color_by == "campo-scores":
                score = self.campo_scores[sequence_position - 1]["campo-score"]
                position_text += " - CAMPO score: %s" % (score)
            elif self.color_by == "secondary-predicted":
                prediction = self.psipred_elements_list[sequence_position - 1]
                pred_text = str(prediction["confidence"]) + " " +pmdt.psipred_element_dict[prediction["sec-str-element"]]
                position_text += " - PSIPRED confidence: %s" % (pred_text)
            elif self.color_by == "dope":
                score = self.dope_items[sequence_position - 1]["dope-score"]
                position_text += " - DOPE score: %s" % (score)

        pymod.residue_bar.helpmessage(position_text)

        # Sequence name bar.
        protein_name = self.full_original_header # self.header_entry.get()
        pymod.sequence_name_bar.helpmessage(protein_name)


    #######################################
    # Methods needed to drag sequences    #
    # and to add/remove gaps to them.     #
    #######################################

    # Stores the X position of an aa when click (useful to calculate the shifting of a sequence
    # when dragging).
    def on_sequence_left_click(self,event):
        self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
        self.sequence_entry.config(state=NORMAL)
        # Sets state to 'NORMAL', so that the sequences can be modified with indels.
        if self.is_child and not self.is_lead_of_collapsed_cluster():
            for sibling in pymod.get_siblings(self):
                sibling.sequence_entry.config(state=NORMAL)

    def on_sequence_left_release(self,event):
        # Sets state to 'DISABLED', so that the sequences can't be modified with keyborad input
        # from the user.
        self.sequence_entry.config(state=DISABLED)
        if self.is_child and not self.is_lead_of_collapsed_cluster():
            for sibling in pymod.get_siblings(self):
                sibling.sequence_entry.config(state=DISABLED)

    def enter_entry(self,event):
        if not self.is_cluster_element():
            self.sequence_entry.bind("<B1-Motion>", self.on_motion)

    # Allows to insert/remove gaps '-' dragging the mouse
    def on_motion(self,event):

        # self.sequence_entry.config(state=NORMAL)

        drag = None

        # If dragging to the right insert an indel '-'.
        if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) > int(self.mypos.split(".")[1]):
            # Change the sequence itself.
            self.sequence_entry.insert(self.mypos, "-",("normal"))
            self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
            # Updates the sequence with new indels.
            # self.sequence_entry.config(width=int(self.sequence_entry['width'])+1)
            # self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END) # This fixes a bug on Ubuntu 14.04.
            self.update_sequence_from_entry()
            drag = "right"

        # If dragging to the left remove the gap '-' (if it exists).
        if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) < int(self.mypos.split(".")[1]) :
            if self.sequence_entry.get(self.sequence_entry.index("@%d,%d" % (event.x, event.y))) == "-":
                self.sequence_entry.delete("%s" % ("@%d,%d" % (event.x, event.y)))
                self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
                self.update_sequence_from_entry()
                drag = "left"

        # self.sequence_entry.config(state=DISABLED)

        # If the sequence is a child, the length of its siblings has to be adjusted and the sequence
        # update the of the mother has to be adjusted.
        if self.is_child and not self.is_lead_of_collapsed_cluster() and drag != None:

            #######################################################################################
            # NOTE:The optimal way to do this would be to rstrip all the sequences, then to ljust #
            # them to the lenght of the "longest" one. However Tkinter is too slow to do this, it #
            # takes too much time to update all the sequences in big clusters at the same time,   #
            # so as long as Tkinter is used the following code has to be applied. This code       #
            # prevents every sequence of a cluster from being updated every time an indel is      #
            # added, and it tries to update only the necessary sequences.                         #
            #######################################################################################

            # Gets the other elements in the cluster.
            mother = pymod.get_mother(self)
            children = pymod.get_children(mother)
            siblings = pymod.get_siblings(self)

            if drag == "right":
                # Removes extra gaps from the sequence being modified.
                self.rstrip_entry()
                rstripped_length = self.get_sequence_entry_length()
                maxlength = self.get_cluster_max_length(children)

                # If after dragging it the rstripped sequence is shorter than the others, adds extra
                # indels to it.
                if rstripped_length < maxlength:
                    self.ljust_entry(maxlength)
                # If the rstripped sequence is longer, adds extra gaps to other sequences to adjust
                # them to the same length.
                else:
                    for s in siblings:
                         s.ljust_entry(rstripped_length)

            elif drag == "left":
                # Removes extra gaps from the sequence being modified.
                self.rstrip_entry()
                rstripped_length = self.get_sequence_entry_length()
                maxlength = self.get_cluster_max_length(children)

                # This happens when, after removing an indel, the rstripped sequence is shorter than
                # the other sequences by just one character. For example
                #
                # before dragging:
                # -AAA-
                # -CCCC <- sequence being dragged
                # -DDD-
                #
                # after dragging:
                # -AAA-
                # CCCC  <- now it's shorter than one character from the maxlength of the cluster
                # -DDD-
                if rstripped_length + 1 == maxlength:
                    # If there are only indels as the last characters in other sequences of the
                    # cluster (such as in the example above) remove them.
                    only_indels = True
                    for s in siblings:
                        if s.get_sequence_entry_last_character() != "-":
                            only_indels = False
                            break
                    if only_indels:
                        for s in siblings:
                            if s.get_sequence_entry_last_character() == "-":
                                s.remove_sequence_entry_last_character()

                    # Adjust the dragged sequence with indels if necessary.
                    maxlength = self.get_cluster_max_length(children)
                    if rstripped_length != maxlength:
                        self.ljust_entry(maxlength)

                # Adjust the dragged sequence with indels.
                else:
                    self.ljust_entry(maxlength)

            # Then updates the mother.
            mother.sequence_entry.config(state=NORMAL)
            pymod.update_stars(mother)
            mother.sequence_entry.delete(1.0,END)
            mother.sequence_entry.insert(1.0, mother.my_sequence,("normal"))
            mother.sequence_entry.config(width=maxlength)
            mother.sequence_entry.config(state=DISABLED)


    # Takes as input a list of children elements and returns as an int the length of the one with
    # the longest entry.
    def get_cluster_max_length(self, children):
        return max([c.get_sequence_entry_length() for c in children])


    def update_sequence_from_entry(self):
        self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END)
        length = self.get_sequence_entry_length()
        self.sequence_entry.config(width=int(length))

    def get_sequence_entry_length(self):
        return len(self.sequence_entry.get("1.0", "%s-1c" % END))
        # return int(self.sequence_entry['width'])

    def get_sequence_entry_last_character(self):
        return self.sequence_entry.get("%s-2c" % END)

    def remove_sequence_entry_last_character(self, update=True):
        self.sequence_entry.delete("%s-2c" % END)
        if update:
            self.update_sequence_from_entry()

    def rstrip_entry(self,maxlength=None,update=True):
        # c.my_sequence = c.my_sequence.rstrip("-")
        found_residue = False
        while not found_residue:
            if maxlength != None and self.get_sequence_entry_length() <= maxlength:
                break
            if self.get_sequence_entry_last_character() == "-":
                self.remove_sequence_entry_last_character(update)
            else:
                found_residue = True

    def ljust_entry(self,maxlength,update=True):
        seql = self.get_sequence_entry_length()
        self.sequence_entry.insert("%s-1c" % END,"-"*(maxlength-seql))
        if update:
            self.update_sequence_from_entry()


    #######################################
    # Other methods needed to interact    #
    # with the sequences loaded into the  #
    # main window.                        #
    #######################################

    def click_residue_with_middle_button(self,event):
        if not self.is_current_position_indel():
            self.center_residue_in_pymol()
            self.select_residue_in_pymol()

    # A popup menu in the right frame to interact with the sequence
    def on_sequence_right_click(self,event):
        if not self.is_current_position_indel():
            try:
                self.popup_menu_right.tk_popup(event.x_root, event.y_root, 0)
            except:
                pass
            #popup_menu2.grab_release()

    ######################################
    # Some methods called in the above   #
    # events.                            #
    ######################################

    def is_current_position_indel(self):
        """
        Check if the current hilighted residue is an indel.
        """
        if self.sequence_entry.get(CURRENT) == "-":
            return True
        else:
            return False


    def get_highlighted_residue_position(self, res_alignment_id=False, pdb_position=False):
        """
        Returns the position of the residue currently highlighted with the mouse in the PyMod main
        window. If 'res_alignment_id' is set to 'True' it will return the id of that residue in the
        aligned sequence. If 'pdb_position' is 'True', it will return the id of that residue in
        the PDB file.
        """
        if res_alignment_id == True and pdb_position == True:
            raise Exception("This is not correct...")
        pos = int(self.sequence_entry.index(CURRENT).split(".")[1]) + 1
        if not res_alignment_id:
            number_gaps = self.sequence_entry.get("1.0", CURRENT).count('-')
            pos -= number_gaps
        if pdb_position:
            pos = self.structure.pdb_chain_sequence[pos -1].pdb_position
        return pos


    def get_pdb_index(self, position_index):
        number_gaps = self.my_sequence[0:position_index].count('-')
        position_index -= number_gaps
        return self.structure.pdb_chain_sequence[position_index].pdb_position


    #################################################################
    # Structure of the right pane popup menu.                       #
    #################################################################

    def build_right_popup_menu(self):
        """
        Builds the popup menu that appears when the user clicks with the left button on the
        sequence in the right pan.
        """
        # Right menu object.
        self.popup_menu_right = Menu(pymod.parent, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.element_type == "primary":
            pass
        elif self.element_type == "structure" or self.element_type == "model":
            self.popup_menu_right.add_command(label="Select Residue in PyMOL", command=self.select_residue_in_pymol)
            self.popup_menu_right.add_command(label="Center Residue in PyMOL", command=self.center_residue_in_pymol)


    attributes_to_store = ["unique_index",
                           "mother_index",
                           "child_index",
                           "is_mother",
                           "is_child",
                           "grid_position",
                           "my_header_fix",
                           "my_header",
                           "my_sequence",
                           "compact_header",
                           "full_original_header",
                           "pir_alignment_string",
                           "annotations",
                           "selected",
                           "is_shown",
                           "show_children",
                           "button_state",
                           "mother_button_color",
                           "color_by",
                           "my_color",
                           "pymol_dss_list",
                           "psipred_elements_list",
                           "campo_scores",
                           "dope_scores",
                           "dope_items",
                           "element_type",
                           "is_blast_query",
                           "is_lead",
                           "is_bridge",
                           "structure",
                           "alignment",
                           "is_model",
                           "models_count",
                           "polymer_type"]
    attributes_to_pickle = []
    lists_to_pickle = {}

    def __getstate__(self):
        return pymod_pickle(self)


###################################################################################################
# STORE PROJECTS' INFORMATION IN PICKLE FILES.                                                    #
###################################################################################################

class PyMod_Pickler(pickle.Pickler):
    pass


class PyMod_Unpickler(pickle.Unpickler):

    def find_class(self, module, name):
        """
        Overridden from the original 'Unpickler' class. Needed to rebuild PyMod object which have
        complex modules names. 'Unpickler' rebuilds objects using the 'fully qualified' name
        reference of their classes (the class name is pickled, along with the name of the module the
        class is defined in). Since PyMOL plugin modules may yield different 'fully qualified' names
        depending on the system, PyMod objects are rebuilt using only the name of their classes.
        """
        try:
            # Try the standard routine of pickle.
            __import__(module)
            mod = sys.modules[module]
            klass = getattr(mod, name)
            return klass
        except:
            # Build object by class name.
            try:
                name = name.rstrip("\n\r") # Needed to fix some old Windows versions behaviour.
            except:
                pass
            return globals()[name]


def pymod_pickle(obj):
    attributes_to_store = obj.attributes_to_store
    lists_to_pickle = obj.lists_to_pickle
    attributes_to_pickle = obj.attributes_to_pickle
    d = dict(obj.__dict__)
    for k in d.keys():
        if k in attributes_to_store:
            pass
        elif k in lists_to_pickle.keys():
            method_name = lists_to_pickle[k][0]
            d[k] = [getattr(list_item, method_name)() for list_item in d[k]]
        else:
            del d[k]
    return d
