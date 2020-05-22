# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing a class to represent the PyMod plugin.
"""

# Python standard library.
import os
import sys
import shutil
import re
import json
import datetime


# PyMOL.
from pymol import cmd

# PyMod modules.
from pymod_lib import pymod_os_specific as pmos # Different OS compatibility-related code.
from pymod_lib import pymod_vars as pmdt # PyMod data used throughout the plugin.
from pymod_lib import pymod_tool as pm_tool # Classes to represent tools used within PyMod.
from pymod_lib import pymod_element as pmel # Classes to represent sequences and alignments in PyMod.

# Path names used shared in various parts of the plugin.
from pymod_lib.pymod_vars import blast_databases_dirname, hmmer_databases_dirname, hmmscan_databases_dirname

# Part of the graphical user interface of PyMod.
from pymod_lib.pymod_gui.main_window.main_window_qt import PyMod_main_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt
from pymod_lib.pymod_gui.specific_gui_components_qt import PyMod_dir_selection_dialog

# Import modules with classes used to extend the 'PyMod' base class.
from ._development import PyMod_development
from ._main_menu_commands import PyMod_main_menu_commands
from ._workspaces import PyMod_workspaces
from ._external import PyMod_external
from ._files_managment import PyMod_files_managment
from ._elements_interactions import PyMod_elements_interactions
from ._elements_loading import PyMod_elements_loading
from ._pymol_interactions import PyMod_pymol_interactions
from ._installer import PyMod_installer


# Function that launches PyMod from the plugin menu of PyMOL.
def pymod_launcher(app, pymod_plugin_name, pymod_version, pymod_revision, parallel_modeller):
    pymod = PyMod(app, pymod_plugin_name, pymod_version, pymod_revision, parallel_modeller)


###################################################################################################
# Class that creates the plugin structure.                                                        #
###################################################################################################

class PyMod(PyMod_main_menu_commands,
            PyMod_workspaces,
            PyMod_external,
            PyMod_files_managment,
            PyMod_elements_interactions,
            PyMod_pymol_interactions,
            PyMod_elements_loading,
            PyMod_installer,
            PyMod_development):
    """
    Class to represent the PyMod plugin.
    """

    def __init__(self, app, pymod_plugin_name, pymod_version, pymod_revision, parallel_modeller):
        """
        Startup of the plugin.
        """

        #---------------------
        # Initializes PyMod. -
        #---------------------

        self.pymod_plugin_name = pymod_plugin_name
        self.pymod_version = pymod_version
        self.pymod_revision = pymod_revision

        self.DEVELOP = False
        self.TEST = False

        self.app = app

        # Initialize PyMod elements and information about the analyses performed by the user.
        self.initialize_pymod_elements_information()


        # --------------------------------------------------------------------------------------
        # Prepare PyMod files and folders that will be created in the project main directory. -
        # --------------------------------------------------------------------------------------

        # PyMod directory. The main folder where all PyMod files (with the exception of the
        # configuration file) will be stored.
        self.pymod_directory_name = "pymod"
        self.pymod_temp_directory_name = "pymod_temp_directory"
        self.projects_directory_name = "projects"
        self.new_project_default_name = "new_pymod_project"
        self.external_tools_dirname = "external_tools"
        self.data_dirname = "data"
        self.temp_directory_name = "temp_dir"
        self.active_session_filename = ".active_session.txt"
        self.inactive_session_filename = ".inactive_session.txt"


        # Structures.
        self.structures_dirname = "structures"  # structures_dirpath

        # Alignments.
        self.alignments_dirname = "alignments"  # alignments_directory
        self.alignments_files_names = "alignment"

        # Images directory
        self.images_dirname = "images"  # images_directory

        # Models.
        self.models_dirname = "models"  # models_directory
        self.models_subdirectory = "modeling_session"

        # PSIPRED.
        self.psipred_dirname = "psipred"

        # BLAST.
        self.similarity_searches_dirname = "similarity_searches"
        self.use_blast_v5_databases = True

        # Domain analysis.
        self.domain_analysis_dirname = 'domain_analysis'

        # File names used in order to store the information of PyMod sessions.
        self.session_temp_dirname = "pymod_project_temp_dir"
        self.session_filename = "pymod_session"

        # Initializes the main paths for the plugin.
        self.initializes_main_paths()


        # -----------------------
        # Prepare PyMod tools. -
        # -----------------------

        self.pymod_tools = []

        # PyMod itself.
        self.pymod_plugin = pm_tool.Tool("pymod", self.pymod_plugin_name)
        self.pymod_plugin.initialize_parameters([pm_tool.Tool_directory("pymod_dir_path", "PyMod Directory", editable=False)])
        self.pymod_tools.append(self.pymod_plugin)

        # ClustalW.
        self.clustalw = pm_tool.Executable_tool("clustalw", "Clustal W",)
        self.clustalw.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File", auto_find=True)])
        self.pymod_tools.append(self.clustalw)

        # Clustal Omega.
        self.clustalo = pm_tool.Executable_tool("clustalo", "Clustal Omega")
        self.clustalo.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File", auto_find=True)])
        self.pymod_tools.append(self.clustalo)

        # MUSCLE.
        self.muscle = pm_tool.Executable_tool("muscle", "MUSCLE")
        self.muscle.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File", auto_find=True)])
        self.pymod_tools.append(self.muscle)

        # BLAST+ suite. Used to run PSI-BLAST and store BLAST sequence databases retrieved from
        # ftp://ftp.ncbi.nlm.nih.gov/blast/db/ .
        self.blast_plus = pm_tool.Executable_tool("blast_plus", "BLAST+ suite")
        self.blast_plus.initialize_parameters([pm_tool.Tool_exec_directory("exe_dir_path", "Executables Directory", auto_find=True),
                                               # A default directory where the database folders available for the
                                               # PSI-BLAST database selection are going to be located.
                                               pm_tool.Tool_directory("database_dir_path", "Database Directory")])
        self.pymod_tools.append(self.blast_plus)

        # HMMER suite.
        self.hmmer_tool = pm_tool.Executable_tool("hmmer", "HMMER suite")
        self.hmmer_tool.initialize_parameters([pm_tool.Tool_exec_directory("exe_dir_path", "Executables Directory", auto_find=True),
                                               pm_tool.Tool_directory("database_dir_path", "HMMER Database Directory"),
                                               pm_tool.Tool_directory("hmmscan_db_dir_path", "HMMSCAN Database Directory")])
        self.pymod_tools.append(self.hmmer_tool)

        # PSIPRED.
        self.psipred = pm_tool.Executable_tool("psipred", "PSIPRED")
        self.psipred.initialize_parameters([pm_tool.Tool_directory("exe_dir_path", "Executables Directory"),
                                            pm_tool.Tool_directory("data_dir_path", "Data Files Directory"),
                                            pm_tool.Tool_directory("database_dir_path", "BLAST Database Directory")])
        self.pymod_tools.append(self.psipred)

        # KSDSSP.
        self.ksdssp = pm_tool.Executable_tool("ksdssp", "KSDSSP")
        self.ksdssp.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.ksdssp)

        # MODELLER.
        self.modeller = pm_tool.Modeller_tool("modeller", "MODELLER")
        self.modeller.initialize_parameters([pm_tool.Use_importable_modeller("use_importable_modeller", "Internal MODELLER"),])
        self.pymod_tools.append(self.modeller)
        self.parallel_modeller = parallel_modeller

        # Initializes tools.
        for tool in self.pymod_tools:
            tool.pymod = self

        # If set to 'True' the most time consuming protocols of PyMod will be run in a thread so
        # that the GUI is not freezed. When developing the code, it is better to set it to 'False',
        # in order to better track exceptions.
        if self.DEVELOP:
            self.use_protocol_threads = False
        else:
            self.use_protocol_threads = True

        # Initializes colors for PyMod and PyMOL.
        self.initialize_colors()


        #----------------------------------------
        # Builds the main window of the plugin. -
        #----------------------------------------

        self.main_window = PyMod_main_window_qt(self)


        #-----------------------
        # Starts up a new job. -
        #-----------------------

        # If it is not found, then treat this session as the first one and asks the user to input
        # the 'PyMod Directory' path before beginning the first PyMod job.
        if not os.path.isfile(self.cfg_file_path):
            self.show_first_time_usage_message()
            self.show_pymod_directory_selection_window()

        # The configuration file is found.
        else:

            # Start an usual PyMod session. Get values options for each PyMod tool and start a
            # new PyMod job.
            try:
                self.get_parameters_from_configuration_file()
                # Checks if the PyMod directory exists.
                if not os.path.isdir(self.pymod_plugin["pymod_dir_path"].get_value()):
                    message = ("The project directory specified in PyMod configuration"
                               " file ('%s') is missing. Please specify a new"
                               " one." % self.pymod_plugin["pymod_dir_path"].get_value())
                    raise Exception(message)
                self.new_job_state()

            except Exception as e:
                if self.DEVELOP:
                    raise e
                self.show_configuration_file_error(e, "read")
                title = 'Configuration file repair'
                message = "Would you like to delete PyMod configuration file and build a new functional copy of it?"
                repair_choice = askyesno_qt(title, message, self.get_qt_parent())
                if repair_choice:
                    self.show_pymod_directory_selection_window()
                else:
                    self.main_window.destroy()


    def show_first_time_usage_message(self):
        title = "PyMod first session"
        message = "This is the first time you run PyMod. Please specify a folder inside which to build the 'PyMod Directory'. "
        message += "All PyMod files (such as its external tools executables, sequence databases and its project files) will be stored in this 'PyMod Directory' on your system."
        self.main_window.show_info_message(title, message)


    def get_pymod_app(self):
        return self.app.root

    def get_qt_parent(self):
        return None


    ###############################################################################################
    # Configuration file and first time usage.                                                    #
    ###############################################################################################

    def get_parameters_from_configuration_file(self):
        """
        Updates the values of the PyMod Tools parameters according to the information in the main
        configuration file.
        """
        cfgfile = open(self.cfg_file_path, 'r')
        # Reads the json configuration file, where PyMod options are stored in a dictionary.
        pymod_config = json.loads(cfgfile.read())
        for tool_object in self.pymod_tools:
            tool_dict = pymod_config[tool_object.name]
            for parameter_name in list(tool_dict.keys()):
                tool_object[parameter_name].value = tool_dict[parameter_name]
        cfgfile.close()


    def show_pymod_directory_selection_window(self):
        """
        Allows to select the 'PyMod Directory' on PyMod first run.
        """
        self.pymod_dir_window = PyMod_dir_selection_dialog(app=self.main_window, pymod=self)
        self.pymod_dir_window.exec_()


    def pymod_directory_selection_state(self):
        """
        This is called when the SUBMIT button on the "PyMod project" window is pressed.
        """
        try:

            # Check if the parent folder of the new PyMod directory exists.
            new_pymod_directory_parent = self.pymod_dir_window.main_entry.text()
            if not os.path.isdir(new_pymod_directory_parent):
                title = 'PyMod directory Error'
                message = "The directory inside which you would like to create your 'PyMod Directory' ('%s') does not exist on your system. Please select an existing path." % (new_pymod_directory_parent)
                self.main_window.show_error_message(title, message)
                return None


            # Check if a PyMod directory already exists in the parent folder.
            new_pymod_dirpath = os.path.join(new_pymod_directory_parent, self.pymod_directory_name)
            if os.path.exists(new_pymod_dirpath):
                title = 'PyMod directory Warning'
                message = "A folder named '%s' already exists on your system. Would you like to overwrite it and all its contents to build a new 'PyMod Directory'?" % (new_pymod_dirpath)
                overwrite = askyesno_qt(title, message, parent=self.get_qt_parent())
                if overwrite:
                    if os.path.isfile(new_pymod_dirpath):
                        os.remove(new_pymod_dirpath)
                    if os.path.isdir(new_pymod_dirpath):
                        shutil.rmtree(new_pymod_dirpath)
                else:
                    return None

            # Check if the configuration file directory exist. If it does not exist, build it.
            if not os.path.exists(self.cfg_directory_path):
                os.mkdir(self.cfg_directory_path)

            # Builds the new PyMod directory with its projects folder.
            os.mkdir(new_pymod_dirpath)
            os.mkdir(os.path.join(new_pymod_dirpath, self.projects_directory_name))

            # Builds empty external tools and data directories.
            os.mkdir(os.path.join(new_pymod_dirpath, self.external_tools_dirname))
            os.mkdir(os.path.join(new_pymod_dirpath, self.data_dirname))

            # Writes the configuration file.
            cfgfile = open(self.cfg_file_path, 'w')
            pymod_config = {}

            # Initializes the parameters for each tool.
            for tool in self.pymod_tools:
                new_tool_parameters = {}
                for parameter in tool.parameters:
                    # new_tool_parameters.update({parameter.name: parameter.get_starting_value()})
                    new_tool_parameters.update({parameter.name: ""})
                pymod_config.update({tool.name: new_tool_parameters})
            pymod_config["pymod"]["pymod_dir_path"] = new_pymod_dirpath
            pymod_config["modeller"]["use_importable_modeller"] = True

            # Define the default database paths for the BLAST and HMMER suites. Users may later change
            # these through the Options window.
            for tool_name, param_name, dirname in [("blast_plus", "database_dir_path", blast_databases_dirname),
                                                   ("hmmer", "database_dir_path", hmmer_databases_dirname),
                                                   ("hmmer", "hmmscan_db_dir_path", hmmscan_databases_dirname)]:
                database_dirpath = os.path.join(new_pymod_dirpath, self.data_dirname, dirname)
                pymod_config[tool_name][param_name] = database_dirpath
                if not os.path.isdir(database_dirpath):
                    os.mkdir(database_dirpath)

            cfgfile.write(json.dumps(pymod_config))
            cfgfile.close()

        except Exception as e:
            if self.DEVELOP:
                raise e
            title = "PyMod Directory Error"
            message = "Unable to write the PyMod configuration directory '%s' because of the following error: %s." % (self.cfg_directory_path, e)
            self.main_window.show_error_message(title, message)
            return None

        # Begin a new PyMod job.
        self.pymod_dir_window.close()
        self.get_parameters_from_configuration_file()
        # tkMessageBox.showinfo("PyMod first run", "You are about to begin your first PyMod session.", parent=self.get_qt_parent())
        self.new_job_state()


    def show_configuration_file_error(self, error, mode):
        if mode == "read":
            action = "reading"
        elif mode == "write":
            action = "writing"
        else:
            action = None
        title = "Configuration file error"
        message = "There was an error while %s the PyMod configuration file (located at '%s'). The error is: '%s'." % (action, self.cfg_file_path, error)
        self.main_window.show_error_message(title, message)


    ###############################################################################################
    # Start a new job.                                                                            #
    ###############################################################################################

    def initialize_pymod_elements_information(self):
        """
        Initializes those attributes storing information on sequences and structures loaded in PyMod.
        Called when starting a new PyMod session.
        """

        # This is the list where are going to be stored all the sequences displayed in the main
        # window represented as objects of the "PyMod_element" class.
        # self.pymod_elements_list = []
        self.root_element = pmel.PyMod_root_element(header="PyMod root")

        # An index that increases by one every time an element is added to the above list by using
        # the .add_element_to_pymod() method.
        self.unique_index = 0
        self.new_objects_index = 0

        # A list of all PDB files loaded in PyMod by the user.
        self.pdb_list = []
        # Count of the PDB files parsed in PyMod.
        self.pdb_counter = 0

        # Alignments counters.
        self.alignment_count = 0
        self.new_clusters_counter = 0

        # Logo images counters.
        self.logo_image_counter = 0

        # Attributes that will keep track of how many models the user builds in a PyMod session.
        self.performed_modeling_count = 0
        # This will keep track of how many multiple chains models the user builds in a PyMod session.
        self.multiple_chain_models_count = 0
        # This will contain ojects of the 'Modeling_session' class in order to build the 'Models'
        # submenu on the plugin's main menu.
        self.modeling_session_list = []

        # BLAST analysis.
        self.blast_cluster_counter = 0

        # Domain analysis.
        self.active_domain_analysis_count = 0

        # This is an index that will bw used to color each structure loaded into PyMod with a
        # different color taken from the list above. It will be used in "color_struct()".
        self.color_index = 0
        self.custom_colors_index = 0 # Each time a user adds a custom color, this will be increased.

        # Attributes to store the "compact headers" of elements loaded in PyMod. They allow to
        # automatically change the name of new elements if there are already some elements with the
        # same name.
        self.compact_header_root_dict = {}
        self.compact_header_set = set()

        # Similar set for full PDB files.
        self.original_pdb_files_set = set()


    def initializes_main_paths(self):

        home_dirpath, cfg_directory_path, pymod_envs_dirpath, pymod_env_dirpath, pymod_pkgs_dirpath = pmos.get_pymod_cfg_paths()
        # Gets the home directory of the user.
        self.home_directory = home_dirpath
        self.current_project_name = None
        self.current_project_dirpath = None

        # Creates the preferences file in an hidden directory in the home directory.
        self.cfg_directory_path = cfg_directory_path
        self.cfg_file_name = "preferences.txt"
        self.cfg_file_path = os.path.join(self.cfg_directory_path, self.cfg_file_name)

        # PyMod conda environment. These directories will only be used if the conda root of PyMOL
        # is located in some directory without writing permissions. In this case, new conda
        # packages will be downloaded and installed in this directory.
        self.pymod_envs_dirpath = pymod_envs_dirpath
        self.pymod_env_dirpath = pymod_env_dirpath
        self.pymod_pkgs_dirpath = pymod_pkgs_dirpath


    # def show_new_job_window(self):
    #     """
    #     Builds a window that let users choose the name of the new projects direcotory at the
    #     beginning of a PyMod session.
    #     """
    #     # self.new_dir_window = New_project_window(self, self.main_window)
    #     self.new_job_state()

    def new_job_state(self, overwrite=False):
        """
        This is called when the SUBMIT button on the "New Job" window is pressed.
        """

        # Checks if the name is valid.
        new_dir_name = self.new_project_default_name # self.new_dir_window.main_entry.get()
        if bool(re.search('[^A-z0-9-_]', new_dir_name)) or "\\" in new_dir_name:
            self.main_window.show_error_message('Name Error', 'Only a-z, 0-9, "-" and "_" are allowed in the project name.')
            return None

        # Writes the directory.
        try:

            # Composes the name of the new projects directory.
            pymod_projects_dir_path = os.path.join(self.pymod_plugin["pymod_dir_path"].get_value(), self.projects_directory_name)
            if not os.path.isdir(pymod_projects_dir_path):
                os.mkdir(pymod_projects_dir_path)

            new_project_dir_path = os.path.join(pymod_projects_dir_path, new_dir_name)

            # Decide whether to build a new directory or overwrite the content
            # of the previous one.
            if os.path.isdir(new_project_dir_path):

                if overwrite:
                    # Just overwrite the previous directory.
                    active_session = False
                else:
                    # Look for a file signaling an active sessions in the project
                    # directory.
                    active_session = self.check_exisisting_session_tempfile(new_project_dir_path)

                if not active_session:
                    # Remove the content of the previous folder.
                    self.remove_previous_project_files(new_project_dir_path)

                else:
                    # Let the user decide whether to remove the content of the
                    # previous directory.
                    title = "Overwrite an existing project?"
                    message = ("An non-finished PyMod project was found in the"
                               " 'projects' folder of your PyMod Directory. This"
                               " either means that PyMod has shut down unexpectedly last"
                               " time or that another PyMod instance is already running"
                               " on your system.\nWould you like to overwrite the"
                               " existing project folder? If another PyMod instance"
                               " is running, its data will be lost (only one PyMod"
                               " instance can be run at any time). If you are not"
                               " running any other PyMod instances, you can safely"
                               " overwrite the old project.")
                    buttons_text = ("Overwrite and continue", "Skip and quit")

                    if self.DEVELOP or self.TEST:
                        overwrite_choice = True
                    else:
                        overwrite_choice = askyesno_qt(title, message, self.get_qt_parent(),
                                                       buttons_text=buttons_text)

                    if overwrite_choice:
                        # Remove the content of the old directory.
                        self.remove_previous_project_files(new_project_dir_path)
                    else:
                        # Quit PyMod.
                        self.main_window.close()
                        return None

            # Simply builds the new directory.
            else:
                os.mkdir(new_project_dir_path)

            # Start the new project.
            self.initialize_session(new_dir_name)

        except Exception as e:
            if self.DEVELOP:
                raise e
            message = "Unable to write directory '%s' because of the following error: %s." % (new_dir_name, e)
            self.main_window.show_error_message("Initialization error", message)
            try: # This raises an exception if the PyMod project has not been initialized.
                self.inactivate_session()
            except:
                pass
            self.main_window.close()
            return None


    def check_exisisting_session_tempfile(self, new_project_dir_path):
        return os.path.isfile(os.path.join(new_project_dir_path, self.active_session_filename))


    def initialize_session(self, new_project_name, saved_session=False):
        """
        Initializes a session and shows the main window, which was previously hidden.
        """
        self.current_pymod_dirpath = self.pymod_plugin["pymod_dir_path"].get_value()
        self.current_project_name = new_project_name
        self.current_projects_dirpath = os.path.join(self.current_pymod_dirpath, self.projects_directory_name)
        self.current_project_dirpath = os.path.join(self.current_projects_dirpath, self.current_project_name)
        self.structures_dirpath = os.path.join(self.current_project_dirpath, self.structures_dirname)
        self.alignments_dirpath = os.path.join(self.current_project_dirpath, self.alignments_dirname)
        self.images_dirpath = os.path.join(self.current_project_dirpath, self.images_dirname)
        self.models_dirpath = os.path.join(self.current_project_dirpath, self.models_dirname)
        self.psipred_dirpath = os.path.join(self.current_project_dirpath, self.psipred_dirname)
        self.similarity_searches_dirpath = os.path.join(self.current_project_dirpath, self.similarity_searches_dirname)
        self.temp_directory_dirpath = os.path.join(self.current_project_dirpath, self.temp_directory_name)
        self.domain_analysis_dirpath = os.path.join(self.current_project_dirpath, self.domain_analysis_dirname)

        # Writes a temporary file to signal that the session is active.
        with open(os.path.join(self.current_project_dirpath, self.active_session_filename), "w") as t_fh:
            date_obj = datetime.datetime.now()
            t_fh.write("Session started: %s\n" % (date_obj.strftime("%d/%m/%Y at %H:%M:%S")))

        # Creates the subdirectories.
        if not saved_session:
            self.create_project_subdirectories()
            self.launch_default()


    def launch_default(self):
        """
        Actions performed when initializing a new PyMod session.
        """
        self._launch_develop()


    def create_project_subdirectories(self):
        for single_dirpath in (self.alignments_dirpath, self.images_dirpath,
                               self.images_dirpath, self.models_dirpath,
                               self.structures_dirpath, self.psipred_dirpath,
                               self.similarity_searches_dirpath, self.temp_directory_dirpath,
                               self.domain_analysis_dirpath):
            try:
                os.mkdir(single_dirpath)
            except:
                pass


    def remove_previous_project_files(self, project_dirpath):
        """
        Removes the previously used files and subdirectories of a project directory
        when users decide to overwrite an existing project's directory.
        """
        # Safely remove all the files in the previously used directory.
        for the_file in os.listdir(project_dirpath):
            file_path = os.path.join(project_dirpath, the_file)
            if os.path.isfile(file_path):
                os.remove(file_path)

        # Remove the directories.
        for single_dir in (self.structures_dirname, self.alignments_dirname, self.models_dirname, self.psipred_dirname,
                           self.similarity_searches_dirname, self.images_dirname, self.temp_directory_name,
                           self.domain_analysis_dirname):
            dir_to_remove_path = os.path.join(project_dirpath, single_dir)
            if os.path.isdir(dir_to_remove_path):
                shutil.rmtree(dir_to_remove_path)


    def inactivate_session(self):
        active_session_filepath = os.path.join(self.current_project_dirpath, self.active_session_filename)
        if os.path.isfile(active_session_filepath):
            inactive_session_filepath = os.path.join(self.current_project_dirpath, self.inactive_session_filename)
            if os.path.isfile(inactive_session_filepath):
                os.remove(inactive_session_filepath)
            with open(active_session_filepath, "a") as t_fh:
                date_obj = datetime.datetime.now()
                t_fh.write("Session finished: %s\n" % (date_obj.strftime("%d/%m/%Y at %H:%M:%S")))
            shutil.move(active_session_filepath, inactive_session_filepath)


    ###############################################################################################
    # Colors.                                                                                     #
    ###############################################################################################

    def initialize_colors(self):
        """
        Prepares colors for PyMod and PyMOL.
        """

        # Stores the color used for PyMod elements.
        self.all_colors_dict_tkinter = {}
        # Import in PyMod some PyMOL colors.
        self.update_pymod_color_dict_with_dict(pmdt.pymol_regular_colors_dict_rgb, update_in_pymol=False)
        # Light PyMOL colors.
        self.update_pymod_color_dict_with_dict(pmdt.pymol_light_colors_dict_rgb)
        # PyMod colors.
        # self.update_pymod_color_dict_with_list(pmdt.pymod_regular_colors_list, update_in_pymol=False)
        self.all_colors_dict_tkinter.update(pmdt.pymod_regular_colors_dict)
        # Prepares other colors for PyMOL and PyMod.
        self.update_pymod_color_dict_with_dict(pmdt.sec_str_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.psipred_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.campo_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.scr_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.entropy_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.pc_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.dope_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.polarity_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.residue_type_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.domain_colors_dict)

    def update_pymod_color_dict_with_dict(self, color_dict, update_in_pymol=True):
        for color_name in list(color_dict.keys()):
            if update_in_pymol:
                cmd.set_color(color_name, color_dict[color_name])
            self.all_colors_dict_tkinter.update({color_name: pmdt.convert_rgb_to_hex(color_dict[color_name])})

    def update_pymod_color_dict_with_list(self, color_list, update_in_pymol=False):
        for color_name in color_list:
            if update_in_pymol:
                cmd.set_color(color_name, color_name)
            self.all_colors_dict_tkinter.update({color_name: color_name})

    def add_new_color(self, new_color_rgb, new_color_hex):
        """
        Adds a new color to PyMod/PyMOL.
        """
        new_color_name = "custom_color_%s" % self.custom_colors_index
        self.all_colors_dict_tkinter[new_color_name] = new_color_hex
        cmd.set_color(new_color_name, new_color_rgb)
        self.custom_colors_index += 1
        return new_color_name
