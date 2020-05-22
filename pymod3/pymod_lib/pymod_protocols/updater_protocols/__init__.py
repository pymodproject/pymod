# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for updating the BLAST and HMMER databases in PyMod by downloading them from the NCBI and
EBI servers.
"""

import os
import subprocess
import json
import time

from pymod_lib.pymod_os_specific import get_exe_file_name, get_formatted_date, check_network_connection
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_vars import data_installer_log_filename
from .updater_internal import all_components_list, Pfam_db_installer, BLAST_db_installer
from .updater_gui import InstallerUpdaterWindow


class UpdaterProtocol(PyMod_protocol):

    protocol_name = "database_updater"

    def __init__(self, pymod):
        PyMod_protocol.__init__(self, pymod)

        # selected components to be updated/installed
        self.selected_comps = []
        self.pending_downloads = 0
        self.pymod_data_dirpath = os.path.join(self.pymod.current_pymod_dirpath, self.pymod.data_dirname)
        # Path of the file where the download dates will be stored.
        self.download_log_filepath = os.path.join(self.pymod_data_dirpath, data_installer_log_filename)


    def launch_from_gui(self):

        # If there is no internet connection, everything is blocked.
        if not check_network_connection('https://www.google.it/'):
            self.pymod.main_window.show_error_message("Connection Error",
                "No internet connection. Can not update PyMod databases.")
            return None

        # Associates the correct installer to each component.
        self.components_dict_installer = {}
        self.components_list = []

        # Builds Qthread objects.
        _all_components_list = all_components_list[:]
        if self.pymod.use_blast_v5_databases:
            _all_components_list.pop(0) # Remove BLAST v4 databases.
        else:
            _all_components_list.pop(1) # Remove BLAST v5 databases.
        for comp, installer_class in zip(_all_components_list, [BLAST_db_installer, Pfam_db_installer, BLAST_db_installer]):
            component_qthread = installer_class(component=comp)
            self.components_dict_installer[comp.name] = component_qthread

            if comp.name == 'hmmscan_databases':
                comp.installer.hmmpress_exe_filepath = os.path.join(self.pymod.hmmer_tool["exe_dir_path"].get_value(),
                                                                    get_exe_file_name("hmmpress"))
            self.components_list.append(comp)
            comp.reinitialize()

        # Showing GUI.
        self.window = InstallerUpdaterWindow(parent=self.pymod.main_window, installer_protocol=self)

    #----------------
    # Installation. -
    #----------------

    def install_selected(self):
        """Called from the GUI slot that responds to the Install Selected button"""

        # Builds the data directory.
        if not os.path.isdir(self.pymod_data_dirpath):
            os.mkdir(self.pymod_data_dirpath)

        # Builds a temporary directory in PyMod configuration folder in /home/user
        # where to unzip the external tools and data components installation directories.
        # The files will later be moved in PyMod directory.
        self.build_temp_directory_in_cfg_dir()

        # Prepares to download the components in this temporary directory.
        # for inst in self.components_dict_installer.values():
        #     inst.download_destination_dirpath = self.temp_files_cfg_directory_path
        for comp in self.selected_comps:
            comp.installer.download_destination_dirpath = self.temp_files_cfg_directory_path

        self.pending_downloads = len(self.selected_comps)

        # Starts the installer QThread for every selected component in the GUI.
        for comp in self.selected_comps:
            comp.installer.install_mode = "download"
            time.sleep(0.1)
            comp.installer.start()

        self.window.activate_progressbar()


    def finish_installation(self):
        """
        Updates the download log and changes the GUI.
        """

        last_downloaded_dict = {}
        installed_count = 0
        for comp in self.selected_comps:
            try:
                print(comp)
                for element in comp.installer.downloaded_filepath_list:
                    print(element)
                    comp.installer.unzip_downloaded_files(comp.target_installation_path,
                                                          in_thread=False)
            except:
                continue
            if comp.installer.installation_status == "success":
                # Set the last downloaded time.
                last_downloaded_dict[comp.name] = comp.last_downloaded
                installed_count += 1

        # Writes the downloads log.
        if last_downloaded_dict:
            # Get the values from the old log, if present.
            if os.path.isfile(self.download_log_filepath):
                old_downloaded_dict = {}
                with open(self.download_log_filepath, "r") as l_fh:
                    old_downloaded_dict = json.loads(l_fh.read())
                # Updates with new values.
                for k in old_downloaded_dict:
                    if not k in last_downloaded_dict:
                        last_downloaded_dict[k] = old_downloaded_dict[k]
            # Actually writes the log.
            with open(self.download_log_filepath, "w") as l_fh:
                l_fh.write(json.dumps(last_downloaded_dict))

        # Close the database updater GUI.
        message = ("Database update terminated successfully for %s of %s"
                   " components." % (installed_count, len(self.selected_comps)))
        self.window.info_slot(message, None)
        # self.window.close()


    def build_temp_directory_in_cfg_dir(self):
        self.temp_files_cfg_directory_path = os.path.join(self.pymod.cfg_directory_path,
                                                          self.pymod.pymod_temp_directory_name)
        # self.temp_files_cfg_directory_path = self.pymod.temp_directory_dirpath
        if not os.path.exists(self.temp_files_cfg_directory_path):
            os.mkdir(self.temp_files_cfg_directory_path)


    def set_component_target_path(self, component):
        """
        Set installation paths inside the 'data' directory inside the PyMod directory. It will produce
        something like: /home/user/pymod/data/blast_databases.
        """

        if component.name == 'blast_databases':
            path_list = [self._get_database_dirpath("blast", "database_dir_path")]
        elif component.name == 'hmmer_databases':
            path_list = [self._get_database_dirpath("hmmer", "database_dir_path")]
        elif component.name == "hmmscan_databases":
            path_list = [self._get_database_dirpath("hmmer", "hmmscan_db_dir_path")]
        else:
            raise KeyError("Unknown component: %s" % (component.name))

        try:
            component.target_installation_path = os.path.join(*path_list)
        except TypeError: # Some path was undefined.
            component.target_installation_path = None

    def _get_database_dirpath(self, tool_name, param_name):
        if tool_name == "blast":
            param_value = self.pymod.blast_plus[param_name].get_value()
        elif tool_name == "hmmer":
            param_value = self.pymod.hmmer_tool[param_name].get_value()
        else:
            raise KeyError("Unknown tool: %s" % (tool_name))
        if param_value.replace(" ", "") == "":
            return None
        else:
            return param_value
