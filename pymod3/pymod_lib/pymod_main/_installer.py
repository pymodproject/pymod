# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for automatically downloading and installing the external PyMod tools (BLAST+, HMMER, Clustal,
MUSCLE, psipred, ksdssp) and the MODELLER conda package.
"""

import os
import sys
import shutil
import re
import subprocess
import zipfile
import urllib.request
import json
import importlib
import time

from pymol import cmd

from pymol.Qt import QtCore, QtWidgets

from pymod_lib.pymod_os_specific import (get_python_architecture, get_os_architecture, get_formatted_date,
                                         get_exe_file_name, get_home_dir,
                                         check_network_connection, check_importable_modeller)
from pymod_lib.pymod_vars import (blast_databases_dirname, hmmer_databases_dirname, hmmscan_databases_dirname,
                                  data_installer_log_filename)
from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt, askopenfile_qt
import warnings

# Installer parameters.
tools_installer_log_filename = "pymod_installer_log.txt"

allowed_platforms = ("linux", "darwin", "win32")

packages_url_dict = {"linux":  "https://github.com/pymodproject/pymod/releases/download/v3.0/linux_pymod_3.0_installer_bundle.zip",
                     "darwin": "https://github.com/pymodproject/pymod/releases/download/v3.0/macos_pymod_3.0_installer_bundle.zip",
                     "win32":  "https://github.com/pymodproject/pymod/releases/download/v3.0/windows_pymod_3.0_installer_bundle.zip"}
use_external_tools_local_install = False

python_minor_version = "%s.%s" % (sys.version_info.major, sys.version_info.minor)


def catch_errors_installer_threads(function):
    """
    Catches errors in the installer threads and closes the dialogs.
    """

    def wrapper(self, *args, **kwargs):
        try:
            return function(self, *args, **kwargs)
        except Exception as e:
            self.emit_exception.emit(e)

    return wrapper


class PyMod_installer:
    """
    Mixin class to extend the main PyMod class with an automatic installer for PyMod components.
    """

    def launch_components_installation(self, event=None):
        """
        Called when selecting the installer option from PyMod main menu.
        """

        #-----------------------------------------------
        # Checks if the installation can be performed. -
        #-----------------------------------------------

        # Check the internet connection, as the installation requires to download packages from the
        # internet.
        if not check_network_connection("https://google.com", timeout=3): # "http://216.58.192.142"
            has_network = False
            if not use_external_tools_local_install:
                message = ("An internet connection is not available. The automatic installation"
                           " of PyMod Components requires a connection to download their"
                           " packages. Please use an internet connection if you want to perform"
                           " the automatic installation.")
                self.main_window.show_error_message("Installation Error", message)
        else:
            has_network = True

        # Check if the user's system supports installation of components.
        if sys.platform not in allowed_platforms:
            message = ("Your operating system ('%s') is not supported by the automatic PyMod"
                       " Installer. You need to install PyMod components manually. Please"
                       " refer to the PyMod documentation for details about manual installation"
                       " of PyMod components." % (sys.platform))
            self.main_window.show_error_message("Installation Error", message)
            return None


        #---------------------------------------------------------------
        # Get some parameters before proceeding with the installation. -
        #---------------------------------------------------------------

        # Get the specs of the user's system.
        python_arch = get_python_architecture()
        os_arch = get_os_architecture()

        # Checks if the OS as a 64 bit architecture.
        if os_arch != "64":
            message = ("Your operating system has a %s bit architecture. The PyMod Installer"
                       " fully supports only systems with a 64 bit architecture. You may be"
                       " able to install MODELLER through Conda (if you have an incentive"
                       " PyMOL build) but you will not be able to download most of the"
                       " external tools of PyMod (BLAST+, HMMER, Clustal, MUSCLE, psipred, ksdssp)."
                       " If you need to use these tools, you will have to install and"
                       " configure them manually. Please refer to the PyMod documentation"
                       " for details on how to do it." % (os_arch))
            self.main_window.show_warning_message("Installation Warning", message)

        # Checks if the external tools directory exists and is empty. On PyMod first session, it
        # should exists and be empty.
        external_tools_dirpath = os.path.join(self.current_pymod_dirpath, self.external_tools_dirname)

        populated_tools_dir = False
        if os.path.isdir(external_tools_dirpath):
            populated_tools_dir = bool(os.listdir(external_tools_dirpath))

        # Check if Conda is integrated in PyMOL.
        try:
            import conda
            has_conda = True
        except ImportError:
            has_conda = False

        # Get the status of MODELLER.
        try:
            import modeller
            # MODELLER is installed correctly.
            modeller_status = "installed"
        except:
            modeller_spec = importlib.util.find_spec("modeller")
            # MODELLER is not installed.
            if modeller_spec is None:
                modeller_status = "missing"
            # License key is probably invalid.
            else: 
                try:
                    import modeller
                except ImportError as e:
                    error = e
                    print(f"Error importing some_library: {e}")
                except Exception as e:
                    error = e
                    print(f"An unexpected error occurred: {e}")
            
                if "DLL load failed" in str(error) or "dll" in str(error):
                    modeller_status = "broken-modeller"
                else:
                    modeller_status = "broken"

        # Show the installation options dialog.
        install_options_dialog = Installation_options_dialog(pymod=self,
                                                             populated_tools_dir=populated_tools_dir,
                                                             has_conda=has_conda,
                                                             modeller_status=modeller_status,
                                                             has_network=has_network,
                                                             os_arch=os_arch)
        install_options_dialog.setModal(True)
        install_options_dialog.exec_()

        install_modeller_flag = install_options_dialog.get_modeller_cb()
        external_tools_flag = install_options_dialog.get_external_tools_state()
        if not install_options_dialog.proceed_flag:
            return None


        #----------------------------------------------
        # Dialogs to actually install the components. -
        #----------------------------------------------

        # External tools installation.
        if external_tools_flag:

            # Show the external components download dialog.
            install_dialog = External_components_dialog(pymod=self,
                                                        url=packages_url_dict[sys.platform],
                                                        os_arch=os_arch,
                                                        local_mode=use_external_tools_local_install)
            install_dialog.setModal(True)
            install_dialog.exec_()

            # Finishes the installation.
            if install_dialog.complete_status:
                self.get_parameters_from_configuration_file()


        # MODELLER installation.
        if install_modeller_flag:

            # Show the MODELLER installation download dialog.
            if modeller_status == "missing":
                install_dialog = Install_MODELLER_dialog(pymod=self)
                install_dialog.setModal(True)
                install_dialog.exec_()

            elif modeller_status == "broken-modeller":
                install_dialog = Install_MODELLER_dialog(pymod=self)
                install_dialog.setModal(True)
                install_dialog.exec_()

            # Asks the MODELLER license key.
            elif modeller_status == "broken":
                modeller_key = get_modeller_license_key_dialog(pymod=self)
                if modeller_key is None:
                    return None

                try:
                    modeller_config_filepath = set_modeller_license_key(modeller_key)
                    message = ("The license key was sucessfully inserted in the MODELLER"
                               " configuration file (at: '%s'). Please restart PyMOL to"
                               " check if the key is valid and if MODELLER was correctly"
                               " installed." % modeller_config_filepath)
                    self.main_window.show_info_message("License Key Set", message)

                except Exception as e:
                    message = ("Could not insert the MODELLER license key because of the"
                               " following error: %s" % e)
                    self.main_window.show_error_message("MODELLER Error", message)


###################################################################################################
# Qt Dialogs.                                                                                     #
###################################################################################################

class Installer_dialog_mixin:
    """
    Mixin class to be incorporated in all the installation process dialogs.
    """

    def keyPressEvent(self, event):
        """
        By overriding this method, the dialog will not close when pressing the "esc" key.
        """
        if event.key() == QtCore.Qt.Key_Escape:
            pass
        else:
            QtWidgets.QDialog.keyPressEvent(self, event)


################################
# Installation options dialog. #
################################

class Installation_options_dialog(Installer_dialog_mixin, QtWidgets.QDialog):
    """
    Dialog to select the components to install. This is shown when launching the installer from
    the PyMod main window.
    """

    is_pymod_window = True

    def __init__(self, pymod,
                 populated_tools_dir,
                 has_conda, modeller_status,
                 has_network,
                 os_arch):

        self.pymod = pymod
        QtWidgets.QDialog.__init__(self, parent=self.pymod.main_window)

        self.populated_tools_dir = populated_tools_dir
        self.external_tools_dirpath = os.path.join(pymod.current_pymod_dirpath,
                                                   pymod.external_tools_dirname)
        self.has_conda = has_conda
        self.modeller_status = modeller_status
        self.has_network = has_network
        self.os_arch = os_arch

        self.proceed_flag = False

        self.initUI()
        # self.resize(QtCore.QSize(500, 180)) # (450, 140), (470, 180)


    def initUI(self):

        self.setWindowTitle('PyMod Tools Installation Options')

        vertical_layout = QtWidgets.QVBoxLayout()

        # Installation options label.
        self.download_info_label = QtWidgets.QLabel("Select components to install:", self)
        # self.download_info_label.setStyleSheet(label_style_1)
        vertical_layout.addWidget(self.download_info_label)

        # Install external tools checkbutton.
        if self.os_arch == "64":
            if self.has_network:
                external_tools_label = "External Tools (BLAST+, HMMER, Clustal, MUSCLE, psipred, ksdssp)"
                external_tools_cb_status = True
                external_tools_cb_checked = not self.populated_tools_dir
            else:
                external_tools_label = "External Tools (can not install, no internet connection)"
                external_tools_cb_status = False
                external_tools_cb_checked = False
        else:
            external_tools_label = "External Tools (can not install, only 64 bit OS are supported)"
            external_tools_cb_status = False
            external_tools_cb_checked = False

        self.external_tools_cb = QtWidgets.QCheckBox(external_tools_label, self)
        self.external_tools_cb.setChecked(external_tools_cb_checked)
        self.external_tools_cb.setEnabled(external_tools_cb_status)
        # self.external_tools_cb.setStyleSheet(label_font_1)
        vertical_layout.addWidget(self.external_tools_cb)

        # Install MODELLER checkbutton.
        install_modeller_label = None
        install_modeller_cb_status = None


        # Decide which actions to perform according to MODELLER status.
        fix_modeller_label = "Fix MODELLER (set licence key)"
        unknown_modeller_status_label = "Unknown MODELLER status"

        if self.modeller_status == "installed":

            install_modeller_label = "MODELLER (already installed)"
            install_modeller_cb_status = False

        else:

            if self.has_network:

                if self.has_conda:

                    if self.modeller_status == "missing":
                        install_modeller_label = "MODELLER (through Conda)"
                        install_modeller_cb_status = True
                    elif self.modeller_status == "broken":
                        install_modeller_label = fix_modeller_label
                        install_modeller_cb_status = True
                    elif self.modeller_status == "broken-modeller":
                        install_modeller_label = "MODELLER (through Conda)"
                        install_modeller_cb_status = True                        
                    else:
                        install_modeller_label = unknown_modeller_status_label
                        install_modeller_cb_status = False
                else:
                    install_modeller_label = "MODELLER (can not install automatically, PyMOL has no Conda)"
                    install_modeller_cb_status = False

            else:
                if self.modeller_status == "missing":
                    install_modeller_label = "MODELLER (can not install automatically, no internet connection)"
                    install_modeller_cb_status = False
                elif self.modeller_status == "broken":
                    install_modeller_label = fix_modeller_label
                    install_modeller_cb_status = True
                else:
                    install_modeller_label = unknown_modeller_status_label
                    install_modeller_cb_status = False


        self.install_modeller_cb = QtWidgets.QCheckBox(install_modeller_label, self)
        self.install_modeller_cb.setChecked(install_modeller_cb_status)
        self.install_modeller_cb.setEnabled(install_modeller_cb_status)
        # self.install_modeller_cb.setStyleSheet(label_font_1)
        # self.install_modeller_cb.stateChanged.connect(lambda:self.click_button(self.install_modeller_cb))

        vertical_layout.addWidget(self.install_modeller_cb)
        vertical_layout.addStretch(1)
        # m = 10
        # vertical_layout.setContentsMargins(m, m, m, m)

        # Button for starting the installation.
        horizontal_layout = QtWidgets.QHBoxLayout()
        self.start_button = QtWidgets.QPushButton('Start', self)
        # self.start_button.setStyleSheet(label_style_2)
        self.start_button.clicked.connect(self.on_button_click)
        horizontal_layout.addWidget(self.start_button)
        if not external_tools_cb_status and not install_modeller_cb_status:
            self.start_button.setEnabled(False)

        # Button for obtaining installation information.
        self.info_button = QtWidgets.QPushButton('Info', self)
        # self.info_button.setStyleSheet(label_style_2)
        self.info_button.clicked.connect(self.on_info_button_click)
        horizontal_layout.addWidget(self.info_button)

        horizontal_layout.addStretch(1)

        # Button for canceling the installation.
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        # self.cancel_button.setStyleSheet(label_style_2)
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)
        self.setLayout(vertical_layout)


    # Interactions with the buttons.
    def on_button_click(self):

        if not (self.external_tools_cb.isChecked() or self.install_modeller_cb.isChecked()):
            title = "Installation Warning"
            message = "Please select on the two components in order to proceed to their installation."
            self.pymod.main_window.show_warning_message(title, message)
            return None

        if self.populated_tools_dir and self.external_tools_cb.isChecked():
            title = 'Installation Warning'
            message = ("Looks like that your PyMod External Tools directory (%s) is not empty."
                       " This might mean that you have already installed some PyMod External"
                       " Tools in that directory (either manually, or by using this same automatic"
                       " installer before). Performing the automatic External Tools installation"
                       " will overwrite that directory and all files within it will be lost. Would you"
                       " like to proceed with the automatic installation?" % self.external_tools_dirpath)
            reply = askyesno_qt(title, message, self.pymod.get_qt_parent())
            if not reply:
                return None

        self.proceed_flag = True
        self.close()

    def on_cancel_button_click(self):
        self.proceed_flag = False
        self.close()

    def on_info_button_click(self):
        message = ("By using this dialog you can automatically install PyMod Components.\n\n"
                   "External Tools installation: if you decide to install the 'External Tools' (the first"
                   " checkbutton in the dialog), PyMod will download a ZIP file (%s) containing all the"
                   " executable and data files of the components and it will configure them so that they can be"
                   " readily used in the plugin.\n\n"
                   "MODELLER installation: if you decide to install MODELLER (the second checkbutton),"
                   " PyMod will use the Conda package manager to download and install MODELLER. This"
                   " function is accessible only if your PyMOL build is integrated with a Conda package"
                   " manager (see PyMod manual for more details). Please note that you need a MODELLER"
                   " licence key to perform this operation (which can be readily obtained from"
                   " https://salilab.org/modeller/registration.html)." % (packages_url_dict[sys.platform]))
        self.pymod.main_window.show_info_message("PyMod Components Installation", message)

    # Get parameters from the GUI.
    def get_external_tools_state(self):
        return self.external_tools_cb.isChecked()

    def get_modeller_cb(self):
        return self.install_modeller_cb.isChecked()


##############################################
# Download progress and installation dialog. #
##############################################

class External_tools_download_thread(QtCore.QThread):
    """
    Runs a download thread to download PyMod components.
    """

    # Signals.
    update_progressbar = QtCore.pyqtSignal(int)
    get_current_size = QtCore.pyqtSignal(int, str)
    get_total_size = QtCore.pyqtSignal(int)
    emit_exception = QtCore.pyqtSignal(Exception)
    signal_start_install = QtCore.pyqtSignal(int)

    def set_params(self, url, download_path):
        self.url = url
        self.download_path = download_path
        self.block_size = 8192
        self.write_mode = "wb"

    @catch_errors_installer_threads
    def run(self):

        # Connects with the web server.
        with urllib.request.urlopen(self.url) as http_resp:

            file_size = http_resp.headers["Content-Length"]

            if file_size is None:
                file_size = 0
            else:
                file_size = int(file_size)
            self.get_total_size.emit(file_size)

            with open(self.download_path, self.write_mode) as o_fh:

                file_size_dl = 0
                chunks_count = 0

                while True:

                    buffer_data = http_resp.read(self.block_size)
                    file_size_dl += len(buffer_data)

                    # Updates the fraction in the progressbar.
                    if file_size != 0:
                        # Can compute a fraction of the total size of the file to download.
                        frac = file_size_dl/file_size*100
                        self.update_progressbar.emit(int(frac))

                    # Updates the MB values in the progressbar.
                    if chunks_count % 100 == 0:
                        self.get_current_size.emit(file_size_dl, "update")
                    else:
                        self.get_current_size.emit(0, "pass")

                    if not buffer_data:
                        break

                    o_fh.write(buffer_data)
                    chunks_count += 1

                # Completes and updates the GUI.
                self.get_current_size.emit(file_size_dl, "completed")
                self.update_progressbar.emit(100)
                # Wait a little bit of time.
                time.sleep(0.5)
                self.signal_start_install.emit(0)


class External_tools_installation_thread(QtCore.QThread):
    """
    Runs a installation thread to install PyMod components.
    """

    # Signals.
    complete_installation = QtCore.pyqtSignal(int)
    emit_exception = QtCore.pyqtSignal(Exception)

    # Name of the PyMod tools present in the installer package.
    tool_names = ["clustalw", "muscle", "clustalo", "blast_plus", "psipred", "ksdssp", "hmmer"]

    # Name of the data directories of the installer.
    data_components_names = [blast_databases_dirname, hmmer_databases_dirname]

    def set_params(self, archive_filepath,
                   backup_archive_filepath,
                   pymod_dirpath,
                   pymod_cfg_file_path,
                   external_tools_dirname,
                   data_dirname,
                   os_arch,
                   use_blast_v5_databases=False,
                   from_main_window=True):

        self.archive_filepath = archive_filepath
        self.backup_archive_filepath = backup_archive_filepath
        self.pymod_dirpath = pymod_dirpath
        self.pymod_cfg_file_path = pymod_cfg_file_path
        self.external_tools_dirname = external_tools_dirname
        self.data_dirname = data_dirname
        self.os_arch = os_arch
        self.use_blast_v5_databases = use_blast_v5_databases
        if self.use_blast_v5_databases:
            self.default_psipred_db = "swissprot"
        else:
            self.default_psipred_db = "swissprot_v4"
        self.from_main_window = from_main_window

        self.pymod_external_tools_dirpath = os.path.join(self.pymod_dirpath, self.external_tools_dirname)
        self.pymod_data_dirpath = os.path.join(self.pymod_dirpath, self.data_dirname)

    @catch_errors_installer_threads
    def run(self):
        """
        Unzips the archive containing the files of PyMod external tools.
        """

        #-----------------------------------------------------------
        # Builds a directory where to extract the temporary files. -
        #-----------------------------------------------------------

        pymod_files_dirpath = os.path.join(os.path.dirname(self.archive_filepath), "installation_files")
        source_pymod_files_dirpath = os.path.join(pymod_files_dirpath, "pymod_files")
        if os.path.isdir(pymod_files_dirpath):
            shutil.rmtree(pymod_files_dirpath)
        os.mkdir(pymod_files_dirpath)

        if self.backup_archive_filepath is not None:
            shutil.copy(self.backup_archive_filepath, self.archive_filepath)

        #-----------------------------------------
        # Check if the file is a valid zip file. -
        #-----------------------------------------

        if not zipfile.is_zipfile(self.archive_filepath):
            raise TypeError("The '%s' file is not a zip file." % self.archive_filepath)


        #---------------------------------------------
        # Extract the file to a temporary directory. -
        #---------------------------------------------

        shutil.unpack_archive(self.archive_filepath, pymod_files_dirpath, format="zip")


        #-------------------------------------------------------------
        # Put PyMod tools executable files in the right directories. -
        #-------------------------------------------------------------

        source_external_tools_dirpath = os.path.join(source_pymod_files_dirpath,
                                                     self.external_tools_dirname,
                                                     str(self.os_arch))

        # Build an 'external_tools' directory.
        if os.path.isdir(self.pymod_external_tools_dirpath):
            shutil.rmtree(self.pymod_external_tools_dirpath)
        os.mkdir(self.pymod_external_tools_dirpath)

        # Copy all tools files.
        installed_tools_list = []
        for tool_name in self.tool_names:
            tool_source_dirpath = os.path.join(source_external_tools_dirpath, tool_name)
            if not os.path.isdir(tool_source_dirpath):
                continue
            installed_tools_list.append(tool_name)
            tool_dest_dirpath = os.path.join(self.pymod_external_tools_dirpath, tool_name)

            # Copy the files in the target PyMod directory.
            shutil.move(tool_source_dirpath, tool_dest_dirpath)

            # Change executable file permissions on UNIX systems.
            if os.name == "posix" and tool_name in ("clustalw", "clustalo", "muscle", "ksdssp", "psipred", "blast_plus", "hmmer"):
                for exec_file in os.listdir(os.path.join(tool_dest_dirpath, "bin")):
                    command = "chmod 777 '%s'" % (os.path.join(tool_dest_dirpath, "bin", exec_file))
                    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()

        # Leave a log file in the data directory to store installation information.
        with open(os.path.join(self.pymod_external_tools_dirpath, tools_installer_log_filename), "w") as l_fh:
            l_fh.write(get_formatted_date() + "\n")


        #-------------------------------------------------
        # Put PyMod data files in the right directories. -
        #-------------------------------------------------

        # Istall PyMod tools (except MODELLER).
        source_data_dirpath = os.path.join(source_pymod_files_dirpath, self.data_dirname)

        # Build a 'data' directory and its subdirectories.
        if not os.path.isdir(self.pymod_data_dirpath):
            os.mkdir(self.pymod_data_dirpath)

        for comp_name in self.data_components_names:
            comp_dirpath = os.path.join(self.pymod_data_dirpath, comp_name)
            if not os.path.isdir(comp_dirpath):
                os.mkdir(comp_dirpath)

        # Moves the files of all components.
        last_downloaded_dict = {}
        for comp_name in self.data_components_names:

            component_source_dirpath = os.path.join(source_data_dirpath, comp_name)
            component_dest_dirpath = os.path.join(self.pymod_data_dirpath, comp_name)

            # Copy the files in the target PyMod directory.
            for element_name in os.listdir(component_source_dirpath):

                if not ("swissprot" in element_name or "pdbaa" in element_name):
                    continue
                element_source_dirpath = os.path.join(component_source_dirpath, element_name)
                element_dest_dirpath = os.path.join(component_dest_dirpath, element_name)

                # Moves directories.
                if os.path.isdir(element_source_dirpath):
                    if os.path.isdir(element_dest_dirpath):
                        shutil.rmtree(element_dest_dirpath)
                    shutil.move(element_source_dirpath, element_dest_dirpath)

                # Moves files.
                if os.path.isfile(element_source_dirpath):
                    if os.path.isfile(element_dest_dirpath):
                        os.remove(element_dest_dirpath)
                    shutil.move(element_source_dirpath, element_dest_dirpath)


        # Leave a log file in the data directory to store the date of installation.
        last_downloaded_dict = {}

        # Gets the download dates from the installer files.
        source_db_downloa_filepath = os.path.join(source_data_dirpath, data_installer_log_filename)
        if os.path.isfile(source_db_downloa_filepath):
            with open(source_db_downloa_filepath, "r") as l_fh:
                last_downloaded_dict = json.loads(l_fh.read())

        # Gets the download dates of additional databases components from the PyMod data directory
        # file.
        db_download_filepath = os.path.join(self.pymod_data_dirpath, data_installer_log_filename)
        if os.path.isfile(db_download_filepath):
            with open(db_download_filepath, "r") as l_fh:
                 old_last_downloaded_dict = json.loads(l_fh.read())
            for k in old_last_downloaded_dict:
                if not k in last_downloaded_dict:
                    last_downloaded_dict[k] = old_last_downloaded_dict[k]

        # Writes the log file.
        with open(db_download_filepath, "w") as l_fh:
            l_fh.write(json.dumps(last_downloaded_dict))


        #-------------------------------------
        # Updating PyMod configuration file. -
        #-------------------------------------

        # Loads the configuration file.
        with open(self.pymod_cfg_file_path, "r") as c_fh:
            config_dict = json.loads(c_fh.read())

        # Tools with only one executable.
        for tool_name in ("clustalw", "clustalo", "muscle", "ksdssp"):
            if not tool_name in installed_tools_list:
                continue
            config_dict[tool_name] = {}
            tool_exe_filepath = os.path.join(self.pymod_external_tools_dirpath, tool_name, "bin", get_exe_file_name(tool_name))
            config_dict[tool_name]["exe_file_path"] = tool_exe_filepath

        # Tools with multiple executables and databases.
        for tool_name in ("blast_plus", "hmmer", "psipred"):

            if not tool_name in installed_tools_list:
                continue

            config_dict[tool_name] = {}
            tool_exe_dirpath = os.path.join(self.pymod_external_tools_dirpath, tool_name, "bin")
            config_dict[tool_name]["exe_dir_path"] = tool_exe_dirpath

            if tool_name == "blast_plus":
                config_dict[tool_name]["database_dir_path"] = os.path.join(self.pymod_data_dirpath, blast_databases_dirname)

            elif tool_name == "hmmer":
                config_dict[tool_name]["database_dir_path"] = os.path.join(self.pymod_data_dirpath, hmmer_databases_dirname)
                config_dict[tool_name]["hmmscan_db_dir_path"] = os.path.join(self.pymod_data_dirpath, hmmscan_databases_dirname)

            elif tool_name == "psipred":
                config_dict[tool_name]["database_dir_path"] = os.path.join(self.pymod_data_dirpath, blast_databases_dirname, self.default_psipred_db)

        # Psipred data directory.
        if "psipred" in installed_tools_list:
            config_dict["psipred"]["data_dir_path"] = os.path.join(self.pymod_external_tools_dirpath, "psipred", "data")

        # MODELLER.
        # config_dict["modeller"] = {"use_importable_modeller": True}

        # Writes the configuration file.
        with open(self.pymod_cfg_file_path, "w") as c_fh:
            c_fh.write(json.dumps(config_dict))


        #----------------------------
        # Cleaning temporary files. -
        #----------------------------

        shutil.rmtree(pymod_files_dirpath)
        os.remove(self.archive_filepath)


        #------------------------------------------------------------------
        # Updates the GUI so that the user can complete the installation. -
        #------------------------------------------------------------------

        if self.from_main_window:
            self.complete_installation.emit(0)


class External_components_dialog(Installer_dialog_mixin, QtWidgets.QDialog):
    """
    A dialog to download from GitHub an archive with binary files for the external tools of PyMod.
    """

    is_pymod_window = True

    def __init__(self, pymod, url, os_arch, local_mode=False):

        self.pymod = pymod
        QtWidgets.QDialog.__init__(self, parent=self.pymod.main_window)

        download_path = os.path.join(pymod.temp_directory_dirpath, "pymod_files.zip")
        self.local_mode = local_mode

        # Performs "local" installation.
        self.backup_download_path = None
        if self.local_mode:
            backup_download_path = askopenfile_qt("Select archive file", name_filter="*.zip", parent=self.get_qt_parent())
            if backup_download_path == "" or not os.path.isfile(backup_download_path):
                self.pymod.main_window.show_warning_message("Input Warning", "Invalid archive file.")
            else:
                self.backup_download_path = backup_download_path

        self.mb_label = "0"
        self.tot_label = "0"
        self.total_size = None
        self.complete_status = False
        self.initUI()


        # Initializes the download thread.
        self.dl_thread = External_tools_download_thread()
        self.dl_thread.update_progressbar.connect(self.on_update_progressbar)
        self.dl_thread.set_params(url, download_path)
        self.dl_thread.get_current_size.connect(self.on_get_current_size)
        self.dl_thread.get_total_size.connect(self.on_get_total_size)
        self.dl_thread.emit_exception.connect(self.on_emit_exception)
        self.dl_thread.signal_start_install.connect(self.start_install_thread)

        # Initializes the installation thread.
        self.inst_thread = External_tools_installation_thread()
        # self.inst_thread.update_progressbar.connect(self.on_update_progressbar)
        self.inst_thread.set_params(archive_filepath=download_path,
                                    pymod_dirpath=pymod.current_pymod_dirpath,
                                    pymod_cfg_file_path=pymod.cfg_file_path,
                                    external_tools_dirname=pymod.external_tools_dirname,
                                    data_dirname=pymod.data_dirname,
                                    os_arch=os_arch,
                                    backup_archive_filepath=self.backup_download_path,
                                    use_blast_v5_databases=self.pymod.use_blast_v5_databases)
        self.inst_thread.complete_installation.connect(self.on_complete_installation)
        self.inst_thread.emit_exception.connect(self.on_emit_exception)


    def initUI(self):

        self.setWindowTitle('Download and Install External Tools')

        vertical_layout = QtWidgets.QVBoxLayout()
        self.download_progressbar = QtWidgets.QProgressBar(self)
        # self.progress.setGeometry(0, 0, 340, 25)
        self.download_progressbar.setMaximum(100)
        if not self.local_mode:
            progressbar_text = "Click the button to start components download and installation"
        else:
            if self.backup_download_path is not None:
                progressbar_text = "Click the button to start components installation"
            else:
                progressbar_text = "Invalid archive file"

        self.download_progressbar.setFormat("")
        self.download_progressbar.setValue(0)
        vertical_layout.addWidget(self.download_progressbar)
        self.progressbar_label = QtWidgets.QLabel(progressbar_text)
        vertical_layout.addWidget(self.progressbar_label)

        # Button for starting the installation.
        horizontal_layout = QtWidgets.QHBoxLayout()
        if not self.local_mode:
            start_button_text = 'Start Download'
        else:
            if self.backup_download_path is not None:
                start_button_text = 'Start Installation'
            else:
                start_button_text = 'Exit Installation'
        self.start_button = QtWidgets.QPushButton(start_button_text, self)
        # self.start_button.setStyleSheet(label_style_2)
        self.start_button.clicked.connect(self.on_button_click)
        horizontal_layout.addWidget(self.start_button)

        horizontal_layout.addStretch(1)

        # Button for canceling the installation.
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        # self.cancel_button.setStyleSheet(label_style_2)
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)

        self.setLayout(vertical_layout)


    # Interactions with the buttons.
    def on_button_click(self):
        # The installation has not been launched yet.
        if not self.complete_status:
            # Remote installation.
            if not self.local_mode:
                self.dl_thread.start()
                self.download_progressbar.setFormat("Connecting...")
                self.progressbar_label.setText("")
                self.start_button.setEnabled(False)
            # Local installation.
            else:
                if self.backup_download_path is not None:
                    self.start_install_thread()
                    self.start_button.setEnabled(False)
                else:
                    self.close()
        # The installation has been already launched.
        else:
            self.close()

    def on_cancel_button_click(self):
        self._terminate_threads()
        self.close()

    def closeEvent(self, evnt):
        self._terminate_threads()

    def _terminate_threads(self):
        if self.dl_thread.isRunning():
            self.dl_thread.terminate()
        if self.inst_thread.isRunning():
            self.inst_thread.terminate()


    # Interactions with the download thread.
    def on_get_total_size(self, value):
        """
        Gets the total size of the archive to download.
        """
        self.total_size = value
        self.tot_label = self._convert_to_mb_string(value)
        if self.total_size == 0:
            pass

    def on_update_progressbar(self, value):
        """
        Updates the percentage in the progressbar.
        """
        self.download_progressbar.setValue(value)

    def on_get_current_size(self, value, state):
        """
        Updates the MB values in the progressbar.
        """
        if not state in ("update", "pass", "completed"):
            raise KeyError(state)

        if state in ("update", "completed"):
            self.mb_label = self._convert_to_mb_string(value)

        if self.total_size != 0:
            if state == "completed":
                self.download_progressbar.setFormat(
                    "Download Completed: %p% (" + self.mb_label + " of " + self.tot_label + " MB)")
            else:
                self.download_progressbar.setFormat(
                    "Download Progress: %p% (" + self.mb_label + " of " + self.tot_label + " MB)")
        else:
            if state == "completed":
                self.download_progressbar.setFormat("Download Completed: " + self.mb_label + " MB")
            else:
                self.download_progressbar.setFormat("Download Progress: " + self.mb_label + " MB of ???")

    def _convert_to_mb_string(self, value):
        if 0 <= value < 100000:
            round_val = 3
        elif 100000 <= value:
            round_val = 1
        else:
            ValueError(value)
        return str(round(value/1000000.0, round_val))


    # Interactions with the installation thread.
    def start_install_thread(self, signal=0):
        self.inst_thread.start()
        # self.download_progressbar.setRange(0, 0)
        # self.download_progressbar.reset()
        self.download_progressbar.setValue(0)
        self.download_progressbar.setFormat("Unzipping (please wait...)")
        self.start_button.setText("Finish Installation")

    def on_complete_installation(self):
        self.download_progressbar.setValue(100)
        self.download_progressbar.setFormat("Installation Completed Successfully")
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.complete_status = True

    def on_emit_exception(self, e):
        """
        Quit the threads and close the installation dialog.
        """
        self._terminate_threads()

        message = "There was an error: %s" % str(e)
        if hasattr(e, "url"):
            message += " (%s)." % e.url
        else:
            message += "."
        message += " Quitting the installation process."
        self.pymod.main_window.show_error_message("Installation Error", message)

        self.close()


############################
# Install MODELLER dialog. #
############################

class MODELLER_conda_install_thread(QtCore.QThread):
    """
    Runs an installation thread to download the MODELLER conda package.
    """

    # Signals.
    complete_installation = QtCore.pyqtSignal(dict)

    def set_params(self, pymod_envs_dirpath, pymod_env_dirpath, pymod_pkgs_dirpath):
        self.pymod_envs_dirpath = pymod_envs_dirpath
        self.pymod_env_dirpath = pymod_env_dirpath
        self.pymod_pkgs_dirpath = pymod_pkgs_dirpath

    def set_modeller_key(self, modeller_key):
        self.modeller_key = modeller_key

    def run(self):
        mod_install_result = perform_modeller_installation(self.modeller_key,
                                                           self.pymod_envs_dirpath,
                                                           self.pymod_env_dirpath,
                                                           self.pymod_pkgs_dirpath)
        self.complete_installation.emit(mod_install_result)


def perform_modeller_installation(modeller_key,
                                  pymod_envs_dirpath,
                                  pymod_env_dirpath,
                                  pymod_pkgs_dirpath,
                                  uninstall=False):

    try:

        print("\n# Installing the MODELLER conda package.")

        # Uses the conda python module.
        import conda.cli.python_api as conda_api

        # Gets environment information.
        info_results = conda_api.run_command(conda_api.Commands.INFO, "--json")
        conda_info_dict = json.loads(info_results[0])
        
        # NOTE: installing modeller on the default PyMOL environment sometimes causes
        # problems in PyMOL < 2.5.7, therefore it is safer to install it in a PyMod "pseudo" conda
        # environment.        
        pymol_version_for_modeller = tuple(map(int, cmd.get_version()[0].split('.')[:3]))
        minimum_pymol_version = (2, 5, 5)
        install_modeller_in_root_conda_env = False
        
        # Adds the salilab channel. This should write a .condarc file in the user's home.
        print("- Adding 'salilab' to the conda channels.")
        conda_api.run_command(conda_api.Commands.CONFIG, "--add", "channels", "salilab")

        # Uninstalls MODELLER if necessary. Not actually used right now in PyMod.
        if uninstall:
            stdout = conda_api.run_command(conda_api.Commands.REMOVE, "modeller")

        # Checks if we are in windows:
        if os.name == 'nt':
        
            # Checks if the conda root has write permissions.
            if "root_writable" in conda_info_dict:
                root_writable = conda_info_dict["root_writable"]
            else:
                root_writable = False
                
            # Checks if the user has Admin privileges.   
            if "is_windows_admin" in conda_info_dict:
                is_windows_admin = conda_info_dict["is_windows_admin"]
            else:
                is_windows_admin = False
                
            # Checks if PyMOL version is > 2.5.5. In this case, Modeller will be installed on the default PyMOL environment
            if pymol_version_for_modeller > minimum_pymol_version:
                print("Your PyMOL version is %s, modeller will be installed in root PyMOL conda env" % str(cmd.get_version()[0]))
                install_modeller_in_root_conda_env = True
            else:
                print("Your PyMOL version is %s, modeller will be installed in a new conda env" % str(cmd.get_version()[0]))
                install_modeller_in_root_conda_env = False
                
        # if we are in other OS, modeller will be installed in new conda env anyway:                
        else:
            install_modeller_in_root_conda_env = False 

        # Sets an environmental variable which will be used by conda to edit the MODELLER
        # configuration file.
        print("- Setting the MODELLER key as an environmental variable.")
        os.environ["KEY_MODELLER"] = modeller_key

        # Edit the condarc file and builds a new conda environment for PyMod.
        if not install_modeller_in_root_conda_env:

            # Add an environment directory.
            print("- Adding a writable envs directory.")
            if os.path.isdir(pymod_envs_dirpath):
                shutil.rmtree(pymod_envs_dirpath)
            os.mkdir(pymod_envs_dirpath)
            conda_api.run_command(conda_api.Commands.CONFIG, "--add", "envs_dirs", pymod_envs_dirpath)

            # Add a packages directory.
            print("- Adding a writable pkgs directory.")
            if os.path.isdir(pymod_pkgs_dirpath):
                shutil.rmtree(pymod_pkgs_dirpath)
            os.mkdir(pymod_pkgs_dirpath)
            conda_api.run_command(conda_api.Commands.CONFIG, "--add", "pkgs_dirs", pymod_pkgs_dirpath)

            # Creates a new environment. Use the same Python version of the Python
            # used in PyMOL.
            print("- Creating a new python %s conda environment for PyMod modules -" % python_minor_version)
            conda_api.run_command(conda_api.Commands.CREATE, "-p", pymod_env_dirpath, "--no-default-packages", "python=%s" % python_minor_version)

        # Installs the MODELLER package.
        
        # Checks if windows:
        if os.name == 'nt':
            if install_modeller_in_root_conda_env:
                if not root_writable and not is_windows_admin:
                    # message 
                    print("You must run PyMOL as Administrator in order to install Modeller!")
                    return {"successful": False, "error": "You must run PyMOL as Administrator in order to install Modeller!\n \
                    Please right-click the PyMOL icon and select 'Run as Administrator'"}
                else:
                    print("- Installing the 'modeller' conda package in root conda env")
                    stdout = conda_api.run_command(conda_api.Commands.INSTALL, "modeller")
            else:
                print("- Installing the 'modeller' conda package in a new conda env %s" % str(pymod_env_dirpath))
                stdout = conda_api.run_command(conda_api.Commands.INSTALL, "-p", pymod_env_dirpath, "modeller")
                
        # other OS:        
        else:
            if not install_modeller_in_root_conda_env:
                print("- Installing the 'modeller' conda package in a new conda env %s" % str(pymod_env_dirpath))
                stdout = conda_api.run_command(conda_api.Commands.INSTALL, "-p", pymod_env_dirpath, "modeller")
            else:    
                print("- Installing the 'modeller' conda package in root conda env")
                stdout = conda_api.run_command(conda_api.Commands.INSTALL, "modeller")
        
        print("- The 'modeller' package was successfully installed.")


        # Removes the temporary directories from the .condarc file.
        if not install_modeller_in_root_conda_env:
            conda_api.run_command(conda_api.Commands.CONFIG, "--remove", "envs_dirs", pymod_envs_dirpath)
            conda_api.run_command(conda_api.Commands.CONFIG, "--remove", "pkgs_dirs", pymod_pkgs_dirpath)

        return {"successful": True, "error": None}


    except Exception as e:

        return {"successful": False, "error": e}


class Install_MODELLER_dialog(Installer_dialog_mixin, QtWidgets.QDialog):
    """
    Dialog to install the MODELLER conda package.
    """

    is_pymod_window = True

    def __init__(self, pymod):

        self.pymod = pymod
        QtWidgets.QDialog.__init__(self, parent=self.pymod.main_window)

        self.complete_status = False

        self.initUI()
        self.adjustSize()

        # Initializes the package installation thread.
        self.mod_conda_thread = MODELLER_conda_install_thread()
        self.mod_conda_thread.complete_installation.connect(self.on_complete_installation)
        self.mod_conda_thread.set_params(pymod_envs_dirpath=pymod.pymod_envs_dirpath,
                                         pymod_env_dirpath=pymod.pymod_env_dirpath,
                                         pymod_pkgs_dirpath=pymod.pymod_pkgs_dirpath)


    def initUI(self):

        self.setWindowTitle('Install MODELLER package')

        # Installation progress bar.
        vertical_layout = QtWidgets.QVBoxLayout()
        self.mod_install_progressbar = QtWidgets.QProgressBar(self)
        self.mod_install_progressbar.setMaximum(100)
        vertical_layout.addWidget(self.mod_install_progressbar)

        # Installation progress label.
        self.mod_install_label = QtWidgets.QLabel("Click the button to start the MODELLER package installation", self)
        # self.mod_install_label.setStyleSheet(label_style_2)
        vertical_layout.addWidget(self.mod_install_label)

        # Button for starting the installation.
        horizontal_layout = QtWidgets.QHBoxLayout()
        self.start_button = QtWidgets.QPushButton('Start Installation', self)
        # self.start_button.setStyleSheet(label_style_2)
        self.start_button.clicked.connect(self.on_button_click)
        horizontal_layout.addWidget(self.start_button)

        horizontal_layout.addStretch(1)

        # Button for canceling the installation.
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        # self.cancel_button.setStyleSheet(label_style_2)
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)

        self.setLayout(vertical_layout)


    # Interactions with the buttons.
    def on_button_click(self):
        if not self.complete_status:

            # Asks the user to provide a MODELLER licence key.
            modeller_key = get_modeller_license_key_dialog(pymod=self.pymod)
            if modeller_key is None:
                return None
            self.modeller_key = modeller_key

            # Starts the installation thread.
            self.mod_conda_thread.set_modeller_key(self.modeller_key)
            self.mod_conda_thread.start()
            self.start_button.setEnabled(False)

            self.mod_install_progressbar.setRange(0, 0)
            self.mod_install_progressbar.setValue(0)
            self.mod_install_progressbar.setTextVisible(True)
            self.mod_install_progressbar.setAlignment(QtCore.Qt.AlignCenter)

            self.mod_install_label.setText("Retrieving and installing the MODELLER conda package...")

        else:
            self.close()

    def on_cancel_button_click(self):
        self._terminate_threads()
        self.close()

    def closeEvent(self, evnt):
        self._terminate_threads()

    def _terminate_threads(self):
        if self.mod_conda_thread.isRunning():
            self.mod_conda_thread.terminate()


    # Interactions with the installation thread.
    def on_complete_installation(self, mod_install_result):
        """
        Launched when the MODELLER installation process is terminated successfully.
        """

        if not mod_install_result["successful"]:
            self.mod_install_failure(mod_install_result)
            return None

        # Finds the MODELLER module to check if the licence key was correct. On the first installation
        # on Windows, this usually fails. It also fails when MODELLER is installed in the PyMod conda
        # environment.
        modeller_spec = importlib.util.find_spec("modeller")


        # Updates the GUI.
        if mod_install_result["successful"]:
            # Completes successfully.
            self.mod_install_progressbar.setRange(0, 100)
            self.mod_install_progressbar.setValue(100)
            self.mod_install_progressbar.setFormat("")

            installation_text = ("MODELLER Package installed. Please restart PyMOL"
                                 " to check if the installation was successful.")
            # if modeller_spec is not None:
            #     installation_text = ("MODELLER Package installed successfully."
            #                          " Please restart PyMOL to use MODELLER.")
            self.mod_install_label.setText(installation_text)

            self.start_button.setText("Finish Installation")
            self.start_button.setEnabled(True)
            self.cancel_button.setEnabled(False)
            self.complete_status = True

        else:
            self.mod_install_failure(mod_install_result)


    def mod_install_failure(self, mod_install_result):
        """
        Completes a failed MODELLER installation.
        """

        self.mod_install_progressbar.setRange(0, 100)
        self.mod_install_progressbar.setValue(0)
        self.mod_install_progressbar.setFormat("")
        self.mod_install_label.setText("MODELLER Conda Package installation failed")

        message = 'MODELLER installation failed because of the following error: \n\n%s.' % str(mod_install_result["error"])
        self.pymod.main_window.show_error_message("Installation Error", message)

        self.start_button.setText("Exit Installation")
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(True) # False
        self.complete_status = True


###########################
# Repair MODELLER dialog. #
###########################

def get_modeller_license_key_dialog(pymod):

    # Asks the user to provide a MODELLER licence key.
    modeller_key, pressed_ok = QtWidgets.QInputDialog.getText(pymod.main_window,
                                                              "Insert MODELLER Key",
                                                              "Insert a valid MODELLER licence key:",
                                                              QtWidgets.QLineEdit.Password, "")
    if not pressed_ok:
        return None

    modeller_key = re.sub("[^a-zA-Z0-9 _-]", "", modeller_key)

    if not modeller_key:
        message = "Please enter a valid MODELLER licence key."
        pymod.main_window.show_warning_message("Input Warning", message)
        return None

    return modeller_key


def set_modeller_license_key(modeller_key):
    """
    Edits the MODELLER configuration file to enter a licence key.
    """

    # Attempts to find the MODELLER config file by importing MODELLER. Usually does not work on
    # Windows.
    modeller_spec = importlib.util.find_spec("modeller")

    if modeller_spec is not None:
        modeller_dirpath = os.path.dirname(modeller_spec.origin)
        modeller_config_filepath = os.path.join(modeller_dirpath, "config.py")

    else:
        modeller_config_filepath = None

    if modeller_config_filepath is None or not os.path.isfile(modeller_config_filepath):
        raise ValueError("MODELLER 'config.py' not found.")

    with open(modeller_config_filepath, "r") as c_fh:
        old_mod_config_file_lines = c_fh.readlines()

    new_mod_config_file_lines = []
    licence_line_found = False
    for line in old_mod_config_file_lines:
        if line.startswith("license"):
            licence_line_found = True
            new_mod_config_file_lines.append("license = r'%s'" % modeller_key)
        else:
            new_mod_config_file_lines.append(line)

    if not licence_line_found:
        raise ValueError("'licence' line not found in %s" % modeller_config_filepath)

    with open(modeller_config_filepath, "w") as c_fh:
        c_fh.writelines(new_mod_config_file_lines)

    return modeller_config_filepath
