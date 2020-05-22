###########################################################################
# Copyright (C) 2020 Giacomo Janson, Alessandro Paiardini
# Copyright (C) 2016-2019 Giacomo Janson, Chengxin Zhang, Alessandro Paiardini
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

"""
PyMod 3: PyMOL Front-end to MODELLER and various other bioinformatics tools.
"""

# Module first imported by the PyMOL plugin system. Contains:
#     - code to check for the presence of Python libraries needed for PyMod to work.
#     - code to initialize PyMod as a PyMOL plugin.


#--------------------------
# Check Python libraries. -
#--------------------------

import os
import sys
from importlib import util
import shutil

from pymol import cmd


# Gets the Python version.
python_version = sys.version_info.major
python_minor_version = "%s.%s" % (sys.version_info.major, sys.version_info.minor)
python_micro_version = "%s.%s.%s" % (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)


# Checks for dependencies.
try:
    # Checks for some Qt bindings in PyMOL.
    from pymol.Qt import QtWidgets

    def showerror(title, message):
        QtWidgets.QMessageBox.critical(None, title, message)

    has_gui = "qt"
    pyqt_found = True

except ImportError:

    # Checks for Tkinter.
    try:
        if python_version == 3: # Python 3.
            from tkinter.messagebox import showerror
        else: # Python 2.
            from tkMessageBox import showerror
        has_gui = "tkinter"
    except ImportError: # On some open source builds, tkinter is missing.
        has_gui = None

    pyqt_found = False

try:
    import numpy
    numpy_found = True
except ImportError:
    numpy_found = False

try:
    import Bio
    biopython_found = True
except ImportError:
    biopython_found = False


# Sets the version of the PyMod plugin.
__pymod_version__ = "3.0"
__revision__ = "0"
__version__ = float(__pymod_version__ + __revision__.replace(".", ""))
pymod_plugin_name = "PyMod " + __pymod_version__


#----------------------------------------
# Initialize the PyMod plugin in PyMOL. -
#----------------------------------------

def __init_plugin__(app):
    """
    Initializes the plugin in the plugin menu of PyMOL 3.
    """

    # Adds a "PyMod" item to the plugin menu of PyMOL.
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt(pymod_plugin_name, startup_pymod)
    # Adds a "PyMod" item to the legacy plugin menu of PyMOL.
    except ImportError:
        app.menuBar.addmenuitem('Plugin', 'command', pymod_plugin_name,
                                label=pymod_plugin_name,
                                command=lambda a=app: startup_pymod(a))


def startup_pymod(app=None):
    """
    Executed when clicking on the 'PyMod' item in PyMOL's plugin menu.
    """

    if has_gui is None:
        print("\n# No GUI library (either Tkinter or Qt bindings) was found. PyMod"
              " can not be launched.")
        return None

    # Check if a PyMod main window is already open.
    if pyqt_found:
        try:
            for widget in QtWidgets.QApplication.instance().topLevelWidgets():
                if hasattr(widget, "is_pymod_main_window") and widget.isVisible():
                    title = "PyMod Error"
                    message = ("PyMod is already running. Please close its main"
                               " window or restart PyMOL in order to launch it again.")
                    showerror(title, message)
                    return None
        except Exception as e:
            pass

    # Checks if Python 3 is available.
    if python_version != 3:
        title = "Python Version Error"
        message = "PyMod %s requires Python 3. Your current Python version is %s." % (__pymod_version__, python_micro_version)
        showerror(title, message)
        return None

    # Checks the PyMOL version.
    pymol_version = float(".".join(cmd.get_version()[0].split(".")[0:2]))
    if pymol_version < 2.3:
        title = "PyMOL Version Error"
        message = "PyMod %s requires a PyMOL version of 2.3 or higher. Your current PyMOL version is %s." % (__pymod_version__, pymol_version)
        showerror(title, message)
        return None

    # Checks for PyQt.
    if not pyqt_found:
        title = "Import Error"
        message = "PyQt5 is not installed on your system. Please install it in order to use PyMod."
        showerror(title, message)
        return None

    # Checks if NumPy and Biopython (PyMod core Python dependencies) are present before launching
    # PyMod.
    if not numpy_found:
        title = "Import Error"
        message = "NumPy is not installed on your system. Please install it in order to use PyMod."
        showerror(title, message)
        return None
    if not biopython_found:
        title = "Import Error"
        message = "Biopython is not installed on your system. Please install it in order to use PyMod."
        showerror(title, message)
        return None


    # Adds to the sys.path the directory where the PyMod module is located.
    pymod_plugin_dirpath = os.path.dirname(__file__)
    if os.path.isdir(pymod_plugin_dirpath):
        sys.path.append(pymod_plugin_dirpath)


    # Attempts to import MODELLER from the default sys.path.
    modeller_spec = util.find_spec("modeller")

    parallel_modeller = False # A flag, when set as 'True' MODELLER parallel jobs can be used.

    # MODELLER can be imported.
    if modeller_spec is not None:
        parallel_modeller = True

    # MODELLER can not be imported. Tries to import it from a PyMod conda
    # pseudo-environment.
    else:

        # Adds to the sys.path the PyMod conda environment directory. Used to import
        # modules fetched by conda within PyMod.
        from pymod_lib import pymod_os_specific

        pymod_env_dirpath = pymod_os_specific.get_pymod_cfg_paths()[3]

        # Checks if a directory for a PyMod conda pseudo-environment exists.
        if os.path.isdir(pymod_env_dirpath):

            if sys.platform in ("linux", "darwin"):
                conda_lib_dirpath_list = [os.path.join(pymod_env_dirpath, "lib", "python%s" % python_minor_version),
                                          os.path.join(pymod_env_dirpath, "lib", "python%s" % python_minor_version, "lib-dynload"),
                                          os.path.join(pymod_env_dirpath, "lib", "python%s" % python_minor_version, "site-packages")]
                path_extend_list = []

            elif sys.platform == "win32":
                dll_dirpath = os.path.join(pymod_env_dirpath, "Lib", "site-packages")
                conda_lib_dirpath_list = [os.path.join(pymod_env_dirpath, "Library", "modeller", "modlib"),
                                          dll_dirpath]

                path_extend_list = [os.path.join(pymod_env_dirpath, "Scripts"),
                                    pymod_env_dirpath,
                                    dll_dirpath]

            else:
                conda_lib_dirpath_list = []
                path_extend_list = []

            # Add the environment directories to the sys.path of Python.
            for conda_lib_dirpath in conda_lib_dirpath_list:
                sys.path.append(conda_lib_dirpath)
            # Also add new directories to the PATH to use the MODELLER DLLs.
            if path_extend_list:
                os.environ["PATH"] += os.pathsep.join(path_extend_list)
            os.environ["HDF5_DISABLE_VERSION_CHECK"] = "1" # Prevents problems with hdf5 libraries versions.

            # Edit the 'modslave.py' file of MODELLER in the PyMod environment directory
            # (the user's original MODELLER files will remain untouched) so that the
            # MODELLER version installed in the PyMod environment directory can be
            # used to run parallel jobs. This operation will be performed only once,
            # that is, the first time PyMod is launched after MODELLER has been
            # installed through the GUI of the plugin.
            try:

                # Gets the installation directory of MODELLER. This can be acquired
                # also if the MODELLER key is not valid.
                from modeller import config

                modeller_bin_dirpath = os.path.join(config.install_dir, "bin")
                modslave_filepath = os.path.join(modeller_bin_dirpath, "modslave.py")
                modslave_bak_filepath = os.path.join(modeller_bin_dirpath, "modslave.py.bak")

                # Checks if the modslave file is present.
                if os.path.isfile(modslave_filepath):

                    # Checks if the backup modslave file is present (if it already
                    # exists, then PyMod has already modified the original modslave
                    # file).
                    if not os.path.isfile(modslave_bak_filepath):

                        print("\n# Configuring the MODELLER modslave.py in the PyMod environment.")
                        print("- modslave.py file found at: %s" % modslave_filepath)

                        # Build a backup copy of the original modslave file.
                        print("- Backing up the modslave.py file.")
                        shutil.copy(modslave_filepath, modslave_bak_filepath)

                        # Gets the content from the old modslave file.
                        with open(modslave_filepath, "r") as m_fh:
                            modfile_lines = m_fh.readlines()

                        insert_idx = 1
                        for line_idx, line in enumerate(modfile_lines):
                            if line.startswith("import sys"):
                                print("- Found insert line")
                                insert_idx = line_idx + 1
                                break

                        # This will add the paths above also to the 'modslave.py'
                        # file, so that MODELLER can be imported in child processes.
                        import_lines = ['\n# Added by the PyMod installer.\n',
                                        'import os\n',
                                        'import sys\n\n',
                                        'conda_lib_dirpath_list = %s\n' % repr(conda_lib_dirpath_list),
                                        'path_extend_list = %s\n\n' % repr(path_extend_list),
                                        'for path in conda_lib_dirpath_list:\n',
                                        '    sys.path.append(path)\n',
                                        'os.environ["PATH"] += os.pathsep.join(path_extend_list)\n',
                                        'os.environ["HDF5_DISABLE_VERSION_CHECK"] = "1"\n',
                                        '# End of the code inserted by PyMod.\n']

                        modfile_lines = modfile_lines[:insert_idx] + import_lines + modfile_lines[insert_idx:]

                        # Edits the modslave file.
                        with open(modslave_filepath, "w") as m_fh:
                            m_fh.writelines(modfile_lines)

                        print("- Finishing.")

                    parallel_modeller = True

            except Exception as e:
                print("- Could not import MODELLER from the PyMod environment: %s." % e)


    # Actually launches PyMod.
    from pymod_lib import pymod_main

    pymod_main.pymod_launcher(app=app,
                              pymod_plugin_name=pymod_plugin_name,
                              pymod_version=__pymod_version__,
                              pymod_revision=__revision__,
                              parallel_modeller=parallel_modeller)
