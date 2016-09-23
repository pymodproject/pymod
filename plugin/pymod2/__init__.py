# GUI.
import tkMessageBox

# Check for systemwide libraries.
import os
import sys
from pymod_lib import pymod_os_specific as pmos

# Cheks for NumPy.
global numpy_found
try:
    # Try to directly import NumPy in PyMOL.
    import numpy
    numpy_found = True
# If NumPy can not be imported, check for a systemwide NumPy installation.
except ImportError, e:
    system_numpy_path = pmos.find_systemwide_lib("numpy")
    if system_numpy_path:
        # If a systemwide NumPy was found, try to import it.
        try:
            sys.path.append(system_numpy_path)
            import numpy
            numpy_found = True
        except:
            numpy_found = False
    else:
        numpy_found = False

# Cheks for Biopython.
global biopython_found
try:
    pmos.check_biopython(raise_exception_on_fail=True)
    biopython_found = True
except ImportError, e:
    system_biopython_path = pmos.find_systemwide_lib("Bio")
    if system_biopython_path:
        try:
            sys.path.append(system_biopython_path)
            pmos.check_biopython(raise_exception_on_fail=True)
            biopython_found = True
        except:
            biopython_found = False
    else:
        biopython_found = False

# Check if there if there is a 'python_libs' directory in PyMod plugin folder and try to import from
# it the missing Pyhon libraries. The 'python_libs' folder may have been created by the PyMod
# Installer Bundle, if some Python libraries were missing.
pymod_plugin_dirpath = os.path.dirname(__file__)
python_libs_dirpath = os.path.join(pymod_plugin_dirpath, "python_libs")
if (not numpy_found or not biopython_found) and os.path.isdir(python_libs_dirpath):
    sys.path.insert(0, python_libs_dirpath)
    # Try to import the NumPy version provided by the PyMod Installer Bundle.
    if not numpy_found:
        try:
            import numpy
            numpy_found = True
        except:
            for m in sys.modules.keys():
                if m.startswith("numpy"):
                    del sys.modules[m]
            try:
                import numpy
                numpy_found = True
            except:
                pass
    # Try to import the Biopython version provided by the PyMod Installer Bundle.
    if not biopython_found:
        try:
            pmos.check_biopython(raise_exception_on_fail=True)
            biopython_found = True
        except:
            # Some old PyMOL builds have old Biopython versions lacking some new modules (such as
            # Phylo). If this old Biopython version was imported while initializing PyMOL, it must be
            # reloaded in order to import the new modules in PyMod.
            for m in sys.modules.keys():
                if m.startswith("Bio"):
                    del sys.modules[m]
            try:
                pmos.check_biopython(raise_exception_on_fail=True)
                biopython_found = True
            except:
                pass

# Sets the version of the PyMod plugin.
global __version__
__version__ = "2.0"
global __revision__
__revision__ = "3"
global pymod_plugin_name
pymod_plugin_name = "PyMod " + __version__

def __init__(self):
    """
    Adds a "PyMod" item to the plugin menu of PyMOL.
    """
    self.menuBar.addmenuitem('Plugin', 'command', pymod_plugin_name, label = pymod_plugin_name,
                             command = lambda s=self : startup_pymod(s))

def startup_pymod(app):
    """
    Checks if NumPy and Biopython (PyMod core Python dependencies) are present before launching
    PyMod.
    """
    if not numpy_found:
        title = "Import Error"
        message = "NumPy is not installed on your system. Please install it in order to use PyMod."
        tkMessageBox.showerror(title, message)
        return False
    if not biopython_found:
        title = "Import Error"
        message = "Biopython is not installed on your system. Please install it in order to use PyMod."
        tkMessageBox.showerror(title, message)
        return False
    import pymod_main
    pymod_main.pymod_launcher(app, pymod_plugin_name, __version__, __revision__)
