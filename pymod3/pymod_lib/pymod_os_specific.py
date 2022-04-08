# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for cross-platform compatibility.
"""

import os, sys
import platform
import re
import subprocess
import shutil
import zipfile
import struct
import urllib.request, urllib.error, urllib.parse
from datetime import datetime

from pymod_lib import pymod_vars


#####################################################################
# Directories managment.                                            #
#####################################################################

def get_home_dir():
    if sys.platform=="win32":
        d=os.path.join(os.environ["HOMEDRIVE"],os.environ["HOMEPATH"])
    else:
        d=os.getenv("HOME")
    if d and os.path.isdir(d):
        return d
    else:
        sys.stderr.write("Error! Cannot locate home folder!\n")
        return os.getcwd()

# def check_user_permissions(file_to_check):
#     if os.access(file_to_check, os.W_OK):
#         return True
#     else:
#         return False

# def is_writable(dirname):
#     """
#     Taken from the 'plugin.installation' module of PyMOL (since version 1.5.0.5).
#     """
#     path = os.path.join(dirname, '__check_writable')
#     try:
#         f = open(path, 'wb')
#         f.close()
#         os.remove(path)
#         return True
#     except (IOError, OSError):
#         return False


def get_pymod_cfg_paths():
    home_dirpath = get_home_dir()
    cfg_directory_path = os.path.join(home_dirpath, pymod_vars.pymod_cfg_dirname)
    pymod_envs_dirpath = os.path.join(cfg_directory_path, pymod_vars.pymod_envs_dirname)
    pymod_env_dirpath = os.path.join(pymod_envs_dirpath, pymod_vars.pymod_env_name)
    pymod_pkgs_dirpath = os.path.join(cfg_directory_path, pymod_vars.pymod_pkgs_dirname)
    return home_dirpath, cfg_directory_path, pymod_envs_dirpath, pymod_env_dirpath, pymod_pkgs_dirpath


#####################################################################
# Find external tools on users' systems.                            #
#####################################################################

def pymod_which(executable, paths=None):
    """ Find 'executable' in the directories listed in 'path'. """

    if sys.platform == "win32":
        ext=".exe"
    else:
        ext=""

    if not paths:
        paths=os.environ["PATH"]
    paths=paths.split(os.pathsep)

    if sys.platform == "darwin" and not "/usr/local/bin" in paths and os.path.isdir("/usr/local/bin"):
        paths.append("/usr/local/bin")

    exe_full_path=None

    # MODELLER.
    if executable == "modeller": # return directory on windows
        # try to locate modX.XX from importable modeller
        try:
            import modeller
            modinstall=os.sep.join(os.path.realpath(modeller.__file__).split(os.sep)[:-3])
        except:
            modinstall=None
        if modinstall:
            version=modinstall.split(os.sep)[-1].lstrip('modeller-')
            exe_full_path=find_executable("mod"+version,paths)

        if not exe_full_path and sys.platform!="win32":
            pattern = re.compile("^mod[0-9]{1,2}[.v][0-9]+$")
            for p in [p for p in set(paths) if p and os.path.isdir(p)]:
                f_list=[os.path.join(p,f) for f in os.listdir(p) if \
                    os.path.isfile(os.path.join(p,f)) and pattern.match(f)]
                if f_list:
                    exe_full_path=os.path.realpath(f_list[-1])
                    break

        # On Windows, use the Python interpreter to run modeller scripts.
        # TODO: locate win32 modeller when architecture is different
        elif not exe_full_path and sys.platform=="win32":
            # pattern=re.compile("^Python[0-9]{2}$")
            # for p in os.listdir(os.environ["HOMEDRIVE"]+'\\'):
            #     if pattern.match(p):
            #         f=os.sep.join([os.environ["HOMEDRIVE"],p,"python.exe"])
            #         if os.path.isfile(f):
            #             exe_full_path=f
            #             break
            pass

    # elif executable in ["TMalign","TMalign.exe","tmalign","tmalign.exe"]:
    #     exe_full_path=find_executable("TMalign"+ext,paths)
    #     if not exe_full_path:
    #         exe_full_path=find_executable("TMscore"+ext,paths)
    #
    # elif executable in ["dssp","dssp.exe"]:
    #     exe_full_path=find_executable(executable,paths)
    #     if not exe_full_path:
    #         exe_full_path=find_executable("mk"+executable,paths)

    elif executable=="clustalw":
        if sys.platform=="win32":
            if os.getenv("ProgramFiles(x86)"):
                p=os.getenv("ProgramFiles(x86)")
            else: # clustalw2 for win is 32bit only
                p=os.getenv("ProgramFiles")
            paths.append(os.path.join(p,"ClustalW2"))
        exe_full_path=find_executable("clustalw2"+ext,paths)
        if not exe_full_path:
            exe_full_path=find_executable("clustalw"+ext,paths)

    elif executable in ["psiblast","blast_plus"]:
        if sys.platform=="win32":
            for p in set([os.getenv("ProgramFiles(x86)"),os.getenv("ProgramFiles"),os.getenv("programw6432")]):
                if p:
                    ncbi_path=os.path.join(p,"NCBI")
                    if os.path.isdir(ncbi_path):
                        for d in os.listdir(ncbi_path):
                            if d.startswith("blast"):
                                blast_path=os.path.join(ncbi_path,d,"bin")
                                if os.path.isdir(blast_path):
                                    paths=[blast_path]+paths
        exe_full_path=find_executable("psiblast"+ext,paths)
        if not exe_full_path:
            exe_full_path=find_executable("blastpgp"+ext,paths)

    elif executable == "hmmer":
        exe_full_path=find_executable("phmmer"+ext,paths)

    else: # executable in [muscle,clustalo,predator], clustalw on posix
        exe_full_path=find_executable(executable+ext,paths)

    return exe_full_path # might be None


def find_executable(executable,paths=[]):
    exe_full_path=None
    for p in [p for p in set(paths) if p and os.path.isdir(p)]:
        f=os.path.join(p,executable)
        if os.path.isfile(f):
            exe_full_path=os.path.abspath(f)
            break
    # If the executble was not found, try to use the 'which' progam.
    if not exe_full_path and os.name == 'posix':
        try:
            exe_full_path = subprocess.check_output(["which",executable]).rstrip("\n")
        except:
            pass
    return exe_full_path


#####################################################################
# Python libraries.                                                 #
#####################################################################

def check_importable_modeller(get_exception=False):
    """
    Checks if systemwide MODELLER can be imported. If it can be imported, returns its version.
    """

    try:
        import modeller
        import _modeller
        import modeller.automodel
        from modeller.scripts import complete_pdb

        if hasattr(modeller, "__version__"):
            hasmodeller = modeller.__version__
        elif hasattr(_modeller, "mod_short_version_get"):
            hasmodeller = _modeller.mod_short_version_get()
        else:
            hasmodeller = "unknown"

        if not get_exception:
            return hasmodeller
        else:
            return None

    except Exception as e:
        if not get_exception:
            return ""
        else:
            return e


def check_valid_pyqt():

    try:
        import PyQt5
        pyqt5_version_str = PyQt5.QtCore.PYQT_VERSION_STR
        pyqt5_version = float(".".join(pyqt5_version_str.split(".")[0:2]))
        pyqt5_version_major = float(pyqt5_version_str.split(".")[0])
        pyqt5_version_minor = float(pyqt5_version_str.split(".")[1])
        valid_version = pyqt5_version_minor < 13
        return valid_version

    except Exception as e:
        return False


#####################################################################
# Archives.                                                         #
#####################################################################

def zip_directory(directory_path, zipfile_path):
    """
    Zips a directory. This function will ignore empty subdirectories inside the directory to zip.
    """
    if not os.path.isdir(directory_path):
        raise Exception("The target path (%s) is not a directory." % directory_path)
    directory_name = os.path.basename(directory_path)
    directory_parent = os.path.dirname(directory_path)
    directory_name_count = directory_path.split(os.path.sep).count(directory_name)
    zipfile_handle = zipfile.ZipFile(zipfile_path, 'w', zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(directory_path):
        for f in files:
            # Build the relative path of the zipped file.
            dir_index_in_root = 0
            for i,dirname in enumerate(root.split(os.path.sep)):
                if dirname == directory_name:
                    dir_index_in_root += 1
                    if dir_index_in_root == directory_name_count:
                        index_to_use = i
                        break
            # absolute_path = os.path.join(root, f)
            # relative_path = absolute_path[len(directory_parent)+len(os.path.sep):]
            # z.write(absolute_path, relative_path)
            relative_path = os.path.sep.join(root.split(os.path.sep)[index_to_use:])
            zipfile_handle.write(os.path.join(root, f), os.path.join(relative_path, f))
    zipfile_handle.close()


#####################################################################
# Files names and document viewing.                                 #
#####################################################################

def clean_file_name(file_name):
    """
    Replaces characters not allowed in Windows. The '/' character is invalid on UNIX and the ':' on
    Mac OS. Used in 'build_header_string()' in 'pymod_main.py'.
    """
    return re.sub("""[\\\\/:*?"<>|']""", "_", file_name)


def get_exe_file_name(program_name, force_not_extended=False):
    if sys.platform == "win32" and not force_not_extended:
        program_name += ".exe"
    return program_name


def open_document_with_default_viewer(document_path):
    # MAC.
    if sys.platform.startswith('darwin'):
        subprocess.call(('open', document_path))
    # WIN.
    elif os.name == 'nt':
        os.startfile(document_path)
    # Other UNIX.
    elif os.name == 'posix':
        subprocess.call(('xdg-open', document_path))


#####################################################################
# Commandline.                                                      #
#####################################################################

def build_commandline_path_string(file_path):
    if is_unix():
        return "'%s'" % (file_path)
    elif sys.platform == "win32":
        return '"%s"' % (file_path)
    else:
        return program_path


def build_commandline_file_argument(argument, path_name, extension=None):
    if extension:
        extension = "." + extension
    else:
        extension = ""
    return """ -%s %s""" % (argument, build_commandline_path_string(path_name+extension))


#####################################################################
# Operating systems and architectures.                              #
#####################################################################

os_names_dict = {"win32": "Windows", "darwin": "Mac OS X", "linux": "Linux"}

def get_python_architecture(): # "32" for x86
    """
    Gets the architecture of the PyMOL built in which PyMod is running. "64" for x86-64.
    """
    return str(8*struct.calcsize("P"))


def get_os_architecture():
    if sys.platform == "win32":
        if 'PROGRAMFILES(X86)' in os.environ:
            return "64"
        else:
            return "32"
    else:
        return platform.architecture()[0][0:2]


def is_unix():
    if sys.platform == "linux2" or sys.platform == "darwin" or os.name == "posix":
        return True
    else:
        return False


def get_unicode_support():
    ucs="0"
    if hasattr(sys,"maxunicode"):
        if sys.maxunicode==65535:
            ucs="2"
        elif sys.maxunicode==1114111:
            ucs="4"
    return ucs


def get_mac_ver():
    """
    Gets the version of Mac OS X.
    """
    mac_ver = platform.mac_ver()[0]
    # On some PyMOL builds from Schrodinger the platform.mac_ver()[0] returns an empty string.
    if not mac_ver:
        r = subprocess.Popen("sw_vers -productVersion", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()
        mac_ver = r[0].rstrip()
    return mac_ver


def get_linux_distribution():
    """
    Gets the Linux distribution name.
    """
    # subprocess.check_call("lsb_release -ic")
    if not float(get_python_version()) >= 2.7:
        return platform.linux_distribution()[0]
    else:
        return platform.dist()[0]

def get_python_version():
    return sys.version[:3]


#####################################################################
# Internet connection.                                              #
#####################################################################

def check_network_connection(remote_address, timeout=None):
    try:
        if timeout is None:
            response = urllib.request.urlopen(remote_address)
        else:
            response = urllib.request.urlopen(remote_address, timeout=timeout)
        return True
    except urllib.error.URLError:
        #as a last resort in countries where Google is banned...
        try:
            response = urllib.request.urlopen("https://www.ncbi.nlm.nih.gov", timeout=5)
            return True
        except:
            return False


#####################################################################
# Dates and times.                                                  #
#####################################################################

def get_formatted_date():
    return datetime.today().strftime('%d/%m/%Y (%H:%M)') #)'%d/%m/%Y (%H:%M:%S)')
