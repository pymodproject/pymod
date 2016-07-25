#!/usr/bin/env python
# os_specific.py: submodules for cross-platform compatibility
# Copyright (C) 2016 Chengxin Zhang, Giacomo Janson
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
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-
# 1301  USA
import os,sys
import platform
import re
import subprocess
import shutil
import struct
import pymol
from pymol import cmd
import platform
try:
    import pwd
except:
    pass

pymod_dir=os.path.dirname(os.path.abspath(__file__))
if (sys.platform in ["win32","darwin"] and \
    not pymod_dir in os.environ["PATH"].lower().split(os.pathsep)) or \
    not pymod_dir in os.environ["PATH"].split(os.pathsep):
    os.environ["PATH"]+=os.pathsep+pymod_dir # add current dir to PATH

# suffix of executables
exe_filetypes=[]
if sys.platform=="win32":
    exe_filetypes=[("Windows Executable",("*.exe","*.cmd","*.bat","*.com"))]

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


def get_working_dir():
    if sys.platform in ["darwin","win32"]:
        # Since MacOSX10.9, $HOME is no longer easily accessible from Finder
        d=os.path.join(get_home_dir(),"Desktop")
        if os.path.isdir(d):
            return d
        else:
            return get_home_dir()
    else: # linux2
        return os.getcwd() # on Win Vista +, pwd is often not writable


def check_user_permissions(file_to_check):
    if os.access(file_to_check, os.W_OK):
        return True
    else:
        return False


def is_writable(dirname):
    '''
    Taken from the 'plugin.installation' module of PyMOL (since version 1.5.0.5).
    '''
    path = os.path.join(dirname, '__check_writable')
    try:
        f = open(path, 'wb')
        f.close()
        os.remove(path)
        return True
    except (IOError, OSError):
        return False


#####################################################################
# Files managment.                                                  #
#####################################################################

def pymod_rm(file_path):
    if os.path.isfile(file_path):
        os.remove(file_path)
    elif os.path.isdir(file_path):
        shutil.rmtree(file_path)


def get_target_path_owner(path):
    stat_info = os.stat(path)
    return pwd.getpwuid(stat_info.st_uid).pw_name


#####################################################################
# Find external tools on users' systems.                            #
#####################################################################

def pymod_which(executable, paths=None):
    ''' Find 'executable' in the directories listed in 'path'. '''

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

# TODO: check systemwide importables. This might give problems on Mac OS.
def check_numpy():
    try: # check numpy
        import numpy
        return numpy.__version__, numpy.__file__
    except:
        return "",""


def check_biopython(raise_exception_on_fail=False):
    # Unpatched Bio.PDB requires md5, which was missing in PyMOL1.2/1.3
    # because its Python2.5 was not linked against OpenSSL libraries
    if not raise_exception_on_fail:
        try:
            import Bio.PDB, Bio, Bio.Phylo # Phylo was missing in PyMOL1.5
            from Bio.Align.Applications import ClustalwCommandline # This was missing in PyMOL 1.4.
            from Bio.Align.Applications import MuscleCommandline
            return Bio.__version__, Bio.__file__
        except:
            return "",""
    else:
        import Bio.PDB, Bio, Bio.Phylo
        from Bio.Align.Applications import ClustalwCommandline
        from Bio.Align.Applications import MuscleCommandline
        return Bio.__version__, Bio.__file__


def find_systemwide_lib(lib_name):
    """
    Check if a systemwide Python library can be imported in PyMOL.
    """
    lib_path = ""
    if sys.platform in ("linux2", "darwin"):
        try:
            # lib_file = subprocess.check_output(["/usr/bin/python2", "-Ec", "import %s; print %s.__file__" % (lib_name, lib_name)]).rstrip()
            r = subprocess.Popen("/usr/bin/python2 -Ec 'import %s; print %s.__file__'" % (lib_name, lib_name),
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
                                 ).communicate()
            lib_file = r[0].rstrip()
            if lib_file:
                lib_path = os.path.split(os.path.dirname(lib_file))[0]
        except:
            pass
    return lib_path


def check_importable_modeller():
    """
    Checks if systemwide MODELLER can be imported. If it can be imported, returns its version.
    """
    try:
        import modeller, _modeller
        import modeller.automodel
        from modeller.scripts import complete_pdb
        if hasattr(_modeller,"mod_short_version_get"):
            hasmodeller=_modeller.mod_short_version_get()
        else:
            hasmodeller=[e.lower()[8:].strip("-") for e in os.path.realpath(
                os.path.dirname(modeller.__file__)).split(os.sep) if e.lower(
                ).startswith("modeller") and len(e)>8][0]
        if not hasmodeller:
            hasmodeller="unknown"
    except:
        hasmodeller=""
    return hasmodeller


def find_modlib_path(modeller_exe_path=None):
    """
    Try to find an importbale modeller package on the user's system.
    """
    modeller_lib_path = ""
    # On Linux official PyMOL builds, the systemwide MODELLER can't be directly imported, so the
    # MODELLER Python library path needs to be appended to 'sys.path'.
    if sys.platform == "linux2":
        try:
            modeller_lib_path = find_systemwide_lib("modeller")
        except:
            pass

    elif sys.platform == "win32":
        # Checks 'modlib' on standard MODELLER installation paths.
        modeller_lib_path = find_modlib_dir_on_windows()
        # If it was not found, checks 'modlib' on the user supplied MODELLER directory.
        if not modeller_lib_path and modeller_exe_path:
            if os.path.isfile(modeller_exe_path):
                custom_modeller_dir_path = os.path.dirname(modeller_exe_path)
                for item in os.listdir(custom_modeller_dir_path):
                    if item == "modlib":
                        modeller_lib_path = os.path.join(custom_modeller_dir_path, item)

    elif sys.platform == "darwin":
        pass

    return modeller_lib_path


def find_modlib_dir_on_windows():
    """
    Check if there is a standard MODELLER installation directory.
    """
    program_files_paths = []
    if get_python_architecture() == "32":
        program_files_paths.append(os.getenv("ProgramFiles(x86)"))
    elif get_python_architecture() == "64":
        program_files_paths.append(os.getenv("ProgramFiles"))
    for program_files_path in program_files_paths:
        if os.path.isdir(program_files_path):
            for program_directory in os.listdir(program_files_path):
                match_object = re.match("^Modeller[0-9]{1,2}[.v][0-9]+$", program_directory)
                if match_object and os.path.isdir(os.path.join(program_files_path, program_directory, "modlib")):
                    return os.path.join(program_files_path, program_directory, "modlib")


#####################################################################
# PyMOL behaviour.                                                  #
#####################################################################

def check_pymol_builtin_cealign():
    """
    Try to see if PyMOL 'cealign' command can be used within PyMod.
    """
    has_builtin_cealign = False
    if hasattr(cmd,"cealign"):
        # Some old PyMOL builds have the 'cealign' method, but not they don't have its 'object'
        # argument, which is needed to make it work in PyMod.
        import inspect
        if "object" in inspect.getargspec(cmd.cealign)[0]:
            has_builtin_cealign = True
    return has_builtin_cealign


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


def get_askopenfilename_tuple(askopenfilename_result):
    """
    Some Tkinter versions reutrn unicode strings when using the 'askopenfilename' dialog to selection
    multiple files. Depending on the presence of '{' or '}' characters in the file name, it they
    will return an unicode string like this:
        'C:/Users/username/pymod/projects/Q6P988.fasta {C:/Users/username/pymod/projects/filename with spaces .fasta} C:/Users/username/pymod/projects/filename\ with\ parenthesis\{.fasta'
    This method will convert these kind of strings in tuples with containing the file names.
    """
    if isinstance(askopenfilename_result, tuple):
        return askopenfilename_result
    elif isinstance(askopenfilename_result, str) or isinstance(askopenfilename_result, unicode):
        substituted_string = re.sub(" \{([A-Z]:/)", " \g<1>", askopenfilename_result).lstrip("{")
        substituted_string = re.sub("\} ([A-Z]:/)", " \g<1>", substituted_string).rstrip("}")
        result_list = re.split("([A-Z]:/)", substituted_string)
        return tuple([disk+path.rstrip(" ").replace("\\","") for disk, path in zip(result_list[1::2], result_list[2::2])])


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
# BLAST databases.                                                  #
#####################################################################

def get_blast_database_prefix(dbpath):
    database_files = filter(lambda f: f != ".DS_Store", os.listdir(dbpath))
    return os.path.commonprefix(database_files)[:-1]


def verify_valid_blast_dbdir(dbpath):
    """
    Checks if the folder specified in 'dbpath' contains a valid set of sequence database files.
    """
    dbprefix = get_blast_database_prefix(dbpath)
    if dbprefix == "":
        return False
    else:
        return True


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

os_names_dict = {"win32": "Windows", "darwin": "Mac OS X", "linux2": "Linux"}

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


def check_dll(executable=''):
    '''Check DLL directly loaded by executable'''
    if not os.path.isfile(executable):
        return []
    exe_txt=open(executable,'rU').read()
    dll_pat=re.compile('\w*\.[Dd][Ll][Ll]\W')
    return sorted(set(dll_pat.findall(exe_txt)))


def check_msvc(executable=''):
    '''Detect Microsoft Visual C++ Runtime of executable'''
    dll_list=check_dll(executable)
    msvc_pat=re.compile('MSVC[MPR]\d*\.DLL')
    msvc_set=set([d[5:][:-5] for d in dll_list if msvc_pat.match(d.upper())])
    if msvc_set:
        return sorted(msvc_set)
    else:
        return []


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


def check_old_macos_version():
    """
    Checks if the Mac OS X version is an old one. Old versions do not support recent MODELLER
    versions, and this will be used to check which version of MODELLER to fetch from MODELLER
    website.
    """
    mac_ver = get_mac_ver()
    for v in ('10.3', '10.4', '10.5'):
        if mac_ver.startswith(v):
            return True
    return False


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

def check_network_connection(remote_address):
    try:
        response=urllib2.urlopen(remote_address)
        return True
    except urllib2.URLError:
        return False
