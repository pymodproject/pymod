# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import gzip
import shutil
import subprocess
import time
import traceback
from ftplib import FTP

from pymol.Qt import QtCore

from pymod_lib.pymod_vars import blast_databases_dirname, hmmer_databases_dirname, hmmscan_databases_dirname
from pymod_lib.pymod_os_specific import get_formatted_date


##################################################################################################
# Supporting classes and functions.                                                              #
##################################################################################################

class PyMod_component:

    def __init__(self, name, full_name, subdir_name=None, remote_source=None, databases_string=None):
        self.name = name
        self.full_name = full_name
        self.target_installation_path = None
        self.remote_source = remote_source
        self.installer = None
        self.subdir_name = subdir_name
        self.can_be_downloaded = False
        self.databases_string = databases_string
        self.last_downloaded = None

    def reinitialize(self):
        self.target_installation_path = None
        self.can_be_downloaded = False
        self.last_downloaded = None


def gunzip_shutil(source_filepath, dest_dirpath, block_size=65536):
    """
    Support for gzipped files.
    """
    dest_filename = os.path.basename(source_filepath).replace('.gz', '')
    dest_filepath = os.path.join(dest_dirpath, dest_filename)
    with gzip.open(source_filepath, 'rb') as s_file, open(dest_filepath, 'wb') as d_file:
        shutil.copyfileobj(s_file, d_file, block_size)

def is_compressed(file_name):
    return file_name.endswith(compressed_files_extensions)


##################################################################################################
# DEFAULT VARIABLES                                                                              #
##################################################################################################

# Installer files information.
compressed_files_extensions = ('.tar.gz', '.gz', '.tgz', '.zip', '.tar', '.tar.xz', '.txz', '.tar.bz2', '.tbz2')


# Lists of PyMod databases that can be installed through the PyMod installer GUI.
all_components_list = [
    # BLAST+ v4. # Old BLAST+ versions (< 2.8.0) only recognize "version 4" databases
    # (found in the "v4" directory of the BLAST ftp site). The default directory,
    # contains "version 5" databases, which are only compatible with newer BLAST+
    # versions.
    PyMod_component(name=blast_databases_dirname, full_name="BLAST+ v4 databases",
                    remote_source="ftp://ftp.ncbi.nlm.nih.gov/blast/db/v4/",
                    subdir_name=["swissprot_v4", "pdbaa_v4"],
                    databases_string="swissprot, pdbaa"),
    # BLAST+ v5.
    PyMod_component(name=blast_databases_dirname, full_name="BLAST+ databases",
                    remote_source="ftp://ftp.ncbi.nlm.nih.gov/blast/db/",
                    subdir_name=["swissprot", "pdbaa"],
                    databases_string="swissprot, pdbaa"),
    # PFAM.
    PyMod_component(name=hmmscan_databases_dirname, full_name="HMMSCAN database",
                    remote_source="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/",
                    databases_string="PFAM"),
    # FASTA files for HMMER programs.
    PyMod_component(name=hmmer_databases_dirname, full_name="HMMER databases",
                    remote_source="ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/",
                    subdir_name=["swissprot", "pdbaa"],
                    databases_string="swissprot, pdbaa")
    ]


###################################################################################################
# PYMOD COMPONENT INSTALLER CLASSES                                                               #
###################################################################################################

class InstallerQThread(QtCore.QThread):
    """
    A QThread wrapper for the installers.
    Contains methods that make sense only when this class is extended
    together with a PyMod_component_installer (or sub) class, creating a subclass
    that inherits both from this class and PyMod_component_installer.
    """

    critical_error = QtCore.pyqtSignal(str, PyMod_component)
    info_message = QtCore.pyqtSignal(str, PyMod_component)
    established_connection = QtCore.pyqtSignal(PyMod_component)

    retrieved_size = QtCore.pyqtSignal(PyMod_component, int)
    set_update_status = QtCore.pyqtSignal(PyMod_component, str, str)

    terminated_installation = QtCore.pyqtSignal(PyMod_component)
    freeze_thread = QtCore.pyqtSignal(PyMod_component, int)

    def run(self):
        """this method is never called alone (see QThread documentation)
        but is executed when the Qthread.start() method is called
        if the flag 'install_mode' is True, the thread will install the databases,
        if not, it will only retrieve the file size from the server."""

        # Only pings the remote source.
        if self.install_mode == "ping":
            # Returns 'False' if there is a connection error.
            time.sleep(0.5)
            connection = self.ping()
            if connection:
                self.component.can_be_downloaded = True

        else:
            try:
                # Actually performs the installation.
                self.download_and_install()

            except TerminateQThread as e:
                pass

            # emette un type error, oltre al gaierror, se non c'e' connessione internet.
            # Emette anche un EOFError e un TimeoutError se non fa in tempo a scaricare.
            except (gaierror, TypeError) as e:
                self.critical_error.emit("Cannot connect to server. Please check Internet connection.", self.component)
            except EOFError as e:
                self.critical_error.emit("Timeout expired for connection. Try again later.", self.component)
            except Exception as e:
                msg = str(e)
                self.critical_error.emit(msg, self.component)
                traceback.print_exc()


    ##########################################################
    # Wrapper for the signals. Used in installer subclasses. #
    ##########################################################

    def emit_signal(self, signal_name, *args, **kwargs):
        # getattr(self, signal_name).emit(*args, **kwargs)
        if signal_name == "set_update_status":
            self.set_update_status.emit(*args, **kwargs)
        elif signal_name == "retrieved_size":
            self.retrieved_size.emit(*args, **kwargs)
        elif signal_name == "critical_error":
            self.critical_error.emit(*args, **kwargs)
        elif signal_name == "terminated_installation":
            self.terminated_installation.emit(*args, **kwargs)
        else:
            raise KeyError("Unknown 'signal_name': %s" % signal_name)


class PyMod_component_installer(InstallerQThread):
    """The installer class must be extended together with the InstallerQThread class. Use the
    'build_installer_qthread' function in 'updater_gui' module. """

    def __init__(self, component, download_destination_dirpath=''):
        InstallerQThread.__init__(self)
        self.component = component
        self.component.installer = self
        self.installation_status = None
        self.download_destination_dirpath = download_destination_dirpath
        self.download_destination_filepath = ''
        self.downloaded_filepath_list = []
        # Initializes the 'install_mode'.
        self.install_mode = None

    def has_zipped_files(self):
        if [fp for fp in self.downloaded_filepath_list if is_compressed(fp)]:
            return True
        else:
            return False

    def download(self):
        '''This method performs the installation and must set the installation_status flag as True or False.
        It has to be overridden in subclasses.'''
        raise NotImplementedError


    def unzip_downloaded_files(self, destination_dirpath,
                               explicit_list_of_files=None, in_thread=True):
        self._unzip_downloaded_files(destination_dirpath=destination_dirpath,
                                     explicit_list_of_files=explicit_list_of_files,
                                     in_thread=in_thread)


    def _unzip_downloaded_files(self, destination_dirpath,
                                explicit_list_of_files=None, in_thread=True):
        # Shutil is the best choice here, bc calling tarfile or zipfile separately for each format is
        # error-prone, while shutil.unpack_archive handles every case on its own.
        # However, shutil does not support simple '.gz' compressed files. I have to register an unpack format,
        # named "bio_gz". Then, I provide a list of extensions corresponding to the format.
        # Cannot put '.gz' directly, because it may overlap with the .tar.gz format and I don't want to mess
        # with builtin modules.
        # Then, the method requires a callable that will be used to unpack archives.
        # The callable must receive the path of the archive, followed by the directory
        # the archive must be extracted to. This callable is created in this module, is the reimplementation of
        # the gunzip command, called gunzip_shutil.
        try:
            shutil.register_unpack_format("bio_gz", [".fasta.gz", ".hmm.gz", ".gz"], gunzip_shutil)
        except shutil.RegistryError: # if it is already registered, ignore it.
            pass

        # selecting compressed files
        if not explicit_list_of_files:
            zip_files = [fp for fp in self.downloaded_filepath_list if is_compressed(fp)]
        else:
            zip_files = explicit_list_of_files

        # decompressing
        for element in zip_files:

            if not os.path.exists(destination_dirpath):
                os.makedirs(destination_dirpath)
            # In Python 3.8, unpacking a tar.gz file in this QThread causes a
            # segmentation fault. They will be upacked later in the main thread.
            if in_thread and element.endswith(".tar.gz"):
                continue

            shutil.unpack_archive(element, destination_dirpath)
            os.remove(element)     # cleaning compressed files
            self.downloaded_filepath_list.remove(element)


    def finish_generic_installation(self):
        """Concludes an installation with facultative additional actions,
        then sets the installation_status flag as True and the
        component installed path."""

        self.additional_actions()
        # installation successful
        self.installation_status = "success"


    def additional_actions(self, *args, **kwargs):
        """Additional actions to perform after the download and the unzipping,
        in order to conclude the installation, e.g. the HMM database preparation with hmmpress"""
        # to be overridden.
        pass


    def build_destination_filepath(self, filename):
        new_filepath = os.path.join(self.download_destination_dirpath, filename)
        self.downloaded_filepath_list.append(new_filepath)
        self.download_destination_filepath = new_filepath
        return new_filepath


class PyMod_FTP_component_installer(PyMod_component_installer):

    def __init__(self, component, destination_dirpath=''):
        PyMod_component_installer.__init__(self, component, destination_dirpath)
        self.hostname, self.subdir = self.get_host_and_subdir_from_component_url()

    def get_host_and_subdir_from_component_url(self):
        '''From the component.remote_source URL, retrieves the host name (es. "ftp.ebi.ac.uk")
        and the subdirectory path of the tool (es. "./pub/databases/Pfam/current_release/database_files/")'''
        remoteurl = self.component.remote_source
        hostname, subdir = remoteurl.replace('ftp://', '').split('/', 1)
        return hostname, subdir


    def ping(self):
        """this method checks for the server and retrieves the size of the compressed file(s)"""
        try:
            ftp = FTP(self.hostname, timeout=10)
            ftp.login()
            if self.subdir:
                ftp.cwd(self.subdir)

            files_on_server = self.get_files_on_server(ftp)
            item_basenames = self.get_item_basenames(files_on_server)
            self.set_total_size(ftp, item_basenames) # sets the self.total_size attribute
            ftp.quit()
            return True

        except Exception as e:
            # traceback.print_exc()
            self.emit_signal("set_update_status", self.component, "Connection Error (%s)" % e, "red")
            return False


    def ftp_connect(self, host, subdirectory_path=None):
        """this method connects to the FTP server and downloads the files, handling signals and
        connection closing. Instruction to download each specific database file must be written in
        'retrieve_files' method."""

        ftp = FTP(host, timeout=30)
        ftp.login()
        if subdirectory_path:
            ftp.cwd(subdirectory_path)

        # emitting the 'set_update_status' signal to show in the table that the donwload has started.
        self.emit_signal("set_update_status", self.component, "Downloading...", "light_green")

        # this method actually downloads the files and must be implemented in subclasses
        self.retrieve_files(ftp)

        ftp.quit()


    def set_total_size(self, ftp_connection, filenames):
        total_size = 0
        ftp_connection.sendcmd("TYPE i") # Switch to Binary mode
        total_size = sum([ftp_connection.size(filename) for filename in filenames])
        ftp_connection.sendcmd("TYPE A") # Switch back to ASCII mode
        self.total_size = total_size
        self.emit_signal("retrieved_size", self.component, total_size)
        return total_size


    def retrieve_files(self, ftp_connection):
        '''this method must provide the code for the download of the correct files and return the
        downloaded file path. It has to be overridden in subclasses.'''
        raise NotImplementedError


    def download_and_install(self):
        """
        This is the method that is called externally. It calls 'ftp_connect' and handles its
        exceptions. After the download has finished, it proceeds to prepare the database files for
        PyMod (that is, it proceeds to the "installation" of the database).
        """

        # Downloads the component files.
        try:
            self.ftp_connect(self.hostname, subdirectory_path=self.subdir)
            # Once the download has been completed, updates the GUI.
            self.installation_status = "downloaded"
            self.emit_signal("set_update_status", self.component, "Downloaded", "light_green")

        except Exception as e:
            # traceback.print_exc()
            messagestr = 'Error during the download of %s via FTP: %s' % (self.component.full_name, e)
            self.emit_signal("critical_error", messagestr, self.component)
            return None

        # "Installs" the component files.
        try:
            # Unzips the files.
            if self.has_zipped_files():
                self.emit_signal("set_update_status", self.component, "Unpacking...", "light_green")
                self.unzip_downloaded_files(self.component.target_installation_path)
            # Performs additional actions.
            if self.installation_status == "stopped":
                return None
            self.finish_generic_installation()

            # Set the last downloaded time.
            self.component.last_downloaded = get_formatted_date()
            self.emit_signal("set_update_status", self.component, "Completed", "green")
            self.emit_signal("terminated_installation", self.component)

        except Exception as e:
            # traceback.print_exc()
            messagestr = 'Error during the installation of %s: %s' % (self.component.full_name, e)
            self.emit_signal("critical_error", messagestr, self.component)


# subclasses
class Pfam_db_installer(PyMod_FTP_component_installer):

    pfam_hmm_filename = "Pfam-A.hmm"

    def __init__(self, component, destination_dirpath=''):
        PyMod_FTP_component_installer.__init__(self, component, destination_dirpath)
        self.hmmpress_exe_filepath = None


    def get_files_on_server(self, ftp_connection):
        listdir = []
        ftp_connection.dir(listdir.append)
        return listdir

    def get_item_basenames(self, listdir):
        """
        Searches for the 'Pfam-A.hmm.gz' file in the ftp directory.
        """
        for item in listdir:
            item_basename = item.split()[-1]
            if item_basename == self.pfam_hmm_filename + '.gz':
                return [item_basename]
        raise ValueError("HMM file not found on PFAM server.")

    def retrieve_files(self, ftp_connection):
        listdir = self.get_files_on_server(ftp_connection)
        item_basename = self.get_item_basenames(listdir)
        download_ftp_with_batches(ftp_connection=ftp_connection,
                                  ftp_source='RETR ' + item_basename[0],
                                  dest_filepath=self.build_destination_filepath(item_basename[0]),
                                  installer=self)

    def additional_actions(self):
        # message = ("The Pfam database has been successfully downloaded. It will"
        #            " now be decompressed using the HMMPRESS tool. Please wait as"
        #            " it might take a little while.")
        self.emit_signal("set_update_status", self.component, "Running hmmpress...", "light_green")

        # Remove the databases files. Their presence will prevent hmmpress to build the new database.
        for ext in (".h3m", ".h3i", ".h3f", ".h3p"):
            db_filepath = os.path.join(self.component.target_installation_path, self.pfam_hmm_filename + ext)
            if os.path.isfile(db_filepath):
                os.remove(db_filepath)

        # Actually launches hmmpress to complete the database installation.
        database_filepath = os.path.join(self.component.target_installation_path, self.pfam_hmm_filename)
        try:
            assert os.path.exists(self.hmmpress_exe_filepath)
            extract_hmmer_databases(database_filepath, self.hmmpress_exe_filepath)
            #little loop to ensure that ALL the files are present in the directory
            #they are big and the following command must not be executed if the extraction is incomplete
            extracted_files = [i for i in os.listdir(self.component.target_installation_path) if not i.endswith('.hmm') and not i.startswith('.')]
            while len(extracted_files) < 4:
                time.sleep(0.5)
            os.remove(database_filepath)

        except:
            # traceback.print_exc()
            message = ("Error in decompressing Pfam databases with HMMPRESS.")
            raise Exception(message)


def extract_hmmer_databases(hmmer_data_filepath, hmmpress_exepath):
    """Indexing HMM databases with hmmpress"""
    cline = [hmmpress_exepath, hmmer_data_filepath]
    print("- Executing the following command:", cline)
    subprocess.check_call(cline)


class BLAST_db_installer(PyMod_FTP_component_installer):

    def __init__(self, component, destination_dirpath=''):
        PyMod_FTP_component_installer.__init__(self, component, destination_dirpath)


    def get_files_on_server(self, ftp_connection):
        listdir = []
        ftp_connection.retrlines('MLSD', listdir.append)
        listdir = [i for i in listdir if not i.endswith('md5')]  # throw away the md5 checksumm files
        self.blastdbs = []
        for name_of_db in self.component.subdir_name:
            self.blastdbs.extend([i for i in listdir if name_of_db in i])
        return self.blastdbs

    def get_item_basenames(self, listdir):
        if self.blastdbs:
            item_basenames = [item.split(';')[-1].strip() for item in self.blastdbs]
            return item_basenames
        else:
            return []

    def retrieve_files(self, ftp_connection):
        self.get_files_on_server(ftp_connection)
        item_basenames = self.get_item_basenames(self.blastdbs)
        if item_basenames:
            for item_basename in item_basenames:
                download_ftp_with_batches(ftp_connection=ftp_connection,
                                          ftp_source='RETR ' + item_basename,
                                          dest_filepath=self.build_destination_filepath(item_basename),
                                          installer=self)

        else:
            raise NotImplementedError
        del self.blastdbs


    def unzip_downloaded_files(self, base_destination_dirpath,
                               explicit_list=None, in_thread=True):

        # selecting compressed files
        zip_files = [fp for fp in self.downloaded_filepath_list if is_compressed(fp)]
        zip_dict = {}
        for subd_name in self.component.subdir_name:
            zip_dict.update({subd_name: [f for f in zip_files if subd_name in f]})

        for k in zip_dict.keys():
            new_subdir = os.path.join(base_destination_dirpath, k)

            self._unzip_downloaded_files(new_subdir, zip_dict[k],
                                         in_thread=in_thread)

        # adjustment of files
        # another for loop bc the unzipping MUST be finished, if included in the previous for loop
        # raises errors and neglects some files
        for k in zip_dict.keys():
            new_subdir = os.path.join(base_destination_dirpath, k)
            content = [l for l in os.listdir(new_subdir) if not l.startswith('.')]
            tab_files = [f for f in content if f.startswith('tax')]
            for tf in tab_files:
                os.remove(os.path.join(new_subdir, tf))

            if self.component.name != 'blast_databases':
                for file_to_move in content:
                    dest_filepath = os.path.join(base_destination_dirpath, file_to_move+'.fasta')
                    if os.path.isfile(dest_filepath):
                        os.remove(dest_filepath)
                    os.rename(os.path.join(new_subdir, file_to_move), dest_filepath)
                os.removedirs(new_subdir)


def download_ftp_with_batches(ftp_connection, ftp_source, dest_filepath, installer):
    """
    Downloads an FTP file, and for each downloaded batch checks whether the thread has stopped if
    it has, the download is also stopped. This system allows to quit a download thread without
    abruptly terminating the thread (which might freeze the entire application).
    """
    with open(dest_filepath, 'wb') as fh:
        def ftp_callback(batch_data):
            if installer.installation_status != "stopped":
                fh.write(batch_data)
            else:
                raise TerminateQThread("Stop download")
        resp = ftp_connection.retrbinary(ftp_source, callback=ftp_callback)


class TerminateQThread(Exception):
    pass
