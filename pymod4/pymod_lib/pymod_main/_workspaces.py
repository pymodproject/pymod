# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
PyMod sessions managment.
"""

import os
import shutil
import pickle
import zipfile

from pymod_lib.pymod_gui.shared_gui_components_qt import askopenfile_qt

from pymol import cmd

import pymod_lib.pymod_vars as pmdt


class PyMod_workspaces:

    def start_new_session(self):
        # Cleans the main window.
        self.clear_main_window()
        # Reinitializes all PyMod elements and counters.
        self.initialize_pymod_elements_information()
        # Build a new project directory.
        self.new_job_state(overwrite=True)


    def clear_main_window(self):
        # Cleans the main window.
        for e in self.get_pymod_elements_list():
            self.main_window.delete_element_widgets(e)


    ##################
    # Save sessions. #
    ##################

    def save_pymod_session(self, project_arc_full_path):

        try:
            # Builds a temporary directory in which to store the project files and which will become zipped.
            project_temp_dirpath = os.path.join(self.current_pymod_dirpath, self.session_temp_dirname)

            if os.path.isdir(project_temp_dirpath):
                shutil.rmtree(project_temp_dirpath)
            os.mkdir(project_temp_dirpath)


            # Saves a pickle file with the information about the PyMod session. This will remove
            # PyQt and Tkinter objects stored in the 'PyMod_element' objects, because they can't be pickled.
            project_pickle_filepath = os.path.join(project_temp_dirpath, "%s.pkl" % self.session_filename)

            # Temporarily removes the PyQt objects from the main PyMod object so that it can be
            # pickled.
            _app = self.app
            _main_window = self.main_window
            self.app = None
            self.main_window = None
            window_dict = {}
            for attr_name in dir(self):
                attr_obj = getattr(self, attr_name)
                # PyMod windows have the 'is_pymod_window' attribute.
                if hasattr(attr_obj, "is_pymod_window"):
                    window_dict[attr_name] = attr_obj
                    setattr(self, attr_name, None)

            with open(project_pickle_filepath, "wb") as a_fh:
                pickle.dump(self, a_fh)

            # Restores the GUI objects.
            self.app = _app
            self.main_window = _main_window
            for attr_name in window_dict:
                setattr(self, attr_name, window_dict[attr_name])


            # Saves a PyMOL session.
            cmd.save(os.path.join(project_temp_dirpath, "%s.pse" % self.session_filename))


            # Copies the current project files in the directory.
            src = self.current_project_dirpath
            dst = os.path.join(project_temp_dirpath, self.session_filename)
            shutil.copytree(src, dst)


            # Builds a .zip file of the temporary directory.
            src = project_temp_dirpath
            zpf = project_arc_full_path # os.path.join(project_temp_dirpath, project_arc_name)
            shutil.make_archive(zpf, 'zip', src)

            # Removes the .zip extension which is added by the 'make_archive' method.
            if os.path.isfile(zpf + ".zip"):
                shutil.move(zpf + ".zip", zpf)
            # Finally removes the temporary directory.
            shutil.rmtree(project_temp_dirpath)

        except Exception as e:
            if os.path.isdir(project_temp_dirpath):
                shutil.rmtree(project_temp_dirpath)
            title = "Save Project Error"
            message = "Could not save the project file to path: %s" % (project_arc_full_path)
            if self.main_window is None:
                self.main_window = _main_window
            self.main_window.show_error_message(title, message)

    ###################
    # Loads sessions. #
    ###################

    def open_pymod_session(self, project_archive_filepath):

        # Where the archive is unpacked, and where the project of this session
        # lives. The second is remembered before the file is read: opening a
        # session replaces everything this object holds, the paths of the
        # machine it was saved on included, and the project has to come back
        # here rather than wherever it was saved.
        project_temp_dirpath = os.path.join(self.current_pymod_dirpath,
                                            self.session_temp_dirname)
        project_dirpath_in_use = self.current_project_dirpath

        # If some errors happen here, continue the PyMod session.
        try:

            # Check if the file is a valid zip file.
            if not zipfile.is_zipfile(project_archive_filepath):
                raise Exception("The file is not a zip file.")

            # Check if the file is a valid PyMod session file.
            zfh = open(project_archive_filepath, 'rb')
            zipfile_obj = zipfile.ZipFile(zfh)
            files_to_check = ["%s.pkl" % self.session_filename, "%s.pse" % self.session_filename]
            if not set(files_to_check) < set(zipfile_obj.namelist()):
                zfh.close()
                raise Exception("The file is not a valid PyMod session file.")
            zfh.close()

            # Builds a temporary directory in which to store project files.
            if os.path.isdir(project_temp_dirpath):
                shutil.rmtree(project_temp_dirpath)
            os.mkdir(project_temp_dirpath)

            # Extract the file to a temporary directory.
            shutil.unpack_archive(project_archive_filepath, project_temp_dirpath, format="zip")

        except Exception as e:
            self.load_project_failure(project_archive_filepath, e, project_temp_dirpath, continue_session=True)
            return None


        # If some errors happens here, close PyMod.
        try:

            # Clears the main window.
            self.clear_main_window()

            # Unpickle the data of the saved PyMod project.
            project_pickle_filepath = os.path.join(project_temp_dirpath, "%s.pkl" % self.session_filename)

            _app = self.app
            _main_window = self.main_window

            a_fh = open(project_pickle_filepath, "rb")
            self.__dict__ = pickle.load(a_fh).__dict__
            a_fh.close()

            self.app = _app
            self.main_window = _main_window

            # Load the PyMOL session.
            pymol_session_filepath = os.path.join(project_temp_dirpath, "%s.pse" % self.session_filename)
            cmd.load(pymol_session_filepath)

            # Ensure af_modeling_session_list exists for compatibility with older sessions
            if not hasattr(self, 'af_modeling_session_list'):
                self.af_modeling_session_list = {}

            # Ensure esm_modeling_session_list exists for compatibility with older sessions
            if not hasattr(self, 'esm_modeling_session_list'):
                self.esm_modeling_session_list = []

            # Ensure omegafold_modeling_session_list exists for compatibility with older sessions
            if not hasattr(self, 'omegafold_modeling_session_list'):
                self.omegafold_modeling_session_list = []

            # Ensure performed_boltz_modeling_count exists for compatibility with older sessions
            if not hasattr(self, 'performed_boltz_modeling_count'):
                self.performed_boltz_modeling_count = 0

            # The project of the archive takes the place of the one in use. This
            # happens before the paths are worked out, because they are worked
            # out from where the project ends up: pointing them at the directory
            # the archive was unpacked into left PyMod with its project inside a
            # temporary directory, which the lines below then removed --
            # "FileNotFoundError: session_temp_dir/pymod_session", and no
            # project directory at all afterwards.
            if os.path.isdir(project_dirpath_in_use):
                shutil.rmtree(project_dirpath_in_use)
            shutil.move(os.path.join(project_temp_dirpath, self.session_filename),
                        project_dirpath_in_use)

            self.current_project_dirpath = project_dirpath_in_use
            self.structures_dirpath = os.path.join(self.current_project_dirpath, self.structures_dirname)
            self.alignments_dirpath = os.path.join(self.current_project_dirpath, self.alignments_dirname)
            self.images_dirpath = os.path.join(self.current_project_dirpath, self.images_dirname)
            self.models_dirpath = os.path.join(self.current_project_dirpath, self.models_dirname)
            self.psipred_dirpath = os.path.join(self.current_project_dirpath, self.psipred_dirname)
            self.similarity_searches_dirpath = os.path.join(self.current_project_dirpath, self.similarity_searches_dirname)
            self.temp_directory_dirpath = os.path.join(self.current_project_dirpath, self.temp_directory_name)
            self.domain_analysis_dirpath = os.path.join(self.current_project_dirpath, self.domain_analysis_dirname)

            # Reinitializes the GUI of the PyMod elements.
            for el in self.get_pymod_elements_list():
                # Initializes in PyMod and add the GUI swidgets.
                el.initialize(self)
                # Updates the structure paths.
                if el.has_structure():
                    el.update_all_structure_paths(self.structures_dirpath)
                # Updates the alignment files (such as trees and dendrograms).
                if el.is_cluster():
                    el.update_alignment_files(self.alignments_dirpath)

################Add ETE code 21/05/2026-SARA################
            # Loads a PyMOL session. The file is still in the temporary
            # directory: only the project folder beside it has been moved.
            pymol_session_filepath = os.path.join(project_temp_dirpath,  "%s.pse" % self.session_filename)
            cmd.reinitialize()
            cmd.load(pymol_session_filepath)
            # Removes the temporary session directory, and works from the
            # project as every other command of PyMod does.
            shutil.rmtree(project_temp_dirpath)
            os.chdir(self.current_project_dirpath)
            # Updates PyMod main window.
            self.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True, update_elements=True)
        except Exception as e:
            # 'raise e' stood here, before the call which tells the user: the
            # exception went to the shell and the window said nothing at all.
            self.load_project_failure(project_archive_filepath, e,
                                      project_temp_dirpath,
                                      continue_session=False)
######################END ETE#######################

    def load_project_failure(self, project_archive_filepath, error, project_temp_dirpath, continue_session=True):
        if os.path.isdir(project_temp_dirpath):
            shutil.rmtree(project_temp_dirpath)
        title = "Open Project Error"
        message = "Could not open the project file '%s': because of the following error: %s." % (project_archive_filepath, error)
        if not continue_session:
            message += " PyMod is now shutting down."
        self.main_window.show_error_message(title, message)
        if not continue_session:
            self.main_window.destroy()
