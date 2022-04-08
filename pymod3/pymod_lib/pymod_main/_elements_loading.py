# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing a series of methods which allow to load PyMod elements in PyMod.
"""

import os
import re

from Bio import SeqIO

from pymol import cmd

from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt, askopenfile_qt
from pymod_lib.pymod_gui.specific_gui_components_qt import (Search_string_window_qt,
                                                            Add_feature_window_qt,
                                                            Edit_sequence_window_qt)

from pymod_lib import pymod_structure
from pymod_lib import pymod_vars
from pymod_lib.pymod_protocols import structural_databases_protocols
from pymod_lib.pymod_seq.seq_manipulation import check_correct_sequence
from pymod_lib.pymod_element_feature import Element_feature

from pymod_lib.pymod_exceptions import PyModInvalidFile


class PyMod_elements_loading:

    ###############################################################################################
    # FILES MANAGMENT.                                                                            #
    ###############################################################################################

    #################################################################
    # Check correct files formats.                                  #
    #################################################################

    def is_sequence_file(self, file_path, file_format, show_error=True):
        """
        Try to open a sequence file using Biopython. Returns 'True' if the file is a valid file of
        the format specified in the 'file_format' argument.
        """
        valid_file = False
        file_handler = None
        try:
            file_handler = open(file_path,"r")
            r = list(SeqIO.parse(file_handler, file_format))
            if len(r) > 0:
                valid_file = True
        except:
            valid_file = False
        if file_handler != None:
            file_handler.close()
        return valid_file


    def is_valid_structure_file(self,file_name, format="pdb", show_error=True):
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
            self.main_window.show_error_message(title,message)
        return valid_pdb


    #################################################################
    # Load sequence files.                                          #
    #################################################################

    def open_sequence_file(self, file_full_path, file_format="fasta"):
        """
        Method for loading in PyMod new sequences parsed from sequence files. It will build new
        PyMod elements, but it will not display its widgets in the main window.
        """
        if not os.path.isfile(file_full_path):
            raise IOError("File does not exist: %s." % file_full_path)
        if not self.is_sequence_file(file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file: %s." % (file_format, file_full_path))

        # Parses a sequence file through Biopython. This will automatically crop headers that have
        # " " (space) characters.
        elements_to_return = []
        for record in SeqIO.parse(file_full_path, file_format):
            # Then builds a PyMod_element object and add it to the 'pymod_elements_list'.
            c = self.build_pymod_element_from_seqrecord(record)
            e = self.add_element_to_pymod(c)
            elements_to_return.append(e)
        return elements_to_return


    def build_cluster_from_alignment_file(self, alignment_file, extension="fasta"):
        """
        Creates a cluster with all the sequences contained in an alignment file.
        """
        # Gets the sequences using Biopython.
        aligned_elements = []
        records = SeqIO.parse(alignment_file, extension)
        for record in records:
            new_child_element = self.build_pymod_element_from_seqrecord(record)
            self.add_element_to_pymod(new_child_element)
            aligned_elements.append(new_child_element)
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="alignment", child_elements=aligned_elements, algorithm="imported")
        return new_cluster


    #################################################################
    # Opening PDB files.                                            #
    #################################################################

    def open_structure_file(self, pdb_file_full_path, file_format="pdb"):
        """
        Opens a PDB file (specified in 'pdb_file_full_path'), reads its content, imports in PyMod
        the sequences of the polypeptide chains and loads in PyMOL their 3D structures.
        """
        if not self.is_valid_structure_file(pdb_file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file." % file_format)
        p = pymod_structure.Parsed_pdb_file(self, pdb_file_full_path, output_directory=self.structures_dirpath)
        elements_to_return = []
        for element in p.get_pymod_elements():
            e = self.add_element_to_pymod(element)
            elements_to_return.append(e)

        # Renames the full PDB file if needed.
        original_pdb_filename = p.structure_file_name
        if original_pdb_filename in self.original_pdb_files_set:
            counter = 1
            while original_pdb_filename in self.original_pdb_files_set:
                original_pdb_filename = "%s_%s" % (counter, original_pdb_filename)
                counter += 1
            for e in elements_to_return:
                e.rename_file_name_root(original_pdb_filename)
        self.original_pdb_files_set.add(original_pdb_filename)

        return elements_to_return


    def color_struct(self):
        color_to_return = pymod_vars.pymod_regular_colors_list[self.color_index % len(pymod_vars.pymod_regular_colors_list)]
        self.color_index += 1
        return color_to_return


    #################################################################
    # Open files dialogs from PyMod.                                #
    #################################################################

    def choose_alignment_file(self):
        """
        Lets users choose an alignment file.
        """
        # Creates a PyQt widget that lets the user select multiple files.
        alignment_file_path = askopenfile_qt("Open an alignment file",
                                             name_filter="*fasta *aln *clu *sto *sth",
                                             parent=self.get_qt_parent())
        if not alignment_file_path:
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(alignment_file_path)[1].replace(".","")

        if extension == "fasta":
            pass
        elif extension in ("aln", "clu"):
            extension = "clustal"
        elif extension in ("sto", "sth"):
            extension = "stockholm"
        # Unknown format.
        else:
            title = "Format Error"
            message = "Unknown alignment file format: %s" % (extension)
            self.main_window.show_error_message(title, message)
            return (None, None)
        return alignment_file_path, extension


    def choose_structure_file(self):
        """
        Lets users choose a strcture file.
        """
        # Creates a PyQt widget that lets the user select multiple files.
        open_filepath = askopenfile_qt("Open an alignment file",
                                       name_filter="*pdb *ent",
                                       parent=self.get_qt_parent())
        if open_filepath == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(os.path.basename(open_filepath))[1].replace(".","")
        return open_filepath, extension


    ###############################################################################################
    # EDIT SEQUENCE AND STRUCTURES.                                                               #
    ###############################################################################################

    def show_edit_sequence_window(self, pymod_element):
        """
        Edit a sequence.
        """
        self.edit_sequence_window = Edit_sequence_window_qt(self.main_window,
                                                            pymod_element=pymod_element,
                                                            title="Edit Sequence",
                                                            upper_frame_title="Edit your Sequence",
                                                            submit_command=self.edit_sequence_window_state)
        self.edit_sequence_window.show()
        # self.edit_sequence_window.resize(700, self.edit_sequence_window.sizeHint().height())


    def edit_sequence_window_state(self):
        """
        Accept the new sequence.
        """
        edited_sequence = self.edit_sequence_window.get_sequence()
        # When editing elements with a structure loaded in PyMOL, only indels can be added/removed.
        if self.edit_sequence_window.pymod_element.has_structure():
            if self.edit_sequence_window.pymod_element.my_sequence.replace("-", "") != edited_sequence.replace("-", ""):
                message = ("The amino acid sequence of an element with a 3D structure"
                           " loaded in PyMOL can not be edited. Only indels can be"
                           " added or removed.")
                self.main_window.show_error_message("Sequence Error", message)
                return None
        if not len(edited_sequence):
            self.main_window.show_error_message("Sequence Error", "Please submit a non empty string.")
            return None
        if not check_correct_sequence(edited_sequence):
            self.main_window.show_error_message("Sequence Error", "Please provide a sequence with only standard amino acid characters.")
            return None
        self.edit_sequence_window.pymod_element.set_sequence(edited_sequence, permissive=True)
        self.main_window.gridder(update_clusters=True, update_elements=True)
        self.edit_sequence_window.destroy()


    def duplicate_sequence(self, element_to_duplicate):
        """
        Make a copy of a certain element.
        """
        if element_to_duplicate.has_structure():
            p = pymod_structure.Parsed_pdb_file(self, element_to_duplicate.get_structure_file(basename_only=False),
                                                output_directory=self.structures_dirpath,
                                                new_file_name=pymod_vars.copied_chain_name % self.new_objects_index)
            self.new_objects_index += 1
            for element in p.get_pymod_elements():
                self.add_element_to_pymod(element, color=element_to_duplicate.my_color) # Add this to use the old color shceme of PyMod: color=self.color_struct()
        else:
            #duplicated_element = self.build_pymod_element(pmel.PyMod_sequence_element, element_to_duplicate.my_sequence, element_to_duplicate.my_header_root)
            duplicated_element = self.build_pymod_element_from_args(element_to_duplicate.my_header_root, element_to_duplicate.my_sequence)
            self.add_element_to_pymod(duplicated_element)
            return duplicated_element


    #################################################################
    # Add features.                                                 #
    #################################################################

    def show_add_feature_window(self, pymod_element, selected_residue):
        """
        Edit a feature to a residue (or series of residues).
        """

        self.add_feature_window = Add_feature_window_qt(self.main_window,
                                                        pymod_element=pymod_element,
                                                        selected_residue=selected_residue,
                                                        title="Add a Feature to %s" % pymod_element.compact_header,
                                                        upper_frame_title="Add a Feature to your Sequence",
                                                        submit_command=self.add_feature_window_state)
        self.add_feature_window.show()


    def add_feature_window_state(self):
        """
        Accept the new feature.
        """

        selected_element = self.add_feature_window.pymod_element

        #-----------------------------------
        # Get the parameters from the GUI. -
        #-----------------------------------

        # Gets the residue range.
        residue_range_str = self.add_feature_window.get_residue_range()

        selection_warning_title = "Selection Warning"
        examples_string = "'34' for a single residue, '12-20' for a residue range"
        selection_warning_message = "Please provide a valid string to select one ore more residues (use the format: %s)." % examples_string

        # Checks for an empty string.
        if not residue_range_str:
            self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
            return None

        # Checks for a single residue.
        sre_match = re.search("^\d+$", residue_range_str)
        if sre_match:
            try:
                residue_range = (int(residue_range_str), int(residue_range_str))
            except ValueError:
                self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
                return None

        # Checks for a residue range.
        else:
            sre_match = re.search("^(\d+\-\d+)$", residue_range_str)
            if not sre_match:
                self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
                return None

            try:
                min_res, max_res = residue_range_str.split("-")
                residue_range = (int(min_res), int(max_res))
            except (ValueError, IndexError):
                self.main_window.show_warning_message(selection_warning_title, selection_warning_message)
                return None

            if residue_range[1] < residue_range[0]:
                self.main_window.show_warning_message(selection_warning_title,
                    "Invalid residue range. The index of the second residue in the range (%s) can not be smaller than the first one (%s)." % (residue_range[1], residue_range[0]))
                return None

        # Convert the 'db_index' values obtained from the GUI into 'seq_index' values.
        try:
            feature_start = selected_element.get_residue_by_db_index(residue_range[0]).seq_index
            feature_end = selected_element.get_residue_by_db_index(residue_range[1]).seq_index
         # The 'db_index' values do not correspond to any residue.
        except KeyError as e:
            self.main_window.show_warning_message(selection_warning_title,
                'The selected sequence does not have a residue with the following id: %s.' % e)
            return None


        # Gets the feature name.
        feature_name = self.add_feature_window.get_feature_name()
        sre_match = re.search("[^ a-zA-Z0-9_-]", feature_name)
        if not feature_name or sre_match:
            self.main_window.show_warning_message(selection_warning_title,
                'Please provide a valid "Feature Name" string (only alphanumeric characters and the "-", "_" and " " characters are allowed).')
            return None
            if len(feature_name) > 15:
                self.main_window.show_warning_message(selection_warning_title,
                    'Please provide a "Feature Name" string shorter than 15 characters.')
                return None


        # Gets the feature color.
        selected_rgb, selected_hex = self.add_feature_window.get_selected_colors()
        feature_color = self.add_new_color(selected_rgb, selected_hex)

        #------------------------------------------------
        # Actually adds the new feature to the element. -
        #------------------------------------------------

        # Adds the domain to the sequence.
        new_feature = Element_feature(id=None, name=feature_name,
                                      start=feature_start, end=feature_end,
                                      description=feature_name,
                                      feature_type='sequence',
                                      color=feature_color)
        selected_element.add_feature(new_feature)

        self.main_window.color_element_by_custom_scheme(selected_element, selected_feature=new_feature)

        if self.add_feature_window.get_select_in_pymol():
            selection_name = "%s_%s_%s" % (feature_name.replace(" ", "_"), selected_element.unique_index, selected_element.features_count)
            cmd.select(selection_name,
                       "object %s and resi %s-%s" % (selected_element.get_pymol_selector(), residue_range[0], residue_range[1]))
            selected_element.features_selectors_list.append(selection_name)

        self.add_feature_window.destroy()


    def delete_features_from_context_menu(self, pymod_element):
        for selection in pymod_element.features_selectors_list:
            try:
                cmd.delete(selection)
            except:
                pass
        pymod_element.clear_features()
        pymod_element.revert_original_color()
        self.main_window.color_element(pymod_element)


    #################################################################
    # Search subsequences.                                          #
    #################################################################

    def show_search_string_window(self, pymod_element):
        """
        Search for a sub-sequence in a protein sequence.
        """
        self.search_string_window = Search_string_window_qt(self.main_window,
                                                            pymod_elements=[pymod_element],
                                                            title="Search in %s" % pymod_element.compact_header,
                                                            upper_frame_title="Search a string in your sequence",
                                                            submit_command=self.search_string_window_state,
                                                            submit_button_text="Search")
        self.search_string_window.show()


    def search_string_window_state(self, event=None):

        selected_element = self.search_string_window.pymod_element

        #-----------------------------------
        # Get the parameters from the GUI. -
        #-----------------------------------

        # Get the string to search for in the selected sequence.
        search_string = self.search_string_window.get_search_string()
        _search_string = "" # " '" + search_string + "'"

        selection_warning_title = "Search Warning"
        # examples_string = "'34' for a single residue, '12-20' for a residue range"
        # selection_warning_message = "Please provide a valid string to select one ore more residues (use the format: %s)." % examples_string

        # Checks for an empty string.
        if not search_string:
            self.main_window.color_element(selected_element)
            self.search_string_window.show_results(self.search_string_window.default_message, state="empty")
            return None

        # Decides whether to use regex to search for a string.
        use_regex = self.search_string_window.get_regex_use()

        # Checks the search string.
        if not use_regex:
            # Search for invalid characters.
            search_sre = re.search("[^A-Za-z]", search_string)
            if search_sre:
                self.main_window.show_warning_message(selection_warning_title,
                    "Please provide a string containing only standard amino acid characters (example: 'ATGV').")
                return None

        else:
            # Test the regular expression.
            try:
                search_sre = re.search(search_string, "test_string")
            except re.error:
                self.main_window.show_warning_message(selection_warning_title,
                    "Please provide a string containing a valid regular expression.")
                return None

        # Gets the highlight color.
        highlight_color = self.search_string_window.get_highlight_color()


        #------------------------------------------------------------
        # Searches for the string in the selected protein sequence. -
        #------------------------------------------------------------

        selection_results_title = "Search Results"
        full_sequence = str(selected_element.my_sequence.replace("-", ""))

        # Actually executes the regular expression to look for sequences.
        finditer_results = list(re.finditer(search_string, full_sequence, flags=re.IGNORECASE))

        if not finditer_results:
            self.main_window.color_element(selected_element)
            results_message = "Pattern" + _search_string + " not found."
            self.search_string_window.show_results(results_message, state="not_found")
            return None

        # Builds a list of residue indices encompassing the matched strings.
        selected_residues = []
        for re_match in finditer_results:
            selected_residues.extend(list(range(re_match.start(), re_match.end())))
        selected_residues = set(selected_residues)

        # Stores the original colors and assign the temporary highlight color to the matched
        # substrings.
        original_color_scheme = selected_element.color_by
        for res in selected_element.get_polymer_residues():
            res.store_current_colors() # Stores the original colors.
            if res.seq_index in selected_residues:
                res.custom_color = highlight_color # Sets the highlight color.
            else:
                res.custom_color = res.get_default_color() # Use the original color.

        # Colors only the sequence and highlights the matching residues.
        self.main_window.color_element_by_custom_scheme(selected_element, use_features=False,
                                                        color_structure=False)

        # Restores the original colors. The next time the sequence is changed (for example
        # when inserting/removing gaps), the highlighted residues will be colored back with
        # their original color.
        for res in selected_element.get_polymer_residues():
            res.revert_original_colors()
        selected_element.color_by = original_color_scheme


        # Show a message with the summary of the results.
        def _get_results_string(n_results):
            if n_results == 1:
                return str(n_results) + " time."
            else:
                return str(n_results) + " times."

        results_message = ("Pattern" + _search_string + " found " + _get_results_string(len(finditer_results)))
        self.search_string_window.show_results(results_message, state="found")


        # Select in PyMOL.
        if self.search_string_window.get_select_in_pymol():
            residues_to_select_list = []
            for re_match in finditer_results:
                try:
                    pymol_start_res = selected_element.get_residue_by_index(re_match.start(), only_polymer=True).db_index
                    pymol_end_res = selected_element.get_residue_by_index(re_match.end(), only_polymer=True).db_index
                    residues_to_select_list.extend(list(range(pymol_start_res, pymol_end_res)))
                except IndexError:
                    pass
            if residues_to_select_list:
                selection_name = "pymod_search_string"
                selector_str = "object %s" % selected_element.get_pymol_selector()
                selector_str += " and resi " + self.main_window._join_residues_list(residues_to_select_list)
                cmd.select(selection_name, selector_str)


    #################################################################
    # Clusters.                                                     #
    #################################################################

    def update_cluster_sequences(self, cluster_element):
        """
        Updates the sequences of a cluster when some sequences are removed or added from the
        cluster.
        """
        children = cluster_element.get_children()
        if len(children) > 1:
            cluster_element.adjust_aligned_children_length()
            cluster_element.update_stars()
        else:
            if len(children) == 1:
                children[0].extract_to_upper_level()
            cluster_element.delete()


    def extract_selection_to_new_cluster(self):
        selected_sequences = self.get_selected_sequences()
        original_cluster_index = self.get_pymod_element_index_in_container(selected_sequences[0].mother) + 1
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="generic", child_elements=selected_sequences, algorithm="extracted")
        self.change_pymod_element_list_index(new_cluster, original_cluster_index)
        self.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)


    #################################################################
    # Transfer alignment files.                                     #
    #################################################################

    def transfer_alignment(self, alignment_element):
        """
        Changes the sequences of the elements contained in a PyMod cluster according to the
        information present in an externally supplied file (chosen by users through a file dialog)
        containing the same sequences aligned in a different way. Right now it supports transfer
        only for sequences having the exactly same sequences in PyMod and in the external alignment.
        """
        # Let users choose the external alignment file.
        openfilename, extension = self.choose_alignment_file()
        if None in (openfilename, extension):
            return False

        # Sequences in the aligment currently loaded into PyMod.
        aligned_elements = alignment_element.get_children()[:]

        # Sequences in the alignment files.
        external_records = list(SeqIO.parse(openfilename, extension))

        if len(external_records) < len(aligned_elements):
            title = "Transfer error"
            message = "'%s' has more sequences (%s) than the alignment in '%s' (%s) and the 'Transfer Alignment' function can't be used in this situation." % (alignment_element.my_header, len(aligned_elements), openfilename, len(external_records))
            self.main_window.show_error_message(title,message)
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
            self.main_window.show_error_message(title, message)
            return False

        # Finally transfer the sequences.
        for match in correspondance_list[:]:
            if match["identity"]:
                match["target-seq"].set_sequence(str(match["external-seq"].seq))
                correspondance_list.remove(match)

        self.main_window.gridder(update_clusters=True)


    def delete_cluster_dialog(self, cluster_element):
        title = "Delete Cluster?"
        message = "Are you sure you want to delete %s?" % (cluster_element.my_header)
        remove_cluster_choice = askyesno_qt(message=message, title=title, parent=self.get_qt_parent())
        if not remove_cluster_choice:
            return None
        title = "Delete Sequences?"
        message = "Would you like to delete all the sequences contained in the %s cluster? By selecting 'No', you will only extract them from the cluster." % (cluster_element.my_header)
        remove_children_choice = askyesno_qt(message=message, title=title, parent=self.get_qt_parent())

        # Delete both the cluster and its children.
        if remove_children_choice:
            cluster_element.delete()
        # Delete only the cluster element and extract the child sequences.
        else:
            children = cluster_element.get_children()
            for c in reversed(children[:]):
                c.extract_to_upper_level()
            cluster_element.delete()

        self.main_window.gridder(update_menus=True)


    #################################################################
    # Import PDB files.                                             #
    #################################################################

    def fetch_pdb_files(self, mode, target_selection):
        fp = structural_databases_protocols.Fetch_structure_file(self)
        fp.initialize_from_gui(mode, target_selection)
        fp.launch_from_gui()


    # def associate_structure_from_popup_menu(self, target_element):
    #     """
    #     Launched when users press the 'Associate 3D Structure' from the leeft popup menu.
    #     """
    #     a = structural_databases_protocols.Associate_structure(self, target_element)
    #     a.launch_from_gui()
