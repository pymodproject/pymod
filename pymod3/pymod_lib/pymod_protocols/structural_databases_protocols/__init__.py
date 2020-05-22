# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import gzip
import urllib.request
import time

from pymod_lib.pymod_os_specific import check_network_connection
from pymod_lib.pymod_structure import Parsed_pdb_file
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import PyModInterruptedProtocol


###################################################################################################
# Fetch PDB files.                                                                                #
###################################################################################################

class Fetch_structure_file(PyMod_protocol):
    """
    Class for downloading a PDB file from the sequences retrieved from BLAST.
    """

    protocol_name = "fetch_pdb"

    def __init__(self, pymod, output_directory=None, mode=None, structures_to_fetch=None, import_mode=None):
        PyMod_protocol.__init__(self, pymod, output_directory)
        self.mode = mode
        self.structures_to_fetch = structures_to_fetch
        self.import_mode = import_mode


    def initialize_from_gui(self, mode, structures_to_fetch):
        # Builds a list of structures to be fetched.
        self.structures_to_fetch = []
        if mode == "single":
            self.structures_to_fetch.append(structures_to_fetch)
        elif mode == "selection":
            self.structures_to_fetch.extend(self.pymod.get_selected_sequences())


    def launch_from_gui(self):
        """
        Let the user choose the way in which to retrieve the structures.
        """
        if not check_network_connection("https://www.google.com"):
            message = ("Could not connect to the Internet to download structure files"
                       " from the PDB. Please check your Internet connection.")
            self.pymod.main_window.show_error_message("Connection Error", message)
            return None

        self.fetch_pdb_dialog = Fetch_pdb_dialog(self.pymod.main_window, pymod=self.pymod)
        self.fetch_pdb_dialog.exec_()

        # Gets the import mode from the GUI and then fetch the files.
        self.import_mode = self.fetch_pdb_dialog.fetch_choice

        # Interrupt the process if users close the dialog window.
        if self.import_mode not in ("single-chain", "multiple-chains"):
            return None

        self.fetch_pdb_files()


    def fetch_pdb_files(self):
        """
        This method will call other methods to tetch the PDB files and load the corresponding
        elements in PyMod.
        """
        # List which will store information for each downloaded file.
        self.fetch_results_list = []

        try:
            if not self.pymod.use_protocol_threads:
                self.download_all_pdb_files()
            else:
                label_text = ("Connecting to the PDB to download %s. Please wait for"
                              " completion..." % self.get_seq_text(self.structures_to_fetch, "PDB file"))
                p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                                function=self.download_all_pdb_files,
                                                args=(), title="Downloading PDB files",
                                                wait_start=0.15, wait_end=0.15,
                                                label_text=label_text)
                p_dialog.exec_()

        except PyModInterruptedProtocol:
            return None

        # Warns the user if some structure files could not be fetched.
        n_failures = len([r for r in self.fetch_results_list if r["failure"]])
        if n_failures != 0:
            title = "Connection Error"
            message = ("Can not access the PDB database to download %s structures out of %s."
                       "\nPlease check your Internet connection or if the PDB ids of the"
                       " structures are valid." % (n_failures, len(self.structures_to_fetch)))
            self.pymod.main_window.show_warning_message(title, message)
            if n_failures == 0:
                return None

        # Actually loads the elements in PyMod.
        for fetch_results in self.fetch_results_list:
            if not fetch_results["failure"]:
                try:
                    self._load_structure(fetch_results)
                except Exception as e:
                    title = "Associate Structure Failure"
                    message = ("The structure association for %s chain %s failed because"
                               " of the following error: %s" % (fetch_results["pdb_code"],
                                                                fetch_results["pdb_chain_id"],
                                                                str(e)))
                    self.pymod.main_window.show_error_message(title, message)

        self.pymod.main_window.gridder(clear_selection=True, update_elements=True, update_clusters=True)


    def download_all_pdb_files(self):
        """
        Downloads all the structure files from the PDB.
        """
        for element in self.structures_to_fetch:
            fetch_results = self._fetch_single_element(element)
            self.fetch_results_list.append(fetch_results)


    def _fetch_single_element(self, old_element):
        """
        Download the PDB file for a single element.
        """
        pdb_code = str(old_element.pdb_id)
        pdb_chain_id = str(old_element.pdb_chain)
        fetch_results = {"pdb_code": pdb_code, "pdb_chain_id": pdb_chain_id,
                         "element": old_element}
        try:
            pdb_file_shortcut = fetch_structure_file(pdb_code=pdb_code,
                output_dir=self.pymod.temp_directory_dirpath)
            fetch_results["failure"] = False
            fetch_results["pdb_file_shortcut"] = pdb_file_shortcut
        except IOError:
            fetch_results["failure"] = True
        return fetch_results


    def _load_structure(self, fetch_results):
        """
        Parses a downloaded PDB file, builds PyMod elements for its chains and loads the corresponding
        elements in PyMod/PyMOL.
        """

        old_element = fetch_results["element"]
        pdb_code = fetch_results["pdb_code"]
        pdb_file_shortcut = fetch_results["pdb_file_shortcut"]
        pdb_chain_id = fetch_results["pdb_chain_id"]

        #----------------------------------------------------------------------------
        # Load in PyMod only the chain corresponding to the hit sequence and adjust -
        # its length to the region identified by BLAST.                             -
        #----------------------------------------------------------------------------
        if self.import_mode == "single-chain":
            a = Associate_structure(self.pymod, old_element)
            a.associate(pdb_file_shortcut, pdb_chain_id)

        #----------------------------------------------------------------------
        # Load each chain found in the PDB file where the 3D structure of the -
        # hit sequence is present.                                            -
        #----------------------------------------------------------------------
        elif self.import_mode == "multiple-chains":
            # Deletes the original hit sequence retrieved by BLAST and replaces it with
            # a new element with an associated structure loaded in PyMOL.
            old_element.delete()
            # Builds 'Pymod_elements' objects for each chain present in the PDB file.
            new_elements = self.pymod.open_structure_file(os.path.abspath(pdb_file_shortcut))
            # Color other chains by gray.
            other_chains_elements = [e for e in new_elements if e.get_chain_id() != pdb_chain_id]
            if other_chains_elements:
                self.pymod.main_window.color_selection("multiple", other_chains_elements, "regular", regular_color="gray")


class Structure_file_fetcher:
    """
    Class to fetch a PDB file from the internet.
    """

    def __init__(self, pdb_code, output_dir, new_name=None):
        self.pdb_code = pdb_code
        self.output_dir = output_dir
        # Form the pdb output name
        if new_name:
            self.output_name = '%s.pdb' % new_name
        else:
            self.output_name = '%s.pdb' % pdb_code

    def fetch(self):
        """
        Actually retrieves the PDB file from the internet.
        """
        pdb_url = "http://www.rcsb.org/pdb/files/%s.pdb.gz" % self.pdb_code
        temp_gzip_file_name = urllib.request.urlretrieve(pdb_url)[0]
        open_gzip_file = gzip.open(temp_gzip_file_name) # Uncompress the file while reading
        output_path = os.path.join(self.output_dir, self.output_name)
        saved_file = open(output_path, 'wb')
        saved_file.write(open_gzip_file.read()) # Write pdb file
        open_gzip_file.close()
        saved_file.close()
        return output_path

def fetch_structure_file(pdb_code, output_dir, new_name=None):
    sf = Structure_file_fetcher(pdb_code, output_dir, new_name=new_name)
    return sf.fetch()


###################################################################################################
# Associate structures to PyMod elements.                                                         #
###################################################################################################

from pymol import cmd

from pymod_lib.pymod_seq.seq_manipulation import global_pairwise_alignment


class Associate_structure(PyMod_protocol):
    """
    Once 'build_structure_objects()' has been used, this will edit the 'PyMod_element' and
    'Structure' objects corresponding to the 'chain_id' according to the original element sequence.
    Usually this is used when fetching a PDB file corresponding to some hit from a BLAST search,
    because hits in HSPs may have a shorter sequence with respect to the full PDB chain.
    This method can be called to crop the full 3D chain according to the hit sequence in the
    HSP.
    """

    temp_full_name = "pymod_full_temp"
    temp_fragment_name = "pymod_fragment_temp"
    protocol_name = "associate_structure"

    def __init__(self, pymod, pymod_element):
        PyMod_protocol.__init__(self, pymod, output_directory=pymod.temp_directory_dirpath)
        # Directory in which to save the temporary files.
        self.temp_directory = self.pymod.temp_directory_dirpath
        self.target_element = pymod_element
        self.associate_pdb_file = None

    def associate(self, pdb_file_path, chain_id):
        """
        Actually associates the structure.
        """
        self._set_options(pdb_file_path, chain_id)

        # Parses the source structure file.
        if not self.associate_pdb_file:
            self.associate_pdb_file = Parsed_pdb_file(self.pymod, self.original_pdb_file_path,
                                                      copy_original_file=False, save_chains_files=False)

        #-----------------------------------------------------------------------
        # Check if the pymod element can be associated to the structure chain. -
        #-----------------------------------------------------------------------
        if not self.chain_id in self.associate_pdb_file.get_chains_ids():
            raise KeyError("The structure file '%s' does not have chain '%s'." % (self.original_pdb_file_path, self.chain_id))

        structure_chain_element = self.associate_pdb_file.get_pymod_element_by_chain(self.chain_id)

        # Check if the the target sequence and the sequence of the structure to associate match by
        # aligning the two sequences using dynamic programming.
        ali = global_pairwise_alignment(self.target_element.my_sequence.replace("-", ""),
                                        structure_chain_element.my_sequence.replace("-", ""),
                                        toss_modres=True)

        # If the sequences do not match, interrupt the process.
        # if ali["id"] < 99.9: # TODO: use a 'PyMod_alignment' class.
        #     raise ValueError("The target sequence does not match with the sequence of the structure to associate (sequence identity percentage = %s)." % ali["id"])

        #-------------------------------------------------------------------------------------
        # Gets information about matching and missing residues in the two aligned sequences. -
        #-------------------------------------------------------------------------------------
        pc = 0 # Alignment position counter.
        hc = 0 # Target residue counter.
        tc = 0 # PDB structure residue counter.
        matching_positions = [] # list of matching positions.
        missing_positions = [] # list of missing residues in the pdb structure with respect to the target sequence.
        for hr, tr in zip(ali["seq1"], ali["seq2"]):
            if hr != "-" and tr != "-" and hr == tr:
                matching_positions.append({"pc": pc, "hc": hc, "tc": tc})
            if tr == "-" and hr != "-":
                missing_positions.append({"pc": pc, "hc": hc, "tc": tc})
            if hr != "-":
                hc += 1
            if tr != "-":
                tc += 1
            pc += 1
        # Gets the starting and ending positions (using the PDB numeration) that will be used to
        # crop the 3D structure.
        start_position = structure_chain_element.get_residue_by_index(matching_positions[0]["tc"]).db_index
        end_position = structure_chain_element.get_residue_by_index(matching_positions[-1]["tc"]).db_index

        #------------------------------------------------
        # Use PyMOL to build the new cropped structure. -
        #------------------------------------------------
        # First loads the full PDB structure of the chain in PyMOL.
        cmd.load(self.original_pdb_file_path, self.temp_full_name)
        # Select amino acid residues ranging from the starting and ending positions which define
        # the fragment to "excise".
        cmd.select(self.temp_fragment_name, "resi %s-%s and chain %s and object %s and not hetatm" % (start_position, end_position, self.chain_id, self.temp_full_name))
        # Join the selections and save a file PDB file of the cropped fragment.
        pdb_filename = os.path.splitext(os.path.basename(self.original_pdb_file_path))[0] # structure_chain_element.get_structure_file_root()
        pdb_basename = "%s_cropped.pdb" % pdb_filename
        cropped_structure_file_shortcut = os.path.join(self.temp_directory, pdb_basename)
        cmd.save(cropped_structure_file_shortcut, self.temp_fragment_name)
        # Clean up the selections.
        cmd.delete(self.temp_full_name)
        cmd.delete(self.temp_fragment_name)

        #----------------------------------------------------------------------------------
        # Builds a 'Parsed_pdb_file' object for the PDB file of the structure just saved. -
        #----------------------------------------------------------------------------------
        p = Parsed_pdb_file(self.pymod, os.path.abspath(cropped_structure_file_shortcut),
                            output_directory=self.pymod.structures_dirpath)
        new_element_with_structure = p.get_pymod_element_by_chain(self.chain_id)
        adjust_sequence = self.target_element.my_sequence
        self.pymod.replace_element(self.target_element, new_element_with_structure)

        #--------------------------------------------------------------------------------------
        # Updates the sequence of the fragment to keep it in frame with the original sequence -
        # by including the target sequence indels.                                            -
        #--------------------------------------------------------------------------------------
        list_of_missing_positions = [p["hc"] for p in missing_positions]
        new_sequence = []
        adc = 0
        for i, p in enumerate(adjust_sequence):
            if adc in list_of_missing_positions:
                new_sequence.append("-")
            else:
                # Right now modified residues are not included in the cropped structures,
                # this prevents them from being included in the chain sequence.
                if p != "X":
                    new_sequence.append(p)
                else:
                    new_sequence.append("-")
            if p != "-":
                adc += 1
        new_sequence = "".join(new_sequence)

        # If the old and new sequences do not match, then align them. This situation does not occur
        # often, and this procedure might introduce alignment errors (but it will not make the
        # plugin crash).
        if new_element_with_structure.my_sequence.replace("-", "") != new_sequence.replace("-", ""):
            ali_new = global_pairwise_alignment(adjust_sequence,
                                                new_element_with_structure.my_sequence,
                                                toss_modres=True)
            new_sequence = ali_new["seq2"]

        new_element_with_structure.set_sequence(new_sequence, permissive=False)


    def _set_options(self, pdb_file_path, chain_id):
        self.original_pdb_file_path = pdb_file_path
        self.chain_id = chain_id


    #################################################################
    # Launch from the GUI.                                          #
    #################################################################

    def launch_from_gui(self):
        # Builds a new window.
        self.associate_structure_window = Associate_structure_window(
            parent = self.pymod.get_pymod_app(),
            protocol = self,
            title = "Associate Structure",
            upper_frame_title = "Associate 3D Structure Options",
            submit_command = self.associate_structure_state)
        # This will be set to 'True' once the users select a valid PDB file and press the 'SUBMIT'
        # button.
        self._select_associate_chain = False


    def associate_structure_state(self):
        """
        Called when users press the 'SUBMIT' window of the 'Associate Structure' protocol. This is
        actually both when selecting the structure file path and when selecting a chain of the file.
        """
        #-----------------------------------------------------------------
        # Checks if a correct structure file has been provided as input. -
        #-----------------------------------------------------------------
        if not self._select_associate_chain:

            if not self.associate_structure_window.check_general_input():
                return False

            self.pdb_file_path_from_gui = self.associate_structure_window.get_structure_file()

            if not os.path.isfile(self.pdb_file_path_from_gui):
                self.pymod.main_window.show_error_message("File Error", "Please select a valid file path.")
                return False

            if not self.pymod.is_valid_structure_file(self.pdb_file_path_from_gui, show_error=False):
                self.pymod.main_window.show_error_message("File Type Error", "Please select a valid PDB file.")
                return False

            # Parses the structure file.
            self.associate_pdb_file = Parsed_pdb_file(self.pymod, self.pdb_file_path_from_gui,
                                                      copy_original_file=False, save_chains_files=False)
            # Gets its chains.
            self.available_chains = self.associate_pdb_file.get_chains_ids()
            self.associate_structure_window.show_chain_selection_frame()

            self._select_associate_chain = True

        #----------------------------------------------------------------------------------------
        # If a valid structure file has been provided, this will try to associate the structure -
        # of the chain specified in the combobox to the target element.                         -
        #----------------------------------------------------------------------------------------
        elif self._select_associate_chain:
            try:
                self.associate(self.pdb_file_path_from_gui, self.associate_structure_window.get_structure_chain())
                self.associate_structure_window.destroy()
                self.pymod.main_window.gridder(update_elements=True, update_clusters=True)
            # except Exception, e:
            except Exception:
                title = "Associate Structure Failure"
                # message = "The structure association failed because of the following error: %s" % e
                message = "The structure association failed because of an error"
                self.pymod.main_window.show_error_message(title, message)


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

from pymol.Qt import QtWidgets

# from pymod_lib.pymod_gui.shared_gui_components_qt import *


class Fetch_pdb_dialog(QtWidgets.QDialog):
    """
    Dialog to select the way in which structure files downloaded from the PDB have to be fetched.
    """

    is_pymod_window = True

    def __init__(self, parent, pymod, *args, **kwargs):
        super(Fetch_pdb_dialog, self).__init__(parent, *args, **kwargs)
        self.initUI()
        self.fetch_choice = None
        self.pymod = pymod


    def initUI(self):

        self.setWindowTitle('Import PDB Options')

        vertical_layout = QtWidgets.QVBoxLayout()

        # Installation options label.
        info_text = "Please select the 3D structure import mode:"
        self.fetch_info_label = QtWidgets.QLabel(info_text, self)
        # self.fetch_info_label.setStyleSheet(label_style_1)
        vertical_layout.addWidget(self.fetch_info_label)

        vertical_layout.addStretch(1)

        # Import options radiobuttons.
        horizontal_layout = QtWidgets.QHBoxLayout()

        self.import_all_radiobutton = QtWidgets.QRadioButton("Import in PyMod the structure of every chain of the PDB files.")
        # self.import_all_radiobutton.setChecked(True)
        # self.import_all_radiobutton.setStyleSheet(label_font_1)
        vertical_layout.addWidget(self.import_all_radiobutton)
        self.import_fragment_radiobutton = QtWidgets.QRadioButton("Import in PyMod only the structure of the hit sequences fragments.")
        # label_font_1
        self.import_fragment_radiobutton.setStyleSheet("margin-bottom: 10px")
        vertical_layout.addWidget(self.import_fragment_radiobutton)


        # Import fragments button.
        self.import_button = QtWidgets.QPushButton("Import 3D Structures", self)
        # self.import_button.setStyleSheet(label_style_2)
        self.import_button.clicked.connect(self.on_import_button_click)
        horizontal_layout.addWidget(self.import_button)

        horizontal_layout.addStretch(1)

        # Cancel button.
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        # self.cancel_button.setStyleSheet(label_style_2)
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)
        self.setLayout(vertical_layout)


    def on_cancel_button_click(self):
        self.fetch_choice = None
        self.close()


    def on_import_button_click(self):
        if self.import_all_radiobutton.isChecked():
            self.fetch_choice = "multiple-chains"
            time.sleep(0.1)
            self.close()
        elif self.import_fragment_radiobutton.isChecked():
            self.fetch_choice = "single-chain"
            time.sleep(0.1)
            self.close()
        else:
            self.fetch_choice = None
            message = ("Please select one out of two available import modes in order"
                       " to download the PDB files.")
            self.pymod.main_window.show_warning_message("Warning", message)


# import pymod_lib.pymod_vars as pmdt
# from pymod_lib.pymod_gui import shared_gui_components
#
# class Associate_structure_window(shared_gui_components.PyMod_tool_window):
#
#     def __init__(self, parent = None, protocol = None, **configs):
#         shared_gui_components.PyMod_tool_window.__init__(self, parent=parent , **configs)
#         self.current_protocol = protocol
#         # An entryfield to select the structure file.
#         self.structure_file_enf = shared_gui_components.PyMod_path_entryfield(self.midframe,
#             label_text = "Select Structure File",
#             label_style = shared_gui_components.label_style_1,
#             path_type = "file",
#             file_types = pmdt.all_structure_file_types_atl,
#             askpath_title = "Select Structure File")
#         self.structure_file_enf.pack(**shared_gui_components.pack_options_1)
#         self.add_widget_to_align(self.structure_file_enf)
#         self.add_widget_to_validate(self.structure_file_enf)
#         self.align_widgets(15)
#
#     def get_structure_file(self):
#         return self.structure_file_enf.getvalue()
#
#
#     def show_chain_selection_frame(self):
#         # Removes the entryfield to select the structure file.
#         self.structure_file_enf.pack_forget()
#
#         # Displays a combobox to select the chain id of corresponind to the structure to be
#         # associated with the target sequence.
#         self.chain_selection_cbx = shared_gui_components.PyMod_combobox(self.midframe,
#             label_text = 'Select Chain to Associate',
#             label_style = shared_gui_components.label_style_1,
#             scrolledlist_items=self.current_protocol.available_chains)
#         self.chain_selection_cbx.pack(**shared_gui_components.pack_options_1)
#         self.chain_selection_cbx.selectitem(0)
#         self.add_widget_to_align(self.chain_selection_cbx)
#         self.align_widgets(15)
#
#     def get_structure_chain(self):
#         return self.chain_selection_cbx.get()
