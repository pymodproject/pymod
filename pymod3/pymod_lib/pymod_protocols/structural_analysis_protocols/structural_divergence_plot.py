# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing a protocol for analyzing the structural divergence between aligned protein
structures in PyMod. Builds a plot with CA-CA distance for each residue or a RMSF plot.
"""

import os

import numpy as np
from pymol import cmd

from pymod_lib import pymod_vars
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_combobox_qt)
from pymod_lib.pymod_protocols.structural_analysis_protocols._base_structural_analysis import get_distance, Structural_analysis_mixin
from pymod_lib.pymod_plot.pymod_plot_qt import PyMod_plot_window_qt
from pymod_lib.pymod_threading import Protocol_exec_dialog


###################################################################################################
# Main class.                                                                                     #
###################################################################################################

class Structural_divergence_plot(PyMod_protocol, Structural_analysis_mixin):
    """
    Class for implementing structural divergence analyses (sda) in PyMod.
    """

    protocol_name = "structural_divergence_plot"

    def __init__(self, pymod):
        PyMod_protocol.__init__(self, pymod)
        self.target_sequences = self.get_pymod_elements()


    def launch_from_gui(self):
        """
        Checks the condition for launching a structural divergence plot analysis.
        """

        error_title = "Selection Error"

        # Check that no nucleic acid structure has been selected.
        if not self.check_protein_selection(error_title):
            return None

        # Check that at least two elements have been selected.
        if len(self.target_sequences) < 2:
            self.pymod.main_window.show_error_message(error_title,
                "Please select at least two structures in the same alignment.")
            return None

        # Check for multiple aligned structures in the same cluster.
        if not self.check_multiple_selection(error_title):
            return None


        if not self.target_sequences[0].mother.algorithm in pymod_vars.structural_alignment_tools:
            message = pymod_vars.structural_alignment_warning % "Structural Divergence Plot"
            self.pymod.main_window.show_warning_message("Alignment Warning", message)

        # Opens the options window.
        self.sda_options_window = SDA_options_window_qt(self.pymod.main_window,
            protocol=self,
            title="Structural Divergence Plot Options",
            upper_frame_title="Options for Structural Divergence Plot",
            submit_command=self.sda_state)
        self.sda_options_window.show()


    def sda_state(self):

        #--------------------------------------
        # Checks the parameters from the GUI. -
        #--------------------------------------

        # Analysis mode.
        analysis_mode, analysis_mode_full_name = self.sda_options_window.get_analysis_mode()

        self.reference_element_id = self.sda_options_window.get_reference_id()
        self.reference_element = self.target_sequences[self.reference_element_id]

        self.sda_options_window.destroy()

        #----------------------------------
        # Actually launches the analysis. -
        #----------------------------------

        if analysis_mode == "rmsf":
            self.launch_rmsf_analysis()

        elif analysis_mode == "distances":
            self.launch_distances_divergence_analysis()

        else:
            KeyError(analysis_mode)


    ##########################################
    # Root mean square fluctuation analysis. #
    ##########################################

    def launch_rmsf_analysis(self):
        """
        Launches the RMSF analysis. A plot containing the RMSF profile of the selected structures
        will be shown.
        """

        # Prepares the atomic coordinates of the structures.
        ali_seqs_dicts = []
        res_dicts = []
        for element in self.target_sequences:
            ali_seq_dict, res_dict = self._get_sda_data(element)
            ali_seqs_dicts.append(ali_seq_dict)
            res_dicts.append(res_dict)


        # Actually computes the RMSF of the residues.
        ali_len = len(self.target_sequences[0].my_sequence)
        rmsf_list = []
        residues_additional_data = []


        for ali_pos in range(0, ali_len):

            # For each alignment position, get the coordinates of the residues.
            coords_list = []
            reference_res = None
            for structure_idx, (ali_seq_dict, res_dict) in enumerate(zip(ali_seqs_dicts, res_dicts)):

                # Gets a residue.
                res = self._get_res(ali_pos, ali_seq_dict)
                # Gets its coordinates.
                coords_i = self._get_coords(res, res_dict)
                if coords_i is None:
                    continue
                coords_list.append(coords_i)
                # Stores the reference residue.
                if structure_idx == self.reference_element_id:
                    reference_res = res

            # Get the reference structure information.
            if reference_res is not None:
                residues_additional_data.append({"residue_name": reference_res.three_letter_code,
                                                 "residue_pdb_position": reference_res.db_index,
                                                 "pymol_selector": reference_res.get_pymol_selector(),
                                                 "export_label": "%s %s" % (reference_res.three_letter_code, reference_res.db_index)})
            else:
                residues_additional_data.append({"residue_name": "none",
                                                 "residue_pdb_position": "-",
                                                 "pymol_selector": None,
                                                 "export_label": "None"})

            # For each alignment position, get the coordinates of the residues.
            if len(coords_list) > 1:
                coords_list = np.array(coords_list)
                # Take the centroid.
                centroid = coords_list.mean(axis=0)
                # Computes the RMSF.
                rmsf = np.sqrt(np.mean([np.linalg.norm(ci-centroid)**2 for ci in coords_list]))
                rmsf_list.append(round(rmsf, 3))

            elif len(coords_list) == 1:
                rmsf_list.append(None)

            else:
                rmsf_list.append(None)


        # Builds the structural divergence plot.
        cp = PyMod_plot_window_qt(self.pymod.main_window)
        cp.initialize(pymod=self.pymod, title="RMSF Plot")
        messagebar_text_on_update_qt = ("Selected: __residue_name__ __residue_pdb_position__ " +
                                        self.reference_element.compact_header +
                                        " (alignment position: __x__), RMSF value: __y__ " + "\u212B")
        cp.build_plotting_area(messagebar_initial_text="Set the 'Interact on click' option to 'Yes' to click on the plot and to highlight residues in PyMOL.",
                               update_messagebar=True,
                               messagebar_text_on_update=messagebar_text_on_update_qt,
                               on_click_action=self._highlight_in_pymol_from_sda_plot_qt,
                               x_label_text="Alignment position",
                               y_label_text="RMSF %s" % "\u212B",)

        # Values on the x-axis.
        x_series = list(range(1, ali_len+1))
        # Adds the RMSF plot.
        cp.add_plot(x_series, rmsf_list, label="RMSF plot", additional_data=residues_additional_data)
        # Actually shows the plot.
        cp.show()


    #################################
    # Distance divergence analysis. #
    #################################

    def launch_distances_divergence_analysis(self):
        """
        Launches the "Distance Divergence" analysis. Once the user has selected a reference structure,
        the distances between its CA atoms and the aligned CA atoms of the other structures will be
        calculated. A plot with the CA-CA distance profiles (one profile for each structure aligned
        to the reference one) will be shown.
        """

        # Prepares the atomic coordinates of the structures. First, get the reference element data.
        ref_ali_seq_dict, ref_res_dict = self._get_sda_data(self.reference_element)

        # Get the data for the rest of the elements.
        ali_seqs_dicts = []
        res_dicts = []

        structures_to_analyze = [e for e in self.target_sequences if e is not self.reference_element]
        for element in structures_to_analyze:
            ali_seq_dict, res_dict = self._get_sda_data(element)
            ali_seqs_dicts.append(ali_seq_dict)
            res_dicts.append(res_dict)


        # Initializes the plot.
        cp = PyMod_plot_window_qt(self.pymod.main_window)
        cp.initialize(pymod=self.pymod,
                      title="CA-CA Distance Plot with reference: %s" % self.reference_element.compact_header)
        messagebar_text_on_update_qt = ("Selected: __residue_name__ __residue_pdb_position__"
                                        " of __plot_name__ (alignment position: __x__), CA-CA"
                                        " distance value: __y__ " + "\u212B")
        cp.build_plotting_area(messagebar_initial_text="Set the 'Interact on click' option to 'Yes' to click on the plot and to highlight residues in PyMOL.",
                               update_messagebar=True,
                               messagebar_text_on_update=messagebar_text_on_update_qt,
                               on_click_action=self._highlight_in_pymol_from_sda_plot_qt,
                               x_label_text="Alignment position",
                               y_label_text="CA distance %s" % "\u212B",)

        # Values on the x-axis.
        ali_len = len(self.target_sequences[0].my_sequence)
        y_series_list = [[] for i in range(0, len(structures_to_analyze))]
        additional_data_list = [[] for i in range(0, len(structures_to_analyze))]

        # For each alignment position, get the coordinates of the reference residue (if present in
        # that position).
        for ali_pos in range(0, ali_len):

            # Gets a reference residue and its coordinates.
            ref_res = self._get_res(ali_pos, ref_ali_seq_dict)
            ref_coords_i = self._get_coords(ref_res, ref_res_dict)

            # Cycle through all the structures aligned to the reference one.
            for ali_seq_dict, res_dict, y_series, additional_data in zip(ali_seqs_dicts,
                                                                         res_dicts,
                                                                         y_series_list,
                                                                         additional_data_list):
                # Gets a residue and its coordinates.
                res = self._get_res(ali_pos, ali_seq_dict)
                coords_i = self._get_coords(res, res_dict)
                # The reference structure does not have coordinates in this alignment position.
                if coords_i is None or ref_coords_i is None:
                    y_series.append(None)
                    additional_data.append({"residue_name": "none",
                                            "residue_pdb_position": "-",
                                            "pymol_selector": None,
                                            "export_label": "None"})
                    continue

                # Gets the distance between the reference residue and the target one.
                dist = get_distance(ref_coords_i, coords_i)
                y_series.append(round(dist, 3))
                additional_data.append({"residue_name": res.three_letter_code,
                                        "residue_pdb_position": res.db_index,
                                        "pymol_selector": res.get_pymol_selector(),
                                        "export_label": "%s %s" % (res.three_letter_code, res.db_index)})

        # Adds the CA-CA distances plot.
        x_series = list(range(1, ali_len+1))
        for element, y_series, additional_data in zip(structures_to_analyze,
                                                      y_series_list,
                                                      additional_data_list):
            cp.add_plot(x_series, y_series, label=element.compact_header,
                        additional_data=additional_data)

        # Actually shows the plot.
        cp.show()


    def _highlight_in_pymol_from_sda_plot_qt(self, point_data):
        sel = point_data["pymol_selector"]
        if sel is not None:
            cmd.select("pymod_selection", sel)
            cmd.center("pymod_selection")


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class SDA_options_window_qt(PyMod_protocol_window_qt):
    """
    Window used in the contact map analysis.
    """

    analysis_modes_list = ["CA-CA Distances", "CA RMSF"]

    def build_protocol_middle_frame(self):

        # Analysis mode.
        self.analysis_mode_rds = PyMod_radioselect_qt(label_text="Analysis Mode",
                                                      buttons=self.analysis_modes_list)
        self.analysis_mode_rds.setvalue(self.analysis_modes_list[0])
        self.middle_formlayout.add_widget_to_align(self.analysis_mode_rds)

        # Reference structure combobox.
        structures_list = [element.my_header for element in self.protocol.target_sequences]
        self.reference_combobox = PyMod_combobox_qt(label_text="Reference Structure",
                                                    items=structures_list)
        self.reference_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.reference_combobox)

        self.middle_formlayout.set_input_widgets_width("auto", padding=10)


    def get_analysis_mode(self):
        feature_type = self.analysis_mode_rds.getvalue()
        if feature_type == self.analysis_modes_list[0]:
            return "distances", feature_type
        elif feature_type == self.analysis_modes_list[1]:
            return "rmsf", feature_type
        else:
            raise KeyError(feature_type)


    def get_reference_id(self):
        try:
            return self.reference_combobox.get_index()
        except:
            print("- WARNING: could not obtain the reference structure id.")
            return 0
