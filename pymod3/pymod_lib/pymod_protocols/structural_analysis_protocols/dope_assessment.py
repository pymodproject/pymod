# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module used to compute the DOPE scores (using the MODELLER package) of elements loaded in PyMod.
"""

import os
import sys

import numpy as np

from pymol import cmd

try:
    import modeller
    from modeller.scripts import complete_pdb
except:
    pass

import pymod_lib.pymod_seq.seq_manipulation as pmsm

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol, MODELLER_common
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import catch_protocol_exception

from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt
from pymod_lib.pymod_plot.pymod_plot_qt import PyMod_plot_window_qt


###################################################################################################
# DOPE PROFILES.                                                                                  #
###################################################################################################

class DOPE_assessment(PyMod_protocol, MODELLER_common):
    """
    Compute the DOPE (Discrete optimized protein energy) of a polypeptidic chain using MODELLER.
    """

    script_file_basename = "dope_profile_script"
    script_file_name = "%s.py" % script_file_basename
    script_temp_output_name = "dope_profile_temp_out.txt"
    protocol_name = "dope_assessment"

    def additional_initialization(self):
        MODELLER_common.__init__(self)
        self.selected_sequences = []
        self.unfiltered_dope_scores_dict = {} # Items will contain DOPE scores for ligands and water molecules.
        self.dope_scores_dict = {} # Items will contain DOPE scores of only polymer residues.
        self.assessed_structures_list = []
        self.global_dope_scores = []
        self.remove_temp_files = True


    @catch_protocol_exception
    def launch_from_gui(self):
        """
        Called when users decide calculate DOPE of a structure loaded in PyMod.
        """

        self.selected_sequences = self.pymod.get_selected_sequences() # self.get_pymod_elements(self.selected_sequences)

        #-----------------------------------------------
        # Checks if the DOPE profiles can be computed. -
        #-----------------------------------------------
        modeller_error = self.pymod.modeller.check_exception()
        if modeller_error is not None:
            message = "In order to compute DOPE scores of a structure, MODELLER must be installed and configured correctly. %s" % modeller_error
            self.pymod.main_window.show_error_message("MODELLER Error", message)
            return None

        if len(self.selected_sequences) == 0:
            self.pymod.main_window.show_error_message("Selection Error", "Please select at least one structure to assess.")
            return None
        if any([e.polymer_type == "nucleic_acid" for e in self.selected_sequences]):
            self.pymod.main_window.show_error_message("Selection Error", "Can not perform DOPE assessment on nucleci acids structures.")
            return None
        if not self.pymod.all_sequences_have_structure(self.selected_sequences):
            self.pymod.main_window.show_error_message("Selection Error", "Please select only elements that have a 3D structure currently loaded in PyMOL.")
            return None
        if len(self.selected_sequences) > 1:
            mothers_set = set([seq.mother for seq in self.selected_sequences])
            if self.pymod.root_element in mothers_set or len(mothers_set) > 1:
                self.pymod.main_window.show_error_message("Selection Error", "You can assess multiple structures DOPE only if they are aligned in the same cluster.")
                return None

        # Ask users if they would like to color the sequences according to their DOPE values.
        title = "Color Option"
        message = "Would you like to color the selected sequences by their DOPE values, once they have been calculated?"
        color_by_dope_choice = askyesno_qt(message=message, title=title, parent=self.pymod.get_qt_parent())

        #----------------------------------------
        # Initializes the MODELLER environment. -
        #----------------------------------------

        if self.run_modeller_internally:
            with open(os.devnull, "w") as n_fh: # Silence some MODELLER output.
                original_stdout = sys.stdout
                sys.stdout = n_fh
                env = self._initialize_env()
                sys.stdout = original_stdout
        else:
            env = None


        #------------------------------------------------------------------------------------
        # Actually computes the DOPE scores of the polypeptide chains in the user selection -
        # and assigns to each residue of a corresponding color according to its DOPE.       -
        #------------------------------------------------------------------------------------

        if not self.pymod.use_protocol_threads:
            self.compute_all_dopes(env=env)
            self.assign_dope_items()

        else:

            label_text = ("Computing DOPE of %s. Please wait for the process to"
                          " complete..." % self.get_seq_text(self.selected_sequences, "structure"))

            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self._launch_dope_thread,
                                            args=(None, ), wait_start=1.25, wait_end=0.25,
                                            title="Computing DOPE scores.",
                                            label_text=label_text,
                                            lock=True, # a thread with MODELLER can not exit safely.
                                            lock_title=self.pymod.modeller_lock_title,
                                            stdout_silence=True,
                                            lock_message=self.pymod.modeller_lock_message)
            p_dialog.exec_()

        #----------------------
        # Color the elements. -
        #----------------------

        if color_by_dope_choice:
            for element in self.selected_sequences:
                self.pymod.main_window.color_element_by_dope(element)

        #--------------------------------
        # Shows the DOPE profiles plot. -
        #--------------------------------

        if len(self.selected_sequences) == 1:
            dope_graph_mode = "single"
        elif len(self.selected_sequences) >= 2:
            dope_graph_mode = "multiple"

        # Prepares the data to show in the plot.
        self.dope_plot_data = self.prepare_dope_plot_data(self.selected_sequences, mode=dope_graph_mode)
        # Shows the plot.
        self.show_plot()

        #-----------------------------
        # Shows an assessment table. -
        #-----------------------------

        column_headers = ["Structure Name", "DOPE score"]

        data_array = list(zip([e.get_pymol_selector() for e in self.assessed_structures_list], self.global_dope_scores))
        assessment_table_args = {"column_headers": column_headers,
                                 # "row_headers": [m["name"] for m in a.outputs],
                                 "data_array": data_array,
                                 "title": "Assessment of Selected Structures",
                                 "width": 850, "height": 420,
                                 "sortable": True,
                                 }
        self.pymod.show_table(**assessment_table_args)



    def add_element(self, pymod_element):
        self.selected_sequences.append(pymod_element)


    def _launch_dope_thread(self, env):
        self.compute_all_dopes(env=env)
        self.assign_dope_items()


    def compute_all_dopes(self, env=None):
        for element in self.selected_sequences:
            self._compute_dope_of_element(element, env=env)


    def _compute_dope_of_element(self, element, env=None):
        # Prepares the input for MODELLER.
        e_file_name = element.get_structure_file(strip_extension=False)
        e_file_shortcut = os.path.join(self.pymod.structures_dirpath, e_file_name)
        e_profile_file_shortcut = os.path.join(self.output_directory, e_file_name+".profile")
        # Computes the DOPE of the 3D structure of the chain of the 'element'.
        global_dope_score = self._compute_dope_of_structure_file(e_file_shortcut, e_profile_file_shortcut, env=env)
        # Reads the output file produced by MODELLER with the DOPE scores of the chain of the
        # 'element'.
        dope_scores = self._get_dope_profile(e_profile_file_shortcut)
        # Removes the temporary profile file.
        if self.remove_temp_files:
            os.remove(e_profile_file_shortcut)
        # Stores the results.
        self.unfiltered_dope_scores_dict.update({element: dope_scores})
        self.assessed_structures_list.append(element)
        self.global_dope_scores.append(global_dope_score)


    def _compute_dope_of_structure_file(self, str_file_path, profile_file_path, env=None):
        return self._compute_dope(str_file_path, profile_file_path, env=env,)


    def _compute_dope(self, str_file_path, profile_file_path, env=None):
        """
        Uses MODELLER to compute the DOPE of a polypeptidic chain, and ouptuts the results in
        'profile_file_path'. When 'env' is set to 'None', MODELLER will be initialized. If
        MODELLER has already been initialized, the its 'env' varibale can be passed in this
        argument so that it is not initialized again.
        """
        if env == None:
            env = self._initialize_env()
        modstr = complete_pdb(env, str(str_file_path))
        # Assess with DOPE.
        s = modeller.selection(modstr).only_std_residues() # only_het_residues, only_std_residues, only_water_residues
        # Gets the DOPE score.
        score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',file=str(profile_file_path), normalize_profile=True, smoothing_window=15)

        return score


    def _get_dope_profile(self, profile_file_name, seq=None):
        """
        Read 'profile_file' into a Python list, and add gaps corresponding to the alignment
        sequence 'seq'.
        """
        profile_file = open(profile_file_name,"r")
        vals = []
        for line in profile_file.readlines():
            res_three_letter_code = line[8:11]
            # Read all non-comment and non-blank lines from the file:
            if not line.startswith('#') and len(line) > 10:
                # Initially do not exclude also water molecules (named 'TIP3') and heteroresidues
                # from the graph.
                spl = line.split()
                vals.append(float(spl[-1]))
        profile_file.close()
        return vals


    def assign_dope_items(self):
        # Retain only the DOPE values for residues of the chain (standard and modified residues).
        for chain_element in self.assessed_structures_list:
            all_chain_dope_scores = self.unfiltered_dope_scores_dict[chain_element]
            filtered_chain_dope_scores = []
            for res, score in zip(chain_element.residues, all_chain_dope_scores):
                if res.is_polymer_residue():
                    filtered_chain_dope_scores.append(score)
            self.dope_scores_dict.update({chain_element: filtered_chain_dope_scores})
        # Builds a list of all DOPE values of the residues in the selection.
        ldope = []
        for chain_element in self.assessed_structures_list:
            ldope.extend(self.dope_scores_dict[chain_element])
        # Takes the min and max values among all the selected residues.
        min_value = min(ldope)
        max_value = max(ldope)
        # An array with the equally sapced limits generated with the list above.
        bins = np.array(np.linspace(min_value, max_value, num=10))
        for chain_element in self.assessed_structures_list:
            # An array with all the DOPE values of a single chain in the selection.
            adope = np.array(self.dope_scores_dict[chain_element])
            # An array with the id of the bins where those values reside.
            inds = np.digitize(adope, bins)
            # Returns a list like:
            # [(-0.052, 4), (-0.03, 3), (-0.04, 5), (-0.04, 6), (-0.041, 7), (-0.042, 8), (-0.043, 10), ...]
            # which contains for all standard residues of a polypeptidic chain a tuple. The
            # first value of the tuple is the DOPE score of that residues, the second is the id
            # (going from 1 to 10) of the bin where that value resides.
            dope_items = []
            for dope_score, bin_id in zip(adope, inds): # zip(ldope, inds):
                dope_items.append({"score": dope_score, "interval": bin_id})
            # Actually assigns the DOPE score to the residues of the PyMod element.
            for res, dope_item in zip(chain_element.get_polymer_residues(), dope_items):
                res.dope_score = dope_item
            chain_element._has_dope_scores = True


    def prepare_dope_plot_data(self, selection, start_from=0, mode="single"):
        """
        Takes a selection of 'PyMod_elemet' objects, takes their DOPE scores and returns the data in
        a dictionary which can be supplied as an argument to the 'show_dope_plot()' in order to
        display it in a plot.
        """
        dope_plot_data = []
        for element in selection:
            # Prepares a list with the PyMOL additional data for each residue of the chain, so that
            # when clicking on some point, the corresponding residue will be highlighted in PyMOL,
            # and the message bar of the plot will be updated.
            residues_names = [res.three_letter_code for res in element.get_polymer_residues()]
            residues_pdb_positions = [res.db_index for res in element.get_polymer_residues()]
            pymol_selectors = [res.get_pymol_selector() for res in element.get_polymer_residues()]
            residue_additional_data = []
            for r_name, r_position, r_selector in zip(residues_names, residues_pdb_positions, pymol_selectors):
                residue_additional_data.append({"residue_name": r_name,
                                                "residue_pdb_position": r_position,
                                                "pymol_selector": r_selector,
                                                "export_label": "%s %s"%(r_name, r_position)})
            element_dope_scores = self.dope_scores_dict[element][:]
            # If the sequences in the selection are aligned, adjust the profiles by inserting 'None'
            # values for gap positions.
            if mode == "multiple":
                # Insert gaps into the profile corresponding to those in seq:
                # R: seq = str(seq).replace("X","-")
                ri = 0
                seq = element.my_sequence
                for i, res in enumerate(seq):
                    if res != "-":
                        # If the first residue is preceeded by some indels.
                        first_residue_with_preceeding_gaps = False
                        if ri == 0 and i != 0:
                            n_gaps = i
                            for gap in range(n_gaps):
                                element_dope_scores.insert(ri, None)
                                residue_additional_data.insert(ri, {"export_label": "Gap"})
                            ri += 1 + n_gaps
                            first_residue_with_preceeding_gaps = True
                        # Applies to the rest of residues in the sequence.
                        if not first_residue_with_preceeding_gaps:
                            n_gaps = pmsm.get_leading_gaps(seq, i)
                            for gap in range(n_gaps):
                                element_dope_scores.insert(ri+1, None)
                                residue_additional_data.insert(ri+1, {"export_label": "Gap"})
                            ri += 1 + n_gaps

            # For DOPE plots of multiple chains models and templates.
            for g in range(start_from):
                element_dope_scores.insert(0, None)
                residue_additional_data.insert(0, {None: None})

            # Prepares the data.
            dope_plot_data.append({"dope_scores": element_dope_scores,
                                   "additional_data": residue_additional_data,
                                   "label": element.compact_header})

        return dope_plot_data


    def show_plot(self):
        show_dope_plot(self.dope_plot_data, self.pymod.get_qt_parent(), pymod=self.pymod)


class DOPE_profile_window_session:
    """
    A class used to show DOPE plots. The constructor takes as first argument the output of the
    'prepare_dope_plot_data' method of the 'DOPE_assessment' class.
    """

    def __init__(self, dope_plot_data, parent_window, pymod=None):
        self.dope_plot_data = dope_plot_data
        self.parent_window = parent_window
        self.pymod = pymod

    def show(self):
        """
        Uses the 'pymod_plot' module to show a DOPE profile.
        """
        x_label_text = None
        messagebar_text_on_update = None
        if len(self.dope_plot_data) > 1:
            x_label_text = "Alignment position"
            messagebar_text_on_update_qt = "Selected: __residue_name__ __residue_pdb_position__ of __plot_name__ (alignment position: __x__), DOPE value: __y__"
        else:
            x_label_text = "Residue position"
            messagebar_text_on_update_qt = "Selected: __residue_name__ __residue_pdb_position__ of __plot_name__, DOPE value: __y__"

        cp = PyMod_plot_window_qt(self.pymod.main_window)
        cp.initialize(pymod=self.pymod, title="DOPE Profile")
        cp.build_plotting_area(messagebar_initial_text="Set the 'Interact on click' option to 'Yes' to click on the plot and to highlight residues in PyMOL.",
                               update_messagebar=True,
                               messagebar_text_on_update=messagebar_text_on_update_qt,
                               on_click_action=self._highlight_in_pymol_from_dope_plot_qt,
                               x_label_text=x_label_text,
                               y_label_text="DOPE score")
        for chain_dope_data in self.dope_plot_data:
            # Adds the new plot corresponding to the current chain.
            x_data = list(range(1, len(chain_dope_data["dope_scores"])+1))
            y_data = chain_dope_data["dope_scores"]
            cp.add_plot(x_data, y_data,
                        label=chain_dope_data["label"],
                        additional_data=chain_dope_data["additional_data"])
        cp.show()


    def _highlight_in_pymol_from_dope_plot_qt(self, point_data):
        cmd.select("pymod_selection", point_data["pymol_selector"])
        cmd.center("pymod_selection")


def show_dope_plot(dope_plot_data, parent_window, pymod=None):
    """
    Shortcut function to use the 'DOPE_profile_window_session' class.
    """
    dpw = DOPE_profile_window_session(dope_plot_data, parent_window, pymod=pymod)
    dpw.show()


def compute_dope_of_structure_file(pymod, str_file_path, profile_file_path, env=None):
    """
    Quickly computes the DOPE energy of the molecules contained in a structure file.
    """
    model_dope_protocol = DOPE_assessment(pymod)
    dope_score = model_dope_protocol._compute_dope_of_structure_file(str_file_path, profile_file_path, env=env)
    return dope_score
