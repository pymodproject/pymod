# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing a base class for "PyMod protocols". A "PyMod protocol" can be many kind of
actions performed in PyMod (for example, running a BLAST search, performing a multiple sequence
alignment or using MODELLER). Each protocol will be represented by a class derived from this base
class.
"""

import os
import sys
import math
import time # For development.

import numpy as np

# This will give a message if MODELLER is installed, but the licence key is not inserted.
try:
    import modeller
except:
    pass
from io import StringIO

from pymol import cmd

from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt


class PyMod_protocol:
    """
    A base class for PyMod protocols.
    """

    def __init__(self, pymod, output_directory=os.path.curdir):
        # 'PyMod' class object, used to access all the information of the plugin.
        self.pymod = pymod
        # self.protocol_name = protocol_name

        # Original stdout.
        self.sys_stdout = sys.stdout
        # Temporary stdout used by some protocols.
        self.my_stdout = None
        # Directory where the output files of the protocol will be built.
        self.output_directory = output_directory
        # Perform an additional initialization, which is protocol specific.
        self.additional_initialization()


    def additional_initialization(self):
        """
        Overridden in child classes.
        """
        pass


    def get_pymod_elements(self, pymod_elements=None):
        if pymod_elements == None:
            pymod_elements = self.pymod.get_selected_sequences()
        return pymod_elements


    def get_seq_text(self, sequences, elements_string="sequence"):
        n_seqs = len(sequences)
        if n_seqs == 1:
            seq_text = "1 %s" % (elements_string)
        else:
            seq_text = "%s %ss" % (n_seqs, elements_string)
        return seq_text


    def launch_from_gui(self):
        """
        Override in children classes. This should contain a list fof methods
        organized like this (but of course variations are possible and actually
        used in PyMod):
            self.check_parameters_from_gui()
            self.show_options_window()
            self.execute_protocol()
            self.remove_temp_files()
            self.import_results_in_pymod()
        """

        pass


    def quit_protocol(self):
        """
        Quit the protocol. Each protocol may extend this method.
        """
        os.chdir(self.pymod.current_project_dirpath)


    ###############################################################################################
    # Build selection of PyMod elements for the protocols.                                        #
    ###############################################################################################

    def extend_selection_to_hidden_children(self):
        selected_elements = self.pymod.get_selected_sequences()
        extend_selection = None
        # Checks if in the selection there are some cluster leads of which mothers are not selected.
        if False in [e.mother.selected for e in selected_elements if self.pymod.main_window.is_lead_of_collapsed_cluster(e)]:
            # Ask to include the hidden children of the collapsed clusters.
            title = "Selection Message"
            message = "Would you like to include in the alignment the hidden sequences of the collapsed clusters?"
            extend_selection = askyesno_qt(title, message, parent=self.pymod.get_qt_parent())
        else:
            extend_selection = False
        if not extend_selection:
            return None
        # Actually selects the hidden children.
        for e in selected_elements:
            if self.pymod.main_window.is_lead_of_collapsed_cluster(e) and not e.mother.selected:
                e.mother.widget_group.select_collapsed_cluster_descendants()


    def build_cluster_lists(self):
        """
        This will build the self.involved_clusters_list, which will contain the elements
        belonging to cluster that were either selected entirely or with at least one selected child.
        """
        # A list that will contain all the clusters that have at least one selected element.
        self.involved_clusters_list = []
        for e in self.pymod.get_selected_elements():
            if not e.mother in self.involved_clusters_list:
                self.involved_clusters_list.append(e.mother)
        if self.pymod.root_element in self.involved_clusters_list:
            self.involved_clusters_list.remove(self.pymod.root_element)

        # A list that will contain all the clusters that have all of their sequences selected.
        self.selected_clusters_list = self.pymod.get_selected_clusters()

        # A list that will contain all the selected sequences in the root level of PyMod.
        self.selected_root_sequences_list = set([s for s in self.pymod.get_selected_sequences() if s.mother == self.pymod.root_element])


    ###############################################################################################
    # PyMOL related methods.                                                                      #
    ###############################################################################################

    def superpose_in_pymol(self, mobile_selection, fixed_selection, save_superposed_structure=True, output_directory=None):
        """
        Superpose 'mobile' to 'fixed' in PyMOL.
        """
        if not output_directory:
            output_directory = self.pymod.structures_dirpath
        if hasattr(cmd, "super"): # 'super' is sequence-independent.
            cmd.super(mobile_selection, fixed_selection)
        else: # PyMOL 0.99 does not have 'cmd.super'.
            cmd.align(mobile_selection, fixed_selection)
        if save_superposed_structure:
            cmd.save(os.path.join(output_directory, mobile_selection+".pdb"), mobile_selection)


    def get_coords_array(self, pymod_element, interaction_center, get_selectors=True):
        """
        Takes a 'pymod_element' as input and returns a list of its PyMod residues and the corresponding
        atomic coordinates of the interaction center (specified 'interaction_center') and selectors.
        Only residues for which atomic coordinates are defined in PyMOL are retrieved.
        """

        objsel = pymod_element.get_pymol_selector()
        residues = pymod_element.get_polymer_residues()

        if interaction_center == "ca":
            atm_types = ("CA", )
        elif interaction_center == "cb":
            atm_types = ("CB", "CA", )
        else:
            raise KeyError("Invalid 'interaction_center': %s" % interaction_center)

        # Gets the coordinates from PyMOL.
        from pymol import stored

        # Builds a dictionary in which the keys are tuples with (residue_id, atom_name) and the values
        # are coordinates.
        stored.pymod_element_coords = {}

        # The 'cmd.iterate_state' function is particularly fast in PyMOL, much faster than iterating to
        # every residue and calling 'cmd.get_coords' for them.
        cmd.iterate_state(1, objsel, "stored.pymod_element_coords[(int(resv), name)] = [x, y, z]")

        # Store the coordinates of residues having an interaction center atom (CA or CB).
        coords = [] # Coordinates list.
        _residues = [] # New residues list, which will substitute the 'residues' one.
        selectors = []
        for res in residues:

            # Attempts to get various atoms.
            for atm_type in atm_types:
                if (res.db_index, atm_type) in stored.pymod_element_coords:
                    coords_i = stored.pymod_element_coords[(res.db_index, atm_type)]
                    coords.append(coords_i)
                    _residues.append(res)
                    if get_selectors:
                        sel = "object %s and n. %s and i. %s" % (objsel, atm_type, res.db_index)
                        selectors.append(sel)
                    break

        coords = np.array(coords)

        if get_selectors:
            return _residues, coords, selectors
        else:
            return _residues, coords


    def get_rmsd(self, element_1, element_2, coords_dict_1=None, coords_dict_2=None,
                 threshold_dist=100.0):
        """
        Takes two 'PyMod_elements' objects and computes a RMSD deviation between their structures
        loaded in PyMOL. The RMSD is computed between a list of residues pairs defined by the
        alignment currently existing in PyMod between the two sequences. Remove aligned pairs with
        a distance above 'threshold_dist'.
        """

        t1 = time.time()

        # Get the list of the ids of the aligned residues for both sequences.
        ali_id = 0
        coords_1 = []
        coords_2 = []
        for id_1, id_2 in zip(element_1.my_sequence, element_2.my_sequence):
            if id_1 != "-" and id_2 != "-":
                db_index_1 = element_1.get_residue_by_index(ali_id, aligned_sequence_index=True).db_index
                db_index_2 = element_2.get_residue_by_index(ali_id, aligned_sequence_index=True).db_index
                # Get their coordinates from the coordinate dictionaries and stores them in form of matrices.
                if db_index_1 in coords_dict_1 and db_index_2 in coords_dict_2:
                    coords_1.append(coords_dict_1[db_index_1])
                    coords_2.append(coords_dict_2[db_index_2])
            ali_id += 1
        coords_1 = np.array(coords_1)
        coords_2 = np.array(coords_2)

        # Perform substraction between the two matrices.
        coords_diff = coords_2 - coords_1
        # Take the square of every element in the matrix and sum along rows.
        square_dist_list = (coords_diff**2).sum(axis=1)
        # Filters the distances > 'threshold_dist' Angstrom.
        square_dist_list = square_dist_list[square_dist_list < threshold_dist**2]

        # Computes the RMSD.
        square_dist_sum = square_dist_list.sum()
        n_aligned = float(square_dist_list.shape[0])
        if n_aligned == 0:
            rmsd = None
        else:
            if square_dist_sum == 0:
                rmsd = 0.0
            else:
                rmsd = math.sqrt(square_dist_sum/n_aligned)

        return rmsd


###################################################################################################
# MODELLER protocol.                                                                              #
###################################################################################################

class MODELLER_common:
    """
    A base class for running MODELLER based scripts. To be used as a parent class along
    'PyMod_protocol' class.
    """

    def __init__(self):
        self.run_modeller_internally = self.pymod.modeller.run_internally()

    def _initialize_env(self, use_hetatm=True, use_water=True):
        env = modeller.environ()
        env.io.atom_files_directory = []
        env.io.atom_files_directory.append(".")
        env.io.hetatm = use_hetatm
        env.io.water = use_water
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        return env
