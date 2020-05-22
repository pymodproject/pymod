# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Common functions and classes for the structural analysis protocols in PyMod.
"""

import math
import numpy as np

from pymol import cmd


class Structural_analysis_mixin:
    """
    Mixin class used in the 'Contact Map' and 'Structural Divergence Plot' protocols.
    """

    def check_protein_selection(self, error_title="Selection Error"):
        if any([e.polymer_type == "nucleic_acid" for e in self.target_sequences]):
            self.pymod.main_window.show_error_message(error_title,
                "Can not perform the analysis for nucleic acids structures.")
            return False
        return True


    def check_multiple_selection(self, error_title="Selection Error"):

        # Check that all elements are currently aligned.
        if not all([e.is_child() for e in self.target_sequences]):
            self.pymod.main_window.show_error_message(error_title,
                "Not all the selected elements are currently in an alignment.")
            return False

        # Check that all elements are currently in the same alignment.
        if len(set([e.mother.unique_index for e in self.target_sequences])) != 1:
            self.pymod.main_window.show_error_message(error_title,
                "Not all the selected elements are currently in the same alignment.")
            return False

        # Check that all elements have a structure loaded in PyMod.
        if not all([e.has_structure() for e in self.target_sequences]):
            self.pymod.main_window.show_error_message(error_title,
                "Please select only elements having a 3D structure loaded in PyMOL.")
            return False

        # Check the alignment length.
        if len(set([len(element.my_sequence) for element in self.target_sequences])) != 1:
            self.pymod.main_window.show_error_message(error_title,
                ("Alignment error: the selected elements do not have aligned sequences of"
                 " the same length. Try to align the elements again or to modify the alignment."))
            return False

        return True


    def _get_sda_data(self, pymod_element, interaction_center="ca"):
        """
        This method takes a 'pymod_element' as argument and returns two dictionaries.
        'ali_seq_dict' associates the alignment ids to the corresponding PyMod residues.
        'res_dict' associates PyMod residues to the atomic coordinates of the corresponding
        interaction center (specified in the 'interaction_center' argument) and their PyMOL selector.
        """
        rc = 0
        ali_seq_dict = {}
        residues = pymod_element.get_polymer_residues()
        for ali_pos, p in enumerate(pymod_element.my_sequence):
            if p != "-":
                ali_seq_dict[ali_pos] = residues[rc]
                rc += 1
        residues, coords_array, selectors = self.get_coords_array(pymod_element, interaction_center)
        res_dict = dict([(res, (c, s)) for (res, c, s) in zip(residues, coords_array, selectors)])
        return ali_seq_dict, res_dict


    def _get_res(self, ali_pos, ali_seq_dict):
        try:
            return ali_seq_dict[ali_pos]
        except KeyError:
            return None

    def _get_coords(self, res, res_dict):
        try:
            return res_dict[res][0]
        except KeyError:
            return None

    def _get_selector(self, res, res_dict):
        try:
            return res_dict[res][1]
        except KeyError:
            return None


def get_distance(coords_1, coords_2):
    return math.sqrt((coords_1[0]-coords_2[0])**2 + (coords_1[1]-coords_2[1])**2 + (coords_1[2]-coords_2[2])**2)
