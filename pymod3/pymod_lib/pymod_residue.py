# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Residues of sequences and structures in PymMod.
"""

from pymod_lib.pymod_vars import feature_types


class PyMod_residue:
    """
    Class for a residue of PyMod element.
    """

    def __init__(self, three_letter_code, one_letter_code, index=None, seq_index=None, db_index=None):

        self.three_letter_code = three_letter_code
        self.one_letter_code = one_letter_code
        self.full_name = three_letter_code

        self.index = index
        self.seq_index = seq_index
        self.db_index = db_index # "Database" index: for elements with a 3D structure in
                                 # PyMOL, it is the PDB residue index.

        self.pymod_element = None
        self.color_seq = None # Color for the sequence in the main window. Contains an RGB code.
        self.color_str = None # Color for the 3D structure loaded in PyMOL. Contains a color name.
        self.custom_color = None # Used to attribute a residue a particular color.
        # These attributes store the colors used at some time for a residue, so that when temporarily
        # changing its color, it can be colored back to its original color.
        self._color_seq = None
        self._color_str = None
        self._custom_color = None

        self.secondary_structure = None
        self.psipred_result = None
        self.campo_score = None
        self.dope_score = None
        self.scr_score = None

        # Features.
        self.features = dict([(ft, []) for ft in feature_types])


    def __repr__(self):
        return self.one_letter_code+'__'+str(self.index)

    def is_polymer_residue(self):
        """
        Check if the residue is part of a polymer chain, or a ligand molecule/atom.
        """
        return not self.is_water() and not self.is_ligand()

    def is_standard_residue(self):
        return isinstance(self, PyMod_standard_residue)

    def is_water(self):
        return isinstance(self, PyMod_water_molecule)

    def is_ligand(self):
        return isinstance(self, PyMod_ligand)

    def is_modified_residue(self):
        return isinstance(self, PyMod_modified_residue)

    def is_heteroresidue(self, exclude_water=True):
        if exclude_water:
            return self.is_non_water_heteroresidue()
        else:
            return issubclass(self.__class__, PyMod_heteroresidue)

    def is_non_water_heteroresidue(self):
        return issubclass(self.__class__, PyMod_heteroresidue) and not self.is_water()


    def get_pymol_selector(self):
        """
        Gets the correspondig selector in PyMOL.
        """
        # Selectors that work:
        #     #     /1UBI_Chain_A//A/LEU`43/CA
        #     #     1UBI_Chain_A and resi 43
        return "%s and resi %s" % (self.pymod_element.get_pymol_selector(), self.db_index)


    def get_id_in_aligned_sequence(self):
        # Only polymer residues can be included in alignments.
        assert self.is_polymer_residue()
        res_counter = 0
        index = self.pymod_element.get_polymer_residues().index(self)
        for i, p in enumerate(self.pymod_element.my_sequence):
            if p != "-":
                if index == res_counter:
                    return i
                res_counter += 1
        return None


    def get_parent_structure_chain_id(self):
        return self.pymod_element.get_chain_id()


    def get_default_color(self):
        if self.color_str is None:
            return self.pymod_element.my_color
        else:
            return self.color_str

    # Used in the "Search subsequence" feature of PyMod.
    def store_current_colors(self):
        self._color_seq = self.color_seq
        self._color_str = self.color_str
        self._custom_color = self.custom_color

    def revert_original_colors(self):
        self.color_seq = self._color_seq
        self.color_str = self._color_str
        self.custom_color = self._custom_color


class PyMod_standard_residue(PyMod_residue):
    hetres_type = None

class PyMod_heteroresidue(PyMod_residue):
    hetres_type = "?"

class PyMod_ligand(PyMod_heteroresidue):
    hetres_type = "ligand"

class PyMod_modified_residue(PyMod_heteroresidue):
    hetres_type = "modified residue"

class PyMod_water_molecule(PyMod_heteroresidue):
    hetres_type = "water"
