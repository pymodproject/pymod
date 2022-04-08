# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Interactions of PyMod with PyMOL.
"""

import os
from pymol import cmd

from pymod_lib.pymod_gui.specific_gui_components_qt import Import_from_pymol_window_qt
from pymod_lib.pymod_vars import sphere_hetres_names


class PyMod_pymol_interactions:

    ###############################################################################################
    # Basic PyMOL interactions.                                                                   #
    ###############################################################################################

    def load_element_in_pymol(self, element, as_hedgehog=True):
        """
        Loads the PDB structure of the chain into PyMol.
        """
        file_to_load = element.get_structure_file(basename_only=False)
        pymol_object_name = element.get_pymol_selector()
        cmd.load(file_to_load, pymol_object_name)
        cmd.color(element.my_color, pymol_object_name)
        cmd.hide("everything", pymol_object_name)
        cmd.show("cartoon", pymol_object_name) # Show the new chain as a cartoon.
        if as_hedgehog:
            cmd.show("sticks", "(%s and name cb) or (%s and name ca)" % (pymol_object_name, pymol_object_name))
        cmd.util.cnc(pymol_object_name) # Colors by atom.
        cmd.zoom(pymol_object_name)
        cmd.center(pymol_object_name)

    def center_chain_in_pymol(self, pymod_element):
        cmd.center(pymod_element.get_pymol_selector())

    def hide_chain_in_pymol(self, pymod_element):
        cmd.disable(pymod_element.get_pymol_selector())

    def show_chain_in_pymol(self, pymod_element):
        cmd.enable(pymod_element.get_pymol_selector())

    def show_hedgehog_in_pymol(self, pymod_element):
        pymol_object_name = pymod_element.get_pymol_selector()
        cmd.hide("everything", pymol_object_name)
        cmd.show("cartoon", pymol_object_name) # Show the new chain as a cartoon.
        cmd.show("sticks", "(%s and name cb) or (%s and name ca)" % (pymol_object_name, pymol_object_name))

    def show_het_in_pymol(self, pymod_element):
        self._set_het_show_state_in_pymol(pymod_element, show=True)

    def hide_het_in_pymol(self, pymod_element):
        self._set_het_show_state_in_pymol(pymod_element, show=False)

    def _set_het_show_state_in_pymol(self, pymod_element, show):
        pymol_object_name = pymod_element.get_pymol_selector()
        het_residues = pymod_element.get_residues(standard=False, ligands=True, modified_residues=True, water=False)
        for het_res in het_residues:
            resname = het_res.three_letter_code
            res_sel = het_res.get_pymol_selector()
            if self.TEST or self.DEVELOP:
                print(resname, res_sel)
            if resname in sphere_hetres_names:
                if show:
                    cmd.show("spheres", res_sel)
                else:
                    cmd.hide("spheres", res_sel)
            else:
                if show:
                    cmd.show("sticks", res_sel)
                else:
                    cmd.hide("sticks", res_sel)


    ###############################################################################################
    # Import structures from PyMOL.                                                               #
    ###############################################################################################

    def import_pymol_selections_from_main_menu(self):
        """
        Method for importing PyMOL Selections into PyMod. It saves PyMOL objects selected by users
        to file, and loads it into PyMOL using 'open_structure_file()'.
        """
        # Find all structures already loaded into PyMod: items in struct_list are excluded from
        # importable PyMOL object list.
        struct_list = [member.get_pymol_selector() for member in self.get_pymod_elements_list() if member.has_structure()]

        # Importable PyMOL objects.
        scrolledlist_items = [str(obj) for obj in cmd.get_names("objects") if not obj in struct_list and cmd.get_type(obj) == "object:molecule"]

        if not scrolledlist_items:
            if struct_list:
                self.main_window.show_error_message("No Importable Object", "All PyMOL objects are already imported into PyMod.")
            else:
                self.main_window.show_error_message("No Importable Object", "No PyMOL object to import.")
            return

        # Builds a new window.
        self.import_from_pymol_window = Import_from_pymol_window_qt(self.main_window,
            title="Import from PyMOL",
            upper_frame_title="Load PyMOL Objects into PyMod",
            submit_command=self.import_selected_pymol_object,
            selections_list=scrolledlist_items)
        self.import_from_pymol_window.show()


    def import_selected_pymol_object(self):
        selections_to_import = self.import_from_pymol_window.get_objects_to_import()
        if len(selections_to_import) > 0:
            for selected_num, sele in enumerate(selections_to_import):
                selected_num+=1
                filename=sele+".pdb"
                pdb_file_shortcut = os.path.join(self.structures_dirpath, filename)
                cmd.save(pdb_file_shortcut,sele)
                cmd.delete(sele)
                self.open_structure_file(os.path.abspath(pdb_file_shortcut))
            self.import_from_pymol_window.destroy()
            self.main_window.gridder(update_elements=True)
        else:
            self.main_window.show_error_message("Selection Error", "Please select at least one object to import.")
