# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Interactions of PyMod with PyMOL.
"""

import os
from pymol import cmd

from pymod_lib.pymod_gui.specific_gui_components_qt import Import_from_pymol_window_qt, Multistate_import_dialog_qt
from pymod_lib.pymod_vars import sphere_hetres_names


class PyMod_pymol_interactions:

    def _safe_pymol_operation(self, operation, *args, **kwargs):
        """
        Safely execute PyMOL operations to prevent crashes when windows overlap.
        
        This wrapper ensures that PyMOL operations don't cause conflicts with the
        PyMod GUI, especially when windows are overlapping or during object loading.
        """
        try:
            # Notify main window that a PyMOL operation is starting
            if hasattr(self, 'main_window') and hasattr(self.main_window, '_handling_pymol_conflict'):
                was_handling = self.main_window._handling_pymol_conflict
                self.main_window._handling_pymol_conflict = True
                
            try:
                # Execute the PyMOL operation
                result = operation(*args, **kwargs)
                
                # Brief pause to let PyMOL settle before returning control
                from pymol.Qt import QtCore
                QtCore.QThread.msleep(100)  # 100ms pause
                
                return result
                
            finally:
                # Restore the handling flag
                if hasattr(self, 'main_window') and hasattr(self.main_window, '_handling_pymol_conflict'):
                    self.main_window._handling_pymol_conflict = was_handling
                    # Ensure window separation after PyMOL operation
                    if hasattr(self.main_window, '_ensure_window_separation'):
                        QtCore.QTimer.singleShot(200, self.main_window._ensure_window_separation)
                        
        except Exception as e:
            if self.DEVELOP:
                print(f"Error in safe PyMOL operation: {e}")
            # Don't let PyMOL errors crash PyMod
            if hasattr(self, 'main_window'):
                self.main_window.show_error_message("PyMOL Operation Error", 
                                                   f"An error occurred during PyMOL operation: {str(e)}")
            return None

    ###############################################################################################
    # Basic PyMOL interactions.                                                                   #
    ###############################################################################################

    def load_element_in_pymol(self, element, as_hedgehog=True):
        """
        Loads the PDB structure of the chain into PyMol.
        Protected against crashes during PyMOL operations.
        """
        self._safe_pymol_operation(self._load_element_in_pymol_impl, element, as_hedgehog)
        
    def _load_element_in_pymol_impl(self, element, as_hedgehog=True):
        """
        Internal implementation of load_element_in_pymol.
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
        self._safe_pymol_operation(cmd.center, pymod_element.get_pymol_selector())

    def hide_chain_in_pymol(self, pymod_element):
        self._safe_pymol_operation(cmd.disable, pymod_element.get_pymol_selector())

    def show_chain_in_pymol(self, pymod_element):
        self._safe_pymol_operation(cmd.enable, pymod_element.get_pymol_selector())

    def show_hedgehog_in_pymol(self, pymod_element):
        self._safe_pymol_operation(self._show_hedgehog_in_pymol_impl, pymod_element)
        
    def _show_hedgehog_in_pymol_impl(self, pymod_element):
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
        molecule_objects = [str(obj) for obj in cmd.get_names("objects")
                            if not obj in struct_list
                            and cmd.get_type(obj) == "object:molecule"]

        # Only the ones which hold a chain. PyMod imports a structure as a
        # sequence with coordinates, and a molecule with no chain in it -- a
        # ligand, an ion, something drawn in the builder of PyMOL -- is not a
        # sequence: it was listed here all the same, among the structures, with
        # nothing to say which was which. Nucleic acids are chains too, and are
        # still offered: PyMod models them.
        scrolledlist_items = [obj for obj in molecule_objects
                              if self.holds_a_chain(obj)]

        if not scrolledlist_items:
            if molecule_objects:
                self.main_window.show_error_message(
                    "No Importable Object",
                    "The PyMOL objects which are not already in PyMod hold no"
                    " chain: %s.\n\nPyMod imports structures -- a sequence with"
                    " its coordinates -- so a ligand, an ion or a molecule drawn"
                    " in PyMOL is not offered here. A ligand comes in with the"
                    " structure it belongs to."
                    % ", ".join(molecule_objects[:6]))
            elif struct_list:
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


    @staticmethod
    def holds_a_chain(pymol_object):
        """
        Whether a PyMOL object has a polymer in it -- a protein or a nucleic
        acid -- rather than being a molecule on its own.

        Asked of PyMOL, which knows: a ligand, an ion and a molecule from the
        builder all answer zero, and a protein answers with its atoms whether or
        not it carries ligands of its own.
        """

        try:
            return bool(cmd.count_atoms("(%s) and polymer" % pymol_object))
        except Exception:
            # An object which cannot be asked is left in the list rather than
            # being hidden: the import itself will say what is wrong with it.
            return True


    def import_selected_pymol_object(self):
        selections_to_import = self.import_from_pymol_window.get_objects_to_import()
        if len(selections_to_import) > 0:
            for selected_num, sele in enumerate(selections_to_import):
                selected_num+=1
                
                # CHECK FOR MULTI-STATE OBJECTS
                num_states = cmd.count_states(sele)
                print(f"Object '{sele}' has {num_states} state(s)")
                
                if num_states > 1:
                    # MULTI-STATE OBJECT: Ask user what to do
                    choice = self.handle_multistate_import(sele, num_states)
                    if choice == "cancel":
                        continue  # Skip this object
                    elif choice == "first_only":
                        # Import only first state
                        filename = sele + ".pdb"
                        pdb_file_shortcut = os.path.join(self.structures_dirpath, filename)
                        cmd.save(pdb_file_shortcut, sele, state=1)
                        cmd.delete(sele)
                        self.open_structure_file(os.path.abspath(pdb_file_shortcut))
                    elif choice == "all_states":
                        # Import all states as separate objects
                        for state_num in range(1, num_states + 1):
                            filename = f"{sele}_state_{state_num:03d}.pdb"
                            pdb_file_shortcut = os.path.join(self.structures_dirpath, filename)
                            cmd.save(pdb_file_shortcut, sele, state=state_num)
                            self.open_structure_file(os.path.abspath(pdb_file_shortcut))
                        cmd.delete(sele)
                    elif isinstance(choice, int):
                        # Import specific state chosen by user
                        filename = f"{sele}_state_{choice:03d}.pdb"
                        pdb_file_shortcut = os.path.join(self.structures_dirpath, filename)
                        cmd.save(pdb_file_shortcut, sele, state=choice)
                        cmd.delete(sele)
                        self.open_structure_file(os.path.abspath(pdb_file_shortcut))
                else:
                    # SINGLE-STATE OBJECT: Import normally
                    filename = sele + ".pdb"
                    pdb_file_shortcut = os.path.join(self.structures_dirpath, filename)
                    cmd.save(pdb_file_shortcut, sele)
                    cmd.delete(sele)
                    self.open_structure_file(os.path.abspath(pdb_file_shortcut))
                    
            self.import_from_pymol_window.destroy()
            self.main_window.gridder(update_elements=True)
        else:
            self.main_window.show_error_message("Selection Error", "Please select at least one object to import.")

    def handle_multistate_import(self, object_name, num_states):
        """
        Handle import of multi-state objects (e.g., NMR structures).
        Shows dialog asking user how to import the object.
        
        Returns:
        - "cancel": User cancelled the operation
        - "first_only": Import only the first state
        - "all_states": Import all states as separate objects 
        - int: Import specific state number chosen by user
        """
        dialog = Multistate_import_dialog_qt(
            parent=self.main_window,
            object_name=object_name,
            num_states=num_states,
            title=f"Multi-State Object: {object_name}"
        )
        
        if dialog.exec_():
            return dialog.get_user_choice()
        else:
            return "cancel"
