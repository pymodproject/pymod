# Copyright 2025 by me. All rights reserved.
# This code is part of the PyMod package and governed by its license.
#
# This file is adapted from clash_refine_gui.py for single-molecule optimization.

"""
A module for the Molecule Optimizer configuration window.
"""

from pymol.Qt import QtWidgets
from pymol import cmd

class MoleculeOptimizerWindow(QtWidgets.QDialog):
    """
    A class for the Molecule Optimizer configuration window.
    """

    def __init__(self, parent=None, pymod_elements=None, has_openbabel=False):
        super(MoleculeOptimizerWindow, self).__init__(parent)

        self.setWindowTitle("Molecule Optimizer")

        layout = QtWidgets.QGridLayout(self)

        # --- SELECTION ---
        self.selection_label = QtWidgets.QLabel("Molecule to Optimize:")
        self.selection_combobox = QtWidgets.QComboBox()
        self.selection_combobox.setToolTip("Select the single molecule you wish to optimize.\nAll other molecules will be treated as a fixed environment.")
        self.update_selection_options(pymod_elements)
        layout.addWidget(self.selection_label, 0, 0)
        layout.addWidget(self.selection_combobox, 0, 1)

        # --- CLASH REFINEMENT CHECKBOX ---
        self.resolve_clashes_checkbox = QtWidgets.QCheckBox("Find and resolve clashes with sculpting")
        self.resolve_clashes_checkbox.setToolTip("If checked, the molecule will be scanned for internal clashes, which will then be repaired using sculpting.")
        self.resolve_clashes_checkbox.setChecked(False) # Default to off
        layout.addWidget(self.resolve_clashes_checkbox, 1, 0, 1, 2)

        # --- ADVANCED SETTINGS ---
        self.advanced_button = QtWidgets.QPushButton("Show Advanced Settings")
        self.advanced_button.setCheckable(True)
        self.advanced_button.toggled.connect(self.toggle_advanced_settings)
        layout.addWidget(self.advanced_button, 2, 0, 1, 2)

        # --- PARAMETERS GROUP (Initially hidden) ---
        self.params_group = QtWidgets.QGroupBox("Clash Detection Parameters")
        params_layout = QtWidgets.QGridLayout(self.params_group)
        layout.addWidget(self.params_group, 3, 0, 1, 2)

        self.threshold_label = QtWidgets.QLabel("Clash Threshold:")
        self.threshold_spinbox = QtWidgets.QDoubleSpinBox()
        self.threshold_spinbox.setToolTip("VDW radii multiplication factor. Clash if distance < (VDW1 + VDW2) * threshold.")
        self.threshold_spinbox.setRange(0.1, 1.0); self.threshold_spinbox.setValue(0.85); self.threshold_spinbox.setSingleStep(0.05)
        params_layout.addWidget(self.threshold_label, 0, 0)
        params_layout.addWidget(self.threshold_spinbox, 0, 1)

        self.min_overlap_label = QtWidgets.QLabel("Min Overlap (Å):")
        self.min_overlap_spinbox = QtWidgets.QDoubleSpinBox()
        self.min_overlap_spinbox.setToolTip("Minimum overlap in Angstroms to be considered a clash.")
        self.min_overlap_spinbox.setRange(0.0, 2.0); self.min_overlap_spinbox.setValue(0.1); self.min_overlap_spinbox.setSingleStep(0.05)
        params_layout.addWidget(self.min_overlap_label, 1, 0)
        params_layout.addWidget(self.min_overlap_spinbox, 1, 1)

        # --- SCULPTING GROUP (Initially hidden) ---
        self.sculpting_group = QtWidgets.QGroupBox("Sculpting Parameters")
        sculpting_layout = QtWidgets.QGridLayout(self.sculpting_group)
        layout.addWidget(self.sculpting_group, 4, 0, 1, 2)

        self.sculpting_radius_label = QtWidgets.QLabel("Sculpting Radius (Å):")
        self.sculpting_radius_spinbox = QtWidgets.QDoubleSpinBox()
        self.sculpting_radius_spinbox.setToolTip("Radius for the sculpting sphere around clashing residues.")
        self.sculpting_radius_spinbox.setRange(1.0, 10.0); self.sculpting_radius_spinbox.setValue(5.0); self.sculpting_radius_spinbox.setSingleStep(0.5)
        sculpting_layout.addWidget(self.sculpting_radius_label, 0, 0)
        sculpting_layout.addWidget(self.sculpting_radius_spinbox, 0, 1)

        self.sculpting_cycles_label = QtWidgets.QLabel("Sculpting Cycles:")
        self.sculpting_cycles_spinbox = QtWidgets.QSpinBox()
        self.sculpting_cycles_spinbox.setToolTip("Number of sculpting cycles to perform if clashes are found.")
        self.sculpting_cycles_spinbox.setRange(1, 1000); self.sculpting_cycles_spinbox.setValue(100)
        sculpting_layout.addWidget(self.sculpting_cycles_label, 1, 0)
        sculpting_layout.addWidget(self.sculpting_cycles_spinbox, 1, 1)

        # --- OPEN BABEL GROUP (Initially hidden) ---
        self.ob_group = QtWidgets.QGroupBox("Minimization Parameters")
        ob_layout = QtWidgets.QGridLayout(self.ob_group)
        layout.addWidget(self.ob_group, 5, 0, 1, 2)

        if not has_openbabel:
            self.ob_group.setTitle("Minimization Parameters (Open Babel not installed)")
            self.ob_group.setEnabled(False)

        self.ob_ff_label = QtWidgets.QLabel("Force Field:")
        self.ob_ff_combo = QtWidgets.QComboBox()
        self.ob_ff_combo.addItems(['MMFF94s', 'GAFF', 'MMFF94', 'UFF', 'Ghemical']); self.ob_ff_combo.setCurrentText("MMFF94s")
        ob_layout.addWidget(self.ob_ff_label, 0, 0)
        ob_layout.addWidget(self.ob_ff_combo, 0, 1)

        self.ob_nsteps_label = QtWidgets.QLabel("Max Steps:")
        self.ob_nsteps_spinbox = QtWidgets.QSpinBox()
        self.ob_nsteps_spinbox.setToolTip("Maximum number of minimization steps.")
        self.ob_nsteps_spinbox.setRange(1, 20000); self.ob_nsteps_spinbox.setValue(300)
        ob_layout.addWidget(self.ob_nsteps_label, 1, 0)
        ob_layout.addWidget(self.ob_nsteps_spinbox, 1, 1)

        self.ob_cutoff_checkbox = QtWidgets.QCheckBox("Enable Cutoff")
        self.ob_cutoff_checkbox.setToolTip("Enable/disable cutoff for non-bonded interactions.")
        self.ob_cutoff_checkbox.setChecked(True)
        ob_layout.addWidget(self.ob_cutoff_checkbox, 2, 0, 1, 2)

        self.ob_vdw_cutoff_label = QtWidgets.QLabel("VDW Cutoff (Å):")
        self.ob_vdw_cutoff_spinbox = QtWidgets.QDoubleSpinBox()
        self.ob_vdw_cutoff_spinbox.setToolTip("Cutoff in Angstroms for Van der Waals interactions.")
        self.ob_vdw_cutoff_spinbox.setRange(1.0, 20.0)
        self.ob_vdw_cutoff_spinbox.setSingleStep(0.5)
        self.ob_vdw_cutoff_spinbox.setValue(6.0)
        ob_layout.addWidget(self.ob_vdw_cutoff_label, 3, 0)
        ob_layout.addWidget(self.ob_vdw_cutoff_spinbox, 3, 1)

        self.ob_elec_cutoff_label = QtWidgets.QLabel("Electrostatic Cutoff (Å):")
        self.ob_elec_cutoff_spinbox = QtWidgets.QDoubleSpinBox()
        self.ob_elec_cutoff_spinbox.setToolTip("Cutoff in Angstroms for electrostatic interactions.")
        self.ob_elec_cutoff_spinbox.setRange(1.0, 20.0)
        self.ob_elec_cutoff_spinbox.setSingleStep(0.5)
        self.ob_elec_cutoff_spinbox.setValue(8.0)
        ob_layout.addWidget(self.ob_elec_cutoff_label, 4, 0)
        ob_layout.addWidget(self.ob_elec_cutoff_spinbox, 4, 1)
        
        self.keep_original_checkbox = QtWidgets.QCheckBox("Keep a copy of the original structure")
        self.keep_original_checkbox.setChecked(True)
        ob_layout.addWidget(self.keep_original_checkbox, 5, 0, 1, 2)

        # --- BUTTONS ---
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.button(QtWidgets.QDialogButtonBox.Ok).setText("Optimize")
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        layout.addWidget(self.button_box, 6, 0, 1, 2)

        # --- INITIAL STATE ---
        self.toggle_advanced_settings(False)
        self.adjustSize()

    def toggle_advanced_settings(self, checked):
        self.params_group.setVisible(checked)
        self.sculpting_group.setVisible(checked)
        self.advanced_button.setText("Hide Advanced Settings" if checked else "Show Advanced Settings")
        # Keep minimization options always visible but grouped
        # self.ob_group.setVisible(checked)
        self.adjustSize()

    def show_forcefield_fallback_dialog(self, failed_forcefield, parent=None):
        """
        Show a dialog to let the user choose an alternative force field when the initial one fails.
        Returns the selected force field name or None if user cancels.
        Similar to the implementation in clash_and_refine.py
        """
        # Available force fields in order of recommendation
        available_forcefields = ['GAFF', 'MMFF94', 'UFF', 'Ghemical']
        
        # Remove the failed force field from the list if it's there
        if failed_forcefield in available_forcefields:
            available_forcefields.remove(failed_forcefield)
        
        if not available_forcefields:
            # No alternatives available
            return None
        
        # Create a simple dialog
        msg_box = QtWidgets.QMessageBox(parent or self)
        msg_box.setWindowTitle("Force Field Setup Failed")
        msg_box.setText(f"Force field setup failed for '{failed_forcefield}'.")
        msg_box.setInformativeText("Would you like to try with a different force field?")
        msg_box.setIcon(QtWidgets.QMessageBox.Warning)
        
        # Add buttons for each available force field
        buttons = {}
        for ff in available_forcefields:
            button = msg_box.addButton(ff, QtWidgets.QMessageBox.ActionRole)
            buttons[button] = ff
        
        # Add cancel button
        cancel_button = msg_box.addButton("Cancel", QtWidgets.QMessageBox.RejectRole)
        msg_box.setDefaultButton(cancel_button)
        
        # Show dialog and get result
        msg_box.exec_()
        clicked_button = msg_box.clickedButton()
        
        if clicked_button in buttons:
            selected_ff = buttons[clicked_button]
            print(f"User selected alternative force field: {selected_ff}")
            return selected_ff
        else:
            print("User cancelled force field selection")
            return None

    def update_selection_options(self, pymod_elements=None):
        """Update selection options similar to clash_refine_gui.py"""
        self.selection_combobox.clear()
        elements = pymod_elements or []

        if len(elements) >= 1:
            # Add PyMod objects
            element_selectors = [e.get_pymol_selector() for e in elements]
            self.selection_combobox.addItems(element_selectors)
            
            # Add relevant PyMOL selections for each object
            for selector in element_selectors:
                all_selections = [name for name in cmd.get_names("selections") if cmd.get_type(name) == 'selection']
                relevant_selections = [s for s in all_selections if cmd.count_atoms(f"({s}) and ({selector})") > 0 and s != selector]
                if relevant_selections:
                    for sel in relevant_selections:
                        # Add with object context to avoid duplicates
                        display_name = f"{sel} (from {selector})"
                        self.selection_combobox.addItem(display_name)
        else:
            # Fallback: add all available selections or objects with context
            all_objects = cmd.get_names("objects")
            all_selections = [name for name in cmd.get_names("selections") if cmd.get_type(name) == 'selection']
            
            items = set(all_objects)
            
            for sel in all_selections:
                is_processed = False
                for obj in all_objects:
                    if sel != obj and cmd.count_atoms(f"({sel}) and ({obj})") > 0:
                        items.add(f"{sel} (from {obj})")
                        is_processed = True
                if not is_processed:
                    items.add(sel)

            if items:
                self.selection_combobox.addItems(sorted(list(items)))
            else:
                self.selection_combobox.addItem("No objects found")
                self.selection_combobox.setEnabled(False)

    def get_values(self):
        # Parse the selection text to extract the actual selection name
        selected_text = self.selection_combobox.currentText()
        mobile_object = selected_text
        
        # If it's a selection with context (e.g., "selection_name (from object_name)"),
        # extract just the selection name -- and keep the object it was listed under, which
        # is the one the protocol has to optimize.
        mobile_parent = None
        if " (from " in selected_text and selected_text.endswith(")"):
            mobile_object = selected_text.split(" (from ")[0]
            mobile_parent = selected_text.split(" (from ")[1][:-1]

        vals = {
            "mobile_object": mobile_object,
            "mobile_parent": mobile_parent,
            "keep_original": self.keep_original_checkbox.isChecked(),
            "resolve_clashes": self.resolve_clashes_checkbox.isChecked(),
            # Clash detection
            "threshold": self.threshold_spinbox.value(),
            "min_overlap": self.min_overlap_spinbox.value(),
            # Sculpting
            "sculpting_radius": self.sculpting_radius_spinbox.value(),
            "sculpting_cycles": self.sculpting_cycles_spinbox.value(),
        }
        if self.ob_group.isEnabled():
            vals["run_ob_minimize"] = True
            vals["ob_forcefield"] = self.ob_ff_combo.currentText()
            vals["ob_nsteps"] = self.ob_nsteps_spinbox.value()
            vals["ob_cutoff"] = self.ob_cutoff_checkbox.isChecked()
            vals["ob_vdw_cutoff"] = self.ob_vdw_cutoff_spinbox.value()
            vals["ob_elec_cutoff"] = self.ob_elec_cutoff_spinbox.value()
        else:
            vals["run_ob_minimize"] = False
            
        return vals
