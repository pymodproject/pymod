# Copyright 2025 by me. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod4 folder.

"""
A module for the Clash & Refine configuration window.
"""

from pymol.Qt import QtWidgets
from pymol import cmd

class ClashRefineWindow(QtWidgets.QDialog):
    """
    A class for the Clash & Refine configuration window.
    """

    def __init__(self, parent=None, selected_elements=None, has_openbabel=False):
        super(ClashRefineWindow, self).__init__(parent)

        self.setWindowTitle("Optimize Interface")
        self.action = None # To store which button was clicked

        layout = QtWidgets.QGridLayout(self)

        # Selection dropdown
        self.selection_label = QtWidgets.QLabel("Selection:")
        self.selection_combobox = QtWidgets.QComboBox()
        self.selection_combobox.setToolTip(
            "When 'Complex' is selected, the interface between all selected objects will be minimized.\n" 
            "When a single object or selection is chosen, only that part will be optimized."
        )
        self.update_selection_options(selected_elements)

        layout.addWidget(self.selection_label, 0, 0)
        layout.addWidget(self.selection_combobox, 0, 1)

        # --- Advanced Settings Button ---
        self.advanced_button = QtWidgets.QPushButton("Show Advanced Settings")
        self.advanced_button.setCheckable(True)
        self.advanced_button.toggled.connect(self.toggle_advanced_settings)
        layout.addWidget(self.advanced_button, 1, 0, 1, 2)


        # --- Parameters Group ---
        self.params_group = QtWidgets.QGroupBox("Detection Parameters")
        params_layout = QtWidgets.QGridLayout(self.params_group)
        layout.addWidget(self.params_group, 2, 0, 1, 2)

        # Clash detection threshold
        self.threshold_label = QtWidgets.QLabel("Clash Detection Threshold:")
        self.threshold_spinbox = QtWidgets.QDoubleSpinBox()
        self.threshold_spinbox.setToolTip("Factor to multiply VDW radii. A clash is detected if distance < (VDW1 + VDW2) * threshold.")
        self.threshold_spinbox.setRange(0.1, 1.0)
        self.threshold_spinbox.setSingleStep(0.05)
        self.threshold_spinbox.setValue(0.85)
        params_layout.addWidget(self.threshold_label, 0, 0)
        params_layout.addWidget(self.threshold_spinbox, 0, 1)

        # Minimum Overlap
        self.min_overlap_label = QtWidgets.QLabel("Minimum Overlap (Å):")
        self.min_overlap_spinbox = QtWidgets.QDoubleSpinBox()
        self.min_overlap_spinbox.setToolTip("Minimum overlap in Angstroms to consider a contact a clash.")
        self.min_overlap_spinbox.setRange(0.0, 2.0)
        self.min_overlap_spinbox.setSingleStep(0.05)
        self.min_overlap_spinbox.setValue(0.1)
        params_layout.addWidget(self.min_overlap_label, 1, 0)
        params_layout.addWidget(self.min_overlap_spinbox, 1, 1)

        # Max Display
        self.max_display_label = QtWidgets.QLabel("Max Bumps to Display:")
        self.max_display_spinbox = QtWidgets.QSpinBox()
        self.max_display_spinbox.setToolTip("Maximum number of clashes to display.")
        self.max_display_spinbox.setRange(1, 1000)
        self.max_display_spinbox.setValue(200)
        params_layout.addWidget(self.max_display_label, 2, 0)
        params_layout.addWidget(self.max_display_spinbox, 2, 1)

        # Ignore H
        self.ignore_h_checkbox = QtWidgets.QCheckBox("Ignore Hydrogen Atoms")
        self.ignore_h_checkbox.setToolTip("If checked, hydrogen atoms are ignored in clash detection.")
        self.ignore_h_checkbox.setChecked(True)
        params_layout.addWidget(self.ignore_h_checkbox, 3, 0, 1, 2)

        # Aggressive Mode
        self.aggressive_checkbox = QtWidgets.QCheckBox("Only Severe Clashes (Aggressive Mode)")
        self.aggressive_checkbox.setToolTip("If checked, uses a more aggressive clash detection preset.")
        self.aggressive_checkbox.toggled.connect(self.toggle_aggressive_mode)
        params_layout.addWidget(self.aggressive_checkbox, 4, 0, 1, 2)

        # --- Refinement Group ---
        self.refine_group = QtWidgets.QGroupBox("Refinement Parameters")
        refine_layout = QtWidgets.QGridLayout(self.refine_group)
        layout.addWidget(self.refine_group, 3, 0, 1, 2)

        # Keep original checkbox
        self.keep_original_checkbox = QtWidgets.QCheckBox("Keep a copy of the original structure")
        self.keep_original_checkbox.setToolTip("If checked, a copy of the original structure is kept.")
        # Checked, as in Optimize Molecule: the refinement replaces the structures it
        # started from, and getting them back means loading them again.
        self.keep_original_checkbox.setChecked(True)
        refine_layout.addWidget(self.keep_original_checkbox, 0, 0, 1, 2)
        
        # --- Sculpting Group ---
        self.sculpting_group = QtWidgets.QGroupBox("Sculpting Parameters")
        self.sculpting_group.setCheckable(True)
        self.sculpting_group.setChecked(True)  # Default: enabled
        self.sculpting_group.setToolTip("Enable/disable PyMOL sculpting optimization.")
        sculpting_layout = QtWidgets.QGridLayout(self.sculpting_group)
        layout.addWidget(self.sculpting_group, 4, 0, 1, 2)

        # Sculpting Radius
        self.sculpting_radius_label = QtWidgets.QLabel("Sculpting Radius (Å):")
        self.sculpting_radius_spinbox = QtWidgets.QDoubleSpinBox()
        self.sculpting_radius_spinbox.setToolTip("Radius in Angstroms for the sculpting sphere around clashing residues.")
        self.sculpting_radius_spinbox.setRange(1.0, 10.0)
        self.sculpting_radius_spinbox.setSingleStep(0.5)
        self.sculpting_radius_spinbox.setValue(6.0)
        sculpting_layout.addWidget(self.sculpting_radius_label, 0, 0)
        sculpting_layout.addWidget(self.sculpting_radius_spinbox, 0, 1)

        # Sculpting Cycles
        self.sculpting_cycles_label = QtWidgets.QLabel("Sculpting Cycles:")
        self.sculpting_cycles_spinbox = QtWidgets.QSpinBox()
        self.sculpting_cycles_spinbox.setToolTip("Number of sculpting cycles to perform.")
        self.sculpting_cycles_spinbox.setRange(1, 1000)
        self.sculpting_cycles_spinbox.setValue(100)
        sculpting_layout.addWidget(self.sculpting_cycles_label, 1, 0)
        sculpting_layout.addWidget(self.sculpting_cycles_spinbox, 1, 1)

        # --- Open Babel Group ---
        self.ob_group = QtWidgets.QGroupBox("Optional: Open Babel Minimization")
        self.ob_group.setToolTip("Enable/disable Open Babel minimization.")
        self.ob_group.setCheckable(True)
        self.ob_group.setChecked(True)
        ob_layout = QtWidgets.QGridLayout(self.ob_group)
        layout.addWidget(self.ob_group, 5, 0, 1, 2)

        if not has_openbabel:
            self.ob_group.setTitle("Open Babel Minimization (not installed)")
            self.ob_group.setEnabled(False)

        # Optimization Mode
        self.ob_mode_label = QtWidgets.QLabel("Optimization Mode:")
        self.ob_mode_label.setToolTip(
            "Auto: Torsion-only for proteins/nucleic acids, Full for ligands.\n"
            "Torsion Only: Preserves bond lengths/angles, good for macromolecules.\n"
            "Full (CG/SD): Optimizes all coordinates, good for small molecules."
        )
        self.ob_mode_combo = QtWidgets.QComboBox()
        self.ob_mode_combo.setToolTip(
            "Auto: Torsion-only for proteins/nucleic acids, Full for ligands.\n"
            "Torsion Only: Preserves bond lengths/angles, good for macromolecules.\n"
            "Full (CG/SD): Optimizes all coordinates, good for small molecules."
        )
        self.ob_mode_combo.addItems(["Auto", "Torsion Only (Monte Carlo)", "Full (CG/SD)"])
        ob_layout.addWidget(self.ob_mode_label, 0, 0)
        ob_layout.addWidget(self.ob_mode_combo, 0, 1)

        # --- NEW: MC Temperature ---
        self.mc_temp_label = QtWidgets.QLabel("MC Temperature (K):")
        self.mc_temp_spinbox = QtWidgets.QDoubleSpinBox()
        self.mc_temp_spinbox.setToolTip("Initial temperature for Monte Carlo simulated annealing.")
        self.mc_temp_spinbox.setRange(50.0, 5000.0)
        self.mc_temp_spinbox.setSingleStep(300)
        self.mc_temp_spinbox.setValue(100.0)
        ob_layout.addWidget(self.mc_temp_label, 1, 0)
        ob_layout.addWidget(self.mc_temp_spinbox, 1, 1)

        # --- NEW: MC Cycles ---
        self.mc_cycles_label = QtWidgets.QLabel("MC Cycles:")
        self.mc_cycles_spinbox = QtWidgets.QSpinBox()
        self.mc_cycles_spinbox.setToolTip("Number of cycles for Monte Carlo optimization.")
        self.mc_cycles_spinbox.setRange(50, 5000)
        self.mc_cycles_spinbox.setSingleStep(100)
        self.mc_cycles_spinbox.setValue(500)
        ob_layout.addWidget(self.mc_cycles_label, 2, 0)
        ob_layout.addWidget(self.mc_cycles_spinbox, 2, 1)
        
        # Force Field
        self.ob_ff_label = QtWidgets.QLabel("Force Field:")
        self.ob_ff_combo = QtWidgets.QComboBox()
        self.ob_ff_combo.setToolTip("Force field to use for Open Babel minimization.")
        # In the order the protocol falls back through, so the list reads as what it
        # is: what will be tried, in which order, if the one chosen cannot take the
        # structure. GAFF is last because Open Babel prints a typing warning for it on
        # a protein, and the van der Waals term it then computes -- the term a clash
        # refinement rests on -- comes out far below the others on the same atoms.
        self.ob_ff_combo.addItems(['MMFF94s', 'MMFF94', 'UFF', 'Ghemical', 'GAFF'])
        self.ob_ff_combo.setCurrentText("MMFF94s")
        ob_layout.addWidget(self.ob_ff_label, 3, 0)
        ob_layout.addWidget(self.ob_ff_combo, 3, 1)

        # Method (for Full minimization)
        self.ob_method_label = QtWidgets.QLabel("Full Minimization Method:")
        self.ob_method_combo = QtWidgets.QComboBox()
        self.ob_method_combo.setToolTip("Minimization method to use with Open Babel (for Full mode).")
        self.ob_method_combo.addItems(['Conjugate Gradients', 'Steepest Descent'])
        ob_layout.addWidget(self.ob_method_label, 4, 0)
        ob_layout.addWidget(self.ob_method_combo, 4, 1)

        # NSteps (for Full minimization)
        self.ob_nsteps_label = QtWidgets.QLabel("Full Minimization Steps:")
        self.ob_nsteps_spinbox = QtWidgets.QSpinBox()
        self.ob_nsteps_spinbox.setToolTip("Maximum number of minimization steps (for Full mode)." )
        self.ob_nsteps_spinbox.setRange(1, 10000)
        self.ob_nsteps_spinbox.setValue(500)
        ob_layout.addWidget(self.ob_nsteps_label, 5, 0)
        ob_layout.addWidget(self.ob_nsteps_spinbox, 5, 1)

        # Convergence
        self.ob_conv_label = QtWidgets.QLabel("Convergence:")
        self.ob_conv_spinbox = QtWidgets.QDoubleSpinBox()
        self.ob_conv_spinbox.setToolTip("Convergence criterion for minimization.")
        self.ob_conv_spinbox.setDecimals(5)
        self.ob_conv_spinbox.setRange(0.00001, 0.1)
        self.ob_conv_spinbox.setSingleStep(0.001)
        self.ob_conv_spinbox.setValue(0.01)
        ob_layout.addWidget(self.ob_conv_label, 6, 0)
        ob_layout.addWidget(self.ob_conv_spinbox, 6, 1)
        
        # Cutoff
        self.ob_cutoff_checkbox = QtWidgets.QCheckBox("Enable Cutoff")
        self.ob_cutoff_checkbox.setToolTip("Enable/disable cutoff for non-bonded interactions.")
        self.ob_cutoff_checkbox.setChecked(True)
        ob_layout.addWidget(self.ob_cutoff_checkbox, 7, 0, 1, 2)

        # VDW Cutoff
        self.ob_vdw_cutoff_label = QtWidgets.QLabel("VDW Cutoff (Å):")
        self.ob_vdw_cutoff_spinbox = QtWidgets.QDoubleSpinBox()
        self.ob_vdw_cutoff_spinbox.setToolTip("Cutoff in Angstroms for Van der Waals interactions.")
        self.ob_vdw_cutoff_spinbox.setRange(1.0, 20.0)
        self.ob_vdw_cutoff_spinbox.setSingleStep(0.5)
        self.ob_vdw_cutoff_spinbox.setValue(6.0)
        ob_layout.addWidget(self.ob_vdw_cutoff_label, 8, 0)
        ob_layout.addWidget(self.ob_vdw_cutoff_spinbox, 8, 1)

        # Elec Cutoff
        self.ob_elec_cutoff_label = QtWidgets.QLabel("Electrostatic Cutoff (Å):")
        self.ob_elec_cutoff_spinbox = QtWidgets.QDoubleSpinBox()
        self.ob_elec_cutoff_spinbox.setToolTip("Cutoff in Angstroms for electrostatic interactions.")
        self.ob_elec_cutoff_spinbox.setRange(1.0, 20.0)
        self.ob_elec_cutoff_spinbox.setSingleStep(0.5)
        self.ob_elec_cutoff_spinbox.setValue(8.0)
        ob_layout.addWidget(self.ob_elec_cutoff_label, 9, 0)
        ob_layout.addWidget(self.ob_elec_cutoff_spinbox, 9, 1)

        # Mobile Region Cutoff
        self.mobile_cutoff_label = QtWidgets.QLabel("Mobile Region Radius (Å):")
        self.mobile_cutoff_spinbox = QtWidgets.QDoubleSpinBox()
        self.mobile_cutoff_spinbox.setToolTip(
            "Radius in Angstroms around clash sites where atoms can move during minimization.\n"
            "Atoms outside this radius will be kept fixed to speed up the calculation.\n"
            "Smaller values = faster but more constrained, larger values = slower but more flexible."
        )
        self.mobile_cutoff_spinbox.setRange(3.0, 15.0)
        self.mobile_cutoff_spinbox.setSingleStep(0.5)
        self.mobile_cutoff_spinbox.setValue(4.0)
        ob_layout.addWidget(self.mobile_cutoff_label, 10, 0)
        ob_layout.addWidget(self.mobile_cutoff_spinbox, 10, 1)
        
        # --- Button Box ---
        self.button_box = QtWidgets.QDialogButtonBox()
        self.refine_button = self.button_box.addButton("Refine Clashes", QtWidgets.QDialogButtonBox.AcceptRole)
        self.show_button = self.button_box.addButton("Show Clashes Only", QtWidgets.QDialogButtonBox.ActionRole)
        self.cancel_button = self.button_box.addButton(QtWidgets.QDialogButtonBox.Cancel)
        layout.addWidget(self.button_box, 6, 0, 1, 2)

        self.refine_button.clicked.connect(self.accept_and_refine)
        self.show_button.clicked.connect(self.accept_and_show)
        self.cancel_button.clicked.connect(self.reject)

        # --- Initial State ---
        self.toggle_advanced_settings(False) # Hide advanced settings initially
        self.update_minimization_controls("Auto-detect") # Set initial control state
        self.adjustSize() # Adjust window size

    def toggle_advanced_settings(self, checked):
        self.params_group.setVisible(checked)
        self.refine_group.setVisible(checked)
        self.sculpting_group.setVisible(checked)
        self.ob_group.setVisible(checked)
        if checked:
            self.advanced_button.setText("Hide Advanced Settings")
        else:
            self.advanced_button.setText("Show Advanced Settings")
        self.adjustSize()

    def update_minimization_controls(self, minimization_type):
        """Enable/disable controls based on minimization type."""
        is_full = minimization_type == "Full"
        is_torsion = minimization_type == "Torsion-only"
        
        # Method is only relevant for full minimization
        self.ob_method_label.setEnabled(is_full or minimization_type == "Auto-detect")
        self.ob_method_combo.setEnabled(is_full or minimization_type == "Auto-detect")
        
        # Convergence is only relevant for full minimization
        self.ob_conv_label.setEnabled(is_full or minimization_type == "Auto-detect")
        self.ob_conv_spinbox.setEnabled(is_full or minimization_type == "Auto-detect")
        
        # Cutoff controls are only for full minimization
        self.ob_cutoff_checkbox.setEnabled(is_full or minimization_type == "Auto-detect")
        self.ob_vdw_cutoff_label.setEnabled(is_full or minimization_type == "Auto-detect")
        self.ob_vdw_cutoff_spinbox.setEnabled(is_full or minimization_type == "Auto-detect")
        self.ob_elec_cutoff_label.setEnabled(is_full or minimization_type == "Auto-detect")
        self.ob_elec_cutoff_spinbox.setEnabled(is_full or minimization_type == "Auto-detect")
        
        # Update tooltip based on type
        if is_torsion:
            self.ob_nsteps_spinbox.setToolTip("Number of Monte Carlo moves for torsion optimization.")
        elif is_full:
            self.ob_nsteps_spinbox.setToolTip("Maximum number of minimization steps.")
        else:  # Auto-detect
            self.ob_nsteps_spinbox.setToolTip("Number of optimization steps (auto-adapts based on molecule type).")

    def accept_and_refine(self):
        # Check that at least one optimization method is selected
        sculpting_enabled = self.sculpting_group.isChecked()
        minimization_enabled = self.ob_group.isEnabled() and self.ob_group.isChecked()
        
        if not sculpting_enabled and not minimization_enabled:
            QtWidgets.QMessageBox.warning(
                self, 
                "No Optimization Method Selected", 
                "Please select at least one optimization method:\n"
                "- Enable Sculpting Parameters, or\n"
                "- Enable Open Babel Minimization\n\n"
                "Both methods can be used together for better results."
            )
            return
        
        self.action = 'refine'
        self.accept()

    def accept_and_show(self):
        self.action = 'show'
        self.accept()

    def toggle_aggressive_mode(self, checked):
        is_enabled = not checked
        self.min_overlap_spinbox.setEnabled(is_enabled)
        self.max_display_spinbox.setEnabled(is_enabled)
        self.ignore_h_checkbox.setEnabled(is_enabled)
        if checked:
            self.min_overlap_spinbox.setValue(0.6)
            self.max_display_spinbox.setValue(500)
            self.ignore_h_checkbox.setChecked(False)
        else:
            self.min_overlap_spinbox.setValue(0.30)
            self.max_display_spinbox.setValue(200)
            self.ignore_h_checkbox.setChecked(True)

    def update_selection_options(self, selected_elements=None):
        self.selection_combobox.clear()
        elements = selected_elements or []

        if len(elements) > 1:
            # --- COMPLEX MODE ---
            element_selectors = [e.get_pymol_selector() for e in elements]
            complex_sele_str = " or ".join(element_selectors)
            cmd.select("clash_refine_complex", complex_sele_str)

            self.selection_combobox.addItem("Complex")
            self.selection_combobox.addItems(element_selectors)

        elif len(elements) == 1:
            # --- SINGLE SELECTION MODE ---
            selector = elements[0].get_pymol_selector()
            self.selection_combobox.addItem(selector)
            all_selections = [name for name in cmd.get_names("selections") if cmd.get_type(name) == 'selection']
            relevant_selections = [s for s in all_selections if cmd.count_atoms(f"({s}) and ({selector})") > 0 and s != selector]
            self.selection_combobox.addItems(relevant_selections)

        else:
            # --- DEFAULT/NO SELECTION MODE ---
            selections = [name for name in cmd.get_names("selections") if cmd.get_type(name) == 'selection']
            if not selections:
                self.selection_combobox.addItem("Full Molecule")
            else:
                self.selection_combobox.addItems(selections)

    def get_values(self):
        vals = {
            "selection": self.selection_combobox.currentText(),
            "threshold": self.threshold_spinbox.value(),
            "min_overlap": self.min_overlap_spinbox.value(),
            "ignore_h": self.ignore_h_checkbox.isChecked(),
            "max_display": self.max_display_spinbox.value(),
            "aggressive": self.aggressive_checkbox.isChecked(),
            "action": self.action,
            "keep_original": self.keep_original_checkbox.isChecked(),
            "run_sculpting": self.sculpting_group.isChecked(),
            "sculpting_radius": self.sculpting_radius_spinbox.value(),
            "sculpting_cycles": self.sculpting_cycles_spinbox.value(),
            
        }
        if self.ob_group.isEnabled():
            vals["run_ob_minimize"] = self.ob_group.isChecked()
            vals["ob_optimization_mode"] = self.ob_mode_combo.currentText()
            vals["mc_temp"] = self.mc_temp_spinbox.value()
            vals["mc_cycles"] = self.mc_cycles_spinbox.value()
            vals["ob_forcefield"] = self.ob_ff_combo.currentText()
            vals["ob_method"] = self.ob_method_combo.currentText()
            vals["ob_nsteps"] = self.ob_nsteps_spinbox.value()
            vals["ob_conv"] = self.ob_conv_spinbox.value()
            vals["ob_cutoff"] = self.ob_cutoff_checkbox.isChecked()
            vals["ob_vdw_cutoff"] = self.ob_vdw_cutoff_spinbox.value()
            vals["ob_elec_cutoff"] = self.ob_elec_cutoff_spinbox.value()
            vals["mobile_cutoff_radius"] = self.mobile_cutoff_spinbox.value()
        return vals
