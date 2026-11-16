# Copyright 2025 by me. All rights reserved.
# This code is part of the PyMod package and governed by its license.
#
# This file is adapted from clash_and_refine.py for single-molecule optimization.

import math
import re
import pymol
from pymol import cmd

# Use temporary GUI for now - replace with actual import when available
try:
    from pymod_lib.pymod_gui.molecule_optimizer_gui import MoleculeOptimizerWindow
except ImportError:
    try:
        from tmp_rovodev_molecule_optimizer_gui import MoleculeOptimizerWindow
    except ImportError:
        MoleculeOptimizerWindow = None

# Check for Open Babel and Qt
try:
    from openbabel import openbabel as ob
    has_openbabel = True
except ImportError:
    has_openbabel = False

try:
    from pymol.Qt import QtWidgets
    has_qt = True
except ImportError:
    has_qt = False

VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
    'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
}

# The window the command runs in front of: what it is doing, how far it has got, and
# everything it has to say. The messages used to go to the PyMOL console only, where a run
# looked like nothing happening at all, and the bar only existed for the minimisation --
# which is the last stage, and never reached when an earlier one failed.
if has_qt:
    class OptimizationProgress(QtWidgets.QDialog):

        def __init__(self, title="Optimize Molecule", parent=None):
            super().__init__(parent)
            self.setWindowTitle(title)
            self.setModal(False)
            self.resize(620, 380)
            self.cancelled = False

            layout = QtWidgets.QVBoxLayout(self)
            self.stage_label = QtWidgets.QLabel("Starting...", self)
            layout.addWidget(self.stage_label)

            self.progress = QtWidgets.QProgressBar(self)
            self.busy()
            layout.addWidget(self.progress)

            self.log_box = QtWidgets.QPlainTextEdit(self)
            self.log_box.setReadOnly(True)
            # Wrapped, not scrolled sideways: a message about a force field is long and
            # the end of it is the part which says why.
            self.log_box.setLineWrapMode(QtWidgets.QPlainTextEdit.WidgetWidth)
            font = self.log_box.font()
            font.setFamily("Courier New")
            self.log_box.setFont(font)
            layout.addWidget(self.log_box)

            buttons = QtWidgets.QHBoxLayout()
            buttons.addStretch()
            self.close_button = QtWidgets.QPushButton("Close", self)
            self.close_button.setEnabled(False)
            self.close_button.clicked.connect(self.close)
            buttons.addWidget(self.close_button)
            layout.addLayout(buttons)

        def busy(self):
            """A bar which moves without knowing how much is left, for a stage which
            cannot count its steps."""
            self.progress.setRange(0, 0)

        def counting(self, total):
            self.progress.setRange(0, max(1, total))
            self.progress.setValue(0)

        def stage(self, text):
            self.stage_label.setText(text)
            self.busy()
            self._breathe()

        def step(self, done, text=None):
            if self.progress.maximum() > 0:
                self.progress.setValue(min(done, self.progress.maximum()))
            if text:
                self.stage_label.setText(text)
            self._breathe()

        def log(self, text):
            self.log_box.appendPlainText(text)
            bar = self.log_box.verticalScrollBar()
            bar.setValue(bar.maximum())
            self._breathe()

        def finished(self, text):
            self.progress.setRange(0, 1)
            self.progress.setValue(1)
            self.stage_label.setText(text)
            self.close_button.setEnabled(True)
            self.close_button.setDefault(True)
            self._breathe()

        def _breathe(self):
            """Let Qt paint: the work runs in this thread, so nothing is drawn unless the
            event loop is given a turn."""
            QtWidgets.QApplication.processEvents()


class MoleculeOptimizer:
    def __init__(self, pymod_instance, selected_elements=None):
        self.pymod = pymod_instance
        self.gui_elements = selected_elements or []
        self.progress_window = None
        # Set when a stage fails, so the window does not end on "Done".
        self.failure_note = None

    # ------------------------------------------------------------------
    # What the user is told while the command runs.
    # ------------------------------------------------------------------

    def _log(self, *parts):
        """Say it in the window of the run, and in the console as before."""
        text = " ".join(str(part) for part in parts)
        print(text)
        if self.progress_window is not None:
            try:
                self.progress_window.log(text)
            except Exception:
                # A window which has been closed must not take the run down with it.
                self.progress_window = None

    def _stage(self, text):
        """A new stage: named, with a bar which moves while it lasts."""
        self._log("\n== %s" % text)
        if self.progress_window is not None:
            try:
                self.progress_window.stage(text)
            except Exception:
                self.progress_window = None

    def _open_progress_window(self, title="Optimize Molecule"):
        if not has_qt:
            return None
        try:
            parent = getattr(self.pymod, "main_window", None)
            window = OptimizationProgress(title=title, parent=parent)
            window.show()
            window.raise_()
            QtWidgets.QApplication.processEvents()
            self.progress_window = window
            return window
        except Exception as e:
            print("Could not open the progress window: %s" % e)
            self.progress_window = None
            return None

    def _close_progress_window(self, text):
        if self.progress_window is None:
            return
        try:
            self.progress_window.finished(text)
        except Exception:
            pass

    def structure_elements(self):
        """
        The elements of the session which have a structure in PyMOL.

        get_pymod_elements_list() returns everything in the main window, clusters and
        alignments included, and those have no PyMOL object: asking one for its selector
        raises AttributeError. A project with an alignment in it -- after an inverse folding
        run, say -- used to bring this command down before it opened.
        """

        elements = []
        for element in self.pymod.get_pymod_elements_list():
            try:
                if element.has_structure():
                    elements.append(element)
            except Exception:
                continue
        return elements

    def run(self):
        if not MoleculeOptimizerWindow:
            self.pymod.main_window.show_error_message("GUI Error", "Molecule Optimizer GUI not available.")
            return
            
        dialog = MoleculeOptimizerWindow(self.pymod.main_window, self.gui_elements, has_openbabel)
        if not dialog.exec_():
            return

        values = dialog.get_values()
        mobile_selection = values["mobile_object"]

        # From here on the run has something to show for itself.
        self._open_progress_window("Optimize Molecule — %s" % mobile_selection)
        self._stage("Reading the selection")
        
        # Determine if we're working with a specific selection or full object
        is_selection = mobile_selection in cmd.get_names('selections')
        
        if is_selection:
            self._log(f"Working with PyMOL selection: {mobile_selection}")

            # Which object the selection belongs to. Asked of PyMOL, which knows: every
            # atom carries the name of its object in 'model'. It used to be worked out by
            # looking for an atom of the same index in each object in turn, and an index is
            # numbered within its own object, so the first object long enough to have that
            # index answered instead -- for a selection on the second of two objects, the
            # first one was returned, and it was then the first one which got replaced.
            parent_objects = set()
            cmd.iterate(mobile_selection, "parent_objects.add(model)",
                        space={"parent_objects": parent_objects})

            # The window says which object the selection was listed under; when it does,
            # that is the answer, and this only confirms it.
            named_parent = values.get("mobile_parent")
            if named_parent and named_parent in parent_objects:
                parent_objects = {named_parent}

            if not parent_objects:
                self.pymod.main_window.show_error_message("Selection Error", f"Could not determine parent objects for selection '{mobile_selection}'.")
                self._close_progress_window("Stopped: the selection has no object.")
                return

            if len(parent_objects) > 1:
                self.pymod.main_window.show_error_message(
                    "Selection Error",
                    "The selection '%s' spans more than one object (%s).\n\nThis tool"
                    " optimizes one molecule: please select a region of a single object,"
                    " or use Optimize Interface."
                    % (mobile_selection, ", ".join(sorted(parent_objects))))
                self._close_progress_window("Stopped: the selection spans more than one"
                                            " object.")
                return

            # Nothing else in the session takes part. This command optimizes one molecule:
            # the residues which were selected move, the rest of their own object holds
            # them in place, and every other object is left out of it. Taking them all as
            # an environment is what it used to do, and it made the command refuse to run
            # on a model: a model sits where its template sits, so it clashes with it by
            # construction, and with anything else already in that frame.
            static_objects = []
            mobile_object_name = sorted(parent_objects)[0]
            self._log(f"  the selection belongs to '{mobile_object_name}'")
            
        else:
            self._log(f"Working with full object: {mobile_selection}")
            mobile_object_name = mobile_selection
            # As above: this object alone.
            static_objects = []

        other_objects = [name for name in cmd.get_names("objects")
                         if name != mobile_object_name]
        if other_objects:
            self._log("  %s other object(s) are open and take no part in this: %s"
                  % (len(other_objects), ", ".join(other_objects[:6])
                     + (", ..." if len(other_objects) > 6 else "")))

        # Create trajectory object for minimization frames
        trajectory_object = f"{mobile_object_name}_minimization_trajectory"
        self._cleanup_objects([trajectory_object])
        
        if values["keep_original"]:
            original_copy_name = f"{mobile_object_name}_original"
            self._cleanup_objects([original_copy_name])
            if is_selection:
                # For selections, save the entire parent object
                cmd.create(original_copy_name, mobile_object_name)
            else:
                cmd.create(original_copy_name, mobile_selection)

        self._stage("Adding hydrogens")
        if is_selection:
            cmd.h_add(mobile_object_name)
        else:
            cmd.h_add(mobile_selection)

        # Optional: Find and resolve internal clashes
        if values.get("resolve_clashes", False):
            self._stage("Looking for clashes inside the molecule")
            clash_selection = mobile_selection if not is_selection else mobile_object_name
            internal_clashes = self._find_clashes(values, selection_name=clash_selection)

            if internal_clashes:
                self._log(f"Found {len(internal_clashes)} internal clashes. Starting refinement...")
                if self.progress_window is not None:
                    self.progress_window.stage("Relaxing %d clashing residues"
                                               % len(internal_clashes))
                self._show_bumps(values, internal_clashes, "internal_clashes_found")
                if is_selection:
                    self._refine_clashes(mobile_object_name, internal_clashes, mobile_selection, values)
                else:
                    self._refine_clashes(mobile_selection, internal_clashes, mobile_selection, values)
            else:
                self._log("No significant internal clashes found.")

        if values["run_ob_minimize"]:
            self._stage("Minimising")
            self._run_minimization(mobile_selection, static_objects, values, is_selection, mobile_object_name)
        
        self._cleanup_objects(["internal_clashes_found", "bumps_*"])
        self._log("Optimization complete.")
        self._close_progress_window(
            self.failure_note or "Done. The messages above are the whole run.")

    def _refine_clashes(self, target_object, clashes, mobile_object_name, values):
        from pymod_lib import pymod_vars as pmdt
        self._log("Refining clashes...")

        # 1. Filter for clashing residues on the mobile object only
        clashing_residues_info = set()
        resn_map = {}
        cmd.iterate(mobile_object_name, "resn_map[(chain, str(resi))] = resn", space={'resn_map': resn_map})

        for clash in clashes:
            # For internal clashes, both atoms are from the mobile object
            atom_mobile1 = clash['a1']
            resn1 = resn_map.get((atom_mobile1.chain, str(atom_mobile1.resi)))
            clashing_residues_info.add((atom_mobile1.chain, atom_mobile1.resi, resn1, mobile_object_name))

            atom_mobile2 = clash['a2']
            resn2 = resn_map.get((atom_mobile2.chain, str(atom_mobile2.resi)))
            clashing_residues_info.add((atom_mobile2.chain, atom_mobile2.resi, resn2, mobile_object_name))

        if not clashing_residues_info:
            self._log("No clashes involving the mobile object found to refine.")
            return

        self._log(f"Refining {len(clashing_residues_info)} clashing residues on mobile object: {mobile_object_name}")

        # 2. Mutagenesis Wizard to relax sidechains
        try:
            cmd.wizard("mutagenesis")
            wiz = cmd.get_wizard()
            if self.progress_window is not None:
                self.progress_window.counting(len(clashing_residues_info))
            for done, (chain, resi, resn, _) in enumerate(clashing_residues_info, start=1):
                sele = f"/{target_object}//{chain}/{resi}"
                if self.progress_window is not None:
                    self.progress_window.step(
                        done, "Relaxing side chains — residue %d of %d (%s %s)"
                        % (done, len(clashing_residues_info), resn or "?", resi))
                if resn in ["ALA", "GLY"]: continue
                relax_target_res = "ALA"
                wiz.set_mode(relax_target_res)
                wiz.do_select(sele)
                wiz.apply()
                cmd.refresh()
                wiz.set_mode(resn)
                wiz.do_select(sele)
                wiz.apply()
            cmd.set_wizard()
            self._log("Residue relaxation with wizard complete.")
        except Exception as e:
            self._log(f"Warning: Mutagenesis wizard failed: {e}")
            cmd.set_wizard()

        # 3. Sculpting on the refined residues
        self._log("Performing local sculpting...")
        core_sele_name = "sculpt_core"
        sculpt_area_name = "sculpt_area"
        core_sele_parts = [f"(chain {chain} and resi {resi})" for chain, resi, _, _ in clashing_residues_info]
        core_sele_string = f"({target_object}) and ({' or '.join(core_sele_parts)})"
        if cmd.select(core_sele_name, core_sele_string) > 0:
            sculpt_sele_string = f"(byres ({core_sele_name} around {values['sculpting_radius']})) and ({target_object})"
            if cmd.select(sculpt_area_name, sculpt_sele_string) > 0:
                try:
                    cmd.sculpt_activate(target_object)
                    cmd.sculpt_iterate(target_object, cycles=values['sculpting_cycles'])
                finally:
                    cmd.sculpt_deactivate(target_object)
            self._cleanup_objects([core_sele_name, sculpt_area_name])
        self._log("Sculpting finished.")

    def _run_minimization(self, mobile_selection, static_object_names, values, is_selection=False, mobile_object_name=None):
        """Usa la logica di clash_and_refine.py che funziona bene per la minimizzazione"""
        if not has_openbabel:
            self._log("OpenBabel not found. Skipping minimization")
            return

        # Determina quali residui devono essere mobili
        if is_selection:
            self._log(f"🔒 SELECTION MODE: Solo '{mobile_selection}' sarà mobile")
            # Per le selezioni, trova i residui della selezione mobile
            mobile_residues = self._get_mobile_residue_specs_from_selection(mobile_selection, mobile_object_name)
            
            # Crea il target object combinato come in clash_and_refine
            trajectory_object = f"{mobile_object_name}_minimization_trajectory"
            self._cleanup_objects([trajectory_object])
            all_objects_sele = " or ".join([mobile_object_name] + static_object_names)
            cmd.create(trajectory_object, all_objects_sele)
            
        else:
            self._log(f"🔒 OBJECT MODE: Solo '{mobile_selection}' sarà mobile")
            # Per oggetti completi, trova tutti i residui dell'oggetto mobile
            mobile_residues = self._get_mobile_residue_specs_from_object(mobile_selection)
            
            trajectory_object = f"{mobile_selection}_minimization_trajectory"  
            self._cleanup_objects([trajectory_object])
            all_objects_sele = " or ".join([mobile_selection] + static_object_names)
            cmd.create(trajectory_object, all_objects_sele)

        # Le specifiche dei residui mobili sono già state calcolate
        self._log(f"Mobile residues identified: {len(mobile_residues)}")
            
        # Cutoffs make the energy cheaper to compute and are on by default in the window.
        # They used to be forced on here, which left the checkbox doing nothing; turning
        # them off is slower and steadier, and is now the user's choice again.
        use_cutoff = values.get("ob_cutoff", True)
        vdw_cutoff = values.get("ob_vdw_cutoff", 6.0)
        elec_cutoff = values.get("ob_elec_cutoff", 8.0)

        self._log(f"Energy calculation settings:")
        self._log(f"  Use cutoff: {use_cutoff} (from the window)")
        self._log(f"  → VDW Cutoff: {vdw_cutoff} Å")
        self._log(f"  → Electrostatic Cutoff: {elec_cutoff} Å")
        
        # CORREZIONE: Salva frame iniziale senza idrogeni per consistenza
        self._log("Saving initial frame without hydrogens...")
        initial_pdb = cmd.get_pdbstr(f"({trajectory_object}) and not hydrogens")
        if initial_pdb:
            cmd.read_pdbstr(initial_pdb, trajectory_object, state=1, quiet=1)
        
        # Chiama il metodo di minimizzazione di clash_and_refine
        result = self._run_openbabel_minimize_adapted(
            selection=trajectory_object,
            forcefield=values["ob_forcefield"],
            method="Steepest Descent",
            nsteps=values["ob_nsteps"],
            conv=0.01,
            core_residue_specs=mobile_residues,
            cutoff=1 if use_cutoff else 0,
            cut_vdw=vdw_cutoff,
            cut_elec=elec_cutoff,
            trajectory_object=trajectory_object
        )
        
        if result == 1:  # Success
            self._log("Minimization successful, proceeding to finalization.")
            # Apply trajectory coloring
            self._apply_trajectory_coloring(trajectory_object, mobile_selection, is_selection)
            
            # Finalize by replacing the original object with the refined one
            original_object = mobile_object_name if is_selection else mobile_selection
            self._finalize_and_replace(trajectory_object, original_object, values)
            
            self._log(f"✅ Trajectory saved with {cmd.count_states(trajectory_object)} states in '{trajectory_object}'")
        else:
            self._log("❌ Minimization failed!")
            self.failure_note = ("The minimisation did not run. What was done before it"
                                 " stands; read the messages above.")

    def _get_mobile_residue_specs_from_selection(self, selection, parent_object):
        """Ottieni le specifiche dei residui per una selezione PyMOL"""
        residue_specs = []
        model = cmd.get_model(selection)
        seen_residues = set()
        
        for atom in model.atom:
            res_key = (atom.chain, atom.resi)
            if res_key not in seen_residues:
                residue_specs.append((atom.chain, atom.resi, atom.resn))
                seen_residues.add(res_key)
                
        self._log(f"Found {len(residue_specs)} mobile residues in selection '{selection}'")
        return residue_specs
    
    def _get_mobile_residue_specs_from_object(self, object_name):
        """Ottieni le specifiche dei residui per un oggetto completo"""
        residue_specs = []
        model = cmd.get_model(object_name)
        seen_residues = set()
        
        for atom in model.atom:
            res_key = (atom.chain, atom.resi)
            if res_key not in seen_residues:
                residue_specs.append((atom.chain, atom.resi, atom.resn))
                seen_residues.add(res_key)
                
        self._log(f"Found {len(residue_specs)} mobile residues in object '{object_name}'")
        return residue_specs
    
    def _get_mobile_atom_indices(self, trajectory_object, mobile_residues):
        """Ottieni gli indici degli atomi mobili basati sui residui mobili"""
        mobile_indices = []
        
        # Crea set per ricerca veloce
        mobile_res_set = {(chain, str(resi)) for chain, resi, _ in mobile_residues}
        
        model = cmd.get_model(trajectory_object)
        for i, atom in enumerate(model.atom):
            if (atom.chain, atom.resi) in mobile_res_set:
                mobile_indices.append(i + 1)  # OpenBabel usa 1-based indexing
                
        return mobile_indices



    def _apply_trajectory_coloring(self, trajectory_object, mobile_selection, is_selection):
        """CORREZIONE BUG 2: Apply coloring to trajectory basato su residui mobili corretti"""
        try:
            self._log("Applying trajectory coloring...")
            
            # Color everything gray first (fixed regions)
            cmd.color("gray", trajectory_object)
            
            if is_selection:
                # CORREZIONE BUG 2: Per le selezioni, identificare correttamente i residui mobili
                self._log(f"SELECTION MODE: Coloring mobile selection '{mobile_selection}' in green")
                
                # Trova i residui della selezione mobile
                mobile_residues = set()
                model = cmd.get_model(mobile_selection)
                for atom in model.atom:
                    mobile_residues.add((atom.chain, atom.resi))
                
                if mobile_residues:
                    # Crea selezione per residui mobili nella traiettoria
                    residue_selectors = []
                    for chain, resi in mobile_residues:
                        residue_selectors.append(f"(chain {chain} and resi {resi})")
                    
                    mobile_traj_sele = f"({trajectory_object}) and ({' or '.join(residue_selectors)})"
                    mobile_count = cmd.count_atoms(mobile_traj_sele)
                    
                    if mobile_count > 0:
                        # Apply green coloring using PyMOL color by atom
                        cmd.util.cbag(mobile_traj_sele)
                        cmd.color("green", f"({mobile_traj_sele}) and element C")
                        self._log(f"✅ Colored {len(mobile_residues)} mobile residues ({mobile_count} atoms) in green")
                    else:
                        self._log(f"⚠️ No atoms found for mobile residues in trajectory")
                else:
                    self._log(f"⚠️ No mobile residues identified from selection")
                    
            else:
                # CORREZIONE BUG 2: Per oggetti completi, colorare basandosi su catene
                self._log(f"OBJECT MODE: Coloring mobile object '{mobile_selection}' in green")
                
                try:
                    # Ottieni le catene dall'oggetto mobile originale
                    mobile_chains = []
                    cmd.iterate(mobile_selection, "mobile_chains.append(chain)", space={'mobile_chains': mobile_chains})
                    mobile_chains = list(set(mobile_chains))
                    
                    if mobile_chains:
                        chain_selector = " or ".join([f"chain {chain}" for chain in mobile_chains])
                        mobile_traj_sele = f"({trajectory_object}) and ({chain_selector})"
                        mobile_count = cmd.count_atoms(mobile_traj_sele)
                        
                        if mobile_count > 0:
                            # Apply green coloring using PyMOL color by atom
                            cmd.util.cbag(mobile_traj_sele)
                            cmd.color("green", f"({mobile_traj_sele}) and element C")
                            self._log(f"✅ Colored mobile chains {mobile_chains} ({mobile_count} atoms) in green")
                        else:
                            self._log(f"⚠️ No atoms found for mobile chains in trajectory")
                    else:
                        self._log(f"⚠️ Could not determine chains from mobile object")
                        
                except Exception as e:
                    self._log(f"⚠️ Chain-based coloring failed: {e}")
                    # Fallback: prova a usare nome dell'oggetto
                    try:
                        cmd.util.cbag(f"({trajectory_object})")
                        cmd.color("green", f"({trajectory_object}) and element C")
                        self._log("Used fallback coloring for entire trajectory")
                    except:
                        self._log("Warning: Could not apply any mobile coloring")
                        
            # Set representation to show the trajectory nicely
            cmd.show("cartoon", trajectory_object)
            cmd.show("sticks", f"({trajectory_object}) and sidechain")
            
            self._log(f"✅ Trajectory coloring complete: Fixed regions = gray, Mobile regions = green")
            
        except Exception as e:
            self._log(f"❌ Warning: Trajectory coloring failed: {e}")
            # Fallback to simple coloring
            try:
                cmd.color("gray", trajectory_object)
                cmd.show("cartoon", trajectory_object)
                self._log("Applied fallback gray coloring")
            except:
                pass

    def _update_pymol_coordinates_from_trajectory(self, target_object, trajectory_object, final_state, mobile_selection=None):
        """CORREZIONE BUG 3: Aggiorna le coordinate di PyMOL dal trajectory object (solo heavy atoms)"""
        try:
            self._log(f"Updating coordinates for '{target_object}' from '{trajectory_object}' state {final_state}")
            
            # IMPORTANTE: Aggiorna solo gli heavy atoms per evitare problemi con idrogeni OpenBabel
            heavy_atom_filter = "not hydrogens"
            
            if mobile_selection:
                # Per selezioni, aggiorna solo gli heavy atoms della selezione mobile
                self._log(f"Selection mode: updating only heavy atoms from selection '{mobile_selection}'")
                
                # Filtra solo heavy atoms della selezione mobile
                heavy_mobile_sele = f"({mobile_selection}) and ({heavy_atom_filter})"
                heavy_count = cmd.count_atoms(heavy_mobile_sele)
                
                if heavy_count > 0:
                    self._log(f"Found {heavy_count} heavy atoms to update")
                    
                    # CORREZIONE: Usa PyMOL per mappare coordinate automaticamente
                    try:
                        self._log(f"Creating temporary object from final trajectory state for coordinate mapping")
                        
                        # Crea oggetto temporaneo dal frame finale del trajectory
                        temp_final_obj = "tmp_final_state_coords"
                        cmd.create(temp_final_obj, f"{trajectory_object} and state {final_state}", 
                                  source_state=final_state, target_state=1)
                        
                        temp_count = cmd.count_atoms(temp_final_obj)
                        self._log(f"  Temporary object created: {temp_count} atoms from trajectory state {final_state}")
                        
                        if temp_count > 0:
                            # Usa PyMOL update per copiare coordinate dei residui mobili
                            # Crea selezione per gli stessi residui nel temp object
                            mobile_residues = set()
                            cmd.iterate(heavy_mobile_sele, "mobile_residues.add((chain, resi))", 
                                       space={'mobile_residues': mobile_residues})
                            
                            if mobile_residues:
                                # Costruisci selezione per residui corrispondenti nel temp object
                                residue_selectors = []
                                for chain, resi in mobile_residues:
                                    residue_selectors.append(f"(chain {chain} and resi {resi})")
                                
                                temp_mobile_sele = f"({temp_final_obj}) and ({' or '.join(residue_selectors)}) and not hydrogens"
                                temp_mobile_count = cmd.count_atoms(temp_mobile_sele)
                                
                                self._log(f"  Found {temp_mobile_count} corresponding heavy atoms in temporary object")
                                
                                if temp_mobile_count > 0:
                                    # Usa PyMOL update per copiare coordinate
                                    cmd.update(heavy_mobile_sele, temp_mobile_sele)
                                    self._log(f"✅ Updated {heavy_count} heavy atoms using PyMOL coordinate mapping")
                                else:
                                    raise Exception("No corresponding atoms found in temporary object")
                            else:
                                raise Exception("No mobile residues identified")
                        else:
                            raise Exception("Temporary object is empty")
                        
                        # Cleanup
                        cmd.delete(temp_final_obj)
                        
                    except Exception as e:
                        self._log(f"⚠️ PyMOL coordinate mapping failed: {e}, using fallback method")
                        
                        # Cleanup in case of error
                        if "tmp_final_state_coords" in cmd.get_names("objects"):
                            cmd.delete("tmp_final_state_coords")
                        
                        # Fallback: usa il metodo standard di PyMOL senza allineamento
                        self._log("Using PyMOL update method as fallback...")
                        
                        try:
                            # Usa cmd.update ma specifica gli heavy atoms della selezione mobile
                            cmd.update(heavy_mobile_sele, 
                                      f"({trajectory_object}) and state {final_state} and ({heavy_mobile_sele})",
                                      source_state=final_state, target_state=1)
                            self._log(f"✅ Updated heavy atoms using PyMOL update method")
                            
                        except Exception as e2:
                            self._log(f"❌ PyMOL update method also failed: {e2}")
                            self._log("⚠️ Coordinate update unsuccessful - structure may not be updated")
                else:
                    self._log(f"⚠️ No heavy atoms found in mobile selection")
            else:
                # Per oggetti completi, usa la stessa logica robusta della modalità selezione
                self._log(f"Object mode: updating heavy atoms of entire object '{target_object}'")
                
                # Usa l'intero oggetto come selezione mobile
                heavy_mobile_sele = f"({target_object}) and ({heavy_atom_filter})"
                heavy_count = cmd.count_atoms(heavy_mobile_sele)
                
                if heavy_count > 0:
                    self._log(f"Found {heavy_count} heavy atoms to update in '{target_object}'")
                    
                    # Usa la logica di mapping robusta
                    try:
                        self._log(f"Creating temporary object from final trajectory state for coordinate mapping")
                        
                        temp_final_obj = "tmp_final_state_coords"
                        cmd.create(temp_final_obj, f"{trajectory_object} and state {final_state}", 
                                  source_state=final_state, target_state=1)
                        
                        temp_count = cmd.count_atoms(temp_final_obj)
                        self._log(f"  Temporary object created: {temp_count} atoms from trajectory state {final_state}")
                        
                        if temp_count > 0:
                            # Usa PyMOL update per copiare le coordinate. PyMOL mapperà gli atomi corrispondenti.
                            # Non è necessario un selettore di residui qui perché stiamo aggiornando l'intero oggetto.
                            cmd.update(heavy_mobile_sele, f"({temp_final_obj}) and not hydrogens")
                            self._log(f"✅ Updated {heavy_count} heavy atoms using PyMOL coordinate mapping")
                        else:
                            raise Exception("Temporary object is empty")
                        
                        # Cleanup
                        cmd.delete(temp_final_obj)
                        
                    except Exception as e:
                        self._log(f"❌ PyMOL coordinate mapping failed: {e}")
                        self._log("⚠️ Coordinate update unsuccessful - structure may not be updated")
                else:
                    self._log(f"⚠️ No heavy atoms found in target object '{target_object}'")
            
            # IMPORTANTE: Rigenera idrogeni di PyMOL dopo l'aggiornamento delle coordinate
            self._log("Regenerating PyMOL hydrogens...")
            cmd.h_add(target_object)
            self._log("✅ PyMOL hydrogens regenerated")
                    
        except Exception as e:
            self._log(f"❌ Error updating coordinates: {e}")
            
    def _update_pymod_elements_coordinates(self, object_names):
        """CORREZIONE BUG 3: Aggiorna le coordinate negli elementi PyMod"""
        try:
            self._log("Updating PyMod element coordinates...")
            
            for obj_name in object_names:
                # Trova l'elemento PyMod corrispondente
                for element in self.structure_elements():
                    element_selector = element.get_pymol_selector()
                    if element_selector == obj_name:
                        try:
                            # Aggiorna le coordinate nell'elemento PyMod
                            if hasattr(element.structure, 'update_coordinates_from_pymol'):
                                element.structure.update_coordinates_from_pymol()
                                self._log(f"✅ Updated PyMod coordinates for '{obj_name}'")
                            elif hasattr(element, 'update_structure_from_pymol'):
                                element.update_structure_from_pymol()
                                self._log(f"✅ Updated PyMod structure for '{obj_name}'")
                            else:
                                self._log(f"⚠️ No suitable method found to update PyMod coordinates for '{obj_name}'. The 3D view may be updated, but the internal state is not.")
                        except Exception as e:
                            self._log(f"⚠️ Could not update PyMod coordinates for '{obj_name}': {e}")
                        break
                        
        except Exception as e:
            self._log(f"❌ Error updating PyMod elements: {e}")

    def _finalize_and_replace(self, trajectory_object, original_object_name, values):
        """
        Finalizes the structure by replacing the original object with the refined one.
        This is adapted from the robust logic in clash_and_refine.py.
        """
        import tempfile
        import os
        
        self._log("\n=== FINALIZING STRUCTURE REPLACEMENT ===")
        
        last_state = cmd.count_states(trajectory_object)
        if last_state <= 1:
            self._log("⚠️ No new states in trajectory, cannot finalize replacement.")
            return

        # 1. Save the final refined state to a temporary file
        temp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        temp_pdb_path = temp_pdb.name
        temp_pdb.close()

        try:
            # Explicitly create a new object from the LAST state of the trajectory
            temp_extract_obj = "temp_extract_final"
            cmd.create(temp_extract_obj, trajectory_object, source_state=last_state, target_state=1)
            
            # Now select the relevant atoms (by chain) from this new single-state object
            original_chains = cmd.get_chains(original_object_name)
            valid_chains = [c for c in original_chains if c] # Filter out empty chain IDs

            if not valid_chains:
                self._log(f"⚠️ Could not determine valid chains for '{original_object_name}'. Using entire object for final save.")
                final_atoms_sele = temp_extract_obj
            else:
                chain_selector = " or ".join([f"chain {c}" for c in valid_chains])
                final_atoms_sele = f"({temp_extract_obj}) and ({chain_selector})"
            
            if cmd.count_atoms(final_atoms_sele) == 0:
                self._log(f"⚠️ Could not isolate original object by chain, using entire final state as fallback.")
                final_atoms_sele = temp_extract_obj

            cmd.save(temp_pdb_path, final_atoms_sele, format="pdb")
            cmd.delete(temp_extract_obj)
            
            if not (os.path.exists(temp_pdb_path) and os.path.getsize(temp_pdb_path) > 0):
                raise Exception("Failed to save final structure to temporary file.")
            
            self._log(f"  ✓ Saved final refined state to {temp_pdb_path}")

            # 2. Find and delete the original PyMod element
            element_to_delete = None
            for element in self.structure_elements():
                if element.get_pymol_selector() == original_object_name:
                    element_to_delete = element
                    break
            
            if element_to_delete:
                self._log(f"  ✓ Found original PyMod element to replace: {element_to_delete}")
                try:
                    element_to_delete.delete()
                    self._log(f"  ✓ Deleted original PyMod element and associated PyMOL object '{original_object_name}'")
                except Exception as e:
                    self._log(f"  ⚠️ Failed to delete original PyMod element: {e}. Deleting PyMOL object as fallback.")
                    if original_object_name in cmd.get_names("objects"):
                        cmd.delete(original_object_name)
            else:
                self._log(f"  ⚠️ Could not find original PyMod element. Deleting PyMOL object '{original_object_name}' directly.")
                if original_object_name in cmd.get_names("objects"):
                    cmd.delete(original_object_name)

            # 3. Load refined structure, now with corrected naming
            base_name = re.sub(r"_chain_[A-Z0-9]+", "", original_object_name)
            refined_object_name = f"{base_name}_refined"

            if refined_object_name in cmd.get_names("objects"):
                cmd.delete(refined_object_name)

            self._log(f"  Loading refined structure into PyMOL as '{refined_object_name}'...")
            cmd.load(temp_pdb_path, refined_object_name)

            # Replicate clash_and_refine workflow: save to PyMod's dir, then open that file.
            self._log(f"  Staging refined structure for PyMod import...")
            pymod_import_filename = refined_object_name + ".pdb"
            pymod_import_path = os.path.join(self.pymod.structures_dirpath, pymod_import_filename)
            
            cmd.save(pymod_import_path, refined_object_name)
            self._log(f"  Saved to PyMod directory: {pymod_import_path}")

            self._log(f"  Importing '{refined_object_name}' into PyMod from its directory...")
            self.pymod.open_structure_file(os.path.abspath(pymod_import_path))
            
            # 4. Final UI updates
            self.pymod.main_window.gridder(clear_selection=True, update_elements=True, update_menus=True)
            self._log(f"  ✓ Imported refined structure and updated PyMod interface.")

        except Exception as e:
            self._log(f"❌ An error occurred during structure replacement: {e}")
        finally:
            # 5. Cleanup
            if os.path.exists(temp_pdb_path):
                os.unlink(temp_pdb_path)
                self._log(f"  ✓ Cleaned up temporary file.")

    def _debug_coordinate_changes(self, trajectory_object, mobile_selection, is_selection, final_state):
        """DEBUG: Verifica che le coordinate siano effettivamente cambiate durante la minimizzazione"""
        try:
            self._log("🔍 DEBUG: Checking coordinate changes in trajectory...")
            
            # Selezione per il debug (solo heavy atoms dal trajectory)
            debug_sele = f"({trajectory_object}) and not hydrogens"
            
            debug_count = cmd.count_atoms(debug_sele)
            if debug_count == 0:
                self._log("⚠️ No atoms found for coordinate debug")
                return
            
            self._log(f"Analyzing coordinate changes for {debug_count} heavy atoms...")
            
            # Ottieni coordinate stato 1 (iniziale)
            coords_initial = []
            cmd.iterate_state(1, debug_sele, "coords_initial.append([x, y, z])", space={'coords_initial': coords_initial})
            
            # Ottieni coordinate stato finale
            coords_final = []
            cmd.iterate_state(final_state, debug_sele, "coords_final.append([x, y, z])", space={'coords_final': coords_final})
            
            if len(coords_initial) != len(coords_final):
                self._log(f"❌ Coordinate count mismatch: {len(coords_initial)} initial vs {len(coords_final)} final")
                return
            
            # Calcola differenze
            total_displacement = 0
            max_displacement = 0
            moved_atoms = 0
            
            for i in range(len(coords_initial)):
                x1, y1, z1 = coords_initial[i]
                x2, y2, z2 = coords_final[i]
                
                displacement = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                total_displacement += displacement
                max_displacement = max(max_displacement, displacement)
                
                if displacement > 0.01:  # Atomo si è mosso di almeno 0.01 Å
                    moved_atoms += 1
            
            avg_displacement = total_displacement / len(coords_initial) if coords_initial else 0
            
            self._log(f"📊 COORDINATE ANALYSIS:")
            self._log(f"  Total atoms analyzed: {len(coords_initial)}")
            self._log(f"  Atoms moved > 0.01Å: {moved_atoms}")
            self._log(f"  Average displacement: {avg_displacement:.3f} Å")
            self._log(f"  Maximum displacement: {max_displacement:.3f} Å")
            
            if avg_displacement < 0.001:
                self._log("⚠️  WARNING: Very small coordinate changes - minimization may not be working")
            elif avg_displacement > 5.0:
                self._log("⚠️  WARNING: Very large coordinate changes - possible distortion")
            else:
                self._log("✅ Coordinate changes look reasonable")
                
        except Exception as e:
            self._log(f"❌ Coordinate debug failed: {e}")

    def _find_clashes(self, values, selection_name, inter_clash_selectors=None):
        min_overlap = values.get("min_overlap", 0.1)
        ignore_h = values.get("ignore_h", True)
        # Aggiunto per coerenza con clash_and_refine
        aggressive = values.get("aggressive", False)
        include_adjacent = aggressive

        if aggressive:
            min_overlap = 0.6
            ignore_h = False

        if cmd.count_atoms(selection_name) == 0:
            return []

        clashes = []

        if inter_clash_selectors and len(inter_clash_selectors) >= 2:
            selector1 = inter_clash_selectors[0]
            selector2 = inter_clash_selectors[1]
            model1 = cmd.get_model(selector1)
            model2 = cmd.get_model(selector2)

            for atom1 in model1.atom:
                sym1 = getattr(atom1, 'symbol', '').strip().upper() or atom1.name.strip()[0].upper()
                if ignore_h and sym1 == 'H': continue

                for atom2 in model2.atom:
                    sym2 = getattr(atom2, 'symbol', '').strip().upper() or atom2.name.strip()[0].upper()
                    if ignore_h and sym2 == 'H': continue

                    dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(atom1.coord, atom2.coord))
                    dist = math.sqrt(dist_sq)

                    if sym1 == 'S' and sym2 == 'S' and 1.9 <= dist <= 2.2: continue

                    vdw_sum = VDW_RADII.get(sym1, 1.7) + VDW_RADII.get(sym2, 1.7)
                    overlap = vdw_sum * values["threshold"] - dist

                    if overlap > min_overlap and dist > 0.0:
                        clashes.append({'a1': atom1, 'a2': atom2, 'overlap': overlap, 'obj1': selector1, 'obj2': selector2})
        else:
            model = cmd.get_model(selection_name)
            natoms = len(model.atom)
            for i in range(natoms):
                atom1 = model.atom[i]
                sym1 = getattr(atom1, 'symbol', '').strip().upper() or atom1.name.strip()[0].upper()
                if ignore_h and sym1 == 'H': continue

                for j in range(i + 1, natoms):
                    atom2 = model.atom[j]
                    sym2 = getattr(atom2, 'symbol', '').strip().upper() or atom2.name.strip()[0].upper()
                    if ignore_h and sym2 == 'H': continue

                    if atom1.chain == atom2.chain and atom1.resi == atom2.resi:
                        continue

                    if not include_adjacent and atom1.chain == atom2.chain and abs(int(atom1.resi) - int(atom2.resi)) == 1:
                        continue

                    dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(atom1.coord, atom2.coord))
                    dist = math.sqrt(dist_sq)

                    if sym1 == 'S' and sym2 == 'S' and 1.9 <= dist <= 2.2: continue

                    vdw_sum = VDW_RADII.get(sym1, 1.7) + VDW_RADII.get(sym2, 1.7)
                    overlap = vdw_sum * values["threshold"] - dist

                    if overlap > min_overlap and dist > 0.0:
                        clashes.append({'a1': atom1, 'a2': atom2, 'overlap': overlap, 'obj1': selection_name, 'obj2': selection_name})

        clashes.sort(key=lambda x: x['overlap'], reverse=True)
        return clashes

    def _show_bumps(self, values, clashes, group_name):
        self._cleanup_objects(["bumps_red", "bumps_orange", "bumps_yellow", group_name])
        if not clashes: return

        for clash in clashes:
            a1, a2, overlap = clash['a1'], clash['a2'], clash['overlap']
            bucket = 'bumps_yellow'
            if overlap >= 1.0: bucket = 'bumps_red'
            elif overlap >= 0.6: bucket = 'bumps_orange'

            mid_pos = [(c1 + c2) / 2.0 for c1, c2 in zip(a1.coord, a2.coord)]
            cmd.pseudoatom(bucket, pos=mid_pos, vdw=max(0.15, min(2.0, overlap * 0.6)), label=f"{overlap:.2f}")

        # ... (resto della logica di show_bumps)

    def _cleanup_objects(self, object_list):
        for obj_name in object_list:
            try:
                if '*' in obj_name:
                    cmd.delete(obj_name)
                elif obj_name in cmd.get_names('objects'):
                    cmd.delete(obj_name)
                elif obj_name in cmd.get_names('selections'):
                    cmd.delete(obj_name)
            except Exception:
                pass

    def _debug_selection_mapping(self, trajectory_object, mobile_residues):
        """Maps mobile residues to OpenBabel and colors the PyMOL object for visual debugging."""
        self._log("1. Setting up OpenBabel molecule...")
        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats('pdb', 'pdb')
        mol_string = cmd.get_pdbstr(trajectory_object)
        mol = ob.OBMol()
        obconversion.ReadString(mol, mol_string)
        mol.AddHydrogens()

        self._log("2. Identifying mobile vs. fixed atoms in OpenBabel...")
        mobile_res_set = {(chain, str(resi)) for chain, resi, _ in mobile_residues}
        mobile_atoms_info = []
        fixed_atoms_info = []

        for ob_atom in ob.OBMolAtomIter(mol):
            ob_res = ob_atom.GetResidue()
            if not ob_res: continue

            ob_chain = ob_res.GetChain()
            ob_resi = str(ob_res.GetNum())
            ob_name = ob_res.GetAtomID(ob_atom).strip()

            # Skip hydrogens for coloring
            if ob_atom.GetAtomicNum() == 1:
                continue

            atom_info = (ob_chain, ob_resi, ob_name)
            if (ob_chain, ob_resi) in mobile_res_set:
                mobile_atoms_info.append(atom_info)
            else:
                fixed_atoms_info.append(atom_info)
        
        self._log(f"  - PyMOL mobile residues: {len(mobile_res_set)}")
        self._log(f"  - OpenBabel mobile heavy atoms: {len(mobile_atoms_info)}")
        self._log(f"  - OpenBabel fixed heavy atoms: {len(fixed_atoms_info)}")

        self._log("3. Applying debug colors to PyMOL object...")
        if not mobile_atoms_info and not fixed_atoms_info:
            self._log("  - ERROR: Could not map any atoms.")
            return

        # Color fixed atoms red
        if fixed_atoms_info:
            # Process in chunks to avoid overly long selection strings
            chunk_size = 100
            for i in range(0, len(fixed_atoms_info), chunk_size):
                chunk = fixed_atoms_info[i:i+chunk_size]
                fixed_sele_parts = [f"(chain {c} and resi {r} and name {n})" for c, r, n in chunk]
                fixed_sele = " or ".join(fixed_sele_parts)
                cmd.color("red", f"({trajectory_object}) and ({fixed_sele})")

        # Color mobile atoms green
        if mobile_atoms_info:
            chunk_size = 100
            for i in range(0, len(mobile_atoms_info), chunk_size):
                chunk = mobile_atoms_info[i:i+chunk_size]
                mobile_sele_parts = [f"(chain {c} and resi {r} and name {n})" for c, r, n in chunk]
                mobile_sele = " or ".join(mobile_sele_parts)
                cmd.color("green", f"({trajectory_object}) and ({mobile_sele})")

        cmd.show("sticks")
        cmd.show("cartoon")
        self._log("✅ Debug coloring complete. Green = Mobile, Red = Fixed.")

    # The force fields the two windows offer, in the order to fall back through.
    FORCE_FIELDS = ("MMFF94s", "MMFF94", "GAFF", "UFF", "Ghemical")

    # What a force field for organic molecules is expected to cover. Anything else in the
    # molecule is worth naming when one of them cannot type it.
    ORGANIC_ELEMENTS = frozenset(("H", "C", "N", "O", "S", "P", "F", "CL", "BR", "I", "SE"))

    def _unusual_atoms(self, mol):
        """
        What in this molecule an organic force field may not know: the elements outside
        the organic set, and the atoms whose element the file does not give at all --
        MODELLER writes 'X' in the element column of a ligand it carried as a rigid body,
        so the magnesium of a template arrives as an atom of no element.

        Returns a list of phrases, each naming something the user can find in the
        structure.
        """

        from openbabel import openbabel as _ob

        def symbol_of(atomic_number):
            if hasattr(_ob, "GetSymbol"):
                return _ob.GetSymbol(atomic_number)
            return _ob.OBElementTable().GetSymbol(atomic_number)

        elements, unnamed = {}, {}
        for atom in _ob.OBMolAtomIter(mol):
            number = atom.GetAtomicNum()
            residue = atom.GetResidue()
            residue_name = residue.GetName().strip() if residue else "?"
            if not number:
                unnamed.setdefault(residue_name, 0)
                unnamed[residue_name] += 1
                continue
            symbol = symbol_of(number)
            if symbol.upper() in self.ORGANIC_ELEMENTS:
                continue
            elements.setdefault(symbol, residue_name)

        phrases = ["%s (in %s)" % (symbol, residue)
                   for symbol, residue in sorted(elements.items())]
        for residue, count in sorted(unnamed.items()):
            phrases.append("%d atom(s) of %s whose element the file does not give"
                           % (count, residue))
        return phrases

    def _setup_force_field(self, chosen, mol, mol_string, obconversion):
        """
        Set up the force field the user asked for, or the first one which can take this
        molecule. Returns (name, force field), or (name, None) when none of them can.
        """

        order = [chosen] + [name for name in self.FORCE_FIELDS if name != chosen]
        unusual = None

        for position, name in enumerate(order):
            ff = ob.OBForceField.FindForceField(name)
            if not ff:
                self._log("  '%s' is not in this Open Babel." % name)
                continue

            # A fresh molecule for each attempt: a force field which failed to type the
            # first one may have left it changed.
            if position:
                mol.Clear()
                obconversion.ReadString(mol, mol_string)
                mol.AddHydrogens()

            try:
                ok = ff.Setup(mol)
            except Exception as e:
                self._log("  '%s' raised while setting up: %s" % (name, e))
                ok = False

            if ok:
                if position:
                    self._log("  '%s' can, so the minimisation uses it instead of '%s'."
                              % (name, chosen))
                else:
                    self._log("Force field: %s." % name)
                return name, ff

            if unusual is None:
                unusual = self._unusual_atoms(mol)
            reason = ""
            if unusual:
                reason = " It contains %s." % ", ".join(unusual)
            self._log("  '%s' cannot assign parameters to this molecule.%s"
                      % (name, reason))

        self.failure_note = ("No force field could take this molecule, so it was not"
                             " minimised.")
        self._log("None of the force fields of this Open Babel can take this molecule,"
                  " so there is nothing to minimise. The clash refinement, if it was asked"
                  " for, has been done.")
        if getattr(self.pymod, "main_window", None) is not None:
            self.pymod.main_window.show_error_message(
                "Minimization Not Possible",
                "No force field available in Open Babel could assign parameters to this"
                " molecule, so it was not minimised.\n\nTried: %s." % ", ".join(order))
        return chosen, None

    def _run_openbabel_minimize_adapted(self, selection, forcefield, method, nsteps, conv, core_residue_specs, cutoff=1, cut_vdw=6.0, cut_elec=8.0, trajectory_object=None):
        """CORREZIONE COMPLETA: Metodo di minimizzazione basato su clash_and_refine.py funzionante con mappatura robusta degli atomi."""
        if not has_openbabel:
            return 0
            
        # Setup OpenBabel
        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats('pdb', 'pdb')
        
        # Get PDB string e crea la molecola
        mol_string = cmd.get_pdbstr(selection)
        mol = ob.OBMol()
        obconversion.ReadString(mol, mol_string)
        
        self._log("Adding explicit hydrogens to OpenBabel molecule...")
        mol.AddHydrogens()
        
        # Which force field can take this molecule.
        #
        # A force field which is present is not a force field which fits: MMFF94 and MMFF94s
        # are parameterised for organic molecules and have no parameters for a metal, so
        # Setup() answers false on a model which carries, say, a magnesium ion from its
        # template -- and the command used to stop there, saying the molecule "may have
        # issues". It says which atoms are unusual now, and goes on with the first force
        # field which can type them, in the order the window offers them.
        forcefield, ff = self._setup_force_field(forcefield, mol, mol_string, obconversion)
        if ff is None:
            return 0

        # ROBUST ATOM MAPPING AND CONSTRAINT SETUP
        self._log("Matching PyMOL mobile residues to OpenBabel atoms via fingerprints...")
        mobile_ob_indices = set()
        mobile_res_set = {(chain, str(resi)) for chain, resi, _ in core_residue_specs}

        for ob_atom in ob.OBMolAtomIter(mol):
            ob_res = ob_atom.GetResidue()
            if not ob_res:
                continue

            ob_chain = ob_res.GetChain()
            ob_resi = str(ob_res.GetNum())
            
            if (ob_chain, ob_resi) in mobile_res_set:
                mobile_ob_indices.add(ob_atom.GetIdx())

        self._log(f"  ✓ Found {len(mobile_ob_indices)} mobile atoms in OpenBabel molecule.")

        if not mobile_ob_indices:
            self._log("ERROR: No mobile atoms could be mapped in OpenBabel! Minimization aborted.")
            return 0

        constraints = ob.OBFFConstraints()
        total_atoms = mol.NumAtoms()
        for i in range(1, total_atoms + 1):
            if i not in mobile_ob_indices:
                constraints.AddAtomConstraint(i)
        ff.SetConstraints(constraints)
        self._log(f"Constrained {total_atoms - len(mobile_ob_indices)} atoms, {len(mobile_ob_indices)} atoms mobile")
        
        # USER-REQUESTED CHANGE: Respect the cutoff parameter.
        if cutoff:
            self._log("ENABLING cutoffs for performance...")
            ff.EnableCutOff(True)
            ff.SetVDWCutOff(cut_vdw)
            ff.SetElectrostaticCutOff(cut_elec)
            # The pair list has to be built before any energy is read. Enabling the cutoffs
            # only sets the radii: Energy() then leaves out every interaction which is not
            # in the list yet, and answers with a figure well below the energy of the
            # molecule -- 38.6 against 110.1 on a test dipeptide. Open Babel builds the
            # list itself during the descent, so the minimisation was never affected; what
            # was affected is the number reported here as the initial energy, which could
            # make an improving structure look as though its energy had risen.
            ff.UpdatePairsSimple()
        else:
            self._log("DISABLING cutoffs for stability...")
            ff.EnableCutOff(False)
        
        # Test energia iniziale per verificare setup
        try:
            initial_energy = ff.Energy()
            if math.isnan(initial_energy) or math.isinf(initial_energy):
                self._log(f"ERROR: Invalid initial energy: {initial_energy}")
                return 0
            self._log(f"Initial energy: {initial_energy:.2f}")
        except Exception as e:
            self._log(f"ERROR: Cannot calculate initial energy: {e}")
            return 0
        
        if self.progress_window is not None:
            self.progress_window.counting(nsteps)
        
        # Minimization loop con controllo errori
        self._log("Starting minimization...")
        for step in range(1, nsteps + 1):
            try:
                ff.SteepestDescent(1)
                current_energy = ff.Energy()

                # Safety check for exploding energy
                if current_energy > 1e9:
                    self._log(f"❌ WARNING: Energy has exploded at step {step} (E = {current_energy:.2f}). Minimization is unstable and will be stopped.")
                    break # Stop the loop but do not crash the script
                
                if math.isnan(current_energy) or math.isinf(current_energy):
                    self._log(f"ERROR: Invalid energy at step {step}: {current_energy}")
                    break
                
                if step % 20 == 0 or step == nsteps:
                    ff.GetCoordinates(mol)
                    pdb_out = obconversion.WriteString(mol)
                    
                    if trajectory_object and pdb_out:
                        temp_obj_name = f"tmp_mol_frame_{step}"
                        cmd.read_pdbstr(pdb_out, temp_obj_name, state=1, quiet=1)
                        heavy_atom_sele = f"({temp_obj_name}) and not hydrogens"
                        heavy_count = cmd.count_atoms(heavy_atom_sele)
                        
                        if heavy_count > 0:
                            next_state = cmd.count_states(trajectory_object) + 1
                            cmd.create(trajectory_object, heavy_atom_sele, 1, next_state)
                            self._log(f"  Saved frame {step}: {heavy_count} heavy atoms to state {next_state}")
                        else:
                            self._log(f"  WARNING: No heavy atoms found in frame {step}")
                        
                        cmd.delete(temp_obj_name)
                
                if self.progress_window is not None and (step % 5 == 0 or step == nsteps):
                    self.progress_window.step(
                        step, "Minimising — step %d of %d, energy %.2f"
                        % (step, nsteps, current_energy))
                    
            except Exception as e:
                self._log(f"ERROR during minimization at step {step}: {e}")
                break
        
        try:
            final_energy = ff.Energy()
            self._log(f"Final energy: {final_energy:.2f}, ΔE: {final_energy - initial_energy:.2f}")
        except Exception as e:
            self._log(f"ERROR calculating final energy: {e}")
            final_energy = initial_energy
        
        
        # Salva frame finale con heavy atoms
        try:
            ff.GetCoordinates(mol)
            self._log("Saving final frame with heavy atoms only...")
            final_pdb = obconversion.WriteString(mol)
            
            if trajectory_object and final_pdb:
                temp_final_obj = f"tmp_final_frame"
                cmd.read_pdbstr(final_pdb, temp_final_obj, state=1, quiet=1)
                heavy_atom_sele = f"({temp_final_obj}) and not hydrogens"
                heavy_count = cmd.count_atoms(heavy_atom_sele)
                
                if heavy_count > 0:
                    next_state = cmd.count_states(trajectory_object) + 1
                    cmd.create(trajectory_object, heavy_atom_sele, 1, next_state)
                    self._log(f"✅ Final frame saved: {heavy_count} heavy atoms to state {next_state}")
                else:
                    self._log("❌ WARNING: No heavy atoms found in final frame")
                
                cmd.delete(temp_final_obj)
            else:
                self._log("❌ ERROR: Cannot save final frame")
                
        except Exception as e:
            self._log(f"❌ ERROR updating final coordinates: {e}")
        
        self.pymod.main_window.show_info_message("Minimization", f"Energy: {initial_energy:.2f} → {final_energy:.2f}")
        return 1

def molecule_optimizer_protocol(pymod_instance, selected_elements=None):
    protocol = MoleculeOptimizer(pymod_instance, selected_elements=selected_elements)
    protocol.run()
