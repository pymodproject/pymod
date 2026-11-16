import math
import os
import re
import shutil
import pymol
from pymol import cmd
from pymod_lib.pymod_gui.clash_refine_gui import ClashRefineWindow
from pymod_lib import pymod_vars as pmdt

# Check for Open Babel
try:
    from openbabel import openbabel as ob
    has_openbabel = True
except ImportError:
    has_openbabel = False

# Try Qt for GUI progress window
try:
    from pymol.Qt import QtWidgets, QtCore
    has_qt = True
except ImportError:
    try:
        from PyQt5 import QtWidgets, QtCore
        has_qt = True
    except ImportError:
        has_qt = False


# Extended VDW radii dictionary including nucleic acid atoms
VDW_RADII = {
    # Standard atoms
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'P': 1.80, 'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
    # Additional atoms commonly found in nucleic acids
    'MG': 1.73, 'CA': 2.31, 'K': 2.75, 'NA': 2.27, 'ZN': 1.39
}

# Nucleic acid residue names for DNA and RNA - expanded for PDB variations
DNA_RESIDUES = {
    'DA', 'DT', 'DG', 'DC',  # Standard DNA
    'A', 'T', 'G', 'C',      # Short DNA names
    'ADE', 'THY', 'GUA', 'CYT',  # Full names
    'dA', 'dT', 'dG', 'dC',  # Alternative naming
}
RNA_RESIDUES = {
    'A', 'U', 'G', 'C',      # Standard RNA
    'RA', 'RU', 'RG', 'RC',  # RNA with R prefix
    'AMP', 'UMP', 'GMP', 'CMP',  # Monophosphates
    'ATP', 'UTP', 'GTP', 'CTP',  # Triphosphates
    'ADP', 'UDP', 'GDP', 'CDP',  # Diphosphates
}
NUCLEIC_RESIDUES = DNA_RESIDUES | RNA_RESIDUES

# Standard protein residues
PROTEIN_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

def is_nucleic_acid_residue(resn):
    """Check if a residue is a nucleic acid."""
    return resn.strip() in NUCLEIC_RESIDUES

def is_protein_residue(resn):
    """Check if a residue is a protein."""
    return resn.strip() in PROTEIN_RESIDUES

def debug_residue_types(model, debug_enabled=True):
    """Debug function to analyze residue types in the model."""
    if not debug_enabled:
        return
    
    residue_counts = {}
    unique_residues = set()
    
    for atom in model.atom:
        resn = atom.resn.strip()
        unique_residues.add(resn)
        
        if resn not in residue_counts:
            residue_counts[resn] = {'nucleic': 0, 'protein': 0, 'other': 0}
        
        if is_nucleic_acid_residue(resn):
            residue_counts[resn]['nucleic'] += 1
        elif is_protein_residue(resn):
            residue_counts[resn]['protein'] += 1
        else:
            residue_counts[resn]['other'] += 1
    
    print("\n=== RESIDUE DEBUG INFO ===")
    print(f"Unique residues found: {sorted(unique_residues)}")
    for resn, counts in residue_counts.items():
        total = sum(counts.values())
        nucleic_flag = "✅ NUCLEIC" if counts['nucleic'] > 0 else ""
        protein_flag = "✅ PROTEIN" if counts['protein'] > 0 else ""
        other_flag = "❓ OTHER" if counts['other'] > 0 else ""
        flags = " ".join([f for f in [nucleic_flag, protein_flag, other_flag] if f])
        print(f"  {resn}: {total} atoms ({flags})")
    print("========================\n")

def are_adjacent_nucleotides(atom1, atom2, include_adjacent=False):
    """
    Check if two atoms are from adjacent nucleotides and should be excluded
    from clash detection due to covalent bonding.
    """
    if not include_adjacent:
        return False
    
    # Only check if both atoms are from nucleic acids
    if not (is_nucleic_acid_residue(atom1.resn) and is_nucleic_acid_residue(atom2.resn)):
        return False
    
    # Must be same chain and adjacent residue numbers
    if atom1.chain != atom2.chain:
        return False
    
    try:
        resi1 = int(atom1.resi)
        resi2 = int(atom2.resi)
        if abs(resi1 - resi2) != 1:
            return False
    except ValueError:
        return False
    
    # Check for phosphodiester bond atoms (P-O bonds between nucleotides)
    phosphodiester_pairs = [
        ('P', "O3'"), ("O3'", 'P'), ('P', "O5'"), ("O5'", 'P'),
        ('P', 'OP1'), ('OP1', 'P'), ('P', 'OP2'), ('OP2', 'P'),
        ('P', 'O1P'), ('O1P', 'P'), ('P', 'O2P'), ('O2P', 'P')
    ]
    
    atom1_name = atom1.name.strip()
    atom2_name = atom2.name.strip()
    
    return (atom1_name, atom2_name) in phosphodiester_pairs

def are_normal_dna_contacts(atom1, atom2, distance_threshold=3.5):
    """
    Check if two atoms are from different DNA chains and form normal 
    double helix contacts that should not be considered clashes.
    """
    # Both must be nucleic acids
    if not (is_nucleic_acid_residue(atom1.resn) and is_nucleic_acid_residue(atom2.resn)):
        return False
    
    # Must be different chains
    if atom1.chain == atom2.chain:
        return False
    
    # Distance-based filter: normal DNA contacts are typically closer than severe clashes
    dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(atom1.coord, atom2.coord))
    dist = (dist_sq) ** 0.5
    
    # If distance is reasonable for DNA structure, likely a normal contact
    if dist < distance_threshold:
        atom1_name = atom1.name.strip()
        atom2_name = atom2.name.strip()
        
        # Base-pairing atoms and sugar-phosphate backbone interactions
        normal_contact_atoms = {
            "C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'",
            "P", "OP1", "OP2", "O1P", "O2P",
            "N1", "N2", "N3", "N4", "N6", "N7", "N9", 
            "O2", "O4", "O6", "C2", "C4", "C5", "C6", "C8"
        }
        
        # If both atoms are typical DNA atoms, likely normal contact
        if atom1_name in normal_contact_atoms and atom2_name in normal_contact_atoms:
            return True
    
    return False

def are_adjacent_amino_acids(atom1, atom2, include_adjacent=False):
    """
    Check if two atoms are from adjacent amino acids and should be excluded
    from clash detection due to peptide bonding.
    """
    if not include_adjacent:
        return False
        
    # Only check if both atoms are from proteins
    if not (is_protein_residue(atom1.resn) and is_protein_residue(atom2.resn)):
        return False
    
    # Must be same chain and adjacent residue numbers
    if atom1.chain != atom2.chain:
        return False
        
    try:
        resi1 = int(atom1.resi)
        resi2 = int(atom2.resi)
        if abs(resi1 - resi2) != 1:
            return False
    except ValueError:
        return False
    
    # Check for peptide bond atoms
    atom1_name = atom1.name.strip()
    atom2_name = atom2.name.strip()
    
    return ((atom1_name == 'C' and atom2_name == 'N') or 
            (atom1_name == 'N' and atom2_name == 'C'))

def is_cysteine_sidechain_clash(atom1, atom2):
    """
    Check if the clash involves sidechain atoms of two cysteine residues.
    Returns True if both atoms are from CYS sidechains (should be ignored).
    """
    # Both residues must be cysteine
    if not (atom1.resn.strip() == 'CYS' and atom2.resn.strip() == 'CYS'):
        return False
    
    # Define cysteine sidechain atoms (excluding backbone atoms)
    cys_sidechain_atoms = {'CB', 'SG'}  # Beta carbon and sulfur gamma
    
    atom1_name = atom1.name.strip()
    atom2_name = atom2.name.strip()
    
    # Check if both atoms are sidechain atoms
    if atom1_name in cys_sidechain_atoms and atom2_name in cys_sidechain_atoms:
        return True
    
    return False


# --- PROGRESS WINDOW CLASS ---
if has_qt:
    class MinimizationProgress(QtWidgets.QDialog):
        def __init__(self, nsteps, parent=None, unit=""):
            super().__init__(parent)
            self.unit = (" " + unit) if unit else ""
            self.setWindowTitle("OpenBabel Minimization")
            self.resize(400, 100)
            layout = QtWidgets.QVBoxLayout(self)
            self.label = QtWidgets.QLabel("Running minimization...", self)
            layout.addWidget(self.label)
            self.progress = QtWidgets.QProgressBar(self)
            self.progress.setMaximum(nsteps)
            self.progress.setValue(0)
            layout.addWidget(self.progress)
            self.initial_energy = None  # Store initial energy for ΔE calculation

        def update_progress(self, step, energy):
            self.progress.setValue(step)
            if self.initial_energy is None:
                # First step - store initial energy and show absolute value
                self.initial_energy = energy
                self.label.setText(f"Step {step} — Initial Energy: {energy:.2f}{self.unit}")
            else:
                # Subsequent steps - show ΔE (change from initial)
                delta_e = energy - self.initial_energy
                sign = "+" if delta_e >= 0 else ""
                self.label.setText(f"Step {step} — ΔE: {sign}{delta_e:.2f}{self.unit}")
            QtWidgets.QApplication.processEvents()


def refine_sculpting(target_object, clashing_residues_info, sphere=8.0, cycles=150, mobile_obj_name=None):
    """
    FIXED: Sculpting step that respects mobile object constraints.
    Only sculpts atoms from the mobile object, keeping fixed objects completely rigid.
    """
    from pymol import cmd

    core_sele_name = "clashing_residues_core"
    sculpt_area_name = "sculpt_area"

    try:
        if mobile_obj_name:
            print(f"🔒 CONSTRAINED SCULPTING: Only {mobile_obj_name} will be sculpted, other objects remain FIXED")
            
            # Filter clashing residues to only include those from mobile object
            mobile_clashing_residues = [res for res in clashing_residues_info if res[3] == mobile_obj_name]
            
            if not mobile_clashing_residues:
                print(f"No clashing residues found in mobile object {mobile_obj_name}. Skipping sculpting.")
                return
            
            print(f"DEBUG: Found {len(mobile_clashing_residues)} clashing residues in mobile object {mobile_obj_name}")
            
            # Group residues by chain for mobile object only
            clashing_by_chain = {}
            for chain, resi, resn, obj_name in mobile_clashing_residues:
                if chain not in clashing_by_chain:
                    clashing_by_chain[chain] = set()
                clashing_by_chain[chain].add(resi)
        else:
            print("🔓 UNCONSTRAINED SCULPTING: All objects may be sculpted")
            # Group residues by chain to build an efficient selection string
            clashing_by_chain = {}
            for chain, resi, resn, obj_name in clashing_residues_info:
                if chain not in clashing_by_chain:
                    clashing_by_chain[chain] = set()
                clashing_by_chain[chain].add(resi)

        core_sele_parts = []
        for chain, resi_set in clashing_by_chain.items():
            resi_str = "+".join(sorted(list(resi_set)))
            core_sele_parts.append(f"(chain {chain} and resi {resi_str})")

        # Build core selection string - NO byobj for combined target
        core_sele_string = f"/{target_object}// and ({' or '.join(core_sele_parts)})"

        print(f"DEBUG: Core selection string: {core_sele_string}")
        core_atom_count = cmd.select(core_sele_name, core_sele_string)
        print(f"DEBUG: Core selection '{core_sele_name}' selected {core_atom_count} atoms.")

        if core_atom_count > 0:
            # Around selection - NO byobj
            sculpt_sele_string = f"({core_sele_name} around {sphere}) and {target_object}"
            
            sculpt_atom_count = cmd.select(sculpt_area_name, sculpt_sele_string)
            print(f"DEBUG: Sculpt area selection '{sculpt_area_name}' selected {sculpt_atom_count} atoms.")

            if sculpt_atom_count > 0:
                try:
                    cmd.sculpt_activate(target_object)
                    
                    # Primary: Try iterate on selection
                    if cmd.count_atoms(sculpt_area_name) > 0:  # Verify exists
                        cmd.sculpt_iterate(sculpt_area_name, cycles=cycles)
                        print(f"✅ Constrained sculpting complete on {sculpt_atom_count} atoms.")
                    else:
                        raise pymol.CmdException("Sculpt area empty")
                        
                except pymol.CmdException as e:
                    if "object sculpt_area not found" in str(e):
                        print("Fallback: Sculpt on copied selection")
                        cmd.select("sculpt_fix_copy", sculpt_area_name)  # Copy if needed
                        cmd.sculpt_iterate(target_object, cycles=cycles // 2)  # Reduced on object
                        cmd.delete("sculpt_fix_copy")
                    else:
                        print(f"Sculpt failed: {e}")
                finally:
                    # Always try to deactivate sculpting on the target object
                    try:
                        cmd.sculpt_deactivate(target_object)
                    except pymol.CmdException:
                        pass
            else:
                print("Warning: Sculpt area selection is empty.")
        else:
            print("Warning: Core clashing selection is empty.")
    finally:
        # Always clean up selections
        if core_sele_name in cmd.get_names("selections"):
            cmd.delete(core_sele_name)
        if sculpt_area_name in cmd.get_names("selections"):
            cmd.delete(sculpt_area_name)
        if "temp_sculpt_mobile" in cmd.get_names("selections"):
            cmd.delete("temp_sculpt_mobile")
        if "sculpt_fix_copy" in cmd.get_names("selections"):
            cmd.delete("sculpt_fix_copy")


class ClashAndRefine:
    def __init__(self, pymod_instance, selected_elements=None):
        self.pymod = pymod_instance
        self.selected_elements = selected_elements or []
        self.objects_to_delete = []
        # What the minimisation did, for the one window at the end of the protocol.
        self.energy_report = None

    # The force fields the window offers, in the order to fall back through. GAFF comes
    # last here, where Optimize Molecule has it third: this command is pointed at
    # proteins and nucleic acids, and on those Open Babel prints "Cannot perform atom
    # type translation" while setting GAFF up, after which its van der Waals term -- the
    # one a clash refinement rests on -- comes out orders of magnitude below what UFF
    # and Ghemical find on the same atoms.
    FORCE_FIELDS = ("MMFF94s", "MMFF94", "UFF", "Ghemical", "GAFF")

    def _setup_force_field(self, chosen, mol, mol_string, obconversion):
        """
        Set up the force field the user asked for, or the first one which can take this
        structure. Returns (name, force field), or (name, None) when none of them can.

        A force field which is installed is not a force field which fits: MMFF94 and
        MMFF94s have no parameters for a magnesium ion or a phosphorylated residue, and
        a complex carrying either of them -- which is to say a good share of the
        structures this command is pointed at -- stopped at a dialog asking the user to
        name another force field, with nothing said about why the first one had failed.
        """

        order = [chosen] + [name for name in self.FORCE_FIELDS if name != chosen]

        for position, name in enumerate(order):
            ff = ob.OBForceField.FindForceField(name)
            if not ff:
                print("  '%s' is not in this Open Babel." % name)
                continue

            # A fresh molecule for each attempt: a force field which failed to type the
            # first one may have left it changed.
            if position:
                mol.Clear()
                obconversion.ReadString(mol, mol_string)

            try:
                ok = ff.Setup(mol)
            except Exception as e:
                print("  '%s' raised while setting up: %s" % (name, e))
                ok = False

            if ok:
                if position:
                    print("  '%s' can take this structure, so the minimisation uses it"
                          " instead of '%s'." % (name, chosen))
                else:
                    print("Force field: %s." % name)
                return name, ff

            print("  '%s' cannot assign parameters to this structure." % name)

        self.pymod.main_window.show_error_message(
            "Minimization Not Possible",
            "No force field available in Open Babel could assign parameters to this"
            " structure, so it was not minimised.\n\nTried: %s." % ", ".join(order))
        return chosen, None

    # 1. DISPATCHER FUNCTION
    def run_clash_refine(self):
        """Displays the dialog and calls the appropriate protocol based on user choice."""
        dialog = ClashRefineWindow(self.pymod.main_window, selected_elements=self.selected_elements, has_openbabel=has_openbabel)
        if not dialog.exec_():
            return

        values = dialog.get_values()
        action = values.get("action")

        if action == 'show':
            self.show_clashes_only_protocol(values)
        elif action == 'refine':
            self.refine_clashes_protocol(values)

    # 2. SHOW CLASHES ONLY PROTOCOL
    def show_clashes_only_protocol(self, values):
        """Runs the 'Show Clashes Only' workflow. It should not affect any other state."""
        print("--- Running Show Clashes Only Protocol ---")
        # --- Setup ---
        interacting_objects = self._get_interacting_objects()
        if not interacting_objects:
            return

        # --- Find and show clashes ---
        temp_clash_obj = "clash_refine_temp"
        self._cleanup_objects([temp_clash_obj])
        cmd.create(temp_clash_obj, " or ".join(interacting_objects))

        clashes = self.find_clashes(temp_clash_obj, values, inter_clash_selectors=interacting_objects)
        if not clashes:
            self.pymod.main_window.show_info_message("No Clashes Found", "No significant clashes were found.")
            self._cleanup_objects([temp_clash_obj])
            return

        self.show_bumps(temp_clash_obj, values, inter_clash_selectors=interacting_objects)
        self.pymod.main_window.show_info_message("Clash Display", f"{len(clashes)} clashes displayed.")

        # --- Final Cleanup for SHOW ---
        self._cleanup_objects([temp_clash_obj])
        print("--- Show Clashes Protocol Finished ---")

    # 3. REFINE CLASHES PROTOCOL
    def refine_clashes_protocol(self, values):
        """Runs the full refinement workflow, ensuring a clean state."""
        import re
        print("--- Running Refine Clashes Protocol ---")
        # --- Setup ---
        interacting_objects = self._get_interacting_objects()
        if not interacting_objects:
            return

        mobile_obj_name = None
        selected_in_dialog = values["selection"]
        if selected_in_dialog != "Complex":
            mobile_obj_name = selected_in_dialog

        # --- Store colors ---
        # The hydrogens used to be stripped from the user's own objects here and added
        # again, which changed the structures the command was pointed at before it had
        # done anything. The copy it works on, 'target_object' below, is protonated on
        # its own, so the originals are only read.
        from pymol import stored
        stored.clash_refine_original_colors = {}
        for obj_name in interacting_objects:
            if obj_name in cmd.get_names("objects"):
                colors = []
                cmd.iterate(obj_name, "colors.append((index, color))", space={'colors': colors})
                stored.clash_refine_original_colors[obj_name] = colors

        # --- Create temporary objects ---
        short_names = sorted(list(set([re.sub(r'_chain_[A-Z0-9]+', '', obj) for obj in interacting_objects])))
        base_name = "_".join(short_names)
        if not base_name:
            base_name = "complex"
        # RESTORED: Always create unified target_object to maintain intermolecular interactions
        target_object = f"{base_name}_optimization_process"
        temp_clash_obj = "clash_refine_temp"
        self._cleanup_objects([target_object, temp_clash_obj])

        source_selection = " or ".join(interacting_objects)
        if cmd.count_atoms(source_selection) == 0:
            self.pymod.main_window.show_error_message("Selection Error", f"The selection for clash detection ('{source_selection}') is empty.")
            return

        cmd.create(target_object, source_selection, source_state=1, target_state=1)
        cmd.h_add(target_object)
        cmd.create(temp_clash_obj, target_object)
        
        if mobile_obj_name:
            print(f"=== CONSTRAINED MODE: Only {mobile_obj_name} will be optimized, others remain FIXED ===")
            print(f"Fixed objects: {[obj for obj in interacting_objects if obj != mobile_obj_name]}")
            print(f"  ✓ All molecules present in {target_object} for proper energy calculations")
            print(f"  ✓ Constraints will be applied to keep fixed objects immobile")
        else:
            print("=== UNCONSTRAINED MODE: All objects will be optimized ===")
        
        # State 1: Initial structures (before any optimization)
        print(f"  ✓ State 1: Initial structures saved in '{target_object}'")
        current_state = 1

        # --- Find clashes ---
        clashes = self.find_clashes(temp_clash_obj, values, inter_clash_selectors=interacting_objects)
        if not clashes:
            # Use PyMod's GUI system for dialog instead of direct QtWidgets
            if hasattr(self.pymod.main_window, 'show_question_message'):
                proceed = self.pymod.main_window.show_question_message('No Clashes Found', "No clashes found. Proceed with contact-based refinement?")
            elif has_qt:
                proceed = QtWidgets.QMessageBox.question(self.pymod.main_window, 'No Clashes Found', "No clashes found. Proceed with contact-based refinement?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes
            else:
                proceed = False
                
            if proceed:
                clashes = self.find_contacts_as_dummy_clashes(interacting_objects)
                if not clashes:
                    self.pymod.main_window.show_info_message("No Contacts Found", "No residues found for contact-based refinement.")
                    self._cleanup_objects([target_object, temp_clash_obj])
                    return
            else:
                self._cleanup_objects([target_object, temp_clash_obj])
                return

        # --- Refine ---
        bump_objects_created = self.show_bumps(temp_clash_obj, values, inter_clash_selectors=interacting_objects)
        self.refine_clashes(target_object, temp_clash_obj, clashes, values, interacting_objects, mobile_obj_name)

        # --- Final Cleanup for REFINE ---
        self._cleanup_objects(bump_objects_created)
        print("--- Refine Clashes Protocol Finished ---")

    def find_clashes(self, selection_name, values, inter_clash_selectors=None):
        """
        Improved clash detection function that handles both proteins and nucleic acids properly.
        """
        min_overlap = values["min_overlap"]
        ignore_h = values["ignore_h"]
        aggressive = values["aggressive"]
        include_adjacent = aggressive

        if aggressive:
            min_overlap = 0.6
            ignore_h = False

        if cmd.count_atoms(selection_name) == 0:
            return []

        print(f"\n=== CLASH DETECTION DEBUG ===")
        print(f"Parameters: min_overlap={min_overlap}, ignore_h={ignore_h}, aggressive={aggressive}")
        print(f"Selection: {selection_name}")
        
        clashes = []

        if inter_clash_selectors and len(inter_clash_selectors) > 1:
            # Inter-object clash detection
            print(f"INTER-OBJECT DETECTION: {len(inter_clash_selectors)} objects")
            print(f"Objects: {inter_clash_selectors}")
            
            clash_count = 0
            skip_count = 0
            
            for i in range(len(inter_clash_selectors)):
                for j in range(i + 1, len(inter_clash_selectors)):
                    selector1 = inter_clash_selectors[i]
                    selector2 = inter_clash_selectors[j]

                    model1 = cmd.get_model(selector1)
                    model2 = cmd.get_model(selector2)
                    
                    print(f"Comparing {selector1} ({len(model1.atom)} atoms) vs {selector2} ({len(model2.atom)} atoms)")
                    
                    # Debug dei tipi di residui in entrambi i modelli
                    print(f"=== DEBUGGING {selector1} ===")
                    debug_residue_types(model1)
                    print(f"=== DEBUGGING {selector2} ===")
                    debug_residue_types(model2)

                    for atom1 in model1.atom:
                        sym1 = getattr(atom1, 'symbol', '').strip().upper() or atom1.name.strip()[0].upper()
                        if ignore_h and sym1 == 'H': 
                            continue

                        for atom2 in model2.atom:
                            sym2 = getattr(atom2, 'symbol', '').strip().upper() or atom2.name.strip()[0].upper()
                            if ignore_h and sym2 == 'H': 
                                continue

                            dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(atom1.coord, atom2.coord))
                            dist = math.sqrt(dist_sq)

                            # Skip disulfide bonds
                            if sym1 == 'S' and sym2 == 'S' and 1.9 <= dist <= 2.2: 
                                continue

                            # Skip cysteine sidechain clashes
                            if is_cysteine_sidechain_clash(atom1, atom2):
                                skip_count += 1
                                if clash_count < 3:  # Debug primi 3 casi
                                    print(f"DEBUG: Skipping cysteine sidechain clash: {atom1.resn}{atom1.resi}:{atom1.name} - {atom2.resn}{atom2.resi}:{atom2.name} (dist: {dist:.2f}Å)")
                                continue

                            # Skip normal DNA double helix contacts
                            if are_normal_dna_contacts(atom1, atom2):
                                skip_count += 1
                                if clash_count < 3:  # Debug primi 3 casi
                                    print(f"DEBUG: Skipping normal DNA contact: {atom1.resn}{atom1.resi}:{atom1.name} - {atom2.resn}{atom2.resi}:{atom2.name} (dist: {dist:.2f}Å)")
                                continue

                            vdw_sum = VDW_RADII.get(sym1, 1.7) + VDW_RADII.get(sym2, 1.7)
                            overlap = vdw_sum * values["threshold"] - dist

                            if overlap > min_overlap and dist > 0.0:
                                clash_count += 1
                                if clash_count <= 5:  # Debug primi 5 clash
                                    nucleic1 = is_nucleic_acid_residue(atom1.resn)
                                    nucleic2 = is_nucleic_acid_residue(atom2.resn)
                                    print(f"INTER-OBJECT CLASH #{clash_count}: {atom1.resn}{atom1.resi}:{atom1.name}({sym1}) - {atom2.resn}{atom2.resi}:{atom2.name}({sym2})")
                                    print(f"  Distance: {dist:.2f}Å, VDW sum: {vdw_sum:.2f}Å, Overlap: {overlap:.2f}Å")
                                    print(f"  Nucleic acid flags: {nucleic1}, {nucleic2}")
                                    print(f"  Objects: {selector1} vs {selector2}")
                                clashes.append({
                                    'a1': atom1, 'a2': atom2, 'overlap': overlap, 
                                    'obj1': selector1, 'obj2': selector2
                                })
        else:
            # Intra-object clash detection
            model = cmd.get_model(selection_name)
            natoms = len(model.atom)
            
            # Debug: analizza i tipi di residui
            debug_residue_types(model)
            
            print(f"Analyzing {natoms} atoms for intra-object clashes...")
            clash_count = 0
            skip_count = 0
            
            for i in range(natoms):
                atom1 = model.atom[i]
                sym1 = getattr(atom1, 'symbol', '').strip().upper() or atom1.name.strip()[0].upper()
                if ignore_h and sym1 == 'H': 
                    continue

                for j in range(i + 1, natoms):
                    atom2 = model.atom[j]
                    sym2 = getattr(atom2, 'symbol', '').strip().upper() or atom2.name.strip()[0].upper()
                    if ignore_h and sym2 == 'H': 
                        continue

                    # Skip atoms from the same residue
                    if atom1.chain == atom2.chain and atom1.resi == atom2.resi:
                        skip_count += 1
                        continue

                    # Handle adjacent residues based on molecule type
                    if atom1.chain == atom2.chain and abs(int(atom1.resi) - int(atom2.resi)) == 1:
                        # Skip covalently bonded atoms in adjacent residues
                        if are_adjacent_amino_acids(atom1, atom2, include_adjacent):
                            skip_count += 1
                            if clash_count < 5:  # Debug primi 5 casi
                                print(f"DEBUG: Skipping adjacent amino acids: {atom1.resn}{atom1.resi}:{atom1.name} - {atom2.resn}{atom2.resi}:{atom2.name}")
                            continue
                        if are_adjacent_nucleotides(atom1, atom2, include_adjacent):
                            skip_count += 1
                            if clash_count < 5:  # Debug primi 5 casi
                                print(f"DEBUG: Skipping adjacent nucleotides: {atom1.resn}{atom1.resi}:{atom1.name} - {atom2.resn}{atom2.resi}:{atom2.name}")
                            continue

                    dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(atom1.coord, atom2.coord))
                    dist = math.sqrt(dist_sq)

                    # Skip disulfide bonds
                    if sym1 == 'S' and sym2 == 'S' and 1.9 <= dist <= 2.2: 
                        continue

                    # Skip cysteine sidechain clashes
                    if is_cysteine_sidechain_clash(atom1, atom2):
                        skip_count += 1
                        if clash_count < 5:  # Debug primi 5 casi
                            print(f"DEBUG: Skipping cysteine sidechain clash: {atom1.resn}{atom1.resi}:{atom1.name} - {atom2.resn}{atom2.resi}:{atom2.name} (dist: {dist:.2f}Å)")
                        continue

                    vdw_sum = VDW_RADII.get(sym1, 1.7) + VDW_RADII.get(sym2, 1.7)
                    overlap = vdw_sum * values["threshold"] - dist

                    if overlap > min_overlap and dist > 0.0:
                        clash_count += 1
                        if clash_count <= 5:  # Debug primi 5 clash
                            nucleic1 = is_nucleic_acid_residue(atom1.resn)
                            nucleic2 = is_nucleic_acid_residue(atom2.resn)
                            print(f"DEBUG CLASH #{clash_count}: {atom1.resn}{atom1.resi}:{atom1.name}({sym1}) - {atom2.resn}{atom2.resi}:{atom2.name}({sym2})")
                            print(f"  Distance: {dist:.2f}Å, VDW sum: {vdw_sum:.2f}Å, Overlap: {overlap:.2f}Å")
                            print(f"  Nucleic acid flags: {nucleic1}, {nucleic2}")
                        clashes.append({'a1': atom1, 'a2': atom2, 'overlap': overlap})

        print(f"CLASH DETECTION SUMMARY:")
        if 'natoms' in locals():
            print(f"  Detection mode: INTRA-OBJECT")
            print(f"  Total atoms analyzed: {natoms}")
            print(f"  Atom pairs skipped: {skip_count}")
        else:
            print(f"  Detection mode: INTER-OBJECT")
            print(f"  Atom pairs skipped (normal contacts): {skip_count if 'skip_count' in locals() else 'N/A'}")
            print(f"  Inter-object clashes found: {clash_count if 'clash_count' in locals() else len(clashes)}")
        print(f"  Final clashes returned: {len(clashes)}")
        print(f"============================\n")

        clashes.sort(key=lambda x: x['overlap'], reverse=True)
        return clashes[:values["max_display"]]

    def show_bumps(self, selection, values, name='bump_check', inter_clash_selectors=None):
        """Show visual representation of clashes."""
        bump_objects = [name, 'bumps_red', 'bumps_orange', 'bumps_yellow']
        self._cleanup_objects(bump_objects)

        clashes = self.find_clashes(selection, values, inter_clash_selectors=inter_clash_selectors)
        if not clashes: 
            return []

        for clash in clashes:
            a1, a2, overlap = clash['a1'], clash['a2'], clash['overlap']
            bucket = 'bumps_yellow'
            if overlap >= 1.0: 
                bucket = 'bumps_red'
            elif overlap >= 0.6: 
                bucket = 'bumps_orange'

            mid_pos = [(c1 + c2) / 2.0 for c1, c2 in zip(a1.coord, a2.coord)]
            cmd.pseudoatom(bucket, pos=mid_pos, vdw=max(0.15, min(2.0, overlap * 0.6)), label=f"{overlap:.2f}")

        created_bump_objects = []
        for bucket_obj in ('bumps_red', 'bumps_orange', 'bumps_yellow'):
            if bucket_obj in cmd.get_names('objects'):
                cmd.color(bucket_obj.split('_')[-1], bucket_obj)
                cmd.show('spheres', bucket_obj)
                cmd.set('sphere_scale', 1.0, bucket_obj)
                try:
                    if cmd.count_atoms(bucket_obj) > 20:
                        cmd.hide('labels', bucket_obj)
                except pymol.CmdException:
                    pass
                created_bump_objects.append(bucket_obj)

        if created_bump_objects:
            cmd.group(name, ' or '.join(created_bump_objects))

        return [name] + created_bump_objects

    def find_contacts_as_dummy_clashes(self, interacting_objects):
        """Find inter-molecular contacts to use as dummy clashes for refinement."""
        if len(interacting_objects) < 2:
            return []
        
        contacts = []
        for i in range(len(interacting_objects)):
            for j in range(i + 1, len(interacting_objects)):
                obj1 = interacting_objects[i]
                obj2 = interacting_objects[j]
                
                # Find atoms within contact distance (typically 4-5 Å)
                contact_distance = 5.0
                
                model1 = cmd.get_model(obj1)
                model2 = cmd.get_model(obj2)
                
                for atom1 in model1.atom:
                    for atom2 in model2.atom:
                        dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(atom1.coord, atom2.coord))
                        dist = math.sqrt(dist_sq)
                        
                        if dist <= contact_distance:
                            # Create dummy clash with small overlap
                            contacts.append({
                                'a1': atom1, 'a2': atom2, 'overlap': 0.1, 
                                'obj1': obj1, 'obj2': obj2
                            })
                            
                if len(contacts) >= 20:  # Limit number of contacts
                    break
        
        return contacts[:20]

    def refine_clashes(self, target_object, temp_object, clashes, values, interacting_objects, mobile_obj_name=None):
        """Main refinement function handling both proteins and nucleic acids."""
        import re
        current_state = 1
        
        # SALVA I COLORI DEGLI ELEMENTI ORIGINALI ALL'INIZIO
        print("Saving original PyMod element colors...")
        self.original_element_colors = {}
        for obj_name in interacting_objects:
            # Find the corresponding PyMod element
            for element in self.pymod.get_pymod_elements_list():
                if hasattr(element, 'structure') and element.structure:
                    element_selector = element.get_pymol_selector()
                    if element_selector == obj_name:
                        if hasattr(element, 'my_color') and element.my_color:
                            self.original_element_colors[obj_name] = element.my_color
                            print(f"  ✓ Saved color '{element.my_color}' for '{obj_name}'")
                        else:
                            print(f"  ⚠ No color found for '{obj_name}', will use default")
                        break
        
        if values["keep_original"]:
            original_copy_name = f"{target_object}_original"
            self._cleanup_objects([original_copy_name])
            cmd.create(original_copy_name, target_object)

        clashing_residues_info = set()
        model = cmd.get_model(temp_object)
        clashing_resi_set = set()
        
        for clash in clashes:
            if 'obj1' in clash and 'obj2' in clash:
                clashing_resi_set.add((clash['a1'].chain, clash['a1'].resi, clash['obj1']))
                clashing_resi_set.add((clash['a2'].chain, clash['a2'].resi, clash['obj2']))
            else:
                clashing_resi_set.add((clash['a1'].chain, clash['a1'].resi, interacting_objects[0]))
                clashing_resi_set.add((clash['a2'].chain, clash['a2'].resi, interacting_objects[1]))

        resi_to_resn_map = {(atom.chain, atom.resi): atom.resn for atom in model.atom}

        for chain, resi, obj_name in clashing_resi_set:
            resn = resi_to_resn_map.get((chain, resi))
            if resn:
                clashing_residues_info.add((chain, resi, resn, obj_name))

        residues_to_refine = clashing_residues_info
        if mobile_obj_name:
            print(f"DEBUG: Filtering clashing residues for mobile object: {mobile_obj_name}")
            mobile_residues = set()
            for c, r, n, obj in clashing_residues_info:
                if obj == mobile_obj_name:
                    mobile_residues.add((c, r, n, obj))
            if mobile_residues:
                residues_to_refine = mobile_residues
            print(f"DEBUG: Found {len(residues_to_refine)} clashing residues for wizard/sculpting.")
        
        # Store clashing residues info for Monte Carlo targeting
        self._temp_clashing_residues = residues_to_refine
        print(f"🎯 Stored {len(residues_to_refine)} clashing residues for Monte Carlo targeting")

        # Create a temporary object for sculpting to preserve state 1
        sculpting_obj_name = f"{target_object}_sculpting"
        cmd.create(sculpting_obj_name, target_object, 1, 1)
        self.objects_to_delete.append(sculpting_obj_name)

        # Refinement with mutagenesis wizard
        try:
            cmd.wizard("mutagenesis")
            wiz = cmd.get_wizard()
            if not wiz:
                self.pymod.main_window.show_error_message("PyMOL Error", "Could not activate mutagenesis wizard.")
                return

            for chain, resi, resn, obj_name in residues_to_refine:
                sele = f"/{sculpting_obj_name}//{chain}/{resi}"
                
                # Handle nucleic acids differently from proteins
                if is_nucleic_acid_residue(resn):
                    relax_target_res = "C"  # Use cytosine for nucleic acids
                    if resn.strip() == "C":
                        continue
                else:
                    relax_target_res = "ALA"  # Use alanine for proteins
                    if resn in ["ALA", "GLY"]:
                        continue
                        
                wiz.set_mode(relax_target_res)
                wiz.do_select(sele)
                wiz.apply()
                cmd.refresh()
                cmd.deselect()
                wiz.set_mode(resn)
                wiz.do_select(sele)
                wiz.apply()

        except Exception as e:
            self.pymod.main_window.show_error_message("PyMOL Refine Error", f"Refinement with wizard failed: {e}")
        finally:
            cmd.set_wizard()

        # Local sculpting (only if enabled)
        if values.get("run_sculpting", True) and residues_to_refine:
            if mobile_obj_name:
                print(f"Performing CONSTRAINED local sculpting around clashing residues (mobile object: {mobile_obj_name})...")
            else:
                print("Performing local sculpting around clashing residues...")
            refine_sculpting(sculpting_obj_name, residues_to_refine, sphere=values.get("sculpting_radius", 6.0), cycles=values.get("sculpting_cycles", 100), mobile_obj_name=mobile_obj_name)
            
            # Update the main object with the result of sculpting
            cmd.delete(target_object)
            cmd.set_name(sculpting_obj_name, target_object)
            self.objects_to_delete.remove(sculpting_obj_name) # No longer need to delete this separately
            print(f"  ✓ Updated '{target_object}' with sculpting result.")
        elif not values.get("run_sculpting", True):
            print("⏭ Sculpting disabled by user, skipping sculpting step...")
            # Still need to clean up the sculpting object if sculpting is disabled
            cmd.delete(sculpting_obj_name)
            self.objects_to_delete.remove(sculpting_obj_name)
        else:
            print("⏭ No residues to refine, skipping sculpting step...")
            # Clean up the sculpting object
            cmd.delete(sculpting_obj_name)
            self.objects_to_delete.remove(sculpting_obj_name)

        # OpenBabel minimization (if requested)
        minimization_method = None
        if values.get("run_ob_minimize"):
            if not has_openbabel:
                self.pymod.main_window.show_error_message("Open Babel Error", "Open Babel is not installed, cannot perform minimization.")
            else:
                # Create a temporary object for minimization from the last state of target_object
                last_state = cmd.count_states(target_object)
                min_temp_obj = f"{target_object}_min_temp"
                cmd.create(min_temp_obj, target_object, source_state=last_state, target_state=1)
                self.objects_to_delete.append(min_temp_obj)

                # Now, add hydrogens with the PyMOL command on the temporary object
                print("Adding hydrogens for OpenBabel minimization...")
                cmd.h_add(min_temp_obj)

                # Determine which residues should be mobile during local minimization
                mobile_cutoff_radius = values.get("mobile_cutoff_radius", 4.0)  # Default 4.0 Å radius
                
                print(f"Using mobile cut-off radius: {mobile_cutoff_radius} Å")
                
                residues_for_ob_min = clashing_residues_info
                if mobile_obj_name:
                    print(f"DEBUG: Constraining local minimization to mobile object: {mobile_obj_name}")
                    residues_for_ob_min = {res for res in clashing_residues_info if res[3] == mobile_obj_name}

                print(f"Core clashing residues involved: {len(residues_for_ob_min)}")
                for c, r, n, obj_name in list(residues_for_ob_min)[:5]:
                    print(f"  - {obj_name} chain {c} residue {r} ({n})")
                if len(residues_for_ob_min) > 5:
                    print(f"  ... and {len(residues_for_ob_min) - 5} more")

                core_clash_selection_parts = []
                processed_residues = set() # Avoid duplicates
                print("  - Expanding mobile selection to include adjacent residues.")
                for c, r, n, obj_name in residues_for_ob_min:
                    try:
                        resi_num = int(r)
                        # Add clashing, previous, and next residues
                        for i in range(resi_num - 1, resi_num + 2):
                            if (c, i) not in processed_residues:
                                core_clash_selection_parts.append(f"(chain {c} and resi {i})")
                                processed_residues.add((c, i))
                    except ValueError:
                        # Fallback for non-integer residue identifiers
                        if (c, r) not in processed_residues:
                            core_clash_selection_parts.append(f"(chain {c} and resi {r})")
                            processed_residues.add((c, r))
                
                core_clash_selection = " or ".join(core_clash_selection_parts) if core_clash_selection_parts else ""
                
                mobile_atom_indices = []
                expanded_atom_count = 0  # Initialize to prevent undefined variable error
                
                if core_clash_selection:
                    final_mobile_selection = f"({min_temp_obj}) and ({core_clash_selection})"
                    
                    print(f"Selection strategy:")
                    print(f"  1. Core clash residues: ALWAYS included")
                    print(f"  2. Expansion: atoms within {mobile_cutoff_radius}Å of core residues")
                    
                    residue_indices = []
                    try:
                        core_atom_count = cmd.count_atoms(f"({min_temp_obj}) and ({core_clash_selection})")
                        cmd.iterate(final_mobile_selection, "residue_indices.append(index)", space={'residue_indices': residue_indices})
                        mobile_atom_indices = list(set(residue_indices))
                        expanded_atom_count = len(mobile_atom_indices)  # Set here to ensure it's defined
                        
                        # ULTRA-STRICT mobile object isolation: ONLY mobile object atoms allowed
                        if mobile_obj_name:
                            print(f"🔒 ULTRA-STRICT MODE: Applying maximum constraints for {mobile_obj_name}")
                            
                            # Get chains from mobile object in original structure
                            mobile_chains = []
                            cmd.iterate(mobile_obj_name, "mobile_chains.append(chain)", space={'mobile_chains': mobile_chains})
                            mobile_chains = list(set(mobile_chains))
                            print(f"  - Chains in original mobile object: {mobile_chains}")

                            if not mobile_chains:
                                print(f"  ❌ Cannot determine chains for mobile object!")
                                return 1

                            # Use chain-based selection instead of byobj
                            chain_selector = " or ".join([f"chain {chain}" for chain in mobile_chains])
                            mobile_only_selection = f"({min_temp_obj}) and ({chain_selector})"
                            mobile_only_count = cmd.count_atoms(mobile_only_selection)
                            print(f"  - Atoms found using chain-based selection: {mobile_only_count}")

                            if mobile_only_count == 0:
                                print(f"  ❌ Cannot identify mobile atoms using chains!")
                                return 1

                            # Intersect with clash region
                            strict_mobile_selection = f"({mobile_only_selection}) and (({core_clash_selection}) or (({core_clash_selection}) around {mobile_cutoff_radius}))"
                            
                            # Re-iterate with ultra-strict selection
                            residue_indices = []
                            cmd.iterate(strict_mobile_selection, "residue_indices.append(index)", space={'residue_indices': residue_indices})
                            mobile_atom_indices = list(set(residue_indices))
                            
                            print(f"  - Mobile-only selection: '{mobile_only_selection}'")
                            print(f"  - Atoms in mobile-only: {mobile_only_count}")
                            print(f"  - Final mobile atoms for minimization: {len(mobile_atom_indices)}")
                            
                            # VERIFICATION: Double-check no fixed object atoms are included
                            fixed_objects = [obj for obj in interacting_objects if obj != mobile_obj_name]
                            if fixed_objects and len(mobile_atom_indices) > 0:
                                print(f"  - Verifying selection against fixed objects: {fixed_objects}")
                                fixed_chains = []
                                for fixed_obj in fixed_objects:
                                    cmd.iterate(fixed_obj, "fixed_chains.append(chain)", space={'fixed_chains': fixed_chains})
                                fixed_chain_selector = " or ".join([f"chain {ch}" for ch in set(fixed_chains)])
                                
                                if fixed_chain_selector:
                                    contaminated_selection = f"index {'+'.join(map(str, mobile_atom_indices))} and ({fixed_chain_selector})"
                                    contaminated_count = cmd.count_atoms(contaminated_selection)
                                    
                                    if contaminated_count > 0:
                                        print(f"  ⚠️  WARNING: Found {contaminated_count} atoms from fixed objects in mobile selection!")
                                        print(f"  🛠️  REMOVING contaminated atoms...")
                                        
                                        clean_selection = f"index {'+'.join(map(str, mobile_atom_indices))} and not ({fixed_chain_selector})"
                                        clean_indices = []
                                        cmd.iterate(clean_selection, "clean_indices.append(index)", space={'clean_indices': clean_indices})
                                        mobile_atom_indices = list(set(clean_indices))
                                        print(f"  ✅ Cleaned mobile atoms: {len(mobile_atom_indices)}")
                                    else:
                                        print(f"  ✅ No contamination detected - all atoms belong to mobile chains.")
                            
                            # Update expanded_atom_count after potential cleanup
                            expanded_atom_count = len(mobile_atom_indices)
                        
                        print(f"Mobile region composition:")
                        print(f"  - Core clashing atoms: {core_atom_count}")
                        print(f"  - Total mobile atoms: {expanded_atom_count}")
                        print(f"  - Expansion factor: {expanded_atom_count/core_atom_count:.1f}x")
                        
                        temp_mobile_vis = "temp_mobile_preview"
                        cmd.select(temp_mobile_vis, final_mobile_selection)
                        print(f"  - Mobile selection created: '{temp_mobile_vis}' (will be deleted after minimization)")
                        
                    except Exception as e:
                        print(f"Warning: Could not determine mobile atoms: {e}")
                        mobile_selection_parts = []
                        for c, r, n, obj_name in residues_for_ob_min:
                            res_id = int(r)
                            mobile_selection_parts.append(f"(chain {c} and resi {res_id-2}-{res_id+2})")

                        full_mobile_selection = " or ".join(mobile_selection_parts)
                        if full_mobile_selection:
                            cmd.iterate(f"({min_temp_obj}) and ({full_mobile_selection})", "residue_indices.append(index)", space={'residue_indices': residue_indices})
                            mobile_atom_indices = list(set(residue_indices))

                if mobile_atom_indices:
                    core_residue_specs = {(c, str(r)) for c, r, n, obj in residues_for_ob_min}
                    minimization_method = self.run_openbabel_minimize(
                        min_temp_obj, mobile_atom_indices,
                        values["ob_forcefield"], values["ob_method"],
                        values["ob_nsteps"], values["ob_conv"], 
                        core_residue_specs,
                        values.get("ob_cutoff", 1), values.get("ob_vdw_cutoff", 6.0), values.get("ob_elec_cutoff", 8.0),
                        original_elements=None, 
                        minimization_type=values.get("ob_optimization_mode", "Auto"),
                        mc_temp=values.get("mc_temp", 300.0),
                        mc_cycles=values.get("mc_cycles", 1000),
                        target_object=target_object
                    )

                    if minimization_method == "error":
                        return

                    num_states_in_temp = cmd.count_states(min_temp_obj)
                    if num_states_in_temp > 1:
                        print(f"  Copying {num_states_in_temp - 1} frames from minimization trajectory...")
                        for i in range(2, num_states_in_temp + 1):
                            new_target_state = cmd.count_states(target_object) + 1
                            cmd.create(target_object, min_temp_obj, source_state=i, target_state=new_target_state)
                            print(f"    ✓ State {new_target_state}: After minimization saved in '{target_object}'")
                    else:
                        # Handle the case where only one state was produced (e.g., Monte Carlo)
                        new_target_state = cmd.count_states(target_object) + 1
                        cmd.create(target_object, min_temp_obj, source_state=1, target_state=new_target_state)
                        print(f"  ✓ State {new_target_state}: After minimization saved in '{target_object}'")
                else:
                    print("⚠ No atoms selected for optimization. Skipping energy minimization.")

        # ================================================================

        self._finalize_visualization(target_object, clashing_residues_info)

        # POST-REFINEMENT CLEANUP AND FINALIZATION 
        # ================================================================
        
        print("\n=== POST-REFINEMENT FINALIZATION ===")
        
        # CREATE AND IMPORT REFINED STRUCTURES (OPTIMIZED WORKFLOW)
        last_state = cmd.count_states(target_object)
        source_state_to_copy = last_state if last_state > 0 else 1
        keep_originals = values.get("keep_original", False)
        
        if keep_originals:
            print(f"Creating refined objects from state {source_state_to_copy} of '{target_object}'.")
            print("Original objects will be preserved.")
        else:
            print(f"Creating refined objects from state {source_state_to_copy} of '{target_object}'.")
            print("Original objects will be removed after creating refined versions.")

        import tempfile
        import os
        
        # Step 1: Save final state as single refined complex
        temp_cif = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
        temp_cif_path = temp_cif.name
        temp_cif.close()
        
        try:
            # Verify target_object exists and has valid states before saving
            if target_object not in cmd.get_names("objects"):
                raise Exception(f"Target object '{target_object}' no longer exists")
            
            states_count = cmd.count_states(target_object)
            if states_count == 0:
                raise Exception(f"Target object '{target_object}' has no states")
            
            atoms_count = cmd.count_atoms(target_object)
            if atoms_count == 0:
                raise Exception(f"Target object '{target_object}' has no atoms")
                
            print(f"  Target object verification: {states_count} states, {atoms_count} atoms")
            
            # Save the final state of the entire complex
            cmd.save(temp_cif_path, target_object, state=source_state_to_copy, format="cif")
            print(f"  ✓ Saved refined coordinates to temporary CIF file")
            
            if os.path.exists(temp_cif_path) and os.path.getsize(temp_cif_path) > 0:
                
                # Step 2: Create refined complex name (clean from duplicates)
                complex_base_name = "_".join(sorted(set(interacting_objects)))
                
                # Remove redundant suffixes and duplicated PDB codes
                # Remove _chain_X suffixes
                complex_base_name = re.sub(r"_chain_[A-Z]", "", complex_base_name)
                
                # Extract unique PDB codes (remove duplicates like "1atc_1atc")
                parts = complex_base_name.split("_")
                unique_parts = []
                seen = set()
                for part in parts:
                    if part not in seen and part:  # Only add if not seen and not empty
                        unique_parts.append(part)
                        seen.add(part)
                
                if unique_parts:
                    complex_base_name = "_".join(unique_parts)
                else:
                    complex_base_name = "complex"  # Fallback name
                
                refined_complex_name = f"{complex_base_name}_refined"
                
                print(f"  Creating refined complex: '{refined_complex_name}'")
                
                # Step 3: Load refined complex into PyMOL
                if refined_complex_name in cmd.get_names("objects"):
                    cmd.delete(refined_complex_name)
                cmd.load(temp_cif_path, refined_complex_name)
                self.objects_to_delete.append(refined_complex_name)
                
                # Verify single state
                final_states = cmd.count_states(refined_complex_name)
                final_atoms = cmd.count_atoms(refined_complex_name)
                print(f"  ✓ '{refined_complex_name}' created: {final_states} state(s), {final_atoms} atoms")
                
                # Step 4: Apply coloring
                cmd.util.cnc(refined_complex_name)  # Element-based coloring
                print(f"  ✓ Applied element-based coloring")
                
                # Step 5: Remove the originals -- unless the user asked to keep them.
                # "Keep a copy of the original structure" is a box in the window, and it
                # used to pick which of two congratulatory messages was shown at the end:
                # the objects and the PyMod elements were deleted either way.
                if keep_originals:
                    print("  Keeping the original objects, as the window asked.")
                else:
                    self._remove_original_elements(interacting_objects)

                # Step 6: Import refined complex into PyMod as single structure
                print(f"  Importing refined complex into PyMod...")
                try:
                    # Save refined complex to PDB file for PyMod import
                    filename = refined_complex_name + ".pdb"
                    pdb_file_shortcut = os.path.join(self.pymod.structures_dirpath, filename)
                    cmd.save(pdb_file_shortcut, refined_complex_name)
                    
                    # Import into PyMod
                    self.pymod.open_structure_file(os.path.abspath(pdb_file_shortcut))
                    print(f"  ✓ Imported '{refined_complex_name}' as single structure")
                    
                    # Update PyMod gridder
                    self.pymod.main_window.gridder(update_elements=True)
                    print(f"  ✓ PyMod interface updated")
                    
                    # Now find the refined objects that were created by PyMod import
                    # They should be named like: 1atc_refined_chain_A, 1atc_refined_chain_B, etc.
                    refined_objects = []
                    for obj_name in cmd.get_names("objects"):
                        if refined_complex_name in obj_name and "_chain_" in obj_name:
                            refined_objects.append(obj_name)
                    
                    if refined_objects:
                        print(f"  Found refined objects: {refined_objects}")
                        
                        # Remove hydrogens from refined objects
                        print(f"  Removing hydrogens from refined objects...")
                        for obj_name in refined_objects:
                            h_count_before = cmd.count_atoms(f"({obj_name}) and elem H")
                            if h_count_before > 0:
                                print(f"    Removing {h_count_before} hydrogen atoms from '{obj_name}'...")
                                cmd.remove(f"({obj_name}) and elem H")
                                h_count_after = cmd.count_atoms(f"({obj_name}) and elem H")
                                print(f"    ✓ Hydrogens removed from '{obj_name}': {h_count_before} → {h_count_after}")
                            else:
                                print(f"    No hydrogens found in '{obj_name}'")
                    else:
                        print(f"  ⚠ No refined chain objects found after import")
                    
                except Exception as e:
                    print(f"  ❌ PyMod import failed: {e}")
                    
            else:
                print(f"  ❌ Failed to create temporary CIF file")
                
        except Exception as e:
            print(f"  ❌ Refined complex creation failed: {e}")
            
        finally:
            # Clean up temporary file
            if os.path.exists(temp_cif_path):
                os.unlink(temp_cif_path)
                print(f"  ✓ Cleaned up temporary file")
        
        # 2. RIMOZIONE ATOMI DI IDROGENO dal complesso refined (fallback)
        # Use same naming logic as above
        complex_base_name = "_".join(sorted(set(interacting_objects)))
        complex_base_name = re.sub(r"_chain_[A-Z]", "", complex_base_name)
        parts = complex_base_name.split("_")
        unique_parts = []
        seen = set()
        for part in parts:
            if part not in seen and part:
                unique_parts.append(part)
                seen.add(part)
        if unique_parts:
            complex_base_name = "_".join(unique_parts)
        else:
            complex_base_name = "complex"
        refined_complex_name = f"{complex_base_name}_refined"
        
        print(f"Checking refined complex '{refined_complex_name}' for hydrogens...")
        if refined_complex_name in cmd.get_names("objects"):
            h_count_before = cmd.count_atoms(f"({refined_complex_name}) and elem H")
            if h_count_before > 0:
                print(f"Removing {h_count_before} hydrogen atoms from refined complex '{refined_complex_name}'...")
                cmd.remove(f"({refined_complex_name}) and elem H")
                h_count_after = cmd.count_atoms(f"({refined_complex_name}) and elem H")
                print(f"  ✓ Hydrogens removed from refined complex: {h_count_before} → {h_count_after}")
            else:
                print(f"  No hydrogens found in refined complex '{refined_complex_name}'")
        else:
            print(f"  ⚠ Refined complex '{refined_complex_name}' not found for hydrogen removal")
        
        # 3. RIPRISTINO COLORI ORIGINALI SEMPLICE
        print("\n--- Restoring original colors ---")
        
        # Find refined objects for color restoration
        complex_base_name = "_".join(sorted(set(interacting_objects)))
        complex_base_name = re.sub(r"_chain_[A-Z]", "", complex_base_name)
        parts = complex_base_name.split("_")
        unique_parts = []
        seen = set()
        for part in parts:
            if part not in seen and part:
                unique_parts.append(part)
                seen.add(part)
        if unique_parts:
            complex_base_name = "_".join(unique_parts)
        else:
            complex_base_name = "complex"
        refined_complex_name = f"{complex_base_name}_refined"
        
        # Find refined chain objects
        refined_objects = []
        for obj_name in cmd.get_names("objects"):
            if refined_complex_name in obj_name and "_chain_" in obj_name:
                refined_objects.append(obj_name)
        
        if refined_objects:
            print(f"Found refined objects for coloring: {refined_objects}")
            
            # Use saved colors to apply to refined objects
            colors_applied = 0
            for original_obj in interacting_objects:
                # Find corresponding refined object
                refined_obj = None
                for ref_obj in refined_objects:
                    if "_chain_" in original_obj and "_chain_" in ref_obj:
                        orig_chain = original_obj.split("_chain_")[-1]
                        ref_chain = ref_obj.split("_chain_")[-1]
                        if orig_chain == ref_chain:
                            refined_obj = ref_obj
                            break
                
                if not refined_obj:
                    print(f"  ⚠ No refined object found for '{original_obj}'")
                    continue
                
                # Use saved color
                if hasattr(self, 'original_element_colors') and original_obj in self.original_element_colors:
                    element_color = self.original_element_colors[original_obj]
                    print(f"  Applying saved color '{element_color}' from '{original_obj}' to '{refined_obj}'")
                    cmd.color(element_color, refined_obj)
                    colors_applied += 1
                else:
                    print(f"  ⚠ No saved color for '{original_obj}', using default")
                    cmd.util.cnc(refined_obj)  # Element-based coloring as fallback
            
            print(f"  ✓ Applied colors to {colors_applied}/{len(refined_objects)} refined objects")
            
            # Also apply colors to the corresponding PyMod elements
            print("  Applying colors to PyMod elements...")
            pymod_colors_applied = 0
            for original_obj in interacting_objects:
                if hasattr(self, 'original_element_colors') and original_obj in self.original_element_colors:
                    saved_color = self.original_element_colors[original_obj]
                    
                    # Find corresponding refined PyMod element
                    for element in self.pymod.get_pymod_elements_list():
                        if hasattr(element, 'structure') and element.structure:
                            element_selector = element.get_pymol_selector()
                            # Check if this element corresponds to a refined object
                            for refined_obj in refined_objects:
                                if "_chain_" in original_obj and "_chain_" in refined_obj:
                                    orig_chain = original_obj.split("_chain_")[-1]
                                    ref_chain = refined_obj.split("_chain_")[-1]
                                    if (orig_chain == ref_chain and 
                                        element_selector == refined_obj):
                                        element.my_color = saved_color
                                        print(f"    ✓ Applied color '{saved_color}' to PyMod element '{element_selector}'")
                                        pymod_colors_applied += 1
                                        break
            
            print(f"  ✓ Applied colors to {pymod_colors_applied} PyMod elements")
            
            # Force PyMod interface update to show new colors
            if pymod_colors_applied > 0:
                try:
                    self.pymod.main_window.gridder(update_elements=True)
                    print("  ✓ Updated PyMod interface to show new colors")
                except Exception as e:
                    print(f"  ⚠ Failed to update PyMod interface: {e}")
        else:
            print("  No refined objects found for color restoration")
        
        # 4. DESELEZIONA OGGETTI PYMOD
        print("\n--- Deselecting PyMod objects ---")
        try:
            # Deseleziona tutti gli elementi selezionati in PyMod
            if hasattr(self, 'selected_elements') and self.selected_elements:
                for element in self.selected_elements:
                    if hasattr(element, 'set_selected'):
                        element.set_selected(False)
                        print(f"  ✓ Deselected PyMod element")
                print(f"  ✓ Deselected {len(self.selected_elements)} PyMod elements")
            
            # Deseleziona anche in PyMOL
            cmd.deselect()
            print("  ✓ Cleared PyMOL selection")
            
        except Exception as e:
            print(f"  ⚠ Warning during deselection: {e}")
        
        # 5. PULIZIA FINALE MIGLIORATA
        print("\n--- Final cleanup ---")
        
        # Lista di oggetti temporanei da rimuovere
        cleanup_objects = [
            temp_object,
            "clash_refine_temp", 
            "mobile_residues_to_minimize",
            "temp_mobile_preview",
            "clashing_residues_core",
            "sculpt_area"
        ]
        
        # Aggiungi oggetti bump se esistono
        bump_objects = ["bump_check", "bumps_red", "bumps_orange", "bumps_yellow"]
        cleanup_objects.extend(bump_objects)
        
        # Rimuovi selezioni temporanee
        cleanup_selections = [
            "clash_refine_complex",
            "mobile_residues_to_minimize", 
            "temp_mobile_preview",
            "clashing_residues_core",
            "sculpt_area"
        ]
        
        self._cleanup_objects(cleanup_objects)
        
        for sel_name in cleanup_selections:
            if sel_name in cmd.get_names("selections"):
                cmd.delete(sel_name)
                print(f"  ✓ Removed selection: {sel_name}")

        self._cleanup_objects(self.objects_to_delete)
        
        # 6. RINOMINA E COLORA OGGETTO TRAIETTORIA (se minimizzazione eseguita)
        if minimization_method == "full" and target_object in cmd.get_names("objects"):
            final_traj_name = "local_minimization_trajectory"
            if final_traj_name in cmd.get_names("objects"):
                cmd.delete(final_traj_name)
            cmd.set_name(target_object, final_traj_name)
            print(f"  ✓ Saved minimization trajectory as '{final_traj_name}'")
            
            # Applica colorazione intelligente alla traiettoria
            try:
                self._apply_trajectory_coloring_enhanced(
                    final_traj_name, 
                    clashes, 
                    clashing_residues_info, 
                    interacting_objects,
                    values.get("mobile_cutoff_radius", 4.0)
                )
            except Exception as e:
                print(f"  ⚠ Enhanced trajectory coloring failed: {e}")
                # Fallback alla colorazione semplice
                try:
                    cmd.color("gray", final_traj_name)
                    print(f"  ✓ Applied fallback coloring to trajectory")
                except:
                    pass
        elif minimization_method is None and target_object in cmd.get_names("objects"):
            # No minimization was run, so delete the temporary object
            cmd.delete(target_object)
            print(f"  ✓ Removed temporary object: {target_object}")
        
        # 7. VERIFICA FINALE SICURA
        print("\n--- Final verification ---")
        remaining_objects = cmd.get_names("objects")
        remaining_selections = cmd.get_names("selections")
        
        print(f"Objects remaining: {len(remaining_objects)}")
        for obj in remaining_objects[:5]:  # Mostra solo i primi 5
            try:
                # Verifica che l'oggetto esista prima di contare gli atomi
                if obj in cmd.get_names("objects"):
                    atom_count = cmd.count_atoms(obj)
                    print(f"  {obj}: {atom_count} atoms")
                else:
                    print(f"  {obj}: ⚠ Object not found")
            except Exception as e:
                print(f"  {obj}: ❌ Error counting atoms: {e}")
        if len(remaining_objects) > 5:
            print(f"  ... and {len(remaining_objects) - 5} more objects")
            
        print(f"Selections remaining: {len(remaining_selections)}")
        if remaining_selections:
            print(f"  {', '.join(remaining_selections[:5])}")
        
        print("\n=== REFINEMENT FINALIZATION COMPLETE ===")
        
        # 8. MESSAGGIO FINALE ALL'UTENTE
        keep_originals = values.get("keep_original", False)
        # Use same naming logic as above
        complex_base_name = "_".join(sorted(set(interacting_objects)))
        complex_base_name = re.sub(r"_chain_[A-Z]", "", complex_base_name)
        parts = complex_base_name.split("_")
        unique_parts = []
        seen = set()
        for part in parts:
            if part not in seen and part:
                unique_parts.append(part)
                seen.add(part)
        if unique_parts:
            complex_base_name = "_".join(unique_parts)
        else:
            complex_base_name = "complex"
        refined_complex_name = f"{complex_base_name}_refined"
        
        # One window at the end, with what was done in it. There used to be four of
        # these messages, one for each combination of sculpting or minimisation and
        # keeping the originals or not, listing what had been done with a tick in front
        # of each line -- and before any of them, the minimisation had opened a window
        # of its own carrying the energies.
        if values.get("run_ob_minimize"):
            lines = ["Refinement with minimization complete.\n"]
        else:
            lines = ["Sculpting complete.\n"]

        lines.append("Refined structure: %s" % refined_complex_name)
        if keep_originals:
            lines.append("Original structures kept: %s"
                         % ", ".join(interacting_objects))
        else:
            lines.append("Original structures removed: %s"
                         % ", ".join(interacting_objects))
        # Named rather than described: the trajectory is kept under a name of its own
        # after a full minimisation, and the object the work happened in stays behind
        # after a torsion one.
        left_behind = [name for name in ("local_minimization_trajectory", target_object)
                       if name in cmd.get_names("objects")]
        if left_behind:
            lines.append("Also in the scene: %s." % ", ".join(left_behind))

        self.pymod.main_window.show_info_message(
            "Clash Refinement Complete",
            "\n".join(lines) + "\n" + self._energy_summary())


    def _get_interacting_objects(self):
        """Determines the target objects for the operation."""
        if len(self.selected_elements) >= 2:
            return [e.get_pymol_selector() for e in self.selected_elements]
        elif len(self.selected_elements) == 1:
            selected_obj_name = self.selected_elements[0].get_pymol_selector()
            all_objects = cmd.get_names("objects")
            undesired_suffixes = [
                "_minimization_process", "_original", "_global_min_debug", "_min_state", 
                "bump_check", "bumps_red", "bumps_orange", "bumps_yellow", "clash_refine_temp"
            ]
            real_objects = [name for name in all_objects if not any(name.endswith(s) for s in undesired_suffixes)]
            if len(real_objects) != 2:
                self.pymod.main_window.show_error_message(
                    "Ambiguous Selection", 
                    f"You selected one object ({selected_obj_name}), but there are {len(real_objects)} total objects in the scene. Please select the two specific interacting objects."
                )
                return None
            partner_obj = [name for name in real_objects if name != selected_obj_name][0]
            return [selected_obj_name, partner_obj]
        else:
            self.pymod.main_window.show_error_message("Selection Error", "Please select one or two objects to analyze.")
            return None

    def should_use_torsion_only_minimization(self, selection):
        """
        Determines if torsion-only minimization should be used based on molecular composition.
        Returns True for proteins/nucleic acids, False for primarily ligands.
        """
        protein_atoms = cmd.count_atoms(f"({selection}) and (polymer.protein)")
        nucleic_atoms = cmd.count_atoms(f"({selection}) and (polymer.nucleic)")
        ligand_atoms = cmd.count_atoms(f"({selection}) and (organic)")

        # If there are significant protein or nucleic acid components, default to torsion-only
        if protein_atoms > 0 or nucleic_atoms > 0:
            return True
        
        # If only ligands (or small organic molecules) are present
        if ligand_atoms > 0:
            return False
        
        # Default to torsion-only if no clear composition (e.g., single atom, or unknown type)
        return True

    def run_openbabel_minimize(self, selection, mobile_atom_indices, forcefield, method, nsteps, conv, core_residue_specs, cutoff=1, cut_vdw=6.0, cut_elec=8.0, original_elements=None, minimization_type="Auto", mc_temp=300.0, mc_cycles=1000, target_object=None):
        """Run OpenBabel minimization (torsion-only or full) based on the selected mode."""
        
        # Determine minimization type based on user selection
        if minimization_type == "Auto":
            use_torsion_only = self.should_use_torsion_only_minimization(selection)
            actual_type = "Torsion-only" if use_torsion_only else "Full"
        elif minimization_type == "Torsion Only (Monte Carlo)":
            use_torsion_only = True
            actual_type = "Torsion-only"
        elif minimization_type == "Full (CG/SD)":
            use_torsion_only = False
            actual_type = "Full"
        else:
            # Default fallback
            use_torsion_only = True
            actual_type = "Torsion-only (Default)"
        
        print(f"\n=== {actual_type.upper()} MINIMIZATION ===")
        
        if use_torsion_only:
            print("Using Monte Carlo torsion optimization to preserve ring geometries and bond lengths.")
            print("Suitable for proteins, nucleic acids, and complex macromolecules.")
            exit_status = self.run_torsion_only_minimize(selection, mobile_atom_indices, forcefield, method, nsteps, conv, mc_temp, mc_cycles, core_residue_specs, target_object=target_object, cut_vdw=cut_vdw, cut_elec=cut_elec)
            if exit_status != 0:
                return "error"
            return "torsion"
        else:
            print("Using full energy minimization for all degrees of freedom.")
            print("Suitable for small molecules, ligands, and simple compounds.")
            exit_status = self.run_full_minimize(selection, mobile_atom_indices, forcefield, method, nsteps, conv, cutoff, cut_vdw, cut_elec, original_elements)
            if exit_status != 0:
                return "error"
            return "full" 

    def run_torsion_only_minimize(self, selection, mobile_atom_indices, forcefield, method, nsteps, conv, mc_temp, mc_cycles, core_residue_specs, target_object=None, cut_vdw=12.0, cut_elec=18.0):
        """Run torsion-only minimization using Monte Carlo optimization."""
        mobile_sele_name = "mobile_residues_to_minimize"
        
        from pymol import stored
        stored.original_colors = {}
        progress_dialog = None
        exit_status = 0
        
        try:
            # Setup visual feedback
            sele_string = f"({selection}) and index { '+'.join(map(str, mobile_atom_indices)) }"
            cmd.select(mobile_sele_name, sele_string)
            cmd.iterate(mobile_sele_name, "stored.original_colors[index] = color")
            
            print(f"Highlighting {len(mobile_atom_indices)} atoms for torsion optimization...")
            cmd.color("yellow", mobile_sele_name)
            cmd.show("spheres", mobile_sele_name)
            cmd.zoom(mobile_sele_name, buffer=3.0)
            cmd.refresh()
            
            # 1. Collect detailed mobile atom fingerprints using iterate_state for coordinate access
            def collect_mobile_fp(chain, resi, name, x, y, z, mobile_atoms_fp):
                """Helper to safely append atom fingerprint."""
                if name.strip():
                    mobile_atoms_fp.append({
                        'chain': chain.strip() if chain else '',
                        'resi': resi.strip() if resi else '',
                        'name': name.strip(),
                        'coord': (float(x), float(y), float(z))
                    })

            mobile_atoms_fp = []
            mobile_sele_str = f"({selection}) and index {'+'.join(map(str, mobile_atom_indices))}"
            iterate_expr = "collect_mobile_fp(chain, resi, name, x, y, z, mobile_atoms_fp)"
            
            try:
                # Use iterate_state on state 1, which is the state we are optimizing
                cmd.iterate_state(1, mobile_sele_str, iterate_expr, space={'collect_mobile_fp': collect_mobile_fp, 'mobile_atoms_fp': mobile_atoms_fp})
                
                print(f"✅ Collected {len(mobile_atoms_fp)} mobile atom fingerprints from PyMOL.")
                
                if len(mobile_atoms_fp) == 0 and len(mobile_atom_indices) > 0:
                     print("🚨 ERROR: Fingerprint collection resulted in 0 atoms. Aborting.")
                     return 1

            except Exception as e:
                print(f"❌ FATAL: iterate_state failed while collecting fingerprints: {e}")
                return 1

            # 2. Add hydrogens and export PDB string for OpenBabel
            print("Adding hydrogens for OpenBabel optimization...")
            cmd.h_add(selection)
            mol_string = cmd.get_pdbstr(f"{selection} and state 1")

            # 3. Setup OpenBabel
            obconversion = ob.OBConversion()
            obconversion.SetInAndOutFormats('pdb', 'pdb')
            mol = ob.OBMol()
            if not obconversion.ReadString(mol, mol_string):
                self.pymod.main_window.show_error_message("Open Babel Error", "Failed to read molecule into OpenBabel.")
                return 1

            forcefield, ff = self._setup_force_field(forcefield, mol, mol_string,
                                                     obconversion)
            if ff is None:
                return 1

            # 4. Robustly match PyMOL fingerprints to OB atoms and apply constraints
            from math import sqrt
            mobile_ob_indices_set = set()
            pymol_fp_map = {(fp['chain'], fp['resi'], fp['name']): fp['coord'] for fp in mobile_atoms_fp}
            
            constraints = ob.OBFFConstraints()
            for ob_atom in ob.OBMolAtomIter(mol):
                ob_res = ob_atom.GetResidue()
                if not ob_res:
                    constraints.AddAtomConstraint(ob_atom.GetIdx())
                    continue

                ob_chain = ob_res.GetChain()
                ob_resi = str(ob_res.GetNum())
                ob_name = ob_res.GetAtomID(ob_atom).strip()
                key = (ob_chain, ob_resi, ob_name)

                is_mobile = False
                if key in pymol_fp_map:
                    pymol_coord = pymol_fp_map[key]
                    if pymol_coord is None:
                        is_mobile = True
                    else:
                        ob_coord_vec = ob_atom.GetVector()
                        dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(pymol_coord, (ob_coord_vec.GetX(), ob_coord_vec.GetY(), ob_coord_vec.GetZ())))
                        if dist_sq < 0.25: # 0.5 Angstrom tolerance
                            is_mobile = True
                
                if is_mobile:
                    mobile_ob_indices_set.add(ob_atom.GetIdx())
                
                if not is_mobile:
                    constraints.AddAtomConstraint(ob_atom.GetIdx())

            ff.SetConstraints(constraints)
            # Cutoffs are deliberately off on this path: the radii are set so that they
            # are ready if it is turned back on, but the energies below are the full sum.
            ff.EnableCutOff(False)
            ff.SetVDWCutOff(float(cut_vdw))
            ff.SetElectrostaticCutOff(float(cut_elec))
            print(f"⚡️ Energy cutoffs disabled; radii kept at VDW={cut_vdw}Å, Elec={cut_elec}Å")
            
            self._temp_mobile_indices = mobile_ob_indices_set # Store for MC

            # 5. Report matching results and verify
            match_rate = (len(mobile_ob_indices_set) / len(mobile_atoms_fp)) * 100 if mobile_atoms_fp else 0
            print(f"Matching results: {len(mobile_ob_indices_set)}/{len(mobile_atoms_fp)} atoms matched ({match_rate:.1f}%)")
            if not mobile_ob_indices_set:
                print("ERROR: No mobile atoms identified in OpenBabel. Skipping optimization.")
                return 0

            initial_energy = ff.Energy()
            # In the unit this force field answers in: GAFF, UFF and Ghemical answer in
            # kJ/mol, MMFF94 and MMFF94s in kcal/mol, and the figures were all labelled
            # kcal/mol whatever they were.
            energy_unit = ff.GetUnit()
            print(f"Initial energy: {initial_energy:.3f} {energy_unit}")

            if has_qt:
                # Parent to the PyMod main window so the dialog inherits PyMod's
                # stylesheet and is properly owned (otherwise it appears unstyled
                # and can misbehave with focus/the event loop).
                progress_dialog = MinimizationProgress(mc_cycles, parent=self.pymod.main_window,
                                                      unit=energy_unit)
                progress_dialog.show()
                progress_dialog.update_progress(0, initial_energy)

            # 6. Find rotatable torsions using the matched mobile OB atom indices
            dihedrals = self.get_rotatable_dihedrals_for_core_residues(mol, core_residue_specs, mobile_ob_indices_set)
            
            print(f"Found {len(dihedrals)} rotatable torsions in mobile atoms.")
            
            # Single Monte Carlo run (iterative logic moved to wrapper)
            optimized_mol, best_energy = self.monte_carlo_torsion_minimize(
                mol, dihedrals, ff, mc_cycles, selection, mobile_sele_name, progress_dialog, mc_temp, target_object
            )
            
            # Update PyMOL with the best structure found
            final_pdb = obconversion.WriteString(optimized_mol)
            if final_pdb.strip():
                cmd.read_pdbstr(final_pdb, selection, state=1) 
                cmd.refresh()
                print(f"  ✓ Updated '{selection}' with the lowest energy structure found.")
            
            print(f"Torsion optimization completed!")
            print(f"Final energy: {best_energy:.3f} {energy_unit}")
            
            self._record_energies("Torsion optimization", initial_energy, best_energy,
                                  energy_unit)

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.pymod.main_window.show_error_message("Open Babel Error", f"Minimization failed: {e}")

        finally:
            if has_qt and progress_dialog:
                progress_dialog.close()
            if mobile_sele_name in cmd.get_names("selections"):
                cmd.delete(mobile_sele_name)

        return exit_status

    def get_mobile_residue_specs(self, selection, mobile_atom_indices):
        """Estrai (chain, resi) da PyMOL selection - stabile, no indici post-h_add."""
        specs = set()
        cmd.iterate(f"({selection}) and index {'+'.join(map(str, mobile_atom_indices))}", 
                    "specs.add((chain.strip(), resi.strip()))", 
                    space={'specs': specs})
        return specs

    def get_rotatable_dihedrals_for_core_residues(self, mol, core_residue_specs, mobile_ob_indices):
        """
        Finds specific, valid side-chain torsions (chi angles) for amino acids
        and glycosidic/backbone torsions for nucleic acids
        to prevent protein/nucleic acid distortion. This version uses predefined dictionaries
        of torsions to ensure only correct torsions are rotated.
        """
        import random
        mobile_set = set(mobile_ob_indices)
        dihedrals = []

        # Predefined chi angles for standard amino acids
        PROTEIN_CHI_ANGLES = {
            'ARG': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'NE'), ('CG', 'CD', 'NE', 'CZ'), ('CD', 'NE', 'CZ', 'NH1')],
            'ASN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')],
            'ASP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')],
            'CYS': [('N', 'CA', 'CB', 'SG')],
            'GLN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'OE1')],
            'GLU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'OE1')],
            'HIS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'ND1')],
            'ILE': [('N', 'CA', 'CB', 'CG1'), ('CA', 'CB', 'CG1', 'CD1')],
            'LEU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
            'LYS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'CE'), ('CG', 'CD', 'CE', 'NZ')],
            'MET': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'SD'), ('CB', 'CG', 'SD', 'CE')],
            'PHE': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
            'PRO': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')],
            'SER': [('N', 'CA', 'CB', 'OG')],
            'THR': [('N', 'CA', 'CB', 'OG1')],
            'TRP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
            'TYR': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
            'VAL': [('N', 'CA', 'CB', 'CG1')],
        }

        # Predefined torsions for nucleic acids (DNA/RNA)
        # Glycosidic bond (χ) and key backbone torsions
        NUCLEIC_TORSIONS = {
            # Purines: A, G (χ: O4'-C1'-N9-C4)
            'A': [('O4\'', 'C1\'', 'N9', 'C4')],
            'G': [('O4\'', 'C1\'', 'N9', 'C4')],
            # Pyrimidines: C, T/U (χ: O4'-C1'-N1-C2)
            'C': [('O4\'', 'C1\'', 'N1', 'C2')],
            'T': [('O4\'', 'C1\'', 'N1', 'C2')],
            'U': [('O4\'', 'C1\'', 'N1', 'C2')],
            # Additional backbone torsions (common for flexibility: γ, δ, ε)
            # These are added for all nucleic residues if present
            'A': [('O4\'', 'C1\'', 'N9', 'C4'), ('C5\'', 'C4\'', 'C3\'', 'O3\'')],  # χ + γ
            'G': [('O4\'', 'C1\'', 'N9', 'C4'), ('C5\'', 'C4\'', 'C3\'', 'O3\'')],  # χ + γ
            'C': [('O4\'', 'C1\'', 'N1', 'C2'), ('C5\'', 'C4\'', 'C3\'', 'O3\'')],  # χ + γ
            'T': [('O4\'', 'C1\'', 'N1', 'C2'), ('C5\'', 'C4\'', 'C3\'', 'O3\'')],  # χ + γ
            'U': [('O4\'', 'C1\'', 'N1', 'C2'), ('C5\'', 'C4\'', 'C3\'', 'O3\'')],  # χ + γ
        }
        
        print(f"Finding specific torsions for {len(core_residue_specs)} core residues (proteins + nucleic acids)...")
        core_resi_nums = {int(resi_str) for _, resi_str in core_residue_specs}

        for res in ob.OBResidueIter(mol):
            if res.GetNum() not in core_resi_nums:
                continue
            
            res_name = res.GetName().strip().upper()
            
            # Create a mapping of atom name to atom index for the current residue
            atom_map = {res.GetAtomID(atom).strip(): atom.GetIdx() for atom in ob.OBResidueAtomIter(res)}

            # Check for protein chi angles
            if res_name in PROTEIN_CHI_ANGLES:
                for chi_atom_names in PROTEIN_CHI_ANGLES[res_name]:
                    try:
                        atom_indices = [atom_map[name] for name in chi_atom_names]
                        if all(idx in mobile_set for idx in atom_indices):
                            dihedrals.append(tuple(atom_indices))
                            print(f"  ✓ Protein chi angle for {res_name}{res.GetNum()}: {chi_atom_names}")
                    except KeyError:
                        continue

            # Check for nucleic acid torsions
            if res_name in NUCLEIC_TORSIONS:
                for torsion_atom_names in NUCLEIC_TORSIONS[res_name]:
                    try:
                        atom_indices = [atom_map[name] for name in torsion_atom_names]
                        if all(idx in mobile_set for idx in atom_indices):
                            dihedrals.append(tuple(atom_indices))
                            print(f"  ✓ Nucleic torsion for {res_name}{res.GetNum()}: {torsion_atom_names}")
                    except KeyError:
                        continue

        dihedrals = list(set(dihedrals))  # Remove duplicates
        print(f"✅ Found {len(dihedrals)} unique torsions (proteins + nucleic acids).")

        if len(dihedrals) > 100:
            print(f"⚠️ Torsion count high ({len(dihedrals)}). Sampling 100 random.")
            dihedrals = random.sample(dihedrals, 100)
        
        return dihedrals

    def get_dihedral_angle(self, mol, dih):
        """Get dihedral angle in degrees."""
        a, b, c, d = dih
        return math.degrees(mol.GetTorsion(a, b, c, d))

    def set_dihedral_angle(self, mol, dih, angle_deg):
        """Set dihedral angle in degrees."""
        a, b, c, d = dih
        mol.SetTorsion(a, b, c, d, math.radians(angle_deg))

    def compute_energy_safe(self, mol, ff):
        """Compute energy with error handling, preserving constraints."""
        try:
            # Don't call Setup again as it would reset constraints
            # The force field should already be set up with the original molecule
            # and we only need to update coordinates
            ff.SetCoordinates(mol)
            return ff.Energy()
        except:
            # Fallback: try full setup if coordinate update fails
            try:
                if not ff.Setup(mol):
                    return float('inf')
                return ff.Energy()
            except:
                return float('inf')

    def check_mobile_clashes(self, mol, mobile_indices, values):
        """
        Check for clashes in mobile atoms vs all.
        """
        if not mobile_indices:
            return []
        
        mobile_atoms = [atom for atom in ob.OBMolAtomIter(mol) if atom.GetIdx() in mobile_indices]
        all_atoms = list(ob.OBMolAtomIter(mol))
        
        clashes = []
        min_overlap = values.get("min_overlap", 0.3)
        ignore_h = values.get("ignore_h", True)
        threshold = values.get("threshold", 0.8)
        
        atomic_num_to_symbol = {
            1: 'H', 6: 'C', 7: 'N', 8: 'O', 15: 'P', 16: 'S', 
            9: 'F', 17: 'CL', 35: 'BR', 53: 'I', 12: 'MG', 
            20: 'CA', 19: 'K', 11: 'NA', 30: 'ZN'
        }
        
        for mobile_atom in mobile_atoms:
            if ignore_h and atomic_num_to_symbol.get(mobile_atom.GetAtomicNum(), 'C') == 'H':
                continue
            
            mobile_pos = mobile_atom.GetVector()
            mobile_res = mobile_atom.GetResidue()
            if not mobile_res:
                continue
            
            for other_atom in all_atoms:
                if other_atom.GetIdx() == mobile_atom.GetIdx():
                    continue
                
                other_res = other_atom.GetResidue()
                if not other_res or mobile_res.GetNum() == other_res.GetNum():
                    continue
                
                if ignore_h and atomic_num_to_symbol.get(other_atom.GetAtomicNum(), 'C') == 'H':
                    continue
                
                other_pos = other_atom.GetVector()
                dist = mobile_pos.distSq(other_pos) ** 0.5
                
                sym1 = atomic_num_to_symbol.get(mobile_atom.GetAtomicNum(), 'C')
                sym2 = atomic_num_to_symbol.get(other_atom.GetAtomicNum(), 'C')
                vdw_sum = VDW_RADII.get(sym1, 1.7) + VDW_RADII.get(sym2, 1.7)
                overlap = vdw_sum * threshold - dist
                
                if overlap > min_overlap and dist > 0.0:
                    clashes.append({'mobile': mobile_atom.GetIdx(), 'other': other_atom.GetIdx(), 'overlap': overlap})
        
        return clashes

    def get_clashing_residues_from_mobile_atoms(self, mol, mobile_atom_indices):
        """
        Identify clashing residues only from mobile atoms that have torsions.
        This focuses energy calculation on the actual optimization region.
        """
        if not mobile_atom_indices:
            print("⚠️  No mobile atoms provided for clash detection")
            return set()
            
        # Get all atoms
        atoms = list(ob.OBMolAtomIter(mol))
        mobile_atoms = [atom for atom in atoms if atom.GetIdx() in mobile_atom_indices]
        
        print(f"🔍 Clash detection in mobile region:")
        print(f"  - Total atoms in molecule: {len(atoms)}")
        print(f"  - Mobile atoms with torsions: {len(mobile_atoms)}")
        
        # Element symbol mapping
        atomic_num_to_symbol = {
            1: 'H', 6: 'C', 7: 'N', 8: 'O', 15: 'P', 16: 'S', 
            9: 'F', 17: 'CL', 35: 'BR', 53: 'I', 12: 'MG', 
            20: 'CA', 19: 'K', 11: 'NA', 30: 'ZN'
        }
        
        clashing_residues = set()
        clash_count = 0
        
        # Check clashes only between mobile atoms and all other atoms
        for mobile_atom in mobile_atoms:
            mobile_pos = mobile_atom.GetVector()
            mobile_res = mobile_atom.GetResidue()
            if not mobile_res:
                continue
                
            for other_atom in atoms:
                # Skip same atom
                if mobile_atom.GetIdx() == other_atom.GetIdx():
                    continue
                    
                other_res = other_atom.GetResidue()
                if not other_res:
                    continue
                    
                # Skip same residue
                if mobile_res.GetNum() == other_res.GetNum():
                    continue
                
                # Calculate distance
                other_pos = other_atom.GetVector()
                dist = mobile_pos.distSq(other_pos) ** 0.5
                
                # Get VDW radii
                sym1 = atomic_num_to_symbol.get(mobile_atom.GetAtomicNum(), 'C')
                sym2 = atomic_num_to_symbol.get(other_atom.GetAtomicNum(), 'C')
                vdw1 = VDW_RADII.get(sym1, 1.7)
                vdw2 = VDW_RADII.get(sym2, 1.7)
                vdw_sum = vdw1 + vdw2
                
                # Check for clash (more strict: 75% of VDW sum)
                if dist < vdw_sum * 0.75:
                    clash_count += 1
                    # Add the mobile residue to clashing residues
                    clashing_residues.add(mobile_res.GetNum())
                    if clash_count <= 10:  # Debug first 10 clashes
                        print(f"    Clash #{clash_count}: mobile res {mobile_res.GetNum()} atom {mobile_atom.GetIdx()} vs res {other_res.GetNum()} atom {other_atom.GetIdx()} (dist: {dist:.2f}Å < {vdw_sum*0.75:.2f}Å)")
        
        # Convert residue numbers to atom indices for those residues
        clashing_atom_indices = set()
        for atom in mobile_atoms:
            if atom.GetResidue() and atom.GetResidue().GetNum() in clashing_residues:
                clashing_atom_indices.add(atom.GetIdx())
        
        print(f"🎯 Clash-based energy optimization:")
        print(f"  - Found {clash_count} actual clashes")
        print(f"  - Clashing mobile residues: {len(clashing_residues)}")
        print(f"  - Mobile atoms in clashing residues: {len(clashing_atom_indices)}")
        print(f"  - Energy computed only for {len(clashing_atom_indices)}/{len(mobile_atoms)} mobile atoms")
        if len(clashing_atom_indices) > 0:
            speedup = len(atoms) / len(clashing_atom_indices)
            print(f"  - Expected speedup: ~{speedup:.1f}x")
        
        return clashing_atom_indices

    def get_atoms_from_clashing_residues(self, mol, mobile_atom_indices, clashing_residues):
        """
        Get atom indices for atoms that belong to specific clashing residues.
        This ensures Monte Carlo operates ONLY on residues that underwent sculpting.
        """
        if not mobile_atom_indices or not clashing_residues:
            return set()
            
        atoms = list(ob.OBMolAtomIter(mol))
        mobile_atoms = [atom for atom in atoms if atom.GetIdx() in mobile_atom_indices]
        
        # Extract residue info from clashing_residues set: (chain, resi, resn, obj_name)
        target_residues = set()
        for residue_info in clashing_residues:
            if len(residue_info) >= 2:  # (chain, resi, ...)
                target_residues.add((residue_info[0], residue_info[1]))  # (chain, resi)
        
        # Get atoms that belong to the target clashing residues
        target_atoms = set()
        for atom in mobile_atoms:
            atom_res = atom.GetResidue()
            if atom_res:
                atom_chain = atom_res.GetChain()
                atom_resi = str(atom_res.GetNum())
                if (atom_chain, atom_resi) in target_residues:
                    target_atoms.add(atom.GetIdx())
        
        print(f"🎯 Targeted Monte Carlo optimization:")
        print(f"  - Target clashing residues: {len(target_residues)}")
        print(f"  - Mobile atoms in target residues: {len(target_atoms)}")
        print(f"  - Optimization focus: {len(target_atoms)}/{len(mobile_atoms)} mobile atoms")
        if len(target_atoms) > 0:
            speedup = len(atoms) / len(target_atoms)
            print(f"  - Expected speedup: ~{speedup:.1f}x")
        
        return target_atoms

    # 'compute_localized_energy' stood here: it was never called from anywhere, and it
    # would not have worked if it had been. To leave the distant atoms out of the sum it
    # moved every one of them to (1000, 1000, 1000) -- all to the same point, on top of
    # each other and far from the atoms they are bonded to -- so the energy it read was
    # the energy of that, not of the molecule. It also read it without handing the moved
    # coordinates to the force field, which keeps its own copy.

    def _update_pymol_coords_from_ob(self, mol, selection, mobile_sele_name):
        """
        Update PyMOL coordinates from OpenBabel molecule.
        This is needed for saving trajectory frames.
        """
        try:
            from pymol import cmd
            
            # Get coordinates from OpenBabel molecule
            coords = []
            for atom in ob.OBMolAtomIter(mol):
                coords.append([atom.GetX(), atom.GetY(), atom.GetZ()])
            
            # Update PyMOL coordinates
            cmd.load_coords(coords, selection)
            
        except Exception as e:
            print(f"  ⚠️ Could not update PyMOL coordinates: {e}")

    def monte_carlo_torsion_minimize(self, mol, dihedrals, ff, max_steps, selection, mobile_sele_name, progress_dialog=None, mc_temp=300.0, target_object=None):
        """Monte Carlo torsion optimization with simulated annealing and fast energy calculation."""
        import random
        import math
        
        current_mol = ob.OBMol(mol)
        best_mol = ob.OBMol(mol)
        
        current_angles = [self.get_dihedral_angle(current_mol, d) for d in dihedrals]
        
        current_energy = self.compute_energy_safe(current_mol, ff)
        best_energy = current_energy
        
        unit = ff.GetUnit()
        print(f"Starting FAST Monte Carlo with {max_steps} steps, {len(dihedrals)} torsions...")
        print(f"Initial energy: {current_energy:.3f} {unit}")
        
        # Parameters for simulated annealing
        temp = mc_temp * 1.5  # Start higher for better exploration
        # Boltzmann's constant in the unit the force field answers in. The energies below
        # are its own, and a constant in another unit moves the temperature of the search
        # by the factor between them: with a kcal/mol constant and GAFF, which answers in
        # kJ/mol, the search ran four times colder than the temperature asked for.
        kB = 0.0019872 if "kcal" in unit.lower() else 0.0083145
        
        accepted_moves = 0
        save_interval = max(1, max_steps // 20)
        
        for step in range(max_steps):
            if not dihedrals: break
            torsion_idx = random.randrange(len(dihedrals))
            trial_angles = list(current_angles)
            
            delta = random.uniform(-5, 5)
            trial_angles[torsion_idx] = (trial_angles[torsion_idx] + delta) % 360.0
            
            trial_mol = ob.OBMol(current_mol)
            self.set_dihedral_angle(trial_mol, dihedrals[torsion_idx], trial_angles[torsion_idx])
            
            trial_energy = self.compute_energy_safe(trial_mol, ff)
            
            accept = False
            dE = trial_energy - current_energy
            
            if dE <= 0:
                accept = True
            else:
                # Metropolis, with nothing in front of it: the cut at 10 which used to
                # stand here refused every larger step whatever the temperature was, and
                # 10 of a unit the constant above was not in.
                exponent = -dE / (kB * temp) if temp > 0 else -float("inf")
                if exponent > -700 and random.random() < math.exp(exponent):
                    accept = True
            
            if accept:
                current_mol = trial_mol
                current_angles = trial_angles
                current_energy = trial_energy
                accepted_moves += 1
                
                if trial_energy < best_energy:
                    best_mol = ob.OBMol(trial_mol)
                    best_energy = trial_energy
            
            # Faster cooling
            temp *= 0.99
            
            if step % save_interval == 0 or step == max_steps - 1:
                acceptance_rate = (accepted_moves / (step + 1)) * 100
                print(f"Step {step+1}/{max_steps}: E = {current_energy:.3f}, Best = {best_energy:.3f}, T = {temp:.1f}, Accept = {acceptance_rate:.1f}%")
                if progress_dialog:
                    progress_dialog.update_progress(step + 1, current_energy)
        
        print(f"Optimization completed! Accepted {accepted_moves}/{max_steps} moves.")
        return best_mol, best_energy

    def _get_ligand_mobile_indices(self, selection, mobile_atom_indices):
        """Identifica gli atomi dei ligandi nella regione mobile."""
        ligand_indices = []
        
        try:
            # Ottieni solo gli atomi di ligandi (non proteine/acidi nucleici)
            ligand_selection = f"({selection}) and (organic) and not (polymer.protein or polymer.nucleic)"
            
            # Interseca con gli atomi mobili
            mobile_set = set(mobile_atom_indices)
            cmd.iterate(ligand_selection, "ligand_indices.append(index)", 
                        space={'ligand_indices': ligand_indices})
            
            # Ritorna solo gli atomi che sono sia ligandi che mobili
            final_ligand_indices = [idx for idx in ligand_indices if idx in mobile_set]
            
            if final_ligand_indices:
                print(f"   → Identified {len(final_ligand_indices)} ligand atoms for energy minimization")
            
            return final_ligand_indices
            
        except Exception as e:
            print(f"   ⚠ Error identifying ligand atoms: {e}")
            return []

    def analyze_clash_composition_and_strategy(self, selection, clashes):
        """Analyze clash composition and determine optimal optimization strategy."""
        try:
            # Get molecular composition
            protein_atoms = cmd.count_atoms(f"({selection}) and (polymer.protein)")
            nucleic_atoms = cmd.count_atoms(f"({selection}) and (polymer.nucleic)")
            ligand_atoms = cmd.count_atoms(f"({selection}) and (organic)")
            water_atoms = cmd.count_atoms(f"({selection}) and (resn HOH)")
            total_atoms = cmd.count_atoms(f"({selection})")
            
            if total_atoms == 0:
                return "torsion_only", False
            
            # Analyze clashes involving ligands
            ligand_clashes = 0
            total_clashes = len(clashes)
            
            for clash in clashes:
                a1_is_ligand = not (is_protein_residue(clash['a1'].resn) or is_nucleic_acid_residue(clash['a1'].resn))
                a2_is_ligand = not (is_protein_residue(clash['a2'].resn) or is_nucleic_acid_residue(clash['a2'].resn))
                
                # Skip water molecules
                if clash['a1'].resn.strip() == 'HOH' or clash['a2'].resn.strip() == 'HOH':
                    continue
                    
                if a1_is_ligand or a2_is_ligand:
                    ligand_clashes += 1
            
            has_ligand_clashes = ligand_clashes > 0
            ligand_clash_percentage = (ligand_clashes / max(1, total_clashes)) * 100
            
            # Determine strategy
            if has_ligand_clashes:
                strategy = "hybrid"  # Monte Carlo + full minimization for ligands
                print(f"\n🎯 OPTIMIZATION STRATEGY: Hybrid approach")
                print(f"   Monte Carlo torsion optimization + ligand energy minimization")
                print(f"   Reason: {ligand_clashes}/{total_clashes} clashes involve ligands ({ligand_clash_percentage:.1f}%)")
            else:
                strategy = "torsion_only"  # Monte Carlo only
                print(f"\n🎯 OPTIMIZATION STRATEGY: Monte Carlo torsion optimization")
                print(f"   Reason: Preserving ring geometries and bond lengths (no ligand clashes)")
            
            return strategy, has_ligand_clashes
            
        except Exception as e:
            print(f"Error in strategy analysis: {e}")
            return "torsion_only", False
    
    def run_full_minimize(self, selection, mobile_atom_indices, forcefield, method, nsteps, conv, cutoff=1, cut_vdw=6.0, cut_elec=8.0, original_elements=None):
        """Run full OpenBabel minimization (original implementation)."""
        from pymol import cmd
        mobile_sele_name = "mobile_residues_to_minimize"

        from pymol import stored
        stored.original_colors = {}
        progress_dialog = None
        exit_status = 0

        try:
            sele_string = f"({selection}) and index { '+'.join(map(str, mobile_atom_indices)) }"
            cmd.select(mobile_sele_name, sele_string)

            cmd.iterate(mobile_sele_name, "stored.original_colors[index] = color")

            print(f"Highlighting {len(mobile_atom_indices)} atoms for full minimization...")
            cmd.color("yellow", mobile_sele_name)
            cmd.show("spheres", mobile_sele_name)
            cmd.zoom(mobile_sele_name, buffer=3.0)
            cmd.refresh()

            import time
            time.sleep(1)

            print("\nStarting Open Babel full minimization...")
            print("CRITICAL: Ensuring molecule is properly protonated before minimization...")
            
            # Ensure hydrogens are present before export
            cmd.h_add(selection)
            
            mol_string = cmd.get_pdbstr(f"{selection} and state 1")
            cmd.delete(selection)
            cmd.read_pdbstr(mol_string, selection)

            obconversion = ob.OBConversion()
            obconversion.SetInAndOutFormats('pdb', 'pdb')
            mol = ob.OBMol()
            obconversion.ReadString(mol, mol_string)

            forcefield, ff = self._setup_force_field(forcefield, mol, mol_string,
                                                     obconversion)
            if ff is None:
                return 1

            # Create mapping between PyMOL indices and OpenBabel indices
            constraints = ob.OBFFConstraints()
            mobile_set = set(mobile_atom_indices)
            
            # Create a proper mapping by getting atom indices from the actual structure
            pymol_to_ob_mapping = {}
            
            # Get the mapping by iterating through the structure in the same order
            atom_counter = 1  # OpenBabel uses 1-based indexing
            for atom in ob.OBMolAtomIter(mol):
                pymol_to_ob_mapping[atom_counter] = atom.GetIdx()
                atom_counter += 1
            
            # Convert PyMOL mobile indices to OpenBabel indices
            mobile_ob_indices = set()
            for pymol_idx in mobile_atom_indices:
                if pymol_idx in pymol_to_ob_mapping:
                    mobile_ob_indices.add(pymol_to_ob_mapping[pymol_idx])
            
            print(f"DEBUG: Mobile atoms in PyMOL indices: {len(mobile_atom_indices)}")
            print(f"DEBUG: Mobile atoms in OpenBabel indices: {len(mobile_ob_indices)}")
            
            # Apply constraints to non-mobile atoms in OpenBabel
            constrained_atoms = 0
            for atom in ob.OBMolAtomIter(mol):
                if atom.GetIdx() not in mobile_ob_indices:
                    constraints.AddAtomConstraint(atom.GetIdx())
                    constrained_atoms += 1
            
            print(f"DEBUG: Applied constraints to {constrained_atoms} fixed atoms")
            ff.SetConstraints(constraints)
            
            # There used to be a "verification" of the constraints here: it moved a
            # constrained atom, read the energy, and concluded that the constraints were
            # weak if the energy had not risen by 1000. Open Babel holds an atom by
            # zeroing its gradient during the descent, not by adding a penalty to the
            # energy, so the figure it looked at never said anything about whether the
            # constraints were on, and the remedy it then applied was to set the same
            # constraints again. It cost three energies of the whole complex.

            if int(cutoff):
                ff.EnableCutOff(True)
                ff.SetVDWCutOff(float(cut_vdw))
                ff.SetElectrostaticCutOff(float(cut_elec))
                # As in molecule_optimizer: the pair list has to be built before the energy
                # below is read, or the figure is short of every interaction not yet in the
                # list. Here it matters twice, because that figure is both what the user is
                # shown and what the abort threshold a few lines down is compared against.
                ff.UpdatePairsSimple()

            initial_energy = self.compute_energy_safe(mol, ff)
            unit = ff.GetUnit()

            # What "much too high to be worth minimising" is depends on how many atoms
            # there are and on the unit the force field answers in. A flat 100000 was
            # neither: a protein complex under UFF, which answers in kJ/mol, starts in
            # the millions and is perfectly ordinary, so the command refused to run on
            # the structures it was written for, while a small molecule which really is
            # broken can sit below the figure. A thousand kcal/mol for every atom is
            # past anything a structure which can be minimised reaches.
            per_atom_limit = 1000.0 * (1.0 if "kcal" in unit.lower() else 4.184)
            energy_limit = per_atom_limit * max(mol.NumAtoms(), 1)

            if not math.isfinite(initial_energy) or initial_energy > energy_limit:
                self.pymod.main_window.show_error_message(
                    "High Initial Energy Warning",
                    f"Initial energy is extremely high ({initial_energy:.2f} {unit} for"
                    f" {mol.NumAtoms()} atoms).\n\n"
                    "This usually means something is wrong with the input structure:"
                    " atoms on top of each other, or a geometry a force field cannot"
                    " make sense of.\n\n"
                    "Minimization aborted. Please check your structure or try a different force field."
                )
                if mobile_sele_name in cmd.get_names("selections"):
                    cmd.delete(mobile_sele_name)
                exit_status = 1
                return exit_status

        except Exception as e:
            self.pymod.main_window.show_error_message("Setup Error", f"An error occurred during minimization setup: {e}")
            return 1

        try:
            if has_qt:
                # Parent to the PyMod main window so the dialog inherits PyMod's
                # stylesheet and is properly owned (otherwise it appears unstyled
                # and can misbehave with focus/the event loop).
                progress_dialog = MinimizationProgress(nsteps, parent=self.pymod.main_window,
                                                      unit=ff.GetUnit())
                progress_dialog.show()

            converged = False
            movie_state = 2 # Start movie from state 2
            for i in range(1, nsteps + 1):
                if method == 'Conjugate Gradients':
                    converged = ff.ConjugateGradients(1, conv)
                else:
                    converged = ff.SteepestDescent(1, conv)

                # Save intermediate states for a "movie" effect
                if (i > 0 and i % 20 == 0) or (i == nsteps) or (converged and i > 1):
                    current_view = cmd.get_view()
                    ff.GetCoordinates(mol)
                    pdb_out = obconversion.WriteString(mol)

                    # Launder the new state through a temporary object
                    temp_state_obj = "temp_local_min_state"
                    cmd.delete(temp_state_obj)
                    if pdb_out.strip():
                        cmd.read_pdbstr(pdb_out, temp_state_obj)
                        cmd.create(selection, temp_state_obj, source_state=1, target_state=movie_state)
                        cmd.delete(temp_state_obj)

                    movie_state += 1
                    cmd.color("yellow", mobile_sele_name)
                    cmd.show("spheres", mobile_sele_name)
                    cmd.set_view(current_view)
                    if has_qt:
                        QtWidgets.QApplication.processEvents()

                energy = self.compute_energy_safe(mol, ff)
                if i % 10 == 0 or converged:
                    print(f"Step {i}/{nsteps}... Energy: {energy:.2f} {ff.GetUnit()}")
                    if progress_dialog:
                        progress_dialog.update_progress(i, energy)

                if converged:
                    print("Minimization converged.")
                    break

            final_energy = self.compute_energy_safe(mol, ff)
            self._apply_trajectory_coloring(selection, mobile_atom_indices)
            
            self._record_energies("Full minimization", initial_energy, final_energy,
                                  ff.GetUnit())

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.pymod.main_window.show_error_message("Open Babel Error", f"Full minimization failed: {e}")

        finally:
            if has_qt and progress_dialog:
                progress_dialog.close()
            if mobile_sele_name in cmd.get_names("selections"):
                print("Cleaning up temporary selections and colors...")
                cmd.hide("spheres", mobile_sele_name)
                cmd.delete(mobile_sele_name)
                if 'original_colors' in stored.__dict__:
                    del stored.original_colors
            
            if "temp_mobile_preview" in cmd.get_names("selections"):
                cmd.delete("temp_mobile_preview")
                print("Cleaned up mobile preview selection")

        return exit_status

    def _apply_trajectory_coloring_enhanced(self, trajectory_name, clashes, clashing_residues_info, interacting_objects, mobile_radius):
        """
        Apply intelligent coloring to the trajectory object:
        - RED: Residues/nucleotides directly involved in clashes
        - GREEN: Mobile area around clashes (minimized region)
        - GRAY: Fixed atoms (constrained during minimization)
        """
        try:
            print(f"\n=== ENHANCED TRAJECTORY COLORING ===")
            print(f"Applying intelligent coloring to '{trajectory_name}'...")
            
            # Define custom colors
            cmd.set_color("clash_red", [0.9, 0.1, 0.1])        # Bright red for clashes
            cmd.set_color("mobile_green", [0.1, 0.8, 0.1])     # Bright green for mobile
            cmd.set_color("fixed_gray", [0.6, 0.6, 0.6])       # Gray for fixed
            
            # Step 1: Color everything gray (fixed) by default
            cmd.color("fixed_gray", trajectory_name)
            print(f"  ✓ Base coloring: Gray (fixed atoms)")
            
            # Step 2: Identify and color clash residues RED
            clash_residues = set()
            for clash in clashes:
                # Aggiungi residui coinvolti nei clash
                a1, a2 = clash['a1'], clash['a2']
                clash_residues.add((a1.chain, a1.resi))
                clash_residues.add((a2.chain, a2.resi))
            
            if clash_residues:
                clash_selection_parts = []
                for chain, resi in clash_residues:
                    clash_selection_parts.append(f"(chain {chain} and resi {resi})")
                
                clash_selection = f"({trajectory_name}) and ({' or '.join(clash_selection_parts)})"
                cmd.color("clash_red", clash_selection)
                print(f"  ✓ Clash residues colored RED: {len(clash_residues)} residues")
                
                # Step 3: Identify and color mobile area GREEN
                if mobile_radius > 0:
                    # Crea selezione mobile attorno ai clash
                    core_clash_selection = " or ".join(clash_selection_parts)
                    mobile_selection = f"({trajectory_name}) and (({core_clash_selection}) around {mobile_radius})"
                    
                    # Rimuovi i residui clash dalla selezione mobile (per non sovrascrivere il rosso)
                    mobile_selection_filtered = f"({mobile_selection}) and not ({core_clash_selection})"
                    
                    try:
                        mobile_count = cmd.count_atoms(mobile_selection_filtered)
                        if mobile_count > 0:
                            cmd.color("mobile_green", mobile_selection_filtered)
                            print(f"  ✓ Mobile area colored GREEN: {mobile_count} atoms within {mobile_radius}Å")
                        else:
                            print(f"  ⚠ No mobile atoms found around clashes")
                    except Exception as e:
                        print(f"  ⚠ Mobile area coloring failed: {e}")
            
            # Step 4: Enhanced visual representation
            try:
                # Show clashing residues as sticks (more visible)
                if clash_residues:
                    clash_stick_selection = f"({trajectory_name}) and ({' or '.join([f'(chain {c} and resi {r})' for c, r in clash_residues])})"
                    cmd.show("sticks", clash_stick_selection)
                    cmd.set("stick_radius", 0.2, clash_stick_selection)
                    print(f"  ✓ Clash residues shown as thick sticks")
                
                # Show mobile area as thinner sticks
                if mobile_radius > 0:
                    mobile_stick_selection = f"({trajectory_name}) and (({' or '.join([f'(chain {c} and resi {r})' for c, r in clash_residues])}) around {mobile_radius}) and not ({' or '.join([f'(chain {c} and resi {r})' for c, r in clash_residues])})"
                    try:
                        mobile_stick_count = cmd.count_atoms(mobile_stick_selection)
                        if mobile_stick_count > 0:
                            cmd.show("sticks", mobile_stick_selection)
                            cmd.set("stick_radius", 0.1, mobile_stick_selection)
                            print(f"  ✓ Mobile area shown as thin sticks")
                    except:
                        pass
                
            except Exception as e:
                print(f"  ⚠ Enhanced visualization failed: {e}")
            
            # Step 5: Create informative selections for user
            try:
                # Crea selezioni permanenti per facile accesso
                if clash_residues:
                    clash_sel_name = "clash_residues"
                    clash_sel_string = f"({trajectory_name}) and ({' or '.join([f'(chain {c} and resi {r})' for c, r in clash_residues])})"
                    cmd.select(clash_sel_name, clash_sel_string)
                    print(f"  ✓ Created selection '{clash_sel_name}' for clash residues")
                
                if mobile_radius > 0:
                    mobile_sel_name = "mobile_area"
                    mobile_sel_string = f"({trajectory_name}) and (({' or '.join([f'(chain {c} and resi {r})' for c, r in clash_residues])}) around {mobile_radius})"
                    cmd.select(mobile_sel_name, mobile_sel_string)
                    print(f"  ✓ Created selection '{mobile_sel_name}' for mobile area")
                
            except Exception as e:
                print(f"  ⚠ Selection creation failed: {e}")
            
            # Step 6: Summary report
            total_atoms = cmd.count_atoms(trajectory_name)
            clash_atoms = sum(cmd.count_atoms(f"({trajectory_name}) and chain {c} and resi {r}") for c, r in clash_residues)
            
            try:
                mobile_atoms = cmd.count_atoms(mobile_selection) if mobile_radius > 0 else 0
                fixed_atoms = total_atoms - mobile_atoms
            except:
                mobile_atoms = 0
                fixed_atoms = total_atoms - clash_atoms
            
            print(f"\n  📊 TRAJECTORY COLORING SUMMARY:")
            print(f"     🔴 RED (clash residues): {clash_atoms} atoms")
            print(f"     🟢 GREEN (mobile area): {mobile_atoms - clash_atoms} atoms")
            print(f"     ⚫ GRAY (fixed region): {fixed_atoms} atoms")
            print(f"     📊 Total trajectory: {total_atoms} atoms")
            print(f"     📏 Mobile radius: {mobile_radius} Å")
            
            print(f"  ✅ Enhanced trajectory coloring complete!")
            
        except Exception as e:
            print(f"❌ Enhanced trajectory coloring failed: {e}")
            # Fallback to simple coloring
            try:
                cmd.color("gray", trajectory_name)
                print(f"  ✓ Applied fallback gray coloring")
            except:
                pass

    def _apply_trajectory_coloring(self, selection_name, mobile_atom_indices):
        """Apply special coloring to the trajectory object: light gray for all atoms, green for mobile atoms."""
        try:
            print("Applying trajectory coloring scheme...")
            
            # Define colors
            light_gray_color = "gray80"
            mobile_color = "green"
            
            cmd.set_color("light_gray_trajectory", [0.8, 0.8, 0.8])
            
            # First, color everything light gray
            cmd.color("light_gray_trajectory", selection_name)
            
            # Then, process mobile atoms in chunks to avoid long selection strings
            if mobile_atom_indices:
                chunk_size = 100  # Adjust chunk size if needed
                for i in range(0, len(mobile_atom_indices), chunk_size):
                    chunk = mobile_atom_indices[i:i + chunk_size]
                    indices_str = "+".join(map(str, chunk))
                    mobile_selection_chunk = f"({selection_name}) and index {indices_str}"
                    
                    # Apply coloring and visualization to the chunk
                    cmd.color(mobile_color, mobile_selection_chunk)
                    cmd.show("sticks", mobile_selection_chunk)
                    cmd.set("stick_radius", 0.15, mobile_selection_chunk)

                print(f"✅ Trajectory coloring applied:")
                print(f"   - All atoms: Light gray")
                print(f"   - Mobile atoms ({len(mobile_atom_indices)}): Green")
                
        except Exception as e:
            print(f"Warning: Could not apply trajectory coloring: {e}")

    def _record_energies(self, what, before, after, unit):
        """
        Keep what the minimisation did, for the window at the end of the protocol.

        Each minimisation used to open a window of its own as it finished, and the
        protocol opened another one a moment later: two windows for one run, the first
        of them holding the pair of figures this is about.
        """

        self.energy_report = {"what": what, "before": before, "after": after,
                              "unit": unit}
        print(f"{what} complete. {before:.2f} -> {after:.2f} {unit}"
              f" ({after - before:+.2f})")

    def _energy_summary(self):
        """
        What the minimisation did, as lines for the window at the end, or nothing when
        no minimisation was run.

        The change is the part which says something. The totals are the energy of the
        whole complex, nearly all of which is held fixed and contributes the same figure
        at both ends, so they come out large and nearly equal however well the
        refinement went: five million kJ/mol at both ends with a few thousand between
        them is an ordinary result, not a failed run.
        """

        report = getattr(self, "energy_report", None)
        if not report:
            return ""
        return ("\n%s: %+.2f %s\n"
                "(whole complex, %.2f → %.2f)\n"
                % (report["what"], report["after"] - report["before"], report["unit"],
                   report["before"], report["after"]))

    def _remove_original_elements(self, interacting_objects):
        """
        Delete the structures the refinement started from, from PyMOL and from PyMod.

        An element is matched to an object by its PyMOL selector and by nothing else.
        It used to be enough for the name of the object to appear anywhere in the text
        of an element, so refining '1ol5_chain_A' could take '1ol5_chain_A_refined', or
        any other element whose header happened to contain those letters, down with it.
        """

        print("  Removing the original objects from PyMod and PyMOL...")

        for obj_name in interacting_objects:
            if obj_name in cmd.get_names("objects"):
                cmd.delete(obj_name)
                print(f"    ✓ Removed '{obj_name}' from PyMOL")

        wanted = set(interacting_objects)
        elements = list(self.selected_elements)
        if not elements:
            elements = self.pymod.get_pymod_elements_list()

        removed_count = 0
        for element in elements:
            if not getattr(element, "structure", None):
                continue
            if element.get_pymol_selector() not in wanted:
                continue
            try:
                element.delete()
                removed_count += 1
                print(f"    ✓ Removed PyMod element: {element}")
            except Exception as e:
                print(f"    ❌ Error removing PyMod element {element}: {e}")

        if removed_count:
            self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                           update_menus=True)
        else:
            print("    ⚠ No PyMod element matched the objects which were refined")


    def _cleanup_objects(self, object_list):
        """
        Improved cleanup method that handles both objects and selections robustly.
        """
        if not object_list:
            return
            
        print(f"Cleaning up: {object_list}")
        
        for obj_name in object_list:
            try:
                # Try to delete as object first
                if obj_name in cmd.get_names("objects"):
                    cmd.delete(obj_name)
                    print(f"  ✓ Deleted object: {obj_name}")
                
                # Try to delete as selection
                elif obj_name in cmd.get_names("selections"):
                    cmd.delete(obj_name)
                    print(f"  ✓ Deleted selection: {obj_name}")
                    
            except Exception as e:
                print(f"  ❌ Error deleting {obj_name}: {e}")

        # A sweep by name stood here: anything whose name ended in '_preview' or
        # '_optimized_final' was deleted, whoever had made it. The one object of this
        # protocol which those patterns would have caught is 'temp_mobile_preview',
        # which is named outright in the lists passed to this method, and nothing here
        # has ever been called '_optimized_final'. What was left was a sweep over the
        # user's own objects.

    def _finalize_visualization(self, target_object, clashing_residues_info):
        """
        Apply final visualization settings to the result object.
        """
        print("\n--- Finalizing visualization ---")

        # 1. Remove hydrogens
        cmd.remove(f"({target_object}) and hydro")
        print(f"  ✓ Removed hydrogens from '{target_object}'")

        # 2. Show as cartoon and sticks for side chains
        cmd.hide("everything", target_object)
        cmd.show("cartoon", target_object)
        cmd.show("sticks", f"({target_object}) and not (name C+N+O)")
        print(f"  ✓ Set representation to cartoon and sticks for '{target_object}'")

        # 3. Coloring
        cmd.color("gray80", target_object)
        cmd.set("cartoon_color", "gray80", target_object)
        cmd.set("stick_color", "gray80", target_object)
        print(f"  ✓ Colored '{target_object}' gray")

        if clashing_residues_info:
            clash_sele = " or ".join([f"(chain {c} and resi {r})" for c, r, n, obj in clashing_residues_info])
            num_states = cmd.count_states(target_object)
            print(f"  ✓ Coloring '{num_states}' states")

            # State 1: Red
            if num_states >= 1:
                sele = f"({target_object} and ({clash_sele})) and state 1"
                cmd.set("cartoon_color", "red", sele)
                cmd.set("stick_color", "red", sele)
                print("  ✓ Colored clashing residues red in state 1")

            # State 2: Orange
            if num_states >= 2:
                sele = f"({target_object} and ({clash_sele})) and state 2"
                cmd.set("cartoon_color", "orange", sele)
                cmd.set("stick_color", "orange", sele)
                print("  ✓ Colored clashing residues orange in state 2")

            # Final state: Green
            if num_states > 2:
                last_state = num_states
                sele = f"({target_object} and ({clash_sele})) and state {last_state}"
                cmd.set("cartoon_color", "green", sele)
                cmd.set("stick_color", "green", sele)
                print(f"  ✓ Colored clashing residues green in state {last_state}")

        # 4. Set final state
        last_state = cmd.count_states(target_object)
        if last_state > 1:
            cmd.frame(last_state)
            print(f"  ✓ Set view to final state ({last_state})")

def clash_and_refine_protocol(pymod_instance, selected_elements=None):
    """Entry point function for the protocol."""
    protocol = ClashAndRefine(pymod_instance, selected_elements)
    protocol.run_clash_refine()
