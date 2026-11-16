# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Key-gated ESM / Forge modeling protocols for PyMod.

Two protocols are provided, both under the Modeling -> ESMFold menu:

* ``ESM_forge_fold`` -- sequence -> structure, via the Forge ``/fold`` endpoint
  (an ESMFold2-class model). This is offered *alongside* the free, keyless
  ESMAtlas "Build model" protocol, not as a replacement.

* ``ESM_inverse_fold`` -- structure -> sequence (protein design), via the Forge
  ``/inverse_fold`` endpoint (an ESM3-class model). This is a new capability.

Both require an API key, configured by the user in the PyMod Options window
(the "ESM service" tool parameters). When no key is set, the protocols
show an explanatory message instead of running. Networking uses ``requests``
(already a PyMod dependency) through the vendored, SDK-free ``ForgeClient``.
"""

import os

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_combobox_qt)
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import PyModInterruptedProtocol
import pymod_lib.pymod_structure as pmstr
import pymod_lib.pymod_vars as pmdt

from .forge_client import ForgeClient, DEFAULT_FORGE_BASE_URL
from . import atom37


# Model identifiers exposed in the GUI. These are the defaults from the
# prototype scripts; they are editable here in one place if the service renames
# them.
FOLD_MODELS = ["esmfold2-fast-2026-05"]
INVERSE_FOLD_MODELS = ["esm3-open-2024-03"]


# ---------------------------------------------------------------------------
# Shared helpers.
# ---------------------------------------------------------------------------

class _Forge_protocol_mixin:
    """
    Common credential handling and error reporting for the Forge protocols.
    """

    def get_forge_credentials(self):
        """
        Returns ``(api_key, base_url)`` from the PyMod preferences, or
        ``(None, base_url)`` if no key is configured.
        """
        try:
            api_key = self.pymod.esm["forge_api_key"].get_value()
            base_url = self.pymod.esm["forge_base_url"].get_value()
        except Exception:
            api_key, base_url = None, DEFAULT_FORGE_BASE_URL
        api_key = (api_key or "").strip() or None
        base_url = (base_url or "").strip() or DEFAULT_FORGE_BASE_URL
        return api_key, base_url

    def require_api_key(self):
        """
        Returns a ready ``ForgeClient`` if a key is set, otherwise shows a
        message pointing the user to the Options window and returns ``None``.
        """
        api_key, base_url = self.get_forge_credentials()
        if not api_key:
            self.pymod.main_window.show_warning_message(
                "ESM API key required",
                "This feature uses the ESM REST service of EvolutionaryScale, "
                "which needs an API key.\n\n"
                "Set your key in PyMod under Tools → Options → "
                "'ESM service' → ESM API Key, then try again.\n\n"
                "The free, keyless ESMFold 'Build model' option does not require "
                "a key.")
            return None
        return ForgeClient(api_key=api_key, base_url=base_url)

    def show_api_error(self, title, exc):
        # Always log the full error to the console for debugging, but show a
        # length-capped message in the dialog: an unbounded string (e.g. a
        # server response echoing the request) would stretch the message box
        # across the whole screen.
        full = str(exc)
        print("[%s] %s" % (title, full))
        message = full if len(full) <= 600 else (full[:600] + " […]")
        self.pymod.main_window.show_error_message(title, message)

    def call_forge_with_progress(self, client, endpoint, payload,
                                 title, label_text, error_title):
        """
        Perform a (slow) Forge REST call on a background thread while showing a
        modal progress dialog, so the PyMod GUI stays responsive.

        Only the network request runs off the main thread; the caller does all
        Qt/element work after this returns (Qt objects must be touched from the
        main thread only). Returns the decoded response on success, or ``None``
        if the user cancelled or the call failed (the error dialog is shown
        here in the latter case).
        """
        # The worker thread cannot return a value, so it stashes the result on
        # the protocol instance.
        self._forge_response = None

        def _worker():
            self._forge_response = client.post_json(endpoint, payload)

        p_dialog = Protocol_exec_dialog(
            app=self.pymod.main_window, pymod=self.pymod,
            function=_worker, args=(),
            title=title, label_text=label_text,
            wait_start=0.15, wait_end=0.15)
        try:
            p_dialog.exec_()
        except PyModInterruptedProtocol:
            # The user pressed Cancel / closed the dialog.
            return None
        except Exception as e:
            self.show_api_error(error_title, e)
            return None

        return self._forge_response


# ---------------------------------------------------------------------------
# Forge fold: sequence -> structure.
# ---------------------------------------------------------------------------

def rescale_plddt(values):
    """
    The per residue confidence of a prediction, on the 0-100 scale of a B-factor.

    The service answers on a scale of 0 to 1 and PyMOL wants the conventional
    one, so the numbers are multiplied by a hundred -- unless they are already on
    that scale, which the service also does.

    A complex has a gap in this array: the chains are sent joined by ':', and the
    separator has no confidence, so its entry is null. Those nulls are kept where
    they are, because the array has to stay in step with the sequence which was
    sent, and they are stepped over when the largest value is looked for. Asking
    for the largest without stepping over them ended a prediction of two chains
    in "'>' not supported between instances of 'NoneType' and 'float'".
    """

    if not values:
        return values
    known = [value for value in values if value is not None]
    if not known or max(known) > 1.0:
        return values
    return [None if value is None else value * 100.0 for value in values]


def mean_of(values):
    """
    The mean of the numbers there are, or None when there are none. The nulls
    which separate the chains of a complex take no part in it: counting them
    would drag the mean of every complex down.
    """

    known = [value for value in (values or []) if value is not None]
    if not known:
        return None
    return sum(known) / float(len(known))


class ESM_forge_fold(_Forge_protocol_mixin, PyMod_protocol):
    """
    Fold one or more selected protein sequences through the Forge ``/fold``
    endpoint and load the predicted structure into PyMod.

    With a single selected sequence a monomer is folded. With several selected
    sequences they are folded together as one complex (multimer): the chains
    are joined with ':' for the request, and the predicted complex is loaded
    back as separate chains in PyMOL and separate elements in PyMod.
    """

    protocol_name = "ESMFold2 prediction"

    def launch_from_gui(self):
        selected_sequences = self.pymod.get_selected_sequences()

        if len(selected_sequences) == 0:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Please select one or more target protein sequences to fold "
                "(select several to fold them together as a complex).")
            return None

        for element in selected_sequences:
            if element.polymer_type in ("dna", "rna"):
                self.pymod.main_window.show_error_message(
                    "Selection Error",
                    "Structure prediction can only be used on protein sequences. ESM does "
                    "not model DNA or RNA, so nucleic-acid chains can not be "
                    "included in the complex.")
                return None
            if not element.can_be_modeled():
                self.pymod.main_window.show_error_message(
                    "Selection Error",
                    "Please select only sequences that do not already have a "
                    "structure loaded in PyMOL.")
                return None

        # Key check before opening the options window.
        if self.require_api_key() is None:
            return None

        self.target_elements = list(selected_sequences)
        # Kept for backward compatibility with code expecting a single element.
        self.target_element = self.target_elements[0]
        self.options_window = Forge_fold_options_window_qt(
            parent=self.pymod.main_window, protocol=self,
            title="Structure prediction with ESMFold2",
            upper_frame_title="Predict the structure of one or more protein chains with ESMFold2",
            submit_command=self.run_from_options_window)
        self.options_window.show()

    def run_from_options_window(self):
        model = self.options_window.model_combobox.get()
        try:
            num_steps = int(self.options_window.steps_entry.getvalue(validate=False))
        except (ValueError, TypeError):
            num_steps = 32
        num_steps = max(1, min(256, num_steps))
        try:
            num_loops = int(self.options_window.loops_entry.getvalue(validate=False))
        except (ValueError, TypeError):
            num_loops = 3
        num_loops = max(1, min(12, num_loops))
        self.options_window.close()
        self.run(model=model, num_sampling_steps=num_steps, num_loops=num_loops)

    def run(self, model, num_sampling_steps=32, num_loops=3):
        client = self.require_api_key()
        if client is None:
            return None

        # One sequence per selected chain. For a complex the chains are joined
        # with ':' (the ESMFold2 multimer convention): the service returns one
        # extra, all-null "separator" coordinate row per ':'.
        chain_sequences = []
        for element in self.target_elements:
            chain_seq = "".join(c for c in element.my_sequence.upper()
                                if c.isalpha())
            if chain_seq:
                chain_sequences.append(chain_seq)

        if not chain_sequences:
            self.pymod.main_window.show_error_message(
                "Structure prediction failed",
                "The selected sequence(s) contain no residues to fold.")
            return None

        is_complex = len(chain_sequences) > 1
        joined_sequence = ":".join(chain_sequences)

        payload = {
            "model": model,
            "potential_sequence_of_concern": False,
            "sequence": joined_sequence,
            "msa": None,
            "num_loops": num_loops,
            "num_sampling_steps": num_sampling_steps,
            # Always request PAE so the model offers a "Predicted Aligned Error"
            # plot from the Models menu, like AlphaFold models do.
            "include_pae": True,
            "include_distogram": False,
            # Inter-chain ipTM is only meaningful for a complex.
            "include_pair_chains_iptm": is_complex,
            "include_embeddings": False,
        }

        what = "complex" if is_complex else "sequence"
        result = self.call_forge_with_progress(
            client, "fold", payload,
            title="Structure prediction with ESMFold2",
            label_text=("Predicting the structure of the %s with ESMFold2.\n"
                        "This can take a while depending on the sequence length "
                        "and the service load. Please wait..." % what),
            error_title="Structure prediction failed")
        if result is None:
            return None

        coords = result.get("coordinates") if isinstance(result, dict) else None
        if coords is None:
            self.pymod.main_window.show_error_message(
                "Structure prediction failed",
                "The service response did not contain structure coordinates.\n"
                "Response keys: %s"
                % (sorted(result.keys()) if isinstance(result, dict) else type(result)))
            return None

        # The service returns pLDDT on a 0-1 scale; PyMOL/structure tools expect
        # the conventional 0-100 B-factor scale, so rescale when needed.
        plddt = rescale_plddt(result.get("plddt"))

        # Reconstruct a PDB from the atom37 coordinates, then load it with the
        # same machinery PyMod uses for any model file.
        models_dir = os.path.join(self.pymod.current_project_dirpath,
                                  self.pymod.models_dirname)
        if not os.path.isdir(models_dir):
            os.mkdir(models_dir)

        if not hasattr(self.pymod, "performed_esm_modeling_count"):
            self.pymod.performed_esm_modeling_count = 0
        self.pymod.performed_esm_modeling_count += 1
        model_name = "esmfold2_model_%d" % self.pymod.performed_esm_modeling_count
        pdb_path = os.path.join(models_dir, model_name + ".pdb")

        # Split the (possibly multi-chain) response into per-chain records. The
        # 'coordinates'/'plddt' arrays align 1:1 with the joined sequence
        # *including* the ':' characters, which map to all-null separator rows;
        # we consume and skip those, starting a new chain at each ':'.
        chain_records = []
        cur_seq, cur_coords, cur_b = [], [], []
        coord_index = 0
        coords_ok = True
        for ch in joined_sequence:
            if ch == ":":
                if cur_seq:
                    chain_records.append(("".join(cur_seq), cur_coords, cur_b))
                cur_seq, cur_coords, cur_b = [], [], []
                coord_index += 1  # consume (skip) the separator coordinate row
                continue
            if coord_index >= len(coords):
                coords_ok = False
                break
            cur_seq.append(ch)
            cur_coords.append(coords[coord_index])
            confidence = (plddt[coord_index]
                          if plddt and coord_index < len(plddt) else None)
            cur_b.append(0.0 if confidence is None else confidence)
            coord_index += 1
        if cur_seq:
            chain_records.append(("".join(cur_seq), cur_coords, cur_b))

        if not coords_ok or not chain_records:
            self.pymod.main_window.show_error_message(
                "Structure prediction failed",
                "The structure coordinates returned by the service did not match "
                "the submitted sequence(s).")
            return None

        try:
            if len(chain_records) > 1:
                atom37.write_multichain_pdb_from_atom37(chain_records, pdb_path)
            else:
                seq0, coords0, b0 = chain_records[0]
                atom37.write_pdb_from_atom37(seq0, coords0, pdb_path, b_factors=b0)
        except Exception as e:
            self.show_api_error("Could not build the PDB from the response of the service", e)
            return None

        parsed_model_file = pmstr.Parsed_model_pdb_file(
            self.pymod, pdb_path,
            output_directory=self.pymod.structures_dirpath,
            new_file_name=model_name, model_root_name=model_name)

        colors = pmdt.pymol_light_colors_list
        loaded_elements = []
        for chain_number, model_element in enumerate(parsed_model_file.get_pymod_elements()):
            self.pymod.add_element_to_pymod(
                model_element, color=colors[chain_number % len(colors)])
            loaded_elements.append(model_element)

        mean_plddt = mean_of(plddt)
        if mean_plddt is not None:
            print("Mean pLDDT of the predicted structure: %.2f" % mean_plddt)

        # Register this fold as an ESM modeling session so it shows up under
        # Models -> ESMFold (with its pLDDT plot), like the keyless ESMFold and
        # AlphaFold models. The pLDDT plot reads per-residue values from the CA
        # B-factors of 'full_output_path', which were written above.
        self.model_name = model_name
        self.full_output_path = pdb_path
        self.new_model_element = loaded_elements[0] if loaded_elements else None
        # All chain elements of the (possibly multi-chain) complex, in chain
        # order, so the PAE plot can map the N x N matrix across every chain.
        self.new_model_elements = loaded_elements
        # Store the Predicted Aligned Error matrix (N x N), if the service
        # returned it, so the Models menu can offer a PAE plot for this model.
        self.pae_matrix = None
        if isinstance(result, dict):
            for pae_key in ("pae", "predicted_aligned_error"):
                pae_values = result.get(pae_key)
                if pae_values:
                    self.pae_matrix = pae_values
                    break
        if not hasattr(self.pymod, "esm_modeling_session_list"):
            self.pymod.esm_modeling_session_list = []
        self.pymod.esm_modeling_session_list.append(self)

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                       update_menus=True, update_elements=True)
        return pdb_path


class Forge_fold_options_window_qt(PyMod_protocol_window_qt):

    def build_protocol_middle_frame(self):
        self.model_combobox = PyMod_combobox_qt(
            label_text="Model", items=FOLD_MODELS)
        self.model_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.model_combobox)

        self.steps_entry = PyMod_entryfield_qt(
            label_text="Sampling steps", value="32")
        self.middle_formlayout.add_widget_to_align(self.steps_entry)

        self.loops_entry = PyMod_entryfield_qt(
            label_text="Recycles (num_loops)", value="3")
        self.middle_formlayout.add_widget_to_align(self.loops_entry)

        self.middle_formlayout.set_input_widgets_width("auto")


# ---------------------------------------------------------------------------
# Inverse folding: structure -> designed sequence.
# ---------------------------------------------------------------------------

class ESM_inverse_fold(_Forge_protocol_mixin, PyMod_protocol):
    """
    Design a sequence for a selected 3D structure through the Forge
    ``/inverse_fold`` endpoint and load the designed sequence into PyMod.
    """

    protocol_name = "ESM3 sequence design"

    def launch_from_gui(self):
        selected = self.pymod.get_selected_elements()
        structures = [e for e in selected
                      if e.has_structure() and e.polymer_type not in ("dna", "rna")]

        if len(structures) == 0:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Please select one protein structure (an element with a 3D "
                "structure loaded in PyMOL) to design a sequence for.")
            return None
        if len(structures) > 1:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Please select only one structure. Designing a sequence works on a "
                "single chain at a time.")
            return None

        if self.require_api_key() is None:
            return None

        self.target_element = structures[0]
        self.options_window = Inverse_fold_options_window_qt(
            parent=self.pymod.main_window, protocol=self,
            title="Design a sequence with ESM3",
            upper_frame_title="Design a sequence which folds into the selected structure",
            submit_command=self.run_from_options_window)
        self.options_window.show()

    def run_from_options_window(self):
        model = self.options_window.model_combobox.get()
        try:
            temperature = float(self.options_window.temp_entry.getvalue(validate=False))
        except (ValueError, TypeError):
            temperature = 0.1
        temperature = max(0.0, min(2.0, temperature))
        try:
            num_steps = int(self.options_window.steps_entry.getvalue(validate=False))
        except (ValueError, TypeError):
            num_steps = 8
        num_steps = max(1, num_steps)
        self.options_window.close()
        self.run(model=model, temperature=temperature, num_steps=num_steps)

    def run(self, model, temperature=0.1, num_steps=8):
        client = self.require_api_key()
        if client is None:
            return None

        # The element's per-chain PDB file is the cleanest input: a single
        # chain already split out by PyMod.
        try:
            pdb_path = self.target_element.get_structure_file(basename_only=False)
        except Exception as e:
            self.show_api_error("Could not locate the structure file", e)
            return None
        chain_id = None
        try:
            chain_id = self.target_element.get_chain_id()
        except Exception:
            chain_id = None

        try:
            seq_from_pdb, coordinates, _residues = atom37.pdb_to_atom37(
                pdb_path, chain_id=chain_id)
        except Exception as e:
            self.show_api_error("Could not read coordinates from the structure", e)
            return None

        if not coordinates:
            self.pymod.main_window.show_error_message(
                "Sequence design failed",
                "No standard protein residues with backbone coordinates were "
                "found in the selected structure.")
            return None

        # Inverse folding is performed through the generic '/generate' endpoint,
        # decoding the *sequence* track from the (fixed) structure track. The
        # input sequence is fully masked ('_'), and the number of decoding steps
        # must not exceed the sequence length.
        n = len(coordinates)
        num_steps = max(1, min(n, num_steps))
        payload = {
            "model": model,
            "track": "sequence",
            "num_steps": num_steps,
            "temperature": temperature,
            "potential_sequence_of_concern": False,
            "inputs": {
                "sequence": "_" * n,
                # Missing atoms must be sent as [null, null, null] (a valid
                # list), not a bare null, which the API rejects.
                "coordinates": atom37.coordinates_with_null_atoms(coordinates),
            },
        }

        result = self.call_forge_with_progress(
            client, "generate", payload,
            title="Design a sequence with ESM3",
            label_text=("Designing a sequence for the selected structure with ESM3 "
                        "the ESM service.\nThis can take a while. Please "
                        "wait..."),
            error_title="Sequence design failed")
        if result is None:
            return None

        designed = _extract_sequence_from_response(result)
        if not designed:
            self.pymod.main_window.show_error_message(
                "Sequence design failed",
                "The service response did not contain a designed sequence.\n"
                "Response keys: %s"
                % (sorted(result.keys()) if isinstance(result, dict) else type(result)))
            return None

        # Load the designed sequence as a new PyMod element.
        if not hasattr(self.pymod, "performed_esm_inverse_fold_count"):
            self.pymod.performed_esm_inverse_fold_count = 0
        self.pymod.performed_esm_inverse_fold_count += 1
        base_name = self.target_element.get_structure_file(
            basename_only=True, strip_extension=True) or "design"
        design_name = "%s_design_%d" % (base_name,
                                        self.pymod.performed_esm_inverse_fold_count)

        new_element = self.pymod.build_pymod_element_from_args(design_name, designed)
        self.pymod.add_element_to_pymod(new_element)
        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                       update_menus=True, update_elements=True)

        print("Inverse-folded design (%d aa): %s" % (len(designed), designed))
        return designed


class Inverse_fold_options_window_qt(PyMod_protocol_window_qt):

    def build_protocol_middle_frame(self):
        self.model_combobox = PyMod_combobox_qt(
            label_text="Model", items=INVERSE_FOLD_MODELS)
        self.model_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.model_combobox)

        self.temp_entry = PyMod_entryfield_qt(
            label_text="Temperature", value="0.1")
        self.middle_formlayout.add_widget_to_align(self.temp_entry)

        self.steps_entry = PyMod_entryfield_qt(
            label_text="Decoding steps (num_steps)", value="8")
        self.middle_formlayout.add_widget_to_align(self.steps_entry)

        self.middle_formlayout.set_input_widgets_width("auto")


def _extract_sequence_from_response(result):
    """
    Defensive extraction of the designed sequence. The '/generate' endpoint
    returns it under ``outputs.sequence``; older/alternative shapes are also
    handled in case the response evolves.
    """
    if isinstance(result, dict):
        for key in ("sequence", "designed_sequence"):
            if isinstance(result.get(key), str):
                return result[key]
        for key in ("outputs", "output", "protein"):
            obj = result.get(key)
            if isinstance(obj, dict):
                for seq_key in ("sequence", "designed_sequence"):
                    if isinstance(obj.get(seq_key), str):
                        return obj[seq_key]
    return None
