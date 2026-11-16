# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Protein design with BoltzProt-1 in PyMod.

Designs binders against the sequences selected in PyMod, by calling the design
pipeline of the Boltz API. 'boltz_design' builds the request and 'boltz_client',
shared with the structure prediction, talks to the service.

Unlike a structure prediction, a design run is a pipeline which produces many
binders and is charged for each of them: the service refuses fewer than ten, so
the smallest run costs around $0.25 and a thousand designs around $25. The cost is
therefore always obtained from the service and shown before anything is started.
"""

import os

from pymol import cmd

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import PyModInterruptedProtocol

from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_design as bd
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_ligands as bligands
from pymod_lib.pymod_protocols.modeling_protocols.Boltz._residue_picking import (
    Residue_picking_mixin)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import (
    boltz_design_results as bdres)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_results as bres
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.design_gui_qt import (
    Boltz_design_window_qt)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.design_results_gui_qt import (
    Boltz_design_results_window_qt)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.boltz_client import (
    BoltzClient, BoltzAPIError, DEFAULT_BOLTZ_BASE_URL)


class BOLTZ_design(Residue_picking_mixin, PyMod_protocol):
    """
    Class to represent a BoltzProt-1 protein design run.
    """

    # A design run is a pipeline which generates and scores many binders, so it
    # takes far longer than a prediction. Reaching the ceiling does not stop it on
    # the service.
    poll_interval = 30
    max_wait = 6 * 3600

    # Set while a run is going on. Class attributes, so they are there however the
    # protocol was built.
    client = None
    _design_run = None
    _design_error = None
    _last_reported_status = None

    # The epitope as the user typed it, kept for the confirmation: the window is
    # closed by then, and the request counts residues differently.
    epitope_as_typed = ""

    # The surfaces to keep away from, as the user typed them.
    avoided_as_typed = ""

    # What detect_ligands() found, or None while it has not run.
    _detected_ligands = None

    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)
        self.design_window_class_qt = Boltz_design_window_qt
        self.selected_sequences = self.pymod.get_selected_sequences()


    def launch_from_gui(self):
        """
        Called when the design entry of the menu is clicked.
        """

        if not self.selected_sequences:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Please select the sequences of the target, that is the molecule"
                " the binder has to be designed against.")
            return None

        if False in [s.can_be_modeled() for s in self.selected_sequences]:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Alignments can not be used as a target. Please select single"
                " sequences, not an alignment.")
            return None

        self.build_design_window()


    def build_design_window(self):
        self.design_window = self.design_window_class_qt(
            parent=self.pymod.main_window,
            protocol=self,
            title="BoltzProt-1 Protein Design",
            upper_frame_title="Design a binder against the selected sequences",
            submit_command=self.run_design)
        self.design_window.show()


    ###############################################################################################
    # Credentials, shared in shape with the structure prediction protocol.                        #
    ###############################################################################################

    def get_boltz_credentials(self):
        try:
            api_key = self.pymod.boltz["api_key"].get_value()
            base_url = self.pymod.boltz["base_url"].get_value()
        except Exception:
            api_key, base_url = None, DEFAULT_BOLTZ_BASE_URL
        api_key = (api_key or "").strip() or None
        base_url = (base_url or "").strip() or DEFAULT_BOLTZ_BASE_URL
        return api_key, base_url


    def require_api_key(self):
        api_key, base_url = self.get_boltz_credentials()
        if not api_key:
            self.pymod.main_window.show_error_message(
                "Boltz API Key Missing",
                "No Boltz API key is configured. Insert one in the PyMod options"
                " window (Tools -> Options, 'Boltz service') before running a"
                " design.\n\nA key can be created in the Boltz API console.")
            return None
        return BoltzClient(api_key, base_url=base_url)


    ###############################################################################################
    # Building the request.                                                                       #
    ###############################################################################################

    def target_chain_assignment(self):
        """
        Which chain of the request each selected sequence is, as [(element, chain)].

        The chain a selected structure already has in PyMOL is kept, so a chain A
        stays chain A rather than being renamed. Sequences with no structure have
        none to keep and get the fallback names. The epitope is read against this
        same assignment, so the letters a user types are the letters the request
        carries.
        """

        preferred = []
        for element in self.selected_sequences:
            chain_id = None
            structure = getattr(element, "structure", None)
            if structure is not None:
                chain_id = getattr(structure, "chain_id", None)
            preferred.append(chain_id)
        chain_ids = bd.target_chain_ids(len(self.selected_sequences),
                                        preferred=preferred)
        return list(zip(self.selected_sequences, chain_ids))


    def detect_ligands(self, look_again=False, ask_dictionary=True, forget=False):
        """
        The ligands PyMod has, described so that the request can carry them: see
        'boltz_ligands'. The same list the prediction window offers, and for the
        same reason -- a ligand on screen should not have to be typed.
        """

        if look_again or self._detected_ligands is None:
            # 'forget' is the button in the window: ask the dictionary again
            # from scratch. A plain 'look_again' keeps what has been asked
            # already, which is what makes the second pass -- the one which
            # fills the names in after the window is up -- cost nothing.
            if forget:
                bligands.clear_ccd_cache()
            self._detected_ligands = bligands.detect_ligands(
                self.pymod, selected_elements=self.selected_sequences,
                ask_dictionary=ask_dictionary)
        return self._detected_ligands


    def ligand_chain_ids(self, count):
        """
        Chains for the ligands of the target, after the chains of the sequences and
        away from the binder's.
        """

        if not count:
            return []
        taken = [chain_id for _, chain_id in self.target_chain_assignment()]
        chain_ids = bd.target_chain_ids(len(taken) + count,
                                        preferred=taken + [None] * count)
        return chain_ids[len(taken):]


    def build_design_payload(self):
        """
        Reads the options window and returns the request of the design run,
        raising 'bd.BoltzDesignError' when it can not be built.
        """

        window = self.design_window

        # Kept because the window is closed before the run is confirmed, and the
        # confirmation shows the epitope in the numbering the user typed.
        self.epitope_as_typed = window.get_epitope_text()
        self.avoided_as_typed = window.get_avoided_text()

        assignment = self.target_chain_assignment()
        entities = []
        for element, chain_id in assignment:
            polymer_type = getattr(element, "polymer_type", "protein")
            if polymer_type == "dna":
                entities.append(bd.target_dna([chain_id], element.my_sequence))
            elif polymer_type == "rna":
                entities.append(bd.target_rna([chain_id], element.my_sequence))
            else:
                entities.append(bd.target_protein([chain_id], element.my_sequence))

        # The ligands ticked in the list join the target, each in a chain of its
        # own, in the notation the detection worked out for it.
        ligands = window.get_checked_ligands()
        for ligand, chain_id in zip(ligands, self.ligand_chain_ids(len(ligands))):
            if ligand.kind == bligands.KIND_SMILES:
                ligand_type = bd.TYPE_LIGAND_SMILES
            else:
                ligand_type = bd.TYPE_LIGAND_CCD
            entities.append(bd.target_ligand([chain_id], ligand.value, ligand_type))

        curated = window.get_curated_binder()
        return bd.build_design_payload(
            entities,
            binder_specification=window.get_binder_specification(),
            num_proteins=window.get_num_designs(),
            modality=window.get_modality() or bd.MODALITY_PEPTIDE,
            curated_binder=curated,
            epitope_residues=self.epitope_indices(self.epitope_as_typed),
            non_binding_residues=self.epitope_indices(self.avoided_as_typed))


    ###############################################################################################
    # The epitope: the residues of the target the binder should bind.                             #
    ###############################################################################################

    def epitope_indices(self, text):
        """
        The epitope of the request, from what the user typed: {chain_id: [0-based
        indices]}. The numbers are looked up among the residues of each target, so
        that the numbering of a structure becomes the indices the API counts.
        """

        return self.residues_to_indices(text, self.target_chain_assignment())


    def epitope_text_from_pymol(self, selection="sele"):
        """The epitope written out from a selection in PyMOL, for the button."""

        return self.residue_text_from_pymol(self.target_chain_assignment(),
                                            selection=selection)




    ###############################################################################################
    # Running the design.                                                                         #
    ###############################################################################################

    def run_design(self):
        """
        Called by the Submit button. Builds and validates the request, prices it
        with the service, and asks for a confirmation before it would be started.
        """

        client = self.require_api_key()
        if client is None:
            return None
        # Kept so the results window can download the structures of the designs.
        self.client = client

        try:
            payload = self.build_design_payload()
        except bd.BoltzDesignError as e:
            self.pymod.main_window.show_error_message("Boltz Design Error", str(e))
            return None

        confirm = self.design_window.get_confirm_cost()
        self.design_window.destroy()

        if not self.confirm_design(client, payload, confirm):
            return None

        print("- Running a BoltzProt-1 design:")
        print(bd.describe_design_payload(payload))

        run = self.submit_and_wait(client, payload)
        if run is None:
            return None

        try:
            self.show_designs(client, run)
        except (bdres.BoltzDesignResultError, BoltzAPIError) as e:
            self.pymod.main_window.show_error_message(
                "Boltz Design Error",
                "The design run finished, but its results could not be read:\n\n%s"
                % e)
            return None
        return None


    def submit_and_wait(self, client, payload):
        """
        Starts the design run and waits for it, off the GUI thread and behind a
        progress dialog which can be cancelled.

        Cancelling stops the waiting and NOT the run: it goes on generating designs
        on the service, and is charged for, so the message says where to find it.
        """

        self._design_run = None
        self._design_error = None

        def worker():
            try:
                started = client.start(bd.DESIGN_RESOURCE, payload)
                run_id = started.get("id")
                print("- The design run was started (id: %s). Waiting for it..."
                      % run_id)
                self._design_run = client.wait(
                    bd.DESIGN_RESOURCE, run_id,
                    poll_interval=self.poll_interval,
                    max_wait=self.max_wait,
                    on_status=self._report_status)
            except Exception as e:
                self._design_error = e

        dialog = Protocol_exec_dialog(
            app=self.pymod.main_window, pymod=self.pymod,
            function=worker, args=(),
            title="BoltzProt-1 design",
            label_text="Designing binders. A design run takes considerably longer"
                       " than a structure prediction.",
            wait_start=0.15, wait_end=0.15)
        try:
            dialog.exec_()
        except PyModInterruptedProtocol:
            print("- Waiting for the design run was cancelled. It keeps running on"
                  " the server, and is charged for: its designs can be fetched from"
                  " the Boltz console.")
            return None
        except Exception as e:
            self.pymod.main_window.show_error_message(
                "Boltz Design Error", "The design run could not be run:\n\n%s" % e)
            return None

        if self._design_error is not None:
            self.pymod.main_window.show_error_message(
                "Boltz Design Error",
                "The design run did not succeed:\n\n%s" % self._design_error)
            return None

        return self._design_run


    def _report_status(self, run):
        status = (run or {}).get("status")
        if status != self._last_reported_status:
            print("- The design run is %s." % status)
            self._last_reported_status = status


    ###############################################################################################
    # The designs.                                                                                #
    ###############################################################################################

    def show_designs(self, client, run):
        """
        Fetches the designs of a finished run and lists them for the user to choose
        from. They are not imported on their own: a run produces at least ten and
        may produce thousands.
        """

        run_id = (run or {}).get("id")
        response = client.list_results(bd.DESIGN_RESOURCE, run_id)
        designs = bdres.rank_results(response, run)

        run_name = "boltzprot_%s" % (self.get_design_run_count() + 1)
        print(bdres.format_report(designs, run_name))

        if not designs:
            self.pymod.main_window.show_warning_message(
                "No design",
                "The design run finished but produced no design.")
            return None

        self.results_window = Boltz_design_results_window_qt(
            self.pymod.main_window, protocol=self, designs=designs,
            run_name=run_name)
        self.results_window.show()
        return designs


    def get_design_run_count(self):
        """
        How many design runs this session has imported, so their names do not
        collide. Kept on the PyMod object, like the other counters.
        """
        if not hasattr(self.pymod, "performed_boltz_design_count"):
            self.pymod.performed_boltz_design_count = 0
        return self.pymod.performed_boltz_design_count


    def make_output_directory(self, run_name):
        """The directory of a design run, under the models of the project."""

        models_directory = os.path.join(self.pymod.current_project_dirpath,
                                        self.pymod.models_dirname)
        directory = os.path.join(models_directory, run_name)
        if not os.path.isdir(directory):
            os.makedirs(directory)
        return directory


    def import_designs(self, designs, run_name):
        """
        Imports the chosen designs into PyMod. Called by the results window, which
        stays open so more can be imported afterwards.

        Every design comes with the predicted structure of its complex with the
        target, so that is what is loaded: the binder and the target as chains of
        one model, which also brings the designed sequence in with it. A design
        whose structure can not be fetched falls back to its sequence alone, so
        something is still imported, and the difference is reported.
        """

        output_directory = self.make_output_directory(run_name)

        imported = []
        with_structure = 0
        sequence_only = []
        for design in designs:
            if not design.get("sequence"):
                continue
            name = "%s_%s" % (run_name, design.get("rank"))
            elements = self.import_one_design(design, name, output_directory)
            if elements is None:
                # No structure: keep the sequence, which is the point of the run.
                sequence_only.append(name)
                element = self.pymod.build_pymod_element_from_args(
                    name, design["sequence"])
                self.pymod.add_element_to_pymod(element)
                imported.append(element)
            else:
                with_structure += 1
                imported.extend(elements)

        if not imported:
            return []

        self.pymod.performed_boltz_design_count = self.get_design_run_count() + 1
        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                       update_menus=True, update_elements=True)

        message = ["%s designs were imported." % (with_structure + len(sequence_only))]
        if with_structure:
            message.append("%s of them as the predicted structure of the complex"
                           " with the target, so the binder and the target are"
                           " loaded as chains of one model." % with_structure)
        if sequence_only:
            message.append("%s could not have their structure downloaded and were"
                           " imported as sequences only: %s. The download links of"
                           " a run expire about half an hour after it finishes."
                           % (len(sequence_only), ", ".join(sequence_only)))
        text = "\n\n".join(message)
        print("- " + text.replace("\n\n", " "))
        self.pymod.main_window.show_info_message("Designs imported", text)
        return imported


    def import_one_design(self, design, name, output_directory):
        """
        Downloads the structure of one design and loads it as PyMod elements, or
        returns None when it carries no structure or it can not be fetched.
        """

        url = design.get("structure_url")
        if not url or self.client is None:
            return None

        cif_path = os.path.join(output_directory, "%s.cif" % name)
        pdb_path = os.path.join(output_directory, "%s.pdb" % name)
        try:
            self.client.download_to_file(url, cif_path)
            bres.cif_to_pdb(cif_path, pdb_path)
            return bres.load_model_in_pymod(self.pymod, pdb_path, name)
        except (BoltzAPIError, bres.BoltzResultError) as e:
            print("- The structure of %s could not be imported (%s). Its sequence"
                  " is imported instead." % (name, e))
            return None


    def confirm_design(self, client, payload, confirm):
        """
        Prices the run with the service and, when 'confirm' is set, shows it for
        approval. Returns False when the run should not go ahead.

        A design run is expensive enough that a request the service refuses stops
        everything: there is no point starting something it has already rejected.
        """

        cost_text = None
        try:
            estimate = client.estimate_cost(bd.DESIGN_RESOURCE, payload)
            cost = estimate.get("estimated_cost_usd")
            if cost is not None:
                cost_text = "$%s" % cost
        except BoltzAPIError as e:
            self.pymod.main_window.show_error_message(
                "Boltz Design Request Rejected",
                "The Boltz service did not accept the design request:\n\n%s" % e)
            return False

        if not confirm:
            return True

        message = ["The following design run will be started:", "",
                   bd.describe_design_payload(payload), ""]

        # The summary counts the residues of the sequence which is sent, as the
        # request does. The user typed the numbering of their structure, so their
        # own words are shown next to it rather than left to be recognised.
        epitope_text = self.epitope_as_typed
        if epitope_text:
            message.insert(3, "  epitope, as you gave it: %s" % epitope_text)
        if self.avoided_as_typed:
            message.insert(3, "  kept away from, as you gave it: %s"
                              % self.avoided_as_typed)
        if cost_text is not None:
            message.append("Estimated cost: %s" % cost_text)
            message.append("(an estimate: the final charge is computed at run time)")
        else:
            message.append("The cost of the run could not be estimated.")
        if client.uses_test_key():
            message.append("")
            message.append("The configured key belongs to the TEST environment:"
                           " the results will be synthetic and nothing will be"
                           " charged.")
        message.append("")
        message.append("Do you want to continue?")

        return askyesno_qt("Confirm the Boltz Design Run", "\n".join(message),
                           parent=self.pymod.get_qt_parent())
