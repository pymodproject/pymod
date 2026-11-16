# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Structure and binding prediction with Boltz-2 in PyMod.

The prediction is run by the Boltz REST API rather than by a local installation of
the model: 'boltz_client' talks to the service and 'boltz_payload' builds the
request. This module holds the protocol which collects the input from PyMod and
drives the two.
"""

import os

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import PyModInterruptedProtocol
import pymod_lib.pymod_structure as pmstr
import pymod_lib.pymod_vars as pmdt
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.__gui_qt import (Boltz_options_window_qt,
                                                                        MSA_LABEL_SERVER,
                                                                        MSA_LABEL_PYMOD,
                                                                        MSA_LABEL_NONE,
                                                                        LIGAND_LABEL_CCD,
                                                                        LIGAND_LABEL_SMILES)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_payload as bpay
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_ligands as bligands
from pymod_lib.pymod_protocols.modeling_protocols.Boltz._residue_picking import (
    Residue_picking_mixin)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_msa as bmsa
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_results as bres
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.results_gui_qt import (
    Boltz_results_window_qt)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz.boltz_client import (BoltzClient,
                                                                            BoltzAPIError,
                                                                            DEFAULT_BOLTZ_BASE_URL,
                                                                            STRUCTURE_RESOURCE)
from pymol.Qt import QtWidgets
from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt


###################################################################################################
# Boltz-2 structure and binding prediction.                                                       #
###################################################################################################


class Boltz_model_session:
    """
    One predicted model, as the "Models" menu needs it.

    A run produces several models and each has its own confidence figures, its own
    chains and its own predicted aligned error, so each gets one of these rather
    than the whole run sharing a single record.
    """

    def __init__(self, model_name, job, sample, sample_index, full_output_path,
                 elements, output_directory, pae_matrix=None):
        self.model_name = model_name
        # The whole response, so the figures of the run can still be read, and the
        # sample this model actually is.
        self.job = job
        self.sample = sample
        self.sample_index = sample_index
        self.full_output_path = full_output_path
        self.new_model_elements = list(elements)
        self.new_model_element = self.new_model_elements[0] if elements else None
        self.output_directory = output_directory
        self.pae_matrix = pae_matrix


class BOLTZ_model(Residue_picking_mixin, PyMod_protocol):

    """
    Class to represent a structure and binding prediction run with Boltz-2.
    """

    # How many homologues of an alignment are sent at most. The alignment travels
    # base64 encoded inside the JSON body, so a very deep one would make the
    # request large. A constant of the class rather than an instance attribute, so
    # that it is there however the protocol was built.
    max_msa_sequences = 4000

    # What detect_ligands() found, or None while it has not run. On the class, so
    # that a protocol built any other way than through __init__ -- as the tests
    # build it -- still has it.
    _detected_ligands = None

    # The pocket as the user typed it, kept for the confirmation: the window is
    # closed by then, and the request counts residues from zero.
    pocket_as_typed = ""

    # How the prediction is waited for. A small job takes under a minute, a large
    # complex with an automatic alignment can take much longer, so the ceiling is
    # generous; reaching it does not stop the job on the service.
    poll_interval = 10
    max_wait = 3600

    # Set while a prediction runs. Class attributes so that they exist however the
    # protocol was built.
    client = None
    _last_reported_status = None


    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)

        self.cut_structure = False

        self.boltz_initialization()

        self.selected_sequences = self.pymod.get_selected_sequences()


    def boltz_initialization(self):
        self.boltz_modeling_window_class_qt = Boltz_options_window_qt


    def launch_from_gui(self):
        """
        Called when the Boltz entry of the menu is clicked.
        """

        self.initialize_modeling_protocol_from_gui()


    ###########################################################################
    # Checks the input selection.                                             #
    ###########################################################################

    def initialize_modeling_protocol_from_gui(self):

        #---------------------------------------------------------------------
        # Get the selected sequences to see if there is a correct selection. -
        #---------------------------------------------------------------------

        # First check if at least one sequence is selected.
        if not len(self.selected_sequences) > 0:
            self.pymod.main_window.show_error_message("Selection Error", "Please select at least one target sequence to run a Boltz-2 prediction.")
            return None

        # Nucleic acids used to be refused here, because the local Boltz-1 path
        # was never finished for them. The API accepts DNA and RNA entities, so
        # protein-DNA and protein-RNA complexes can be predicted and the check
        # was removed.

        # Checks if all the selected sequences can be used to build a model.
        # 'can_be_modeled' only excludes cluster elements, that is alignments and
        # search results; the message used to claim it refused sequences with a
        # structure loaded, which it never did.
        if False in [s.can_be_modeled() for s in self.selected_sequences]:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Alignments can not be predicted. Please select single sequences,"
                " not an alignment.")
            return None

        self.build_boltz_modeling_window()


    def build_boltz_modeling_window(self):
        """
        Builds the window with the options of the prediction.
        """

        self.boltz_modeling_window = self.boltz_modeling_window_class_qt(
            parent=self.pymod.main_window,
            protocol=self,
            title="Boltz-2 Structure Prediction",
            upper_frame_title="Options for the Boltz-2 prediction",
            submit_command=self.run_boltz_modeling)
        self.boltz_modeling_window.show()


    ###############################################################################################
    # Importing a prediction which already ran.                                                   #
    ###############################################################################################

    def launch_import_by_id_from_gui(self):
        """
        Imports a prediction which the service has already finished, given its
        identifier.

        The results of a job are kept for a week while its download links last only
        about half an hour, so a run whose import failed, or which was cancelled
        while being waited for, is not lost: it can be fetched again here instead of
        being paid for a second time. No selection is needed, and nothing is charged.
        """

        client = self.require_api_key()
        if client is None:
            return None
        self.client = client

        job_id, confirmed = QtWidgets.QInputDialog.getText(
            self.pymod.get_qt_parent(),
            "Import a Boltz prediction",
            "Identifier of the prediction to import.\n\n"
            "The results of a prediction are kept for a week, so one whose import"
            " failed, or which was cancelled while being waited for, can be fetched"
            " here without running it again. Identifiers look like"
            " 'sab_pred_...' and are shown in the console when a prediction starts,"
            " and in the Boltz console.")
        job_id = (job_id or "").strip()
        if not confirmed or not job_id:
            return None

        try:
            job = client.retrieve(STRUCTURE_RESOURCE, job_id)
        except BoltzAPIError as e:
            self.pymod.main_window.show_error_message(
                "Boltz Prediction Not Found",
                "The prediction %r could not be fetched:\n\n%s" % (job_id, e))
            return None

        status = (job or {}).get("status")
        if status != "succeeded":
            if status == "failed":
                message = ("That prediction did not succeed: %s"
                           % client.job_error_message(job))
            elif status in ("pending", "running"):
                message = ("That prediction is still %s. Try again once it has"
                           " finished." % status)
            else:
                message = "That prediction is in the state %r." % status
            self.pymod.main_window.show_error_message("Boltz Prediction Not Ready",
                                                      message)
            return None

        if job.get("data_deleted_at"):
            self.pymod.main_window.show_error_message(
                "Boltz Results Deleted",
                "The results of that prediction were deleted on %s. The results of"
                " a prediction are only kept for a week, so it has to be run"
                " again." % job["data_deleted_at"])
            return None

        try:
            self.import_prediction(job)
        except (bres.BoltzResultError, BoltzAPIError) as e:
            self.pymod.main_window.show_error_message(
                "Boltz Import Error",
                "The prediction was fetched, but its results could not be"
                " imported:\n\n%s" % e)
            return None
        return None


    ###############################################################################################
    # Information the options window asks for.                                                    #
    ###############################################################################################

    def detect_ligands(self, look_again=False, ask_dictionary=True, forget=False):
        """
        The ligands PyMod has, described so that a request can carry them: see
        'boltz_ligands'. The result is kept, because a scan asks the chemical
        component dictionary over the network; 'look_again' throws it away, which
        is what the button in the window does.

        Every structure open in PyMod is looked at, not only the selected
        sequences: a ligand is usually loaded with the structure it came from,
        while the sequence being predicted has no structure at all. The ligands of
        the selected sequences come first in the list.
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


    def get_available_heteroresidues(self):
        """
        The CCD codes of the ligands found, kept for the tests and for anything
        which asked for codes alone. What the window uses is 'detect_ligands',
        which also says how each ligand should be sent.
        """

        return [ligand.code for ligand in self.detect_ligands()
                if ligand.kind in (bligands.KIND_CCD, bligands.KIND_CCD_UNVERIFIED)]


    def get_msa_option_labels(self):
        """
        The alignment options to offer. The alignment of PyMod is only proposed
        when at least one of the selected proteins actually belongs to one, so the
        choice is never offered when it could not be honoured.
        """

        labels = [MSA_LABEL_SERVER]
        if self.get_elements_with_an_alignment():
            labels.append(MSA_LABEL_PYMOD)
        labels.append(MSA_LABEL_NONE)
        return labels


    def get_elements_with_an_alignment(self):
        """
        The selected proteins which belong to an alignment in PyMod, and can
        therefore have their own alignment sent to the service.
        """

        elements = []
        for element in self.selected_sequences:
            if getattr(element, "polymer_type", "protein") != "protein":
                continue
            if bmsa.alignment_of_element(element):
                elements.append(element)
        return elements


    def get_boltz_credentials(self):
        """
        Returns '(api_key, base_url)' from the PyMod preferences. The key is None
        when it has not been configured; the base URL always falls back to the
        default, since 'get_value' returns None until a parameter is set.
        """

        try:
            api_key = self.pymod.boltz["api_key"].get_value()
            base_url = self.pymod.boltz["base_url"].get_value()
        except Exception:
            api_key, base_url = None, DEFAULT_BOLTZ_BASE_URL
        api_key = (api_key or "").strip() or None
        base_url = (base_url or "").strip() or DEFAULT_BOLTZ_BASE_URL
        return api_key, base_url


    def require_api_key(self):
        """
        Returns a client, or None after telling the user that the key is missing.
        """

        api_key, base_url = self.get_boltz_credentials()
        if not api_key:
            self.pymod.main_window.show_error_message(
                "Boltz API Key Missing",
                "No Boltz API key is configured. Insert one in the PyMod options"
                " window (Tools -> Options, 'Boltz service') before running a"
                " prediction.\n\nA key can be created in the Boltz API console.")
            return None
        return BoltzClient(api_key, base_url=base_url)


    ###############################################################################################
    # Building the request.                                                                       #
    ###############################################################################################

    def build_prediction_payload(self):
        """
        Reads the options window and returns the request body of the prediction,
        raising 'bpay.BoltzPayloadError' with a message meant for the user when
        the selection can not be turned into a valid request.
        """

        window = self.boltz_modeling_window

        # Chain identifiers are assigned here, in the order of the selection: they
        # can not be taken from the loaded structures, since a sequence being
        # predicted usually has none.
        n_ligands = len(window.get_checked_ligands())
        if window.get_typed_ligand():
            n_ligands += 1
        chain_ids = bpay.chain_id_sequence(len(self.selected_sequences) + n_ligands)
        next_chain = iter(chain_ids)

        msa_label = window.get_msa_label()

        entities = []
        polymer_assignment = []
        for element in self.selected_sequences:
            chain_id = next(next_chain)
            polymer_assignment.append((element, chain_id))
            polymer_type = getattr(element, "polymer_type", "protein")
            if polymer_type == "dna":
                entities.append(bpay.polymer_entity(bpay.TYPE_DNA, [chain_id],
                                                    element.my_sequence))
            elif polymer_type == "rna":
                entities.append(bpay.polymer_entity(bpay.TYPE_RNA, [chain_id],
                                                    element.my_sequence))
            else:
                msa = self.build_msa_for_element(element, msa_label, chain_id)
                entities.append(bpay.polymer_entity(bpay.TYPE_PROTEIN, [chain_id],
                                                    element.my_sequence, msa=msa))

        # Ligands: the ones ticked in the list are sent in the notation the
        # detection worked out for each -- its code where the dictionary has one,
        # the SMILES read from its coordinates where it does not -- and the typed
        # one may be either.
        ligand_chain_ids = []
        for ligand in window.get_checked_ligands():
            chain_id = next(next_chain)
            ligand_chain_ids.append(chain_id)
            if ligand.kind == bligands.KIND_SMILES:
                ligand_type = bpay.TYPE_LIGAND_SMILES
            else:
                ligand_type = bpay.TYPE_LIGAND_CCD
            entities.append(bpay.ligand_entity([chain_id], ligand.value,
                                               ligand_type))

        typed_ligand = window.get_typed_ligand()
        if typed_ligand:
            chain_id = next(next_chain)
            ligand_chain_ids.append(chain_id)
            kind_label = window.get_ligand_kind_label()
            if kind_label == LIGAND_LABEL_CCD:
                ligand_type = bpay.TYPE_LIGAND_CCD
            elif kind_label == LIGAND_LABEL_SMILES:
                ligand_type = bpay.TYPE_LIGAND_SMILES
            else:
                ligand_type = None  # guessed from the value
            entities.append(bpay.ligand_entity([chain_id], typed_ligand,
                                               ligand_type))

        binding = None
        if window.get_affinity_requested():
            if not ligand_chain_ids:
                raise bpay.BoltzPayloadError(
                    "The binding affinity can only be predicted for a ligand, but"
                    " no ligand was given. Add a ligand, or set 'Predict binding"
                    " affinity' to No.")
            # The affinity is predicted for one binder. With more than one ligand
            # the first is used, and the user is told which.
            binding = bpay.ligand_protein_binding(ligand_chain_ids[0])
            if len(ligand_chain_ids) > 1:
                print("- More than one ligand was given: the binding affinity"
                      " will be predicted for the one in chain %s."
                      % ligand_chain_ids[0])

        # A pocket: where the ligand should sit. Without one the prediction is
        # free to place it anywhere, which for a protein with several clefts is
        # often not where the chemistry says it goes.
        constraints = []
        self.pocket_as_typed = window.get_pocket_text()
        if self.pocket_as_typed:
            if not ligand_chain_ids:
                raise bpay.BoltzPayloadError(
                    "A pocket says where a ligand should sit, but no ligand was"
                    " given. Add one, or leave the pocket empty.")
            contact_residues = self.residues_to_indices(self.pocket_as_typed,
                                                        polymer_assignment)
            if contact_residues:
                constraints.append(bpay.pocket_constraint(
                    ligand_chain_ids[0], contact_residues, entities,
                    max_distance=window.get_pocket_distance(),
                    force=window.get_pocket_forced()))

        return bpay.build_payload(entities,
                                  num_samples=window.get_num_samples(),
                                  binding=binding,
                                  constraints=constraints)


    def polymer_chain_assignment(self):
        """
        Which chain of the request each selected sequence will be, for the window
        to name and for a pocket to be read against.

        The identifiers are assigned in the order of the selection, as
        build_prediction_payload assigns them; ligands come after the sequences, so
        the sequences keep their letters whatever is ticked.
        """

        chain_ids = bpay.chain_id_sequence(len(self.selected_sequences))
        return list(zip(self.selected_sequences, chain_ids))


    def pocket_text_from_pymol(self, selection="sele"):
        """The pocket written out from a selection in PyMOL, for the button."""

        return self.residue_text_from_pymol(self.polymer_chain_assignment(),
                                            selection=selection)


    def build_msa_for_element(self, element, msa_label, chain_id):
        """
        Returns the MSA descriptor of one protein chain, following what was chosen
        in the options window.

        With the alignment of PyMod, a protein which belongs to no alignment gets
        an empty MSA rather than an automatic one: the service refuses a request
        which mixes automatic and custom alignments, while custom and empty may be
        combined. This is reported, so the user knows that chain went in without
        one.
        """

        if msa_label == MSA_LABEL_NONE:
            return bpay.msa_empty()

        if msa_label != MSA_LABEL_PYMOD:
            return bpay.msa_automatic()

        try:
            a3m_text = bmsa.a3m_for_element(element, max_sequences=self.max_msa_sequences)
        except bmsa.BoltzMSAError as e:
            raise bpay.BoltzPayloadError(
                "The alignment of chain %s can not be used: %s" % (chain_id, e))

        if a3m_text is None:
            print("- Chain %s does not belong to an alignment in PyMod: it is sent"
                  " without one." % chain_id)
            return bpay.msa_empty()

        print("- Chain %s: sending the alignment of PyMod, %s sequences."
              % (chain_id, bmsa.count_sequences(a3m_text)))
        return bpay.msa_custom(a3m_text)


    def run_boltz_modeling(self):
        """
        Called by the Submit button. Builds and validates the request, shows what
        it will cost, and asks for a confirmation before it is sent.
        """

        client = self.require_api_key()
        if client is None:
            return None
        # Kept so the import can download the results with the same credentials.
        self.client = client

        try:
            payload = self.build_prediction_payload()
        except bpay.BoltzPayloadError as e:
            self.pymod.main_window.show_error_message("Boltz Input Error", str(e))
            return None

        confirm = self.boltz_modeling_window.get_confirm_cost()
        self.boltz_modeling_window.destroy()

        if not self.confirm_prediction(client, payload, confirm):
            return None

        print("- Running a Boltz-2 prediction:")
        print(bpay.describe_payload(payload))

        job = self.submit_and_wait(client, payload)
        if job is None:
            return None

        try:
            self.import_prediction(job)
        except (bres.BoltzResultError, BoltzAPIError) as e:
            self.pymod.main_window.show_error_message(
                "Boltz Import Error",
                "The prediction finished, but its results could not be imported:"
                "\n\n%s\n\nThe prediction itself was not lost: it succeeded and"
                " its results are kept on the service for a week (identifier %s)."
                % (e, (job or {}).get("id")))
            return None
        return None


    ###############################################################################################
    # Running the prediction.                                                                     #
    ###############################################################################################

    def submit_and_wait(self, client, payload):
        """
        Submits the prediction and waits for it, off the GUI thread and behind a
        progress dialog which can be cancelled, so PyMod stays usable while the
        service works.

        Returns the finished job, or None when it failed or was cancelled. Note
        that cancelling stops the waiting, not the job: it keeps running on the
        service and can be fetched from the Boltz console.
        """

        # The worker can not return a value, so it leaves the result here.
        self._boltz_job = None
        self._boltz_error = None

        def worker():
            try:
                started = client.start(STRUCTURE_RESOURCE, payload)
                job_id = started.get("id")
                print("- The prediction was submitted (id: %s). Waiting for it..."
                      % job_id)
                self._boltz_job = client.wait(
                    STRUCTURE_RESOURCE, job_id,
                    poll_interval=self.poll_interval,
                    max_wait=self.max_wait,
                    on_status=self._report_status)
            except Exception as e:
                # Reported on the main thread: a dialog can not be raised here.
                self._boltz_error = e

        dialog = Protocol_exec_dialog(
            app=self.pymod.main_window, pymod=self.pymod,
            function=worker, args=(),
            title="Boltz-2 prediction",
            label_text="Waiting for the Boltz service. This may take a few"
                       " minutes.",
            wait_start=0.15, wait_end=0.15)
        try:
            dialog.exec_()
        except PyModInterruptedProtocol:
            print("- Waiting for the Boltz prediction was cancelled. The job may"
                  " still be running on the server.")
            return None
        except Exception as e:
            self.pymod.main_window.show_error_message(
                "Boltz Error", "The prediction could not be run:\n\n%s" % e)
            return None

        if self._boltz_error is not None:
            self.pymod.main_window.show_error_message(
                "Boltz Error", "The prediction did not succeed:\n\n%s"
                % self._boltz_error)
            return None

        return self._boltz_job


    def _report_status(self, job):
        status = (job or {}).get("status")
        if status != self._last_reported_status:
            print("- The prediction is %s." % status)
            self._last_reported_status = status


    ###############################################################################################
    # Importing the results.                                                                      #
    ###############################################################################################

    def import_prediction(self, job):
        """
        Downloads every predicted model, loads them in PyMOL and PyMod, and reports
        the figures.

        Each model becomes a session of its own in the "Models" menu, not one
        session for the whole run. A run asks for several samples by default, and
        each has its own confidence figures, its own chains and its own predicted
        aligned error; a single session for all of them meant one profile drawn over
        every model at once and, worse, the matrix of one sample being labelled with
        the residues of whichever model happened to be loaded first.
        """

        output_directory = self.make_output_directory()
        samples = bres.sample_results(job)
        run_number = self.pymod.performed_boltz_modeling_count + 1
        # "boltz_model_1" rather than "boltz_1": the latter reads as Boltz-1, the
        # model, which is not what it means.
        run_root = "boltz_model_%s" % run_number

        # The archive is the only place the predicted aligned error is, and it holds
        # one matrix per sample, so it is fetched before the models are registered.
        archive_path = self.download_archive(job, output_directory, run_root)

        loaded_elements = []
        sessions = []
        failed = []
        for index, sample in enumerate(samples):
            model_name = (run_root if len(samples) == 1
                          else "%s_%s" % (run_root, index + 1))
            cif_path = os.path.join(output_directory, "%s.cif" % model_name)
            pdb_path = os.path.join(output_directory, "%s.pdb" % model_name)

            print("- Downloading %s..." % model_name)
            try:
                self.client.download_to_file(bres.structure_url_of(sample), cif_path)
                bres.cif_to_pdb(cif_path, pdb_path)
                elements = self.load_model_in_pymod(pdb_path, model_name)
            except (BoltzAPIError, bres.BoltzResultError) as e:
                # One model failing must not lose the others, nor the prediction:
                # it has been paid for and stays on the service for a week.
                print("- %s could not be imported: %s" % (model_name, e))
                failed.append(model_name)
                continue

            loaded_elements.extend(elements)
            sessions.append(self.register_model(
                job=job, sample=sample, sample_index=index, model_name=model_name,
                pdb_path=pdb_path, elements=elements,
                output_directory=output_directory, archive_path=archive_path))

        if not loaded_elements:
            raise bres.BoltzResultError(
                "None of the %s models of the prediction could be downloaded.\n\n"
                "The prediction itself succeeded and is kept on the service for a"
                " week, so it does not have to be run again: its identifier is %s."
                % (len(samples), job.get("id")))

        if failed:
            print("- %s of the %s models could not be imported (%s), the others were."
                  % (len(failed), len(samples), ", ".join(failed)))

        self.pymod.performed_boltz_modeling_count += 1

        print(bres.format_report(job, run_root))

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                       update_menus=True, update_elements=True)

        # The window shows the best model, which is the one the figures of the run
        # refer to.
        best = sessions[0] if sessions else None
        self.results_window = Boltz_results_window_qt(
            self.pymod.main_window, job=job, model_name=best.model_name if best
            else run_root, output_directory=output_directory,
            sample=best.sample if best else None)
        self.results_window.show()

        return loaded_elements


    def download_archive(self, job, output_directory, run_root):
        """
        Downloads the results archive, which holds the predicted aligned error of
        every sample. Returns its path, or None: a missing archive costs the plots,
        not the models.
        """

        archive_url = bres.archive_url_of(job)
        if not archive_url:
            return None
        archive_path = os.path.join(output_directory,
                                    "%s_all_results.tar.gz" % run_root)
        try:
            self.client.download_to_file(archive_url, archive_path)
            print("- The full results of the service were saved in %s."
                  % archive_path)
            return archive_path
        except BoltzAPIError as e:
            print("- The results archive could not be downloaded, so the predicted"
                  " aligned error will not be available: %s." % e)
            return None


    def register_model(self, job, sample, sample_index, model_name, pdb_path,
                       elements, output_directory, archive_path):
        """
        Registers one model for the "Models" menu, carrying only what belongs to it:
        its own chains, its own figures and its own predicted aligned error.
        """

        session = Boltz_model_session(
            model_name=model_name,
            job=job,
            sample=sample,
            sample_index=sample_index,
            full_output_path=pdb_path,
            elements=elements,
            output_directory=output_directory,
            pae_matrix=bres.pae_matrix_from_archive(archive_path,
                                                    sample=sample_index))
        if not hasattr(self.pymod, "boltz_modeling_session_list"):
            self.pymod.boltz_modeling_session_list = []
        self.pymod.boltz_modeling_session_list.append(session)
        return session


    def make_output_directory(self):
        """The directory of this prediction, under the models of the project."""

        models_directory = os.path.join(self.pymod.current_project_dirpath,
                                        self.pymod.models_dirname)
        name = "%s_%s_boltz" % (self.pymod.models_subdirectory,
                                self.pymod.performed_boltz_modeling_count + 1)
        directory = os.path.join(models_directory, name)
        if not os.path.isdir(directory):
            os.makedirs(directory)
        return directory


    def load_model_in_pymod(self, pdb_path, model_name):
        """
        Loads a predicted model as PyMod elements, one per polymer chain, and
        colours them by the per-residue confidence the service put in the B-factor
        column of the file. A ligand is not an element of its own: it stays in the
        PyMOL object as heteroatoms, as it does for any structure PyMod loads.
        """

        parsed = pmstr.Parsed_model_pdb_file(
            self.pymod, pdb_path,
            output_directory=self.pymod.structures_dirpath,
            new_file_name=model_name, model_root_name=model_name)

        colors = pmdt.pymol_light_colors_list
        elements = []
        for chain_number, element in enumerate(parsed.get_pymod_elements()):
            self.pymod.add_element_to_pymod(
                element, color=colors[chain_number % len(colors)])
            elements.append(element)

        return elements


    def confirm_prediction(self, client, payload, confirm):
        """
        Asks the service what the request would cost and, when 'confirm' is set,
        shows it to the user for approval. Returns False when the prediction
        should not be run.

        The cost estimate is free and does not run the prediction. A failure to
        obtain it is reported but does not by itself stop the run, since the
        estimate is only informative.
        """

        cost_text = None
        try:
            estimate = client.estimate_cost(STRUCTURE_RESOURCE, payload)
            cost = estimate.get("estimated_cost_usd")
            if cost is not None:
                cost_text = "$%s" % cost
        except BoltzAPIError as e:
            # A rejected request is worth stopping for: the service has refused
            # the input, so running it would fail as well.
            self.pymod.main_window.show_error_message(
                "Boltz Request Rejected",
                "The Boltz service did not accept the request:\n\n%s" % e)
            return False

        if not confirm:
            return True

        message = ["The following prediction will be run:", "",
                   bpay.describe_payload(payload), ""]
        if cost_text is not None:
            message.append("Estimated cost: %s" % cost_text)
            message.append("(an estimate: the final charge is computed from the"
                           " exact size of the input)")
        else:
            message.append("The cost of the prediction could not be estimated.")
        if client.uses_test_key():
            message.append("")
            message.append("The configured key belongs to the TEST environment:"
                           " the result will be synthetic and nothing will be"
                           " charged.")
        message.append("")
        message.append("Do you want to continue?")

        return askyesno_qt("Confirm the Boltz Prediction", "\n".join(message),
                           parent=self.pymod.get_qt_parent())
