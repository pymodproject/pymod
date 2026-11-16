# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Retrieval of a structure from the ESM Atlas: the protocol, the fetcher and the plot.

This is the second thing PyMod does with ESM, beside the prediction which forge_protocols
runs: the Atlas is the database of what ESMFold2 has already folded, so it is read rather
than asked to compute. Both halves live in this package, each with a client of its own --
forge_client for the API which predicts, atlas_client for the one which is read.

The AlphaFold retrieval was the model for the shape of this one, and the shape is worth
keeping: a sequence is looked up in a database of predicted structures, the entry which was
found becomes an element, it is aligned with the query, and its structure is imported and
optionally cut to the part which matched. Nothing here is inherited from it, though. The
code is this module's own, so that either retrieval can change without moving the other,
and so that what is true of the Atlas -- that it is addressed by the sequence itself, that a
miss is an ordinary answer, that it reports a confidence and features of its own -- can be
said plainly instead of worked around.

What it does share is PyMod's own machinery, which belongs to no predictor: the protocol
base class, ClustalW, the structure association, and the BLAST search which supplies the
sequence similarity the Atlas cannot.
"""

import os
import re

import numpy as np
import requests
from Bio import SeqIO
from Bio.PDB import PDBParser
import Bio.PDB

from pymol import cmd

from pymod_lib.pymod_exceptions import PyModInterruptedProtocol
from pymod_lib.pymod_os_specific import check_network_connection
from pymod_lib.pymod_plot.pymod_plot_qt import PyMod_plot_window_qt
from pymod_lib.pymod_protocols.alignment_protocols.clustalw import Clustalw_regular_alignment
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.modeling_protocols.ESMFold import atlas_client
from pymod_lib.pymod_protocols.modeling_protocols.ESMFold._atlas_gui_qt import (
    ESM_atlas_search_window_qt)
from pymod_lib.pymod_protocols.modeling_protocols.modeling_protocols_common import (
    Associate_structure)
from pymod_lib.pymod_protocols.similarity_searches_protocols.local_blast import (
    LOC_BLAST_search)
from pymod_lib.pymod_seq.seq_conservation import get_conservation_symbol
from pymod_lib.pymod_threading import Protocol_exec_dialog


###################################################################################################
# The structures the Atlas answers with.                                                          #
###################################################################################################

class ESM_atlas_models_file(PyMod_protocol):
    """
    Puts the coordinates the Atlas returns into PyMod, as an element with a structure.

    There is no download here, which is why this is not the fetcher of a database of files:
    the Atlas answers with the coordinates inside its record, so they are written to a file
    and that file is parsed and associated with the element.
    """

    protocol_name = "fetch_esm_atlas_models"

    # The database this reads, as the user hears it named.
    models_database_name = "ESM Atlas"

    def __init__(self, pymod, output_directory=None, structures_to_fetch=None,
                 import_mode=None):
        PyMod_protocol.__init__(self, pymod, output_directory)
        self.structures_to_fetch = structures_to_fetch or []
        self.import_mode = import_mode

    def initialize_from_gui(self, mode, structures_to_fetch):
        self.structures_to_fetch = []
        if mode == "single":
            self.structures_to_fetch.append(structures_to_fetch)
        elif mode == "selection":
            self.structures_to_fetch.extend(self.pymod.get_selected_sequences())

    ###########################################################################
    # Reading the Atlas.                                                      #
    ###########################################################################

    def header_of(self, element):
        """The name of an element, however it can be read."""

        try:
            return self.pymod.get_header(element)
        except Exception:
            return getattr(element, "my_header", "") or ""

    def atlas_record_for(self, element):
        """
        The record of an element in the Atlas, or None.

        Two sequences are worth trying. The one the element carries is the obvious one, and
        for a hit of a BLAST search it is the aligned fragment of the subject, which is
        almost never a sequence anybody has deposited: the Atlas is keyed on whole
        sequences. So when that misses, the entry the hit belongs to is read from UniProt
        and asked about instead, which is the sequence the Atlas actually holds.
        """

        header = self.header_of(element)

        # The accession of the entry this element is of, kept for the name of the model:
        # it is the protein which was searched for, whatever the Atlas calls the sequence.
        accession = atlas_client.accession_in(header)
        self.followed_accession = accession

        sequence = atlas_client.clean_sequence(getattr(element, "my_sequence", "") or "")
        if sequence:
            try:
                return atlas_client.fetch_record(sequence)
            except atlas_client.Atlas_miss:
                pass
            except atlas_client.Atlas_unavailable as e:
                print("[ESM Atlas] %s." % e)
                return None

        if not accession:
            print("[ESM Atlas] %s: not in the Atlas, and its name carries no accession to"
                  " follow." % (header or "this hit"))
            return None

        whole = atlas_client.canonical_sequence(accession)
        if not whole:
            print("[ESM Atlas] %s: the sequence of %s could not be read from UniProt."
                  % (header, accession))
            return None
        try:
            record = atlas_client.fetch_record(whole)
        except (atlas_client.Atlas_miss, atlas_client.Atlas_unavailable) as e:
            print("[ESM Atlas] %s: %s." % (accession, e))
            return None
        print("[ESM Atlas] %s: the aligned fragment is not in the Atlas, the entry %s is."
              % (header, accession))
        return record

    ###########################################################################
    # Importing them.                                                         #
    ###########################################################################

    def fetch_models_files(self):
        """
        Reads the Atlas for every element to import, then loads what came back.
        """

        self.fetch_results_list = []

        try:
            if not self.pymod.use_protocol_threads:
                self.read_all_records()
            else:
                label_text = ("Reading the %s for %s. Please wait for completion..."
                              % (self.models_database_name,
                                 self.get_seq_text(self.structures_to_fetch, "structure")))
                p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                                function=self.read_all_records,
                                                args=(), title="Reading the ESM Atlas",
                                                wait_start=0.15, wait_end=0.15,
                                                label_text=label_text)
                p_dialog.exec_()

        except PyModInterruptedProtocol:
            return None

        # A sequence which is not in the Atlas is an ordinary answer, not a failure of the
        # connection, so it is said that way.
        n_misses = len([r for r in self.fetch_results_list if r["failure"]])
        if n_misses != 0:
            self.pymod.main_window.show_warning_message(
                "Not in the ESM Atlas",
                "The ESM Atlas holds no structure for %s of the %s sequences asked about."
                "\n\nThe Atlas is read by the exact sequence, so a protein which has never"
                " been deposited in that exact form is a miss: a designed sequence, a"
                " mutant, or a fragment cut at boundaries nobody has used. The console says"
                " which ones, and why."
                % (n_misses, len(self.structures_to_fetch)))
            if n_misses == len(self.fetch_results_list):
                return None

        for fetch_results in self.fetch_results_list:
            if not fetch_results["failure"]:
                try:
                    self._load_structure(fetch_results)
                except Exception as e:
                    self.pymod.main_window.show_error_message(
                        "Associate Structure Failure",
                        "The structure association for %s chain %s failed because of the"
                        " following error: %s" % (fetch_results["model_id"],
                                                  fetch_results.get("model_chain"), str(e)))

        self.pymod.main_window.gridder(clear_selection=True, update_elements=True,
                                       update_clusters=True)

    def read_all_records(self):
        for element in self.structures_to_fetch:
            self.fetch_results_list.append(self._fetch_single_element(element))

    def _fetch_single_element(self, old_element):
        results = {"element": old_element, "failure": True,
                   "model_id": atlas_client.atlas_hash(
                       getattr(old_element, "my_sequence", "") or "")[:12]}

        record = self.atlas_record_for(old_element)
        if record is None:
            return results

        # The element is named after the protein the user was looking for: the accession of
        # the hit, which is what the results window listed it as. The accession the Atlas
        # returns is its own -- it is keyed on the sequence, so for a sequence which also
        # sits in a metagenome it answers with an IMG/M gene number, which names nothing
        # which was searched for. It is kept among the figures instead.
        accession = str(getattr(self, "followed_accession", None)
                        or record.get("accession") or results["model_id"])
        name = "ESM-%s" % "".join(c if c.isalnum() or c in "-_." else "_" for c in accession)
        path = os.path.join(self.pymod.temp_directory_dirpath, "%s.pdb" % name)
        with open(path, "w") as structure_file:
            structure_file.write(atlas_client.pdb_text_of(record))

        print("[ESM Atlas] %s" % atlas_client.summary_of(record))

        # Kept so that the figures of this model can be shown again from the Models menu.
        # A list rather than a dictionary: the order is the order they were imported in,
        # which is the order the user saw them in the results window.
        if not hasattr(self.pymod, "esm_atlas_model_list"):
            self.pymod.esm_atlas_model_list = []
        self.pymod.esm_atlas_model_list.append({"model_name": name, "record": record,
                                                "element": old_element})
        results["model_id"] = name
        results["model_file_shortcut"] = path
        results["failure"] = False
        results["atlas_record"] = record
        return results

    def get_model_chain(self, model_file_shortcut, old_element):

        with open(model_file_shortcut, "r") as parsed_file_handle:
            parsed_biopython_structure = PDBParser(
                PERMISSIVE=1, QUIET=True).get_structure("temp", parsed_file_handle)

        for model in parsed_biopython_structure.get_list():
            for chain in model.get_list():
                old_element.model_chain = str(chain.id)

    def _load_structure(self, fetch_results):
        """
        Parses the file the coordinates were written to and associates it with the element.
        """

        old_element = fetch_results["element"]
        model_file_shortcut = fetch_results["model_file_shortcut"]

        self.get_model_chain(model_file_shortcut, old_element)
        fetch_results["model_chain"] = str(old_element.model_chain)
        model_chain = fetch_results["model_chain"]

        # Only the chain of the sequence which was looked up, or every chain of the file.
        # A model of the Atlas has one chain, so the two agree, but the choice is the
        # user's and it is honoured.
        #
        # This is the Associate_structure of modeling_protocols_common, not the one of
        # structural_databases_protocols: there are two classes of that name, and only this
        # one associates a model. It takes ASSOCIATE_ALL, it returns the new element, and it
        # files the coordinates under a session directory of their own in Models. The other
        # imports without complaint and fails here.
        if self.import_mode == "single-chain":
            a = Associate_structure(self.pymod, old_element)
        else:
            a = Associate_structure(self.pymod, old_element, ASSOCIATE_ALL=True)
        a.associate(model_file_shortcut, model_chain)

        self.new_element_with_structure = a.return_new_element()

        # Which element this model became. The record was written when the coordinates
        # arrived, and at that point the element which would carry them did not exist.
        for atlas_model in reversed(getattr(self.pymod, "esm_atlas_model_list", [])):
            if atlas_model.get("model_name") == fetch_results.get("model_id"):
                atlas_model["model_element"] = self.new_element_with_structure
                break


###################################################################################################
# The protocol.                                                                                   #
###################################################################################################

class ESM_atlas_search(PyMod_protocol):
    """
    Looks a sequence up in the ESM Atlas and imports the structure it finds.

    Two ways in, as the window offers them. By identifier, the accession in the name of the
    selected sequence is read from UniProt, the entry becomes an element, it is aligned with
    the query, and the Atlas is asked about it. By sequence, BLAST finds what resembles the
    selection and the Atlas is asked about each hit the user ticks -- because the Atlas
    cannot be asked for something similar: it is addressed by the exact sequence, and the
    search it offers of its own matches the embeddings of ESMC, which returns proteins of
    the same fold rather than of the same protein. Asked about the immunoglobulin domain of
    PD-L1 it answers with butyrophilins and a myelin glycoprotein, and never with PD-L1.
    """

    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)

        self.cut_structure = False
        self.atlas_search_window_class_qt = ESM_atlas_search_window_qt
        self.selected_sequences = self.pymod.get_selected_sequences()

    def launch_from_gui(self):
        """
        Called when the ESM Atlas lookup is clicked in the menu.
        """

        self.initialize_protocol_from_gui()

    ###########################################################################
    # Checks the input selection.                                             #
    ###########################################################################

    def initialize_protocol_from_gui(self):

        if not len(self.selected_sequences) > 0:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Please select at least one target sequence to look up in the ESM Atlas.")
            return None

        if len(self.selected_sequences) > 1:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Please select only one target sequence to look up in the ESM Atlas.")
            return None

        # The Atlas is a database of proteins, and ESMFold2 folds nothing else.
        if any(e.polymer_type in ["dna", "rna"] for e in self.selected_sequences):
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Cannot look up nucleic acids (RNA/DNA) in the ESM Atlas.")
            return None

        if False in [s.can_be_modeled() for s in self.selected_sequences]:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Please select only sequences that do not have a structure loaded in PyMOL.")
            return None

        self.build_atlas_search_window()

    def build_atlas_search_window(self):
        """
        Builds the ESM Atlas lookup window.
        """

        self.atlas_search_window = self.atlas_search_window_class_qt(
            parent=self.pymod.main_window, protocol=self)
        self.atlas_search_window.show()

    ###########################################################################
    # Running it.                                                             #
    ###########################################################################

    def launch_atlas_search(self):

        self.get_input_from_atlas_search_window()

        if self.search_choice == "by_id":
            self.check_existing_uniprot_id()

            if self.uniprot_header is not None:

                if self.new_structure_import_mode == "hit_fragment":
                    self.cut_structure = True

                self.atlas_search_window.close()
                self.by_id_mode_import_results_in_pymod()

            else:
                self.atlas_search_window.close()
                self.pymod.main_window.show_warning_message(
                    "ESM Atlas search Warning",
                    "The provided Sequence-ID does not match with a valid Uniprot-ID.\n"
                    "Please try with a 'Search by Sequence' instead.")

        else:
            self.atlas_search_window.close()
            self.query_element = self.selected_sequences[0]
            self.run_sequence_similarity_search()

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                       update_menus=True, update_elements=True)

    def get_input_from_atlas_search_window(self):

        # If the selected sequence is part of a cluster, PyMod asks whether to build a new
        # alignment or to extend the existing one.
        if self.selected_sequences[0].is_child():
            if self.atlas_search_window.al_mode_cb.currentText() == "Build new alignment":
                self.new_sequences_import_mode = "build-new"
            else:
                self.new_sequences_import_mode = "expand"
        else:
            self.new_sequences_import_mode = "build-new"

        if self.atlas_search_window.sequence_id_rb.isChecked():
            self.search_choice = "by_id"
        else:
            self.search_choice = "by_sequence"

        if self.atlas_search_window.import_mode_cb.currentText() == "Whole Structure":
            self.new_structure_import_mode = "all_structure"
        else:
            self.new_structure_import_mode = "hit_fragment"

    def check_existing_uniprot_id(self):
        """
        The accession in the name of the selection, and the entry it stands for.

        The Atlas has no index of accessions -- it takes the MD5 of a sequence and nothing
        else -- so an identifier is followed to UniProt first, and the sequence of that
        entry is what the Atlas is asked about.
        """

        for seq in self.selected_sequences:
            header = self.pymod.get_header(seq)

        uniprot_id_search = re.search(
            "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",
            header, re.IGNORECASE)

        if not check_network_connection("https://google.com", timeout=3):
            self.pymod.main_window.show_error_message(
                "Connection Error", "An internet connection is not available.")
            return None

        if not uniprot_id_search:
            self.uniprot_header = None
        else:
            self.uniprot_header = uniprot_id_search.group(0)
            path = os.path.join(self.pymod.temp_directory_dirpath,
                                self.uniprot_header + ".fasta")
            fastaurl = "https://rest.uniprot.org/uniprotkb/" + self.uniprot_header + ".fasta"
            response = requests.get(fastaurl)
            with open(path, "wb") as handle:
                handle.write(response.content)
            with open(path) as fasta_handle:
                fasta_records = list(SeqIO.parse(fasta_handle, "fasta"))
            for record in fasta_records:
                self.new_model_element = self.pymod.build_pymod_element_from_seqrecord(record)

    def fetch_structure_for(self, element, import_mode):
        """
        Reads the Atlas for an element and associates the structure with it.

        Returns the fetcher, whose 'new_element_with_structure' is the element which now
        carries the structure. When the Atlas had nothing, no element carries one, and the
        sequence of the entry stays in the main window, aligned with the query, which is
        worth keeping even without coordinates.
        """

        fetcher = ESM_atlas_models_file(self.pymod)
        fetcher.initialize_from_gui("single", element)
        fetcher.import_mode = import_mode
        fetcher.fetch_models_files()

        if not hasattr(fetcher, "new_element_with_structure"):
            fetcher.new_element_with_structure = element
            self.pymod.main_window.show_warning_message(
                "Not in the ESM Atlas",
                "The sequence was found and aligned, but the ESM Atlas holds no structure"
                " for it.\n\n"
                "To obtain a structure for it, predict one with 'Predict with ESMFold2"
                " (API key)' or with Boltz-2, or fetch the AlphaFold model of the entry"
                " with AlphaFold -> Fetch from database.")
        return fetcher

    ###########################################################################
    # The entry, aligned with the query, and its structure.                   #
    ###########################################################################

    def by_id_mode_import_results_in_pymod(self):
        """
        Builds a cluster with the query sequence as a mother and the entry as its child.
        """

        self.query_element = self.selected_sequences[0]

        if self.new_sequences_import_mode == "build-new":

            self.pymod.add_element_to_pymod(self.new_model_element)
            self.query_element.selected = True
            self.new_model_element.selected = True

            a = Clustalw_regular_alignment(self.pymod,
                                           output_directory=self.pymod.alignments_dirpath)
            a.rebuild_single_alignment_choice = False
            a.involved_clusters_list = []
            a.alignment_mode = "build-new-alignment"
            a.target_cluster_index = 0
            a.params_from_gui = {"selected_matrix": "blosum",
                                 "gapopen_value": 10.0,
                                 "gapextension_value": 0.2}
            a.elements_to_align_dict = {}
            # The cluster of a retrieved model, named after the database it was read from.
            a.AF = True
            a.models_label = "ESM Atlas"
            a.launch_no_gui()

        elif self.new_sequences_import_mode == "expand":

            selected_cluster = self.query_element.get_ancestor()

            self.pymod.add_element_to_pymod(self.new_model_element)
            self.query_element.selected = True
            self.new_model_element.selected = True

            a = Clustalw_regular_alignment(self.pymod,
                                           output_directory=self.pymod.alignments_dirpath)
            a.rebuild_single_alignment_choice = False
            a.involved_clusters_list = [selected_cluster]
            a.alignment_mode = "keep-previous-alignment"
            a.target_cluster_index = 0
            a.params_from_gui = {"selected_matrix": "blosum",
                                 "gapopen_value": 10.0,
                                 "gapextension_value": 0.2}
            a.elements_to_align_dict = {}
            a.AF = True
            a.models_label = "ESM Atlas"
            a.launch_no_gui()

        if self.cut_structure:
            self.cut_entry_to_the_matching_region()
            fp = self.fetch_structure_for(self.new_model_element, "single-chain")
        else:
            fp = self.fetch_structure_for(self.new_model_element, "multiple-chains")

        # It is a model, not an experimental structure, whichever database predicted it.
        fp.new_element_with_structure._is_af_model = True

        self.color_by_agreement_with_the_query(fp.new_element_with_structure)

    def cut_entry_to_the_matching_region(self):
        """
        Trims the entry, and the query with it, to the region the two share.

        Worked out from the alignment rather than from numbering: the entry of a database is
        a whole protein and the selection is often a domain of it, so where one begins
        inside the other is what the alignment says.
        """

        pc = 0                  # Alignment position counter.
        hc = 0                  # Query residue counter.
        tc = 0                  # Entry residue counter.
        matching_positions = []
        for hr, tr in zip(self.query_element.my_sequence,
                          self.new_model_element.my_sequence):
            if hr != "-" and tr != "-" and hr == tr:
                matching_positions.append({"pc": pc, "hc": hc, "tc": tc})
            if hr != "-":
                hc += 1
            if tr != "-":
                tc += 1
            pc += 1

        if not matching_positions:
            print("[ESM Atlas] the entry and the selection share no residue, so the whole"
                  " structure is imported.")
            return None

        starting_pos = matching_positions[0]["pc"]
        ending_pos = matching_positions[-1]["pc"]

        for element in (self.new_model_element, self.query_element):
            matching_seq = element.my_sequence[starting_pos:ending_pos + 1]
            gaps = "-" * (ending_pos - (starting_pos + len(matching_seq)))
            if self.new_sequences_import_mode == "expand":
                # Inside an existing alignment the columns have to stay where they are.
                element.my_sequence = "-" * starting_pos + matching_seq + gaps
            else:
                element.my_sequence = matching_seq + gaps

    def color_by_agreement_with_the_query(self, model_element):
        """
        Colours the imported model by how its residues agree with the query.

        Three colours: the same residue, a conservative difference, and a position the
        query does not have. It says at a glance which part of the model is the protein
        which was searched for and which part the entry brought with it.
        """

        scores_dict = {"-": 1, "*": 3, ":": 3, ".": 2, "missing": None}

        sequences = [self.query_element, model_element]
        alignment_length = len(sequences[0].my_sequence)
        list_of_alignment_columns = [[s.my_sequence[num] for s in sequences]
                                     for num in range(alignment_length)]

        residues = model_element.get_polymer_residues()
        rc = 0
        for seq_res, column in zip(model_element.my_sequence, list_of_alignment_columns):
            if column[0] == "-" and column[1] != "-":
                residues[rc].af_score = scores_dict[get_conservation_symbol(column)]
                rc += 1
            elif column[0] != "-" and column[1] != "-":
                if column[0] == column[1]:
                    residues[rc].af_score = scores_dict["*"]
                else:
                    residues[rc].af_score = scores_dict["."]
                rc += 1

        self.new_model_element._has_af_scores = True
        self.pymod.main_window.color_element_by_af_scores(model_element)

    ###########################################################################
    # Searching by sequence.                                                  #
    ###########################################################################

    def run_sequence_similarity_search(self):
        """
        BLAST finds what resembles the selection, and the Atlas supplies each structure.
        """

        self.pymod.esm_atlas_blast_search = ESM_atlas_blast_search(
            self.pymod, output_directory=self.pymod.similarity_searches_dirpath)
        self.pymod.esm_atlas_blast_search.launch_from_gui(AF=True)


class ESM_atlas_blast_search(LOC_BLAST_search):
    """
    PyMod's BLAST search, taking the structure of each hit from the ESM Atlas.

    This is how a search by sequence finds anything in the Atlas, which cannot be asked for
    something similar. BLAST is PyMod's own search tool and belongs to no predictor; what
    this says is only where the structures of the hits come from, and what to call them.
    Every hit of a UniProt database is in UniParc, and therefore in the Atlas.
    """

    models_fetcher_class = ESM_atlas_models_file
    models_database_name = "ESM Atlas"
    models_label = "ESM Atlas"

    def resolve_model_for_hit(self, hsp):
        """
        What to call a hit: ESM- and the accession of the entry it belongs to.

        No database is asked anything here. The accession comes from the name BLAST gave
        the subject, and it is the accession the Atlas record will be named after too.
        """

        try:
            title = self.get_subject_name(hsp)
        except Exception:
            title = ""
        accession = atlas_client.accession_in(title)
        if accession:
            return "ESM-%s" % accession
        sequence = atlas_client.clean_sequence(self.get_hsp_subj_seq(hsp))
        return "ESM-%s" % atlas_client.atlas_hash(sequence)[:12] if sequence else None


###################################################################################################
# The confidence of a model, along its sequence.                                                  #
###################################################################################################

def _highlight_in_pymol_from_plddt_plot(point_data):

    if point_data is not None:

        cmd.select("pymod_selection", point_data["pymol_selector"])
        cmd.center("pymod_selection")


def build_plddt_plot_window(pymod):
    """The window the confidence of an Atlas model is drawn in."""

    cp = PyMod_plot_window_qt(pymod.main_window)
    cp.initialize(pymod=pymod,
                  title="Predicted Local Distance Difference Test (pLDDT) per Position")
    cp.build_plotting_area(
        messagebar_initial_text=("Set the 'Interact on click' option to 'Yes' to click on"
                                 " the plot and to highlight residues in PyMOL."),
        update_messagebar=True,
        messagebar_text_on_update=("Selected: __residue_name__ __residue_pdb_position__ of"
                                   " __plot_name__, pLDDT value: __y__"),
        on_click_action=_highlight_in_pymol_from_plddt_plot,
        x_label_text="Positions",
        y_label_text="pLDDT")
    return cp


def plot_plddt_of_model(pymod, cp, model_name, model_file_path, model_element):
    """
    Adds to a plot the confidence of one model, and colours the model by it.

    ESMFold2 reports a confidence for every residue and the Atlas writes it in the B-factor
    column of the coordinates, which is where this reads it. The record carries the same
    numbers, but the file is what PyMOL has loaded, so the plot and the structure agree.
    """

    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(model_name, model_file_path)

    plddt_data = []
    res_numbers = []
    res_names = []
    pymol_selectors = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" not in residue:
                    continue
                pymol_selectors.append("%s and resi %s"
                                       % (structure.get_id(), residue.get_id()[1]))
                plddt_data.append(residue["CA"].get_bfactor())
                res_numbers.append(residue.get_id()[1])
                res_names.append(residue.get_resname())

    residue_additional_data = [{"residue_name": r_name,
                                "residue_pdb_position": r_position,
                                "pymol_selector": r_selector,
                                "export_label": "%s %s" % (r_name, r_position)}
                               for r_name, r_position, r_selector
                               in zip(res_names, res_numbers, pymol_selectors)]

    # The same values on the residues of the element, so that the main window can colour
    # the sequence by them, in ten bands between the lowest and the highest.
    plddt_scores = [score for res, score in zip(model_element.residues, plddt_data)
                    if res.is_polymer_residue()]
    if plddt_scores:
        bins = np.array(np.linspace(min(plddt_scores), max(plddt_scores), num=10))
        inds = np.digitize(np.array(plddt_scores), bins)
        for res, (score, bin_id) in zip(model_element.get_polymer_residues(),
                                        zip(plddt_scores, inds)):
            res.plddt_score = {"score": score, "interval": bin_id}
        model_element._has_plddt_scores = True
        pymod.main_window.color_element_by_plddt(model_element)

    cp.add_plot(res_numbers, plddt_data, label=model_name,
                additional_data=residue_additional_data)
    cp.show()
