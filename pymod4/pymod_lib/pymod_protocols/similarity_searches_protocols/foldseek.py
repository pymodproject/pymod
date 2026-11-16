# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Searching for structures which look like a structure of the project.

Every other search in PyMod starts from a sequence; this one starts from
coordinates, and finds proteins whose fold resembles the query even when nothing
in their sequence says so. The server is the one behind search.foldseek.com.

What comes back is made into ordinary PyMod elements. A hit of the PDB is named
'1crn_A' and a hit of an AlphaFold database keeps its 'AF-P01542-F1-model_v6'
name, and those are names PyMod itself reads: the structure of a chosen hit is
then fetched and associated by the same protocols which serve the BLAST results
and the 'Fetch from the PDB' command, so a hit behaves like every other
structural element of the project.

**One chain per search.** The server searches the first chain of what it is sent
and ignores the rest -- a two chain file gave exactly the hits of its chain A --
so PyMod sends one chain at a time and says which one it is sending.
"""

import os

from pymol import cmd

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import PyModInterruptedProtocol
from pymod_lib.pymod_seq.seq_star_alignment import save_cstar_alignment

from pymod_lib.pymod_protocols.similarity_searches_protocols import foldseek_client as fsc
from pymod_lib.pymod_protocols.similarity_searches_protocols import foldseek_results as fsr
from pymod_lib.pymod_protocols.similarity_searches_protocols.foldseek_gui_qt import (
    Foldseek_options_window_qt, Foldseek_results_window_qt)

from pymod_lib.pymod_protocols.structural_databases_protocols import Fetch_structure_file
# The association of a model, and the download of one, are taken from the module
# which owns them. Note that two classes are called Associate_structure in PyMod:
# this is the one which goes with the models.
from pymod_lib.pymod_protocols.modeling_protocols.modeling_protocols_common import (
    Associate_structure as Associate_model_structure, fetch_models_file)


class Foldseek_search(PyMod_protocol):
    """
    A structure search: one chain of the project against the structure databases.
    """

    protocol_name = "foldseek"
    search_label = "Foldseek"

    # What the server offered, filled in before the window is built.
    databases = ()

    # Set while a search is going on.
    ticket = None
    payload = None
    results = None
    query_element = None
    query_chain = None
    databases_searched = ()

    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)
        self.options_window_class_qt = Foldseek_options_window_qt
        self.results_window_class_qt = Foldseek_results_window_qt
        self.selected_sequences = self.pymod.get_selected_sequences()


    ###############################################################################################
    # Starting.                                                                                   #
    ###############################################################################################

    def launch_from_gui(self):
        """
        Called from the menu. One structure is searched at a time, because the
        server searches one chain.
        """

        if not fsc.has_requests:
            self.pymod.main_window.show_error_message(
                "Component Missing",
                "The 'requests' module is needed to reach the Foldseek server,"
                " and it is not installed in this PyMOL.")
            return None

        if len(self.selected_sequences) != 1:
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "A structure search starts from one structure: please select"
                " exactly one, rather than %s." % len(self.selected_sequences))
            return None

        self.query_element = self.selected_sequences[0]

        if not self.query_element.has_structure():
            self.pymod.main_window.show_error_message(
                "Selection Error",
                "Foldseek searches with coordinates, and '%s' has none. Load or"
                " predict a structure for it first, or use one of the sequence"
                " searches in this menu."
                % self.element_name(self.query_element))
            return None

        try:
            self.databases = fsc.available_databases()
        except fsc.FoldseekError as e:
            self.pymod.main_window.show_error_message("Foldseek Error", str(e))
            return None

        self.build_options_window()


    def build_options_window(self):
        self.options_window = self.options_window_class_qt(
            parent=self.pymod.main_window,
            protocol=self,
            title="Foldseek Structure Search",
            upper_frame_title="Search the structure databases with %s"
                              % self.element_name(self.query_element),
            submit_command=self.run_search)
        self.options_window.show()


    @staticmethod
    def element_name(element):
        return (getattr(element, "my_header", None)
                or getattr(element, "compact_header", None) or "?")


    def query_chains(self):
        """The chains the object of the query has in PyMOL."""

        try:
            return [chain for chain in cmd.get_chains(
                self.query_element.get_pymol_selector()) if chain]
        except Exception:
            return []


    ###############################################################################################
    # Running the search.                                                                         #
    ###############################################################################################

    def run_search(self):
        """
        Called by the Submit button: writes the chain, sends it, waits, and shows
        what came back.
        """

        window = self.options_window
        databases = window.get_databases()
        if not databases:
            self.pymod.main_window.show_error_message(
                "Foldseek Search", "Choose at least one database to search.")
            return None

        self.query_chain = window.get_chain()
        self.mode = window.get_mode()
        window.destroy()

        try:
            structure_path = self.write_query_structure()
        except Exception as e:
            self.pymod.main_window.show_error_message(
                "Foldseek Search", "The structure could not be written: %s" % e)
            return None

        print("- Foldseek: searching %s with %s, chain %s, against %s"
              % (self.element_name(self.query_element), self.mode,
                 self.query_chain or "the only one", ", ".join(databases)))

        self.databases_searched = list(databases)

        self._error = None
        try:
            if not self.pymod.use_protocol_threads:
                self._search(structure_path, databases)
            else:
                dialog = Protocol_exec_dialog(
                    app=self.pymod.main_window, pymod=self.pymod,
                    function=self._search, args=(structure_path, databases),
                    title="Foldseek Search",
                    # Short: what is going on in detail -- the stage, and how
                    # much of the answer has arrived -- is printed in the shell.
                    label_text=("Searching the selected structural databases"
                                " with Foldseek"),
                    wait_start=0.2, wait_end=0.2)
                dialog.exec_()
        except PyModInterruptedProtocol:
            return None

        if self._error is not None:
            self.pymod.main_window.show_error_message("Foldseek Error",
                                                      str(self._error))
            return None

        if not self.results or not fsr.hit_count(self.results):
            self.pymod.main_window.show_warning_message(
                "Foldseek Search",
                "The search finished and found nothing in the databases which"
                " were chosen.")
            return None

        self.show_results()
        return None


    def _search(self, structure_path, databases):
        """Submit, wait, fetch. Runs off the GUI thread."""

        try:
            self.ticket = fsc.submit(structure_path, databases, mode=self.mode)
            state = fsc.wait(self.ticket, report=self._report_status)
            if state != fsc.STATUS_COMPLETE:
                raise fsc.FoldseekError(self.failed_search_message(state))
            print("  Foldseek: the search is done, fetching its results")
            self.payload = fsc.results(self.ticket, report=self._report_download)
            print("  Foldseek: reading the results")
            self.results = fsr.read_results(self.payload)
            print("  Foldseek: %s hits" % fsr.hit_count(self.results))
        except (fsc.FoldseekError, fsr.FoldseekResultError) as e:
            self._error = e
        except Exception as e:
            self._error = fsc.FoldseekError("The search failed: %s" % e)


    def failed_search_message(self, state):
        """
        What to say when the server ends a search as ERROR.

        Worth saying plainly, because the server names a search by what was
        asked of it: the same chain sent again with the same mode and the same
        databases comes back as the same search, with the same ERROR, however
        many times it is tried. Something has to change for a new search to be
        run at all.
        """

        return ("The search finished as %s rather than completing.\n\n"
                "The server gives the same answer to the same question: this"
                " chain, searched in the same way against the same databases,"
                " is the search which has just failed, so trying it again"
                " unchanged gives this same message. Change something and it"
                " becomes a new search:\n\n"
                "  \u2022 search fewer databases, or other ones -- one database"
                " which cannot answer ends the whole search, not only its own"
                " part of it;\n"
                "  \u2022 choose another alignment (3Di and amino acids, TM-align,"
                " 3Di only);\n"
                "  \u2022 search another chain, if the one chosen is very short"
                " or has few residues with coordinates.\n\n"
                "The search was %s, chain %s, against %s.\n"
                "Its ticket is %s, and it can be looked at on %s."
                % (state, self.mode, self.query_chain or "the only one",
                   ", ".join(self.databases_searched or []), self.ticket,
                   fsc.DEFAULT_BASE_URL))


    def _report_download(self, megabytes):
        """
        Said while the results are being fetched. A search of many databases
        answers with tens of megabytes, and on a slow day the download is the
        longest part of the whole thing -- with nothing to show for it, the
        window which says "please wait" looks stuck.
        """

        print("  Foldseek: %.0f MB of results so far" % megabytes)


    def _report_status(self, state):
        if state != getattr(self, "_last_status", None):
            print("  Foldseek: %s" % state.lower())
            self._last_status = state


    def write_query_structure(self):
        """
        Writes the chain which is searched, and returns its path.

        A PDB written by PyMOL is what the server is given: it accepts it, which
        was measured, so nothing has to be converted.
        """

        directory = self.pymod.temp_directory_dirpath
        path = os.path.join(directory, "foldseek_query.pdb")

        selector = self.query_element.get_pymol_selector()
        if self.query_chain:
            selection = "(%s) and chain %s and polymer" % (selector, self.query_chain)
        else:
            selection = "(%s) and polymer" % selector

        if not cmd.count_atoms(selection):
            raise ValueError("there are no atoms in %s" % selection)
        cmd.save(path, selection)
        return path


    ###############################################################################################
    # The results.                                                                                #
    ###############################################################################################

    def show_results(self):
        self.results_window = self.results_window_class_qt(
            parent=self.pymod.main_window, protocol=self)
        self.results_window.show()


    def import_hits(self, hits, import_mode="single-chain"):
        """
        Builds an element for each chosen hit, gathers them in a cluster with the
        query, and fetches the structure of those which have one.

        Nothing here is special to Foldseek: the elements are named the way PyMod
        names elements, and the structures are fetched and associated by the
        protocols which serve the other searches.
        """

        if not hits:
            return None

        query = self.query_element
        query_index = self.pymod.get_pymod_element_index_in_container(query)

        elements = []
        for hit in hits:
            sequence = hit.sequence or hit.hit_alignment.replace("-", "")
            if not sequence:
                continue
            element = self.pymod.build_pymod_element_from_hsp(
                sequence, hit.pymod_header)
            self.pymod.add_element_to_pymod(element)
            element.foldseek_hit = hit
            elements.append((element, hit))

        if not elements:
            self.pymod.main_window.show_warning_message(
                "Foldseek Search", "None of the chosen hits carried a sequence.")
            return None

        cluster = self.pymod.add_new_cluster_to_pymod(
            cluster_type="blast-cluster",
            cluster_label="Foldseek",
            query=query,
            child_elements=[query] + [element for element, _ in elements],
            algorithm="foldseek",
            update_stars=False)
        self.pymod.change_pymod_element_list_index(cluster, query_index)

        self.align_to_query(query, elements)

        # The structures, through the ordinary protocols of PyMod. What is in
        # PyMOL before they are fetched is noted, because what is there
        # afterwards and was not before is exactly what has just been loaded.
        objects_before = set(cmd.get_object_list())

        self.fetch_structures([element for element, hit in elements
                               if hit.structure_source == "pdb"], "pdb",
                              import_mode)
        self.fetch_structures([element for element, hit in elements
                               if hit.structure_source == "afdb"], "afdb",
                              import_mode)

        # And on top of the query, where the search found them.
        self.superpose_to_query(objects_before)

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True,
                                       update_elements=True, update_menus=True)
        return None


    def align_to_query(self, query, elements):
        """
        Puts the hits in register with the query, using the alignment the search
        already made -- the same centre star construction the BLAST results use.
        """

        query_sequence = query.my_sequence.replace("-", "")
        sequences = [query_sequence]
        identifiers = [query.get_unique_index_header()]
        pairs = []
        to_update = {query.get_unique_index_header(): query}

        for element, hit in elements:
            aligned_query = hit.query_alignment
            aligned_hit = hit.hit_alignment
            if not aligned_query or not aligned_hit:
                continue

            start = (hit.query_start or 1) - 1
            length = len(aligned_query.replace("-", ""))
            updated_query = (query_sequence[:start] + aligned_query
                             + query_sequence[start + length:])
            updated_hit = ("-" * len(query_sequence[:start]) + aligned_hit
                           + "-" * len(query_sequence[start + length:]))

            if updated_query.replace("-", "") != query_sequence:
                # The alignment does not describe this query: leave the hit
                # unaligned rather than corrupting the cluster.
                continue

            pairs.append((updated_query, updated_hit))
            sequences.append(aligned_hit.replace("-", ""))
            identifiers.append(element.get_unique_index_header())
            to_update[element.get_unique_index_header()] = element

        if len(sequences) < 2:
            return None

        pairs.extend([None] * int((len(sequences) - 1) * (len(sequences) - 2) / 2))
        records = save_cstar_alignment(seqs=sequences, all_ids=identifiers,
                                       pairwise_alis=pairs)
        for record in records:
            if record.id in to_update:
                to_update[record.id].my_sequence = str(record.seq)
        return None


    def fetch_structures(self, elements, source, import_mode):
        """
        Fetches the structures of these elements and associates them, each from
        the database it came from.
        """

        if not elements:
            return None

        try:
            if source == "pdb":
                # The protocol which fetches from the PDB takes the elements and
                # the mode directly; its initialize_from_gui() would read the
                # current selection of PyMod instead, which is not what was
                # ticked in the results window.
                fetcher = Fetch_structure_file(self.pymod)
                fetcher.structures_to_fetch = list(elements)
                fetcher.import_mode = import_mode
                fetcher.fetch_pdb_files()
            else:
                self.fetch_models(elements, import_mode)
        except Exception as e:
            self.pymod.main_window.show_error_message(
                "Structure Download Error",
                "The sequences were imported, but their structures could not be"
                " fetched: %s" % e)
        return None


    def fetch_models(self, elements, import_mode):
        """
        Downloads the AlphaFold model of each hit and associates it, one at a
        time, as PyMod associates every other structure.

        Fetch_models_file is not used for this. It does the same download and the
        same association, and then writes what it did into the session of the
        AlphaFold command -- pymod.af_modeling_session_list[pymod.current_af_session]
        -- which exists only while that command is running. Driven from here it
        ended in "'PyMod' object has no attribute 'current_af_session'", reported
        as a failed association although the structure had in fact been
        associated a moment earlier.
        """

        failures = []
        for element in elements:
            hit = getattr(element, "foldseek_hit", None)
            accession = ((hit.uniprot if hit is not None else None)
                         or str(getattr(element, "model_id", "") or ""))
            if not accession:
                failures.append((self.element_name(element),
                                 "no accession to ask the database for"))
                continue

            try:
                path = fetch_models_file(
                    model_id=accession,
                    output_dir=self.pymod.temp_directory_dirpath)
            except Exception as e:
                failures.append((self.element_name(element),
                                 "the model could not be downloaded: %s" % e))
                continue

            if not path:
                failures.append((self.element_name(element),
                                 "the AlphaFold database has no model for %s"
                                 % accession))
                continue

            try:
                associate_all = import_mode == "multiple-chains"
                associator = Associate_model_structure(self.pymod, element,
                                                       ASSOCIATE_ALL=associate_all)
                # A model of the AlphaFold database is one chain, and it is A.
                associator.associate(path, "A")
            except Exception as e:
                failures.append((self.element_name(element),
                                 "the structure could not be associated: %s" % e))

        if failures:
            self.pymod.main_window.show_warning_message(
                "AlphaFold Models",
                "The sequences were imported. The model of %s could not be"
                " loaded:\n\n%s"
                % ("one hit" if len(failures) == 1 else "%s hits" % len(failures),
                   "\n".join("%s: %s" % pair for pair in failures[:6])))
        return None


    def superpose_to_query(self, objects_before):
        """
        Puts every structure which was just loaded on top of the query in PyMOL.

        The alignment which the search returned puts the sequences in register;
        this is the other half of it, and the half the eye reads. Without it the
        hits sit wherever the PDB and the AlphaFold database left them, which for
        a list of them is nowhere near each other.

        What is superposed is read from PyMOL -- the objects which were not there
        before the structures were fetched -- and not from the elements this
        protocol built. Associating a structure does not fill in the element: it
        builds a new one carrying the structure and puts it in the place of the
        old one (pymod.replace_element, and in 'multiple-chains' mode the old
        element is deleted outright), so by the time there is anything to
        superpose, the elements which were passed here no longer exist. Looking
        for them was the reason nothing moved.

        'super' is used rather than 'align' because it superposes on the fold and
        not on the sequence -- which is the whole point of this search: a hit may
        resemble the query in shape while sharing almost nothing of its sequence,
        and a sequence alignment of the two would have little to work with.
        """

        query = self.query_element
        if not query.has_structure():
            return None

        query_object = self.object_of(query.get_pymol_selector())
        fixed = query.get_pymol_selector()
        if self.query_chain:
            fixed = "(%s) and chain %s" % (fixed, self.query_chain)

        arrived = [name for name in cmd.get_object_list()
                   if name not in objects_before and name != query_object]
        if not arrived:
            return None

        # The query does not move: everything is moved onto it. It is left out
        # of what is superposed, and where it stands is written down here so
        # that it can be put back if anything moves it all the same -- a hit
        # which is the same entry as the query, loaded into the same object, is
        # the way that happens.
        standing = self.where_it_stands(query_object)

        superposed, failures, put_back = [], [], False
        for group in self.group_by_entry(arrived):

            # A structure imported with all its chains has each chain in an
            # object of its own. The one which is superposed is the one which has
            # most to say about the query -- the chain the search matched -- and
            # the rest of the entry is moved by the same matrix, or superposing
            # one chain of a complex would pull the complex apart.
            mobile, outcome = self.best_of_the_entry(group, fixed)
            if mobile is None:
                failures.extend(group)
                continue

            try:
                cmd.super(mobile, fixed)
            except Exception:
                failures.extend(group)
                continue

            for other in group:
                if other != mobile:
                    try:
                        cmd.matrix_copy(mobile, other)
                    except Exception:
                        pass
            # Before anything is written down, and before the next hit is
            # superposed: every hit is put on the query where it stands, not on
            # a query which the hit before it has nudged.
            put_back = self.put_it_back(query_object, standing) or put_back

            for name in group:
                self.save_superposed(name)

            superposed.append((mobile, outcome[0], outcome[1],
                               len(group) - 1))

        for name, rmsd, atoms, followers in superposed:
            print("- Foldseek: %s superposed on the query, RMSD %.2f A over %s"
                  " atoms%s" % (name, rmsd, atoms,
                                ", with %s other chains of its structure"
                                % followers if followers else ""))
        if failures:
            print("- Foldseek: nothing could be superposed for %s"
                  % ", ".join(failures[:6]))

        # What happened to the query is said out loud, every time. It is the one
        # thing about this which is not visible on the screen -- a structure
        # which has moved looks exactly like one which has not -- and it is the
        # thing which must not happen.
        moved_by = self.how_far_it_moved(query_object, standing)
        if self.put_it_back(query_object, standing) or put_back:
            print("- Foldseek: %s is the structure which was searched with, and"
                  " it has been put back where it stood (it had been moved %s):"
                  " the hits are moved onto it, never it onto them."
                  % (query_object, "%.2f A" % moved_by
                     if moved_by is not None else "by a matrix"))
        else:
            print("- Foldseek: %s, the structure which was searched with, has not"
                  " moved (%s)."
                  % (query_object, "%.3f A" % moved_by if moved_by is not None
                     else "its coordinates could not be read"))
        return None


    @staticmethod
    def how_far_it_moved(object_name, standing):
        """How far an object is from where it stood, in Angstroms, or None."""

        if standing is None or standing[0] is None:
            return None
        try:
            now = cmd.get_coords(object_name, 1)
            if now is None or now.shape != standing[0].shape:
                return None
            return float(abs(now - standing[0]).max())
        except Exception:
            return None


    @staticmethod
    def where_it_stands(object_name):
        """
        Where an object is now: its coordinates and the matrix it is drawn with.
        Both are needed, because a superposition can move either.
        """

        try:
            return (cmd.get_coords(object_name, 1),
                    cmd.get_object_matrix(object_name))
        except Exception:
            return None


    def put_it_back(self, object_name, standing):
        """
        Puts an object back where it stood, and says whether it had to. Nothing
        should move the query; this is what makes sure of it rather than trusting
        that nothing does.
        """

        if standing is None:
            return False
        coordinates, matrix = standing
        now = self.where_it_stands(object_name)
        if now is None:
            return False

        moved = False
        try:
            if coordinates is not None and now[0] is not None:
                if abs(coordinates - now[0]).max() > 1e-4:
                    cmd.load_coords(coordinates, object_name, 1)
                    moved = True
            if matrix is not None and now[1] is not None:
                if max(abs(a - b) for a, b in zip(matrix, now[1])) > 1e-6:
                    cmd.set_object_ttt(object_name, matrix)
                    moved = True
        except Exception:
            return False
        return moved


    @staticmethod
    def object_of(selector):
        """The name of the PyMOL object a selector of PyMod points at."""

        return str(selector).split(" and ")[0].strip("() ")


    def group_by_entry(self, object_names):
        """
        The objects which came out of the same structure file, grouped together:
        the chains of one entry belong to each other and must move together.

        The file each object came from is read from the PyMod element which owns
        it. An object with no element of its own stands alone.
        """

        elements = {}
        for element in self.pymod.get_pymod_elements_list():
            if element.has_structure():
                elements[self.object_of(element.get_pymol_selector())] = element

        groups, by_entry = [], {}
        for name in object_names:
            element = elements.get(name)
            entry = getattr(getattr(element, "structure", None),
                            "original_structure_file_path", None)
            if not entry:
                groups.append([name])
                continue
            by_entry.setdefault(entry, []).append(name)
        return groups + list(by_entry.values())


    @staticmethod
    def best_of_the_entry(group, fixed):
        """
        Which chain of an entry to superpose, and what its superposition is.

        'super' is asked what it would do to each of them without doing it, and
        the one which matches most of the query wins. This is the chain the
        search found; asking rather than guessing means the answer does not
        depend on how the chains happen to be named.
        """

        best_name, best_outcome, best_score = None, None, None
        for name in group:
            try:
                outcome = cmd.super(name, fixed, transform=0)
            except Exception:
                continue
            if not outcome or not outcome[1]:
                continue
            score = (outcome[1], -outcome[0])  # most atoms, then lowest RMSD.
            if best_score is None or score > best_score:
                best_name, best_outcome, best_score = name, outcome, score
        return best_name, best_outcome


    def save_superposed(self, object_name):
        """
        Writes the moved coordinates back to the file PyMod keeps for this chain,
        so that the structure is where it was put after the project is saved and
        opened again. This is what the structural alignments of PyMod do too.
        """

        for element in self.pymod.get_pymod_elements_list():
            if not element.has_structure():
                continue
            if self.object_of(element.get_pymol_selector()) != object_name:
                continue
            path = getattr(element.structure, "current_chain_file_path", None)
            if path:
                try:
                    cmd.save(path, object_name)
                except Exception:
                    pass
            return


    def show_backbone(self, hit):
        """
        Draws the backbone of a hit in PyMOL from the coordinates the search
        returned, without fetching anything, and superposes it on the query.
        """

        text = fsr.as_pdb_string(hit)
        if not text:
            self.pymod.main_window.show_warning_message(
                "Foldseek", "This hit came back without coordinates.")
            return None

        name = "foldseek_%s" % hit.pymod_header.replace("-", "_")
        try:
            cmd.delete(name)
            cmd.read_pdbstr(text, name)
            cmd.hide("everything", name)
            cmd.show("ribbon", name)
            cmd.color("salmon", name)
            try:
                cmd.align("%s////CA" % name,
                          "(%s)////CA" % self.query_element.get_pymol_selector())
            except Exception:
                pass
        except Exception as e:
            self.pymod.main_window.show_error_message(
                "Foldseek", "The backbone could not be shown: %s" % e)
        return None
