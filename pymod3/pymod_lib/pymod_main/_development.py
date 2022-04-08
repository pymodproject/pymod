"""
Module for development.
"""

import os
import sys
import urllib.request
import random
import shutil
import importlib

from pymod_lib.pymod_seq.seq_star_alignment import global_pairwise_alignment_cs


class PyMod_development:
    """
    Methods used when developing or testing PyMod.
    """

    def _launch_develop(self):
        """
        For development only. The 'open_sequence_file', 'open_structure_file' and
        'build_cluster_from_alignment_file' methods are used to import sequences
        when PyMod starts.
        The code in this method illustrates how to the API of PyMod to developers.
        """

        return None

        #------------------------------------------
        # Set up the paths of some example files. -
        #------------------------------------------

        spec = importlib.util.find_spec("pymod_lib")
        if spec is None:
            raise Exception("'pymod_lib' module not found.")
        pymod_lib_data_dirpath = os.path.join(os.path.dirname(spec.origin), "pymod_data")

        alignment_filepath = os.path.join(pymod_lib_data_dirpath, "sequences", "PF00069_seed_min.fasta")
        sequence_filepath = os.path.join(pymod_lib_data_dirpath, "sequences", "sequence.fasta")
        pdb_filepath = os.path.join(pymod_lib_data_dirpath, "pdb", "1NDD.pdb")

        #------------------------------
        # Examples of sequence files. -
        #------------------------------

        # Loads every sequence present in a FASTA file in PyMod.
        elements = self.open_sequence_file(sequence_filepath)
        # Selects in PyMod the first element.
        elements[0].widget_group.toggle_element()

        # Loads a multiple sequence alignment in PyMod by putting all the
        # sequences in an alignment object.
        cluster = self.build_cluster_from_alignment_file(alignment_filepath)

        # Add to the alignment object the previously loaded elements.
        add_child = False
        if add_child:
            for e in elements:
                cluster.add_child(e)


        #-------------------------------
        # Example of a structure file. -
        #-------------------------------

        # Loads in PyMod/PyMOL all the chains present in a PDB file.
        chain_elements = self.open_structure_file(pdb_filepath)


        #-------------------------------
        # Build elements from scratch. -
        #-------------------------------

        # Takes the sequence of the first PDB chain loaded before, and initializes
        # a new PyMod element with a randomized copy of it.
        template_element = chain_elements[0]
        randomized_seq = self.randomize_sequence(template_element.my_sequence, n_gaps=1, n_insertions=1)
        new_element = self.build_pymod_element_from_args("test_sequence", randomized_seq)
        self.add_element_to_pymod(new_element)

        # Add a new alignment object to PyMod in which contains the randomized
        # sequence and the original PDB chain element.
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="alignment",
                                                    child_elements=[new_element, template_element],
                                                    algorithm="imported")

        # Perform a sequence alignment between the two elements and updated their
        # sequence in PyMod.
        sq, st = global_pairwise_alignment_cs(new_element.my_sequence, template_element.my_sequence, 20, 5)
        new_element.my_sequence = sq
        template_element.my_sequence = st

        #-----------------------
        # Launching protocols. -
        #-----------------------

        # Load sessions.
        if 0:
            # Loads a PyMod session from the API. Useful when developing some
            # complex functionality which requires a lot of steps to test.
            project_archive_filepath = "/home/user/pymod_session.pmse"
            self.open_pymod_session(project_archive_filepath)

        # Automatically launches a protocol (input sequences must be correctly selected
        # in the code above).
        if 0:
            # Homology modeling.
            from pymod_lib.pymod_protocols.modeling_protocols.homology_modeling import MODELLER_homology_modeling
            modeller_session = MODELLER_homology_modeling(self)
            modeller_session.launch_from_gui()

        if 0:
            # Contact map.
            from pymod_lib.pymod_protocols.structural_analysis_protocols.contact_map_analysis import Contact_map_analysis
            Contact_map_analysis(self).launch_from_gui()

        if 0:
            # Multiple alignments.
            from pymod_lib.pymod_protocols.alignment_protocols.clustalo import Clustalomega_regular_alignment
            Clustalomega_regular_alignment(self).launch_from_gui()

        if 0:
            # BLAST.
            from pymod_lib.pymod_protocols.similarity_searches_protocols.psiblast import PSI_BLAST_search
            PSI_BLAST_search(self).launch_from_gui()

        if 0:
            # HMMER.
            from pymod_lib.pymod_protocols.similarity_searches_protocols.phmmer import PHMMER_search
            PHMMER_search(self).launch_from_gui()

        if 0:
            # HMMSCAN
            from pymod_lib.pymod_protocols.domain_analysis_protocols.domain_analysis import Domain_Analysis_Protocol
            Domain_Analysis_Protocol(self, "local").launch_from_gui()


        self.main_window.gridder(update_clusters=True, update_menus=True, update_elements=True)


    def load_uniprot_random(self, reviewed=False, grid=True):
        try:
            if reviewed:
                rev_string = "yes"
            else:
                rev_string = "no"
            temp_fasta_path = urllib.request.urlretrieve("http://www.uniprot.org/uniprot/?query=reviewed:%s+AND+organism:9606&random=yes&format=fasta" % rev_string)[0]
            self.open_sequence_file(temp_fasta_path)
            if grid:
                self.main_window.gridder(update_clusters=True, update_menus=True, update_elements=True)
        except Exception as e:
            if grid:
                self.main_window.show_error_message("UniProt Error",
                    "Could not obtain a sequence from the UniProt server beacuse of the following reason: %s" % str(e))


    def load_pdb_random(self, code=None, grid=True):
        try:
            spec = importlib.util.find_spec("pymod_lib")
            if spec is None:
                raise Exception("'pymod_lib' module not found.")
            pdb_list_filepath = os.path.join(os.path.dirname(spec.origin), "pymod_data", "pdb", "all_proteins.txt")

            if code is None:
                with open(pdb_list_filepath, "r") as p_fh:
                    ids = [i.replace(" ", "") for i in p_fh.read().split("\n")]
                    code = ids[random.randint(0, len(ids)-1)]

            print("\n# Fetching PDB: %s" % code)
            file_url = "https://files.rcsb.org/download/%s.pdb" % code
            temp_path = urllib.request.urlretrieve(file_url)[0]
            temp_pdb_path = os.path.join(os.path.dirname(temp_path), "%s.pdb" % code)
            shutil.move(temp_path, temp_pdb_path)

            elements = self.open_structure_file(temp_pdb_path)
            if grid:
                self.main_window.gridder(update_clusters=True, update_menus=True, update_elements=True)
            return elements

        except Exception as e:
            if grid:
                self.main_window.show_error_message("RCSB PDB Error",
                    "Could not obtain a sequence from the RCSB server beacuse of the following reason: %s" % str(e))


    def randomize_sequence(self, sequence, seqid=30.0, n_gaps=0, n_insertions=2):
        amino_acids = "QWERTYPASDFGHKLCVNM" # + "X"
        list_seq = list(sequence)

        n_substitutions = int(seqid*len(list_seq)/100.0)
        for s in range(0, n_substitutions):
            seq_len = len(list_seq)
            random_index = random.randint(0, seq_len-1)
            list_seq.pop(random_index)
            list_seq.insert(random_index, random.choice(amino_acids))


        max_gap_length = 10
        for g in range(0, n_gaps):
            seq_len = len(list_seq)
            random_index = random.randint(0, seq_len-1)
            gap_length = random.randint(0, max_gap_length)
            for l in range(0, gap_length):
                list_seq.pop(random_index-l)


        max_insertion_length = 18
        for i in range(0, n_gaps):
            seq_len = len(list_seq)
            random_index = random.randint(0, seq_len-1)
            insertion_length = random.randint(0, max_insertion_length)
            for l in range(0, insertion_length):
                list_seq.insert(random_index, random.choice(amino_acids))

        return "".join(list_seq)
