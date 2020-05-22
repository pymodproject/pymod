# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing classes to represent "PyMod elements", that is each sequence or alignment loaded
in the PyMod main window.
"""

import os

from . import pymod_vars as pmdt
from .pymod_seq import seq_manipulation
from .pymod_seq import seq_conservation
from .pymod_seq.seq_headers import check_fetchable_pdb_from_header
from .pymod_residue import PyMod_standard_residue
from .pymod_exceptions import PyModMissingStructure, PyModSequenceConflict


###################################################################################################
# PYMOD ELEMENTS.                                                                                 #
###################################################################################################

class PyMod_element:
    """
    A base that stores all the informations of a sequence or sequence cluster.
    """
    cluster = False

    def __init__(self,
                 header=None,
                 description=None,
                 color="white"
                 ):

        #----------------------------------------
        # Indices and PyMod_element properties. -
        #----------------------------------------

        # It is a unique id to identify each element. It is given by the "unique_index" of the
        # "pymod" object. This values will be assigned when the element is added to the list of
        # PyMod_element objects by the '.add_element_to_pymod' method of the PyMod class.
        self.unique_index = None

        self.mother = None
        self.list_of_children = []

        self.blast_query = False
        self.lead = False
        self.bridge = False

        #--------------------------
        # Sequences and residues. -
        #--------------------------

        self.my_sequence = None

        #--------------------------------------------------------------------------------
        # Headers (name of the sequences and clusters displayed throughout the plugin). -
        #--------------------------------------------------------------------------------

        # The full original header.
        self.original_header = header

        self.my_header = header
        self.my_header_root = header

        self.compact_header = header
        self.compact_header_root = header
        self.compact_header_prefix = ""


        # Check if a PDB can be fetched for this element from its name.
        can_be_fetched_results = check_fetchable_pdb_from_header(self.original_header)
        if can_be_fetched_results == None:
            self.pdb_can_be_fetched = False
            self.pdb_id = None
            self.pdb_chain = None
        else:
            self.pdb_can_be_fetched = True
            self.pdb_id = can_be_fetched_results[0]
            self.pdb_chain = can_be_fetched_results[1]


        # Sets the 'my_sequence' attribute. The primary sequence of an element. If the sequence is
        # changed or aligned by the user, it will be modified to include indels.
        # self.set_sequence(record_seq, adjust_sequence)

        #--------------------------------
        # Descriptions and annotations. -
        #--------------------------------

        self.description = description
        self.annotations = {}

        #----------------------------------------------
        # Appearance and intercations with the users. -
        #----------------------------------------------

        # Its value is False when the sequence header is not selected, it becomes True when the user
        # selects the sequence by left-clicking on its header.
        self.selected = False

        # Defines the way the sequence has to be colored.
        self.color_by = "regular"
        self._color_by = "regular"

        # Name of the color of the sequence when its "color_by" attribute is set to "regular".
        self.my_color = color

        #-------------------
        # Alignment files. -
        #-------------------

        self.tree_file_path = None
        self.tree_file_name = None


    #################################################################
    # Methods for managing PyMod clusters.                          #
    #################################################################

    def __repr__(self):
        return self.my_header + ' ' + self.__class__.__name__


    def is_cluster(self):
        return self.cluster


    def is_mother(self):
        if self.is_cluster() and self.list_of_children:
            return True
        else:
            return False


    def is_child(self, exclude_root_element=True):
        if self.mother:
            if exclude_root_element and self.is_root_child():
                return False
            else:
                return True
        else:
            return False


    def is_root(self):
        return isinstance(self, PyMod_root_element)


    def is_root_child(self):
        return self.mother.is_root()


    def is_root_sequence(self):
        return self.is_root_child() and not self.is_cluster()


    def get_ancestor(self, exclude_root_element=True):
        if self.is_child(exclude_root_element):
            if self.mother.get_ancestor(exclude_root_element):
                return self.mother.get_ancestor(exclude_root_element)
            else:
                return self.mother
        return None


    def get_descendants(self):
        descendants = []
        for d in self.get_children():
            descendants.append(d)
            if d.is_mother():
                descendants.extend(d.get_descendants())
        return descendants


    def filter_for_sequences(method):
        def filterer(self, sequences_only=False):
            if sequences_only:
                return [e for e in method(self, sequences_only) if not e.is_cluster()]
            else:
                return method(self, sequences_only)
        return filterer


    @filter_for_sequences
    def get_siblings(self, sequences_only=False):
        if not self.mother:
            return []
        else:
            return [c for c in self.mother.list_of_children if c != self]


    def extract_to_upper_level(self):
        """
        Extracts an element to the upper level. Overridden by subclasses.
        """
        new_mother = self.mother.mother
        self.remove_all_lead_statuses()
        new_mother.add_child(self)


    def delete(self):
        mother = self.mother
        mother.remove_child(self)


    #################################################################
    # Check element properties.                                     #
    #################################################################

    def has_structure(self):
        """
        This will be overridden in 'PyMod_sequence_element' classes.
        """
        return False


    #################################################################
    # Cluster leads.                                                #
    #################################################################

    def set_as_lead(self):
        self.remove_all_lead_statuses()
        self.lead = True

    def set_as_query(self):
        self.remove_all_lead_statuses()
        self.blast_query = True

    def remove_all_lead_statuses(self):
        self.lead = False
        self.blast_query = False


    def is_blast_query(self):
        return self.blast_query

    def is_lead(self):
        return self.lead or self.blast_query

    def is_bridge(self):
        return self.bridge


    #################################################################
    # Headers formatting.                                           #
    #################################################################

    def get_unique_index_header(self):
        return pmdt.unique_index_header_formatted % self.unique_index


    #################################################################
    # Colors.                                                       #
    #################################################################

    def store_current_color(self):
        self._color_by = self.color_by

    def revert_original_color(self):
        self.color_by = self._color_by


###################################################################################################
# CLUSTERS.                                                                                       #
###################################################################################################

class _PyMod_cluster_element(PyMod_element):
    """
    Cluster elements represent the 'alignment' elements in the PyMod main window. They take one row
    in the main window exactly like sequence elements and behave similarly.
    """

    cluster = True

    def __init__(self, sequence=None, header=None, algorithm=None, cluster_type="generic", cluster_id=None, description=None, color=""):
        PyMod_element.__init__(self, header, description=description, color=color)
        self.algorithm = algorithm
        self.cluster_type = cluster_type
        self.cluster_id = cluster_id
        self.initial_number_of_sequences = None
        self.my_sequence = sequence
        self.rmsd_dict = None

    def add_children(self, children): # TODO: use the 'set_initial_number_of_sequences' argument.
        if not hasattr(children, "__iter__"):
            children = [children]
        for child in children:
            self.add_child(child)
        if self.initial_number_of_sequences == None:
            self.initial_number_of_sequences = len(children)

    def add_child(self, child):
        # Remove from the old mother.
        if child.is_child(exclude_root_element=False):
            old_mother = child.mother
            old_mother.remove_child(child)
        # Add to new mother.
        child.mother = self
        self.list_of_children.append(child)


    def remove_child(self, child):
        self.list_of_children.remove(child)
        child.remove_all_lead_statuses()

    def get_children(self):
        return self.list_of_children

    def get_lead(self):
        for child in self.get_children():
            if child.is_lead():
                return child
        return None


    def update_stars(self, adjust_elements=False):
        self.my_sequence = seq_conservation.compute_stars(self.get_children(), adjust_elements=adjust_elements)


    def remove_gap_only_columns(self):
        seq_manipulation.remove_gap_only_columns(self.get_children())

    def adjust_aligned_children_length(self):
        seq_manipulation.adjust_aligned_elements_length(self.get_children())


    #################################################################
    # Alignment files.                                              #
    #################################################################

    def get_tree_file_path(self):
        return self.tree_file_path

    def set_tree_file_path(self, tree_file_path):
        self.tree_file_path = tree_file_path
        self.tree_file_name = os.path.basename(tree_file_path)

    def update_alignment_files(self, alignments_dirpath):
        if self.tree_file_path is not None and self.tree_file_name is not None:
            self.tree_file_path = os.path.join(alignments_dirpath, self.tree_file_name)


###################################################################################################
# SEQUENCES.                                                                                      #
###################################################################################################

class _PyMod_sequence_element(PyMod_element):
    """
    The objects of this class are the sequences (both from sequences and structure files) that
    appear on the left column or alignments elements.
    """

    def __init__(self, sequence=None, header="", residues=None, description=None, color="white"):

        PyMod_element.__init__(self, header=header, description=description, color=color)

        # TODO: cleans up the sequence.
        pass

        #----------------------------------------
        # Builds the residues and the sequence. -
        #----------------------------------------

        self.my_sequence = sequence
        self.residues = residues
        self._polymer_residues = None

        if self.residues == None and self.my_sequence == None:
            raise Exception("Either a sequence or a residues list must be provided.")
        elif self.residues != None and self.my_sequence != None:
            raise Exception("Can not accept both a sequence and a residues list.")
        # If the residues list is not provided, build it through the sequence provided in the
        # constructor.
        elif self.residues == None and self.my_sequence != None:
            self.set_residues_from_sequence()
        elif self.residues != None and self.my_sequence == None:
            self.set_sequence_from_residues()

        # Update residues information with indices.
        self.update_residues_information()

        #----------------------------------------------------
        # Other sequence- or structure-related information. -
        #----------------------------------------------------

        self.initialize_additional_information()


    def initialize_additional_information(self):
        # A PyMod element with an associated structure will have a 'PyMod_structure' object assigned
        # to this attribute.
        self.structure = None
        # Boolean flags to define if a sequence has some secondary structure assigned to it.
        self.assigned_secondary_structure = False
        self.predicted_secondary_structure = False

        # Store various type of data for an element.
        self._has_campo_scores = False
        self._has_entropy_scores = False
        self._has_pc_scores = False
        self._has_dope_scores = False

        # Contain the features (e.g.: domains) of an element.
        self.features = dict([(ft, []) for ft in pmdt.feature_types])
        self.features_count = 0
        self.features_selectors_list = []

        # Stores the id of the domain analysis session in which an element was analysed.
        self.domain_analysis_id = None

        # When a query element is split in its domains, the derived PyMod elements will be stored
        # in this list.
        self.derived_domains_list = []

        # When a query element is split in its domains, the PyMod element derived from the
        # splitting will receive the corresponding 'Domain_feature' object which will be assigned to
        # this attribute.
        self.splitted_domain = None

        # Checks if the element is a protein or a nucleic acid by checking the names of its residues.
        self.polymer_type = "protein"
        nucleotides = [nt for nt in list(pmdt.nucleic_acids_dictionary.keys())]
        for r in self.residues:
            if r.is_standard_residue() and r.three_letter_code in nucleotides:
                self.polymer_type = "nucleic_acid"
                break

        self.color_by = "regular"


    def add_feature(self, new_feature):
        """
        Used when adding a feature to a sequence.
        """
        self.features[new_feature.feature_type].append(new_feature)
        feature_start = int(new_feature.start)
        feature_end = int(new_feature.end)
        for r_index, res in enumerate(self.get_polymer_residues()):
            if feature_start <= r_index <= feature_end:
                res.features[new_feature.feature_type].append(new_feature)
        new_feature.element = self
        self.features_count += 1


    def add_domain_feature(self, domain_feature):
        """
        Used when adding a domain feature to a sequence.
        """

        self.features["domains"].append(domain_feature)
        domain_start = int(domain_feature.start)
        domain_end = int(domain_feature.end)
        for r_index, res in enumerate(self.get_polymer_residues()):
            if domain_start <= r_index < domain_end:
                res.features["domains"].append(domain_feature)
        domain_feature.element = self

        if domain_feature.is_derived:
            self.splitted_domain = domain_feature


    def get_domains_features(self):
        return [feat for feat in self.features["domains"] if not feat.is_derived]


    def clear_domains(self):
        self.features["domains"] = []
        for r in self.residues:
            r.features["domains"] = []
        self.splitted_domain = None
        self.derived_domains_list = []


    def clear_features(self):
        self.features["sequence"] = []
        for r in self.residues:
            r.features["sequence"] = []
        self.features_selectors_list = []


    def clear_derived_domains(self):
        self.derived_domains_list = []


    def set_residues_from_sequence(self):
        self.residues = []
        self._polymer_residues = None
        for letter in self.my_sequence:
            if letter != "-":
                self.residues.append(PyMod_standard_residue(one_letter_code=letter,
                                                            three_letter_code=pmdt.get_prot_one_to_three(letter)))


    def set_sequence_from_residues(self):
        my_sequence = ""
        for res in self.residues:
            if res.is_polymer_residue():
                my_sequence += res.one_letter_code
        self.my_sequence = my_sequence


    def update_residues_information(self):
        polymer_residue_count = 0
        for i,res in enumerate(self.residues):
            res.index = i # Index considering also heteroresidues in the sequences.
            res.seq_index =  polymer_residue_count # Index considering only the polymer sequence.
            if res.is_polymer_residue():
                polymer_residue_count += 1
            if res.db_index == None:
                res.db_index = i + 1
            res.pymod_element = self


    ###############################################################################################
    # Sequence related.                                                                           #
    ###############################################################################################

    def set_sequence(self, new_sequence, permissive=False):
        """
        Updates the sequence of a PyMod element with a string provided in the 'new_sequence' argument.
        The gapless versions of the new and old sequences will always be compared. If 'permissive'
        is set to 'False' and if the old and new sequences do not match, an exception will be raised.
        If 'permissive' is set to 'True', the new sequence will overwrite the old one and all the
        attributes and features of the PyMod element (residue scores, domains, etc...) will be lost.
        """

        current_sequence_ungapped = self.my_sequence.replace("-","")
        new_sequence_ungapped = new_sequence.replace("-","")

        if new_sequence_ungapped != current_sequence_ungapped:
            if not permissive:
                print("# Sequence:", self.my_header)
                print("# New:", new_sequence_ungapped)
                print("# Old:", current_sequence_ungapped)
                raise PyModSequenceConflict("The new sequence does not match with the previous one.")
            else:
                self.set_new_sequence(new_sequence)
        else:
            self.my_sequence = new_sequence


    def set_new_sequence(self, new_sequence):
        """
        Overwrites the old sequence with a new one. Also sets up again the residues of the element,
        thus removing any previous information.
        """
        self.my_sequence = new_sequence
        self.set_residues_from_sequence()
        self.update_residues_information()
        self.initialize_additional_information()


    def trackback_sequence(self, sequence_to_align):
        ali = seq_manipulation.global_pairwise_alignment(self.my_sequence.replace("-", ""), sequence_to_align)
        self.set_sequence(ali["seq1"])
        return ali["seq1"], ali["seq2"]


    ###############################################################################################
    # Residues related.                                                                           #
    ###############################################################################################

    def get_polymer_residues(self):
        """
        Get the "polymer residues" of a PyMod element. These residues are the ones actually showed
        in the sequence loaded in PyMod main window an comprise standard amino acid residues and
        modified residues (such as a phospho-serine), represented by an X.
        """
        if self._polymer_residues is None: # Stores the list, for fast access when doing multiple calls.
            self._polymer_residues = [r for r in self.residues if r.is_polymer_residue()]
        return self._polymer_residues

    def get_polymer_sequence_string(self):
        return "".join([res.one_letter_code for res in self.get_polymer_residues()])

    def get_standard_residues(self):
        return [r for r in self.residues if r.is_polymer_residue() and not r.is_modified_residue()]

    def get_heteroresidues(self, exclude_water=True):
        return [r for r in self.residues if r.is_heteroresidue(exclude_water=exclude_water)]

    def get_waters(self):
        return [r for r in self.residues if r.is_water()]

    def has_waters(self):
        return len(self.get_waters()) > 0

    def get_residues(self, standard=True, ligands=False, modified_residues=True, water=False):
        list_of_residues = []
        for res in self.residues:
            if res.is_standard_residue() and standard:
                list_of_residues.append(res)
            elif res.is_ligand() and ligands:
                list_of_residues.append(res)
            elif res.is_modified_residue() and modified_residues:
                list_of_residues.append(res)
            elif res.is_water() and water:
                list_of_residues.append(res)
        return list_of_residues


    def get_residue_by_index(self, index, aligned_sequence_index=False, only_polymer=True):
        """
        Returns a residue having the index provided in the 'index' argument in the sequence.
        """
        if aligned_sequence_index:
            index = seq_manipulation.get_residue_id_in_gapless_sequence(self.my_sequence, index)
            if index is None:
                raise KeyError(index)
        if only_polymer:
            return self.get_polymer_residues()[index]
        else:
            return self.residues[index]

    def get_residue_by_db_index(self, db_index):
        for res in self.residues:
            if res.db_index == db_index:
                return res
        raise KeyError(db_index)


    ###############################################################################################
    # Structure related.                                                                          #
    ###############################################################################################

    def set_structure(self, structure_object):
        self.structure = structure_object
        self.structure.set_pymod_element(self)


    def has_structure(self):
        return self.structure != None


    def check_structure(method):
        """
        Decorator method to check if PyMod_element object hasn an associated 3D structure.
        """
        def checker(self, *args, **config):
            if not self.has_structure():
                raise PyModMissingStructure("The element does not have a structure.")
            return method(self, *args, **config)
        return checker


    @check_structure
    def remove_structure(self):
        self.structure = None

    @check_structure
    def rename_chain_structure_files(self, chain_structure_file):
        self.set_initial_chain_file(chain_structure_file)
        self.set_current_chain_file(chain_structure_file)

    @check_structure
    def set_initial_chain_file(self, initial_chain_file):
        self.structure.initial_chain_file_path = initial_chain_file
        self.structure.initial_chain_file_name = os.path.basename(self.structure.initial_chain_file_path)

    @check_structure
    def set_current_chain_file(self, current_chain_file):
        self.structure.current_chain_file_path = current_chain_file
        self.structure.current_chain_file_name = os.path.basename(self.structure.current_chain_file_path)

    @check_structure
    def update_all_structure_paths(self, structure_dirpath):
        try:
            self.structure.initial_chain_file_path = os.path.join(structure_dirpath, self.structure.initial_chain_file_name)
            self.structure.current_chain_file_path = os.path.join(structure_dirpath, self.structure.current_chain_file_name)
            self.structure.current_full_file_path = os.path.join(structure_dirpath, self.structure.current_full_file_name)
        # Compatibility with old PyMod versions which did not have the "*_file_name" attributes.
        except AttributeError:
            pass


    @check_structure
    def rename_file_name_root(self, file_name_root):
        self.structure.file_name_root = file_name_root

    @check_structure
    def get_structure_file(self, basename_only=True, strip_extension=False, initial_chain_file=False, full_file=False):
        """
        Returns the name of the structure file of the element.
        'basename_only': if 'True' returns only the basename of the structure file.
        'strip_extension': if 'True' removes the extension to the from the returned file path.
        'initial_chain_file': if 'True' returns the structure file of the original chain (with
                                   its original atomic coordinates). If 'False' returns the current
                                   structure file of the chain (with modified atomic coordinates if
                                   the structure has been superposed).
        'full_file': if 'True' returns the file path of the original PDB file of the structure (the
                     PDB file containing the chain of the element along any other chain belonging to
                     other elements).
        """
        assert(not (initial_chain_file and full_file))
        if initial_chain_file:
            result = self.structure.initial_chain_file_path
        elif full_file:
            result = self.structure.current_full_file_path
        else:
            result = self.structure.current_chain_file_path
        if basename_only:
            result = os.path.basename(result)
        if strip_extension:
            result = os.path.splitext(result)[0]
        return result

    @check_structure
    def get_chain_id(self):
        return self.structure.chain_id

    @check_structure
    def get_chain_numeric_id(self):
        return self.structure.numeric_chain_id


    #################################################################
    # PyMOL related.                                                #
    #################################################################

    @check_structure
    def get_pymol_selector(self):
        return os.path.splitext(os.path.basename(self.structure.initial_chain_file_path))[0]


    #################################################################
    # Structural features.                                          #
    #################################################################

    @check_structure
    def get_disulfides(self):
        return self.structure.disulfides_list

    @check_structure
    def has_disulfides(self):
        return self.structure.disulfides_list != []

    @check_structure
    def add_disulfide(self, disulfide=None):
        self.structure.disulfides_list.append(disulfide)


    ###############################################################################################
    # Annotation related.                                                                         #
    ###############################################################################################

    def has_assigned_secondary_structure(self):
        return bool(self.assigned_secondary_structure)

    def has_predicted_secondary_structure(self):
        return bool(self.predicted_secondary_structure)

    def has_campo_scores(self):
        return self._has_campo_scores

    def has_entropy_scores(self):
        return self._has_entropy_scores

    def has_pc_scores(self):
        return self._has_pc_scores

    def has_dope_scores(self):
        return self._has_dope_scores

    def has_domains(self, only_original=False):
        if only_original:
            return self.features["domains"] and self.splitted_domain == None
        else:
            return self.features["domains"]

    def has_features(self):
        return self.features["sequence"]


    def pdb_is_fetchable(self):
        return self.pdb_can_be_fetched

    def can_be_modeled(self):
        r = False
        # Exlude cluster elements (alignments and BLAST-searches).
        if not self.is_cluster():
            r = True

        # if self.has_structure():
        #     r = False
        # # The element is a primary sequence imported by the user.
        # else:
        #     r = True

        return r

    def is_suitable_template(self):
        return self.has_structure() # TODO: and sequence.is_model != True:


    def is_model(self):
        return isinstance(self, PyMod_model_element)


class PyMod_root_element(_PyMod_cluster_element):
    pass


class _PyMod_model_element(_PyMod_sequence_element):

    def __init__(self, model_root_name, *args, **configs):
        self.model_root_name = model_root_name
        _PyMod_sequence_element.__init__(self, *args, **configs)


###################################################################################################
# CLASSES FOR EXTENDING PYMOD ELEMENTS TO CONTROL DATA OF THE PLUGIN.                             #
###################################################################################################

class Added_PyMod_element:
    """
    A mixin class for PyMod elements storing information about the whole PyMod plugin. It extends
    their methods so that they will also control other elements of PyMod.
    """

    def initialize(self, pymod):
        self.pymod = pymod

    def set_as_lead(self):
        """
        Also remove other lead statuses from siblings, since there can be only one lead per cluster.
        """
        for sibling in self.get_siblings():
            sibling.remove_all_lead_statuses()
        super(Added_PyMod_element, self).set_as_lead()

    def set_as_query(self):
        for sibling in self.get_siblings():
            sibling.remove_all_lead_statuses()
        super(Added_PyMod_element, self).set_as_query()

    def extract_to_upper_level(self, place_below_mother=True):
        """
        Extract elements from their clusters. Also changes the position of the item in the list of
        PyMod elements.
        """
        old_mother_index = self.pymod.get_pymod_element_index_in_container(self.mother) + 1
        super(Added_PyMod_element, self).extract_to_upper_level()
        if place_below_mother:
            self.pymod.change_pymod_element_list_index(self, old_mother_index)
        # Allows to show the widgets if the child element was in collapsed cluster.
        self.widget_group.show = True
        # Hides the child signes.
        self.widget_group.cluster_button.setVisible(False)

    def delete(self):
        """
        Used to remove definitively an element from PyMod.
        """
        # Delete its structure in PyMOL.
        if self.has_structure():
            self.pymod.delete_pdb_file_in_pymol(self)
        # Delete its children.
        if self.is_mother():
            children = self.get_children()
            for c in children[:]:
                c.delete()
        # Actually delete the element.
        super(Added_PyMod_element, self).delete()


###################################################################################################
# CLASSES FOR EXTENDING PYMOD ELEMENTS BEHAVIOUR TO CONTROL PYMOD GUI.                            #
###################################################################################################

class PyMod_element_GUI(Added_PyMod_element):
    """
    A mixin class extending the behaviour of 'PyMod_element' objects, so that their methods will
    also control the elements' widgets in PyMod main window.
    """

    def initialize(self, *args, **configs):
        Added_PyMod_element.initialize(self, *args, **configs)
        # Builds for it some GUI widgets to show in PyMod window. They will be gridded later.
        self.pymod.main_window.add_pymod_element_widgets(self)

    def remove_all_lead_statuses(self):
        # self.show_collapsed_mother_widgets()
        super(PyMod_element_GUI, self).remove_all_lead_statuses()
        # Removes the cluster button of leads of collapsed clusters.
        self.widget_group.cluster_button.set_visibility()

    def set_as_lead(self):
        super(PyMod_element_GUI, self).set_as_lead()
        self.widget_group.cluster_button.set_visibility()


    def extract_to_upper_level(self, *args, **configs):
        """
        Extract elements from their clusters. Also modifies its widgets and those of their parents.
        """
        # self.show_collapsed_mother_widgets()
        super(PyMod_element_GUI, self).extract_to_upper_level(*args, **configs)

    def delete(self, *args, **configs):
        """
        Delete an element.
        """
        # self.show_collapsed_mother_widgets()
        super(PyMod_element_GUI, self).delete(*args, **configs)
        # Remove the element widgets.
        self.pymod.main_window.delete_element_widgets(self)

    def show_collapsed_mother_widgets(self):
        """
        If the element is a lead of a collapsed cluster, then show the widgets of the cluster (which
        were hidden) after the element lead element is extracted.
        """
        if self.is_lead() and self.pymod.main_window.is_collapsed_cluster(self.mother):
            self.mother.widget_group.show = True

    def add_child(self, child, *args, **configs):
        super(PyMod_element_GUI, self).add_child(child, *args, **configs)
        # If the cluster is not collapsed, show the widgets of the children.
        if not self.pymod.main_window.is_collapsed_cluster(self):
            child.widget_group.show = True
        # If the cluster is collapsed hide it.
        else:
            child.widget_group.hide_widgets()


    def __getstate__(self):
        """
        Called by pickle when dumping an object to a file.
        """
        return pymod_getstate(self)


###################################################################################################
# CLASSES ACTUALLY USED IN THE PLUGIN.                                                            #
###################################################################################################

class PyMod_cluster_element(PyMod_element_GUI, _PyMod_cluster_element):
    pass

class PyMod_sequence_element(PyMod_element_GUI, _PyMod_sequence_element):
    pass

class PyMod_model_element(PyMod_element_GUI, _PyMod_model_element):
    pass


###################################################################################################
# Pickling behaviour.                                                                             #
###################################################################################################

# Attributes which can not be pickled because they contain PyQt or Tkinter objects. Used when saving
# Python objects in order to save a PyMod session.
pymod_element_unpickable_elements = ["widget_group", "pymod"]
pymod_unpickable_elements = []
unpickable_attributes = pymod_element_unpickable_elements + pymod_unpickable_elements

def pymod_getstate(obj):
    d = dict(obj.__dict__)
    for k in unpickable_attributes:
        if k in d:
            del d[k]
    return d
