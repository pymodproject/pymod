# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for performing homology modeling with MODELLER in PyMod.
"""

import os
import sys
import shutil
import importlib
import pickle

from pymol import cmd

try:
    import modeller
    import modeller.automodel
    from modeller.scripts import complete_pdb
    from modeller import soap_pp
    from modeller import parallel
except:
    pass

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol, MODELLER_common
from pymod_lib.pymod_protocols.structural_analysis_protocols.dope_assessment import DOPE_assessment, show_dope_plot, compute_dope_of_structure_file
from pymod_lib.pymod_element_feature import Element_feature
from pymod_lib.pymod_threading import Protocol_exec_dialog

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_seq.seq_manipulation as pmsm
import pymod_lib.pymod_structure as pmstr
from pymod_lib.pymod_protocols.modeling_protocols.homology_modeling._gui_qt import Modeling_window_qt


###################################################################################################
# HOMOLOGY MODELING.                                                                              #
###################################################################################################

class MODELLER_homology_modeling(PyMod_protocol, MODELLER_common):
    """
    Class to represent an homology model building session with MODELLER.
    """

    # The maximum number of models that Modeler can produce at the same time.
    max_models_per_session = 1000
    default_modeller_random_seed = -8123
    min_modeller_random_seed = -50000
    max_modeller_random_seed = -2

    # If set to True, creates a file "my_model.py" that can be used by command line MODELLER to
    # perform the modellization. It is necessary when using MODELLER as an external coomand line
    # tool, that is, when using MODELLER on PyMOL version which can't import the systemwide
    # 'modeller' library. Users may also take the script and modify it according to their needs.
    write_modeller_script_option = True

    modeling_dirpath = ""
    modeling_files_name = "my_model"
    modeling_script_name = "%s.py" % modeling_files_name
    modeling_script_lib_name = "%s_lib.py" % modeling_files_name
    modeling_lib_name = modeling_script_lib_name.split(".")[0]
    modeling_log_name = "%s.log" % modeling_files_name
    pir_file_name = "align-multiple.ali"
    modeller_temp_output_name = "modeller_saved_outputs.txt"
    modeller_vars_filename = "modeller_vars.pkl"

    # single_chain_models_formatted_string = "m%s_%s" # (decoy_number, model_name)
    single_chain_models_name = "model"
    single_chain_loop_models_name = "loop"
    multiple_chains_models_name = "multi_model" # "MyMultiModel"
    multiple_chains_loop_models_name = "multi_loop"

    tc_temp_pymol_name = "template_complex_temp"
    mc_temp_pymol_name = "model_complex_temp"

    single_chain_model_color = "white"
    list_of_model_chains_colors = pmdt.pymol_light_colors_list
    loop_default_color = "orange" # "red"

    protocol_name = "modeller_homology_modeling"


    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)
        MODELLER_common.__init__(self)
        self.use_hetatm_in_session = False
        self.use_water_in_session = False
        self.models_count = 0

        self.modeling_protocol_initialization()

    def modeling_protocol_initialization(self):
        self.modeling_mode = "homology_modeling"
        self.modeling_window_class_qt = Modeling_window_qt


    def round_assessment_value(self, assessment_value, digits_after_point=3):
        return round(assessment_value, digits_after_point)


    def launch_from_gui(self):
        """
        This method is called when the "MODELLER" option is clicked in the "Tools" menu.
        """

        # Try to find if Modeller is installed on the user's computer.
        modeller_error = self.pymod.modeller.check_exception()
        if modeller_error is not None:
            message = "In order to perform 3D modeling in PyMod, MODELLER must be installed and configured correctly. %s" % modeller_error
            self.pymod.main_window.show_error_message("MODELLER Error", message)
            return None

        self.initialize_modeling_protocol_from_gui()


    ###########################################################################
    # Checks the input selection.                                             #
    ###########################################################################

    def initialize_modeling_protocol_from_gui(self):

        #---------------------------------------------------------------------
        # Get the selected sequences to see if there is a correct selection. -
        #---------------------------------------------------------------------

        selected_sequences = self.pymod.get_selected_sequences()

        # First check if at least one sequence is selected.
        if not len(selected_sequences) > 0:
            self.pymod.main_window.show_error_message("Selection Error", "Please select at least one target sequence to use MODELLER.")
            return None

        # Checks if all the selected sequences can be used to build a model.
        if False in [s.can_be_modeled() for s in selected_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only sequences that do not have a structure loaded in PyMOL.")
            return None

        # Checks that all the selected sequences are currently aligned to some other sequence
        # (aligned sequences are always 'children'). Only sequences aligned to some template can be
        # modeled.
        if False in [e.is_child() for e in selected_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Please select only target sequences that are currently aligned to some structure.")
            return None

        #----------------------------------------------------------------------------------------
        # Builds the modeling clusters which will store the information needed to run MODELLER. -
        #----------------------------------------------------------------------------------------

        # Using the 'build_cluster_list()' method build the 'self.involved_cluster_elements_list'
        # just like when performing an alignment.
        self.build_cluster_lists()

        # This will contain a list of 'Modeling_cluster' objects.
        self.modeling_clusters_list = []

        # Checks other conditions in each cluster (that is, checks if there is only one target
        # sequence selected per cluster and if it is aligned to some suitable template).
        for cluster_element in self.involved_clusters_list:

            # Checks that only one sequence per cluster is selected.
            if not self.pymod.check_only_one_selected_child_per_cluster(cluster_element):
                self.pymod.main_window.show_error_message("Selection Error", "Please select only one target sequence in the following cluster: %s" % (cluster_element.my_header))
                return False

            # Look if there is at least one suitable template aligned to the target sequence.
            templates_temp_list = []
            target_name = None
            found_nucleic_acid_template = False
            for e in cluster_element.get_children():
                 if not e.selected and e.is_suitable_template():
                     if e.polymer_type == "nucleic_acid":
                         found_nucleic_acid_template = True
                     templates_temp_list.append(e)
                 if e.selected:
                     target_name = e.my_header

            # Checks if some templates have been found.
            if not len(templates_temp_list) > 0:
                self.pymod.main_window.show_error_message("Selection Error", "The target sequence '%s' is currently not aligned to any suitable template." % (target_name))
                return False

            # Checks for the absence of nucleic acid templates.
            if found_nucleic_acid_template:
                self.pymod.main_window.show_error_message("Selection Error", "The target sequence '%s' is currently aligned to a nucleic acid structure, which can not be used as a template." % (target_name))
                return False

            #-------------------------------------------------
            # Actually builds the modeling clusters objects. -
            #-------------------------------------------------

            self.modeling_clusters_list.append(Modeling_cluster(self, cluster_element))

        #--------------------------------------
        # Build the homology modeling window. -
        #--------------------------------------

        # Define the modeling mode.
        self.multiple_chain_mode = len(self.modeling_clusters_list) > 1

        # Single chain homology modeling mode.
        if not self.multiple_chain_mode:
            self.build_modeling_window()

        # Multiple chains modeling requires the identification of "template complexes" and
        # additional controls.
        else:
            # This will build the 'self.available_template_complex_list'.
            self.initialize_multichain_modeling()
            # Proceeds only if there is at least one suitable "template complex".
            if len(self.available_template_complex_list) > 0:
                self.build_modeling_window()
            else:
                self.pymod.main_window.show_error_message("Selection Error", "There aren't any suitable 'Template Complexes' to perform multiple chain homology modeling.")


    def initialize_multichain_modeling(self):
        """
        This method will prepare data needed to perform multichain modeling. It will:
            - identify suitable template complexes
            - check if there are target sequences with the same sequence, so that symmetry restraints
              can be applied to them when using Modeller.
        """

        if self.modeling_mode == "homology_modeling":

            #--------------------------------------------
            # Generates modeling clusters dictionaries. -
            #--------------------------------------------

            for mc in self.modeling_clusters_list:
                mc.build_template_complex_chains_dict()

            #--------------------------------------------------
            # Builds a list of suitable "template complexes". -
            #--------------------------------------------------

            # A "teplate complex" is available only if in each selected cluster there is at least ONE
            # chain coming from the same original PDB file. For example, with these two cluster:
            #     - cluster 1: <1HHO_Chain:A>, 2DN2_Chain:A
            #     - cluster 2: <1HHO_Chain:B>, 3EOK_Chain:A
            # the "template complex" is 1HHO.

            # self.available_template_complex_list = []
            # self.available_template_complex_dict = {}

            codes_per_cluster = [set(mc.full_structure_files_dict.keys()) for mc in self.modeling_clusters_list]
            self.available_template_complex_list = list(set.intersection(*codes_per_cluster))
            self.available_template_complex_list.sort() # Sorts the list alphabetically.
            self.all_full_structure_fn_dict = {}
            for mc in self.modeling_clusters_list:
                self.all_full_structure_fn_dict.update(mc.full_structure_fn_dict)

        #--------------------------------------------------------------------------------------
        # Checks if there are some target chains with the same sequence, so that the user may -
        # apply symmetry restraints to them.                                                  -
        #--------------------------------------------------------------------------------------

        # Builds a 'Symmetry_restraints_groups' object that is going to be used to keep track of
        # modeling clusters that have a target sequence with the same sequence.
        self.symmetry_restraints_groups = Symmetry_restraints_groups_list()
        for mc in self.modeling_clusters_list:
            seq = str(mc.target.my_sequence).replace("-","")
            # Adds a "symmetry restraints group" for each group of target sequences that share the
            # exact same sequence.
            if not seq in [g.id for g in self.symmetry_restraints_groups.get_groups()]:
                self.symmetry_restraints_groups.add_group(seq)
            self.symmetry_restraints_groups.get_group_by_id(seq).add_cluster(mc)
        # Also assigns "symmetry ids" to each modeling cluster.
        for mc in self.modeling_clusters_list:
            seq = str(mc.target.my_sequence).replace("-","")
            if seq in [g.id for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2)]:
                mc.symmetry_restraints_id = seq
            else:
                mc.symmetry_restraints_id = None


    def build_modeling_window(self):
        """
        Builds the modeling window with all the options for MODELLER.
        """
        # self.modeling_window = self.modeling_window_class(self.pymod.get_pymod_app(), self)

        self.modeling_window = self.modeling_window_class_qt(parent=self.pymod.main_window,
                                                             protocol=self)
        self.modeling_window.show()
        # self.modeling_window.resize(700, 800)


    def launch_modelization(self):
        """
        This method is called when the 'SUBMIT' button in the modelization window is pressed. It
        contains the code to instruct Modeller on how to perform the modelization.
        """

        if not self.prepare_modelization():
            return None
        self.perform_modelization()


    ###########################################################################
    # Checks the modeling window parameters.                                  #
    ###########################################################################

    def get_modeling_options_from_gui(self):
        """
        Add to the 'Modeling_clusters' objects information about which templates to use according to
        the parameters supplied by users.
        """

        #--------------------------------------
        # Get options from the 'Options' tab. -
        #--------------------------------------

        self.starting_model_number = 1

        # Checks if a correct value in the max models entry has been supplied.
        self.ending_model_number = self.modeling_window.max_models_enf.getvalue(validate=True)

        # Get heteroatoms options.
        if self.modeling_mode == "homology_modeling":
            self.exclude_hetatms = pmdt.yesno_dict[self.modeling_window.exclude_heteroatoms_rds.getvalue()]
        elif self.modeling_mode == "loop_refinement":
            self.include_hetatms = self.modeling_window.exclude_heteroatoms_rds.getvalue()

        # Check if a correct random seed value has been supplied.
        self.modeller_random_seed = self.modeling_window.random_seed_enf.getvalue(validate=True)

        # Check that a correct number of parallel jobs has been supplied.
        self.modeller_n_jobs = self.modeling_window.n_jobs_enf.getvalue(validate=True)
        if self.modeller_n_jobs != 1 and self.modeller_n_jobs > self.ending_model_number:
            message = ("You can not use more MODELLER parallel jobs (%s) than the"
                       " number of models you want to build (%s). Please decrease"
                       " the number of parallel jobs." % (self.modeller_n_jobs, self.ending_model_number))
            raise ValueError(message)

        # Checks the optimization level parameters.
        self.optimization_level = self.modeling_window.optimization_level_rds.getvalue()
        if self.optimization_level == "Custom":
            self.check_custom_optimization_level_parameters()

        # Checks if the custom objective function parameters are correct.
        self.use_custom_obj_func = self.modeling_window.get_custom_obj_func_option()
        if self.use_custom_obj_func == "customize":
            self.check_custom_obj_func_parameters()
        elif self.use_custom_obj_func == "altmod":
            self.set_altmod_obj_func_parameters()
        else:
            if self.modeling_mode == "loop_refinement":
                # By default use the 'loopmodel' class.
                self.loop_refinement_class_name = "loopmodel"

        if self.modeling_mode == "homology_modeling":
            self.superpose_to_templates = pmdt.yesno_dict[self.modeling_window.superpose_models_to_templates_rds.getvalue()]
        if self.modeling_mode == "loop_refinement":
            self.loop_initial_conformation = self.modeling_window.loop_initial_conformation_rds.getvalue()

        self.color_models_by_choice = self.modeling_window.color_models_rds.getvalue()

        #------------------------------------------
        # Gets options for each modeling cluster. -
        #------------------------------------------

        for modeling_cluster in self.modeling_clusters_list:
            # Begins a for cycle that is going to get the structures to be used as templates.
            # When using loop refinement, it will get the loop intervals.
            modeling_cluster.set_options_from_gui()

        #------------------------------------
        # Set the template complex options. -
        #------------------------------------

        if self.multiple_chain_mode:
            self._get_multichain_options_from_gui()

        #------------------------------------
        # Check if hetatms have to be used. -
        #------------------------------------

        # Homology modeling.
        if self.modeling_mode == "homology_modeling":
            if self.exclude_hetatms:
                # The 'hetatm' flag of MODELLER will be set to 'True', but each template heteroatom
                # will be aligned to a gap character in the target sequence (that is, the heteroatom
                # will not be included in the model).
                self.use_hetatm_in_session = False
                self.use_water_in_session = False
            else:
                self.use_hetatm_in_session = True
                # Check if water molecules have to be included in the modeling session.
                self.use_water_in_session = True in [mc.use_water_in_cluster() for mc in self.modeling_clusters_list]

        # Loop refinement.
        else:
            if self.include_hetatms == "Het + Water":
                self.use_hetatm_in_session = True
                self.use_water_in_session = True
            elif self.include_hetatms == "Het":
                self.use_hetatm_in_session = True
                self.use_water_in_session = False
            elif self.include_hetatms == "No":
                self.use_hetatm_in_session = False
                self.use_water_in_session = False


        #------------------------------------------------------------------------------------------
        # Builds a list with the "knowns" for MODELLER and sets the name of the target sequences. -
        #------------------------------------------------------------------------------------------

        # Homology modeling.
        if self.modeling_mode == "homology_modeling":

            self.all_templates_namelist = []
            self.modeller_target_name = ""
            self.pir_sequences_dict = {}

            if not self.multiple_chain_mode:
                self.all_templates_namelist = self.modeling_clusters_list[0].get_template_nameslist()
                self.modeller_target_name = self.single_chain_models_name # self.modeling_clusters_list[0].target_name
            elif self.multiple_chain_mode:
                for mc in self.modeling_clusters_list:
                    for t in mc.templates_list:
                        # Includes the "template complex" name only once.
                        if mc.is_template_complex_chain(t):
                            if not self.template_complex_modeller_name in self.all_templates_namelist:
                                self.all_templates_namelist.append(self.template_complex_modeller_name)
                        else:
                            self.all_templates_namelist.append(mc.template_options_dict[t]["modeller_name"])
                self.modeller_target_name = self.multiple_chains_models_name

        # Loop modeling.
        else:
            if not self.multiple_chain_mode:
                self.modeller_target_name = self.single_chain_loop_models_name # self.modeling_clusters_list[0].target_name
            elif self.multiple_chain_mode:
                self.modeller_target_name = self.multiple_chains_loop_models_name


        #-------------------------------------
        # Get options for disulfide bridges. -
        #-------------------------------------

        if self.modeling_mode == "homology_modeling":
            self.use_template_dsb = self.modeling_window.get_use_template_dsb_var()
            self.use_auto_dsb = self.modeling_window.get_auto_dsb_var()
        elif self.modeling_mode == "loop_refinement":
            self.use_template_dsb = 0
            self.use_auto_dsb = 0
        self.use_user_defined_dsb = self.modeling_window.get_use_user_defined_dsb_var()
        if self.use_user_defined_dsb:
            self.all_user_defined_dsb = self.modeling_window.get_user_dsb_list()

        #---------------------------------------
        # Get options for symmetry restraints. -
        #---------------------------------------

        if self.multiple_chain_mode:
            self.symmetry_restraints_groups.get_symmetry_restraints_from_gui()


    def _get_multichain_options_from_gui(self):
        """
        Finds the PDB_file object of the "template complex" selected by the user.
        """
        # Name of the original full structure file of the template complex.
        self.template_complex_filename = self.available_template_complex_list[self.modeling_window.get_template_complex_var()]
        # Name of the template complex.
        self.template_complex_name = self.all_full_structure_fn_dict[self.template_complex_filename]
        self.template_complex_modeller_name = self.template_complex_name[:-4]
        for mc in self.modeling_clusters_list:
            for t in mc.templates_list:
                if self.chain_is_from_template_complex(t):
                    mc.set_template_complex_chain(t)


    def set_modeling_options(self):
        """
        Set additional modeling options after some initial paramaters from the GUI have been
        checked.
        """
        if self.multiple_chain_mode:
            self.list_of_symmetry_restraints = self.symmetry_restraints_groups.get_symmetry_restraints_list()


    def chain_is_from_template_complex(self, pymod_element):
        return pymod_element.get_structure_file(full_file=True) == self.template_complex_filename


    def check_all_modeling_parameters(self):
        """
        This will be used before launching Modeller to check:
            - if the parameters of each modeling clusters are correct
            - when performing multichain modeling
                - if there is exactly 1 template complex chain selected in each cluster
                - if symmetry restraints buttons are selected properly
        if some parameters are not correct, an exception will be raised.
        """

        # Checks if the parameters of all the "modeling clusters" are correct.
        for mc in self.modeling_clusters_list:
            self.check_modeling_cluster_parameters(mc)

        # Check if there are only correct sequences.
        for mc in self.modeling_clusters_list:
            if not pmsm.check_correct_sequence(mc.target.my_sequence):
                raise ValueError("The target sequence '%s' contains an invalid character in its sequence (%s) and MODELLER can't modelize it." % (mc.target.my_header, pmsm.get_invalid_characters_list(mc.target.my_sequence)[0]))
            if self.modeling_mode == "homology_modeling":
                for t in mc.templates_list:
                    if not pmsm.check_correct_sequence(t.my_sequence):
                        raise ValueError("The template '%s' contains an invalid character in its sequence (%s) and MODELLER can't use it as a template." % (t.my_header, pmsm.get_invalid_characters_list(t.my_sequence)[0]))

        # If each "modeling cluster" has correct parameters, when performing multiple chain modeling,
        # there are other conditions that must be satisfied.
        if self.multiple_chain_mode:

            if self.modeling_mode == "homology_modeling":
                # Then perform additional controls for each modeling cluster and also get the list of
                # the "target complex" chains selected by the user.
                for mc in self.modeling_clusters_list:
                    # Gets the "template complex" chains selected in the current modeling cluster.
                    template_complex_selected_chains_in_cluster = [t for t in mc.templates_list if mc.is_template_complex_chain(t)]
                    # Check if the current cluster has a selected chain from the "target complex".
                    if len(template_complex_selected_chains_in_cluster) == 0:
                        raise ValueError("Please select AT LEAST one chain from the 'Template Complex' (%s) as a template for %s!" % (self.template_complex_name, mc.target_name))
                    # Checks if in some cluster there is more than one selected template belonging to the
                    # "template complex". This is needed for because ONLY one chain belonging to the
                    # "template complex" can be selected by ther user in each cluster.
                    if len(template_complex_selected_chains_in_cluster) > 1:
                        raise ValueError("Please select ONLY one chain from the 'Template Complex' (%s) as template for %s!" % (self.template_complex_name, mc.target_name))
                    # Sets the template complex in the modeling cluster. This will first set the
                    # 'block_index' attribute.
                    mc.set_template_complex_chain_to_use(template_complex_selected_chains_in_cluster[0])

                # Fixes the 'block_index' attribute so that each block index will correspond to the
                # numeric ID of the chain in the 3D model.
                for mc_idx, mc in enumerate(self.get_modeling_clusters_list(sorted_by_id=True)):
                    mc.block_index = mc_idx

            # Finally checks if the symmetries checkbuttons are selected properly.
            self.check_symmetry_restraints_vars()

        # Check if the alignments can be given as correct input to MODELLER.
        if self.modeling_mode == "homology_modeling":
            self.check_alignments()

        # Warns the user if parallel MODELLER jobs can not be used.
        if self.modeller_n_jobs > 1 and not self.pymod.parallel_modeller:
            self.modeller_n_jobs = 1
            fix_pythonpath_modeller_ref = "the PyMod user's guide"
            message = ("Can not use multiple MODELLER jobs. MODELLER will be launched"
                       " serially, using only one job. If you want to use MODELLER in"
                       " a parallel way, re-install MODELLER in a way that it can be"
                       " fully compatible with PyMod/PyMOL (see %s for more"
                       " information)." % fix_pythonpath_modeller_ref)
            self.pymod.main_window.show_warning_message("MODELLER warning", message)


    def check_custom_optimization_level_parameters(self):
        """
        Checks if the custom optimization level parameters are correct.
        """

        # Max CG iterations.
        if self.modeling_window.optimization_level_frame.use_max_cg_iterations_enf:
            self.custom_max_cg_iterations = self.modeling_window.optimization_level_frame.max_cg_iterations_enf.getvalue(validate=True)

        # VTFM schedule.
        if self.modeling_window.optimization_level_frame.use_vtfm_schedule_rds:
            self.custom_vtfm_schedule = self.modeling_window.optimization_level_frame.vtfm_schedule_rds.getvalue()

        # MDSA protocol.
        if self.modeling_window.optimization_level_frame.use_md_schedule_rds:
            self.custom_md_schedule = self.modeling_window.optimization_level_frame.md_schedule_rds.getvalue()

        # Repeat optimization.
        if self.modeling_window.optimization_level_frame.use_repeat_optimization_enf:
            self.custom_repeat_optimization = self.modeling_window.optimization_level_frame.repeat_optimization_enf.getvalue(validate=True)

        # Max objective function value.
        if self.modeling_window.optimization_level_frame.use_max_obj_func_value_enf:
            self.custom_max_obj_func_value = self.modeling_window.optimization_level_frame.max_obj_func_value_enf.getvalue(validate=True)

        return True


    def check_custom_obj_func_parameters(self):
        """
        Check the paramaters for customizing the objective function.
        """

        # Automodel class for loop refinement.
        if self.modeling_window.custom_obj_func_frame.use_loop_stat_pot_rds:
            self.loop_refinement_class_name = self.modeling_window.custom_obj_func_frame.loop_stat_pot_rds.getvalue()

        # Homology-derived distance restraints terms.
        if self.modeling_window.custom_obj_func_frame.use_hddr_frame:
            self.custom_w_hddr = self.modeling_window.custom_obj_func_frame.hddr_frame.getvalue(validate=True)

        # Dihedral angle restraints terms.
        if self.modeling_window.custom_obj_func_frame.use_hdar_frame:
            self.custom_w_hdar = self.modeling_window.custom_obj_func_frame.hdar_frame.getvalue(validate=True)

        # CHARM22 terms.
        if self.modeling_window.custom_obj_func_frame.use_charmm22_frame:
            self.custom_w_charmm22 = self.modeling_window.custom_obj_func_frame.charmm22_frame.getvalue(validate=True)

        # Soft sphere overlap terms.
        if self.modeling_window.custom_obj_func_frame.use_soft_sphere_frame:

            # DOPE loopmodel does not use soft sphere terms.
            get_custom_w_soft_sphere = True if self.modeling_mode == "homology_modeling" else self.loop_refinement_class_name != "DOPE loopmodel"

            if get_custom_w_soft_sphere:
                self.custom_w_soft_sphere = self.modeling_window.custom_obj_func_frame.soft_sphere_frame.getvalue(validate=True)
            else:
                self.custom_w_soft_sphere = 1.0

        # Use DOPE in homology modeling.
        if self.modeling_window.custom_obj_func_frame.use_dope_frame:

            self.use_dope_in_obj_func = self.modeling_window.custom_obj_func_frame.get_use_dope_var()

            if self.use_dope_in_obj_func:
                self.custom_w_dope = self.modeling_window.custom_obj_func_frame.dope_frame.getvalue(validate=True)
            else:
                self.custom_w_dope = 0.0

        # Statistical potentials terms in loop modeling.
        if self.modeling_window.custom_obj_func_frame.use_stat_pot_frame:
            self.custom_w_stat_pot_loop = self.modeling_window.custom_obj_func_frame.stat_pot_frame.getvalue(validate=True)

        # Lennard-Jones terms.
        if self.modeling_window.custom_obj_func_frame.use_lj_frame:

            get_custom_w_lj = False if self.modeling_mode == "homology_modeling" else self.loop_refinement_class_name == "DOPE loopmodel"

            if get_custom_w_lj:
                self.custom_w_lj = self.modeling_window.custom_obj_func_frame.lj_frame.getvalue(validate=True)
            else:
                self.custom_w_lj = 1.0

        # GBSA terms.
        if self.modeling_window.custom_obj_func_frame.use_gbsa_frame:

            get_custom_w_gbsa = False if self.modeling_mode == "homology_modeling" else self.loop_refinement_class_name == "DOPE loopmodel"

            if get_custom_w_gbsa:
                self.custom_w_gbsa = self.modeling_window.custom_obj_func_frame.gbsa_frame.getvalue(validate=True)
            else:
                self.custom_w_gbsa = 1.0

        # Non-bonded interactions cutoff.
        self.custom_nb_cutoff = self.modeling_window.custom_obj_func_frame.nb_cutoff_frame.getvalue(validate=True)

        return True

    def set_altmod_obj_func_parameters(self):
        self.custom_w_hddr = 1.0
        self.custom_w_hdar = 1.0
        self.custom_w_charmm22 = 1.0
        self.custom_w_soft_sphere = 1.0
        self.use_dope_in_obj_func = True
        self.custom_w_dope = 0.5
        self.custom_nb_cutoff = 8.0


    def check_modeling_cluster_parameters(self, modeling_cluster):
        """
        Checks the if there are any problems with the user-supplied parameters of a "modeling cluster"
        before starting the modeling process.
        """
        # Checks if there are some templates that have been selected.
        if modeling_cluster.templates_list == []:
            raise ValueError("You have to select at least one template for target '%s' in order to build a model!" % (modeling_cluster.target_name))


    def check_symmetry_restraints_vars(self):
        """
        Check if symmetry restraints for multiple chain modeling can be applied.
        """
        for srg in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2):
            if srg.use == 1:
                raise ValueError("In order to impose symmetry restraints you need select the 'Apply symmetry restraints' option for at least two targets with the same sequence (you selected this option only for target '%s')." % (srg.list_of_clusters[0].target_name))


    def check_targets_with_cys(self):
        """
        Check if there is at least one modeling cluster with a target sequence with at least two CYS
        residues.
        """
        return True in [mc.has_target_with_multiple_cys() for mc in self.modeling_clusters_list]

    def check_structures_with_disulfides(self):
        """
        Checks in all modeling clusters if there is at least one structure with a disulfide bridge.
        """
        return True in [mc.has_structures_with_disulfides() for mc in self.modeling_clusters_list]


    def get_template_complex_chains(self, sorted_by_id=False):
        return [mc.get_template_complex_chain() for mc in self.get_modeling_clusters_list(sorted_by_id=sorted_by_id)]


    def get_targets_list(self, sorted_by_id=False):
        return [mc.target for mc in self.get_modeling_clusters_list(sorted_by_id=sorted_by_id)]


    def get_modeling_clusters_list(self, sorted_by_id=False):
        if sorted_by_id:
            return sorted(self.modeling_clusters_list, key = lambda mc: mc.block_index)
        else:
            return self.modeling_clusters_list


    def check_alignments(self):
        """
        Checks if correct alignments are being provided.
        """
        for modeling_cluster in self.modeling_clusters_list:
            for char_id, target_char in enumerate(modeling_cluster.target.my_sequence):
                if target_char == "X":
                    template_gaps = [template.my_sequence[char_id] == "-" for template in modeling_cluster.templates_list]
                    # The alignment is not correct if there are any 'X' character of the target not
                    # aligned to at least one residue of the templates.
                    if not False in template_gaps:
                        message = ("An X residue of the '%s' target sequence is currently"
                                   " not aligned to a template residue. Make sure that"
                                   " each X residue in your target sequences is aligned"
                                   " with some residue of your templates." % (modeling_cluster.target_name))
                        raise ValueError(message)


    def get_modeller_schedule_obj(self, phase, modeller_exec, name):

        # MODELLER 3D Model building parameters.
        modeller_3d_building_params_dict = {"cg": {"internal": {"fastest": modeller.automodel.autosched.fastest,
                                                                "very fast": modeller.automodel.autosched.very_fast,
                                                                "fast": modeller.automodel.autosched.fast,
                                                                "normal": modeller.automodel.autosched.normal,
                                                                "slow": modeller.automodel.autosched.slow},
                                                   "external": {"fastest": "modeller.automodel.autosched.fastest",
                                                                "very fast": "modeller.automodel.autosched.very_fast",
                                                                "fast": "modeller.automodel.autosched.fast",
                                                                "normal": "modeller.automodel.autosched.normal",
                                                                "slow": "modeller.automodel.autosched.slow"}},
                                            "md": {"internal": {"none": None,
                                                                "very fast": modeller.automodel.refine.very_fast,
                                                                "fast": modeller.automodel.refine.fast,
                                                                "slow": modeller.automodel.refine.slow,
                                                                "very slow": modeller.automodel.refine.very_slow},
                                                   "external": {"none": "None",
                                                                "very fast": "modeller.automodel.refine.very_fast",
                                                                "fast": "modeller.automodel.refine.fast",
                                                                "slow": "modeller.automodel.refine.slow",
                                                                "very slow": "modeller.automodel.refine.very_slow"}}}
        return modeller_3d_building_params_dict[phase][modeller_exec][name]


    ###########################################################################
    # Prepares input files for MODELLER.                                      #
    ###########################################################################

    def prepare_modeling_session_files(self):
        """
        Prepares the directory where MODELLER's output will be generated and moves into it.
        """
        #--------------------------------------------------------------------------
        # Build a directory where all the modeling session files will be located. -
        #--------------------------------------------------------------------------

        # The absolute path of the models directory.
        models_dir = os.path.join(self.pymod.current_project_dirpath, self.pymod.models_dirname)
        # Name of the model subdirectory where Modeller output files are going to be placed.
        model_subdir_name = "%s_%s_%s" % (self.pymod.models_subdirectory, self.pymod.performed_modeling_count + 1, self.modeller_target_name)
        # The absolute path of the model subdirectory.
        modeller_output_dir_path = os.path.join(models_dir, model_subdir_name)

        # Stores the path of the modeling directory.
        self.modeling_dirpath = modeller_output_dir_path
        if not os.path.isdir(self.modeling_dirpath):
            os.mkdir(self.modeling_dirpath)

        self.prepare_modeling_protocol_session_files()

        #-------------------------------------------------------------------
        # Changes the current working directory to the modeling directory. -
        #-------------------------------------------------------------------
        # The current directory has to be changed beacause in Modeller the user can't change the
        # output directory, it has to be the current directory.
        os.chdir(self.modeling_dirpath)


    def prepare_modeling_protocol_session_files(self):

        #-------------------------------------------------------------------------
        # Prepares the structure files of the templates in the output directory. -
        #-------------------------------------------------------------------------

        # Prepares the single chain templates files.
        for mc in self.modeling_clusters_list:
            mc.prepare_single_chains_template_files()

        # Prepares the template complex file.
        if self.multiple_chain_mode:
            list_of_template_complex_files = []
            for t in self.get_template_complex_chains(sorted_by_id=True):
                list_of_template_complex_files.append(os.path.join(self.modeling_dirpath, t.get_structure_file()))
            pmstr.join_pdb_files(list_of_template_complex_files, os.path.join(self.modeling_dirpath, self.template_complex_name))

        #---------------------------------------
        # Prepares input and ouput file paths. -
        #---------------------------------------

        self.pir_file_path = os.path.join(self.modeling_dirpath, self.pir_file_name)


    #--------------------------------------------------------------------------
    # Creates a file with the target-template(s )alignment in the PIR format. -
    #--------------------------------------------------------------------------

    def build_pir_align_file(self):
        """
        This function creates alignments in a PIR format: this is entirely rewritten from the
        original PyMod version. This method should write the sequences as seen by MODELLER.
        """

        # Starts to write the PIR sequence file needed by MODELLER for input.
        pir_align_file_handle = open(self.pir_file_path, "w")


        for modeling_cluster in self.modeling_clusters_list:

            #--------------------------------------------------------------------------------------
            # Prepares the full sequences (with both standard residues and heteroresidues) of the -
            # target and templates.                                                               -
            #--------------------------------------------------------------------------------------

            het_to_use_list = []
            np_hetatm_portions = {}

            for template in modeling_cluster.templates_list:

                # Builds the polymer sequence of the template. It comprises standard amino acid
                # residues and modified residues making part of the polypeptide chain.
                res_count = 0
                res_list = template.get_polymer_residues()
                polymer_seq = ""

                for ali_id, ali_pos in enumerate(template.my_sequence):

                    if ali_pos != "-":

                        res = res_list[res_count]
                        if not res.is_standard_residue():
                            if modeling_cluster.template_options_dict[template]["hetres_dict"][res]:
                                het_to_use_list.append(ali_id)
                            polymer_seq += "."

                        else:
                            polymer_seq += ali_pos

                        res_count += 1

                    else:
                        polymer_seq += ali_pos


                # Adds non-polymer heteroatoms (ligands and water molecules).
                np_hetres_list = template.get_residues(standard=False, ligands=True,
                                                       modified_residues=False,
                                                       water=self.use_water_in_session)
                np_hetatm_portions[template] = np_hetres_list

                self.pir_sequences_dict[template] = {"pir_seq": polymer_seq,
                                                     "modeling_cluster": modeling_cluster}


            # Builds the PIR string for the target sequence.
            polymer_seq = ""

            for ali_id, ali_pos in enumerate(modeling_cluster.target.my_sequence):

                if ali_pos != "-":
                    if ali_id in het_to_use_list and self.use_hetatm_in_session:
                        polymer_seq += "."
                    else:
                        if ali_pos != "X":
                            polymer_seq += ali_pos
                        else:
                            polymer_seq += "A" # Adds an alanine for an undefined residue.
                else:
                    polymer_seq += ali_pos
            self.pir_sequences_dict[modeling_cluster.target] = {"pir_seq": polymer_seq, "modeling_cluster": modeling_cluster}


            #-------------------------------------------------------------------------------
            # Update the sequences by inserting in them ligands and water molecules in the -
            # sequences aligned in PyMod (which miss ligands and water molecules).         -
            #-------------------------------------------------------------------------------

            # Add the non-polymer heteratomic portions (containing ligands and water molecules)
            # for the templates.
            for template_i in modeling_cluster.templates_list:
                for template_j in modeling_cluster.templates_list:
                    np_hetres_list = np_hetatm_portions[template_j]
                    # Add to each template its ligand portions.
                    if template_i is template_j:
                        self.pir_sequences_dict[template_i]["pir_seq"] += self._get_hetres_seq(np_hetres_list)
                    else:
                        self.pir_sequences_dict[template_i]["pir_seq"] += len(self._get_hetres_seq(np_hetres_list))*"-"

            # Add the non-polymer heteratomic portion for the target.
            for template_i in modeling_cluster.templates_list:
                np_hetres_list = np_hetatm_portions[template_i]
                for hetres in np_hetres_list:
                    if modeling_cluster.template_options_dict[template_i]["hetres_dict"][hetres] and self.use_hetatm_in_session:
                        if not hetres.is_water():
                            self.pir_sequences_dict[modeling_cluster.target]["pir_seq"] += "."
                        else:
                            self.pir_sequences_dict[modeling_cluster.target]["pir_seq"] += "w"
                    else:
                        self.pir_sequences_dict[modeling_cluster.target]["pir_seq"] += "-"


        #--------------------------------
        # Write the template sequences. -
        #--------------------------------

        # First write the template complex block (only for multiple chain mode).
        if self.multiple_chain_mode:
            print(">P1;%s" % self.template_complex_modeller_name, file=pir_align_file_handle)
            print("structure:%s:.:.:.:.::::" % self.template_complex_modeller_name, file=pir_align_file_handle) # TODO: (template_code,template_chain,template_chain)
            tc_pir_string = ""
            for template_complex_chain in self.get_template_complex_chains(sorted_by_id=True): # TODO: order them.
                tc_pir_string += self.pir_sequences_dict[template_complex_chain]["pir_seq"] + self.get_chain_separator()
            print(self.get_pir_formatted_sequence(tc_pir_string), file=pir_align_file_handle)

        # Then write the single chain template blocks.
        for modeling_cluster in self.get_modeling_clusters_list(sorted_by_id=True):
            for template in modeling_cluster.get_single_chain_templates():
                # Writes the first line of the template.
                template_code = modeling_cluster.template_options_dict[template]["modeller_name"]
                template_chain = template.get_chain_id()
                print(">P1;%s" % template_code, file=pir_align_file_handle)
                print("structure:%s:.:%s:.:%s::::" % (template_code, template_chain, template_chain), file=pir_align_file_handle)
                sct_pir_string = ""
                for target_chain in self.get_targets_list(sorted_by_id=True):
                    if target_chain != modeling_cluster.target:
                        sct_pir_string += len(self.pir_sequences_dict[target_chain]["pir_seq"])*"-" + self.get_chain_separator() # TODO: adjust separator for single chain modeling.
                    else:
                        sct_pir_string += self.pir_sequences_dict[template]["pir_seq"] + self.get_chain_separator()
                print(self.get_pir_formatted_sequence(sct_pir_string), file=pir_align_file_handle)

        # Finally write the target block.
        print(">P1;%s" % self.modeller_target_name, file=pir_align_file_handle)
        print("sequence:%s:.:.:.:.::::" % self.modeller_target_name, file=pir_align_file_handle)
        targets_pir_string = ""
        for target in self.get_targets_list(sorted_by_id=True):
            targets_pir_string += self.pir_sequences_dict[target]["pir_seq"] + self.get_chain_separator()
        print(self.get_pir_formatted_sequence(targets_pir_string), file=pir_align_file_handle)

        pir_align_file_handle.close()


    def _get_hetres_seq(self, hetres_list):
        hetres_seq = ""
        for hetres in hetres_list:
            if not hetres.is_water():
                hetres_seq += "."
            else:
                hetres_seq += "w"
        return hetres_seq


    def get_pir_formatted_sequence(self,sequence,multi=False):
        """
        Print one block of the PIR alignment file using 60 characters-long lines.
        """
        formatted_sequence = ""
        for s in range(0,len(sequence),60):
            # For all the lines except the last one.
            if (len(sequence) - s) > 60:
                formatted_sequence += sequence[s:s+60] + "\n"
            # For the last line.
            else:
                if not multi:
                    formatted_sequence += sequence[s:]+"*"+"\n"
                else:
                    formatted_sequence += sequence[s:]+"/*"+"\n"
        return formatted_sequence

    def get_chain_separator(self):
        if self.multiple_chain_mode:
            return "/"
        else:
            return ""


    ###########################################################################
    # Actually launches the full 3D modeling protocol.                        #
    ###########################################################################

    def prepare_modelization(self):

        #-----------------------------------------------------------------------------------
        # Takes input supplied by users though the GUI and sets the names of sequences and -
        # template which will be used by MODELLER.                                         -
        #-----------------------------------------------------------------------------------

        try:
            self.get_modeling_options_from_gui()
        except Exception as e:
            title = "Input Error"
            self.pymod.main_window.show_error_message(title, str(e))
            return False

        # Starts the modeling process only if the user has supplied correct parameters.
        try:
            self.check_all_modeling_parameters()
        except Exception as e:
            title = "Input Error"
            self.pymod.main_window.show_error_message(title, str(e))
            return False

        self.set_modeling_options()

        # The modeling window can be destroyed.
        self.modeling_window.destroy()

        # Prepares the directory where MODELLER's output will be generated and moves into it.
        self.prepare_modeling_session_files()

        return True


    def perform_modelization(self):
        try:
            self._perform_modelization()
        except Exception as e:
            self.finish_modeling_session(successful=False, error_message=str(e))


    def launch_modeller(self, in_thread=False):
        """
        Sets the MODELLER environment and script, then launches 3D model building.
        """

        #************************
        # Sets the environment. *
        #************************

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            # Use multiple CPUs, each one will be used to build a model copy.
            if self.modeller_n_jobs > 1:
                j = parallel.job()
                for _job_idx in range(0, self.modeller_n_jobs):
                    j.append(parallel.local_slave())

            modeller.log.verbose()
            if self.modeller_random_seed == self.default_modeller_random_seed:
                env = modeller.environ()
            else:
                env = modeller.environ(rand_seed=self.modeller_random_seed)
            env.io.atom_files_directory = []
            env.io.atom_files_directory.append(".")
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            # This dictionary will store string which will be used to write the modeling script
            # files.
            self.mod_script_dict = {"environment": "",
                                    "hetres": "",
                                    "automodel": {"definition": "",
                                                  "default_patches": "",
                                                  "special_patches": {"definition": "",
                                                                      "multichain": "",
                                                                      "disulfides": ""},
                                                  "special_restraints": "",
                                                  "loop_selection": "",},
                                    "automodel_init": "",
                                    "refinement": "",
                                    "models_indices": "",
                                    "make": "",
                                    "external_post_make": ""}
            self.mod_script_dict["environment"] += "import modeller\n"
            self.mod_script_dict["environment"] += "import modeller.automodel\n"
            if self.modeller_n_jobs > 1:
                self.mod_script_dict["environment"] += "from modeller import parallel\n\n"
                self.mod_script_dict["environment"] += "j = parallel.job()\n"
                self.mod_script_dict["environment"] += "for _job_idx in range(0, %s):\n" % self.modeller_n_jobs
                self.mod_script_dict["environment"] += "    j.append(parallel.local_slave())\n\n"
            else:
                self.mod_script_dict["environment"] += "\n"

            self.mod_script_dict["environment"] += "modeller.log.verbose()\n"
            if self.modeller_random_seed == self.default_modeller_random_seed:
                self.mod_script_dict["environment"] += "env = modeller.environ()\n"
            else:
                self.mod_script_dict["environment"] += "env = modeller.environ(rand_seed=%s)\n" % self.modeller_random_seed
            if not self.run_modeller_internally:
                env = None
            self.mod_script_dict["environment"] += "env.io.atom_files_directory = []\n"
            self.mod_script_dict["environment"] += "env.io.atom_files_directory.append('.')\n"
        #------------------------------------------------------------


        #**************************************
        # Sets heteroatoms and water options. *
        #**************************************

        # If the user wants to include hetero-atoms and water molecules. The 'hetatm' flag is always
        # set to 'True' in homology modeling.
        if self.modeling_mode == "homology_modeling" or self.use_hetatm_in_session:

            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                env.io.hetatm = True
            #--------------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["hetres"] += "env.io.hetatm = True\n"
            #--------------------------------------------------------

            # Use water only if the user chose to include water molecules from some template.
            if self.use_water_in_session:

                # Internal ------------------------------------------
                if self.run_modeller_internally:
                    env.io.water = True
                #----------------------------------------------------

                # External ------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["hetres"] += "env.io.water = True\n"
                #----------------------------------------------------


        #**********************************************************
        # Define the name of the base automodel class to be used. *
        #**********************************************************

        if self.modeling_mode == "homology_modeling":
            session_automodel_class = modeller.automodel.automodel
            session_automodel_class_name = "automodel"

        elif self.modeling_mode == "loop_refinement":
            if self.loop_refinement_class_name == "loopmodel":
                session_automodel_class = modeller.automodel.loopmodel
                session_automodel_class_name = "loopmodel"
            elif self.loop_refinement_class_name == "DOPE loopmodel":
                session_automodel_class = modeller.automodel.dope_loopmodel
                session_automodel_class_name = "dope_loopmodel"


        #*******************************************************
        # Creates a file with the alignment in the PIR format. *
        #*******************************************************

        if self.modeling_mode == "homology_modeling":
            self.build_pir_align_file()


        #*******************************************************************
        # Defines a custom class to use some additional MODELLER features. *
        #*******************************************************************
        # This class is going to be used to build the "a" object used to perform the actual
        # homology modelization. It is going to inherit everything from the automodel class
        # but is going to have dynamically redifined routines to make it possible to:
        #   - include user defined disulfide bridges in the model
        #   - exclude template disulfide bridges in the model
        #   - build multichain models with symmetries restraints
        #   - rename the chains in multichain models

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            self.mod_script_dict["automodel"]["definition"] += "class MyModel(modeller.automodel.%s):\n" % session_automodel_class_name
        #------------------------------------------------------------


        #*******************************
        # Template-derived disulfides. *
        #*******************************

        if self.modeling_mode == "homology_modeling":

            # If there are some targets with at least two CYS residues and also some templates with
            # disulfides, it decides whether to use template disulfides or to let the CYS residues
            # in a "reduced" state.
            if self.check_targets_with_cys() and self.check_structures_with_disulfides():

                # If the user chose to use templates disulfides bridges MODELLER will automatically use
                # the patch_ss_templates() method of the automodel class.
                if self.use_template_dsb:
                    pass

                # Don't use template dsbs: this will not create any dsbs in the model by not disulfide
                # patching and leave the model CYS residues that in the template are engaged in a dsb in
                # a "reduced" state.
                else:

                    # External ------------------------------------------
                    if self.write_modeller_script_option:
                        self.mod_script_dict["automodel"]["default_patches"] += "    def default_patches(self, aln):\n"
                        self.mod_script_dict["automodel"]["default_patches"] += "        pass\n"
                    #----------------------------------------------------


        #***********************************************************************************
        # Part for multichain models and user defined disulfide bridges, which requires to *
        # override the special_patches() method.                                           *
        #***********************************************************************************

        # External --------------------------------------------------
        if self.write_modeller_script_option:

            self.mod_script_dict["automodel"]["special_patches"]["definition"] += "    def special_patches(self, aln):\n"

            #--------------------------------------------------------------------------------
            # Renumber the residues in the new chains starting from 1. When Modeller builds -
            # a multichain model it doesn't restart to count residues from 1 when changing  -
            # chain. The following code renumbers the residues in the correct way.          -
            #--------------------------------------------------------------------------------

            if self.multiple_chain_mode:

                if self.modeling_mode == "homology_modeling":

                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        # Renumber the residues of the chains of multichain models starting from 1.\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        count_dictionary = dict((chain.name, 1) for chain in self.chains)\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        for chain in self.chains:\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "            for residue in chain.residues:\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "                residue.num = '%d' % (count_dictionary[chain.name])\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "                count_dictionary[chain.name] += 1\n"


            #-------------------------------------------------------------
            # Informs Modeller on how to build custom disulfide bridges. -
            #-------------------------------------------------------------

            if self.check_targets_with_cys():

                # If the user wants to use some custom dsb.
                if self.use_user_defined_dsb:

                    # Gets the list of user defined dsb for each modeling cluster (if the
                    # target of the modeling cluster doesn't have at least two cys residues
                    # it will have an [] empty list).
                    for (mci, mc) in enumerate(self.modeling_clusters_list):
                        for dsb in self.all_user_defined_dsb[mci]:
                            # For example "CYS 321".
                            cys1 = dsb[0][4:] # dsb[0][3:]
                            cys2 = dsb[1][4:] # dsb[1][3:]

                            # Redefine the routine to include user defined dsb.
                            # NOTE. If a bridge has the same cys the following error will ensue:
                            # <class '_modeller.ModellerError'>: unqang__247E> Internal error
                            if self.modeling_mode == "homology_modeling":
                                if self.multiple_chain_mode:
                                    chain = mc.block_index
                                    self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.chains[%s].residues['%s'], self.chains[%s].residues['%s']))\n" % (chain, cys1, chain, cys2)
                                else:
                                    self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.residues['%s'], self.residues['%s']))\n" % (cys1, cys2)
                            elif self.modeling_mode == "loop_refinement":
                                chain = mc.target.get_chain_id()
                                self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.chains['%s'].residues['%s'], self.chains['%s'].residues['%s']))\n" % (chain, cys1, chain, cys2)

                # If the user wants MODELLER to build automatically the dsb.
                if self.use_auto_dsb:
                    # Adds disulfides bridges for cys that are sufficently close.
                    self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch_ss()\n"
        #------------------------------------------------------------


        #**************************************************************************
        # Apply symmetry restraints to target chains that have the same sequence. *
        #**************************************************************************

        if self.multiple_chain_mode and len(self.list_of_symmetry_restraints) > 0:

            # External ----------------------------------------------
            if self.write_modeller_script_option:

                # Constrain chains to be identical (but only restrain the C-alpha atoms, to reduce
                # the number of interatomic distances that need to be calculated).
                self.mod_script_dict["automodel"]["special_restraints"] += "    def special_restraints(self, aln):\n"
                for si, symmetry_restraints_group in enumerate(self.list_of_symmetry_restraints):
                    # MODELLER's API only allows to apply symmetry restraints to couple of chains.
                    # Therefore if we want to apply symmetry restraints to chains ("A", "B", "C"),
                    # we will first apply symmetry restraints to ("A", "B") and then to ("B", "C").
                    self.mod_script_dict["automodel"]["special_restraints"] += "        # Symmetry restraints group n. %d.\n" % (si+1)
                    for s in symmetry_restraints_group:
                        # In homology modeling, the chain ids for symmetry restraints are integers
                        # starting from 0.
                        if self.modeling_mode == "homology_modeling":
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s1 = modeller.selection(self.chains[%s]).only_atom_types('CA')\n" % (s[0])
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s2 = modeller.selection(self.chains[%s]).only_atom_types('CA')\n" % (s[1])
                        # In loop refinement, the chain ids for symmetry restraints are the chain ID
                        # letters.
                        elif self.modeling_mode == "loop_refinement":
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s1 = modeller.selection(self.chains['%s']).only_atom_types('CA')\n" % (s[0])
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s2 = modeller.selection(self.chains['%s']).only_atom_types('CA')\n" % (s[1])
                        self.mod_script_dict["automodel"]["special_restraints"] += "        self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "\n"

                # Report on symmetry violations greater than 1A after building each model.
                self.mod_script_dict["automodel"]["special_restraints"] += "    def user_after_single_model(self):\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "        self.restraints.symmetry.report(1.0)\n\n"
            #--------------------------------------------------------


        #******************
        # Loop selection. *
        #******************

        if self.modeling_mode == "loop_refinement":

            # Definition of loop residues.
            loop_intervals = []
            for mcl in self.modeling_clusters_list:
                loop_intervals.extend(mcl.loops_list)


            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["automodel"]["special_restraints"] += "    def select_loop_atoms(self):\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "        return modeller.selection(\n"
                for lt in loop_intervals:
                    self.mod_script_dict["automodel"]["special_restraints"] += "            self.residue_range('%s:%s', '%s:%s'),\n" % (lt[0], lt[2], lt[1], lt[2])
                self.mod_script_dict["automodel"]["special_restraints"] += "        )\n\n"
            #--------------------------------------------------------


            # Loop starting conformation.
            if self.loop_initial_conformation == "Linearized":
                pass

            elif self.loop_initial_conformation == "Original":

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["automodel"]["special_restraints"] += "    def build_ini_loop(self, atmsel):\n"
                    self.mod_script_dict["automodel"]["special_restraints"] += "        pass\n"
                #--------------------------------------------------------


        #**************************************************************
        # Creates the "automodel" object to perform the modelization. *
        #**************************************************************

        # Internal --------------------------------------------------
        if self.run_modeller_internally:

            # Writes a module in the current modeling directory.
            self.write_modeller_script(write_script=False, write_lib=True)
            # Imports the module.
            sys.path.append(self.modeling_dirpath)
            _session_mod_lib = importlib.import_module(self.modeling_lib_name)
            # Reloads it, so that the module from previous modeling sessions (if present) is
            # refreshed with the new one.
            _session_mod_lib = importlib.reload(_session_mod_lib)
            MyModel = _session_mod_lib.MyModel
            sys.path.remove(self.modeling_dirpath)

            # Standard homology modeling.
            if self.modeling_mode == "homology_modeling":
                if not self.multiple_chain_mode:
                    assess_methods = (modeller.automodel.assess.GA341, modeller.automodel.assess.DOPE)
                else:
                    assess_methods = (modeller.automodel.assess.DOPE, )
                a = MyModel(env,
                            alnfile=self.pir_file_name,                                        # alignment filename
                            knowns=tuple([str(tmpn) for tmpn in self.all_templates_namelist]), # tuple(self.all_templates_namelist),                # codes of the templates
                            sequence=str(self.modeller_target_name),                           # code of the target
                            assess_methods=assess_methods)

            # Loop modeling.
            elif self.modeling_mode == "loop_refinement":
                loop_assess_methods = (modeller.automodel.assess.DOPE, ) # soap_loop.Scorer()
                a = MyModel(env,
                            sequence=str(self.modeller_target_name),
                            inimodel=self.loop_refinement_starting_model,       # initial model of the target
                            loop_assess_methods=loop_assess_methods)            # assess loops with statistical potentials

        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            if self.modeling_mode == "homology_modeling":
                self.mod_script_dict["automodel_init"] += "a = MyModel(\n"
                self.mod_script_dict["automodel_init"] += "    env,\n"
                self.mod_script_dict["automodel_init"] += "    alnfile='%s',\n" % self.pir_file_name
                self.mod_script_dict["automodel_init"] += "    knowns=%s,\n" % repr(tuple([str(tmpn) for tmpn in self.all_templates_namelist]))
                self.mod_script_dict["automodel_init"] += "    sequence='%s',\n" % (str(self.modeller_target_name))
                if not self.multiple_chain_mode:
                    self.mod_script_dict["automodel_init"] += "    assess_methods=(modeller.automodel.assess.GA341, modeller.automodel.assess.DOPE))\n"
                else:
                    self.mod_script_dict["automodel_init"] += "    assess_methods=(modeller.automodel.assess.DOPE, ))\n"

            elif self.modeling_mode == "loop_refinement":
                self.mod_script_dict["automodel_init"] += "a = MyModel(\n"
                self.mod_script_dict["automodel_init"] += "    env,\n"
                self.mod_script_dict["automodel_init"] += "    sequence='%s',\n" % self.modeller_target_name
                self.mod_script_dict["automodel_init"] += "    inimodel='%s',\n" % self.loop_refinement_starting_model
                self.mod_script_dict["automodel_init"] += "    loop_assess_methods=(modeller.automodel.assess.DOPE,))\n"
        #------------------------------------------------------------


        #************************************************
        # Sets the level of refinment and optimization. *
        #************************************************

        #---------------------
        # Homology modeling. -
        #---------------------

        if self.modeling_mode == "homology_modeling":

            if self.optimization_level == "Low":
                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    # Low VTFM optimization:
                    a.library_schedule = modeller.automodel.autosched.very_fast
                    # Low MD optimization:
                    a.md_level = modeller.automodel.refine.very_fast
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.very_fast\n"
                    self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.very_fast\n"
                #--------------------------------------------------------

            elif self.optimization_level == "Default":
                # a.library_schedule = modeller.automodel.autosched.normal
                # a.max_var_iterations = 200
                # a.md_level = modeller.automodel.refine.very_fast
                # a.repeat_optimization = 1
                # a.max_molpdf = 1e7
                pass

            elif self.optimization_level == "Mid":
                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    # Thorough VTFM optimization:
                    a.library_schedule = modeller.automodel.autosched.normal
                    a.max_var_iterations = 300
                    # Mid MD optimization:
                    a.md_level = modeller.automodel.refine.fast
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.fast\n"
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = 300\n"
                    self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.fast\n"
                #--------------------------------------------------------

            elif self.optimization_level == "High":
                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    # Very thorough VTFM optimization:
                    a.library_schedule = modeller.automodel.autosched.slow
                    a.max_var_iterations = 300
                    # Thorough MD optimization:
                    a.md_level = modeller.automodel.refine.slow
                    # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
                    a.repeat_optimization = 2
                    a.max_molpdf = 1e8 # 1e6
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.slow\n"
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = 300\n"
                    self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.slow\n"
                    self.mod_script_dict["refinement"] += "a.repeat_optimization = 2\n"
                    self.mod_script_dict["refinement"] += "a.max_molpdf = 1e8\n"
                #--------------------------------------------------------

            elif self.optimization_level == "Approximate":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    # Reduces the amount of homology-derived distance restraints to speed up the
                    # optimization.
                    a.max_ca_ca_distance = 10.0
                    a.max_n_o_distance = 6.0
                    a.max_sc_mc_distance = 5.0
                    a.max_sc_sc_distance = 4.5
                    # Few CG iterations.
                    a.max_var_iterations = 50
                    a.library_schedule = modeller.automodel.autosched.fastest
                    # No MD phase.
                    a.md_level = None
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.max_ca_ca_distance = 10.0\n"
                    self.mod_script_dict["refinement"] += "a.max_n_o_distance = 6.0\n"
                    self.mod_script_dict["refinement"] += "a.max_sc_mc_distance = 5.0\n"
                    self.mod_script_dict["refinement"] += "a.max_sc_sc_distance = 4.5\n"
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = 50\n"
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.fastest\n"
                    self.mod_script_dict["refinement"] += "a.md_level = None\n"
                #--------------------------------------------------------

            elif self.optimization_level == "Custom":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    a.max_var_iterations = self.custom_max_cg_iterations
                    a.library_schedule = self.get_modeller_schedule_obj("cg", "internal", self.custom_vtfm_schedule)
                    a.md_level = self.get_modeller_schedule_obj("md", "internal", self.custom_md_schedule)
                    a.repeat_optimization = self.custom_repeat_optimization
                    a.max_molpdf = 10**self.custom_max_obj_func_value
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = %s\n" % self.custom_max_cg_iterations
                    self.mod_script_dict["refinement"] += "a.library_schedule = %s\n" % self.get_modeller_schedule_obj("cg", "external", self.custom_vtfm_schedule)
                    self.mod_script_dict["refinement"] += "a.md_level = %s\n" % self.get_modeller_schedule_obj("md", "external", self.custom_md_schedule)
                    self.mod_script_dict["refinement"] += "a.repeat_optimization = %s\n" % self.custom_repeat_optimization
                    self.mod_script_dict["refinement"] += "a.max_molpdf = 10**%s\n" % self.custom_max_obj_func_value
                #--------------------------------------------------------


        #-------------------
        # Loop refinement. -
        #-------------------

        elif self.modeling_mode == "loop_refinement":

            if self.optimization_level == "Low":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    a.loop.md_level = modeller.automodel.refine.very_fast
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = modeller.automodel.refine.very_fast\n"
                #--------------------------------------------------------

            elif self.optimization_level == "Mid":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    a.loop.md_level = modeller.automodel.refine.fast
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = modeller.automodel.refine.fast\n"
                #--------------------------------------------------------

            elif self.optimization_level == "Default":
                # a.loop.md_level = modeller.automodel.refine.slow
                pass

            elif self.optimization_level == "High":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    a.loop.md_level = modeller.automodel.refine.very_slow
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = modeller.automodel.refine.very_slow\n"
                #--------------------------------------------------------

            elif self.optimization_level == "Approximate":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    a.loop.md_level = None
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = None\n"
                #--------------------------------------------------------

            elif self.optimization_level == "Custom":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:
                    a.loop.max_var_iterations = self.custom_max_cg_iterations
                    a.loop.md_level = self.get_modeller_schedule_obj("md", "internal", self.custom_md_schedule)
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.max_var_iterations = %s\n" % self.custom_max_cg_iterations
                    self.mod_script_dict["refinement"] += "a.loop.md_level = %s\n" % self.get_modeller_schedule_obj("md", "external", self.custom_md_schedule)
                #--------------------------------------------------------


        #************************************
        # Customize the objective function. *
        #************************************

        if self.use_custom_obj_func in ("customize", "altmod"):

            #---------------------
            # Homology modeling. -
            #---------------------

            if self.modeling_mode == "homology_modeling":

                w_dope = self.custom_w_dope if self.use_dope_in_obj_func else 0

                # Internal ----------------------------------------------
                if self.run_modeller_internally:

                    # Sets the weights of the objective function.
                    a.env.schedule_scale = modeller.physical.values(default=1.0,
                                                                    # Statistical potential terms.
                                                                    nonbond_spline=w_dope,
                                                                    # Distance restraints terms.
                                                                    ca_distance=self.custom_w_hddr,
                                                                    n_o_distance=self.custom_w_hddr,
                                                                    sd_mn_distance=self.custom_w_hddr,
                                                                    sd_sd_distance=self.custom_w_hddr,
                                                                    # Dihedral restraints terms.
                                                                    phi_dihedral=self.custom_w_hdar,
                                                                    psi_dihedral=self.custom_w_hdar,
                                                                    chi1_dihedral=self.custom_w_hdar,
                                                                    chi2_dihedral=self.custom_w_hdar,
                                                                    chi3_dihedral=self.custom_w_hdar,
                                                                    chi4_dihedral=self.custom_w_hdar,
                                                                    phi_psi_dihedral=self.custom_w_hdar,
                                                                    # CHARM22 terms.
                                                                    bond=self.custom_w_charmm22,
                                                                    angle=self.custom_w_charmm22,
                                                                    dihedral=self.custom_w_charmm22,
                                                                    improper=self.custom_w_charmm22,
                                                                    # Soft-sphere terms.
                                                                    soft_sphere=self.custom_w_soft_sphere)

                    edat = a.env.edat
                    edat.contact_shell = self.custom_nb_cutoff

                    # Allow calculation of statistical (dynamic_modeller) potential.
                    if self.use_dope_in_obj_func:
                        edat.dynamic_modeller = True
                        # Group restraints.
                        gprsr = modeller.group_restraints(a.env, classes='$(LIB)/atmcls-mf.lib', parameters='$(LIB)/dist-mf.lib')
                        a.group_restraints = gprsr
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += """

# Sets the weights of the objective function.
a.env.schedule_scale = modeller.physical.values(default=1.0,
                                                # Statistical potential terms.
                                                nonbond_spline=%s,
                                                # Distance restraints terms.
                                                ca_distance=%s,
                                                n_o_distance=%s,
                                                sd_mn_distance=%s,
                                                sd_sd_distance=%s,
                                                # Dihedral restraints terms.
                                                phi_dihedral=%s,
                                                psi_dihedral=%s,
                                                chi1_dihedral=%s,
                                                chi2_dihedral=%s,
                                                chi3_dihedral=%s,
                                                chi4_dihedral=%s,
                                                phi_psi_dihedral=%s,
                                                # CHARM22 terms.
                                                bond=%s,
                                                angle=%s,
                                                dihedral=%s,
                                                improper=%s,
                                                # Soft-sphere terms.
                                                soft_sphere=%s)

# Non-bonded interactions cutoff.
edat = a.env.edat
edat.contact_shell = %s
""" % (w_dope, self.custom_w_hddr, self.custom_w_hddr, self.custom_w_hddr, self.custom_w_hddr, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_soft_sphere, self.custom_nb_cutoff)

                    if self.use_dope_in_obj_func:
                         self.mod_script_dict["refinement"] += "# Adds DOPE terms to the objective function.\n"
                         self.mod_script_dict["refinement"] += "edat.dynamic_modeller = True\n"
                         self.mod_script_dict["refinement"] += "gprsr = modeller.group_restraints(a.env, classes='$(LIB)/atmcls-mf.lib', parameters='$(LIB)/dist-mf.lib')\n"
                         self.mod_script_dict["refinement"] += "a.group_restraints = gprsr\n"
                #--------------------------------------------------------


            #-------------------
            # Loop refinement. -
            #-------------------

            elif self.modeling_mode == "loop_refinement":

                # Internal ----------------------------------------------
                if self.run_modeller_internally:

                    # Sets the weights of the objective function.
                    a.loop.env.schedule_scale = modeller.physical.values(default=1.0,
                                                                         # Statistical potential terms.
                                                                         nonbond_spline=self.custom_w_stat_pot_loop,
                                                                         # Dihedral restraints terms.
                                                                         phi_dihedral=self.custom_w_hdar,
                                                                         psi_dihedral=self.custom_w_hdar,
                                                                         chi1_dihedral=self.custom_w_hdar,
                                                                         chi2_dihedral=self.custom_w_hdar,
                                                                         chi3_dihedral=self.custom_w_hdar,
                                                                         chi4_dihedral=self.custom_w_hdar,
                                                                         phi_psi_dihedral=self.custom_w_hdar,
                                                                         # CHARM22 terms.
                                                                         bond=self.custom_w_charmm22,
                                                                         angle=self.custom_w_charmm22,
                                                                         dihedral=self.custom_w_charmm22,
                                                                         improper=self.custom_w_charmm22,
                                                                         # Soft-sphere terms.
                                                                         soft_sphere=self.custom_w_soft_sphere,
                                                                         # Lennard-Jones terms.
                                                                         lennard_jones=self.custom_w_lj,
                                                                         # GBSA terms.
                                                                         gbsa=self.custom_w_gbsa,
                                                                         )

                    edat = a.loop.env.edat
                    edat.contact_shell = self.custom_nb_cutoff
                #--------------------------------------------------------

                # External ----------------------------------------------
                if self.write_modeller_script_option:

                    self.mod_script_dict["refinement"] += """

# Sets the weights of the objective function.
a.loop.env.schedule_scale = modeller.physical.values(default=1.0,
                                                     # Statistical potential terms.
                                                     nonbond_spline=%s,
                                                     # Dihedral restraints terms.
                                                     phi_dihedral=%s,
                                                     psi_dihedral=%s,
                                                     chi1_dihedral=%s,
                                                     chi2_dihedral=%s,
                                                     chi3_dihedral=%s,
                                                     chi4_dihedral=%s,
                                                     phi_psi_dihedral=%s,
                                                     # CHARM22 terms.
                                                     bond=%s,
                                                     angle=%s,
                                                     dihedral=%s,
                                                     improper=%s,
                                                     # Soft-sphere terms.
                                                     soft_sphere=%s,
                                                     # Lennard-Jones terms.
                                                     lennard_jones=%s,
                                                     # GBSA terms.
                                                     gbsa=%s,
                                                     )

# Non-bonded interactions cutoff.
edat = a.loop.env.edat
edat.contact_shell = %s
""" % (self.custom_w_stat_pot_loop, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_soft_sphere, self.custom_w_lj, self.custom_w_gbsa, self.custom_nb_cutoff)
                #--------------------------------------------------------


        #*********************************************************************
        # Determines how many models to build and actually build the models. *
        #*********************************************************************

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            if self.modeling_mode == "homology_modeling":
                self.mod_script_dict["models_indices"] += "a.starting_model = %s\n" % self.starting_model_number
                self.mod_script_dict["models_indices"] += "a.ending_model = %s\n" % self.ending_model_number

            elif self.modeling_mode == "loop_refinement":
                self.mod_script_dict["models_indices"] += "a.loop.starting_model = %s\n" % self.starting_model_number
                self.mod_script_dict["models_indices"] += "a.loop.ending_model = %s\n" % self.ending_model_number

            if self.modeller_n_jobs > 1:
                self.mod_script_dict["make"] += "a.use_parallel_job(j)\n"
            self.mod_script_dict["make"] += "a.make()\n"

            # Saves an output file that will be read by PyMod when MODELLER is executed externally.
            self.write_modeller_script(write_script=True, write_lib=False)

        if not self.run_modeller_internally:
            raise Exception("Not implemented anymore.")
        #------------------------------------------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            if self.modeling_mode == "homology_modeling":
                a.starting_model = self.starting_model_number # index of the first model
                a.ending_model = self.ending_model_number # index of the last model
            elif self.modeling_mode == "loop_refinement":
                a.loop.starting_model = self.starting_model_number
                a.loop.ending_model = self.ending_model_number

            if self.modeller_n_jobs > 1:
                # Use a parallel job for model building.
                a.use_parallel_job(j)

            # This is the method that launches the model building phase.
            a.make()
        #------------------------------------------------------------

        # Saves the 'automodel' object in a pickle file that can be loaded in the main thread.
        if in_thread:
            with open(self.modeller_vars_filename, "wb") as v_fh:
                pickle.dump(a, v_fh)
        # Directly returns the 'automodel' object.
        else:
            return env, a


    def load_modeller_vars(self):
        with open(self.modeller_vars_filename, "rb") as v_fh:
            a = pickle.load(v_fh)
        return None, a


    def _perform_modelization(self):

        ###########################################################################################
        # Sets the MODELLER run and build 3D models.                                              #
        ###########################################################################################

        # Launches in the main thread the method in which the 3D models are actually built.
        if not self.pymod.use_protocol_threads:
            env, a = self.launch_modeller(in_thread=False)

        # Launches in a separate thread the method.
        else:
            label_text = "Running 3D model building. Please wait for the process to complete..."
            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.launch_modeller,
                                            args=(True, ),
                                            title="Building 3D models with MODELLER.",
                                            label_text=label_text,
                                            # A thread with MODELLER can not exit safely.
                                            lock=True,
                                            wait_start=1.0, wait_end=1.0, wait_close=1.0,
                                            lock_title=self.pymod.modeller_lock_title,
                                            lock_message=self.pymod.modeller_lock_message,
                                            stdout_filepath=os.path.join(self.modeling_dirpath, "modeller_log.txt"))
            p_dialog.exec_()

            # Loads the MODELLER 'automodel' object saved on a pickle file.
            with open(os.devnull, "w") as n_fh: # Silence some MODELLER output.
                original_stdout = sys.stdout
                sys.stdout = n_fh
                env, a = self.load_modeller_vars()
                sys.stdout = original_stdout


        ###########################################################################################
        # Finishes to set options for MODELLER and returns back to the PyMod projects directory.  #
        ###########################################################################################

        os.chdir(self.pymod.current_project_dirpath)

        #-----------------------------------------------------------------------------------
        # Cycles through all models built by MODELLER to import them into PyMod and PyMOL. -
        #-----------------------------------------------------------------------------------

        # Gets the ouput of MODELLER in a list of dictionaries.
        automodel_output = self._get_modeller_outputs(a)
        # Filters the output by removing models that were not built (when using the standard
        # MODELLER algorithm, a modeling failure is very rare).
        automodel_output = [m for m in automodel_output if m["failure"] is None]

        # Insert here code to perform additional actions on the models.
        # ...

        # If the list of correctly built models is empty, terminates the modeling session.
        if not automodel_output:
            error_message = ("No models could be correctly built by MODELLER."
                             " Terminating the modeling session.")
            self.finish_modeling_session(successful=False, error_message=error_message)
            return None


        for model_file_number, model in enumerate(automodel_output):

            #--------------------------------------------------------
            # Builds PyMod elements for each of the model's chains. -
            #--------------------------------------------------------

            # Gets the file name generated by MODELLER (stored in model['name']).
            model_file_full_path = os.path.join(self.modeling_dirpath, model['name'])

            # Builds a new file name for the model.
            pymod_model_name = self._get_model_filename()
            # Parses the PDB file of the model.
            parsed_model_file = pmstr.Parsed_model_pdb_file(self.pymod, model_file_full_path,
                                                            output_directory=self.pymod.structures_dirpath,
                                                            new_file_name=pymod_model_name,
                                                            model_root_name=self.modeller_target_name)
            current_model_chains_elements = []
            # Multiple chain models will yield multiple PyMod elements (one for each chain, thus one
            # for each modeling cluster).
            for chain_number, (model_element, modeling_cluster) in enumerate(zip(parsed_model_file.get_pymod_elements(), self.get_modeling_clusters_list(sorted_by_id=True))):
                # Add the new element to PyMod.
                self.pymod.add_element_to_pymod(model_element, color=self.get_model_color(chain_number, self.multiple_chain_mode))
                modeling_cluster.model_elements_list.append(model_element)
                current_model_chains_elements.append(model_element)
                # Gets the aligned sequence of the original target element.
                original_target_aligned_sequence = modeling_cluster.target.my_sequence
                '''
                # Substitute the first model with the target element.
                if model_file_number == 0 and modeling_cluster.target.models_count == 0 and not modeling_cluster.target.has_structure():
                    # self.original_target_aligned_sequence = modeling_cluster.target.my_sequence
                    self.pymod.replace_element(old_element=modeling_cluster.target, new_element=model_element)
                    modeling_cluster.target = model_element
                '''
                # Adds gaps to the various copies of the target sequencs.
                model_element.trackback_sequence(original_target_aligned_sequence)


            #------------------------------------------
            # Superpose models to templates in PyMOL. -
            #------------------------------------------

            if self.modeling_mode == "homology_modeling":

                # Just superpose the model's chain to the first template.
                if self.superpose_to_templates and not self.multiple_chain_mode:
                    super_template_selector = self.modeling_clusters_list[0].templates_list[0].get_pymol_selector()
                    # Builds only a selector for the first and only chain models.
                    for mod_e in current_model_chains_elements:
                        self.superpose_in_pymol(mod_e.get_pymol_selector(), super_template_selector)

                # Superposing for multichain modeling is more complex, and follows a different strategy.
                elif self.superpose_to_templates and self.multiple_chain_mode:
                    # Loads the full template complex file in PyMOL when the first model is loaded.
                    if model_file_number == 0:
                        cmd.load(os.path.join(self.modeling_dirpath, self.template_complex_name), self.tc_temp_pymol_name)
                        # Superpose each separated chain of the template complex to the corresponding
                        # chains of the full template complex.
                        for mc in self.get_modeling_clusters_list(sorted_by_id=True):
                            self.superpose_in_pymol(mc.get_template_complex_chain().get_pymol_selector(),          # Single template complex chain selector.
                                                    "%s and chain %s" % (self.tc_temp_pymol_name, mc.tc_chain_id), # Same chain of the full template complex structure.
                                                    save_superposed_structure=False)
                    # Loads the full model complex file.
                    cmd.load(model_file_full_path, self.mc_temp_pymol_name)
                    # Superpose the full model complex file on the template complex using PyMOL.
                    self.superpose_in_pymol(self.mc_temp_pymol_name, self.tc_temp_pymol_name, save_superposed_structure=False)
                    # Saves the new superposed file in the structures directory.
                    # cmd.save(os.path.join(self.pymod.structures_dirpath, model['name']), self.mc_temp_pymol_name)
                    # Superpose single model chains to the correspondig one of the full model complex.
                    for mod_e in current_model_chains_elements:
                        self.superpose_in_pymol(mod_e.get_pymol_selector(),
                                                "%s and chain %s" % (self.mc_temp_pymol_name, mod_e.get_chain_id()),
                                                save_superposed_structure=False)
                    # Cleans up after having superposed the last multiple chain model.
                    if model_file_number == len(a.outputs) - 1:
                        cmd.delete(self.mc_temp_pymol_name)
                        cmd.delete(self.tc_temp_pymol_name)

            # Increases the models count.
            self.increase_model_number()

        #------------------------------------------------------------
        # Quality assessment of the models.                         -
        #------------------------------------------------------------

        # Launches in the main thread the quality assessment.
        if not self.pymod.use_protocol_threads:
            self.perform_quality_assessment(env, automodel_output)

        # Launches in a separate thread the assessement.
        else:
            label_text = "Running 3D model quality assessment. Please wait for the process to complete..."
            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.perform_quality_assessment,
                                            args=(env, automodel_output),
                                            wait_start=2, wait_end=1.0, wait_close=1.0,
                                            title="Assessing the quality of 3D models.",
                                            label_text=label_text,
                                            lock=True,
                                            lock_title=self.pymod.modeller_lock_title,
                                            lock_message=self.pymod.modeller_lock_message,
                                            stdout_filepath=os.path.join(self.modeling_dirpath, "modeller_assessment_log.txt"))
            p_dialog.exec_()


        #------------------------------------------------------------
        # Finishes the modeling process.                            -
        #------------------------------------------------------------

        # Adds the information of this new modeling session to PyMod.
        self.pymod.modeling_session_list.append(self.current_modeling_session)

        # Finally shows the table and the previously built DOPE profile comprising DOPE curves of
        # every model and templates.
        self.pymod.show_table(**self.current_modeling_session.assessment_table_data)
        show_dope_plot(dope_plot_data=self.current_modeling_session.dope_profile_data,
                       parent_window=self.pymod.get_qt_parent(),
                       pymod=self.pymod)

        # Completes the process.
        self.finish_modeling_session(successful=True)


    def perform_quality_assessment(self, env, automodel_output):

        # Starts to build the 'current_modeling_session' which will be used to build a new item on
        # the 'Models' submenu on the main window.
        self.current_modeling_session = Modeling_session_information(session_id=self.pymod.performed_modeling_count + 1,
                                                                     automodel_output=automodel_output,
                                                                     modeling_directory_path=self.modeling_dirpath)

        #------------------------------------------------------------------------------------------
        # Create a DOPE profile plot for all models built in this by computing their DOPE scores. -
        #------------------------------------------------------------------------------------------

        self.all_assessed_structures_list = []
        session_dope_protocol = DOPE_assessment(self.pymod, output_directory=self.modeling_dirpath)
        session_dope_protocol.remove_temp_files = False

        # Actually computes the DOPE profiles the templates and of the models.
        for mc in self.get_modeling_clusters_list(sorted_by_id=True):
            # Templates.
            if self.modeling_mode == "homology_modeling":
                for template in mc.templates_list:
                    session_dope_protocol.add_element(template)
                    self.all_assessed_structures_list.append(template)
            # Initial models for loop refinement.
            elif self.modeling_mode == "loop_refinement":
                session_dope_protocol.add_element(mc.target)
                self.all_assessed_structures_list.append(mc.target)
            # Models.
            for model_element in mc.model_elements_list:
                session_dope_protocol.add_element(model_element)
                self.all_assessed_structures_list.append(model_element)
        session_dope_protocol.compute_all_dopes(env=env)
        session_dope_protocol.assign_dope_items()

        # This for cycle is used to add extra 'None' values in multiple chains profiles. In this way
        # if, for example, there is a model with chains 'A' and 'B', in the plot the profile of
        # chain 'B' will be put just after the end of the profile of chain 'A'.
        alignment_lenght = 0
        for mc in self.get_modeling_clusters_list(sorted_by_id=True):
            # Computes the DOPE profile of the templates.
            if self.modeling_mode == "homology_modeling":
                for template in mc.templates_list:
                    template_dope_data = session_dope_protocol.prepare_dope_plot_data([template], start_from=alignment_lenght, mode="multiple")
                    self.current_modeling_session.add_template_dope_data(template_dope_data[0])
            elif self.modeling_mode == "loop_refinement":
                target_dope_data = session_dope_protocol.prepare_dope_plot_data([mc.target], start_from=alignment_lenght, mode="multiple")
                self.current_modeling_session.add_template_dope_data(target_dope_data[0])

            # Computes the DOPE profile of the models.
            for model_index, model_element in enumerate(mc.model_elements_list):
                model_dope_data = session_dope_protocol.prepare_dope_plot_data([model_element], start_from=alignment_lenght, mode="multiple")
                self.current_modeling_session.add_model_dope_data(model_dope_data[0], model_index)
            alignment_lenght += len(mc.target.my_sequence)

        #------------------------------------------------------------------------------------
        # Gets the objective function and DOPE scores values for each full model (the model -
        # comprising all the chains) built.                                                 -
        #------------------------------------------------------------------------------------

        # This will also save in the modeling directory the DOPE profile files of the models and
        # templates.
        session_assessment_data = []

        # Initialize a new environment for SOAP scoring.
        if self.multiple_chain_mode:
            soap_env = modeller.environ()
            soap_env.libs.topology.read(file='$(LIB)/top_heav.lib')
            soap_env.libs.parameters.read(file='$(LIB)/par.lib')

        # Score each model copy.
        for decoy_index in range(0, len(automodel_output)):
            model = automodel_output[decoy_index]
            fmo = self.current_modeling_session.full_models[decoy_index]

            obj_funct_value = self.get_model_objective_function_value(model)
            dope_score = self.get_model_dope_score_value(model, env)
            assessment_values = [model["name"], self.round_assessment_value(obj_funct_value), self.round_assessment_value(dope_score)]

            # Score with SOAP-PP.
            if self.multiple_chain_mode:
                soap_pp_assessor = soap_pp.Assessor()
                complete_mdl = complete_pdb(soap_env, os.path.join(self.modeling_dirpath, model['name']))
                atmsel = modeller.selection(complete_mdl)
                soap_pp_score = atmsel.assess(soap_pp_assessor)
                assessment_values.append(self.round_assessment_value(soap_pp_score))

            # Score with GA341 (only for single chain modeling).
            if self.modeling_mode == "homology_modeling" and not self.multiple_chain_mode:
                assessment_values.append(self.round_assessment_value(model["GA341 score"][0], digits_after_point=4))

            session_assessment_data.append(assessment_values)
            fmo.assessment_data = assessment_values


        # Prepares data to show a table with assessment values for each model.
        column_headers = ["Model Filename", "Objective Function", "DOPE score"]
        if self.multiple_chain_mode:
            column_headers.append("SOAP-PP")
        if self.modeling_mode == "homology_modeling" and not self.multiple_chain_mode:
            column_headers.append("GA341")

        assessment_table_args = {"column_headers": column_headers,
                                 "data_array": session_assessment_data,
                                 "title": "Assessment of Models of Modeling Session %s" % self.current_modeling_session.session_id,
                                 # "number_of_tabs": 4, "rowheader_width": 25,
                                 "width": 850, "height": 420,
                                 "sortable": True,
                                 }
        self.current_modeling_session.assessment_table_data = assessment_table_args


    def finish_modeling_session(self, successful=False, error_message=""):
        """
        Finishes the modeling session, both when models where sucessully built and when some error
        was encountered.
        """

        # Displays the models in PyMod main window, if some were built.
        self.pymod.main_window.gridder(update_menus=True, clear_selection=True, update_elements=successful, update_clusters=successful)

        if successful:

            # Colors the models and templates according to their DOPE values.
            if self.color_models_by_choice == "DOPE Score":
                for element in self.all_assessed_structures_list:
                    self.pymod.main_window.color_element_by_dope(element)

            else:
                # Color loops when performing loop refinement.
                if self.modeling_mode == "loop_refinement":
                    for mcl in self.modeling_clusters_list: # Cycle each target.
                        for model_element in mcl.model_elements_list: # Cycle each model.
                            for start_res, end_res, _chain in mcl.loops_list: # Cycle each loop.

                                # Builds a loop feature.
                                start_res = int(start_res)
                                end_res = int(end_res)
                                feature_name = "loop %s-%s" % (start_res, end_res)
                                # Convert the 'db_index' values into 'seq_index' values.
                                feature_start = model_element.get_residue_by_db_index(start_res).seq_index
                                feature_end = model_element.get_residue_by_db_index(end_res).seq_index
                                new_feature = Element_feature(id=None, name=feature_name,
                                                              start=feature_start, end=feature_end,
                                                              description=feature_name,
                                                              feature_type='sequence',
                                                              color=self.loop_default_color)
                                # Adds the loop feature to the model element.
                                model_element.add_feature(new_feature)
                            self.pymod.main_window.color_element_by_custom_scheme(model_element)

            # Increases modeling count.
            self.pymod.performed_modeling_count += 1

        # Reverts the stdout to the system one, and removes the modeling files.
        else:
            try:
                if os.path.isdir(self.modeling_dirpath):
                    shutil.rmtree(self.modeling_dirpath)
                self.pymod.main_window.show_error_message("Modeling Session Error",
                    "PyMod has encountered the following error while running MODELLER: %s" % error_message)
            except:
                self.pymod.main_window.show_error_message("Modeling Session Error",
                    "PyMod has encountered an unknown error in the modeling session: %s" % error_message)

        # Moves back to the current project directory.
        os.chdir(self.pymod.current_project_dirpath)


    def _get_modeller_outputs(self, a):
        if self.modeling_mode == "homology_modeling":
            return a.outputs
        elif self.modeling_mode == "loop_refinement":
            return a.loop.outputs

    def _get_model_filename(self):
        return "%s_%s_%s" % (self.modeller_target_name, self.get_model_number()+1, self.pymod.performed_modeling_count+1)


    ###############################################################################################
    # Prepares the modeling script.                                                               #
    ###############################################################################################

    def write_modeller_script(self, write_script=True, write_lib=True):
        """
        This methods will write two files: a script file for running MODELLER externally from PyMod
        (if the 'write_script' argument is set to 'True') and a module file storing the custom
        'automodel' class (if the 'write_lib' argument 'True'). In order to work externally from
        PyMod, both files must be present.
        """

        if not write_script and not write_lib:
            raise ValueError("Either one of 'write_script' or 'write_lib' must be set to 'True'.")

        if write_script:
            mod_script_fh = open(self.modeling_script_name, "w")
            # Environment.
            print("from %s import MyModel\n" % self.modeling_lib_name, file=mod_script_fh)
            print(self.mod_script_dict["environment"], file=mod_script_fh)
            print(self.mod_script_dict["hetres"], file=mod_script_fh)

        if write_lib:
            # Automodel derived class.
            mod_script_lib_fh = open(self.modeling_script_lib_name, "w")
            print("import modeller", file=mod_script_lib_fh)
            print("import modeller.automodel\n", file=mod_script_lib_fh)
            print(self.mod_script_dict["automodel"]["definition"], file=mod_script_lib_fh)
            if self.check_non_empty_automodel_dict():
                # Default patches.
                if self.mod_script_dict["automodel"]["default_patches"] != "":
                    print(self.mod_script_dict["automodel"]["default_patches"], file=mod_script_lib_fh)
                # Special patches.
                if self.check_non_empty_special_patches_dict():
                    print(self.mod_script_dict["automodel"]["special_patches"]["definition"], file=mod_script_lib_fh)
                    if self.mod_script_dict["automodel"]["special_patches"]["multichain"] != "":
                        print(self.mod_script_dict["automodel"]["special_patches"]["multichain"], file=mod_script_lib_fh)
                    if self.mod_script_dict["automodel"]["special_patches"]["disulfides"] != "":
                        print(self.mod_script_dict["automodel"]["special_patches"]["disulfides"], file=mod_script_lib_fh)
                # Special restraints.
                if self.mod_script_dict["automodel"]["special_restraints"] != "":
                    print(self.mod_script_dict["automodel"]["special_restraints"], file=mod_script_lib_fh)
            else:
                print("    pass\n", file=mod_script_lib_fh)
            mod_script_lib_fh.close()

        if write_script:
            # Model building options.
            print(self.mod_script_dict["automodel_init"], file=mod_script_fh)
            print(self.mod_script_dict["refinement"], file=mod_script_fh)
            print(self.mod_script_dict["models_indices"], file=mod_script_fh)
            print(self.mod_script_dict["make"], file=mod_script_fh)
            if not self.run_modeller_internally:
                 print(self.mod_script_dict["external_post_make"], file=mod_script_fh)
            mod_script_fh.close()


    def check_non_empty_automodel_dict(self):
        if self.mod_script_dict["automodel"]["default_patches"] != "":
            return True
        elif self.check_non_empty_special_patches_dict():
            return True
        elif self.mod_script_dict["automodel"]["special_restraints"] != "":
            return True
        else:
            return False

    def check_non_empty_special_patches_dict(self):
        if self.mod_script_dict["automodel"]["special_patches"]["multichain"] != "":
            return True
        elif self.mod_script_dict["automodel"]["special_patches"]["disulfides"] != "":
            return True
        else:
            return False


    ###########################################################################
    # Prepares the output of the modeling session.                            #
    ###########################################################################

    def get_model_number(self):
        return self.models_count

    def increase_model_number(self):
        self.models_count += 1


    def get_model_objective_function_value(self, model):
        """
        Gets the Objective function values of model. The model argument must be an item of the
        'outputs' list of the automodel class of MODELLER.
        """
        model_file = open(os.path.join(self.modeling_dirpath, model['name']), "r")
        obj_funct_value = float(model_file.readlines()[1][39:].replace(" ",""))
        model_file.close()
        return obj_funct_value

    def get_model_dope_score_value(self, model, env=None, use_automodel_assessment=True):
        """
        Gets the total DOPE score of a model.
        """
        if use_automodel_assessment:
            dope_score = model["DOPE score"]
        else:
            dope_score = compute_dope_of_structure_file(self.pymod,
                os.path.join(self.modeling_dirpath, model['name']),
                os.path.join(self.modeling_dirpath, "%s.profile" % model['name'][:-4]),
                env=env)
        return dope_score

    def get_model_color(self, chain_number, multiple_chain_mode):
        if not multiple_chain_mode:
            return self.single_chain_model_color
        else:
            return self.list_of_model_chains_colors[chain_number % len(self.list_of_model_chains_colors)]


###################################################################################################
# Modeling clusters and additional classes for the modeling protocol.                             #
###################################################################################################

class Modeling_cluster:
    """
    Class representing a "modeling cluster", that is a target-template(s) alignment
    used in a homology modeling protocol. When performing multiple chain homology
    modeling, there will be as many modeling clusters as the number of target
    chains.
    """

    def __init__(self, modeling_protocol, cluster):
        self.modeling_protocol = modeling_protocol

        self.cluster_element = cluster
        # This is actually the target sequence that is selected by the user.
        self.target = [c for c in self.cluster_element.get_children() if c.selected][0]
        # This is used to load the model into PyMOL when Modeller has done its job.
        self.target_name = self.get_target_name()

        # self.model_color=target.my_color
        self.aligned_elements_list = self.target.get_siblings(sequences_only=True)

        # Another for cycle to look for templates aligned to the target sequence. These will be
        # displayed in the modeling window.
        self.suitable_templates_list = [e for e in self.aligned_elements_list if e.is_suitable_template()]

        # This will contain a list objects from the Structure_frame class.
        self.structure_frame_list = []
        self.water_molecules_count = 0
        # Dictionary for additional options widgets.
        self.restraints_options_dict = {}

        # Chain ID of the template complex chain for this cluster.
        self.tc_chain_id = " "

        # Names of the template complex files.
        self.full_structure_files_dict = {}
        self.full_structure_fn_dict = {}

        # Index of the modeling cluster (corresponds to the numerical ID of the chain in a model,
        # the first chain has ID = 1, the second ID = 2, etc...).
        self.block_index = 0
        self.symmetry_restraints_id = None

        # List of the elements representing the models built in a session.
        self.model_elements_list = []


    def get_target_name(self):
        if not self.target.is_model():
            return pmos.clean_file_name(self.target.compact_header)
        else:
            return self.target.model_root_name # re.match("m\d+_(.*)_chain_[A-Za-z]","m1_test_chain_x").groups()

    def build_template_complex_chains_dict(self):
        """
        Generates modeling clusters dictionaries.
        They will be needed to check if a suitable "template complex" can be used. A cluster with
        the following templates:
            - 1HHO_Chain:A, 2DN2_Chain:A, 2DN2_Chain:B
        will generate a dictionary with the following structure:
            - {"1HHO.pdb":1, "2DN2.pdb":2}
        The keys are the original PDB files, and the values are the number of chains in the cluster
        which belong to that PDB structure.
        """
        for t in self.suitable_templates_list:
            template_original_file = t.get_structure_file(full_file=True)
            if template_original_file in list(self.full_structure_files_dict.keys()):
                self.full_structure_files_dict[template_original_file] += 1
            else:
                self.full_structure_files_dict.update({template_original_file: 1})
            self.full_structure_fn_dict[template_original_file] = t.structure.file_name_root + ".pdb"


    def set_options_from_gui(self):

        self.templates_list = []
        # Important dictionary that is going to contain informations about the templates.
        self.template_options_dict = {}

        template_count = 0
        for suitable_template, structure_frame in zip(self.suitable_templates_list, self.structure_frame_list):
            # Gets the values of each template checkbutton (it will be 0 if the structure was
            # not selected or 1 if it was selected): selects only structures that were selected
            # by the user to be used as templates.
            if structure_frame.is_selected():
                # Adds some information about the modeling options to the elements.
                template_options_dict = {
                    "id": template_count,
                    # For every selected structure takes the values of HETRES checkbutton.
                    "hetres_dict": structure_frame.get_template_hetres_dict(),
                    # Do the same with the water checkbutton.
                    "water_state": structure_frame.get_water_state_var(),
                    "structure_file": None, "sequence_file": None,
                    "modeller_name": None,
                    "template_complex": False,
                    "selected_template_complex":False}

                # Populate each modeling_cluster "template_list" with the elements selected by
                # the user from the "suitable_templates_list".
                self.add_new_template(pymod_element= suitable_template, template_options = template_options_dict)
                template_count += 1


    def add_new_template(self, pymod_element, template_options):
        self.templates_list.append(pymod_element)
        # This list will be used to inform Modeller about which are the "known" sequences.
        # It will contain the headers of the templates.
        self.template_options_dict.update({pymod_element: template_options})
        self.template_options_dict[pymod_element]["modeller_name"] = self.get_template_modeller_name(pymod_element)


    def get_template_modeller_name(self, pymod_element):
        return pymod_element.get_structure_file().replace(":","_")[:-4]
        # IF THE ORIGINAL PDB FILES ARE TO BE USED:
        #     - self.struct_list[a].structure.original_chain_pdb_file_name.replace(":","_")
        # In the original Pymod it was:
        #     - "1UBI_Chain_A" for non ce-aligned seqs
        #     - "1UBI_Chain_A_aligned.pdb" for aligned seqs
        # These codes must be the same in the .ali file and when assigning the "knowns".
        # If it the names don't contain the .pdb extension, Modeller will still find the
        # right files.


    def get_template_nameslist(self):
        ordered_keys = sorted(list(self.template_options_dict.keys()), key=lambda k:self.template_options_dict[k]["id"])
        return [self.template_options_dict[k]["modeller_name"] for k in ordered_keys]


    def use_water_in_cluster(self):
        return 1 in [self.template_options_dict[t]["water_state"] for t in self.templates_list]

    def has_target_with_multiple_cys(self):
        return self.target.my_sequence.count("C") >= 2


    def set_template_complex_chain(self, template):
        self.template_options_dict[template]["template_complex"] = True

    def set_template_complex_chain_to_use(self, template):
        self.template_options_dict[template]["selected_template_complex"] = True
        self.tc_chain_id = template.get_chain_id()
        # Temporarily sets the 'block_index' of a multiple chain model.
        self.block_index = template.get_chain_numeric_id()

    def is_template_complex_chain(self, template):
        """
        Check if a template chain is part of the 'template complex'.
        """
        return self.template_options_dict[template]["template_complex"]

    def get_template_complex_chain(self):
        for template in self.templates_list:
            if self.is_template_complex_chain(template):
                return template
        return None

    def get_single_chain_templates(self):
        return [t for t in self.templates_list if not self.is_template_complex_chain(t)]

    def get_modeling_elements(self):
        """
        Returns a list with the 'PyMod_elements' objects of the target and the templates of the
        modeling cluster.
        """
        modeling_elements = self.templates_list[:]
        modeling_elements.insert(0, self.target)
        return modeling_elements


    def prepare_single_chains_template_files(self):
        for template in self.templates_list:
            # if not self.is_template_complex_chain(template):
                self.prepare_template_files(template)

    def prepare_template_files(self, template):
        # Copy the templates structure files in the modeling directory.
        template_str_file = template.get_structure_file(basename_only=False)
        copied_template_str_file = os.path.basename(template_str_file)
        shutil.copy(template_str_file, os.path.join(self.modeling_protocol.modeling_dirpath, copied_template_str_file))
        self.template_options_dict[template]["structure_file"] = copied_template_str_file
        # Build a sequence file for the templates.
        # self.build_modeller_sequence_file(template)

    def build_modeller_sequence_file(self, template):
        # From point 17 of https://salilab.org/modeller/manual/node38.html.
        env = modeller.environ()
        modeller.log.none()
        env.io.hetatm = True
        if self.modeling_protocol.use_water_in_session:
            env.io.water = True
        structure_file_name = self.template_options_dict[template]["structure_file"]
        structure_file_code = os.path.splitext(structure_file_name)[0]
        mdl = modeller.model(env, file=os.path.join(self.modeling_protocol.modeling_dirpath, structure_file_name))
        aln = modeller.alignment(env)
        aln.append_model(mdl, align_codes=structure_file_code)
        output_sequence_file = structure_file_code+'_aln.chn'
        aln.write(file=os.path.join(self.modeling_protocol.modeling_dirpath, output_sequence_file))
        self.template_options_dict[template]["sequence_file"] = output_sequence_file

    def has_structures_with_disulfides(self):
        return True in [e.has_disulfides() for e in self.suitable_templates_list]

    def get_symmetry_chain_id(self):
        return self.block_index

    def get_symmetry_restraints_var(self):
        return int(self.restraints_options_dict["symmetry"].isChecked())


class Symmetry_restraints_groups_list:
    """
    This class will be used to store a list of Symmetry_restraint_group objects.
    """
    def __init__(self):
        self.list_of_groups = []

    def get_groups(self, min_number_of_sequences=0):
        """
        Returns a list of Symmetry_restraints_group objects that contain at least the number of
        sequences specified in the argument.
        """
        return [g for g in self.list_of_groups if len(g.list_of_clusters) >= min_number_of_sequences]

    def get_groups_to_use(self):
        return [g for g in self.get_groups(min_number_of_sequences=2) if g.use > 1]

    def add_group(self,symmetry_id):
        srg = Symmetry_restraints_group(symmetry_id)
        self.list_of_groups.append(srg)

    def get_group_by_id(self,symmetry_id):
        for g in self.list_of_groups:
            if g.id == symmetry_id:
                return g

    def get_symmetry_restraints_from_gui(self):
        for srg in self.get_groups(min_number_of_sequences=2):
            si = len([mc for mc in srg.list_of_clusters if mc.get_symmetry_restraints_var() == 1])
            srg.use = si

    def get_symmetry_restraints_list(self):
        """
        Builds a list of symmetry restraints to be used by MODELLER.
        """
        list_of_symmetry_restraints = []
        symmetry_restraints_groups_to_use = self.get_groups_to_use()
        if len(symmetry_restraints_groups_to_use) > 0:
            # Define group of chains on which symmetry restraints have to be applied.
            list_of_groups = []
            for srg in symmetry_restraints_groups_to_use:
                list_of_chains = []
                for mcl in srg.list_of_clusters:
                    if mcl.get_symmetry_restraints_var() == 1:
                        list_of_chains.append(mcl.get_symmetry_chain_id())
                list_of_groups.append(list_of_chains)

            for list_of_chains in list_of_groups:
                s = []
                for c in range(len(list_of_chains)):
                    i1 = list_of_chains[c]
                    i2 = None
                    if c < len(list_of_chains) - 1:
                        i2 = list_of_chains[c+1]

                    if i2 != None:
                        s.append([i1, i2])
                list_of_symmetry_restraints.append(s)
        return list_of_symmetry_restraints


class Symmetry_restraints_group:
    """
    When performing multichain modeling, this will be used to identify a "symmetry restraints
    group", a group of target sequences that share the exact same sequence. By keeping track of
    these groups, PyMod can let the user apply symmetry restraints to those chains when using
    MODELLER.
    """
    def __init__(self, symmetry_id):
        # The "id" is just the target sequence stripped of all indels.
        self.id = symmetry_id
        # This will contain a list of Modeling_cluster objects that contain a target sequence
        # with the same sequence as the "id".
        self.list_of_clusters = []
        # This will be set to True if the user decides to apply symmetry restraints to this group
        # of target sequences.
        self.use = False

    def add_cluster(self, modeling_cluster):
        """
        Adds a Modeling_cluster object to the group list_of_clusters.
        """
        self.list_of_clusters.append(modeling_cluster)


class Modeling_session_information:
    """
    Class for containing information on modeling sessions.
    """

    def __init__(self, session_id, automodel_output, modeling_directory_path):
        self.session_id = session_id
        self.modeling_directory_path = modeling_directory_path
        self.assessment_table_data = None
        self.dope_profile_data = []

        self.full_models = []
        for model in automodel_output:
            self.full_models.append(Full_model(os.path.join(self.modeling_directory_path, model['name'])))

    def add_template_dope_data(self, template_dope_data):
        """
        Stores the templates also profiles to each 'Full_model' object, so that the profile of the
        templates can be inspected by accessing the 'Models' menu.
        """
        for fmo in self.full_models: # All templates' data will be available in each model plot.
            fmo.dope_profile_data.append(template_dope_data)
        self.dope_profile_data.append(template_dope_data)


    def add_model_dope_data(self, model_dope_data, model_index):
        """
        Adds the DOPE profiles of the models.
        """
        self.full_models[model_index].dope_profile_data.append(model_dope_data)
        self.dope_profile_data.append(model_dope_data)


class Full_model:
    """
    Class for containing information on models built in a modeling session. Object of this class
    will be contained in the '.full_models' attribute of 'Modeling_session_information' class
    objects.
    """

    def __init__(self, original_file_path):
        self.original_file_path = original_file_path
        self.model_name = os.path.basename(self.original_file_path)[:-4]
        self.dope_profile_data = []
        self.assessment_data = None
