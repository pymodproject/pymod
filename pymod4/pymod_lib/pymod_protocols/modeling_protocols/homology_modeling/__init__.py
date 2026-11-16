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
from pymod_lib.pymod_protocols.modeling_protocols.modeller_residue_codes import (
    is_modelled_residue_code, is_nucleic_residue_name, split_nucleic_acids_into_chains,
    file_has_nucleic_residues)
from pymod_lib.pymod_protocols.modeling_protocols.homology_modeling._gui_qt import Modeling_window_qt
from pymod_lib.pymod_gui.shared_gui_components_qt import askyesno_qt


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
            ##### MODIFIED #####
            message = ("In order to perform 3D modeling in PyMod, MODELLER must be installed and configured correctly. %s" % modeller_error
                        + "\n" + "Open the 'PyMod Tools Installation Options' window to solve the issue (Help > Install PyMod Components)")
            ##### END #####
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
        
        # Checks that no nucleic acids have been selected.
        ##############MODIFIED on 18/02/2025####################################
        ##############warnings for Homology Modeling###############
        if any(e.polymer_type in ["dna", "rna"] for e in selected_sequences):
            self.pymod.main_window.show_error_message("Selection Error", "Cannot perform Homology Modeling on nucleic acids (RNA/DNA).")
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

        # If on Windows, warns before using parallel MODELLER jobs.
        if self.modeller_n_jobs > 1 and sys.platform == "win32":

            message = ("You decided to use multiple parallel MODELLER jobs on a"
                       " Windows machine. The use of parallel MODELLER jobs does"
                       " not work on certain Windows machines (the modeling task"
                       " could crash or freeze unexpectedly). Are you sure you"
                       " want to use parallel jobs? You may try to do it, and if"
                       " the task fails, you can attempt to run it again with"
                       " only one job."
                       " Note that even if your Windows machine can use parallel"
                       " MODELLER jobs, this will trigger a message from Windows"
                       " Defender (this is because PyMOL will spawn a separate"
                       " child process for each parallel MODELLER job)."
                       " However, you can safely suppress this Windows Defender"
                       " message, in this way you will be able to use all the"
                       " cores you want also on Windows.")
            choice = askyesno_qt("Use parallel MODELLER jobs?", message,
                                 parent=self.pymod.get_qt_parent())

            if not choice:
                raise ValueError("Please set in the 'Options' tab of the MODELLER"
                                 " window the 'N. of Parallel Jobs' option to 1"
                                 " in order to use only one MODELLER job.")


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

    #--------------------------------------------------------------------------
    # Creates a file with the target-template(s )alignment in the PIR format. -
    #--------------------------------------------------------------------------

    ###################MODIFICATION 25/05/2026-SARA#############################
    def build_pir_align_file(self):
        """
        Generates MODELLER PIR alignments by extracting the exact sequence 
        directly from the final filtered PDB files and mirroring it perfectly 
        into the Target sequence to prevent alignment errors.
        """
        import os
        import modeller

        # Open the PIR file for writing
        pir_align_file_handle = open(self.pir_file_path, "w")

        # 1. CHECK WATER STATE FROM THE GUI
        water_requested_by_gui = False
        for modeling_cluster in self.modeling_clusters_list:
            if hasattr(modeling_cluster, 'structure_frame_list'):
                for sf in modeling_cluster.structure_frame_list:
                    if hasattr(sf, 'get_water_state_var') and sf.get_water_state_var() == 1:
                        water_requested_by_gui = True
                        break

        global_water_count = 0
        if water_requested_by_gui:
            for modeling_cluster in self.modeling_clusters_list:
                for template in modeling_cluster.templates_list:
                    if hasattr(template, 'has_waters') and template.has_waters():
                        water_in_this_template = len(template.get_waters())
                        if water_in_this_template > global_water_count:
                            global_water_count = water_in_this_template

        # SINGLE SOURCE OF TRUTH FOR THE READING FLAGS.
        # The PIR alignment row for each template is just len(mdl.residues) for a
        # model MODELLER reads from the template file; at make() time
        # launch_modeller() reads the SAME file again and check_alignment()
        # compares the two counts. If the two reads use different env.io flags the
        # counts diverge and MODELLER aborts with:
        #   read_te_290E> Number of residues in the alignment and pdb files are different
        # (typically triggered by HETATM/nucleic acids or the water option). To make
        # a divergence impossible for ANY protein/DNA/RNA/ligand/water mix, we here
        # mirror EXACTLY the flag decision launch_modeller() makes (see the
        # "Sets heteroatoms and water options" block there): hetatm is forced on for
        # homology modeling / het session / any nucleic-acid template, and water
        # follows the session/GUI state.
        has_forced_nucleic_acids = any(
            getattr(t, 'chain_len_na_info_dict', {})
            for mc in self.modeling_clusters_list
            for t in getattr(mc, 'templates_list', []))
        force_read_hetatm = (self.modeling_mode == "homology_modeling"
                             or getattr(self, 'use_hetatm_in_session', False)
                             or has_forced_nucleic_acids)
        read_water = bool(getattr(self, 'use_water_in_session', False)
                          or water_requested_by_gui)

        # 2. BUILD PIR SEQUENCES DIRECTLY FROM PDB FILES
        for modeling_cluster in self.modeling_clusters_list:
            cluster_max_na_dots = 0
            
            for template in modeling_cluster.templates_list:
                env_check = modeller.environ()
                modeller.log.none()
                
                opts = modeling_cluster.template_options_dict.get(template, {})
                hetres_state = opts.get("hetres_state", 1)
                dict_info_na = getattr(template, 'chain_len_na_info_dict', {})

                # Configure MODELLER reading options to MIRROR launch_modeller()
                # exactly (see force_read_hetatm / read_water above), so this row's
                # length equals what check_alignment() will read from the same file.
                env_check.io.hetatm = True if force_read_hetatm else bool(dict_info_na and hetres_state != 3)
                env_check.io.water = read_water
                try:
                    env_check.io.allow_alternates = True
                except Exception:
                    pass

                structure_file_name = modeling_cluster.template_options_dict[template]["structure_file"]
                pdb_path = os.path.join(modeling_cluster.modeling_protocol.modeling_dirpath, structure_file_name)
                
                # Load the real model to get the exact sequence seen by MODELLER
                mdl_check = modeller.model(env_check, file=pdb_path)
                
                real_pdb_sequence = ""
                for residue in mdl_check.residues:
                    real_pdb_sequence += residue.code
                
                # Calculate additional inherited nucleic acid dots
                template_na_dots = 0
                selected_na_chains = opts.get("selected_na_chains", [])
                if dict_info_na:
                    for chain_id, (count, na_type_str) in dict_info_na.items():
                        if hetres_state == 1:
                            template_na_dots += count
                        elif hetres_state == 2:
                            if chain_id in selected_na_chains:
                                template_na_dots += count
                
                if template_na_dots > cluster_max_na_dots:
                    cluster_max_na_dots = template_na_dots

                # The nucleic acids are already in 'real_pdb_sequence': MODELLER reads
                # them from the file whatever io.hetatm says, and gives them their own
                # lowercase codes. So nothing is appended for them here -- appending
                # 'template_na_dots' dots would count every base twice. The count is
                # kept only to report how many bases this template carries.
                polymer_seq = real_pdb_sequence
                if global_water_count > 0 and not hasattr(template, 'has_waters'):
                    polymer_seq += "w" * global_water_count
                
                self.pir_sequences_dict[template] = {"pir_seq": polymer_seq, "modeling_cluster": modeling_cluster}

            # 3. CREATE TARGET WITH PERFECT MIRROR SYMMETRY
            # We mirror the template sequence directly into the target sequence.
            # An amino acid of the template is matched by a residue of the target,
            # while every position the target has nothing to say about -- a nucleic
            # acid, a ligand, a metal, a water -- is matched by a '.', which is how
            # a PIR alignment says "keep this residue from the template". Both rows
            # come out the same length either way, so this is the step that decides
            # whether the target is read in frame: see is_modelled_residue_code().
            target_poly_seq = ""

            # Use the reference template for this cluster
            if modeling_cluster.templates_list:
                reference_template = modeling_cluster.templates_list[0]
                template_seq_built = self.pir_sequences_dict[reference_template]["pir_seq"]

                target_amino_index = 0
                clean_target_seq = "".join([c for c in modeling_cluster.target.my_sequence if c not in (".", "w", "-")])

                for char in template_seq_built:
                    if not is_modelled_residue_code(char):
                        # A block residue: the target carries a dot in front of it,
                        # and none of its own residues are spent here.
                        target_poly_seq += "."
                    else:
                        # If it is an amino acid position, use the target sequence
                        if target_amino_index < len(clean_target_seq):
                            t_char = clean_target_seq[target_amino_index]
                            target_poly_seq += "A" if t_char == "X" else t_char
                            target_amino_index += 1
                        else:
                            # If target sequence ends but template has more amino acids, add a gap
                            target_poly_seq += "-"
            else:
                # Safety fallback if the cluster is empty
                for ali_pos in modeling_cluster.target.my_sequence:
                    if ali_pos in (".", "w"): target_poly_seq += "-"
                    else: target_poly_seq += "A" if ali_pos == "X" else ali_pos

            self.pir_sequences_dict[modeling_cluster.target] = {"pir_seq": target_poly_seq, "modeling_cluster": modeling_cluster}

        # 4. WRITE THE FORMATTED PIR FILE
        if self.multiple_chain_mode:
            print(">P1;%s" % self.template_complex_modeller_name, file=pir_align_file_handle)
            print("structure:%s:.:.:.:.::::" % self.template_complex_modeller_name, file=pir_align_file_handle)
            tc_pir_string = ""
            for template_complex_chain in self.get_template_complex_chains(sorted_by_id=True):
                tc_pir_string += self.pir_sequences_dict[template_complex_chain]["pir_seq"] + self.get_chain_separator()
            print(self.get_pir_formatted_sequence(tc_pir_string), file=pir_align_file_handle)

        for modeling_cluster in self.get_modeling_clusters_list(sorted_by_id=True):
            for template in modeling_cluster.get_single_chain_templates():
                template_code = modeling_cluster.template_options_dict[template]["modeller_name"]
                template_chain = template.get_chain_id()
                print(">P1;%s" % template_code, file=pir_align_file_handle)

                # Use a WHOLE-MODEL range ('.:.:.:.') for the structure entry rather
                # than an explicit first_res.num:last_res.num range.
                #
                # When nucleic acids / heteroatoms are merged into a protein chain,
                # merge_nucleic_acids_into_single_chain() renumbers the nucleic
                # residues but leaves other HETATM (e.g. Zn ions) on their original
                # numbers, so chain A ends up with NON-MONOTONIC and even COLLIDING
                # residue numbers (in 1AAY: protein 103-187, DNA 188-209, and Zn
                # 201-203 which both overlap the DNA and are not the max). In that
                # case residues[-1].num is not the largest number, so a numeric
                # "first:last" range fails to bracket every aligned residue and
                # MODELLER aborts with:
                #   read_te_290E> Number of residues in the alignment and pdb files
                #   are different
                # '.:.:.:.' tells MODELLER to read the entire model, which matches
                # the alignment row (also built from the whole model) exactly, for
                # any numbering. This mirrors what the template-complex and target
                # entries already do.
                print("structure:%s:.:.:.:.::::" % (template_code), file=pir_align_file_handle)
                
                sct_pir_string = ""
                for target_chain in self.get_targets_list(sorted_by_id=True):
                    if target_chain != modeling_cluster.target:
                        base_seq = self.pir_sequences_dict[target_chain]["pir_seq"]
                        gap_string = "".join([c if c in (".", "w") else "-" for c in base_seq])
                        sct_pir_string += gap_string + self.get_chain_separator()
                    else:
                        sct_pir_string += self.pir_sequences_dict[template]["pir_seq"] + self.get_chain_separator()
                print(self.get_pir_formatted_sequence(sct_pir_string), file=pir_align_file_handle)

        # Final block for target sequences
        print(">P1;%s" % self.modeller_target_name, file=pir_align_file_handle)
        print("sequence:%s:.:.:.:.::::" % self.modeller_target_name, file=pir_align_file_handle)
        targets_pir_string = ""
        for target in self.get_targets_list(sorted_by_id=True):
            targets_pir_string += self.pir_sequences_dict[target]["pir_seq"] + self.get_chain_separator()
        print(self.get_pir_formatted_sequence(targets_pir_string), file=pir_align_file_handle)

        pir_align_file_handle.close()
        
    ###################################END###############################

    """
    ################### MODIFICATION 19/05/2026 - SARA ##############################
   
    def build_pir_align_file(self):
        
        #This function creates alignments in a PIR format.
        #GUI-RESPECTING NA & WATER FIX: Forces 1:1 synchronization of BLK residues (.) 
        #and Water residues (w) between target and templates ONLY if selected in the GUI.
        
        import os

        # Starts to write the PIR sequence file needed by MODELLER for input.
        pir_align_file_handle = open(self.pir_file_path, "w")

        # 1. Determine how many positions to globally force for nucleic acids within the cluster
        global_na_dots = 0
        for modeling_cluster in self.modeling_clusters_list:
            for template in modeling_cluster.templates_list:
                dict_info_na = getattr(template, 'chain_len_na_info_dict', {})
                if dict_info_na:
                    for chain_id, (count, na_type_str) in dict_info_na.items():
                        global_na_dots += count

        # -------------------------------------------------------------------------------------
        # GUI STATE VERIFICATION: Verify if water inclusion is requested by the user
        # -------------------------------------------------------------------------------------
        water_requested_by_gui = False
        
        # Scan the modeling clusters and active frames to check the water checkbox state
        for modeling_cluster in self.modeling_clusters_list:
            if hasattr(modeling_cluster, 'structure_frame_list'):
                for sf in modeling_cluster.structure_frame_list:
                    if hasattr(sf, 'get_water_state_var') and sf.get_water_state_var() == 1:
                        water_requested_by_gui = True
                        break
            elif hasattr(self, 'structure_frame_list'):
                for sf in self.structure_frame_list:
                    if hasattr(sf, 'get_water_state_var') and sf.get_water_state_var() == 1:
                        water_requested_by_gui = True
                        break

        # 2. DIRECT WATER MOLECULE COUNT FROM INTERFACE DATA
        global_water_count = 0
        if water_requested_by_gui:
            for modeling_cluster in self.modeling_clusters_list:
                for template in modeling_cluster.templates_list:
                    # Query the native object components to count loaded water entities
                    if hasattr(template, 'has_waters') and template.has_waters():
                        water_in_this_template = len(template.get_waters())
                        if water_in_this_template > global_water_count:
                            global_water_count = water_in_this_template
            
            print(f">>> HARD-FIX PIR: Water ENABLED by GUI. Forced {global_water_count} 'w' residues in PIR file.")
        else:
            print(">>> HARD-FIX PIR: Water DISABLED by GUI. No 'w' residues added to PIR file.")

        for modeling_cluster in self.modeling_clusters_list:

            # 3. Construct the sanitized polymeric sequence for Template structures
            for template in modeling_cluster.templates_list:
                polymer_seq = ""
                for ali_pos in template.my_sequence:
                    if ali_pos == "." or ali_pos == "w":
                        polymer_seq += "-" 
                    else:
                        polymer_seq += ali_pos
                
                # Append required nucleic acid placeholders and water tokens dynamically
                if global_na_dots > 0:
                    polymer_seq += "." * global_na_dots
                if global_water_count > 0:
                    polymer_seq += "w" * global_water_count
                
                self.pir_sequences_dict[template] = {"pir_seq": polymer_seq, "modeling_cluster": modeling_cluster}

            # 4. Construct the sanitized polymeric sequence for the Target structure
            polymer_seq = ""
            for ali_pos in modeling_cluster.target.my_sequence:
                if ali_pos == "." or ali_pos == "w":
                    polymer_seq += "-"
                elif ali_pos == "X":
                    polymer_seq += "A"
                else:
                    polymer_seq += ali_pos
            
            # Append identical placeholders to ensure strict 1:1 sequence alignment size
            if global_na_dots > 0:
                polymer_seq += "." * global_na_dots
            if global_water_count > 0:
                polymer_seq += "w" * global_water_count
                
            self.pir_sequences_dict[modeling_cluster.target] = {"pir_seq": polymer_seq, "modeling_cluster": modeling_cluster}

        #--------------------------------
        # Write operations for the output PIR file
        #--------------------------------
        if self.multiple_chain_mode:
            print(">P1;%s" % self.template_complex_modeller_name, file=pir_align_file_handle)
            print("structure:%s:.:.:.:.::::" % self.template_complex_modeller_name, file=pir_align_file_handle)
            tc_pir_string = ""
            for template_complex_chain in self.get_template_complex_chains(sorted_by_id=True):
                tc_pir_string += self.pir_sequences_dict[template_complex_chain]["pir_seq"] + self.get_chain_separator()
            print(self.get_pir_formatted_sequence(tc_pir_string), file=pir_align_file_handle)

        for modeling_cluster in self.get_modeling_clusters_list(sorted_by_id=True):
            for template in modeling_cluster.get_single_chain_templates():
                template_code = modeling_cluster.template_options_dict[template]["modeller_name"]
                template_chain = template.get_chain_id()
                print(">P1;%s" % template_code, file=pir_align_file_handle)
                print("structure:%s:.:%s:.:%s::::" % (template_code, template_chain, template_chain), file=pir_align_file_handle)
                
                sct_pir_string = ""
                for target_chain in self.get_targets_list(sorted_by_id=True):
                    if target_chain != modeling_cluster.target:
                        base_seq = self.pir_sequences_dict[target_chain]["pir_seq"]
                        gap_string = ""
                        for char in base_seq:
                            if char == "." or char == "w":
                                gap_string += char
                            else:
                                gap_string += "-"
                        sct_pir_string += gap_string + self.get_chain_separator()
                    else:
                        sct_pir_string += self.pir_sequences_dict[template]["pir_seq"] + self.get_chain_separator()
                print(self.get_pir_formatted_sequence(sct_pir_string), file=pir_align_file_handle)

        print(">P1;%s" % self.modeller_target_name, file=pir_align_file_handle)
        print("sequence:%s:.:.:.:.::::" % self.modeller_target_name, file=pir_align_file_handle)
        targets_pir_string = ""
        for target in self.get_targets_list(sorted_by_id=True):
            targets_pir_string += self.pir_sequences_dict[target]["pir_seq"] + self.get_chain_separator()
        print(self.get_pir_formatted_sequence(targets_pir_string), file=pir_align_file_handle)

        pir_align_file_handle.close()
    """
    def _get_hetres_seq(self, hetres_list):

        """
        Generates the 1-letter code string for non-polymer hetero-residues (ligands/water).
        
        NUCLEIC ACIDS BEHAVIOR NOTE:
        To prevent structural redundancy, this function is bypassed or strictly filtered 
        when processing nucleic acid chains. Because RNA/DNA nucleotides are explicitly 
        handled as custom structural blocks appended at the end of the alignment matrix, 
        this method focuses solely on standard chemical ligands and water molecules ('w') 
        to prevent double-counting sequence characters.
        """

        hetres_seq = ""
        if not hetres_list:
            return hetres_seq
        for hetres in hetres_list:
            if not hetres.is_water():
                hetres_seq += "."
            else:
                hetres_seq += "w"
        return hetres_seq

    #########################END#############################

    """OLD
    def build_pir_align_file(self):
        
        #This function creates alignments in a PIR format: this is entirely rewritten from the
        #original PyMod version. This method should write the sequences as seen by MODELLER.
        

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
                            if modeling_cluster.template_options_dict[template]["hetres_dict"].get(res, False):
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
                    if modeling_cluster.template_options_dict[template_i]["hetres_dict"].get(hetres, False) and self.use_hetatm_in_session:
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

    """

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
            import traceback
            self.finish_modeling_session(successful=False, error_message=traceback.format_exc())

    ######################MODIFICATION on 17/05/2026-SARA########################################
    def launch_modeller(self, in_thread=False):
        """
        Sets the MODELLER environment and script, then launches 3D model building.
        Modded Version: Forces env.io.hetatm, env.io.allow_alternates and contact shells 
        when nucleic acids are present to bypass SequenceMismatchError.

    
        Initializes the MODELLER runtime environment and compiles the execution script.
        
        NUCLEIC ACIDS ENVIRONMENT PATCH:
        When forced nucleic acid structures are detected in the active cluster, the default 
        modeling restrictions are programmatically overridden to accommodate BLK residues:
        1. Forces 'env.io.hetatm = True' to allow the parsing of non-protein atoms.
        2. Forces 'env.io.allow_alternates = True' to ease structural reading tolerance, 
           bypassing 'SequenceMismatchError' caused by spatial coordinate offsets between 
           the edited text sequence and the raw PDB records.
        3. Configures specialized environmental restraints ('env.edat.contact_shell = 4.0' 
           and 'env.edat.dynamic_sphere = True') under default optimization levels. This ensures 
           proper electrostatic and van der Waals contact calculations around the heavy 
           arbitrary BLK blocks representing the nucleic acid backbones.
        """
        import os
        import sys
        import pickle
        import importlib
        from modeller import parallel

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
        # --- FORZATURA STRUTTURALE DI PRIMA FASE ---
        # Verifichiamo se almeno un template nel cluster contiene acidi nucleici mappati
        has_forced_nucleic_acids = False
        for mc in self.modeling_clusters_list:
            # Loop_modeling_cluster (loop refinement) has no templates_list,
            # since loop modeling refines an existing model rather than using
            # templates. Treat the missing attribute as an empty list.
            for template in getattr(mc, 'templates_list', []):
                if getattr(template, 'chain_len_na_info_dict', {}):
                    has_forced_nucleic_acids = True
                    break

        if has_forced_nucleic_acids:
            print(">>> Modeller Engine: Nucleic Acids detected! Forcing env.io.hetatm and allow_alternates = True")

       ##################### MODIFICATION On 19/05/2026 - SARA ######################
        # -------------------------------------------------------------------------------------
        # WATER & HETATM PATCH SECTION (Synchronous Enforcement with the GUI State)
        # -------------------------------------------------------------------------------------
        # Force a direct check on the GUI clusters to see if the user ticked the water checkbox.
        # This bypasses any internal PyMod communication delay or desynchronization issues.
        water_requested_by_gui = False
        if hasattr(self, 'modeling_clusters_list'):
            for modeling_cluster in self.modeling_clusters_list:
                if hasattr(modeling_cluster, 'structure_frame_list'):
                    for sf in modeling_cluster.structure_frame_list:
                        if hasattr(sf, 'get_water_state_var') and sf.get_water_state_var() == 1:
                            water_requested_by_gui = True
                            break
        
        # If the GUI has the water inclusion option active, explicitly force the session flag to True.
        # This guarantees that the core MODELLER environment will initialize its water-parsing routine.
        if water_requested_by_gui:
            self.use_water_in_session = True
            print(">>> LAUNCH MODELLER PATCH: Detected active water option from GUI. Forced 'use_water_in_session = True'")

        # Evaluate if the run requires hetero-atoms (HETATM records) or water molecules.
        # In homology modeling protocols, the 'hetatm' flag is strictly enforced to preserve structural ligands.
        if self.modeling_mode == "homology_modeling" or self.use_hetatm_in_session or has_forced_nucleic_acids:
            
            # --- Internal Execution Mode --------------------------------------------
            if self.run_modeller_internally:
                # Enable general heteroresidue and non-protein parsing (ligands, ions)
                env.io.hetatm = True
                # PATCH: Enable alternate residue matching to handle hybrid setups and custom nucleic acid definitions,
                # effectively preventing fatal SequenceMismatchError exceptions during sequence-structure mapping.
                env.io.allow_alternates = True  
            
            # --- External Script Generation Mode ------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["hetres"] += "env.io.hetatm = True\n"
                # Mirror the patch configuration in the standalone Python executable script
                self.mod_script_dict["hetres"] += "env.io.allow_alternates = True\n"  
            
            # ------------------------------------------------------------------------
            # Explicit Water Molecule Handling
            # ------------------------------------------------------------------------
            # Process crystallographic water molecules ('HOH' / 'w') only if confirmed by the session environment.
            if self.use_water_in_session:
                print(">>> LAUNCH MODELLER: Activating 'env.io.water = True' for the MODELLER engine.")
                
                # --- Internal Execution ---------------------------------------------
                if self.run_modeller_internally:
                    # Instruct the internal MODELLER io handler to map and transfer 'w' residues into the target structure
                    env.io.water = True
                
                # --- External Script ------------------------------------------------
                if self.write_modeller_script_option:
                    # Append the environment flag to the generated external execution pipeline
                    self.mod_script_dict["hetres"] += "env.io.water = True\n"
                    
        ################################# END ##########################################

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
        # External --------------------------------------------------
        if self.write_modeller_script_option:
            self.mod_script_dict["automodel"]["definition"] += "class MyModel(modeller.automodel.%s):\n" % session_automodel_class_name
        #------------------------------------------------------------

        #*******************************
        # Template-derived disulfides. *
        #*******************************
        if self.modeling_mode == "homology_modeling":
            if self.check_targets_with_cys() and self.check_structures_with_disulfides():
                if self.use_template_dsb:
                    pass
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
            if self.multiple_chain_mode:
                if self.modeling_mode == "homology_modeling":
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        # Renumber the residues of the chains of multichain models starting from 1.\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        count_dictionary = dict((chain.name, 1) for chain in self.chains)\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        for chain in self.chains:\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "            for residue in chain.residues:\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "                residue.num = '%d' % (count_dictionary[chain.name])\n"
                    self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "                count_dictionary[chain.name] += 1\n"

            if self.check_targets_with_cys():
                if self.use_user_defined_dsb:
                    for (mci, mc) in enumerate(self.modeling_clusters_list):
                        for dsb in self.all_user_defined_dsb[mci]:
                            cys1 = dsb[0][4:]
                            cys2 = dsb[1][4:]
                            if self.modeling_mode == "homology_modeling":
                                if self.multiple_chain_mode:
                                    chain = mc.block_index
                                    self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.chains[%s].residues['%s'], self.chains[%s].residues['%s']))\n" % (chain, cys1, chain, cys2)
                                else:
                                    self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.residues['%s'], self.residues['%s']))\n" % (cys1, cys2)
                            elif self.modeling_mode == "loop_refinement":
                                chain = mc.target.get_chain_id()
                                self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.chains['%s'].residues['%s'], self.chains['%s'].residues['%s']))\n" % (chain, cys1, chain, cys2)
                if self.use_auto_dsb:
                    self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch_ss()\n"
        #------------------------------------------------------------

        #**************************************************************************
        # Apply symmetry restraints to target chains that have the same sequence. *
        #**************************************************************************
        if self.multiple_chain_mode and len(self.list_of_symmetry_restraints) > 0:
            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["automodel"]["special_restraints"] += "    def special_restraints(self, aln):\n"
                for si, symmetry_restraints_group in enumerate(self.list_of_symmetry_restraints):
                    self.mod_script_dict["automodel"]["special_restraints"] += "        # Symmetry restraints group n. %d.\n" % (si+1)
                    for s in symmetry_restraints_group:
                        if self.modeling_mode == "homology_modeling":
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s1 = modeller.selection(self.chains[%s]).only_atom_types('CA')\n" % (s[0])
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s2 = modeller.selection(self.chains[%s]).only_atom_types('CA')\n" % (s[1])
                        elif self.modeling_mode == "loop_refinement":
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s1 = modeller.selection(self.chains['%s']).only_atom_types('CA')\n" % (s[0])
                            self.mod_script_dict["automodel"]["special_restraints"] += "        s2 = modeller.selection(self.chains['%s']).only_atom_types('CA')\n" % (s[1])
                        self.mod_script_dict["automodel"]["special_restraints"] += "        self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "    def user_after_single_model(self):\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "        self.restraints.symmetry.report(1.0)\n\n"
            #--------------------------------------------------------

        #******************
        # Loop selection. *
        #******************
        if self.modeling_mode == "loop_refinement":
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
            self.write_modeller_script(write_script=False, write_lib=True)
            sys.path.append(self.modeling_dirpath)
            _session_mod_lib = importlib.import_module(self.modeling_lib_name)
            _session_mod_lib = importlib.reload(_session_mod_lib)
            MyModel = _session_mod_lib.MyModel
            sys.path.remove(self.modeling_dirpath)
            if self.modeling_mode == "homology_modeling":
                if not self.multiple_chain_mode:
                    assess_methods = (modeller.automodel.assess.GA341, modeller.automodel.assess.DOPE)
                else:
                    assess_methods = (modeller.automodel.assess.DOPE, )
                a = MyModel(env,
                            alnfile=self.pir_file_name,
                            knowns=tuple([str(tmpn) for tmpn in self.all_templates_namelist]),
                            sequence=str(self.modeller_target_name),
                            assess_methods=assess_methods)
            elif self.modeling_mode == "loop_refinement":
                loop_assess_methods = (modeller.automodel.assess.DOPE, )
                a = MyModel(env,
                            sequence=str(self.modeller_target_name),
                            inimodel=self.loop_refinement_starting_model,
                            loop_assess_methods=loop_assess_methods)
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
        # Sets the level of refinement and optimization. *
        #************************************************
        if self.modeling_mode == "homology_modeling":
            if self.optimization_level == "Low":
                if self.run_modeller_internally:
                    a.library_schedule = modeller.automodel.autosched.very_fast
                    a.md_level = modeller.automodel.refine.very_fast
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.very_fast\n"
                    self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.very_fast\n"
            
            elif self.optimization_level == "Default":
                if has_forced_nucleic_acids:
                    if self.run_modeller_internally:
                        a.env.edat.contact_shell = 4.0
                        a.env.edat.dynamic_sphere = True
                    if self.write_modeller_script_option:
                        self.mod_script_dict["refinement"] += "\n# Support for Nucleic Acids BLK restraints\na.env.edat.contact_shell = 4.0\na.env.edat.dynamic_sphere = True\n"
                else:
                    pass
            
            elif self.optimization_level == "Mid":
                if self.run_modeller_internally:
                    a.library_schedule = modeller.automodel.autosched.normal
                    a.max_var_iterations = 300
                    a.md_level = modeller.automodel.refine.slow
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.fast\n"
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = 300\n"
                    self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.slow\n"
            
            elif self.optimization_level == "High":
                if self.run_modeller_internally:
                    a.library_schedule = modeller.automodel.autosched.slow
                    a.max_var_iterations = 300
                    a.md_level = modeller.automodel.refine.slow
                    a.repeat_optimization = 2
                    a.max_molpdf = 1e8
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.slow\n"
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = 300\n"
                    self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.slow\n"
                    self.mod_script_dict["refinement"] += "a.repeat_optimization = 2\n"
                    self.mod_script_dict["refinement"] += "a.max_molpdf = 1e8\n"
            
            elif self.optimization_level == "Approximate":
                if self.run_modeller_internally:
                    a.max_ca_ca_distance = 10.0
                    a.max_n_o_distance = 6.0
                    a.max_sc_mc_distance = 5.0
                    a.max_sc_sc_distance = 4.5
                    a.max_var_iterations = 50
                    a.library_schedule = modeller.automodel.autosched.fastest
                    a.md_level = None
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.max_ca_ca_distance = 10.0\n"
                    self.mod_script_dict["refinement"] += "a.max_n_o_distance = 6.0\n"
                    self.mod_script_dict["refinement"] += "a.max_sc_mc_distance = 5.0\n"
                    self.mod_script_dict["refinement"] += "a.max_sc_sc_distance = 4.5\n"
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = 50\n"
                    self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.fastest\n"
                    self.mod_script_dict["refinement"] += "a.md_level = None\n"
            
            elif self.optimization_level == "Custom":
                if self.run_modeller_internally:
                    a.max_var_iterations = self.custom_max_cg_iterations
                    a.library_schedule = self.get_modeller_schedule_obj("cg", "internal", self.custom_vtfm_schedule)
                    a.md_level = self.get_modeller_schedule_obj("md", "internal", self.custom_md_schedule)
                    a.repeat_optimization = self.custom_repeat_optimization
                    a.max_molpdf = 10**self.custom_max_obj_func_value
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.max_var_iterations = %s\n" % self.custom_max_cg_iterations
                    self.mod_script_dict["refinement"] += "a.library_schedule = %s\n" % self.get_modeller_schedule_obj("cg", "external", self.custom_vtfm_schedule)
                    self.mod_script_dict["refinement"] += "a.md_level = %s\n" % self.get_modeller_schedule_obj("md", "external", self.custom_md_schedule)
                    self.mod_script_dict["refinement"] += "a.repeat_optimization = %s\n" % self.custom_repeat_optimization
                    self.mod_script_dict["refinement"] += "a.max_molpdf = 10**%s\n" % self.custom_max_obj_func_value

        elif self.modeling_mode == "loop_refinement":
            if self.optimization_level == "Low":
                if self.run_modeller_internally:
                    a.loop.md_level = modeller.automodel.refine.very_fast
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = modeller.automodel.refine.very_fast\n"
            elif self.optimization_level == "Mid":
                if self.run_modeller_internally:
                    a.loop.md_level = modeller.automodel.refine.fast
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = modeller.automodel.refine.fast\n"
            elif self.optimization_level == "Default":
                pass
            elif self.optimization_level == "High":
                if self.run_modeller_internally:
                    a.loop.md_level = modeller.automodel.refine.very_slow
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = modeller.automodel.refine.very_slow\n"
            elif self.optimization_level == "Approximate":
                if self.run_modeller_internally:
                    a.loop.md_level = None
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.md_level = None\n"
            elif self.optimization_level == "Custom":
                if self.run_modeller_internally:
                    a.loop.max_var_iterations = self.custom_max_cg_iterations
                    a.loop.md_level = self.get_modeller_schedule_obj("md", "internal", self.custom_md_schedule)
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += "a.loop.max_var_iterations = %s\n" % self.custom_max_cg_iterations
                    self.mod_script_dict["refinement"] += "a.loop.md_level = %s\n" % self.get_modeller_schedule_obj("md", "external", self.custom_md_schedule)

        #************************************
        # Customize the objective function. *
        #************************************
        if self.use_custom_obj_func in ("customize", "altmod"):
            if self.modeling_mode == "homology_modeling":
                w_dope = self.custom_w_dope if self.use_dope_in_obj_func else 0
                if self.run_modeller_internally:
                    a.env.schedule_scale = modeller.physical.values(default=1.0,
                                                                    nonbond_spline=w_dope,
                                                                    ca_distance=self.custom_w_hddr,
                                                                    n_o_distance=self.custom_w_hddr,
                                                                    sd_mn_distance=self.custom_w_hddr,
                                                                    sd_sd_distance=self.custom_w_hddr,
                                                                    phi_dihedral=self.custom_w_hdar,
                                                                    psi_dihedral=self.custom_w_hdar,
                                                                    chi1_dihedral=self.custom_w_hdar,
                                                                    chi2_dihedral=self.custom_w_hdar,
                                                                    chi3_dihedral=self.custom_w_hdar,
                                                                    chi4_dihedral=self.custom_w_hdar,
                                                                    phi_psi_dihedral=self.custom_w_hdar,
                                                                    bond=self.custom_w_charmm22,
                                                                    angle=self.custom_w_charmm22,
                                                                    dihedral=self.custom_w_charmm22,
                                                                    improper=self.custom_w_charmm22,
                                                                    soft_sphere=self.custom_w_soft_sphere)
                    edat = a.env.edat
                    edat.contact_shell = self.custom_nb_cutoff
                    if self.use_dope_in_obj_func:
                        edat.dynamic_modeller = True
                        gprsr = modeller.group_restraints(a.env, classes='$(LIB)/atmcls-mf.lib', parameters='$(LIB)/dist-mf.lib')
                        a.group_restraints = gprsr
                
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += """
# Sets the weights of the objective function.
a.env.schedule_scale = modeller.physical.values(default=1.0,
                                                nonbond_spline=%s,
                                                ca_distance=%s,
                                                n_o_distance=%s,
                                                sd_mn_distance=%s,
                                                sd_sd_distance=%s,
                                                phi_dihedral=%s,
                                                psi_dihedral=%s,
                                                chi1_dihedral=%s,
                                                chi2_dihedral=%s,
                                                chi3_dihedral=%s,
                                                chi4_dihedral=%s,
                                                phi_psi_dihedral=%s,
                                                bond=%s,
                                                angle=%s,
                                                dihedral=%s,
                                                improper=%s,
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

            elif self.modeling_mode == "loop_refinement":
                if self.run_modeller_internally:
                    a.loop.env.schedule_scale = modeller.physical.values(default=1.0,
                                                                         nonbond_spline=self.custom_w_stat_pot_loop,
                                                                         phi_dihedral=self.custom_w_hdar,
                                                                         psi_dihedral=self.custom_w_hdar,
                                                                         chi1_dihedral=self.custom_w_hdar,
                                                                         chi2_dihedral=self.custom_w_hdar,
                                                                         chi3_dihedral=self.custom_w_hdar,
                                                                         chi4_dihedral=self.custom_w_hdar,
                                                                         phi_psi_dihedral=self.custom_w_hdar,
                                                                         bond=self.custom_w_charmm22,
                                                                         angle=self.custom_w_charmm22,
                                                                         dihedral=self.custom_w_charmm22,
                                                                         improper=self.custom_w_charmm22,
                                                                         soft_sphere=self.custom_w_soft_sphere,
                                                                         lennard_jones=self.custom_w_lj,
                                                                         gbsa=self.custom_w_gbsa,
                                                                         )
                    edat = a.loop.env.edat
                    edat.contact_shell = self.custom_nb_cutoff
                if self.write_modeller_script_option:
                    self.mod_script_dict["refinement"] += """
# Sets the weights of the objective function.
a.loop.env.schedule_scale = modeller.physical.values(default=1.0,
                                                     nonbond_spline=%s,
                                                     phi_dihedral=%s,
                                                     psi_dihedral=%s,
                                                     chi1_dihedral=%s,
                                                     chi2_dihedral=%s,
                                                     chi3_dihedral=%s,
                                                     chi4_dihedral=%s,
                                                     phi_psi_dihedral=%s,
                                                     bond=%s,
                                                     angle=%s,
                                                     dihedral=%s,
                                                     improper=%s,
                                                     soft_sphere=%s,
                                                     lennard_jones=%s,
                                                     gbsa=%s,
                                                     )
# Non-bonded interactions cutoff.
edat = a.loop.env.edat
edat.contact_shell = %s
""" % (self.custom_w_stat_pot_loop, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_hdar, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_charmm22, self.custom_w_soft_sphere, self.custom_w_lj, self.custom_w_gbsa, self.custom_nb_cutoff)

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
            self.write_modeller_script(write_script=True, write_lib=False)
        if not self.run_modeller_internally:
            raise Exception("Not implemented anymore.")
        #------------------------------------------------------------
        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            if self.modeling_mode == "homology_modeling":
                a.starting_model = self.starting_model_number
                a.ending_model = self.ending_model_number
            elif self.modeling_mode == "loop_refinement":
                a.loop.starting_model = self.starting_model_number
                a.loop.ending_model = self.ending_model_number
            if self.modeller_n_jobs > 1:
                a.use_parallel_job(j)
            
            # This is the method that launches the model building phase.
            a.make()
        #------------------------------------------------------------
        if in_thread:
            with open(self.modeller_vars_filename, "wb") as v_fh:
                pickle.dump(a, v_fh)
        else:
            return env, a
    
    ########################END#############################################


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
        print(f"List of dict modeller output: {automodel_output}")
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
                # Give the nucleic acids their own chains back before the element is
                # loaded, so that PyMOL shows one filament per chain.
                self.restore_nucleic_acid_chains(model_element, modeling_cluster)
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


    def restore_nucleic_acid_chains(self, model_element, modeling_cluster):
        """
        Puts the nucleic acids of a finished model back in chains of their own.

        A template is given to MODELLER with its DNA or RNA merged into the protein
        chain, because that is the only way MODELLER reads it, and the model comes
        back the same way: one chain holding the protein and every base. That is not
        what anyone expects to see in PyMOL, where a double helix is two chains, one
        per filament. So the file of the model is rewritten before it is loaded,
        handing the bases back to the chains they came from.

        Does nothing at all when there were no nucleic acids, or when the model does
        not have the number of bases the templates say it should.
        """

        # The filaments as they were in the structure the user opened, kept only if
        # they were asked for: the model carries the chains the template carried, and
        # the template was filtered by what was chosen in the modeling window.
        chain_lengths = []
        for template in getattr(modeling_cluster, 'templates_list', []):
            info = getattr(template, 'chain_len_na_info_dict', {})
            if not info:
                continue
            options = {}
            try:
                options = modeling_cluster.template_options_dict.get(template, {})
            except Exception:
                pass
            hetres_state = options.get("hetres_state", 1)
            if hetres_state == 3:
                # Nothing was kept, so there is nothing to give back.
                return
            # What the template actually contributed, when it is known: in a model of
            # several chains each filament comes from one template only, so what was
            # ticked here is not always what is in this chain of the model.
            contributed = options.get("contributed_na_chains")
            if contributed is None:
                selected = options.get("selected_na_chains", [])
                contributed = [chain_id for chain_id in info
                               if hetres_state != 2 or chain_id in selected]
            # The order of this dictionary is the order the chains had in the file,
            # which is the order the merge kept.
            chain_lengths = [(chain_id, count)
                             for chain_id, (count, na_type) in info.items()
                             if chain_id in contributed]
            break
        if not chain_lengths:
            return

        try:
            model_file_path = model_element.get_structure_file(basename_only=False)
            with open(model_file_path) as model_file:
                lines = model_file.readlines()

            new_lines = split_nucleic_acids_into_chains(lines, chain_lengths)
            if new_lines is None:
                print(">>> The model does not carry the %s nucleic acid residues of"
                      " its templates, so its chains were left as MODELLER wrote"
                      " them." % sum(count for _, count in chain_lengths))
                return

            with open(model_file_path, "w") as model_file:
                model_file.writelines(new_lines)
            print(">>> Nucleic acids of the model put back in their own chains: %s"
                  % ", ".join("%s (%s residues)" % (chain_id, count)
                              for chain_id, count in chain_lengths))

        except Exception as e:
            # Cosmetic work on a model which is otherwise finished: worth a message,
            # not worth losing the model over.
            print(">>> Could not give the nucleic acids of the model their own"
                  " chains: %s" % e)



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

    ###############MODIFICATION 25/05/2026-SARA######################
    def set_options_from_gui(self):

        self.templates_list = []
        # Important dictionary that is going to contain information about the templates.
        self.template_options_dict = {}

        template_count = 0
        for suitable_template, structure_frame in zip(self.suitable_templates_list, self.structure_frame_list):
            # Gets the values of each template checkbutton (it will be 0 if the structure was
            # not selected or 1 if it was selected): selects only structures that were selected
            # by the user to be used as templates.
            if structure_frame.is_selected():
                
                # --- NEW LOGIC FOR NUCLEIC ACID CHAIN SELECTION ---
                # Get the state of the hetero residues selection from your modified GUI
                hetres_state = structure_frame.get_hetres_options_var()
                
                # Retrieve the selected NA chains if the user clicked "Select nucleic acid ligand"
                selected_na_chains = []
                if hetres_state == 2:
                    selected_na_chains = structure_frame.get_selected_na_chain_ids()
                # --------------------------------------------------

                # Adds some information about the modeling options to the elements.
                template_options_dict = {
                    "id": template_count,
                    # Save the high-level state (1=All, 2=Selective, 3=None)
                    "hetres_state": hetres_state,
                    # For every selected structure takes the values of HETRES checkbutton.
                    "hetres_dict": structure_frame.get_template_hetres_dict(),
                    # Store the custom nucleic acid chain IDs selected by the user
                    "selected_na_chains": selected_na_chains,
                    # Do the same with the water checkbutton.
                    "water_state": structure_frame.get_water_state_var(),
                    "structure_file": None, 
                    "sequence_file": None,
                    "modeller_name": None,
                    "template_complex": False,
                    "selected_template_complex": False
                }

                # Populate each modeling_cluster "template_list" with the elements selected by
                # the user from the "suitable_templates_list".
                self.add_new_template(pymod_element=suitable_template, template_options=template_options_dict)
                template_count += 1
        
        ###########################END#########################
    """
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
    """

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
        # Which template of this chain the nucleic acids come from. Every template file of
        # a structure with nucleic acids carries all of them, merged into its own protein
        # chain, so several templates of one chain would each bring a copy of the filaments,
        # and every copy in the frame of its own template. Taking one filament from one
        # template and another from another is no better: they arrive in unrelated frames.
        # So the nucleic acids of a chain come from a single template, the first which asks
        # for any, and the others give none.
        #
        # Per chain of the model, deliberately. A model of several chains gives each chain
        # its own template and its own filaments, and that is right: nothing here reaches
        # across to another chain.
        self.na_donor_template = None
        for template in self.templates_list:
            # if not self.is_template_complex_chain(template):
                self.prepare_template_files(template)
    
    #####################MODIFICATION 25/05/2026 - SARA ########################
    def prepare_template_files(self, template):
        import os
        import shutil
        
        template_str_file = template.get_structure_file(basename_only=False)
        copied_template_str_file = os.path.basename(template_str_file)
        destination_path = os.path.join(self.modeling_protocol.modeling_dirpath, copied_template_str_file)

        opts = self.template_options_dict.get(template, {})
        hetres_state = opts.get("hetres_state", 1)
        selected_na_chains = opts.get("selected_na_chains", [])
        dict_info_na = getattr(template, 'chain_len_na_info_dict', {})

        # Which file the copy is made from. The ordinary chain file of a template
        # holds its protein alone: when the structure carries nucleic acids they were
        # put in a second file, merged into the protein chain, because a chain of
        # bases of its own is not something MODELLER reads. Without that file the
        # bases never reach the model at all, which is what used to happen. The name
        # of the copy does not change, so the alignment and everything downstream
        # still call this template what they called it before.
        source_str_file = template_str_file
        if dict_info_na:
            merged_str_file = getattr(template, 'structure_file_path', None)
            if (merged_str_file and os.path.isfile(merged_str_file)
                    and not file_has_nucleic_residues(template_str_file)
                    and file_has_nucleic_residues(merged_str_file)):
                source_str_file = merged_str_file
                print("[TEMPLATE] '%s': taking the nucleic acids from '%s'."
                      % (copied_template_str_file, os.path.basename(merged_str_file)))
        
        # Which filaments this template is the one to contribute. What the user asked for
        # here, less whatever an earlier template of this model has already put in.
        contributes_na_chains = set()
        if dict_info_na and hetres_state != 3:
            asked_for = (set(selected_na_chains) if hetres_state == 2
                         else set(dict_info_na))
            donor = getattr(self, 'na_donor_template', None)
            if donor is None and asked_for:
                self.na_donor_template = template
                donor = template
            if donor is template:
                contributes_na_chains = asked_for
                taken_by_another = set()
            else:
                # Another template of this same chain is the one giving the nucleic acids.
                contributes_na_chains = set()
                taken_by_another = asked_for
            # What this template really gave. The restoration measures the model against
            # this, not against what was ticked: a filament which another template brought
            # is not in this chain of the model, and counting it would make the numbers
            # disagree and leave the bases unsplit.
            self.template_options_dict.setdefault(template, {})
            self.template_options_dict[template]["contributed_na_chains"] = sorted(
                contributes_na_chains)
            if taken_by_another:
                donor_name = "another template"
                try:
                    donor_name = os.path.basename(
                        self.na_donor_template.get_structure_file(basename_only=False))
                except Exception:
                    pass
                print("[TEMPLATE] '%s': the nucleic acids of this chain are taken from"
                      " '%s', so the chain(s) %s asked for here are left out. Taking one"
                      " filament from one template and another from another would put them"
                      " in unrelated frames."
                      % (copied_template_str_file, donor_name,
                         ", ".join(sorted(taken_by_another))))

        # The file has to be filtered whenever it carries more than this template gives:
        # the chains the user did not ask for, and the ones another template has brought.
        filter_the_nucleic_acids = bool(dict_info_na) and (
            hetres_state == 3 or contributes_na_chains != set(dict_info_na))

        if filter_the_nucleic_acids and contributes_na_chains:
            print(f"\n[FILTER DEBUG] Applying selective filter to nucleic file: {copied_template_str_file}")
            print(f"[FILTER DEBUG] Chains this template contributes: {sorted(contributes_na_chains)}")

            gui_chains_order = sorted(list(dict_info_na.keys()))
            allowed_gui_chains = set(contributes_na_chains)
            
            current_nucleic_index = -1
            last_res_identifier = None
            lines_deleted = 0
            lines_written = 0
            
            with open(source_str_file, "r") as f_in, open(destination_path, "w") as f_out:
                for line in f_in:
                    if line.startswith(("ATOM", "HETATM")):
                        res_num = line[22:26].strip()
                        chain_id = line[21].strip()

                        if is_nucleic_residue_name(line[17:20]):
                            res_identifier = f"{chain_id}_{res_num}"
                            if res_identifier != last_res_identifier:
                                current_nucleic_index += 1
                                last_res_identifier = res_identifier
                            
                            assigned_gui_chain = None
                            accumulated_len = 0
                            for g_chain in gui_chains_order:
                                chain_length = dict_info_na[g_chain][0]
                                if current_nucleic_index < (accumulated_len + chain_length):
                                    assigned_gui_chain = g_chain
                                    break
                                accumulated_len += chain_length
                            
                            # A base which cannot be attributed to any of the chains
                            # goes too. The user asked for these chains and no
                            # others, and keeping what we cannot account for is how
                            # a selected single filament used to come back double.
                            if assigned_gui_chain not in allowed_gui_chains:
                                lines_deleted += 1
                                continue
                                
                    f_out.write(line)
                    lines_written += 1
            print(f"[FILTER DEBUG] Done! Wrote {lines_written} lines. Removed {lines_deleted} lines.\n")
        elif filter_the_nucleic_acids and not contributes_na_chains:
            # Either "Do not use any" was chosen, or every filament this template was
            # asked for is in the model already: nothing of its nucleic acids goes in.
            # Leaving env.io.hetatm off does not do it:
            # MODELLER reads the nucleic acids from the file whatever that flag says,
            # so the model came out with the DNA in it anyway. The only way to leave
            # them out is not to give them to MODELLER at all.
            lines_deleted = 0
            with open(source_str_file, "r") as f_in, open(destination_path, "w") as f_out:
                for line in f_in:
                    if line.startswith(("ATOM", "HETATM")):
                        if is_nucleic_residue_name(line[17:20]):
                            lines_deleted += 1
                            continue
                    f_out.write(line)
            print(f"[FILTER DEBUG] '{copied_template_str_file}': nucleic acids excluded"
                  f" on request, {lines_deleted} lines removed.")

        else:
            # A pure protein, or a template which contributes every filament of the
            # structure: the file goes in as it is.
            if hetres_state == 2:
                print(f"[FILTER DEBUG] Skipping filter for {copied_template_str_file} (pure protein chain).")
            shutil.copy(source_str_file, destination_path)
            
        self.template_options_dict[template]["structure_file"] = copied_template_str_file
        ##########################END#########################

    """
    def prepare_template_files(self, template):
        # Copy the templates structure files in the modeling directory.
        template_str_file = template.get_structure_file(basename_only=False)
        copied_template_str_file = os.path.basename(template_str_file)
        shutil.copy(template_str_file, os.path.join(self.modeling_protocol.modeling_dirpath, copied_template_str_file))
        self.template_options_dict[template]["structure_file"] = copied_template_str_file
        # Build a sequence file for the templates.
        # self.build_modeller_sequence_file(template)
    """
    """
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
    
    """
    ###################MODIFICATION 25/05/2026 - SARA#########################à
    def build_modeller_sequence_file(self, template):
        import os
        import modeller
        
        # Initialize the Modeller environment
        env = modeller.environ()
        modeller.log.none()
        
        # Check the state of heteroatoms and nucleic residues
        opts = self.template_options_dict.get(template, {})
        hetres_state = opts.get("hetres_state", 1)
        
        # If the user chose NOT to use any heteroatoms (state 3), turn off reading
        if hetres_state == 3:
            env.io.hetatm = False
            env.io.water = False
        else:
            env.io.hetatm = True
            if self.modeling_protocol.use_water_in_session:
                env.io.water = True

        structure_file_name = self.template_options_dict[template]["structure_file"]
        structure_file_code = os.path.splitext(structure_file_name)[0]
        
        # Load the PDB file (which has already been filtered by prepare_template_files)
        pdb_path = os.path.join(self.modeling_protocol.modeling_dirpath, structure_file_name)
        mdl = modeller.model(env, file=pdb_path)
        aln = modeller.alignment(env)
        
        # MISMATCH FIX: Find the first and last REAL residues in the modified PDB file
        first_res = mdl.residues[0]
        last_res = mdl.residues[-1]
        
        # Build the exact indices in the format required by Modeller (e.g., "97:D" or "1:X")
        sr_index = f"{first_res.num}:{first_res.chain.name}"
        end_index = f"{last_res.num}:{last_res.chain.name}"
        
        print(f"[MODELLER ALN] Forcing structure boundaries for {structure_file_code}: Start={sr_index}, End={end_index}")
        
        # Define the alignment by explicitly passing sr_index and end_index
        aln.append_model(mdl, 
                         align_codes=structure_file_code, 
                         atom_files=pdb_path,
                         sr_index=sr_index, 
                         end_index=end_index)
                         
        output_sequence_file = structure_file_code + '_aln.chn'
        aln.write(file=os.path.join(self.modeling_protocol.modeling_dirpath, output_sequence_file))
        self.template_options_dict[template]["sequence_file"] = output_sequence_file
    #############################END#####################################

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
