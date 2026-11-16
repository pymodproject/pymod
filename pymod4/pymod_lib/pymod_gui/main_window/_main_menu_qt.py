# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
A module for building the Qt PyMod main window menu.
"""

import sys
import os

from pymol.Qt import QtWidgets, QtGui

from pymod_lib import pymod_vars
from pymod_lib.pymod_gui.shared_gui_components_qt import add_qt_menu_command


class PyMod_main_window_main_menu:
    """
    A class for the Qt PyMod main window menu.
    """

    def make_main_menu(self):
        """
        This method is called when initializing the main window in order to build its main menu.
        """

        self.menubar = self.menuBar()
        self.menubar.setNativeMenuBar(False)

        #---------------
        # "File" menu. -
        #---------------

        self.file_menu = self.menubar.addMenu('File')

        # Sequences submenu.
        self.sequences_submenu = QtWidgets.QMenu('Sequences and Structures', self)
        self.file_menu.addMenu(self.sequences_submenu)
        add_qt_menu_command(self.sequences_submenu, 'Open from File', self.pymod.open_file_from_the_main_menu)
        add_qt_menu_command(self.sequences_submenu, 'Add Raw Sequence', self.pymod.show_raw_seq_input_window)
        add_qt_menu_command(self.sequences_submenu, 'Import PyMOL Objects', self.pymod.import_pymol_selections_from_main_menu)
        add_qt_menu_command(self.sequences_submenu, 'Fetch by UniProt ID', self.pymod.show_uniprot_id_input_window)
        add_qt_menu_command(self.sequences_submenu, 'Fetch by PDB code', self.pymod.fetch_by_pdb_code_from_main_menu)

        # self.file_menu.addSeparator()

        # Submenu to open alignments.
        self.alignment_files_submenu = QtWidgets.QMenu('Alignments', self)
        self.file_menu.addMenu(self.alignment_files_submenu)
        add_qt_menu_command(self.alignment_files_submenu, 'Open from File', self.pymod.open_alignment_from_main_menu)

        # Submenu to open results files from external programs.
        self.external_results_submenu = QtWidgets.QMenu('External Programs Results', self)
        self.file_menu.addMenu(self.external_results_submenu)

        self.hhsuite_results_submenu = QtWidgets.QMenu('HH-suite', self)
        self.external_results_submenu.addMenu(self.hhsuite_results_submenu)
        add_qt_menu_command(self.hhsuite_results_submenu, 'Template Search (.hhr)',
                            self.pymod.open_hhsuite_hhr_from_main_menu)
        add_qt_menu_command(self.hhsuite_results_submenu, 'Multiple Sequence Alignment (.a3m)',
                            self.pymod.open_hhsuite_a3m_from_main_menu)


        self.file_menu.addSeparator()

        # Workspace submenu.
        self.sessions_submenu = QtWidgets.QMenu('Sessions', self)
        self.file_menu.addMenu(self.sessions_submenu)
        add_qt_menu_command(self.sessions_submenu, 'New', self.pymod.start_new_session_from_main_menu)
        add_qt_menu_command(self.sessions_submenu, 'Save', self.pymod.save_session_from_main_menu)
        add_qt_menu_command(self.sessions_submenu, 'Open', self.pymod.open_session_from_main_menu)

        self.file_menu.addSeparator()

        # Exit command.
        add_qt_menu_command(self.file_menu, 'Exit', self.pymod.exit_from_main_menu)


        #----------------
        # "Tools" menu. -
        #----------------

        self.tools_menu = self.menubar.addMenu('Tools')


        # Database search for homologous sequences.
        self.database_search_submenu = QtWidgets.QMenu('Database Search', self)
        self.tools_menu.addMenu(self.database_search_submenu)

        self.blast_tools_submenu = QtWidgets.QMenu('BLAST', self)
        self.database_search_submenu.addMenu(self.blast_tools_submenu)
        add_qt_menu_command(self.blast_tools_submenu, 'Local BLAST', lambda: self.pymod.launch_blast_algorithm("blastp"))
        add_qt_menu_command(self.blast_tools_submenu, 'Remote BLAST', lambda: self.pymod.launch_blast_algorithm("blast"))
        add_qt_menu_command(self.database_search_submenu, 'PSI-BLAST', lambda: self.pymod.launch_blast_algorithm("psi-blast"))
        self.database_search_submenu.addSeparator()
        add_qt_menu_command(self.database_search_submenu, 'phmmer', lambda: self.pymod.launch_hmmer_algorithm("phmmer"))
        add_qt_menu_command(self.database_search_submenu, 'jackhmmer', lambda: self.pymod.launch_hmmer_algorithm("jackhmmer"))
        add_qt_menu_command(self.database_search_submenu, 'hmmsearch', lambda: self.pymod.launch_hmmer_algorithm("hmmsearch"))
        self.database_search_submenu.addSeparator()
        # A search by structure rather than by sequence.
        add_qt_menu_command(self.database_search_submenu, 'Foldseek (structure search)', self.pymod.launch_foldseek_search)


        # Alignment tools
        self.alignments_tools_submenu = QtWidgets.QMenu('Alignment Tools', self)
        self.tools_menu.addMenu(self.alignments_tools_submenu)

        # Sequence alignment tools.
        self.sequence_alignments_submenu = QtWidgets.QMenu('Sequence Alignment', self)
        self.alignments_tools_submenu.addMenu(self.sequence_alignments_submenu)
        add_qt_menu_command(self.sequence_alignments_submenu, 'ClustalW', lambda: self.pymod.launch_alignment_from_the_main_menu("clustalw", "regular"))
        add_qt_menu_command(self.sequence_alignments_submenu, 'Clustal Omega', lambda: self.pymod.launch_alignment_from_the_main_menu("clustalo", "regular"))
        add_qt_menu_command(self.sequence_alignments_submenu, 'MUSCLE', lambda: self.pymod.launch_alignment_from_the_main_menu("muscle", "regular"))
        add_qt_menu_command(self.sequence_alignments_submenu, 'SALIGN (Sequence Alignment)', lambda: self.pymod.launch_alignment_from_the_main_menu("salign-seq", "regular"))

        # Profile alignment tools.
        self.sequence_profile_alignments_submenu = QtWidgets.QMenu('Profile Alignment', self)
        self.alignments_tools_submenu.addMenu(self.sequence_profile_alignments_submenu)
        add_qt_menu_command(self.sequence_profile_alignments_submenu, 'ClustalW', lambda: self.pymod.launch_alignment_from_the_main_menu("clustalw", "profile"))
        add_qt_menu_command(self.sequence_profile_alignments_submenu, 'Clustal Omega', lambda: self.pymod.launch_alignment_from_the_main_menu("clustalo", "profile"))
        add_qt_menu_command(self.sequence_profile_alignments_submenu, 'SALIGN (Sequence Alignment)', lambda: self.pymod.launch_alignment_from_the_main_menu("salign-seq", "profile"))

        # Structural alignment tools.
        self.structural_alignment_submenu = QtWidgets.QMenu('Structural Alignment', self)
        self.alignments_tools_submenu.addMenu(self.structural_alignment_submenu)
        # add_qt_menu_command(self.structural_alignment_submenu, 'Superpose', self.pymod.superpose_from_main_menu)
        add_qt_menu_command(self.structural_alignment_submenu, 'CE Alignment', lambda: self.pymod.launch_alignment_from_the_main_menu("ce", "regular"))
        add_qt_menu_command(self.structural_alignment_submenu, 'SALIGN (Structure Alignment)', lambda: self.pymod.launch_alignment_from_the_main_menu("salign-str", "regular"))


        # Domain tools.
        self.domain_analysis_submenu = QtWidgets.QMenu('Domains Analysis', self)
        self.tools_menu.addMenu(self.domain_analysis_submenu)
        self.hmmscan_tools_submenu = QtWidgets.QMenu('hmmscan', self)
        self.domain_analysis_submenu.addMenu(self.hmmscan_tools_submenu)
        add_qt_menu_command(self.hmmscan_tools_submenu, 'Local hmmscan', lambda: self.pymod.launch_domain_analysis("local"))
        add_qt_menu_command(self.hmmscan_tools_submenu, 'Remote hmmscan', lambda: self.pymod.launch_domain_analysis("remote"))


        # Structural analysis.
        self.structural_analysis_submenu = QtWidgets.QMenu('Structural Analysis', self)
        self.tools_menu.addMenu(self.structural_analysis_submenu)
        add_qt_menu_command(self.structural_analysis_submenu, 'Ramachandran Plot', self.pymod.ramachandran_plot_from_main_menu)
        add_qt_menu_command(self.structural_analysis_submenu, 'Contact Map', self.pymod.contact_map_from_main_menu)
        add_qt_menu_command(self.structural_analysis_submenu, 'Structural Divergence Plot', self.pymod.sda_from_main_menu)
        self.structural_analysis_submenu.addSeparator()
        add_qt_menu_command(self.structural_analysis_submenu, 'Assess with DOPE', self.pymod.dope_from_main_menu)
        self.structural_analysis_submenu.addSeparator()
        add_qt_menu_command(self.structural_analysis_submenu, 'PSIPRED', self.pymod.launch_psipred_from_main_menu)


        # Structure prediction and modeling: everything which produces a structure
        # for a sequence one already has.
        self.modeling_submenu = QtWidgets.QMenu('Structure Modeling', self)
        self.tools_menu.addMenu(self.modeling_submenu)

        # Modeling - MODELLER
        self.modeller_submenu = QtWidgets.QMenu('MODELLER', self)
        self.modeling_submenu.addMenu(self.modeller_submenu)
        add_qt_menu_command(self.modeller_submenu, "Homology modeling", self.pymod.launch_modeller_hm_from_main_menu)
        add_qt_menu_command(self.modeller_submenu, "Loop modeling", self.pymod.launch_modeller_lr_from_main_menu)

        # Modeling - AlphaFold
        self.alphafold_submenu = QtWidgets.QMenu('AlphaFold', self)
        self.modeling_submenu.addMenu(self.alphafold_submenu)
        add_qt_menu_command(self.alphafold_submenu, "Fetch from database", self.pymod.launch_afdb_from_main_menu)

        # Modeling - the ESM models. The submenu is not called ESMFold any more: it now
        # holds a database as well as two predictors, and a reader who saw "ESMFold" with
        # two entries under it saying "no key needed" had no way to tell what the
        # difference between them was. Each entry now says what it does first -- look a
        # structure up, or compute one -- and which model or database it uses.
        self.esm_submenu = QtWidgets.QMenu('ESM', self)
        self.modeling_submenu.addMenu(self.esm_submenu)

        # Looking a structure up: free, immediate, and only for a sequence which is in
        # the database already.
        add_qt_menu_command(self.esm_submenu, "Look up in the ESM Atlas database (no key)",
                            self.pymod.launch_esm_atlas_from_main_menu)

        # Computing one, for a sequence which is in no database: a designed binder, a
        # mutant, a fragment nobody has deposited.
        add_qt_menu_command(self.esm_submenu, "Predict with ESMFold2 (API key)",
                            self.pymod.launch_forge_fold_from_main_menu)

        # The keyless prediction is deliberately not offered any more. It folded the
        # sequence at api.esmatlas.com, and that service has been answering 504 for the
        # whole of a working week: the gateway is up and the front end of esmatlas.com
        # still calls the same endpoint, so the route is broken rather than retired, but a
        # command which always fails is worse than no command. What it was for -- a
        # structure without a key -- is now served by the lookup above, since the Atlas
        # holds the ESMFold2 prediction of anything which has been deposited.
        #
        # The protocol itself is kept, in ESMFold/__init__.py, and it takes only one line
        # here to offer the prediction again if the service comes back.
        # add_qt_menu_command(self.esm_submenu, "Predict with ESMFold (no key)",
        #                     self.pymod.launch_esm_from_main_menu)

        # Modeling - Boltz-2
        self.boltz_submenu = QtWidgets.QMenu('Boltz-2', self)
        self.modeling_submenu.addMenu(self.boltz_submenu)
        add_qt_menu_command(self.boltz_submenu, "Predict structure and binding (API key)", self.pymod.launch_boltz_from_main_menu)
        self.boltz_submenu.addSeparator()
        add_qt_menu_command(self.boltz_submenu, "Import a finished prediction by id...", self.pymod.import_boltz_prediction_by_id_from_main_menu)

        # Modeling - refining a structure one already has, rather than building one.
        self.modeling_submenu.addSeparator()
        add_qt_menu_command(self.modeling_submenu, "Optimize Interface", self.pymod.optimize_from_main_menu)
        add_qt_menu_command(self.modeling_submenu, "Optimize Molecule", self.pymod.launch_molecule_optimizer)


        # Design: everything which invents a sequence, rather than predicting the
        # structure of one which already exists. Kept apart from the modeling
        # commands, since the two answer different questions.
        self.design_submenu = QtWidgets.QMenu('Protein Design', self)
        self.tools_menu.addMenu(self.design_submenu)

        # Design - a sequence for a structure one already has.
        add_qt_menu_command(self.design_submenu, "Design a sequence with ESM3 (API key)", self.pymod.launch_inverse_fold_from_main_menu)
        # Design - a new binder against a target.
        add_qt_menu_command(self.design_submenu, "Design a binder with BoltzProt-1 (API key)", self.pymod.launch_boltz_design_from_main_menu)

        # Options.
        self.tools_menu.addSeparator()
        add_qt_menu_command(self.tools_menu, 'Options', self.pymod.show_pymod_options_window)


        #---------------------
        # "Alignments" menu. -
        #---------------------

        self.alignments_menu = self.menubar.addMenu('Alignments')
        self.build_alignment_submenu()

        #####################Add ETE code 21/05/2026-SARA##############
        #---------------------
        # "Trees" menu. -
        #---------------------

        self.trees_menu = self.menubar.addMenu('Trees')
        self.build_trees_submenu()
        ##################END ETE######################
        #-----------------
        # "Models" menu. -
        #-----------------

        # When the plugin is started there are no models.
        self.models_menu = self.menubar.addMenu('Models')
        self.build_models_submenu()


        #--------------------
        # "Selection" menu. -
        #--------------------

        self.main_selection_menu = self.menubar.addMenu('Selection')
        add_qt_menu_command(self.main_selection_menu, 'Select All [Ctrl+a]', self.pymod.select_all_from_main_menu)
        add_qt_menu_command(self.main_selection_menu, 'Deselect All [Esc]', self.pymod.deselect_all_from_main_menu)

        # Structures selection submenu.
        self.selection_structures_menu = QtWidgets.QMenu('Structures', self)
        self.main_selection_menu.addMenu(self.selection_structures_menu)
        add_qt_menu_command(self.selection_structures_menu, 'Show All in PyMOL [Ctrl+s]', self.pymod.show_all_structures_from_main_menu)
        add_qt_menu_command(self.selection_structures_menu, 'Hide All in PyMOL [Ctrl+h]', self.pymod.hide_all_structures_from_main_menu)
        self.selection_structures_menu.addSeparator()
        add_qt_menu_command(self.selection_structures_menu, 'Select All', self.pymod.select_all_structures_from_main_menu)
        add_qt_menu_command(self.selection_structures_menu, 'Deselect All', self.pymod.deselect_all_structures_from_main_menu)
        # Clusters selection submenu.
        self.selection_clusters_menu = QtWidgets.QMenu('Clusters', self)
        self.main_selection_menu.addMenu(self.selection_clusters_menu)
        add_qt_menu_command(self.selection_clusters_menu, 'Expand All', self.pymod.expand_all_clusters_from_main_menu)
        add_qt_menu_command(self.selection_clusters_menu, 'Collapse All', self.pymod.collapse_all_clusters_from_main_menu)


        #------------------
        # "Display" menu. -
        #------------------

        self.display_menu = self.menubar.addMenu('Display')
        if self.pymod.DEVELOP:
            add_qt_menu_command(self.display_menu, 'Print Selected Sequences', command=self.print_selected)

        # Color menu.
        self.main_color_menu = QtWidgets.QMenu('Color all Sequences', self)
        self.display_menu.addMenu(self.main_color_menu)
        add_qt_menu_command(self.main_color_menu, 'By Default Color', command=lambda: self.color_selection("all", None, "regular"))
        self.main_residues_colors_menu = QtWidgets.QMenu('By Residue Properties', self)
        self.main_color_menu.addMenu(self.main_residues_colors_menu)
        add_qt_menu_command(self.main_residues_colors_menu, 'Residue Type', command=lambda: self.color_selection("all", None, "residue_type"))
        add_qt_menu_command(self.main_residues_colors_menu, 'Polarity', command=lambda: self.color_selection("all", None, "polarity"))
        add_qt_menu_command(self.main_color_menu, 'By Secondary Structure', command=lambda: self.color_selection("all", None, "secondary-auto"))

        # Font selection.
        self.font_selection_menu = QtWidgets.QMenu('Font Selection', self)
        self.display_menu.addMenu(self.font_selection_menu)
        if self.pymod.DEVELOP:
            add_qt_menu_command(self.display_menu, 'Font Type and Size', self.pymod.change_font_from_main_menu)

        # Font size selection.
        self.font_size_selection_menu = QtWidgets.QMenu('Font Size', self)
        self.font_selection_menu.addMenu(self.font_size_selection_menu)
        font_size_action_group = QtWidgets.QActionGroup(self.font_size_selection_menu)
        font_size_action_group.setExclusive(True)
        for font_size in ("6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"):
            action = font_size_action_group.addAction(QtWidgets.QAction(font_size,
                                                                        self.font_size_selection_menu,
                                                                        checkable=True))
            self.font_size_selection_menu.addAction(action)
            if font_size == str(self.font_size):
                action.setChecked(True)
            action.triggered.connect(lambda a=None, t=None, s=font_size: self.pymod.change_font_from_action(t, s))

        # Font type selection.
        self.font_type_selection_menu = QtWidgets.QMenu('Font Type', self)
        self.font_selection_menu.addMenu(self.font_type_selection_menu)
        font_type_action_group = QtWidgets.QActionGroup(self.font_type_selection_menu)
        font_type_action_group.setExclusive(True)
        for font_type in self.get_available_fonts():
            action = font_type_action_group.addAction(QtWidgets.QAction(font_type,
                                                                        self.font_type_selection_menu,
                                                                        checkable=True))
            self.font_type_selection_menu.addAction(action)
            if font_type == str(self.font):
                action.setChecked(True)
            action.triggered.connect(lambda a=None, t=font_type, s=None: self.pymod.change_font_from_action(t, s))


        #---------------
        # "Help" menu. -
        #---------------

        self.help_menu = self.menubar.addMenu('Help')
        if self.pymod.DEVELOP:
            add_qt_menu_command(self.help_menu, 'Test', self.pymod.launch_default)
        add_qt_menu_command(self.help_menu, 'Online Documentation', self.pymod.open_online_documentation)
        add_qt_menu_command(self.help_menu, 'About', self.pymod.show_about_dialog)
        self.help_menu.addSeparator()
        add_qt_menu_command(self.help_menu, "Install PyMod Components", self.pymod.launch_components_installation)
        add_qt_menu_command(self.help_menu, "Check for Database Updates", self.pymod.launch_databases_update)

        if self.pymod.TEST:
            self.help_menu.addSeparator()
            self.examples_menu = QtWidgets.QMenu('Examples', self)
            self.help_menu.addMenu(self.examples_menu)
            add_qt_menu_command(self.examples_menu, "Load random sequence from UniProt", lambda a=None: self.pymod.load_uniprot_random())
            add_qt_menu_command(self.examples_menu, "Load random PDB entry", lambda a=None: self.pymod.load_pdb_random())


    selected_fonts_dict = {"linux": ["Andale Mono",
                                     "Courier",
                                     "DejaVu Sans Mono",
                                     "Liberation Mono",
                                     "Monospace",
                                     "Ubuntu Mono",],
                           "win32": ["Consolas",
                                     # "Courier", # Doesn't seem to get the character size correctly.
                                     "Lucida Sans Typewriter"],
                           "darwin": [
                                      # "consolas",
                                      "courier",
                                      "liberation mono",
                                      "lucida console",
                                      "menlo",
                                      # "monaco",
                                      "monofur",
                                      "prestige elite",
                                      "roboto mono"]}
    def get_available_fonts(self):
        # Get the system fonts.
        font_family_list = [font_name.lower() for font_name in QtGui.QFontDatabase().families()]
        # Get the fonts available in PyMod.
        selected_fonts = [font_name.lower() for font_name in self.selected_fonts_dict.get(sys.platform, [])]
        # Filter the system fonts list.
        font_family_list = [font_name for font_name in font_family_list if font_name in selected_fonts]
        font_family_list.insert(0, "courier new")
        font_family_list.sort()
        return font_family_list


    #########################################################################################
    # Build submenus which get populated with elements when performing operations in PyMod. #
    #########################################################################################

    def build_alignment_submenu(self):
        """
        Build an "Alignment N" voice in the "Alignments" submenu when alignment N is performed.
        """

        # Delete the old alignment submenu.
        self.alignments_menu.clear()

        # Then rebuilds it with the new alignments.
        alignment_list = self.pymod.get_cluster_elements()

        if alignment_list != []:

            for alignment_element in alignment_list:

                # Adds the alignment submenu for each cluster loaded in PyMod to the PyMod main menu.
                label_text = alignment_element.my_header
                single_alignment_submenu = QtWidgets.QMenu(label_text, self)
                self.alignments_menu.addMenu(single_alignment_submenu)

                # Save to a file dialog.
                add_qt_menu_command(single_alignment_submenu, "Save to File",
                                    # The first argument in Qt is a 'False' value, therefore we need to add a dummy 'a' variable
                                    # in order to pass the 'alignment_element' object.
                                    lambda a=None, e=alignment_element: self.pymod.save_alignment_to_file_from_ali_menu(e))
                single_alignment_submenu.addSeparator()


                # Matrices submenu.
                single_alignment_matrices_submenu = QtWidgets.QMenu('Matrices', self)
                single_alignment_submenu.addMenu(single_alignment_matrices_submenu)

                add_qt_menu_command(single_alignment_matrices_submenu, "Sequence Identity Matrix", lambda a=None, e=alignment_element: self.pymod.display_identity_matrix(e))
                add_qt_menu_command(single_alignment_matrices_submenu, "RMSD Matrix", lambda a=None, e=alignment_element: self.pymod.display_rmsd_matrix_from_alignments_menu(e))

                ##############Add ETE code 21/05/2026-SARA#################################
                # Trees.
                
                if len(alignment_element.get_children()) > 2:
                    # single_alignment_trees_submenu = QtWidgets.QMenu('Trees', self)
                    # single_alignment_submenu.addMenu(single_alignment_trees_submenu)

                    # if alignment_element.algorithm in pymod_vars.can_show_guide_tree:
                    #     add_qt_menu_command(single_alignment_trees_submenu, "Show Guide Tree", lambda a=None, e=alignment_element: self.pymod.show_guide_tree_from_alignments_menu(e))
                    # if alignment_element.algorithm in pymod_vars.can_show_dendrogram and alignment_element.tree_file_path is not None:
                    #     add_qt_menu_command(single_alignment_trees_submenu, "Show SALIGN Dendrogram", lambda a=None, e=alignment_element: self.pymod.show_dendrogram_from_alignments_menu(e))
                    # if len(alignment_element.get_children()) >= 2:
                    #     add_qt_menu_command(single_alignment_trees_submenu, "Build Tree from Alignment", lambda a=None, e=alignment_element: self.pymod.build_tree_from_alignments_menu(e))
                    if alignment_element.algorithm in pymod_vars.can_show_guide_tree:
                        add_qt_menu_command(single_alignment_submenu, "Show Guide Tree", lambda a=None, e=alignment_element: self.pymod.show_guide_tree_from_alignments_menu(e))
                    if alignment_element.algorithm in pymod_vars.can_show_dendrogram and alignment_element.tree_file_path is not None:
                        add_qt_menu_command(single_alignment_submenu, "Show SALIGN Dendrogram", lambda a=None, e=alignment_element: self.pymod.show_dendrogram_from_alignments_menu(e))
                    if len(alignment_element.get_children()) >= 2:
                        add_qt_menu_command(single_alignment_submenu, "Build Tree from Alignment", lambda a=None, e=alignment_element: self.pymod.build_tree_from_alignments_menu(e))
                ############################ETE end#######################

                #########original#######
                """
                # Trees.
                if alignment_element.initial_number_of_sequences > 2:
                    single_alignment_trees_submenu = QtWidgets.QMenu('Trees', self)
                    single_alignment_submenu.addMenu(single_alignment_trees_submenu)

                    if alignment_element.algorithm in pymod_vars.can_show_guide_tree:
                        add_qt_menu_command(single_alignment_trees_submenu, "Show Guide Tree", lambda a=None, e=alignment_element: self.pymod.show_guide_tree_from_alignments_menu(e))
                    if alignment_element.algorithm in pymod_vars.can_show_dendrogram and alignment_element.tree_file_path is not None:
                        add_qt_menu_command(single_alignment_trees_submenu, "Show SALIGN Dendrogram", lambda a=None, e=alignment_element: self.pymod.show_dendrogram_from_alignments_menu(e))
                    if len(alignment_element.get_children()) >= 2:
                        add_qt_menu_command(single_alignment_trees_submenu, "Build Tree from Alignment", lambda a=None, e=alignment_element: self.pymod.build_tree_from_alignments_menu(e))
                """
                ###################END################

                # Evolutionary conservation.
                single_alignment_evolutionary_submenu = QtWidgets.QMenu('Evolutionary Conservation', self)
                single_alignment_submenu.addMenu(single_alignment_evolutionary_submenu)

                single_seq_evolutionary_submenu = QtWidgets.QMenu('Sequence Conservation', self)
                single_alignment_evolutionary_submenu.addMenu(single_seq_evolutionary_submenu)
                add_qt_menu_command(single_seq_evolutionary_submenu, "Entropy", lambda a=None, e=alignment_element: self.pymod.launch_entropy_scorer_from_main_menu(e))
                add_qt_menu_command(single_seq_evolutionary_submenu, "CAMPO", lambda a=None, e=alignment_element: self.pymod.launch_campo_from_main_menu(e))
                add_qt_menu_command(single_seq_evolutionary_submenu, "Pair Conservation", lambda a=None, e=alignment_element: self.pymod.launch_pc_from_main_menu(e))

                single_str_evolutionary_submenu = QtWidgets.QMenu('Structural Conservation', self)
                single_alignment_evolutionary_submenu.addMenu(single_str_evolutionary_submenu)
                add_qt_menu_command(single_str_evolutionary_submenu, "SCR_FIND", lambda a=None, e=alignment_element: self.pymod.launch_scr_find_from_main_menu(e))


                # Render alignment.
                single_alignment_render_submenu = QtWidgets.QMenu('Render Alignment', self)
                single_alignment_submenu.addMenu(single_alignment_render_submenu)

                add_qt_menu_command(single_alignment_render_submenu, "Render Figure (local)", lambda a=None, e=alignment_element: self.pymod.launch_alignment_figure_from_main_menu(e))
                add_qt_menu_command(single_alignment_render_submenu, "Generate Logo through WebLogo 3", lambda a=None, e=alignment_element: self.pymod.launch_weblogo_from_main_menu(e))

                # The ESPript entry was removed from this submenu: that protocol uploads the
                # alignment to a relay on the schubert server and then opens the ESPript CGI
                # in a web browser, a chain which breaks whenever one of those servers
                # changes and which never returns the figure to PyMod. "Render Figure
                # (local)" above draws the alignment in PyMod instead, with no server
                # involved. The protocol itself is kept in 'espript.py', which explains this
                # in detail: restoring the command below is all it takes to offer it again.
                # add_qt_menu_command(single_alignment_render_submenu, "Launch ESPript in Web Browser", lambda a=None, e=alignment_element: self.pymod.launch_espript_from_main_menu(e))

        else:

            add_qt_menu_command(self.alignments_menu, "There aren't any alignments")


    def build_models_submenu(self):
        """
        Build an "Modeling Session n" voice in the "Models" submenu once some models have been
        built.
        """

        # Delete the old models submenu.
        self.models_menu.clear()

        self.modeller_submenu_models = QtWidgets.QMenu('MODELLER', self)
        self.models_menu.addMenu(self.modeller_submenu_models)

        self.alphafold_submenu_models = QtWidgets.QMenu('AlphaFold', self)
        self.models_menu.addMenu(self.alphafold_submenu_models)

        self.esmfold_submenu_models = QtWidgets.QMenu('ESMFold', self)
        self.models_menu.addMenu(self.esmfold_submenu_models)

        self.boltz_submenu_models = QtWidgets.QMenu('Boltz-2', self)
        self.models_menu.addMenu(self.boltz_submenu_models)

        # Structures retrieved from the ESM Atlas. They are models like the others, and
        # they carry figures of their own -- what ESMFold2 was confident of, and what ESMC
        # recognises in the protein -- so they are listed under their own heading.
        self.esm_atlas_submenu_models = QtWidgets.QMenu('ESM Atlas', self)
        self.models_menu.addMenu(self.esm_atlas_submenu_models)
        
        if self.pymod.modeling_session_list != []:
            for modeling_session in self.pymod.modeling_session_list:
                # Adds a modeling session submenu to the PyMod main menu.
                label_text = "Modeling Session %s" % (modeling_session.session_id)
                modeling_session_submenu = QtWidgets.QMenu(label_text, self)
                self.modeller_submenu_models.addMenu(modeling_session_submenu)
                add_qt_menu_command(modeling_session_submenu, "Export to File", lambda a=None, ms=modeling_session: self.pymod.save_modeling_session(ms))
                add_qt_menu_command(modeling_session_submenu, "DOPE Profile", lambda a=None, ms=modeling_session: self.pymod.show_session_profile(ms))
                add_qt_menu_command(modeling_session_submenu, "Assessment Table", lambda a=None, ms=modeling_session: self.pymod.show_assessment_table(ms))

        else:
            add_qt_menu_command(self.modeller_submenu_models, "There aren't any models")

    
        ##########################MODIFIED on 28/02/2025########################################
        ######The name of the model in section Models<AlphaFold is the real name visible in pymod 
        ######The plDDT plot working properly
        if self.pymod.af_modeling_session_list:
            for af_session, models in self.pymod.af_modeling_session_list.items():
                if models:
                    for model_name in models:
                        label_text = "%s" % model_name
                        af_session_submenu = QtWidgets.QMenu(label_text, self)
                        self.alphafold_submenu_models.addMenu(af_session_submenu)
                        add_qt_menu_command(af_session_submenu, "pLDDT", lambda a=None, ms=af_session: self.pymod.show_plddt_plot(ms))
                        add_qt_menu_command(af_session_submenu, "Predicted Aligned Error", lambda a=None, ms=af_session: self.pymod.show_pae_plot(ms))
        else:
            add_qt_menu_command(self.alphafold_submenu_models, "There aren't any models")


        """"OLD
        if self.pymod.af_modeling_session_list:
            for af_session, models in self.pymod.af_modeling_session_list.items():
                if models:
                    for model_name in models:
                        ##### MODIFIED #####
                        if self.pymod.af_modeling_session_list != {}:
                        ##### END #####
                            for af_sessions in self.pymod.af_modeling_session_list:
                                ##### MODIFIED #####
                                if self.pymod.af_modeling_session_list[af_sessions] != {}:
                                ##### END #####
                                    label_text = "%s" % model_name
                                    af_session_submenu = QtWidgets.QMenu(label_text, self)
                                    self.alphafold_submenu_models.addMenu(af_session_submenu)
                                    add_qt_menu_command(af_session_submenu, "pLDDT", lambda a=None, ms=af_sessions: self.pymod.show_plddt_plot(ms))
                                    add_qt_menu_command(af_session_submenu, "Predicted Aligned Error", lambda a=None, ms=af_sessions: self.pymod.show_pae_plot(ms))
                        else:
                            add_qt_menu_command(self.alphafold_submenu_models, "There aren't any models")
        """
        """OLD CODE --> name model is the af_session
        ##### MODIFIED #####
        if self.pymod.af_modeling_session_list != {}:
        ##### END #####
            for af_sessions in self.pymod.af_modeling_session_list:
                ##### MODIFIED #####
                if self.pymod.af_modeling_session_list[af_sessions] != {}:
                ##### END #####
                    label_text = "AF model: %s" % af_sessions
                    af_session_submenu = QtWidgets.QMenu(label_text, self)
                    self.alphafold_submenu_models.addMenu(af_session_submenu)
                    add_qt_menu_command(af_session_submenu, "pLDDT", lambda a=None, ms=af_sessions: self.pymod.show_plddt_plot(ms))
                    add_qt_menu_command(af_session_submenu, "Predicted Aligned Error", lambda a=None, ms=af_sessions: self.pymod.show_pae_plot(ms))
        """
        

        # if self.pymod.af_modeling_session_list != {}:
        #     for af_sessions in self.pymod.af_modeling_session_list:
        #         ##### MODIFIED #####
        #         if self.pymod.af_modeling_session_list[af_sessions] != {}:
        #             for af_models in self.pymod.af_modeling_session_list[af_sessions]:
        #                 label_text = "AF model: " + af_models
        #                 af_models_submenu = QtWidgets.QMenu(label_text, self)
        #                 self.alphafold_submenu_models.addMenu(af_models_submenu)
        #                 add_qt_menu_command(af_models_submenu, "pLDDT", lambda a=None, ms=af_models: self.pymod.show_plddt_plot(ms))
        #                 add_qt_menu_command(af_models_submenu, "Predicted Aligned Error", lambda a=None, ms=af_models: self.pymod.show_pae_plot(ms))
        #         ##### END #####
        #
        # else:
        #     add_qt_menu_command(self.alphafold_submenu_models, "There aren't any models")



        #######################MODIFICATION on 05/07/2025######################################
        ##################ESMFold models in the Models option######################
        """
        if self.pymod.esm_modeling_session_list:
            for af_session, models in self.pymod.esm_modeling_session_list.items():
                if models:
                    for model_name in models:
                        ##### MODIFIED #####
                        if self.pymod.esm_modeling_session_list != []:
                        ##### END #####
        
                            for af_sessions in self.pymod.af_modeling_session_list:
                                ##### MODIFIED #####
        """

        #########MODELLER CASE############
        """
        if self.pymod.modeling_session_list != []:
            for modeling_session in self.pymod.modeling_session_list:
                # Adds a modeling session submenu to the PyMod main menu.
                label_text = "Modeling Session %s" % (modeling_session.session_id)
                modeling_session_submenu = QtWidgets.QMenu(label_text, self)
                self.modeller_submenu_models.addMenu(modeling_session_submenu)
                add_qt_menu_command(modeling_session_submenu, "Export to File", lambda a=None, ms=modeling_session: self.pymod.save_modeling_session(ms))
                add_qt_menu_command(modeling_session_submenu, "DOPE Profile", lambda a=None, ms=modeling_session: self.pymod.show_session_profile(ms))
                add_qt_menu_command(modeling_session_submenu, "Assessment Table", lambda a=None, ms=modeling_session: self.pymod.show_assessment_table(ms))

        else:
            add_qt_menu_command(self.modeller_submenu_models, "There aren't any models")
        """
        ##############MODIFICATIO on 11/07/2025################à
        if self.pymod.esm_modeling_session_list:
            for i, esm_session in enumerate(self.pymod.esm_modeling_session_list, start=1):
                base_name = getattr(esm_session, "model_name", "Unnamed_ESM_Model")
                label_text = base_name  

                esm_session_submenu = QtWidgets.QMenu(label_text, self)
                self.esmfold_submenu_models.addMenu(esm_session_submenu)
                add_qt_menu_command(esm_session_submenu, "pLDDT", lambda a=None, ms=esm_session: self.pymod.show_plddt_plot_ESM(ms))
                # Only models that carry Predicted Aligned Error data (Forge fold
                # with include_pae) get a PAE plot entry.
                if getattr(esm_session, "pae_matrix", None) is not None:
                    add_qt_menu_command(esm_session_submenu, "Predicted Aligned Error", lambda a=None, ms=esm_session: self.pymod.show_pae_plot_ESM(ms))
        else:
            add_qt_menu_command(self.esmfold_submenu_models, "There aren't any models")

        ##########END###########

        # Boltz-2 predictions. Each carries the figures the service reported, the
        # per-residue confidence of its model, and the predicted aligned error when
        # the results archive was downloaded with it.
        boltz_sessions = getattr(self.pymod, "boltz_modeling_session_list", [])
        if boltz_sessions:
            for boltz_session in boltz_sessions:
                label_text = getattr(boltz_session, "model_name", "Boltz model")
                boltz_session_submenu = QtWidgets.QMenu(label_text, self)
                self.boltz_submenu_models.addMenu(boltz_session_submenu)
                add_qt_menu_command(boltz_session_submenu, "Confidence figures",
                                    lambda a=None, ms=boltz_session: self.pymod.show_boltz_confidence_figures(ms))
                add_qt_menu_command(boltz_session_submenu, "pLDDT",
                                    lambda a=None, ms=boltz_session: self.pymod.show_plddt_plot_boltz(ms))
                if getattr(boltz_session, "pae_matrix", None) is not None:
                    add_qt_menu_command(boltz_session_submenu, "Predicted Aligned Error",
                                        lambda a=None, ms=boltz_session: self.pymod.show_pae_plot_boltz(ms))
        else:
            add_qt_menu_command(self.boltz_submenu_models, "There aren't any models")

        atlas_models = getattr(self.pymod, "esm_atlas_model_list", [])
        if atlas_models:
            for atlas_model in atlas_models:
                label_text = atlas_model.get("model_name", "ESM Atlas model")
                atlas_submenu = QtWidgets.QMenu(label_text, self)
                self.esm_atlas_submenu_models.addMenu(atlas_submenu)
                # The figures the Atlas reports -- the mean confidence, the pTM and the
                # interpretable features of ESMC -- are not offered yet. The command which
                # shows them is written and tested, in show_esm_atlas_figures(), so it takes
                # this one line to offer them when we come back to how they should look.
                # add_qt_menu_command(atlas_submenu, "Figures of the Atlas",
                #                     lambda a=None, am=atlas_model: self.pymod.show_esm_atlas_figures(am))
                add_qt_menu_command(atlas_submenu, "pLDDT",
                                    lambda a=None, am=atlas_model: self.pymod.show_esm_atlas_plddt(am))
        else:
            add_qt_menu_command(self.esm_atlas_submenu_models, "There aren't any models")

        """
        if self.pymod.esm_modeling_session_list != []:
        ##### END #####
            for esm_session in self.pymod.esm_modeling_session_list:
                label_text = "%s" % model_name
                #label_text = getattr(esm_session, "model_name", "Unnamed_ESM_Model")
                esm_session_submenu = QtWidgets.QMenu(label_text, self)
                self.esmfold_submenu_models.addMenu(esm_session_submenu)
                add_qt_menu_command(esm_session_submenu, "pLDDT", lambda a=None, ms=esm_session: self.pymod.show_plddt_plot(ms))
                add_qt_menu_command(esm_session_submenu, "Predicted Aligned Error", lambda a=None, ms=esm_session: self.pymod.show_pae_plot(ms))
        else:
            add_qt_menu_command(self.esmfold_submenu_models, "There aren't any models")
        """



        ##############OLD CODE ########################
        """
        if self.pymod.esm_modeling_session_list != []:
            pass
        else:
            add_qt_menu_command(self.esmfold_submenu_models, "There aren't any models")
        """




        # The Boltz submenu is filled in above, from 'boltz_modeling_session_list'.
        # A block here used to add "There aren't any models" to it whatever those
        # sessions were, keyed on the OmegaFold list from which this submenu had been
        # copied, so the entry appeared next to the models themselves. It is removed.


###############Add ETE code 21/05/2026-SARA###################
    def build_trees_submenu(self):

            """
            Build a "Trees N" voice in the "Tree" submenu when tree building N is performed.
            """

            # Delete the old trees submenu.
            self.trees_menu.clear()

            add_qt_menu_command(self.trees_menu, "Open from File", self.pymod.open_tree_from_main_menu)

            self.trees_menu.addSeparator()

            # Then rebuilds it with the new alignments.
            alignment_list = self.pymod.get_cluster_elements()
            for alignment_element in alignment_list:

                if not hasattr(alignment_element, "trees_list"):
                    continue
                if alignment_element.trees_list:
                    # Adds the tree submenu for each cluster loaded in PyMod to the PyMod main menu.
                    label_text_aln = "Trees from %s" % alignment_element.my_header
                    single_aln_trees_submenu = QtWidgets.QMenu(label_text_aln, self)
                    self.trees_menu.addMenu(single_aln_trees_submenu)

                    label_list = []
                    for tree_j in alignment_element.trees_list:
                        if tree_j['algorithm'] not in label_list:
                            label_list.append(tree_j['algorithm'])

                            # Parsimony has not calculator option.
                            if next(iter(tree_j.values()))[0] == "NJ" or next(iter(tree_j.values()))[0] == "UPGMA":
                                label_text_tree = "Tree built with %s and %s matrix" % (next(iter(tree_j.values()))[0],next(iter(tree_j.values()))[1])
                                #label_text_tree = "Tree built with %s" % ' and '.join(map(str,tree_j['algorithm']))
                            else:
                                label_text_tree = "Tree built with %s" % tree_j['algorithm']

                            # Add single tree submenu.
                            single_tree_submenu = QtWidgets.QMenu(label_text_tree, self)
                            single_aln_trees_submenu.addMenu(single_tree_submenu)

                            add_qt_menu_command(single_tree_submenu, "View Plot", lambda a = None, e= tree_j['filepath']: self.pymod.visualize_tree(e))
                            ############### MODIFICATION SARA on 31/07/2025 ##############
                            # Added the "Remove Tree" menu command to allow users to remove a specific tree
                            # identified by its filepath (tree_j['filepath']) from the loaded trees list.
                            # This triggers the remove_tree function in pymod to update the internal data
                            # and refresh the menu accordingly.
                            add_qt_menu_command(single_tree_submenu, "Remove Tree", lambda a=None, e=tree_j['filepath']: self.pymod.remove_tree(e))
                            ############### END ###################
                            # Save to a file dialog.
                            add_qt_menu_command(single_tree_submenu, "Save Newick to File",
                                                # The first argument in Qt is a 'False' value, therefore we need to add a dummy 'a' variable
                                                # in order to pass the 'alignment_element' object.
                                                lambda a=None, e=tree_j['filepath']: self.pymod.save_tree_from_trees_menu(e))

            # Imported trees.
            if self.pymod.loaded_trees_list:
                for ltrees in range(len(self.pymod.loaded_trees_list)):
                    label_text = "Imported tree: %s" % (os.path.basename(self.pymod.loaded_trees_list[ltrees]))
                    single_tree_submenu = QtWidgets.QMenu(label_text, self)
                    self.trees_menu.addMenu(single_tree_submenu)
                    
                    
                    """old code , message moved to main_menu_commands 11/04/2026
                    title = "Importing tree"
                    message = "The tree has been imported successfully."
                    self.pymod.main_window.show_info_message(title, message)

                    """
                    add_qt_menu_command(single_tree_submenu, "View Plot", lambda a=None, e=self.pymod.loaded_trees_list[ltrees]: self.pymod.visualize_imported_tree(e))
                    ############### MODIFICATION SARA on 31/07/2025 ##############
                    # Added the "Remove Tree" option for imported trees submenu entries.
                    # This connects the menu action to the remove_tree function, passing the
                    # loaded tree path from loaded_trees_list by index ltrees, enabling removal
                    # of imported trees and automatic update of the tree menu.
                    add_qt_menu_command(
                        single_tree_submenu,
                        "Remove Tree",
                        lambda a=None, e=self.pymod.loaded_trees_list[ltrees]: self.pymod.remove_tree(e)
                    )
                    ###################### END ##################
###############Add ETE code END###################