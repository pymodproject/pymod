# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
A module for building the Qt PyMod main window menu.
"""

import sys

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

        # self.file_menu.addSeparator()

        # Submenu to open alignments.
        self.alignment_files_submenu = QtWidgets.QMenu('Alignments', self)
        self.file_menu.addMenu(self.alignment_files_submenu)
        add_qt_menu_command(self.alignment_files_submenu, 'Open from File', self.pymod.open_alignment_from_main_menu)

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


        # Modeling.
        self.modeling_submenu = QtWidgets.QMenu('Modeling', self)
        self.tools_menu.addMenu(self.modeling_submenu)
        add_qt_menu_command(self.modeling_submenu, "MODELLER (Homology Modeling)", self.pymod.launch_modeller_hm_from_main_menu)
        add_qt_menu_command(self.modeling_submenu, "MODELLER (Loop Modeling)", self.pymod.launch_modeller_lr_from_main_menu)


        # Options.
        self.tools_menu.addSeparator()
        add_qt_menu_command(self.tools_menu, 'Options', self.pymod.show_pymod_options_window)


        #---------------------
        # "Alignments" menu. -
        #---------------------

        self.alignments_menu = self.menubar.addMenu('Alignments')
        self.build_alignment_submenu()


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
        for font_size in ("6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"):
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

                add_qt_menu_command(single_alignment_render_submenu, "Generate Logo through WebLogo 3", lambda a=None, e=alignment_element: self.pymod.launch_weblogo_from_main_menu(e))
                add_qt_menu_command(single_alignment_render_submenu, "Launch ESPript in Web Browser", lambda a=None, e=alignment_element: self.pymod.launch_espript_from_main_menu(e))

        else:

            add_qt_menu_command(self.alignments_menu, "There aren't any alignments")


    def build_models_submenu(self):
        """
        Build an "Modeling Session n" voice in the "Models" submenu once some models have been
        built.
        """

        # Delete the old models submenu.
        self.models_menu.clear()

        if self.pymod.modeling_session_list != []:
            for modeling_session in self.pymod.modeling_session_list:
                # Adds a modeling session submenu to the PyMod main menu.
                label_text = "Modeling Session %s" % (modeling_session.session_id)
                modeling_session_submenu = QtWidgets.QMenu(label_text, self)
                self.models_menu.addMenu(modeling_session_submenu)
                add_qt_menu_command(modeling_session_submenu, "Export to File", lambda a=None, ms=modeling_session: self.pymod.save_modeling_session(ms))
                add_qt_menu_command(modeling_session_submenu, "DOPE Profile", lambda a=None, ms=modeling_session: self.pymod.show_session_profile(ms))
                add_qt_menu_command(modeling_session_submenu, "Assessment Table", lambda a=None, ms=modeling_session: self.pymod.show_assessment_table(ms))

        else:
            add_qt_menu_command(self.models_menu, "There aren't any models")
