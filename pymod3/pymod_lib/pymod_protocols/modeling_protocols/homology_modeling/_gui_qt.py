# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
GUI for MODELLER in PyMod.
"""

import multiprocessing

from pymol.Qt import QtWidgets, QtCore

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_QFormLayout,
                                                          PyMod_entryfield_qt, PyMod_entry_qt,
                                                          PyMod_radioselect_qt, PyMod_spinbox_entry_qt,
                                                          active_entry_style, inactive_entry_style,
                                                          options_title_style,
                                                          large_font_style, small_font_style,
                                                          highlight_color)

import pymod_lib.pymod_seq.seq_manipulation as pmsm


###############################################################################################
# Class for MODELLER options window.                                                          #
###############################################################################################

modeling_window_title_style = options_title_style # "%s; font-weight: bold" % options_title_style

modeling_options_sections_style = "font-weight: bold; color: #BBBBBB" # "color: %s" % highlight_color # "font-weight: bold"
modeling_options_sections_style = "color: #efefef" # "color: #ededed"

modeling_options_subsections_style = small_font_style + "; color: %s" % highlight_color
modeling_options_subsections_style = small_font_style + "; font-weight: bold"
modeling_options_subsections_style = small_font_style + "; " + modeling_options_sections_style

modeling_options_inactive_label_color = "#8c8c8c"

# modeling_window_explanation


class Modeling_window_qt(QtWidgets.QMainWindow):
    """
    A class to represent the 'Homology Modeling Window' of PyMod.
    """

    is_pymod_window = True
    # options_frame_grid_options = {"padx": 10, "pady": 10, "sticky": "w"}
    optimization_level_choices = ("Approximate", "Low", "Default", "Mid", "High", "Custom")

    def __init__(self, parent, protocol):

        super(Modeling_window_qt, self).__init__(parent)
        self.protocol = protocol
        self.initUI()


    def initUI(self):

        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle("MODELLER Options")

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        ################
        # Upper frame. #
        ################

        self.upper_frame_title = QtWidgets.QLabel("Here you can modify options for MODELLER")
        self.main_vbox.addWidget(self.upper_frame_title)


        #################
        # Middle frame. #
        #################

        # Builds the middle widget, which will contain a Notebook and its tabs.
        # This will contain most of the content of the modeling window.
        self.middle_widget = QtWidgets.QWidget()
        self.main_vbox.addWidget(self.middle_widget)

        # Layout of the 'middle_widget'.
        self.middle_gridlayout = QtWidgets.QGridLayout()
        self.middle_widget.setLayout(self.middle_gridlayout)

        # Configure the 'notebook' with its tabs.
        self.notebook = QtWidgets.QTabWidget()

        # Builds the different Notebook pages.
        self.build_main_page()
        self.build_disulfides_page()
        self.build_options_page()

        self.middle_gridlayout.addWidget(self.notebook)


        #################
        # Bottom frame. #
        #################

        # This is the "Submit" button on the modellization window.
        self.main_button = QtWidgets.QPushButton("Submit")
        self.main_button.clicked.connect(lambda a=None: self.protocol.launch_modelization())
        self.main_vbox.addWidget(self.main_button)
        self.main_button.setFixedWidth(self.main_button.sizeHint().width())


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


    ###########################################################################
    # Main tab (template or loop selection).                                  #
    ###########################################################################

    def build_main_page(self):
        """
        Add and configure the 'Main' page and tab fo the modeling window. In the
        homology modeling mode, it will contain options for choosing templates
        and some restraints option. In the loop refinement mode, it will contain
        options for choosing the loop ranges to refine.
        """

        # Main page widget.
        self.main_page = QtWidgets.QWidget()
        self.notebook.addTab(self.main_page, "Main")

        # Main page widget layout.
        self.main_page_gridlayout = QtWidgets.QGridLayout()
        self.main_page.setLayout(self.main_page_gridlayout)

        # Scroll area for the main page.
        self.main_page_scroll = QtWidgets.QScrollArea()

        # Build the main page interior and set properties of the scrollbar.
        self.main_page_interior = QtWidgets.QWidget()
        self.main_page_scroll.setWidgetResizable(True)
        self.main_page_scroll.setWidget(self.main_page_interior)
        self.main_page_gridlayout.addWidget(self.main_page_scroll)

        # self.main_page_interior_layout = QtWidgets.QVBoxLayout()
        self.main_page_interior_layout = QtWidgets.QFormLayout()
        self.main_page_interior.setLayout(self.main_page_interior_layout)

        self.build_modeling_protocol_main_page()


    def build_modeling_protocol_main_page(self):
        """
        Starts to insert content in the "Main" page.
        """

        # If the user choose to build a multiple chain model, it displays an additional option to
        # let the user choose his/her "template complex".
        if self.protocol.multiple_chain_mode:

            # Add widgets for the "template complex" selection.
            self.template_complex_selection_label = QtWidgets.QLabel("Template Complex selection")
            self.template_complex_selection_label.setStyleSheet(modeling_window_title_style)
            self.main_page_interior_layout.addRow(self.template_complex_selection_label)

            # The user can choose the "template complex" with some Radiobuttons.
            self.template_complex_rad_group = QtWidgets.QButtonGroup()

            # Display some information to explain what is a "template complex".
            information = (
            "Select the PDB file containing the complex on which you would like to base the building\n"+
            "of your multiple chain model. The relative orientation in space and the interfaces of\n"+
            "your model's chains will be based on the architecture of the Template Complex.")

            # self.template_complex_message = Label(self.template_complex_selection_frame, text=information, **shared_gui_components.modeling_window_explanation)
            # self.template_complex_message.grid(row=1, column=0, sticky = "w")

            for (tc_i, tc) in enumerate(self.protocol.available_template_complex_list):
                tcb = QtWidgets.QRadioButton(self.protocol.all_full_structure_fn_dict[tc])
                # tcb.setStyleSheet(shared_gui_components.modeling_window_rb_big)
                self.main_page_interior_layout.addRow(tcb)
                self.template_complex_rad_group.addButton(tcb)
                # Initialize by default with the first PDB in the list.
                if tc_i == 0:
                    tcb.setChecked(True)


        # Builds a frame for each modeling_cluster.
        for (mc_i, modeling_cluster) in enumerate(self.protocol.modeling_clusters_list):

            # Add a spacer to separate the sections for each modeling cluster.
            if self.protocol.multiple_chain_mode:
                spacer_frame = QtWidgets.QFrame()
                spacer_frame.setFrameShape(QtWidgets.QFrame.HLine)
                spacer_frame.setFrameShadow(QtWidgets.QFrame.Sunken)
                self.main_page_interior_layout.addRow(spacer_frame)


            #----------------------------------------------------
            # This part should contain also other options like: -
            #     - secondary structure assignment to the model -
            #     - others...                                   -
            #----------------------------------------------------

            show_symmetry_restraints_option = False

            modeling_option_label = QtWidgets.QLabel("Modeling options for target: %s" % (modeling_cluster.target_name))
            modeling_option_label.setStyleSheet(modeling_window_title_style)
            self.main_page_interior_layout.addRow(modeling_option_label)

            # Multiple chain options.
            if self.protocol.multiple_chain_mode:

                # Use symmetry restraints option.
                if modeling_cluster.symmetry_restraints_id != None:
                    show_symmetry_restraints_option = True
                    symmetry_checkbox, symmetry_info = self.build_symmetry_restraints_option(modeling_cluster)


            if any((show_symmetry_restraints_option, )): # This might include other flags in future releases.
                additional_options_label = QtWidgets.QLabel("Restraints options")
                additional_options_label.setStyleSheet(modeling_options_sections_style)
                self.main_page_interior_layout.addRow(additional_options_label)
                if show_symmetry_restraints_option:
                    self.main_page_interior_layout.addRow(symmetry_checkbox, symmetry_info)


            #-------------------------------------------------------------------
            # Builds a frame for each structure aligned to the target sequence -
            # of the current modeling cluster.                                 -
            #-------------------------------------------------------------------

            modeling_option_label = QtWidgets.QLabel("Template selection")
            modeling_option_label.setStyleSheet(modeling_options_sections_style)
            self.main_page_interior_layout.addRow(modeling_option_label)

            # Builds a frame for selecting all templates.
            if len(modeling_cluster.suitable_templates_list) > 1:
                control_all_templates_frame = QtWidgets.QFrame()
                control_all_templates_frame_layout = QtWidgets.QFormLayout()
                control_all_templates_frame.setLayout(control_all_templates_frame_layout)

                control_all_templates_label = QtWidgets.QLabel("Quick template selection")
                control_all_templates_label.setStyleSheet(modeling_options_subsections_style)
                control_all_templates_frame_layout.addRow(control_all_templates_label)

                select_all_templates_button = QtWidgets.QPushButton("Select all templates")
                select_all_templates_button.setStyleSheet(small_font_style)
                select_all_templates_button.setFixedWidth(select_all_templates_button.sizeHint().width()+20)
                select_all_templates_button.clicked.connect(lambda a=None, i=mc_i: self.select_all_templates(i))
                select_no_templates_button = QtWidgets.QPushButton("Deselect all templates")
                select_no_templates_button.setStyleSheet(small_font_style)
                select_no_templates_button.setFixedWidth(select_no_templates_button.sizeHint().width()+20)
                select_no_templates_button.clicked.connect(lambda a=None, i=mc_i: self.deselect_all_templates(i))
                control_all_templates_frame_layout.addRow(select_all_templates_button, select_no_templates_button)

                self.main_page_interior_layout.addRow(control_all_templates_frame)
                control_all_templates_frame_layout.setAlignment(QtCore.Qt.AlignLeft)
                control_all_templates_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)

            # Builds a frame for each template structure.
            modeling_cluster.structure_frame_list = []

            for (si, structure) in enumerate(modeling_cluster.suitable_templates_list):

                # This object is a PyQt one, and contains many other PyQt widgets.
                structure_frame = Structure_frame_qt(parent=None,
                                                     parent_window=self,
                                                     template_element=structure, # target_element=modeling_cluster.target
                                                     template_idx=si,
                                                     modeling_cluster=modeling_cluster,
                                                     modeling_cluster_idx=mc_i)

                self.main_page_interior_layout.addRow(structure_frame)

                # Append the current "structure_frame" to the list of the current modeling cluster.
                # Storing this object will also store the value of each Checkbox, Radiobutton and
                # entry found inside it.
                modeling_cluster.structure_frame_list.append(structure_frame)


    def build_symmetry_restraints_option(self, modeling_cluster):

        symmetry_checkbox = QtWidgets.QCheckBox("Use symmetry restraints for this chain")
        symmetry_checkbox.setStyleSheet(small_font_style)

        symmetry_information = "Show Info"
        symmetry_info = QtWidgets.QPushButton(symmetry_information)
        symmetry_info.setStyleSheet(small_font_style)
        symmetry_info.clicked.connect(lambda a=None, mc=modeling_cluster: self.show_symmetry_info(modeling_cluster))
        symmetry_info.setFixedWidth(symmetry_info.sizeHint().width()+20)

        # Store this widget for this modeling cluster, so that its value may be
        # retrieved later.
        modeling_cluster.restraints_options_dict["symmetry"] = symmetry_checkbox

        return symmetry_checkbox, symmetry_info


    def select_all_templates(self, mc_i):
        modeling_cluster = self.protocol.modeling_clusters_list[mc_i]
        for structure_frame in modeling_cluster.structure_frame_list:
            structure_frame.template_checkbox.setChecked(True)
            structure_frame.click_on_structure_checkbutton()


    def deselect_all_templates(self, mc_i):
        modeling_cluster = self.protocol.modeling_clusters_list[mc_i]
        for structure_frame in modeling_cluster.structure_frame_list:
            structure_frame.template_checkbox.setChecked(False)
            structure_frame.click_on_structure_checkbutton()


    def switch_all_hetres_checkbutton_states(self, het_radio_button_state):
        """
        Launched when the user activates/inactivates the "Include HetAtoms" in the Options page in
        the modeling window.
        """

        for mc in self.protocol.modeling_clusters_list:
            self.switch_hetres_checkbutton_states(mc, het_radio_button_state)


    def switch_hetres_checkbutton_states(self, modeling_cluster, het_radio_button_state):
        for sf in modeling_cluster.structure_frame_list:
            # Activate.
            if het_radio_button_state == 1:
                sf.hetres_radiocluster_button_state = 1
                sf.activate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.activate_het_checkbuttons()
            # Inactivate.
            elif het_radio_button_state == 0:
                sf.hetres_radiocluster_button_state = 0
                sf.inactivate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.inactivate_het_checkbuttons()
            else:
                raise KeyError("Unknown 'het_radio_button_state': %s" % het_radio_button_state)


    def show_symmetry_info(self, modeling_cluster):
        """
        Displays informations about which target sequence shares the same sequence of other targets.
        """
        mc_list = self.protocol.symmetry_restraints_groups.get_group_by_id(modeling_cluster.symmetry_restraints_id).list_of_clusters
        mc_list = [x for x in mc_list if not x is modeling_cluster]
        message1 = "The target '%s' shares the same sequence with these other targets:" % (modeling_cluster.target_name)
        seqs = "\n".join([mc.target_name for mc in mc_list])
        message2 = "so you may apply symmetry restraints for them."
        message = message1 + "\n\n" + seqs + "\n\n" + message2
        self.protocol.pymod.main_window.show_warning_message("Symmetry restraints information",
                                                             message)


    def get_template_complex_var(self):
        sel_id = None
        for rb_idx, rb in enumerate(self.template_complex_rad_group.buttons()):
            if rb.isChecked():
                sel_id = rb_idx
                break
        if sel_id is None:
            raise ValueError("No template complex was selected.")
        print(sel_id)
        return sel_id


    ###########################################################################
    # Disulfides tab.                                                         #
    ###########################################################################

    def build_disulfides_page(self):
        """
        Add the "Disulfides" page to the modeling window notebook.
        """

        #----------------------------------
        # Configures the the "Disulfides" page. -
        #----------------------------------

        # Disulfides page widget.
        self.disulfides_page = QtWidgets.QWidget()
        self.notebook.addTab(self.disulfides_page, "Disulfides")

        # Disulfides page widget layout.
        self.disulfides_page_gridlayout = QtWidgets.QGridLayout()
        self.disulfides_page.setLayout(self.disulfides_page_gridlayout)

        # Scroll area for the disulfides page.
        self.disulfides_page_scroll = QtWidgets.QScrollArea()

        # Build the main page interior and set properties of the scrollbar.
        self.disulfides_page_interior = QtWidgets.QWidget()
        self.disulfides_page_scroll.setWidgetResizable(True)
        self.disulfides_page_scroll.setWidget(self.disulfides_page_interior)
        self.disulfides_page_gridlayout.addWidget(self.disulfides_page_scroll)

        # self.disulfides_page_interior_layout = QtWidgets.QVBoxLayout()
        self.disulfides_page_interior_layout = QtWidgets.QFormLayout()
        self.disulfides_page_interior.setLayout(self.disulfides_page_interior_layout)


        #-------------------------------------------
        # Build widgets for the "Disulfides" page. -
        #-------------------------------------------

        # If at least one cluster has a target with at least two CYS residues, then build the
        # disulfide page with all its options.
        if self.protocol.check_targets_with_cys():

            # Use template disulfides.
            if self.protocol.modeling_mode == "homology_modeling":
                self.build_template_dsb_frame()

            # User defined dsb. Each target is going to have a frame to define additional dsb.
            self.build_user_defined_dsb_frame()

            # Automatically build disulfides.
            if self.protocol.modeling_mode == "homology_modeling":
                self.build_auto_dsb_frame()

        else:
            self.build_no_dsb_frame()


    def build_template_dsb_frame(self):
        """
        Builds the top frame, for the templates disulfides.
        """
        self.templates_dsb_label = QtWidgets.QLabel("Use template disulfides")
        self.templates_dsb_label.setStyleSheet(modeling_options_sections_style)
        self.disulfides_page_interior_layout.addRow(self.templates_dsb_label)

        # If there are some templates with disulfide bridges.
        if self.protocol.check_structures_with_disulfides():

            # Label for the information about the use of this feature.
            information = "Include disulfide bridges from the templates selected in the 'Main' page."
            self.template_disulfides_information = QtWidgets.QLabel(information)
            self.template_disulfides_information.setStyleSheet(small_font_style)
            self.disulfides_page_interior_layout.addRow(self.template_disulfides_information)

            # Initialize the radiobutton: if the are some structures with disulfides use this option
            # by default.
            self.use_template_dsb_rad_group = QtWidgets.QButtonGroup()
            self.use_template_dsb_rad1 = QtWidgets.QRadioButton("Yes")
            self.use_template_dsb_rad1.setChecked(True)
            # self.use_template_dsb_rad1.setStyleSheet(modeling_window_rb_big)
            self.use_template_dsb_rad1.clicked.connect(self.activate_template_dsb_frame)
            self.use_template_dsb_rad_group.addButton(self.use_template_dsb_rad1)

            self.use_template_dsb_rad2 = QtWidgets.QRadioButton("No")
            # self.use_template_dsb_rad2.setStyleSheet(modeling_window_rb_big)
            self.use_template_dsb_rad2.clicked.connect(self.inactivate_template_dsb_frame)
            self.use_template_dsb_rad_group.addButton(self.use_template_dsb_rad2)

            self.disulfides_page_interior_layout.addRow(self.use_template_dsb_rad1, self.use_template_dsb_rad2)

            # Button for displaying the list of disulfide bridges found in the templates.
            toggle_template_dsb_text = "List of templates' disulfides (white: conserved in target, gray: not conserved):"
            # self.toggle_template_frame = QtWidgets.QHBoxLayout()
            # self.toggle_template_dsb_label = QtWidgets.QLabel(toggle_template_dsb_text)
            # self.toggle_template_frame.addWidget(self.toggle_template_dsb_label)

            self.toggle_template_dsb_button = QtWidgets.QPushButton("Show template disulfides")
            self.toggle_template_dsb_button.setStyleSheet(small_font_style)
            self.toggle_template_dsb_button.clicked.connect(self.toggle_template_dsb)
            self.disulfides_page_interior_layout.addRow(self.toggle_template_dsb_button)
            self.toggle_template_dsb_button.setFixedWidth(self.toggle_template_dsb_button.sizeHint().width()+30)
            self.template_dsb_show_state = False

            # Build a frame with the list of all template disulfides (for each target).
            self.build_templates_disulfides_frame()

        # If there aren't templates with disulfide bridges.
        else:
            # Label for the information about the use of this feature.
            information = "There aren't any templates with disulfide bridges."
            self.template_disulfides_information = QtWidgets.QLabel(information)
            self.template_disulfides_information.setStyleSheet(small_font_style)
            self.disulfides_page_interior_layout.addRow(self.template_disulfides_information)


    def activate_template_dsb_frame(self):
        """
        Called when the "Yes" radiobutton of the "Use template disulfide" option is pressed.
        """
        self.toggle_template_dsb_button.show()

    def inactivate_template_dsb_frame(self):
        """
        Called when the "No" radiobutton of the "Use template disulfide" option is pressed.
        This is also called when the "Yes" radiobutton of the "Automatically build disulfides" is
        pressed.
        """
        self.toggle_template_dsb_button.hide()
        self._hide_template_dsb()

    def toggle_template_dsb(self):
        """
        Called when the "Show" button is pressed to show the list of the dsb of the templates.
        """
        if not self.template_dsb_show_state:
            self._show_template_dsb()
        else:
            self._hide_template_dsb()

    def _show_template_dsb(self):
        self.template_disulfides_frame.show()
        self.toggle_template_dsb_button.setText("Hide template disulfides")
        self.template_dsb_show_state = True

    def _hide_template_dsb(self):
        self.template_disulfides_frame.hide()
        self.toggle_template_dsb_button.setText("Show template disulfides")
        self.template_dsb_show_state = False


    def build_templates_disulfides_frame(self):
        """
        Builds the frame for displaying disulfide bridges found in the templates.
        """

        # Frame for template disulfides and its layout.
        self.template_disulfides_frame = QtWidgets.QFrame()
        # self.template_disulfides_frame.setStyleSheet(background='black', bd=1, relief = GROOVE, padx = 15, pady = 10)
        self.disulfides_page_interior_layout.addRow(self.template_disulfides_frame)
        self.template_disulfides_frame_layout = QtWidgets.QGridLayout()
        self.template_disulfides_frame.setLayout(self.template_disulfides_frame_layout)
        self.template_disulfides_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)

        # Show information for every modeling cluster which have templates with disulfides.
        row_counter = 0
        for mci, mc in enumerate([mc for mc in self.protocol.modeling_clusters_list if mc.has_structures_with_disulfides()]):

            target_label = QtWidgets.QLabel("Template disulfides for target: " + mc.target_name)
            target_label.setStyleSheet(small_font_style + "; color: %s" % highlight_color)
            self.template_disulfides_frame_layout.addWidget(target_label, row_counter, 0, 1, 2)
            row_counter += 1

            # Iterate through all the templates for this target.
            for ei, element in enumerate([x for x in mc.suitable_templates_list if x.has_disulfides()]):

                # Label with the name of the template.
                disulfides_label = QtWidgets.QLabel("Template "+ element.my_header + "     ")
                disulfides_label.setStyleSheet(modeling_options_subsections_style)
                self.template_disulfides_frame_layout.addWidget(disulfides_label, row_counter, 0)

                tar_disulfides_label = QtWidgets.QLabel("Corresponding target residues")
                tar_disulfides_label.setStyleSheet(modeling_options_subsections_style)
                self.template_disulfides_frame_layout.addWidget(tar_disulfides_label, row_counter, 1)
                row_counter += 1

                # Begins a for cycle that is going to examine all disulfides bridges of the chain.
                for dsb in element.get_disulfides():

                    # For now, display only intrachain bridges.
                    if dsb.type == "intrachain":

                        # Check if there are homologous CYS in the target according to the alignment.
                        # Take the target sequence.
                        target = mc.target.my_sequence

                        # CYS 1.
                        cys1_alignment_position = pmsm.get_residue_id_in_aligned_sequence(element.my_sequence, dsb.cys1_seq_index)
                        cys1_target_id = pmsm.get_residue_id_in_gapless_sequence(target, cys1_alignment_position)
                        if cys1_target_id is not None:
                            cys1_target_position = cys1_target_id + 1
                            cys1_is_conserved = pmsm.find_residue_conservation(element.my_sequence, target, dsb.cys1_seq_index)
                            cys1_homolog_residue = target[cys1_alignment_position] # The corresponding residue in the target.
                        else:
                            cys1_target_position = ""
                            cys1_is_conserved = False
                            cys1_homolog_residue = "gap"

                        # CYS 2.
                        cys2_alignment_position = pmsm.get_residue_id_in_aligned_sequence(element.my_sequence, dsb.cys2_seq_index)
                        cys2_target_id = pmsm.get_residue_id_in_gapless_sequence(target, cys2_alignment_position)
                        if cys2_target_id is not None:
                            cys2_target_position = cys2_target_id + 1
                            cys2_is_conserved = pmsm.find_residue_conservation(element.my_sequence, target, dsb.cys2_seq_index)
                            cys2_homolog_residue = target[cys2_alignment_position] # The corresponding residue in the target.
                        else:
                            cys2_target_position = ""
                            cys2_is_conserved = False
                            cys2_homolog_residue = "gap"


                        tem_label_text = "C%s - C%s" % (dsb.cys1_pdb_index, dsb.cys2_pdb_index)

                        # If both CYS that form the disulfide in the template are conserved in the target.
                        if cys1_is_conserved and cys2_is_conserved:
                            # Prints also if the CYS are conserved in the target according to the
                            # alignment.
                            tar_label_text = "C%s - C%s" % (cys1_target_position, cys2_target_position)
                            # tar_label_text += " (conserved)"
                            text_color = ""
                        else:
                            tar_label_text = "%s%s - %s%s" % (cys1_homolog_residue, cys1_target_position, cys2_homolog_residue, cys2_target_position)
                            # label_text += " (not conserved)"
                            text_color = "; color: %s" % modeling_options_inactive_label_color

                        # Labels for a single template disulfide bridge.
                        tem_disulfide_label = QtWidgets.QLabel(tem_label_text)
                        tem_disulfide_label.setStyleSheet(small_font_style)
                        self.template_disulfides_frame_layout.addWidget(tem_disulfide_label, row_counter, 0)
                        tar_disulfide_label = QtWidgets.QLabel(tar_label_text)
                        tar_disulfide_label.setStyleSheet(small_font_style + text_color)
                        self.template_disulfides_frame_layout.addWidget(tar_disulfide_label, row_counter, 1)
                        row_counter += 1

        # Align to the left the widgets inside the frame.
        self.template_disulfides_frame_layout.setAlignment(QtCore.Qt.AlignLeft)
        # Don't show the frame initially.
        self.template_disulfides_frame.hide()


    def build_user_defined_dsb_frame(self):
        """
        Builds the bottom frame, for the user-defined disulfides.
        """

        # Main label for this feature.
        self.user_dsb_label = QtWidgets.QLabel("Create new disulfides")
        self.user_dsb_label.setStyleSheet(modeling_options_sections_style)
        self.disulfides_page_interior_layout.addRow(self.user_dsb_label)

        # Label for the information about the use of this feature.
        information = "Define custom disulfide bridges to be included in the model. "
        # information += ("NOTE: if the S atoms of\n" +
        #                 "the two cysteines you selected are going to be located more than 2.5A apart in the\n" +
        #                 "model, MODELLER will not build the bridge." )

        self.user_disulfides_information = QtWidgets.QLabel(information)
        self.user_disulfides_information.setStyleSheet(small_font_style)
        self.disulfides_page_interior_layout.addRow(self.user_disulfides_information)

        # Radiobuttons.
        self.use_user_defined_dsb_rad_group = QtWidgets.QButtonGroup()
        self.use_user_defined_dsb_rad1 = QtWidgets.QRadioButton("Yes")
        # self.use_user_defined_dsb_rad1.setStyleSheet(modeling_window_rb_big)
        self.use_user_defined_dsb_rad1.clicked.connect(self.activate_combo_box_frame)
        self.use_user_defined_dsb_rad_group.addButton(self.use_user_defined_dsb_rad1)

        self.use_user_defined_dsb_rad2 = QtWidgets.QRadioButton("No")
        self.use_user_defined_dsb_rad2.setChecked(True)
        # self.use_user_defined_dsb_rad2.setStyleSheet(modeling_window_rb_big)
        self.use_user_defined_dsb_rad2.clicked.connect(self.inactivate_combo_box_frame)
        self.use_user_defined_dsb_rad_group.addButton(self.use_user_defined_dsb_rad2)

        self.disulfides_page_interior_layout.addRow(self.use_user_defined_dsb_rad1, self.use_user_defined_dsb_rad2)

        # Frame where comboboxes and buttons for user defined disulfides are going to be placed.
        # This is going to be gridded by the "activate_combo_box_frame()" method below.
        self.user_defined_dsb_input_frame = QtWidgets.QFrame()
        # self.user_defined_dsb_input_frame.setStyleSheet(background='black',pady = 5)
        self.user_defined_dsb_input_frame_layout = QtWidgets.QFormLayout()
        self.user_defined_dsb_input_frame.setLayout(self.user_defined_dsb_input_frame_layout)
        self.disulfides_page_interior_layout.addRow(self.user_defined_dsb_input_frame)

        # This will contain a list of User_dsb_selector objects that will store the information
        # about user defined dsb.
        self.user_dsb_selector_list = []

        # Builds a frame where are going to be gridded a series of frames (one for each
        # modeling cluster) in order to let the user define additional disulfide bridges
        # for each target.
        for (mci, mc) in enumerate(self.protocol.modeling_clusters_list):
            uds = User_dsb_selector_frame_qt(parent=None, parent_window=self,
                                             modeling_cluster=mc)
            self.user_dsb_selector_list.append(uds)
            self.user_defined_dsb_input_frame_layout.addRow(uds)

        self.user_defined_dsb_input_frame.hide()


    def activate_combo_box_frame(self):
        self.user_defined_dsb_input_frame.show()

    def inactivate_combo_box_frame(self):
        self.user_defined_dsb_input_frame.hide()


    def build_auto_dsb_frame(self):
        """
        Builds a frame to display the option to make Modeller automatically create all dsb of the
        model.
        """

        # Main label for this feature.
        self.auto_dsb_label = QtWidgets.QLabel("Automatically build disulfides")
        self.auto_dsb_label.setStyleSheet(modeling_options_sections_style)
        self.disulfides_page_interior_layout.addRow(self.auto_dsb_label)

        # Label for the information about the use of this feature.
        information = ("MODELLER will build a disulfide for every pair of cysteines if they are\n"
                       "sufficently close in the model. ")
        # information += "NOTE: by using this option you will not be able to use the two options above."

        self.auto_disulfides_information = QtWidgets.QLabel(information)
        self.auto_disulfides_information.setStyleSheet(small_font_style)
        self.disulfides_page_interior_layout.addRow(self.auto_disulfides_information)

        # Radiobuttons.
        self.auto_dsb_rad_group = QtWidgets.QButtonGroup()
        self.auto_dsb_rad1 = QtWidgets.QRadioButton("Yes")
        # self.auto_dsb_rad1.setStyleSheet(modeling_window_rb_big)
        self.auto_dsb_rad1.clicked.connect(self.activate_auto_dsb)
        self.auto_dsb_rad_group.addButton(self.auto_dsb_rad1)

        self.auto_dsb_rad2 = QtWidgets.QRadioButton("No")
        self.auto_dsb_rad2.setChecked(True)
        # self.auto_dsb_rad2.setStyleSheet(modeling_window_rb_big)
        self.auto_dsb_rad2.clicked.connect(self.inactivate_auto_dsb)
        self.auto_dsb_rad_group.addButton(self.auto_dsb_rad2)

        self.disulfides_page_interior_layout.addRow(self.auto_dsb_rad1, self.auto_dsb_rad2)


    def activate_auto_dsb(self):
        # Inactivates the "use template dsb" radiobuttons and selects the "No" radiobutton.
        if self.protocol.check_structures_with_disulfides():
            self.use_template_dsb_rad2.setChecked(True)
            self.use_template_dsb_rad1.setEnabled(False)
            self.use_template_dsb_rad2.setEnabled(False)
            self.inactivate_template_dsb_frame()

        # Inactivates the "create new dsb" radiobuttons and selects the "No" radiobutton.
        self.use_user_defined_dsb_rad2.setChecked(True)
        self.use_user_defined_dsb_rad1.setEnabled(False)
        self.use_user_defined_dsb_rad2.setEnabled(False)
        self.user_defined_dsb_input_frame.hide()


    def inactivate_auto_dsb(self):
        # Reactivates the "use template dsb" and the "create new dsb" radiobuttons.
        if self.protocol.check_structures_with_disulfides():
            self.use_template_dsb_rad1.setEnabled(True)
            self.use_template_dsb_rad2.setEnabled(True)

        self.use_user_defined_dsb_rad1.setEnabled(True)
        self.use_user_defined_dsb_rad2.setEnabled(True)


    def build_no_dsb_frame(self):
        """
        Builds a frame that is displayed if the target sequence has less than 2 cys.
        """

        self.no_dsb_label = QtWidgets.QLabel("No disulfide bridge can be built.")
        self.no_dsb_label.setStyleSheet(modeling_options_sections_style)
        self.disulfides_page_interior_layout.addRow(self.no_dsb_label)

        information = "No target has at least two CYS residues needed to form a bridge."
        self.no_disulfides_information = QtWidgets.QLabel(information)
        self.no_disulfides_information.setStyleSheet(small_font_style)
        self.disulfides_page_interior_layout.addRow(self.no_disulfides_information)


    def get_use_template_dsb_var(self):
        if self.protocol.check_targets_with_cys():
            if self.protocol.check_structures_with_disulfides():
                if self.use_template_dsb_rad1.isChecked():
                    return 1
                elif self.use_template_dsb_rad2.isChecked():
                    return 0
                else:
                    raise ValueError("No template disulfide option was selected.")
            else:
                return 0
        else:
            return 0

    def get_use_user_defined_dsb_var(self):
        if self.protocol.check_targets_with_cys():
            if self.use_user_defined_dsb_rad1.isChecked():
                return 1
            elif self.use_user_defined_dsb_rad2.isChecked():
                return 0
            else:
                raise ValueError("No custom disulfide option was selected.")
        else:
            return 0

    def get_user_dsb_list(self):
        """
        If some user-defined disulfide bridges have been built, get that information from GUI.
        """
        return [sel.user_defined_disulfide_bridges for sel in self.user_dsb_selector_list]

    def get_auto_dsb_var(self):
        if self.protocol.check_targets_with_cys():
            if self.auto_dsb_rad1.isChecked():
                return 1
            elif self.auto_dsb_rad2.isChecked():
                return 0
            else:
                raise ValueError("No auto disulfide option was selected.")
        else:
            return 0


    ###########################################################################
    # Options tab.                                                            #
    ###########################################################################

    def build_options_page(self):
        """
        Add the "Options" page on modeling window notebook.
        """

        #----------------------------------
        # Configures the the "Options" page. -
        #----------------------------------

        # Options page widget.
        self.options_page = QtWidgets.QWidget()
        self.notebook.addTab(self.options_page, "Options")

        # Options page widget layout.
        self.options_page_gridlayout = QtWidgets.QGridLayout()
        self.options_page.setLayout(self.options_page_gridlayout)

        # Scroll area for the options page.
        self.options_page_scroll = QtWidgets.QScrollArea()

        # Build the options page interior and set properties of the scrollbar.
        self.options_page_interior = QtWidgets.QWidget()
        self.options_page_scroll.setWidgetResizable(True)
        self.options_page_scroll.setWidget(self.options_page_interior)
        self.options_page_gridlayout.addWidget(self.options_page_scroll)

        self.options_page_interior_layout = PyMod_QFormLayout()
        self.options_page_interior.setLayout(self.options_page_interior_layout)

        #--------------------------------------------
        # Start to insert modeling options widgets. -
        #--------------------------------------------

        # Option to choose the number of models that Modeller has to produce.
        if self.protocol.modeling_mode == "homology_modeling":
            n_mod_default = 1
        else:
            n_mod_default = 10
        self.max_models_enf = PyMod_entryfield_qt(label_text="Models to Build",
                                                  value=str(n_mod_default),
                                                  validate={'validator': 'integer',
                                                            'min': 1,
                                                            'max': self.protocol.max_models_per_session})
        self.options_page_interior_layout.add_widget_to_align(self.max_models_enf)


        # Option to choose if Modeller is going to include HETATMs.
        if self.protocol.modeling_mode == "homology_modeling":

            self.exclude_heteroatoms_rds = PyMod_radioselect_qt(label_text='Exclude Heteroatoms',
                                                                buttons=("Yes", "No"))
            self.exclude_heteroatoms_rds.setvalue("No")
            self.options_page_interior_layout.add_widget_to_align(self.exclude_heteroatoms_rds)
            self.exclude_heteroatoms_rds.get_button_at(0).clicked.connect(lambda a=None: self.switch_all_hetres_checkbutton_states(0)) # Yes, inactivate.
            self.exclude_heteroatoms_rds.get_button_at(1).clicked.connect(lambda a=None: self.switch_all_hetres_checkbutton_states(1)) # No, activate.

        elif self.protocol.modeling_mode == "loop_refinement":

            self.exclude_heteroatoms_rds = PyMod_radioselect_qt(label_text='Include Heteroatoms',
                                                                buttons=("Het + Water", "Het", "No"))
            self.exclude_heteroatoms_rds.setvalue("No")
            self.options_page_interior_layout.add_widget_to_align(self.exclude_heteroatoms_rds)


        # Option to choose the level of optimization for Modeller.
        self.optimization_level_rds = PyMod_radioselect_qt(label_text='Optimization Level',
                                                           buttons=self.optimization_level_choices)
        self.optimization_level_rds.setvalue("Default")
        self.options_page_interior_layout.add_widget_to_align(self.optimization_level_rds)

        for button_i in range(0, 5):
            button = self.optimization_level_rds.get_button_at(button_i)
            button.clicked.connect(self.hide_optimization_level_frame)
        self.optimization_level_rds.get_button_at(5).clicked.connect(self.show_optimization_level_frame)

        optimization_level_frame_class = self.get_optimization_level_class()
        self.optimization_level_frame = optimization_level_frame_class()
        self.options_page_interior_layout.addRow(self.optimization_level_frame)
        self.optimization_level_frame.hide()


        # Option to customize the MODELLER objective function.
        self.obj_func_choices = ("Default", "altMOD (v0.2)", "Customize")
        self.obj_func_choices_vals = ("default", "altmod", "customize")
        self.obj_func_choices_dict = dict([(k, v) for (k, v) in zip(self.obj_func_choices, self.obj_func_choices_vals)])

        if self.protocol.modeling_mode == "homology_modeling":
            custom_obj_func_label = 'Objective Function'
            available_obj_func_choices = self.obj_func_choices
            default_obj_func_choice = self.obj_func_choices[0]
        elif self.protocol.modeling_mode == "loop_refinement":
            custom_obj_func_label = "Customize Objective Function"
            available_obj_func_choices = ("Yes", "No")
            default_obj_func_choice = "No"

        self.custom_obj_func_rds = PyMod_radioselect_qt(label_text=custom_obj_func_label,
                                                        buttons=available_obj_func_choices)
        self.custom_obj_func_rds.setvalue(default_obj_func_choice)
        self.options_page_interior_layout.add_widget_to_align(self.custom_obj_func_rds)

        if self.protocol.modeling_mode == "homology_modeling":
            self.custom_obj_func_rds.get_button_at(0).clicked.connect(self.hide_custom_obj_func_frame)
            self.custom_obj_func_rds.get_button_at(1).clicked.connect(self.hide_custom_obj_func_frame)
            self.custom_obj_func_rds.get_button_at(2).clicked.connect(self.show_custom_obj_func_frame)
        elif self.protocol.modeling_mode == "loop_refinement":
            self.custom_obj_func_rds.get_button_at(0).clicked.connect(self.show_custom_obj_func_frame)
            self.custom_obj_func_rds.get_button_at(1).clicked.connect(self.hide_custom_obj_func_frame)

        custom_obj_func_frame_class = self.get_custom_obj_func_frame_class()
        self.custom_obj_func_frame = custom_obj_func_frame_class()
        self.options_page_interior_layout.addRow(self.custom_obj_func_frame)
        self.custom_obj_func_frame.hide()


        # Option to choose the way to color the models.
        self.color_models_choices = ("Default", "DOPE Score")

        self.color_models_rds = PyMod_radioselect_qt(label_text='Color models by',
                                                     buttons=self.color_models_choices)
        self.color_models_rds.setvalue("Default")
        self.options_page_interior_layout.add_widget_to_align(self.color_models_rds)


        # Option to choose whether to super models to template.
        if self.protocol.modeling_mode == "homology_modeling":
            self.superpose_models_to_templates_rds = PyMod_radioselect_qt(label_text='Superpose Models to Templates',
                                                                          buttons=("Yes", "No"))
            self.superpose_models_to_templates_rds.setvalue("Yes")
            self.options_page_interior_layout.add_widget_to_align(self.superpose_models_to_templates_rds)


        # Option to choose which should be the initial conformation of loops.
        if self.protocol.modeling_mode == "loop_refinement":
            self.loop_initial_conformation_rds = PyMod_radioselect_qt(label_text='Initial Conformation for Loops',
                                                                      buttons=("Linearized", "Original"))
            self.loop_initial_conformation_rds.setvalue("Linearized")
            self.options_page_interior_layout.add_widget_to_align(self.loop_initial_conformation_rds)


        # Option to choose the random seed of MODELLER.
        self.random_seed_enf = PyMod_entryfield_qt(label_text="Random Seed (from -50000 to -2)",
                                                   value=str(self.protocol.default_modeller_random_seed),
                                                   validate={'validator': 'integer',
                                                             'min': self.protocol.min_modeller_random_seed,
                                                             'max': self.protocol.max_modeller_random_seed})
        self.options_page_interior_layout.add_widget_to_align(self.random_seed_enf)


        # Option to choose the number of parallel jobs to be used for MODELLER. When using a
        # MODELLER installed in the PyMod conda directory, subprocesses with MODELLER (and thus
        # parallele jobs) can not be used.
        max_n_jobs = multiprocessing.cpu_count()
        self.n_jobs_enf = PyMod_spinbox_entry_qt(label_text="N. of Parallel Jobs",
                                                 value=1,
                                                 # validate={'validator': 'integer', 'min': 1, 'max': max_n_jobs},
                                                 spinbox_min=1, spinbox_max=max_n_jobs
                                                 )
        self.options_page_interior_layout.add_widget_to_align(self.n_jobs_enf)


        self.options_page_interior_layout.set_input_widgets_width("auto", padding=20)


    def show_optimization_level_frame(self):
        self.optimization_level_frame.show()

    def hide_optimization_level_frame(self):
        self.optimization_level_frame.hide()

    def get_optimization_level_class(self):
        return Optimization_level_frame_qt


    def hide_custom_obj_func_frame(self):
        self.custom_obj_func_frame.hide()

    def show_custom_obj_func_frame(self):
        self.custom_obj_func_frame.show()

    def get_custom_obj_func_frame_class(self):
        return Custom_obj_func_frame_qt


    def get_custom_obj_func_option(self):
        if self.protocol.modeling_mode == "homology_modeling":
            return self.obj_func_choices_dict[self.custom_obj_func_rds.getvalue()]
        elif self.protocol.modeling_mode == "loop_refinement":
            input_val = self.custom_obj_func_rds.getvalue()
            if input_val == "Yes":
                return self.obj_func_choices_vals[2]
            elif input_val == "No":
                return self.obj_func_choices_vals[0]


###################################################################################################
# Other classes for the homology modeling window GUI.                                             #
###################################################################################################

#####################################################################
# Template selection window classes.                                #
#####################################################################

class Structure_frame_qt(QtWidgets.QFrame):
    """
    A class to construct the template selection frame and to store all their Qt widgets and
    information.
    """

    labels_width = 16
    # template_options_style = shared_gui_components.modeling_window_option_style.copy()
    # template_options_style.update({"width": labels_width})
    frames_padding = 7

    def __init__(self, parent, parent_window,
                 template_element, template_idx,
                 modeling_cluster, modeling_cluster_idx,
                 *args, **configs):

        super(Structure_frame_qt, self).__init__(parent, *args, **configs)

        self.parent_window = parent_window

        # Modeling cluster.
        self.modeling_cluster = modeling_cluster
        # This is the id of the modeling cluster containing a structure frame.
        self.mc_id = modeling_cluster_idx

        # These will contain a 'PyMod_sequence_element' type object.
        self.template_element = template_element # self.structure_pymod_element
        self.target_element = self.modeling_cluster.target
        # self.target_pymod_element = target_pymod_element

        # The int value that is passed in the for cycle in which the Structure_frame_qt objects are
        # constructed. Identifies different Structure_frame_qt objects.
        self.template_idx = template_idx # self.id

        # This is needed to check what is the state of radiobutton for using hetres. If it is on,
        # then this value should be 1 (by default it is 1 because of the default state of the
        # radiobutton), when it is off, this vaule should be 0.
        self.hetres_radiocluster_button_state = 1


        # Sets the layout of the frame.
        self.structure_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.structure_frame_layout)
        self.tot_rows = 0


        # Builds a frame for each template structure and all its options.
        self.build_use_structure_frame()
        # self.build_sequence_limit_frame()
        self.build_hetres_frame()
        self.build_water_frame()

        self.structure_frame_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.setFrameShape(QtWidgets.QFrame.StyledPanel)


    def build_use_structure_frame(self):
        """
        Builds a Frame which will contain the the checkbox for using the structure as a template.
        """

        # Initial label for the template.
        template_name = self.template_element.my_header[0:-8]
        # Avoids some problems templates with very long names.
        # For really long names it will print something like: sp_P62987...Chain:A
        if len(template_name) < 10:
            template_name = self.template_element.my_header
        else:
            template_name = self.template_element.my_header[0:10] + "..." + self.template_element.my_header[-7:]
        self.template_title_lab = QtWidgets.QLabel("Options for template: " + template_name)
        self.template_title_lab.setStyleSheet(modeling_options_subsections_style)
        self.structure_frame_layout.addWidget(self.template_title_lab, self.tot_rows, 0, 1, 2)
        self.tot_rows += 1


        # Label and checkbox for selecting a template.
        self.template_label = QtWidgets.QLabel("Use as Template: ")
        self.template_label.setStyleSheet(small_font_style)

        # Shows the identity % between the two aligned sequences.
        identity = pmsm.compute_sequence_identity(self.target_element.my_sequence, self.template_element.my_sequence)
        # Checkbox for using the structure as a template.
        checkbox_text = template_name + " (id: " + str(identity) + "%)"
        self.template_checkbox = QtWidgets.QCheckBox(checkbox_text)
        self.template_checkbox.setStyleSheet(small_font_style)
        self.template_checkbox.clicked.connect(self.click_on_structure_checkbutton)

        self.structure_frame_layout.addWidget(self.template_label, self.tot_rows, 0)
        self.structure_frame_layout.addWidget(self.template_checkbox, self.tot_rows, 1)
        self.tot_rows += 1


    def click_on_structure_checkbutton(self):
        """
        This is called when the checkbutton to use the structure as a template is pressed. If users
        want to use hetero-atoms this method will activate the hetres and water checkbuttons, that
        by default are disabled.
        """

        # This is under the influence of the state of the "use hetatm" radiobutton in the
        # options page.
        if self.hetres_radiocluster_button_state == 1:
            # The template is "activated", and so also its checkbuttons are.
            if self.get_use_as_template_var():
                self.activate_water_checkbutton()
                if self.number_of_hetres > 0:
                    self.activate_het_checkbuttons()
            # The template is "inactivated", also its checkbuttons are.
            else:
                self.inactivate_water_checkbutton()
                if self.number_of_hetres > 0:
                    self.inactivate_het_checkbuttons()


    def get_use_as_template_var(self):
        return self.template_checkbox.isChecked()


    def inactivate_het_checkbuttons(self):
        """
        This is launched when the hetres radiobutton state is changed to "NO".
        """
        self.use_all_hetres.setEnabled(False)
        self.select_single_hetres.setEnabled(False)
        self.do_not_use_hetres.setEnabled(False)
        for c in self.structure_hetres_checkbuttons:
            c.setEnabled(False)


    def activate_het_checkbuttons(self):
        """
        This is launched when the hetres radiobutton state is changed to "YES".
        """
        if self.get_use_as_template_var():
            self.use_all_hetres.setEnabled(True)
            self.select_single_hetres.setEnabled(True)
            self.do_not_use_hetres.setEnabled(True)
            for c in self.structure_hetres_checkbuttons:
                c.setEnabled(True)


    def activate_water_checkbutton(self):
        if self.template_element.has_waters():
            if self.get_use_as_template_var():
                self.water_checkbox.setEnabled(True)


    def inactivate_water_checkbutton(self):
        if self.template_element.has_waters():
            self.water_checkbox.setEnabled(False)


    def build_sequence_limit_frame(self):
        """
        Frame for the sequence limits.
        """
        pass


    def build_hetres_frame(self):
        """
        Builds a frame for the Hetero-residues selection.
        """

        # This is going to contain the checkbox states of the HETRES of the structure.
        self.hetres_options_var = 0
        # self.structure_hetres_states = []
        self.structure_hetres_checkbuttons = []
        self.structure_hetres_dict = {}


        # Label.
        self.hetres_label = QtWidgets.QLabel("Hetero Residues: ")
        self.hetres_label.setStyleSheet(small_font_style)

        self.structure_frame_layout.addWidget(self.hetres_label, self.tot_rows, 0)

        # Counts the hetres of this chain.
        self.list_of_hetres = self.template_element.get_heteroresidues()
        self.number_of_hetres = len(self.list_of_hetres)

        # Radiobuttons for hetres options.
        if self.number_of_hetres > 0:

            self.hetres_options_var = 1

            # Use all hetres checkbox.
            use_all_hetres_text = "Use all heteroatomic residues (%s)" % (self.number_of_hetres)
            self.use_all_hetres = QtWidgets.QRadioButton(use_all_hetres_text)
            self.use_all_hetres.setStyleSheet(small_font_style)
            self.use_all_hetres.clicked.connect(self.hide_select_single_hetres_frame)
            self.use_all_hetres.setChecked(True)
            self.use_all_hetres.setEnabled(False)
            self.structure_frame_layout.addWidget(self.use_all_hetres, self.tot_rows, 1)
            self.tot_rows += 1

            # Select single hetres manually.
            self.select_single_hetres = QtWidgets.QRadioButton("Select single heteroatomic residues")
            self.select_single_hetres.setStyleSheet(small_font_style)
            self.select_single_hetres.clicked.connect(self.show_select_single_hetres_frame)
            self.select_single_hetres.setEnabled(False)
            self.structure_frame_layout.addWidget(self.select_single_hetres, self.tot_rows, 1)
            self.tot_rows += 1

            # This is needed to count the "rows" used to grid HETRES checkboxes.
            for hetres in self.list_of_hetres:
                # Get the full name of the hetres.
                checkbox_text = "%s (%s) %s" % (hetres.three_letter_code, hetres.hetres_type, hetres.db_index)
                # Checkbox for each HETRES.
                hetres_checkbutton = QtWidgets.QCheckBox(checkbox_text)
                hetres_checkbutton.setStyleSheet(small_font_style + "; padding-left: 45px")
                hetres_checkbutton.setEnabled(False)

                self.structure_frame_layout.addWidget(hetres_checkbutton, self.tot_rows, 1)
                self.tot_rows += 1
                hetres_checkbutton.hide()

                # Adds the single HETRES state to a list that contains the ones of the structure.
                # self.structure_hetres_states.append(single_hetres_state)
                self.structure_hetres_checkbuttons.append(hetres_checkbutton)
                self.structure_hetres_dict.update({hetres: hetres_checkbutton})

            # Do not use any hetres.
            self.do_not_use_hetres = QtWidgets.QRadioButton("Do not use any heteroatomic residue")
            self.do_not_use_hetres.setStyleSheet(small_font_style)
            self.do_not_use_hetres.clicked.connect(self.hide_select_single_hetres_frame)
            self.do_not_use_hetres.setEnabled(False)
            self.structure_frame_layout.addWidget(self.do_not_use_hetres, self.tot_rows, 1)
            self.tot_rows += 1


        else:
            self.no_hetres_label = QtWidgets.QLabel("No heteroatomic residue found")
            self.no_hetres_label.setStyleSheet(small_font_style + "; color: %s" % modeling_options_inactive_label_color)
            self.structure_frame_layout.addWidget(self.no_hetres_label, self.tot_rows, 1)
            self.tot_rows += 1
            self.hetres_options_var = 3


    def get_hetres_options_var(self):
        if self.number_of_hetres > 0:
            if self.use_all_hetres.isChecked():
                return 1
            elif self.select_single_hetres.isChecked():
                return 2
            elif self.do_not_use_hetres.isChecked():
                return 3
            else:
                raise ValueError("No hetres radiobuttun was pressed.")
        else:
            return 3

    def show_select_single_hetres_frame(self):
        for checkbox in self.structure_hetres_checkbuttons:
            checkbox.show()

    def hide_select_single_hetres_frame(self):
        for checkbox in self.structure_hetres_checkbuttons:
            checkbox.hide()


    def build_water_frame(self):
        """
        Builds a frame for letting the user choose to include water molecules in the model.
        """

        # Label for water.
        self.water_label = QtWidgets.QLabel("Include Water: ")
        self.water_label.setStyleSheet(small_font_style)
        self.structure_frame_layout.addWidget(self.water_label, self.tot_rows, 0)

        # Checkbox for water
        if self.template_element.has_waters():
            n_water = len(self.template_element.get_waters())
            text_for_water_checkbox = "%s water molecules" % (n_water)
            self.water_checkbox = QtWidgets.QCheckBox(text_for_water_checkbox)
            self.water_checkbox.setStyleSheet(small_font_style)
            self.water_checkbox.setEnabled(False)
            self.water_checkbox.clicked.connect(lambda x=self.template_idx: self.click_on_water_checkbutton(x))
            self.structure_frame_layout.addWidget(self.water_checkbox, self.tot_rows, 1)

        else:
            self.no_water_label = QtWidgets.QLabel("This structure has no water molecules")
            self.no_water_label.setStyleSheet(small_font_style + "; color: %s" % modeling_options_inactive_label_color)
            self.structure_frame_layout.addWidget(self.no_water_label, self.tot_rows, 1)


    def click_on_water_checkbutton(self,x):
        """
        When a structure water checkbutton is pressed, this method deselects the water checkbutton of
        all the other structures, because only water from one structure can be used to build the
        model.
        """
        for sf in self.parent_window.protocol.modeling_clusters_list[self.mc_id].structure_frame_list:
            if sf.template_idx != self.template_idx and sf.template_element.has_waters():
                sf.water_checkbox.setChecked(False)


    def get_template_hetres_dict(self):
        """
        Gets a dictionary that indicates which heteroresidue to include in the models.
        """
        template_hetres_dict = {}
        for h in self.structure_hetres_dict:
            if self.get_hetres_options_var() == 1:
                template_hetres_dict.update({h: 1})
            elif self.get_hetres_options_var() == 2:
                template_hetres_dict.update({h: int(self.structure_hetres_dict[h].isChecked())})
            elif self.get_hetres_options_var() == 3:
                template_hetres_dict.update({h: 0})

        for water in self.template_element.get_waters():
            template_hetres_dict.update({water: self.get_water_state_var()})

        return template_hetres_dict

    def is_selected(self):
        return self.template_checkbox.isChecked()

    def get_water_state_var(self):
        if self.template_element.has_waters():
            return int(self.water_checkbox.isChecked())
        else:
            return 0


class User_dsb_selector_frame_qt(QtWidgets.QFrame):
    """
    Each modeling cluster will be used to build an object of this class. It will be used to let
    users define custom disulfides bridges in the model chains.
    """

    def __init__(self, parent, parent_window, modeling_cluster, *args, **configs):
        self.parent_window = parent_window
        self.protocol = self.parent_window.protocol
        self.modeling_cluster = modeling_cluster
        super(QtWidgets.QFrame, self).__init__(parent, *args, **configs)
        self.initUI()
        self.initialize_list()

        self.use_dsb_selector_frame_layout.setAlignment(QtCore.Qt.AlignLeft)
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)


    def initUI(self):

        self.use_dsb_selector_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.use_dsb_selector_frame_layout)
        label_text = ""
        if self.modeling_cluster.has_target_with_multiple_cys():
            label_text = "Select two CYS for target " + self.modeling_cluster.target_name
        else:
            label_text = "Target " + self.modeling_cluster.target_name + " does not have at least two CYS residues."
        self.modeling_cluster_custom_dsb_label = QtWidgets.QLabel(label_text)
        self.modeling_cluster_custom_dsb_label.setStyleSheet(modeling_options_subsections_style)
        self.use_dsb_selector_frame_layout.addWidget(self.modeling_cluster_custom_dsb_label, 0, 0, 1, 3)


    def initialize_list(self):
        """
        Build the initial row in the user-defined disulfide bridges frame.
        """

        self.target_list_of_cysteines = []
        # The target chain sequence.
        self.target_seq = self.modeling_cluster.target.my_sequence
        # This is going to contain User_disulfide_combo_qt objects.
        self.list_of_disulfide_combos = []
        # This list is going to contain info about disulfide bridges defined by the user through the
        # GUI. It is going to contain elements like this [[41,xx],[58,yy]] (the numbers are the
        # position of the target cysteines in both the sequence and the alignment).
        self.user_defined_disulfide_bridges = []

        # Builds an interface to let the user define additional dsb only for targets which have at
        # least two CYS residues.
        if self.modeling_cluster.has_target_with_multiple_cys():

            for (k, r) in enumerate(str(self.target_seq).replace("-", "")):
                if r == "C":
                    cys = {"position": k + 1,
                           "alignment-position": pmsm.get_residue_id_in_aligned_sequence(self.target_seq, k),
                           "state": "free"}
                    self.target_list_of_cysteines.append(cys)

            # If the target sequence has at least two cys, then creates the comboboxes.
            first = User_disulfide_combo_qt(self, self.target_list_of_cysteines)
            self.list_of_disulfide_combos.append(first)


    def add_new_user_disulfide(self):
        """
        This is called when the "Add" button to add a user-defined disulfide is pressed.
        """
        cys_1_sel = self.list_of_disulfide_combos[-1].get_cys1_val()
        cys_2_sel = self.list_of_disulfide_combos[-1].get_cys2_val()

        # Checks that both the comboboxes have been used to select a cys.
        if (cys_1_sel == "" or cys_2_sel == ""):
            txt = "You have to select two cysteines residue to define a disulfide bridge!"
            self.protocol.pymod.main_window.show_warning_message("Warning", txt)
            return None

        # Checks that the same cys has not been selected in both comboboxes.
        elif cys_1_sel == cys_2_sel:
            txt = "You cannot select the same cysteine to form a disulfide bridge!"
            self.protocol.pymod.main_window.show_warning_message("Warning", txt)
            return None

        # Checks that the selected cys are not engaged in other bridges.
        # ...


        # If the two cys are free to form a bridge, then adds the new bridge and updates the
        # frame with a new combobox row.

        # Adds the new row with comboboxes and an "Add" button.
        new_ds_combo = User_disulfide_combo_qt(self, self.target_list_of_cysteines)

        # Activates the previous row and returns the name of the 2 selected cys.
        cysteines = self.list_of_disulfide_combos[-1].activate()

        # Finishes and adds the new row.
        self.list_of_disulfide_combos.append(new_ds_combo)

        # Adds the cys pair to the self.user_defined_disulfide_bridges, which is going to be
        # used in the launch_modelization() method.
        self.user_defined_disulfide_bridges.append(cysteines)


    def remove_user_disulfide(self, udc_to_remove):
        """
        This is called when the "Remove" button is pressed.
        """
        # Deactivate and get the right bridge to remove.
        dsb_to_remove = udc_to_remove.deactivate()
        # Finishes to adds the new row.
        self.list_of_disulfide_combos.remove(udc_to_remove)
        # Also removes the bridge from the self.user_defined_disulfide_bridges.
        self.user_defined_disulfide_bridges.remove(dsb_to_remove)


class User_disulfide_combo_qt:
    """
    Class for building in the 'Disulfide' page in the modeling window a "row" with two comboboxes and
    a button to add or remove a user defined disulfide bridge to be included in the model.
    """

    # This is used in the constructor when a new combobox row is created.
    id_counter = 0

    def __init__(self, parent, cys_list):
        self.parent = parent

        # Selected have the "Add" button, unselected have the "Remove" button.
        self.selected = False

        User_disulfide_combo_qt.id_counter += 1
        self.row = User_disulfide_combo_qt.id_counter

        # The list of cysteines residues of the target sequence.
        self.cys_list = cys_list

        # The list of strings that is going to appear on the scrollable menus of the comboboxes.
        self.scrollable_cys_list = []
        for cys in self.cys_list:
            self.scrollable_cys_list.append("CYS " + str(cys["position"]))

        # Creates the first row with two comboboxes.
        self.create_combobox_row()


    combobox_padding = 30
    button_padding = 40

    def create_combobox_row(self):
        """
        Builds a row with two comboboxes an "Add" button.
        """

        # First CYS combobox.
        self.cys1_combobox = QtWidgets.QComboBox() # Select the first CYS:
        for item in self.scrollable_cys_list:
            self.cys1_combobox.addItem(item)
        self.cys1_combobox.setEditable(False)
        self.cys1_combobox.setStyleSheet(small_font_style)
        self.cys1_combobox.setFixedWidth(self.cys1_combobox.sizeHint().width() + self.combobox_padding)

        # Second CYS combobox.
        self.cys2_combobox = QtWidgets.QComboBox() # Select the second CYS:
        for item in self.scrollable_cys_list:
            self.cys2_combobox.addItem(item)
        self.cys2_combobox.setEditable(False)
        self.cys2_combobox.setStyleSheet(small_font_style)
        self.cys2_combobox.setFixedWidth(self.cys2_combobox.sizeHint().width() + self.combobox_padding)

        self.update_scrollable_cys_list()

        # "Add" button.
        self.new_disulfides_button = QtWidgets.QPushButton("Add")
        self.new_disulfides_button.setStyleSheet(small_font_style)
        self.new_disulfides_button.clicked.connect(self.press_add_button)
        self.new_disulfides_button.setFixedWidth(self.new_disulfides_button.sizeHint().width() + self.button_padding)

        self.parent.use_dsb_selector_frame_layout.addWidget(self.cys1_combobox, self.row, 0)
        self.parent.use_dsb_selector_frame_layout.addWidget(self.cys2_combobox, self.row, 1)
        self.parent.use_dsb_selector_frame_layout.addWidget(self.new_disulfides_button, self.row, 2)

        User_disulfide_combo_qt.id_counter += 1


    def select_cys(self,i):
        """
        This is launched by the combobox when some cys is selected.
        It should also be used to change the color of the cys according to their state.
        """
        pass


    def update_scrollable_cys_list(self):
        """
        Adjust the cysteine list to be displayed on the combobox with the right colors.
        """
        return None

        # # cys = {"position": seq_counter, "alignment-position":k, "state":"free"}
        # for (i,cys) in enumerate(self.cys_list):
        #     if cys["state"] == "free":
        #         self.cys1_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="black")
        #         self.cys2_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="black")
        #     elif cys["state"] == "engaged":
        #         self.cys1_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="gray")
        #         self.cys2_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="gray")
        #     # There must also be a condition used to mark cys residues engaged in disulfides
        #     # present in the templates.


    def press_add_button(self):
        """
        This is going to call a method of the 'User_dsb_selector_frame_qt' class.
        """
        self.parent.add_new_user_disulfide()


    def activate(self):
        """
        This is going to return the information about which cys have been selected when the "Add"
        button is pressed.
        """
        self.selected = True

        # Removes the "Add" button
        self.remove_widget(self.new_disulfides_button)
        self.new_disulfides_button = None

        # Removes both comboboxes, but before get their values.
        self.text1 = self.get_cys1_val()
        self.remove_widget(self.cys1_combobox)
        self.cys1_combobox = None

        self.text2 = self.get_cys2_val()
        self.remove_widget(self.cys2_combobox)
        self.cys2_combobox = None

        self.create_label_row()

        return self.text1, self.text2


    def remove_widget(self, widget):
        self.parent.use_dsb_selector_frame_layout.removeWidget(widget)
        widget.deleteLater()


    def create_label_row(self):
        """
        Build a row with two labels and a "Remove" button. Called when the "Add"
        button is pressed.
        """

        # Create dirst CYS label that tells which cys has been selected.
        self.cys1_label = QtWidgets.QLabel(self.text1)
        self.cys1_label.setStyleSheet(small_font_style)
        self.cys1_label.setAlignment(QtCore.Qt.AlignCenter)
        self.parent.use_dsb_selector_frame_layout.addWidget(self.cys1_label, self.row, 0)

        # Second CYS label.
        self.cys2_label = QtWidgets.QLabel(self.text2)
        self.cys2_label.setStyleSheet(small_font_style)
        self.cys2_label.setAlignment(QtCore.Qt.AlignCenter)
        self.parent.use_dsb_selector_frame_layout.addWidget(self.cys2_label, self.row, 1)

        # Adds the "Remove" button.
        self.remove_disulfides_button = QtWidgets.QPushButton("Remove")
        self.remove_disulfides_button.setStyleSheet(small_font_style)
        self.remove_disulfides_button.clicked.connect(self.press_remove_button)
        self.parent.use_dsb_selector_frame_layout.addWidget(self.remove_disulfides_button, self.row, 2)


    def press_remove_button(self):
        """
        This is going to call a method of the 'User_dsb_selector_frame_qt' class.
        """
        self.parent.remove_user_disulfide(self)


    def deactivate(self):
        """
        This is going to return the information about which bridge has been removed when "Remove"
        button is pressed.
        """

        self.selected = False

        self.remove_widget(self.cys1_label)
        self.cys1_label = None

        self.remove_widget(self.cys2_label)
        self.cys2_label = None

        self.remove_widget(self.remove_disulfides_button)
        self.remove_disulfides_button = None

        return self.text1, self.text2


    def get_cys1_val(self):
        return self.cys1_combobox.currentText()

    def get_cys2_val(self):
        return self.cys2_combobox.currentText()


###############################################################################
# Custom optimization level options.                                          #
###############################################################################

class Optimization_level_frame_qt(QtWidgets.QFrame):

    use_max_cg_iterations_enf = True
    use_vtfm_schedule_rds = True
    use_md_schedule_rds = True
    md_schedule_rds_default_val = "very fast"
    use_repeat_optimization_enf = True
    use_max_obj_func_value_enf = True

    # borderwidth=1
    # relief='groove'
    # pady=5

    def __init__(self, parent=None, *args, **configs):
        super(Optimization_level_frame_qt, self).__init__(parent, *args, **configs)
        self.initUI()
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)


    def initUI(self):

        self.optimization_frame_layout = QtWidgets.QFormLayout()
        # self.optimization_frame_layout.setVerticalSpacing(0)
        self.setLayout(self.optimization_frame_layout)

        self.title_label = QtWidgets.QLabel("Custom Optimization Options")
        self.title_label.setStyleSheet(modeling_options_subsections_style)
        self.optimization_frame_layout.addRow(self.title_label)


        # Max CG iterations.
        if self.use_max_cg_iterations_enf:
            self.max_cg_iterations_enf = Modeling_option_entryfield_qt(None,
                label_text="Max CG iterations",
                # labelmargin=self.enf_labelmargin,
                value="200", validate={'validator': 'integer', 'min': 1, 'max': 2000})
            self.optimization_frame_layout.addRow(self.max_cg_iterations_enf)


        # VTFM schedule radioselect.
        if self.use_vtfm_schedule_rds:
            self.vtfm_schedule_rds = Modeling_option_radioselect_qt(None,
                label_text="VTFM optimization schedule",
                buttons=("fastest", "very fast", "fast", "normal", "slow"))
            self.vtfm_schedule_rds.setvalue("normal")
            self.optimization_frame_layout.addRow(self.vtfm_schedule_rds)


        # MD schedule radioselect.
        if self.use_md_schedule_rds:
            self.md_schedule_rds = Modeling_option_radioselect_qt(None,
                label_text="MD refinement schedule",
                buttons=("none", "very fast", "fast", "slow", "very slow"))
            self.md_schedule_rds.setvalue(self.md_schedule_rds_default_val)
            self.optimization_frame_layout.addRow(self.md_schedule_rds)


        # Repeat optimization.
        if self.use_repeat_optimization_enf:
            self.repeat_optimization_enf = Modeling_option_entryfield_qt(None,
                label_text="Number of times to repeat optimization",
                value="1", validate={'validator': 'integer', 'min': 1, 'max': 10})
            self.optimization_frame_layout.addRow(self.repeat_optimization_enf)


        # Max objective function value.
        if self.use_max_obj_func_value_enf:
            self.max_obj_func_value_enf = Modeling_option_entryfield_qt(None,
                label_text="Max objective function value: 10^",
                value="7", validate={'validator': 'integer', 'min': 5, 'max': 10})
            self.optimization_frame_layout.addRow(self.max_obj_func_value_enf)


class Modeling_option_widget_qt(QtWidgets.QWidget):

    def initUI(self):
        self.setStyleSheet("margin-top: 0px; margin-bottom: 0px; padding-top: 0px; padding-bottom: 0px")


class Modeling_option_entryfield_qt(Modeling_option_widget_qt):

    def __init__(self, parent, label_text, value="", validate={},
                 use_checkbox=False, use_checkbox_command=None, active=True):

        super(Modeling_option_entryfield_qt, self).__init__(parent)
        self.initUI()

        # Layout.
        self.option_layout = QtWidgets.QFormLayout()
        # self.option_layout.setVerticalSpacing(0)
        self.setLayout(self.option_layout)

        # Widgets.
        self.use_checkbox = use_checkbox
        if not self.use_checkbox:
            self.label = QtWidgets.QLabel(label_text)
        else:
            self.checkbox = QtWidgets.QCheckBox(label_text)
            self.use_checkbox_command = use_checkbox_command
            if self.use_checkbox_command is not None:
                self.checkbox.clicked.connect(self.use_checkbox_command)
            self.label = self.checkbox
        self.label.setStyleSheet(small_font_style)

        self.entry = PyMod_entry_qt(value)
        self.entry.set_pmw_validator(validate)
        if active:
            self.entry.setStyleSheet(active_entry_style + "; " + small_font_style)
        else:
            self.entry.setEnabled(False)
            self.entry.setStyleSheet(inactive_entry_style + "; " + small_font_style)
        self.entry.setFixedWidth(60)

        self.option_layout.addRow(self.label, self.entry)

    def getvalue(self, validate=False):
        return self.entry.getvalue(validate=validate, option_name=self.label.text())

    def setvalue(self, value):
        return self.entry.setvalue(value)

    def activate_entry(self):
        self.entry.setEnabled(True)
        self.entry.setStyleSheet(active_entry_style + "; " + small_font_style)

    def inactivate_entry(self):
        self.entry.setEnabled(False)
        self.entry.setStyleSheet(inactive_entry_style + "; " + small_font_style)


class Modeling_option_radioselect_qt(Modeling_option_widget_qt):

    def __init__(self, parent, label_text, buttons=[]):

        super(Modeling_option_radioselect_qt, self).__init__(parent)
        self.initUI()

        # Layout.
        self.option_layout = QtWidgets.QGridLayout()
        # self.option_layout.setVerticalSpacing(0)
        self.setLayout(self.option_layout)

        # Widgets.
        self.label = QtWidgets.QLabel(label_text)
        self.label.setStyleSheet(small_font_style)
        self.option_layout.addWidget(self.label, 0, 0)

        self.button_group = QtWidgets.QButtonGroup()
        self.buttons_names = []
        self.buttons_dict = {}

        for button_idx, button_name in enumerate(buttons):
            radiobutton = QtWidgets.QRadioButton(button_name)
            radiobutton.setStyleSheet(small_font_style)
            self.buttons_names.append(button_name)
            self.buttons_dict[button_name] = radiobutton
            self.button_group.addButton(radiobutton)
            self.option_layout.addWidget(radiobutton, button_idx, 1)

        self.option_layout.setAlignment(QtCore.Qt.AlignLeft)

    def get_buttons(self):
        return self.button_group.buttons()

    def get_button_at(self, index):
        buttons = [b for b in self.get_buttons()]
        return buttons[index]

    def setvalue(self, value):
        self.buttons_dict[value].setChecked(True)

    def getvalue(self):
        checked_button = self.button_group.checkedButton()
        if checked_button is None:
            return None
        else:
            return checked_button.text()


###############################################################################
# Customize objective function options.                                       #
###############################################################################

class Custom_obj_func_frame_qt(QtWidgets.QFrame):

    use_loop_stat_pot_rds = False
    use_hddr_frame = True
    use_hdar_frame = True
    hdar_frame_text = "Homology-derived dihedral restraints. Weight: "
    use_charmm22_frame = True
    use_soft_sphere_frame = True
    use_dope_frame = True
    use_stat_pot_frame = False
    use_lj_frame = False
    use_gbsa_frame = False

    max_weight_val = 3.0

    initial_nb_cutoff_val = 4.0

    # borderwidth=1
    # relief='groove'
    # pady=5

    def __init__(self, parent=None, *args, **configs):
        super(Custom_obj_func_frame_qt, self).__init__(parent, *args, **configs)
        self.initUI()
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)


    def initUI(self):

        self.custom_obj_func_frame_layout = QtWidgets.QFormLayout()
        self.setLayout(self.custom_obj_func_frame_layout)

        self.title_label_1 = QtWidgets.QLabel("Customize Objective Function Terms")
        self.title_label_1.setStyleSheet(modeling_options_subsections_style)
        self.custom_obj_func_frame_layout.addRow(self.title_label_1)


        if self.use_loop_stat_pot_rds:
            self.loop_stat_pot_rds = Modeling_option_radioselect_qt(None,
                label_text="Loop modeling class:",
                buttons=("loopmodel", "DOPE loopmodel"))
            self.loop_stat_pot_rds.setvalue("loopmodel")
            self.custom_obj_func_frame_layout.addRow(self.loop_stat_pot_rds)
            self.loop_stat_pot_rds.get_button_at(0).clicked.connect(self.activate_fm_loop_class)
            self.loop_stat_pot_rds.get_button_at(1).clicked.connect(self.activate_dope_loop_class)

        if self.use_hddr_frame:
            self.hddr_frame = Modeling_option_entryfield_qt(parent=None,
                label_text="Homology-derived distance restraints. Weight: ",
                value="1.0", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val})
            self.custom_obj_func_frame_layout.addRow(self.hddr_frame)


        if self.use_hdar_frame:
            self.hdar_frame = Modeling_option_entryfield_qt(parent=None,
                label_text=self.hdar_frame_text,
                value="1.0", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val})
            self.custom_obj_func_frame_layout.addRow(self.hdar_frame)


        if self.use_charmm22_frame:
            self.charmm22_frame = Modeling_option_entryfield_qt(parent=None,
                label_text="CHARMM22 stereochemical terms. Weight: ",
                value="1.0", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val})
            self.custom_obj_func_frame_layout.addRow(self.charmm22_frame)


        if self.use_soft_sphere_frame:
            self.soft_sphere_frame = Modeling_option_entryfield_qt(parent=None,
                label_text="Sotf-sphere terms. Weight: ",
                value="1.0", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val})
            self.custom_obj_func_frame_layout.addRow(self.soft_sphere_frame)


        if self.use_dope_frame:
            self.dope_frame = Modeling_option_entryfield_qt(parent=None,
                label_text="Include DOPE in the objective function. Weight: ",
                value="0.5", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val},
                use_checkbox=True, use_checkbox_command=self.set_nb_cutoff_on_activate_dope,
                active=False)
            self.custom_obj_func_frame_layout.addRow(self.dope_frame)


        if self.use_stat_pot_frame:
            self.stat_pot_frame = Modeling_option_entryfield_qt(parent=None,
                label_text="Loop statistial potential. Weight: ",
                value="1.0", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val})
            self.custom_obj_func_frame_layout.addRow(self.stat_pot_frame)


        if self.use_lj_frame:
            self.lj_frame = Modeling_option_entryfield_qt(parent=None,
                label_text="Lennard-Jones terms. Weight: ",
                value="1.0", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val},
                active=False)
            self.custom_obj_func_frame_layout.addRow(self.lj_frame)


        if self.use_gbsa_frame:
            self.gbsa_frame = Modeling_option_entryfield_qt(parent=None,
                label_text="GBSA terms. Weight: ",
                value="1.0", validate={'validator': 'real', 'min': 0.0, 'max': self.max_weight_val},
                active=False)
            self.custom_obj_func_frame_layout.addRow(self.gbsa_frame)


        self.title_2_label = QtWidgets.QLabel("Customize Objective Function Parameters")
        self.title_2_label.setStyleSheet(modeling_options_subsections_style)
        self.custom_obj_func_frame_layout.addRow(self.title_2_label)


        self.nb_cutoff_frame = Modeling_option_entryfield_qt(parent=None,
            label_text="Non bonded cutoff (%s): " % ("\u212B"),
            value=str(self.initial_nb_cutoff_val), validate={'validator': 'real', 'min': 3.0, 'max': 14.0})
        self.custom_obj_func_frame_layout.addRow(self.nb_cutoff_frame)


    def set_nb_cutoff_on_activate_dope(self):
        """
        Automatically change the nb interactions cutoff when including DOPE terms in the objective
        function.
        """
        if self.dope_frame.checkbox.isChecked():
            self.nb_cutoff_frame.setvalue("8.0")
            self.dope_frame.activate_entry()
        else:
            self.nb_cutoff_frame.setvalue("4.0")
            self.dope_frame.inactivate_entry()


    def get_use_dope_var(self):
        if self.use_dope_frame:
            return self.dope_frame.checkbox.isChecked()
        else:
            raise ValueError("A DOPE weight selection entry was not built.")


    def activate_fm_loop_class(self):
        """
        Activates lj and gbsa, inactivates soft sphere.
        """
        self.lj_frame.inactivate_entry()

        self.gbsa_frame.inactivate_entry()

        self.soft_sphere_frame.activate_entry()
        self.soft_sphere_frame.entry.setvalue("1.0")

        self.stat_pot_frame.entry.setvalue("1.0")
        self.nb_cutoff_frame.entry.setvalue("7.0")


    def activate_dope_loop_class(self):
        """
        Inactivates lj and gbsa, activates soft sphere.
        """
        self.lj_frame.activate_entry()
        self.lj_frame.entry.setvalue("1.0")

        self.gbsa_frame.activate_entry()
        self.gbsa_frame.entry.setvalue("1.0")

        self.soft_sphere_frame.inactivate_entry()

        self.stat_pot_frame.entry.setvalue("0.6")
        self.nb_cutoff_frame.entry.setvalue("8.0")
