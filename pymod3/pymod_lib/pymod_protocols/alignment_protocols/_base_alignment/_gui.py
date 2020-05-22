# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import sys

import pymod_lib.pymod_vars as pmdt

from pymol.Qt import QtWidgets, QtGui, QtCore
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_radioselect_qt)


class Alignment_window_qt(PyMod_protocol_window_qt):

    def add_middle_frame_widgets(self):
        """
        The middle frame of the window will contain a frame with widgets to choose the alignment
        mode and a frame with widgets to change the alignment algorithm parameters.
        """
        self.build_alignment_mode_frame()
        self.build_algorithm_options_frame()


    def build_alignment_mode_frame(self):
        """
        Builds a frame with some options to choose the alignment mode.
        """

        # Vbox which will store all the widgets for the alignment mode options.
        self.alignment_mode_vbox = QtWidgets.QVBoxLayout()
        self.alignment_mode_label = QtWidgets.QLabel("Alignment Mode")
        self.alignment_mode_vbox.addWidget(self.alignment_mode_label)

        self.alignment_mode_button_group = QtWidgets.QButtonGroup()
        self.build_strategy_specific_modes_frames() # Defined in child classes.

        self.middle_formlayout.addRow(self.alignment_mode_vbox)


    def build_algorithm_options_frame(self):
        """
        Options to choose the parameters of the alignment algoirthm being used.
        """
        self.build_algorithm_options_widgets()


    def build_strategy_specific_modes_frames(self):
        """
        Build components of the GUI to show the alignment options.
        """
        pass


    def build_algorithm_options_widgets(self):
        pass


    def get_alignment_mode(self):
        for radiobutton in self.alignment_mode_button_group.buttons():
            if radiobutton.isChecked():
                return radiobutton._value
        raise ValueError("No alignment mode was selected.")


###################################################################################################
# ALIGNMENT STRATEGIES.                                                                           #
###################################################################################################

class Regular_alignment_window_qt(Alignment_window_qt):
    """
    Base class to build alignment windows for regular alignments.
    """

    def build_strategy_specific_modes_frames(self):

        #----------------------------
        # Rebuild an old alignment. -
        #----------------------------

        if self.protocol.rebuild_single_alignment_choice:
            new_alignment_rb_text = "Rebuild alignment"
            new_alignment_rb_help = "Rebuild the alignment with all its sequences."
            self.new_alignment_radiobutton = QtWidgets.QRadioButton(new_alignment_rb_text)
            self.new_alignment_radiobutton.clicked.connect(self.click_on_build_new_alignment_radio)
            self.new_alignment_radiobutton._value = "rebuild-old-alignment"
            self.new_alignment_radiobutton.setChecked(True)
            self.alignment_mode_vbox.addWidget(self.new_alignment_radiobutton)
            self.alignment_mode_button_group.addButton(self.new_alignment_radiobutton)
            return None

        #------------------------------------------------------
        # Build a new alignment using the selected sequences. -
        #------------------------------------------------------

        new_alignment_rb_text = "Build a new alignment"
        new_alignment_rb_help = "Build a new alignment from scratch using the selected sequences."
        self.new_alignment_radiobutton = QtWidgets.QRadioButton(new_alignment_rb_text)
        self.new_alignment_radiobutton.clicked.connect(self.click_on_build_new_alignment_radio)
        self.new_alignment_radiobutton._value = "build-new-alignment"
        self.new_alignment_radiobutton.setChecked(True)
        self.alignment_mode_vbox.addWidget(self.new_alignment_radiobutton)
        self.alignment_mode_button_group.addButton(self.new_alignment_radiobutton)

        #--------------------
        # Alignment joiner. -
        #--------------------

        # This can be performed only if there is one selected child per cluster.
        if len(self.protocol.involved_clusters_list) > 1 and self.protocol.check_alignment_joining_selection():
            # alignment_joiner_rb_text = "Join the alignments using the selected sequences as bridges (see 'Alignment Joining')."
            self.join_alignments_radiobutton = QtWidgets.QRadioButton("Join Alignments")
            self.join_alignments_radiobutton.clicked.connect(self.click_on_alignment_joiner_radio)
            self.join_alignments_radiobutton._value = "alignment-joining"
            self.alignment_mode_vbox.addWidget(self.join_alignments_radiobutton)
            self.alignment_mode_button_group.addButton(self.join_alignments_radiobutton)


        #---------------------------
        # Keep previous alignment. -
        #---------------------------

        # Right now it can be used only when the user has selected one sequence in a cluster
        # and one sequence outside a cluster.
        if (# Only one selected cluster.
            len(self.protocol.involved_clusters_list) == 1 and
            # Only one selected sequence in the selected cluster.
            self.protocol.pymod.check_only_one_selected_child_per_cluster(self.protocol.involved_clusters_list[0]) and
            # Only one selected sequence outside any cluster.
            len(self.protocol.selected_root_sequences_list) == 1):

            self.keep_previous_alignment_radiobutton = QtWidgets.QRadioButton("Keep previous alignment")
            self.keep_previous_alignment_radiobutton.clicked.connect(self.click_on_keep_previous_alignment_radio)
            self.keep_previous_alignment_radiobutton._value = "keep-previous-alignment"
            self.alignment_mode_vbox.addWidget(self.keep_previous_alignment_radiobutton)
            self.alignment_mode_button_group.addButton(self.keep_previous_alignment_radiobutton)


    def click_on_build_new_alignment_radio(self):
        pass

    def click_on_alignment_joiner_radio(self):
        pass

    def click_on_keep_previous_alignment_radio(self):
        pass


class Profile_alignment_window_qt(Alignment_window_qt):
    """
    Base class to build windows of profile alignment protocols.
    """

    def build_strategy_specific_modes_frames(self):
        """
        Build components of the GUI to show the alignment options.
        """

        #------------------------------------------
        # Perform a profile to profile alignment. -
        #------------------------------------------

        if self.protocol.can_perform_ptp_alignment:
            # profile_profile_rb_text = "Profile to profile: perform a profile to profile alignment."
            profile_profile_rb_text = "Profile to profile"

            self.profile_to_profile_radiobutton = QtWidgets.QRadioButton(profile_profile_rb_text)
            self.profile_to_profile_radiobutton.clicked.connect(self.click_on_profile_to_profile_radio)
            self.profile_to_profile_radiobutton._value = "profile-to-profile"
            self.profile_to_profile_radiobutton.setChecked(True)
            self.alignment_mode_vbox.addWidget(self.profile_to_profile_radiobutton)
            self.alignment_mode_button_group.addButton(self.profile_to_profile_radiobutton)


        #-----------------------------------------
        # Perform sequence to profile alignment. -
        #-----------------------------------------

        sequence_profile_rb_text = None
        build_target_profile_frame = False
        # Shows a different label for the checkbutton if there is one or more clusters involved.
        if len(self.protocol.selected_clusters_list) > 1:
            # sequence_profile_rb_text = "Sequence to profile: align to a target profile the rest of the selected sequences."
            sequence_profile_rb_text = "Sequence to profile"
            build_target_profile_frame = True
        elif len(self.protocol.selected_clusters_list) == 1:
            profile_cluster_name = self.protocol.involved_clusters_list[0].my_header
            # sequence_profile_rb_text = "Sequence to profile: align the selected sequence to the target profile '%s'." % (profile_cluster_name)
            sequence_profile_rb_text = "Sequence to profile"


        # Radiobutton.
        self.sequence_to_profile_radiobutton = QtWidgets.QRadioButton(sequence_profile_rb_text)
        self.sequence_to_profile_radiobutton.clicked.connect(self.click_on_sequence_to_profile_radio)
        self.sequence_to_profile_radiobutton._value = "sequence-to-profile"
        if not self.protocol.can_perform_ptp_alignment:
            self.sequence_to_profile_radiobutton.setChecked(True)
        self.alignment_mode_vbox.addWidget(self.sequence_to_profile_radiobutton)
        self.alignment_mode_button_group.addButton(self.sequence_to_profile_radiobutton)

        # If there is more than one selected cluster, then build a frame to let the user choose
        # which is going to be the target profile.
        if build_target_profile_frame:

            # Frame with the options to choose which is going to be the target profile.
            self.target_profile_frame = QtWidgets.QFormLayout()
            self.alignment_mode_vbox.addLayout(self.target_profile_frame)

            # Label.
            self.target_alignment_label = QtWidgets.QLabel("Target profile:")
            self.target_alignment_label.setStyleSheet("margin-left: 35px")

            # Combobox.
            self.target_alignment_combobox = QtWidgets.QComboBox()
            for cluster in self.protocol.involved_clusters_list:
                self.target_alignment_combobox.addItem(cluster.my_header)
            self.target_alignment_combobox.setEditable(False)

            self.target_profile_frame.addRow(self.target_alignment_label, self.target_alignment_combobox)
            self.target_alignment_combobox.setFixedWidth(self.target_alignment_combobox.sizeHint().width())


    def click_on_profile_to_profile_radio(self):
        if hasattr(self, "target_profile_frame"):
            self.target_alignment_combobox.hide()
            self.target_alignment_label.hide()

    def click_on_sequence_to_profile_radio(self):
        if self.protocol.can_perform_ptp_alignment:
            self.target_alignment_combobox.show()
            self.target_alignment_label.show()

    def get_selected_cluster_index(self):
        return self.target_alignment_combobox.currentIndex()

    def show(self):
        Alignment_window_qt.show(self)
        # If the profile to profile option is available, the 'target_profile_frame' will be
        # hidden until the user clicks on the "sequence_to_profile_radiobutton". This is
        # performed here, because the window has to resize correctly by considering also
        # the space occupied by the widgets in the 'target_profile_frame' (which will
        # not be considered if they are already hidden).
        if self.protocol.can_perform_ptp_alignment:
            self.target_alignment_combobox.hide()
            self.target_alignment_label.hide()


###################################################################################################
# ALGORITHMS SPECIFIC CLASSES.                                                                    #
###################################################################################################

class Structural_alignment_base_window_qt:

    def build_rmsd_option(self):
        # Decide whether to compute the RMSD matrix if the structural alignment.
        self.compute_rmsd_rds = PyMod_radioselect_qt(label_text="Compute RMSD Matrix",
                                                     buttons=('Yes', 'No'))
        self.compute_rmsd_rds.setvalue("Yes")
        self.middle_formlayout.add_widget_to_align(self.compute_rmsd_rds)

    def get_compute_rmsd_option_value(self):
        return pmdt.yesno_dict[self.compute_rmsd_rds.getvalue()]
