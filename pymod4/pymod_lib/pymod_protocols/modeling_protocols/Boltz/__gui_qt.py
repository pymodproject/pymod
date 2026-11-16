# Copyright 2025 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
GUI for the Boltz-2 structure and binding prediction in PyMod.

The window is built on 'PyMod_protocol_window_qt', like the other protocol
windows of PyMod, so that it follows the style of the plugin, gets its Submit
button and its "Show Advanced Options" mechanism for free. The previous version
was a bare QMainWindow which reimplemented part of that by hand.
"""

from pymol.Qt import QtWidgets, QtCore

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_combobox_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_hbox_option_qt)
from pymod_lib.pymod_protocols.modeling_protocols.Boltz._ligands_gui_qt import (
    Ligands_in_pymod_mixin)


# Labels of the alignment options, and the values the protocol reads back.
# The PyMod one is only offered when the selection actually belongs to an
# alignment, see 'get_msa_option_labels' in the protocol.
MSA_LABEL_SERVER = "Built by the server"
MSA_LABEL_PYMOD = "The alignment in PyMod"
MSA_LABEL_NONE = "None (single sequence)"

# Labels of the ligand kinds. "Detect" lets the protocol guess from what was
# typed, which is enough for a CCD code or an obvious SMILES string but can not
# always be right, so the choice stays visible.
LIGAND_LABEL_DETECT = "Detect"
LIGAND_LABEL_CCD = "CCD code"
LIGAND_LABEL_SMILES = "SMILES"

# How close the ligand is kept to the residues of a pocket. The guide of the API
# calls 4 to 8 A typical.
DEFAULT_POCKET_DISTANCE = "6.0"
POCKET_DISTANCE_CHOICES = ["4.0", "5.0", "6.0", "7.0", "8.0", "10.0"]

# Number of predicted structures. The service allows 1 to 10 and charges by
# groups, so 5 samples cost the same as 2: that makes 5 a good default.
SAMPLE_CHOICES = [str(n) for n in range(1, 11)]
DEFAULT_SAMPLES = "5"


class Boltz_options_window_qt(Ligands_in_pymod_mixin, PyMod_protocol_window_qt):
    """
    Window with the options of a Boltz-2 prediction.
    """

    def additional_initialization(self):
        # The ligands PyMod already has, each with the notation a request should
        # carry it in: the code where the chemical component dictionary has one,
        # a SMILES string read from the coordinates where it does not. They are
        # offered as checkboxes so a ligand on screen need not be typed.
        # Without the dictionary, which is over the network: the window opens on
        # what PyMod can see at once, and the names are filled in a moment later
        # (see 'look_up_the_names_later').
        self.detected_ligands = self.protocol.detect_ligands(ask_dictionary=False)


    def build_protocol_middle_frame(self):

        #---------------------------------------------------------------
        # Alignment used for the protein chains.                       -
        #---------------------------------------------------------------

        self.msa_combobox = PyMod_combobox_qt(
            label_text="Alignment",
            items=self.protocol.get_msa_option_labels())
        self.msa_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.msa_combobox)

        #---------------------------------------------------------------
        # Ligand.                                                      -
        #---------------------------------------------------------------

        # A ligand typed by the user, as a CCD code or a SMILES string.
        self.ligand_subframe = PyMod_hbox_option_qt(label_text="Ligand")
        self.ligand_entry = QtWidgets.QLineEdit()
        self.ligand_entry.setPlaceholderText("optional")
        self.ligand_subframe.hbox.addWidget(self.ligand_entry)
        self.ligand_type_combobox = QtWidgets.QComboBox()
        self.ligand_type_combobox.addItems([LIGAND_LABEL_DETECT,
                                            LIGAND_LABEL_CCD,
                                            LIGAND_LABEL_SMILES])
        self.ligand_subframe.hbox.addWidget(self.ligand_type_combobox)
        self.ligand_subframe.set_auto_input_widget_width()
        self.middle_formlayout.add_widget_to_align(self.ligand_subframe)

        #---------------------------------------------------------------
        # The pocket: where the ligand should sit.                     -
        #---------------------------------------------------------------

        # Without one, the prediction puts the ligand where it likes, which on a
        # protein with several clefts is often not where the chemistry says.
        self.pocket_subframe = PyMod_hbox_option_qt(label_text="Pocket residues")
        self.pocket_entry = QtWidgets.QLineEdit()
        self.pocket_entry.setPlaceholderText(self.pocket_example())
        self.pocket_entry.setToolTip(self.pocket_hint_text())
        self.pocket_subframe.hbox.addWidget(self.pocket_entry)
        self.pocket_button = QtWidgets.QPushButton("From the PyMOL selection")
        self.pocket_button.clicked.connect(self.fill_pocket_from_pymol)
        self.pocket_subframe.hbox.addWidget(self.pocket_button)
        self.pocket_subframe.set_auto_input_widget_width()
        # After the automatic width, which is narrower than the example the
        # field shows and than the list of residues which replaces it.
        self.pocket_entry.setMinimumWidth(240)
        self.middle_formlayout.add_widget_to_align(self.pocket_subframe)

        #---------------------------------------------------------------
        # Binding affinity.                                            -
        #---------------------------------------------------------------

        # Only meaningful with a ligand, and the service only predicts the
        # affinity of a ligand for a protein: protein-protein binding is not
        # available with this model.
        self.affinity_rds = PyMod_radioselect_qt(
            label_text="Predict binding affinity", buttons=("Yes", "No"))
        self.affinity_rds.setvalue("No")
        self.middle_formlayout.add_widget_to_align(self.affinity_rds)

        #---------------------------------------------------------------
        # The ligands PyMod and PyMOL already have.                    -
        #---------------------------------------------------------------

        # After the fields which every prediction fills in, and in a box of its
        # own: it is a list, and a list in the middle of a form breaks the form
        # in two.
        self.build_ligands_row()

        #---------------------------------------------------------------
        # Advanced options.                                            -
        #---------------------------------------------------------------

        self.show_advanced_button()

        self.samples_combobox = PyMod_combobox_qt(
            label_text="Number of models", items=SAMPLE_CHOICES)
        self.samples_combobox.combobox.setCurrentIndex(
            SAMPLE_CHOICES.index(DEFAULT_SAMPLES))
        self.middle_formlayout.add_widget_to_align(self.samples_combobox,
                                                   advanced_option=True)

        self.pocket_distance_combobox = PyMod_combobox_qt(
            label_text="Pocket distance (A)", items=POCKET_DISTANCE_CHOICES)
        self.pocket_distance_combobox.combobox.setCurrentIndex(
            POCKET_DISTANCE_CHOICES.index(DEFAULT_POCKET_DISTANCE))
        self.middle_formlayout.add_widget_to_align(self.pocket_distance_combobox,
                                                   advanced_option=True)

        self.pocket_force_rds = PyMod_radioselect_qt(
            label_text="Require the pocket", buttons=("Yes", "No"))
        self.pocket_force_rds.setvalue("No")
        self.middle_formlayout.add_widget_to_align(self.pocket_force_rds,
                                                   advanced_option=True)

        self.confirm_cost_rds = PyMod_radioselect_qt(
            label_text="Confirm the cost before running", buttons=("Yes", "No"))
        self.confirm_cost_rds.setvalue("Yes")
        self.middle_formlayout.add_widget_to_align(self.confirm_cost_rds,
                                                   advanced_option=True)

        self.middle_formlayout.set_input_widgets_width(190)


    ###############################################################################################
    # The ligands PyMod has.                                                                      #
    ###############################################################################################

    ###############################################################################################
    # The pocket.                                                                                 #
    ###############################################################################################

    def pocket_example(self):
        """
        How the residues are written, shown in the field itself: "45-52, 60", or
        "A:45-52, B:88" when the prediction has more than one chain.
        """

        try:
            chains = [chain for _, chain in self.protocol.polymer_chain_assignment()]
        except Exception:
            chains = []
        if len(chains) > 1:
            return "optional, e.g. %s:45-52, %s:88" % (chains[0], chains[1])
        return "optional, e.g. 45-52, 60"


    def pocket_hint_text(self):
        """
        What the field is for. It used to be written out under the field, where
        it was a paragraph of the window for a thing most predictions leave
        empty; it is a tooltip now.
        """

        return ("The residues which line the pocket, numbered as PyMod and PyMOL"
                " number them. The ligand is kept among them; leave it empty to"
                " let Boltz place it where it likes.")


    def fill_pocket_from_pymol(self):
        """Write the current selection of PyMOL into the field."""

        text, message = self.protocol.pocket_text_from_pymol()
        if text:
            self.pocket_entry.setText(text)
        if message:
            self.protocol.pymod.main_window.show_warning_message(
                "No Pocket Read", message)


    ###############################################################################################
    # Reading the options.                                                                        #
    ###############################################################################################

    def get_pocket_text(self):
        return self.pocket_entry.text().strip()

    def get_pocket_distance(self):
        return float(self.pocket_distance_combobox.get())

    def get_pocket_forced(self):
        return self.pocket_force_rds.getvalue() == "Yes"

    def get_msa_label(self):
        return self.msa_combobox.get()

    def get_typed_ligand(self):
        return self.ligand_entry.text().strip()

    def get_ligand_kind_label(self):
        return self.ligand_type_combobox.currentText()

    def get_affinity_requested(self):
        return self.affinity_rds.getvalue() == "Yes"

    def get_num_samples(self):
        # The combobox only offers valid numbers, but the protocol validates the
        # value again when it builds the request.
        return int(self.samples_combobox.get())

    def get_confirm_cost(self):
        # The advanced options may never have been opened, in which case the
        # widget still holds its default of "Yes".
        return self.confirm_cost_rds.getvalue() == "Yes"
