# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
GUI for the BoltzProt-1 protein design in PyMod.

Built on 'PyMod_protocol_window_qt' like the other protocol windows. It covers the
template free case: design a binder against the sequences selected in PyMod,
optionally against a named epitope of them.

The epitope is typed as the residue numbers the user reads in PyMod and in PyMOL,
not as the indices the API counts from zero: the protocol converts them.
"""

from pymol.Qt import QtWidgets

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_combobox_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_hbox_option_qt,
                                                          PyMod_radioselect_qt)

from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_design as bd
from pymod_lib.pymod_protocols.modeling_protocols.Boltz._ligands_gui_qt import (
    Ligands_in_pymod_mixin)


# What can be designed. Each entry is (label, modality, curated binder): a binder is
# either described by hand, giving a modality, or taken from the libraries Boltz
# maintains, giving a curated name instead.
#
# Nanobodies and antibodies go through the curated route. Describing them by hand
# would mean asking for a single chain of a length typed by the user, which is not
# what either of them is: an antibody has a heavy and a light chain. The curated
# route lets the service pick from template libraries it keeps up to date, and needs
# no length at all.
BINDER_CHOICES = [
    ("Peptide", bd.MODALITY_PEPTIDE, None),
    ("Custom protein", bd.MODALITY_CUSTOM, None),
    ("Nanobody (curated by Boltz)", None, bd.CURATED_NANOBODY),
    ("Antibody (curated by Boltz)", None, bd.CURATED_ANTIBODY),
]

# The service refuses fewer than ten designs, and charges for each of them, so the
# smallest run is offered first and the larger ones are named with their price.
DESIGN_COUNT_CHOICES = ["10", "25", "50", "100", "250", "500", "1000"]
DEFAULT_DESIGN_COUNT = "10"

DEFAULT_BINDER = "30"


class Boltz_design_window_qt(Ligands_in_pymod_mixin, PyMod_protocol_window_qt):
    """
    Window with the options of a BoltzProt-1 design run.
    """

    def additional_initialization(self):
        # The ligands of the structures open in PyMod, so that a binder can be
        # designed against the protein as it really is -- with its cofactor, its
        # metal, its drug -- and not against a bare chain.
        # Without the dictionary, which is over the network: the window opens on
        # what PyMod can see at once, and the names are filled in a moment later
        # (see 'look_up_the_names_later').
        self.detected_ligands = self.protocol.detect_ligands(ask_dictionary=False)


    def build_protocol_middle_frame(self):

        # What kind of binder.
        self.modality_combobox = PyMod_combobox_qt(
            label_text="Design a",
            items=[label for label, _, _ in BINDER_CHOICES])
        self.modality_combobox.combobox.setCurrentIndex(0)
        self.modality_combobox.combobox.currentIndexChanged.connect(
            self.on_binder_choice_changed)
        self.middle_formlayout.add_widget_to_align(self.modality_combobox)

        # Its composition. This is a small language, not a sequence, so the window
        # says what it accepts rather than leaving the user to guess.
        self.binder_entryfield = PyMod_entryfield_qt(
            label_text="Binder", value=DEFAULT_BINDER)
        self.middle_formlayout.add_widget_to_align(self.binder_entryfield)
        self.binder_hint = QtWidgets.QLabel(
            "A length (30), a range (10..20), or fixed residues\n"
            "around a designed part (MKTAYI5..10VKSHFSRQ).")
        self.binder_hint.setWordWrap(True)
        self.middle_formlayout.addRow(QtWidgets.QLabel(""), self.binder_hint)

        # Ligands to put in the target beside the selected sequences.
        self.build_ligands_row("Ligands in the target")

        # Which residues of the target to bind. Left empty, the service binds where
        # it likes, which is what PyMod did before this field existed.
        self.targets_label = QtWidgets.QLabel(self.describe_targets())
        self.targets_label.setWordWrap(True)
        self.middle_formlayout.addRow(QtWidgets.QLabel("Target chains"),
                                      self.targets_label)

        # The field and the button which fills it share a row, as the ligand of the
        # prediction window does; the hint goes underneath, across the width.
        self.epitope_subframe = PyMod_hbox_option_qt(label_text="Epitope residues")
        self.epitope_entry = QtWidgets.QLineEdit()
        self.epitope_entry.setPlaceholderText("optional")
        self.epitope_subframe.hbox.addWidget(self.epitope_entry)
        self.epitope_button = QtWidgets.QPushButton("From the PyMOL selection")
        self.epitope_button.clicked.connect(self.fill_epitope_from_pymol)
        self.epitope_subframe.hbox.addWidget(self.epitope_button)
        self.epitope_subframe.set_auto_input_widget_width()
        self.middle_formlayout.add_widget_to_align(self.epitope_subframe)

        self.epitope_hint = QtWidgets.QLabel(self.epitope_hint_text())
        self.epitope_hint.setWordWrap(True)
        self.middle_formlayout.addRow(QtWidgets.QLabel(""), self.epitope_hint)

        # The opposite of an epitope: a surface the binder must keep away from --
        # a membrane face, a crystal contact, the site of another partner.
        self.avoided_subframe = PyMod_hbox_option_qt(label_text="Avoid residues")
        self.avoided_entry = QtWidgets.QLineEdit()
        self.avoided_entry.setPlaceholderText("optional")
        self.avoided_subframe.hbox.addWidget(self.avoided_entry)
        self.avoided_button = QtWidgets.QPushButton("From the PyMOL selection")
        self.avoided_button.clicked.connect(self.fill_avoided_from_pymol)
        self.avoided_subframe.hbox.addWidget(self.avoided_button)
        self.avoided_subframe.set_auto_input_widget_width()
        self.middle_formlayout.add_widget_to_align(self.avoided_subframe,
                                                   advanced_option=True)

        # How many. Every design is charged, so the number is a main option rather
        # than an advanced one.
        self.count_combobox = PyMod_combobox_qt(
            label_text="Number of designs", items=DESIGN_COUNT_CHOICES)
        self.count_combobox.combobox.setCurrentIndex(
            DESIGN_COUNT_CHOICES.index(DEFAULT_DESIGN_COUNT))
        self.middle_formlayout.add_widget_to_align(self.count_combobox)
        self.middle_formlayout.addRow(
            QtWidgets.QLabel(""),
            QtWidgets.QLabel("The service designs at least 10 at a time and\n"
                             "charges for each: the cost is shown before running."))

        self.show_advanced_button()

        self.confirm_cost_rds = PyMod_radioselect_qt(
            label_text="Confirm the cost before running", buttons=("Yes", "No"))
        self.confirm_cost_rds.setvalue("Yes")
        self.middle_formlayout.add_widget_to_align(self.confirm_cost_rds,
                                                   advanced_option=True)

        self.middle_formlayout.set_input_widgets_width(190)

        # Put the Binder field in the state the current choice calls for.
        self.on_binder_choice_changed()


    ###############################################################################################
    # The epitope.                                                                                #
    ###############################################################################################

    def describe_targets(self):
        """Which chain of the request each selected sequence is."""

        try:
            assignment = self.protocol.target_chain_assignment()
        except Exception:
            return "-"
        return "\n".join(
            "%s → chain %s" % (getattr(element, "my_header", "?"), chain_id)
            for element, chain_id in assignment)


    def epitope_hint_text(self):
        try:
            chains = [chain for _, chain in self.protocol.target_chain_assignment()]
        except Exception:
            chains = []
        if len(chains) > 1:
            example = "%s:45-52, %s:88" % (chains[0], chains[1])
        else:
            example = "45-52, 60"
        return ("The residues the binder should bind, numbered as PyMod and PyMOL"
                " number them (%s). Leave it empty to let Boltz bind anywhere."
                % example)


    def fill_epitope_from_pymol(self):
        """Write the current selection of PyMOL into the field."""

        text, message = self.protocol.epitope_text_from_pymol()
        if text:
            self.epitope_entry.setText(text)
        if message:
            self.protocol.pymod.main_window.show_warning_message(
                "No Epitope Read", message)


    def fill_avoided_from_pymol(self):
        """The same, for the surface to keep away from."""

        text, message = self.protocol.epitope_text_from_pymol()
        if text:
            self.avoided_entry.setText(text)
        if message:
            self.protocol.pymod.main_window.show_warning_message(
                "Nothing Read", message)


    def on_binder_choice_changed(self, *args):
        """
        A curated binder has no length: the service chooses from its own templates.
        The field is therefore disabled rather than hidden, so it is clear that it
        exists and why it does not apply.
        """

        curated = self.get_curated_binder()
        self.binder_entryfield.setEnabled(curated is None)
        if curated is None:
            self.binder_hint.setText(
                "A length (30), a range (10..20), or fixed residues\n"
                "around a designed part (MKTAYI5..10VKSHFSRQ).")
        else:
            self.binder_hint.setText(
                "Not used for a curated binder: Boltz picks from the template\n"
                "libraries it maintains, so no length is given.")


    ###############################################################################################
    # Reading the options.                                                                        #
    ###############################################################################################

    def get_modality(self):
        """The modality of a binder described by hand, or None for a curated one."""
        return BINDER_CHOICES[self.modality_combobox.get_index()][1]

    def get_curated_binder(self):
        """The curated binder chosen, or None when one is described by hand."""
        return BINDER_CHOICES[self.modality_combobox.get_index()][2]

    def get_binder_specification(self):
        return self.binder_entryfield.getvalue().strip()

    def get_epitope_text(self):
        return self.epitope_entry.text().strip()

    def get_avoided_text(self):
        return self.avoided_entry.text().strip()

    def get_num_designs(self):
        return int(self.count_combobox.get())

    def get_confirm_cost(self):
        return self.confirm_cost_rds.getvalue() == "Yes"
