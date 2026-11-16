# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
The list of the ligands PyMod has, shared by the two Boltz windows.

A prediction puts them in the complex it predicts; a design puts them in the
target the binder is designed against, so that the binder meets the protein as it
really is -- with its cofactor, its metal, its drug -- rather than a bare chain.
The list itself is the same in both, so it is written once here.

What each ligand is and how it should be sent is decided by 'boltz_ligands'; this
module only shows it and remembers what was ticked.
"""

from pymol.Qt import QtWidgets, QtCore

from pymod_lib.pymod_gui.shared_gui_components_qt import small_font_style
from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_ligands as bligands


class Ligands_in_pymod_mixin:
    """
    Adds a "Ligands in PyMOL" row to a protocol window.

    The window it is mixed into must have a 'protocol' with 'detect_ligands()' and
    a 'middle_formlayout', which every PyMod protocol window has.
    """

    # The first two columns are kept to a width of their own so that the rows
    # read as columns. A name is short; where it came from can be long, and is
    # cut rather than allowed to push the SMILES off the window.
    NAME_WIDTH = 150
    WHERE_WIDTH = 150

    # How tall the list is allowed to be before it scrolls. A structure with a
    # dozen ions would otherwise push the rest of the window off the screen.
    LIST_HEIGHT = 150

    def build_ligands_row(self, label="Ligands in PyMOL"):
        """
        The row, filled with what PyMod has when the window opens: the ligands of
        the structures of the project, and the molecules which are open in PyMOL
        without being part of it -- one drawn in the builder, one opened there by
        hand, what is selected right now. Nothing is ticked to begin with.
        """

        self.ligand_checkboxes = []
        self.ligand_rows = []
        self.ligands_container = QtWidgets.QWidget()
        self.ligands_layout = QtWidgets.QVBoxLayout(self.ligands_container)
        self.ligands_layout.setContentsMargins(0, 0, 0, 0)
        self.ligands_layout.setSpacing(2)

        # After a structure has been loaded, or after a network which was not there
        # when the window opened has come back.
        self.detect_button = QtWidgets.QPushButton("Look again")
        self.detect_button.clicked.connect(self.rebuild_detected_ligands)

        self.ligands_header_label = QtWidgets.QLabel("")
        self.ligands_header_label.setStyleSheet(small_font_style)

        header = QtWidgets.QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        header.addWidget(self.ligands_header_label)
        header.addStretch()
        header.addWidget(self.detect_button)

        # The list scrolls. It holds whatever the session holds, which can be a
        # great deal, and a window which grows with it is a window which walks
        # off the screen.
        self.ligands_scroll = QtWidgets.QScrollArea()
        self.ligands_scroll.setWidgetResizable(True)
        self.ligands_scroll.setWidget(self.ligands_container)
        self.ligands_scroll.setFixedHeight(self.LIST_HEIGHT)
        self.ligands_scroll.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarAsNeeded)

        # A box of its own, with what it holds written on it. The label used to
        # stand in the left column of the form with the list beside it, which
        # left a column of empty space down the window and a label with nothing
        # under it. The box takes the whole width instead, and keeps its height:
        # the list inside it scrolls rather than pushing the rest of the window
        # off the screen.
        self.ligands_group = QtWidgets.QGroupBox(label)
        group_layout = QtWidgets.QVBoxLayout(self.ligands_group)
        group_layout.setContentsMargins(8, 4, 8, 8)
        group_layout.setSpacing(4)
        group_layout.addLayout(header)
        group_layout.addWidget(self.ligands_scroll)
        self.ligands_group.setSizePolicy(QtWidgets.QSizePolicy.Preferred,
                                         QtWidgets.QSizePolicy.Fixed)

        self.middle_formlayout.addRow(self.ligands_group)

        self.fill_detected_ligands()

        # The names of the codes are asked of the chemical component dictionary
        # over the network, which took four and a half seconds in a session of
        # five structures -- four and a half seconds of nothing at all between
        # the menu being clicked and this window appearing. What is shown above
        # was found without asking; the asking happens now, with the window
        # already on the screen, and the list is filled in again when the
        # answers arrive.
        self.look_up_the_names_later()


    def fill_detected_ligands(self):
        """
        One row per ligand, in three columns: what it is called, where it came
        from, and what will be sent -- the SMILES in a field which can be
        corrected, or the name the dictionary gives a code.

        Nothing is ticked. The columns are aligned and each is kept to its own
        width, because a list of ligands is read down a column: everything on
        one long line of its own, which is what this was, cannot be.
        """

        while self.ligands_layout.count():
            item = self.ligands_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
        self.ligand_checkboxes = []
        self.ligand_rows = []

        if not self.detected_ligands:
            self.ligands_header_label.setText("None found")
            return

        usable = [ligand for ligand in self.detected_ligands if ligand.usable]
        self.ligands_header_label.setText(
            "%s found. Tick the ones to include:" % len(self.detected_ligands)
            if usable else "None which can be sent:")

        for ligand in self.detected_ligands:
            row = QtWidgets.QWidget()
            layout = QtWidgets.QHBoxLayout(row)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(6)

            checkbox = QtWidgets.QCheckBox(ligand.label())
            checkbox.setToolTip(ligand.tooltip())
            checkbox.setMinimumWidth(self.NAME_WIDTH)
            checkbox.setMaximumWidth(self.NAME_WIDTH)
            if not ligand.usable:
                checkbox.setEnabled(False)
            layout.addWidget(checkbox)

            where = QtWidgets.QLabel(ligand.where())
            where.setStyleSheet(small_font_style)
            where.setToolTip(ligand.where())
            where.setMinimumWidth(self.WHERE_WIDTH)
            where.setMaximumWidth(self.WHERE_WIDTH)
            layout.addWidget(where)

            if ligand.usable and ligand.kind == bligands.KIND_SMILES:
                field = QtWidgets.QLineEdit(ligand.value)
                field.setToolTip("The SMILES which will be sent. It was read"
                                 " from the coordinates in PyMOL and can be"
                                 " corrected here.")
                field.setCursorPosition(0)   # the start of it, not the end
                field.textChanged.connect(
                    lambda text, l=ligand: self.on_smiles_edited(l, text))
                layout.addWidget(field, 1)
            else:
                describes = QtWidgets.QLabel(ligand.describes())
                describes.setStyleSheet(small_font_style)
                describes.setToolTip(ligand.tooltip())
                layout.addWidget(describes, 1)

            self.ligands_layout.addWidget(row)
            self.ligand_rows.append(row)
            if ligand.usable:
                self.ligand_checkboxes.append((checkbox, ligand))


    @staticmethod
    def on_smiles_edited(ligand, text):
        """What is sent follows what is written in the field."""

        ligand.edited(text)


    def rebuild_detected_ligands(self):
        """Look again, and rebuild the list from what is found now."""

        self.detect_button.setEnabled(False)
        try:
            self.detected_ligands = self.protocol.detect_ligands(look_again=True,
                                                                 forget=True)
            self.fill_detected_ligands()
        finally:
            self.detect_button.setEnabled(True)


    ###############################################################################################
    # The names, which are asked of the network after the window is on the screen.                #
    ###############################################################################################

    def look_up_the_names_later(self):
        """
        Asks the chemical component dictionary about the codes which were found,
        off the window's own thread, and fills the list in again afterwards.

        Nothing of PyMOL and nothing of Qt is touched there: the thread only
        fetches, which fills the cache of 'boltz_ligands', and the list is then
        rebuilt here, in the thread which owns the widgets, where a second scan
        costs a tenth of a second because every answer is already in hand.
        """

        codes = sorted(set(
            ligand.code for ligand in (self.detected_ligands or [])
            if bligands.is_a_real_code(ligand.code)))
        if not codes:
            return

        self.ligands_header_label.setText(
            "%s found. Looking up their names..." % len(self.detected_ligands))

        self._name_thread = _Dictionary_thread(codes)
        self._name_thread.finished_looking.connect(self.on_names_looked_up)
        self._name_thread.start()


    def on_names_looked_up(self):
        """The answers are in the cache: describe everything again, and show it."""

        try:
            self.detected_ligands = self.protocol.detect_ligands(look_again=True)
            self.fill_detected_ligands()
        except RuntimeError:
            # The window was closed while the names were being fetched.
            pass


    def get_checked_ligands(self):
        """The detected ligands the user ticked, as they were described."""

        return [ligand for checkbox, ligand in self.ligand_checkboxes
                if checkbox.isChecked()]


class _Dictionary_thread(QtCore.QThread):
    """
    Asks the chemical component dictionary about a handful of codes.

    It answers by filling the cache of 'boltz_ligands', not by returning
    anything: what the window does with the answers is the window's business,
    and doing it here would mean touching widgets from the wrong thread.
    """

    finished_looking = QtCore.pyqtSignal()

    def __init__(self, codes):
        QtCore.QThread.__init__(self)
        self.codes = list(codes)

    def run(self):
        for code in self.codes:
            try:
                bligands.ccd_entry(code)
            except Exception:
                # A code which cannot be asked about is described without the
                # dictionary, exactly as it is when there is no network.
                pass
        self.finished_looking.emit()
