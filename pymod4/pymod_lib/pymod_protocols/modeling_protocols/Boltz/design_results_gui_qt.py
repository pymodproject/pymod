# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Window listing the designs of a BoltzProt-1 run.

A design run produces at least ten binders and may produce thousands, so they are
not poured into PyMod: they are listed, best first, and the user picks the ones
worth keeping. That is also why the sequences are imported as sequences and not as
structures, with the structure of a design fetched only if it is asked for.
"""

from pymol.Qt import QtCore, QtWidgets

from pymod_lib.pymod_protocols.modeling_protocols.Boltz import (
    boltz_design_results as bdres)


class Boltz_design_results_window_qt(QtWidgets.QMainWindow):
    """
    Shows the designs of a run and imports the chosen ones into PyMod.
    """

    is_pymod_window = True

    def __init__(self, parent, protocol, designs, run_name,
                 title="BoltzProt-1 designs", *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle(title)
        self.protocol = protocol
        self.designs = designs
        self.run_name = run_name

        central = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(central)

        header = QtWidgets.QLabel("BoltzProt-1 design run: %s" % run_name)
        font = header.font()
        font.setBold(True)
        header.setFont(font)
        layout.addWidget(header)

        usable = [d for d in designs if d["sequence"]]
        summary = ("%s designs, ordered by binding confidence. Check the ones to"
                   " keep and import them as sequences in PyMod." % len(designs))
        if len(usable) != len(designs):
            summary += (" %s of them have no readable sequence and can not be"
                        " imported." % (len(designs) - len(usable)))
        note = QtWidgets.QLabel(summary)
        note.setWordWrap(True)
        layout.addWidget(note)

        self.table = self._build_table(designs)
        layout.addWidget(self.table)

        # Selecting many rows by hand is tedious, so the usual choices are buttons.
        selection_bar = QtWidgets.QHBoxLayout()
        for label, count in (("Check the best", 1), ("Check the best 5", 5),
                             ("Check the best 10", 10)):
            button = QtWidgets.QPushButton(label)
            button.clicked.connect(lambda _=False, n=count: self.check_best(n))
            selection_bar.addWidget(button)
        clear_button = QtWidgets.QPushButton("Uncheck all")
        clear_button.clicked.connect(lambda: self.check_best(0))
        selection_bar.addWidget(clear_button)
        selection_bar.addStretch()
        layout.addLayout(selection_bar)

        button_bar = QtWidgets.QHBoxLayout()
        button_bar.addStretch()
        self.import_button = QtWidgets.QPushButton("Import the checked sequences")
        self.import_button.clicked.connect(self.import_checked)
        button_bar.addWidget(self.import_button)
        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        button_bar.addWidget(close_button)
        layout.addLayout(button_bar)

        self.setCentralWidget(central)
        self.resize(900, 560)

        # The best design is a sensible starting point.
        self.check_best(1)


    ###############################################################################################
    # The table.                                                                                  #
    ###############################################################################################

    def _build_table(self, designs):
        columns = ["", "Rank", "Length"] + [label for _, label, _ in bdres.DESIGN_METRICS]
        table = QtWidgets.QTableWidget(len(designs), len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers)
        table.setSelectionBehavior(QtWidgets.QTableWidget.SelectRows)

        # The description of every figure, so the header explains itself.
        for index, (_, label, description) in enumerate(bdres.DESIGN_METRICS):
            table.horizontalHeaderItem(index + 3).setToolTip(description)

        self.checkboxes = []
        for row, design in enumerate(designs):
            checkbox = QtWidgets.QTableWidgetItem()
            checkbox.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            checkbox.setCheckState(QtCore.Qt.Unchecked)
            if not design["sequence"]:
                # Nothing to import: the row is shown but can not be checked.
                checkbox.setFlags(QtCore.Qt.NoItemFlags)
            table.setItem(row, 0, checkbox)
            self.checkboxes.append(checkbox)

            table.setItem(row, 1, self._number_item(str(design["rank"])))
            sequence = design["sequence"] or ""
            table.setItem(row, 2, self._number_item(str(len(sequence)) if sequence
                                                    else "?"))

            for index, (name, _, _) in enumerate(bdres.DESIGN_METRICS):
                value = bdres.format_value(name, design["metrics"].get(name))
                table.setItem(row, index + 3, self._number_item(value))

            if design["warnings"]:
                tip = "\n".join(design["warnings"])
                for column in range(table.columnCount()):
                    item = table.item(row, column)
                    if item is not None:
                        item.setToolTip(tip)

        table.resizeColumnsToContents()
        return table


    @staticmethod
    def _number_item(text):
        item = QtWidgets.QTableWidgetItem(text)
        item.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        return item


    ###############################################################################################
    # Actions.                                                                                    #
    ###############################################################################################

    def check_best(self, count):
        """Checks the first 'count' importable designs, unchecking the rest."""
        checked = 0
        for checkbox, design in zip(self.checkboxes, self.designs):
            if design["sequence"] and checked < count:
                checkbox.setCheckState(QtCore.Qt.Checked)
                checked += 1
            elif design["sequence"]:
                checkbox.setCheckState(QtCore.Qt.Unchecked)


    def get_checked_designs(self):
        return [design for checkbox, design in zip(self.checkboxes, self.designs)
                if checkbox.checkState() == QtCore.Qt.Checked and design["sequence"]]


    def import_checked(self):
        checked = self.get_checked_designs()
        if not checked:
            QtWidgets.QMessageBox.information(
                self, "No design checked",
                "Check the designs to import first.")
            return
        self.protocol.import_designs(checked, self.run_name)
