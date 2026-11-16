# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
The two windows of the Foldseek structure search.

The options window asks what to search and where; the results window shows what
came back, one row per hit with its figures and a box to tick, as the windows of
the other searches of PyMod do.
"""

from pymol.Qt import QtWidgets

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_combobox_qt)
# The table which every search of PyMod shows its hits in. It was written here
# and moved out when the other searches came to use it.
from pymod_lib.pymod_gui.hits_table_qt import Hits_table_qt, coverage_span

from pymod_lib.pymod_protocols.similarity_searches_protocols import foldseek_client as fsc


class Foldseek_options_window_qt(PyMod_protocol_window_qt):
    """
    What to search with, and where.
    """

    def build_protocol_middle_frame(self):

        chains = self.protocol.query_chains()

        # Which chain. The server searches one chain and ignores the others, so
        # the choice is made here rather than being left to it.
        self.chain_combobox = None
        if len(chains) > 1:
            self.chain_combobox = PyMod_combobox_qt(label_text="Chain to search",
                                                    items=chains)
            self.chain_combobox.combobox.setCurrentIndex(0)
            self.middle_formlayout.add_widget_to_align(self.chain_combobox)
            note = QtWidgets.QLabel(
                "This structure has %s chains and Foldseek searches one at a"
                " time: the others take no part in this search." % len(chains))
            note.setWordWrap(True)
            self.middle_formlayout.addRow(QtWidgets.QLabel(""), note)
        elif chains:
            self.single_chain = chains[0]
            self.middle_formlayout.addRow(
                QtWidgets.QLabel("Chain to search"),
                QtWidgets.QLabel("%s, the only one" % chains[0]))

        # How the alignment is computed.
        self.mode_combobox = PyMod_combobox_qt(
            label_text="Alignment", items=[label for _, label in fsc.MODES])
        self.mode_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.mode_combobox)

        # Where to search. The list comes from the server, not from here.
        self.database_checkboxes = []
        container = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        for database in self.protocol.databases:
            checkbox = QtWidgets.QCheckBox(
                "%s%s" % (database["name"],
                          " %s" % database["version"] if database["version"] else ""))
            checkbox.setChecked(database["value"] in fsc.DEFAULT_DATABASES)
            checkbox.setToolTip("Sent as '%s'%s" % (database["value"],
                                                    ", with taxonomy"
                                                    if database["taxonomy"] else ""))
            layout.addWidget(checkbox)
            self.database_checkboxes.append((checkbox, database["value"]))

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(container)
        scroll.setMinimumHeight(170)
        self.middle_formlayout.addRow(QtWidgets.QLabel("Databases"), scroll)

        # Where the structure goes. It leaves the machine, and that is worth
        # saying plainly rather than in a manual.
        warning = QtWidgets.QLabel(
            "The chain is uploaded to the Foldseek server at %s, a free service"
            " shared by everyone who uses it. A published structure is no matter;"
            " an unpublished one leaves your computer." % fsc.DEFAULT_BASE_URL)
        warning.setWordWrap(True)
        self.middle_formlayout.addRow(QtWidgets.QLabel(""), warning)

        self.middle_formlayout.set_input_widgets_width(220)


    ###############################################################################################
    # Reading the options.                                                                        #
    ###############################################################################################

    def get_chain(self):
        if self.chain_combobox is not None:
            return self.chain_combobox.get()
        return getattr(self, "single_chain", None)

    def get_mode(self):
        index = self.mode_combobox.combobox.currentIndex()
        return fsc.MODES[index][0]

    def get_databases(self):
        return [value for checkbox, value in self.database_checkboxes
                if checkbox.isChecked()]


class Foldseek_results_window_qt(QtWidgets.QMainWindow):
    """
    The hits, one row each, with the figures which say how good they are and a
    box to tick for importing them into PyMod.
    """

    is_pymod_window = True

    # (what the header says, how wide it is kept, 0 for as wide as it needs).
    # The description is the one which needs a width of its own: it is the
    # longest thing in the table by far and would push everything else off the
    # window.
    DESCRIPTION_WIDTH = 320
    COVERAGE_WIDTH = 140

    # The colours of a selected row belong to the table which every search
    # shares, and are named here so that the window can be asked for them.
    SELECTION_COLOR = Hits_table_qt.SELECTION_COLOR
    SELECTION_TEXT = Hits_table_qt.SELECTION_TEXT

    COLUMNS = (("", 0),
               ("Hit", 0),
               ("Description", DESCRIPTION_WIDTH),
               ("Database", 0),
               ("E-value", 0),
               ("Probability", 0),
               ("Identity", 0),
               ("Query cov.", COVERAGE_WIDTH),
               ("Hit cov.", COVERAGE_WIDTH),
               ("Length", 0),
               ("Organism", 0))

    def __init__(self, parent, protocol):
        super(Foldseek_results_window_qt, self).__init__(parent)
        self.protocol = protocol

        self.setWindowTitle("Foldseek Search Results")
        self.setGeometry(240, 160, 1420, 620)

        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_vbox = QtWidgets.QVBoxLayout()

        found = sum(len(hits) for hits in self.protocol.results.values())
        databases = ", ".join(sorted(self.protocol.results))
        title = QtWidgets.QLabel(
            "%s structures resemble %s (chain %s) in %s.\nTick the ones to import:"
            " their sequences are added to a new cluster with the query, and the"
            " structures of those which have one are fetched and associated as"
            " usual."
            % (found, self.protocol.element_name(self.protocol.query_element),
               self.protocol.query_chain or "the only one", databases))
        title.setWordWrap(True)
        self.main_vbox.addWidget(title)

        # The table which every search of PyMod shows its hits in: the search
        # bar, the rows with a box to tick, and the strip which writes out the
        # description of the selected row in full.
        self.hits_table = Hits_table_qt(
            parent=self, columns=self.COLUMNS, detail_column="Description")
        self.main_vbox.addWidget(self.hits_table)

        self.fill_table()

        # The buttons.
        buttons = QtWidgets.QHBoxLayout()
        for label, command in (("Select all", self.select_all),
                               ("Select none", self.select_none),
                               ("Select the first 10", self.select_first_ten)):
            button = QtWidgets.QPushButton(label)
            button.clicked.connect(command)
            buttons.addWidget(button)

        buttons.addStretch()

        self.import_mode_label = QtWidgets.QLabel("Structure import mode:")
        buttons.addWidget(self.import_mode_label)
        self.import_mode_combobox = QtWidgets.QComboBox()
        self.import_mode_combobox.addItems(["Only the hit chain",
                                            "The whole structure"])
        buttons.addWidget(self.import_mode_combobox)

        self.show_button = QtWidgets.QPushButton("Show backbone of the selected row")
        self.show_button.clicked.connect(self.show_backbone)
        buttons.addWidget(self.show_button)

        self.import_button = QtWidgets.QPushButton("Import the ticked hits")
        self.import_button.clicked.connect(self.import_hits)
        buttons.addWidget(self.import_button)

        self.main_vbox.addLayout(buttons)
        self.central_widget.setLayout(self.main_vbox)


    def fill_table(self):
        hits = []
        for database in sorted(self.protocol.results):
            hits.extend(self.protocol.results[database])
        # Best first, whatever database it came from.
        hits.sort(key=lambda hit: (hit.evalue if hit.evalue is not None else 1e9))

        for hit in hits:
            name = hit.identifier
            if hit.structure_source is None:
                name += "  (sequence only)"
            values = [name,
                      hit.description or "",
                      hit.database,
                      _format_evalue(hit.evalue),
                      _format_number(hit.probability),
                      _format_percent(hit.identity_percent()),
                      _format_percent(hit.coverage_percent()),
                      _format_percent(hit.hit_coverage_percent()),
                      str(hit.hit_length or ""),
                      hit.taxon_name or ""]

            # The two coverages are drawn rather than written: where the aligned
            # part sits says what kind of match it is.
            spans = {"Query cov.": coverage_span(hit.query_length,
                                                 hit.query_start, hit.query_end),
                     "Hit cov.": coverage_span(hit.hit_length,
                                               hit.hit_start, hit.hit_end)}

            self.hits_table.add_row(
                values, payload=hit, spans=spans,
                detail="%s \u2014 %s" % (hit.identifier,
                                         hit.description or hit.target or ""),
                searchable=" ".join([hit.identifier or "", hit.description or "",
                                     hit.taxon_name or "", hit.database or ""]))

        self.hits_table.finish()


    # What the window was reading from itself before the table was shared. They
    # are kept so that the rest of the window, and the checks, read as they did.

    @property
    def table(self):
        return self.hits_table.table

    @property
    def rows(self):
        return [(checkbox, payload) for checkbox, payload, _, _
                in self.hits_table.rows]

    @property
    def searchbar(self):
        return self.hits_table.searchbar

    @property
    def search_count_label(self):
        return self.hits_table.search_count_label

    @property
    def description_area(self):
        return self.hits_table.detail_area

    @property
    def description_label(self):
        return self.hits_table.detail_label

    @property
    def coverage_delegate(self):
        return self.hits_table.coverage_delegate

    def column_of(self, label):
        return self.hits_table.column_of(label)

    def update_display(self):
        self.hits_table.update_display()

    def shown_rows(self):
        return [(checkbox, payload) for checkbox, payload, _, _
                in self.hits_table.shown_rows()]


    ###############################################################################################
    # The buttons.                                                                                #
    ###############################################################################################

    def select_all(self):
        self.hits_table.select_all()

    def select_none(self):
        self.hits_table.select_none()

    def select_first_ten(self):
        self.hits_table.select_first(10)

    def get_checked_hits(self):
        return self.hits_table.checked_payloads()

    def show_backbone(self):
        row = self.table.currentRow()
        if row < 0 or row >= len(self.rows):
            self.protocol.pymod.main_window.show_warning_message(
                "Foldseek", "Select a row first: its backbone is then drawn in"
                " PyMOL, without anything being downloaded.")
            return
        self.protocol.show_backbone(self.rows[row][1])

    def import_hits(self):
        hits = self.get_checked_hits()
        if not hits:
            self.protocol.pymod.main_window.show_warning_message(
                "Foldseek", "Tick the hits to import first.")
            return
        mode = ("single-chain"
                if self.import_mode_combobox.currentText() == "Only the hit chain"
                else "multiple-chains")
        self.close()
        self.protocol.import_hits(hits, import_mode=mode)


def _format_evalue(value):
    if value is None:
        return ""
    try:
        return "%.1e" % float(value)
    except (TypeError, ValueError):
        return str(value)


def _format_number(value):
    if value is None:
        return ""
    try:
        return "%.3f" % float(value)
    except (TypeError, ValueError):
        return str(value)


def _format_percent(value):
    if value is None:
        return ""
    try:
        return "%.1f%%" % float(value)
    except (TypeError, ValueError):
        return str(value)
