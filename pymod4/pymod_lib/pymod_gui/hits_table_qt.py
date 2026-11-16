# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
The table which shows the hits of a search, and which every search of PyMod uses.

It was written for the Foldseek results and there is nothing about it which is
Foldseek's: a search returns a list of things, each with a name, some figures
saying how good it is and a box to tick, and every search of PyMod returns
exactly that. What the table adds to a list of rows is what a long list needs:

  * a search bar which hides the rows which do not carry what is typed, and says
    how many are left;
  * buttons which select what is shown rather than everything, so that a search
    followed by "Select all" cannot tick hits which are not on the screen;
  * a column which is kept to a fixed width for the long text -- a name, a
    description -- with the whole of it written under the table, where it can be
    read and scrolled instead of being hovered over;
  * a coverage drawn as a bar, the way the Foldseek website draws it: the length
    of the sequence with the aligned part filled in, which says where the
    alignment sits and not only how much of it there is.

A window builds one of these, gives it its columns and its rows, and reads back
which rows are ticked.
"""

from pymol.Qt import QtWidgets, QtCore, QtGui

from pymod_lib.pymod_gui.shared_gui_components_qt import highlight_color


###################################################################################################
# The coverage bars.                                                                              #
###################################################################################################

def coverage_span(length, start, end):
    """
    (start, end, percentage) of an alignment along a sequence, the first two as
    fractions of its length, or None when there is not enough to say.
    """

    try:
        length = float(length)
        first, last = sorted((float(start), float(end)))
    except (TypeError, ValueError):
        return None
    if not length:
        return None
    # A search counts residues from 1; a bar is drawn from 0.
    return ((first - 1) / length, last / length,
            100.0 * (last - first + 1) / length)


def span_from_fraction(start, end, fraction):
    """
    The same, for a search which reports the covered fraction instead of the
    length of the query -- which is what the BLAST and HMMER results carry.
    """

    try:
        start, end, fraction = float(start), float(end), float(fraction)
    except (TypeError, ValueError):
        return None
    if fraction <= 0:
        return None
    covered = abs(end - start) + 1
    return coverage_span(covered / fraction, start, end)


class Coverage_bar_delegate(QtWidgets.QStyledItemDelegate):
    """
    Draws a coverage as the Foldseek website draws it: a bar the length of the
    whole sequence, with the part which takes part in the alignment filled in.

    Where the filled part sits matters as much as how long it is. Two hits
    covering half the query are not the same hit if one aligns to its beginning
    and the other to its end, and a number alone does not say which is which.

    The three numbers are put on the item as (start, end, percentage), the first
    two being fractions of the sequence, and the delegate reads them from there.
    """

    ROLE = QtCore.Qt.UserRole + 1

    TRACK = QtGui.QColor("#3a3f44")
    TRACK_EDGE = QtGui.QColor("#55595e")
    FILL = QtGui.QColor(highlight_color)
    TEXT_ON_DARK = QtGui.QColor("#ffffff")

    # On a selected row. Set by the table from its own colours, so that the two
    # cannot drift apart; a stylesheet does not reach the palette, which is why
    # they are not read from there.
    SELECTED_FILL = QtGui.QColor("#cfe8f7")
    SELECTED_TRACK = QtGui.QColor("#24384a")

    BAR_HEIGHT = 9
    TEXT_WIDTH = 46

    def paint(self, painter, option, index):
        span = index.data(self.ROLE)
        if not span:
            return QtWidgets.QStyledItemDelegate.paint(self, painter, option, index)

        start, end, percent = span

        # The row keeps whatever the table would paint behind it, so a selected
        # row still looks selected.
        style_option = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(style_option, index)
        style_option.text = ""
        widget = style_option.widget
        style = widget.style() if widget is not None else QtWidgets.QApplication.style()
        style.drawControl(QtWidgets.QStyle.CE_ItemViewItem, style_option, painter, widget)

        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)

        selected = bool(option.state & QtWidgets.QStyle.State_Selected)
        fill_color = self.SELECTED_FILL if selected else self.FILL
        track_color = self.SELECTED_TRACK if selected else self.TRACK

        rect = option.rect.adjusted(4, 0, -4, 0)
        text_rect = QtCore.QRect(rect.right() - self.TEXT_WIDTH, rect.top(),
                                 self.TEXT_WIDTH, rect.height())
        track = QtCore.QRect(rect.left(),
                             rect.center().y() - self.BAR_HEIGHT // 2,
                             max(rect.width() - self.TEXT_WIDTH - 6, 1),
                             self.BAR_HEIGHT)

        painter.setPen(QtGui.QPen(self.TRACK_EDGE, 1))
        painter.setBrush(track_color)
        painter.drawRoundedRect(track, 2, 2)

        # The aligned part. It is never drawn thinner than a pixel, so that a
        # hit which covers very little of a long sequence is still visible.
        left = track.left() + int(round(track.width() * max(0.0, min(1.0, start))))
        right = track.left() + int(round(track.width() * max(0.0, min(1.0, end))))
        filled = QtCore.QRect(left, track.top(), max(right - left, 1), track.height())
        painter.setPen(QtCore.Qt.NoPen)
        painter.setBrush(fill_color)
        painter.drawRoundedRect(filled, 2, 2)

        painter.setPen(QtGui.QPen(self.TEXT_ON_DARK))
        painter.setFont(option.font)
        painter.drawText(text_rect, QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight,
                         "%.0f%%" % percent)
        painter.restore()

    def sizeHint(self, option, index):
        size = QtWidgets.QStyledItemDelegate.sizeHint(self, option, index)
        size.setHeight(max(size.height(), self.BAR_HEIGHT + 8))
        return size


###################################################################################################
# The table.                                                                                      #
###################################################################################################

class Hits_table_qt(QtWidgets.QWidget):
    """
    A search bar, a table of hits with a box to tick on each row, and a strip
    underneath which writes out in full whatever the row is too narrow to show.

    Used as:

        table = Hits_table_qt(columns=(("Hit", 0), ("Description", 320),
                                       ("E-value", 0)),
                              detail_column="Description")
        for hit in hits:
            table.add_row(["1crn_A", "CRAMBIN", "1e-30"], payload=hit)
        table.finish()
        ...
        [hit for hit in table.checked_payloads()]

    A column given a width keeps it; one given 0 is as wide as it needs to be.
    A span passed for a column is drawn as a coverage bar instead of written.
    """

    # The row which is selected. The blue of the theme is strong enough to bury
    # what is written on the row, so a quieter one is used, with the same hue.
    SELECTION_COLOR = "#33556b"
    SELECTION_TEXT = "#ffffff"

    def __init__(self, parent=None, columns=(), detail_column=None,
                 search_placeholder="a name, a description or an organism...",
                 searchable_columns=None, with_checkboxes=True):

        super(Hits_table_qt, self).__init__(parent)

        self.columns = tuple(columns)
        self.detail_column = detail_column
        self.searchable_columns = searchable_columns
        self.with_checkboxes = with_checkboxes

        # (checkbox or None, payload, the text the search bar looks through).
        self.rows = []

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # The search bar.
        search_row = QtWidgets.QHBoxLayout()
        search_row.addWidget(QtWidgets.QLabel("Search:"))
        self.searchbar = QtWidgets.QLineEdit()
        self.searchbar.setPlaceholderText(search_placeholder)
        self.searchbar.setToolTip("Only the rows which carry this text are shown,"
                                  " and the buttons which select rows work on"
                                  " what is shown.")
        self.searchbar.textChanged.connect(self.update_display)
        search_row.addWidget(self.searchbar)
        self.search_count_label = QtWidgets.QLabel("")
        search_row.addWidget(self.search_count_label)
        layout.addLayout(search_row)

        # The table.
        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(len(self.columns))
        self.table.setHorizontalHeaderLabels([label for label, _ in self.columns])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QTableWidget.SelectRows)
        self.table.setSortingEnabled(False)
        # One line per row: a long name is cut short with an ellipsis and read in
        # full underneath, which keeps the rows a uniform height and the list
        # short enough to run an eye down.
        self.table.setWordWrap(False)
        self.table.setTextElideMode(QtCore.Qt.ElideRight)
        self.table.itemSelectionChanged.connect(self.show_detail_of_selection)
        self.table.setStyleSheet(
            "QTableWidget::item:selected { background-color: %s; color: %s; }"
            % (self.SELECTION_COLOR, self.SELECTION_TEXT))
        layout.addWidget(self.table)

        self.coverage_delegate = Coverage_bar_delegate(self.table)
        self.coverage_delegate.SELECTED_TRACK = QtGui.QColor(
            self.SELECTION_COLOR).darker(135)

        # The strip which writes out what a column is too narrow to show.
        self.detail_area = QtWidgets.QScrollArea()
        self.detail_area.setWidgetResizable(True)
        self.detail_area.setFixedHeight(48)
        self.detail_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.detail_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.detail_label = QtWidgets.QLabel(
            "Select a row to read the whole of it here.")
        self.detail_label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.detail_label.setContentsMargins(4, 4, 4, 4)
        self.detail_area.setWidget(self.detail_label)
        layout.addWidget(self.detail_area)
        if detail_column is None:
            self.detail_area.setVisible(False)


    ###############################################################################################
    # Filling it.                                                                                 #
    ###############################################################################################

    def column_of(self, label):
        return [name for name, _ in self.columns].index(label)

    def add_row(self, values, payload=None, spans=None, detail=None,
                searchable=None):
        """
        One row. 'values' is one string per column, the box to tick aside;
        'spans' is {column label: (start, end, percent)} for the columns which
        are drawn as a bar; 'detail' is what the strip writes out.
        """

        row = self.table.rowCount()
        self.table.insertRow(row)

        first = 0
        checkbox = None
        if self.with_checkboxes:
            checkbox = QtWidgets.QCheckBox()
            holder = QtWidgets.QWidget()
            holder_layout = QtWidgets.QHBoxLayout(holder)
            holder_layout.addWidget(checkbox)
            holder_layout.setAlignment(QtCore.Qt.AlignCenter)
            holder_layout.setContentsMargins(0, 0, 0, 0)
            self.table.setCellWidget(row, 0, holder)
            first = 1

        for column, value in enumerate(values, start=first):
            self.table.setItem(row, column, QtWidgets.QTableWidgetItem(str(value)))

        for label, span in (spans or {}).items():
            item = self.table.item(row, self.column_of(label))
            if item is not None and span is not None:
                item.setData(Coverage_bar_delegate.ROLE, span)

        if searchable is None:
            searchable = " ".join(str(value) for value in values)
        self.rows.append((checkbox, payload, searchable.lower(),
                          detail if detail is not None else searchable))
        return row


    def finish(self):
        """Called once the rows are in: sets the widths and the delegates."""

        self.table.resizeColumnsToContents()
        for column, (label, width) in enumerate(self.columns):
            if width:
                self.table.setColumnWidth(column, width)
                self.table.horizontalHeader().setSectionResizeMode(
                    column, QtWidgets.QHeaderView.Interactive)

        for column, (label, width) in enumerate(self.columns):
            if self.column_has_bars(column):
                self.table.setItemDelegateForColumn(column, self.coverage_delegate)
        self.table.resizeRowsToContents()


    def column_has_bars(self, column):
        for row in range(self.table.rowCount()):
            item = self.table.item(row, column)
            if item is not None and item.data(Coverage_bar_delegate.ROLE):
                return True
        return False


    ###############################################################################################
    # The search bar.                                                                             #
    ###############################################################################################

    def update_display(self):
        """Shows the rows which carry what was typed, and hides the others."""

        text = self.searchbar.text().strip().lower()
        shown = 0
        for row, (checkbox, payload, searchable, detail) in enumerate(self.rows):
            matches = not text or text in searchable
            self.table.setRowHidden(row, not matches)
            shown += 1 if matches else 0

        if not text:
            self.search_count_label.setText("")
        elif shown:
            self.search_count_label.setText("%s of %s" % (shown, len(self.rows)))
        else:
            self.search_count_label.setText("nothing carries '%s'"
                                            % self.searchbar.text().strip())


    ###############################################################################################
    # Reading it.                                                                                 #
    ###############################################################################################

    def shown_rows(self):
        """The rows the search bar is letting through, in the order they are drawn."""

        return [row for position, row in enumerate(self.rows)
                if not self.table.isRowHidden(position)]

    def checkboxes(self):
        return [checkbox for checkbox, _, _, _ in self.rows]

    def shown_checkboxes(self):
        return [checkbox for checkbox, _, _, _ in self.shown_rows()]

    def checked_payloads(self):
        return [payload for checkbox, payload, _, _ in self.rows
                if checkbox is not None and checkbox.isChecked()]

    def payload_of_row(self, row):
        if row < 0 or row >= len(self.rows):
            return None
        return self.rows[row][1]

    def current_payload(self):
        return self.payload_of_row(self.table.currentRow())


    # The buttons act on what is shown: with an empty search bar that is every
    # row, and with something typed in it, what was searched for. Ticking rows
    # which are not on the screen, and then importing them, is a surprise.

    def select_all(self):
        for checkbox in self.shown_checkboxes():
            if checkbox is not None:
                checkbox.setChecked(True)

    def select_none(self):
        for checkbox in self.shown_checkboxes():
            if checkbox is not None:
                checkbox.setChecked(False)

    def select_first(self, count):
        self.select_none()
        for checkbox in self.shown_checkboxes()[:count]:
            if checkbox is not None:
                checkbox.setChecked(True)


    ###############################################################################################
    # The strip underneath.                                                                       #
    ###############################################################################################

    def show_detail_of_selection(self):
        row = self.table.currentRow()
        if row < 0 or row >= len(self.rows) or self.detail_column is None:
            return
        self.detail_label.setText(self.rows[row][3])


    def make_buttons(self, first_count=10):
        """
        The three buttons every one of these windows has, wired to this table.
        Returned as widgets, for the window to put where it likes.
        """

        buttons = []
        for label, command in (("Select all", self.select_all),
                               ("Select none", self.select_none),
                               ("Select the first %s" % first_count,
                                lambda: self.select_first(first_count))):
            button = QtWidgets.QPushButton(label)
            button.clicked.connect(command)
            buttons.append(button)
        return buttons
