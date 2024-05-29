# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import csv
import warnings

from pymol.Qt import QtCore, QtWidgets, QtGui

from pymod_lib.pymod_gui.shared_gui_components_qt import asksaveasfile_qt


class TableView(QtWidgets.QTableWidget):
    """
    Custom class derived from 'QTableWidget'. See: https://doc.qt.io/qt-5/qtableview.html.
    """

    def __init__(self, data, parent, sortable=False, column_labels=None, row_labels=None, row_labels_height=25, *args):

        self.data = data
        self.column_labels = column_labels
        self.row_labels = row_labels
        self.row_labels_height = row_labels_height

        # Get a set with the length of each column.
        rows_number = len(self.data)
        columns_number = len(self.data[0])

        QtWidgets.QTableWidget.__init__(self, parent=parent,
                                        columnCount=columns_number, rowCount=rows_number,
                                        sortingEnabled=sortable, *args)

        self.set_data()
        self.resizeColumnsToContents()
        # self.resizeRowsToContents()

        self.setStyleSheet("QTableView::item {border: 0px; padding: 0px; margin: 0px;}")
        # self.itemDoubleClicked.connect(self.on_click)

        default_font = QtGui.QFont()
        default_font.setPointSize(default_font.pointSize()-1)
        self.setFont(default_font)


    def set_data(self):

        verticalHeader = self.verticalHeader()
        verticalHeader.setSectionResizeMode(QtWidgets.QHeaderView.Fixed)
        verticalHeader.setDefaultSectionSize(self.row_labels_height)

        for i, row in enumerate(self.data):
            for j, value in enumerate(row):

                # Gets the value to show in the cell.
                if value is not None:
                    # Attempts to convert the value in a float.
                    try:
                        _value = float(value)
                    except ValueError:
                        _value = value
                else:
                    _value = "-"
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    newitem = QtWidgets.QTableWidgetItem(str(_value))
                newitem.setData(QtCore.Qt.DisplayRole, _value)
                newitem.setFlags(QtCore.Qt.ItemIsEnabled)
                newitem.setTextAlignment(QtCore.Qt.AlignCenter)
                self.setItem(i, j, newitem)
        self.setHorizontalHeaderLabels(self.column_labels)
        if self.row_labels != None:
            self.setVerticalHeaderLabels(self.row_labels)


    def on_click(self):
        pass


class QtTableWindow(QtWidgets.QMainWindow):
    """
    Window containing the 'TableView' widget.
    """

    is_pymod_window = True

    def __init__(self, parent, title, data, sortable=False, column_labels=None, row_labels=None, row_labels_height=25, width=800, height=400):
        super(QtTableWindow, self).__init__(parent)
        self.setWindowTitle(title)
        self.setGeometry(50, 50, width, height)

        self.table = TableView(data=data, parent=None, sortable=sortable, column_labels=column_labels, row_labels=row_labels, row_labels_height=row_labels_height)
        self.setCentralWidget(self.table)
        self.table.show()

        save_to_action = QtWidgets.QAction('Save to File', self)
        save_to_action.triggered.connect(lambda a=None: self.save_to_event())

        menubar = self.menuBar()
        file_menu = menubar.addMenu('File')
        file_menu.addAction(save_to_action)


    def save_to_event(self):
        """
        Saves the table data to a .csv file.
        """

        # Let the user select the filepath.
        filepath = asksaveasfile_qt("Save CSV file", name_filter="*.csv")

        if not filepath:
            return None

        try:
            # Writes a .csv file on that path.
            with open(filepath, 'w') as csv_fh:

                writer = csv.writer(csv_fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)

                if self.table.row_labels is not None:
                    writer.writerow([" "] + self.table.column_labels)
                else:
                    writer.writerow(self.table.column_labels)
                for row_idx, row in enumerate(self.table.data):
                    if self.table.row_labels is not None:
                        writer.writerow([self.table.row_labels[row_idx]] + [str(v) for v in row])
                    else:
                        writer.writerow([str(v) for v in row])

        except Exception as e:
            print("- WARNING: could not write a csv file: %s" % str(e))
