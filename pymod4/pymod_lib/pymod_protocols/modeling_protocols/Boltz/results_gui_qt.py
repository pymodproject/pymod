# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Window reporting the confidence and the binding figures of a Boltz-2 prediction.

A plain table: the figures are numbers, and what a user needs is to read them next
to what they mean. The descriptions come from 'boltz_results', so this window knows
nothing about the API.
"""

import os

from pymol.Qt import QtCore, QtGui, QtWidgets

from pymod_lib.pymod_protocols.modeling_protocols.Boltz import boltz_results as bres


class Boltz_results_window_qt(QtWidgets.QMainWindow):
    """
    Shows the figures of a finished prediction, and where its files were written.
    """

    is_pymod_window = True

    def __init__(self, parent, job, model_name, output_directory, sample=None,
                 title="Boltz-2: confidence of the prediction", *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle(title)
        self.job = job
        # The model being shown. A run produces several, each with its own figures,
        # so the window reports the one it was opened for rather than always the
        # best. None means the best one.
        self.sample = sample if sample is not None else bres.best_sample(job)
        self.output_directory = output_directory

        central = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(central)

        # One line of heading. Everything else the window used to spell out is
        # available without taking up room: the meaning of each figure is already a
        # column of the table, the rest lives in tooltips and in the status bar.
        samples = bres.sample_results(job)
        heading = model_name
        if len(samples) > 1:
            heading = "%s  (%s of %s models)" % (model_name,
                                                 self._position_of(job), len(samples))
        header = QtWidgets.QLabel(heading)
        header_font = header.font()
        header_font.setBold(True)
        header.setFont(header_font)
        header.setToolTip(
            "How far the service itself expects this model to be correct: its own"
            " estimate of its reliability, not a comparison with an experimental"
            " structure. Except where stated the figures run from 0 to 1 and higher"
            " is better; above about 0.7 is usually good, below 0.5 doubtful."
            + ("\n\nEach of the %s models of this run has its own figures, under"
               " Models -> Boltz-2." % len(samples) if len(samples) > 1 else ""))
        layout.addWidget(header)

        self.table = self._build_table(job)
        layout.addWidget(self.table)

        # The status bar carries what is worth having at hand but not worth a
        # paragraph: where the files went, and why some rows are dashes.
        status = self.statusBar()
        message = "Models written in %s" % os.path.basename(output_directory)
        if not bres.has_interfaces(job):
            message = ("'%s' = no interface to score, this prediction has a single"
                       " chain   |   " % bres.NOT_APPLICABLE) + message
        status.showMessage(message)
        status.setToolTip(
            "The models of this run were written in:\n%s\n\nFigures shown as '%s'"
            " measure how correctly the chains of a complex are placed with respect"
            " to each other, which does not apply to a prediction of a single chain."
            % (output_directory, bres.NOT_APPLICABLE))

        button_bar = QtWidgets.QHBoxLayout()
        button_bar.addStretch()
        copy_button = QtWidgets.QPushButton("Copy to clipboard")
        copy_button.setToolTip("Copies the figures below as text, to paste them in a"
                               " notebook, a message or a paper.")
        copy_button.clicked.connect(self.copy_report)
        button_bar.addWidget(copy_button)
        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        button_bar.addWidget(close_button)
        layout.addLayout(button_bar)

        self.setCentralWidget(central)

        # Sized from the rows the table actually has, since the descriptions wrap
        # and their height is not known in advance. Capped, so a long report
        # scrolls rather than growing past the screen.
        rows_height = (self.table.verticalHeader().length()
                       + self.table.horizontalHeader().height() + 4)
        self.table.setMinimumHeight(min(rows_height, 620))
        self.resize(640, min(820, rows_height + 190))


    def _position_of(self, job):
        """Which of the models of the run this one is, counting from one."""
        try:
            for index, sample in enumerate(bres.sample_results(job)):
                if sample is self.sample:
                    return index + 1
        except Exception:
            pass
        return 1


    def _build_table(self, job):
        metrics = bres.metrics_of(self.sample)
        binding = bres.binding_metrics_of(job)

        # A figure which does not apply to this prediction is shown as a dash with
        # the reason in its description, instead of the zero the service returns,
        # which reads as a very bad score.
        rows = []
        for name, label, description in bres.METRICS:
            if name not in metrics:
                continue
            rows.append((label, bres.format_value(name, metrics.get(name), job),
                         description))
        if binding:
            rows.append(None)   # a separator before the binding figures
            for name, label, description in bres.BINDING_METRICS:
                if name in binding:
                    rows.append((label,
                                 bres.format_value(name, binding.get(name), job),
                                 description))

        table = QtWidgets.QTableWidget(len(rows), 3)
        table.setHorizontalHeaderLabels(["Figure", "Value", "What it means"])
        table.horizontalHeaderItem(1).setToolTip(
            "Higher is better, except for the distance errors, which are in"
            " Angstrom and where lower is better.")
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers)
        table.setWordWrap(True)

        for index, row in enumerate(rows):
            if row is None:
                item = QtWidgets.QTableWidgetItem("Binding")
                font = item.font()
                font.setBold(True)
                item.setFont(font)
                table.setItem(index, 0, item)
                table.setItem(index, 1, QtWidgets.QTableWidgetItem(""))
                table.setItem(index, 2, QtWidgets.QTableWidgetItem(""))
                continue
            label, value, description = row
            table.setItem(index, 0, QtWidgets.QTableWidgetItem(label))
            value_item = QtWidgets.QTableWidgetItem(value)
            value_item.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            table.setItem(index, 1, value_item)
            table.setItem(index, 2, QtWidgets.QTableWidgetItem(description))

        table.horizontalHeader().setSectionResizeMode(
            2, QtWidgets.QHeaderView.Stretch)
        table.resizeColumnToContents(0)
        table.resizeColumnToContents(1)
        table.resizeRowsToContents()
        return table


    def copy_report(self):
        """
        Puts the figures on the clipboard as text, and says so: a button which does
        its work invisibly looks like a button which does nothing.
        """

        report = bres.format_report(self.job, sample=self.sample)
        QtWidgets.QApplication.clipboard().setText(report)
        QtWidgets.QMessageBox.information(
            self, "Copied",
            "The figures were copied to the clipboard as text, ready to be pasted.")
