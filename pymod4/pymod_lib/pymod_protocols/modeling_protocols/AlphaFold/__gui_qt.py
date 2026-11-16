# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

from pymol.Qt import QtWidgets, QtCore

"""
GUI for AlphaFold in PyMod.
"""

class AFDB_search_window_qt(QtWidgets.QMainWindow):
    """
    A class to represent the type of AFDB search (either 'by ID' or 'by sequence') of PyMod.
    """

    is_pymod_window = True

    def __init__(self, parent, protocol):

        super(AFDB_search_window_qt, self).__init__(parent)
        self.protocol = protocol
        self.initUI()

    def initUI(self):

        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle("AFDB search")
        ##### MODIFIED #####
        self.setGeometry(500, 400, 250, 300)
        ##### END #####

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        ################
        # Upper frame. #
        ################

        self.upper_frame_title = QtWidgets.QLabel("Choose an option below")
        self.main_vbox.addWidget(self.upper_frame_title)


        #################
        # Middle frame. #
        #################

        # Builds the middle widget.
        self.middle_widget = QtWidgets.QWidget()
        self.main_vbox.addWidget(self.middle_widget)

        # Layout of the 'middle_widget'.
        self.middle_gridlayout = QtWidgets.QGridLayout()
        self.middle_widget.setLayout(self.middle_gridlayout)

        # Configure the middle layout with main widget.
        self.middle_scroll_area = QtWidgets.QScrollArea()
        self.middle_gridlayout.addWidget(self.middle_scroll_area)

        self.build_middle_widgets()

        #################
        # Bottom frame. #
        #################

        # This is the "Submit" button on the modellization window.
        self.main_button = QtWidgets.QPushButton("Submit")
        self.main_button.clicked.connect(lambda a=None: self.protocol.launch_afdb_search())
        self.main_vbox.addWidget(self.main_button)
        self.main_button.setFixedWidth(self.main_button.sizeHint().width())


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


    def build_middle_widgets(self):

        self.middle_scroll_layout = QtWidgets.QVBoxLayout()
        self.middle_scroll_area.setLayout(self.middle_scroll_layout)

        self.search_by_label = QtWidgets.QLabel("Search By: ")
        self.sequence_id_rb = QtWidgets.QRadioButton("Sequence ID")
        self.sequence_id_rb.setChecked(True)
        self.sequence_rb = QtWidgets.QRadioButton("Sequence")
        self.sequence_id_rb.clicked.connect(self.show_widgets)
        self.sequence_rb.clicked.connect(self.hide_widgets)

        self.middle_scroll_layout.addWidget(self.search_by_label)
        self.middle_scroll_layout.addWidget(self.sequence_id_rb)
        self.middle_scroll_layout.addWidget(self.sequence_rb)

        self.al_mode_label = QtWidgets.QLabel("Alignment mode: ")
        self.al_mode_cb = QtWidgets.QComboBox()
        self.al_mode_cb.addItems(["Build new alignment", "Expand_existing_alignment"])

        self.middle_scroll_layout.addWidget(self.al_mode_label)
        self.middle_scroll_layout.addWidget(self.al_mode_cb)

        if not self.protocol.selected_sequences[0].is_child():
            self.al_mode_label.hide()
            self.al_mode_cb.hide()

        self.import_mode_label = QtWidgets.QLabel("Structure import mode: ")
        self.import_mode_cb = QtWidgets.QComboBox()
        self.import_mode_cb.addItems(["Whole Structure", "Only hit fragment"])

        self.middle_scroll_layout.addWidget(self.import_mode_label)
        self.middle_scroll_layout.addWidget(self.import_mode_cb)

    def hide_widgets(self):

        self.al_mode_label.setEnabled(False)
        self.al_mode_cb.setEnabled(False)
        self.import_mode_label.setEnabled(False)
        self.import_mode_cb.setEnabled(False)

    def show_widgets(self):

        self.al_mode_label.setEnabled(True)
        self.al_mode_cb.setEnabled(True)
        self.import_mode_label.setEnabled(True)
        self.import_mode_cb.setEnabled(True)
