# Copyright 2025 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

from pymol.Qt import QtWidgets, QtCore

"""
GUI for Boltz-1 in PyMod.
"""

class Boltz_search_window_qt(QtWidgets.QMainWindow):

    """
    A class to represent the Boltz-1 window in PyMod
    """

    is_pymod_window = True

    def __init__(self, parent, protocol, selected_sequences):

        super(Boltz_search_window_qt, self).__init__(parent)
        self.protocol = protocol
        self.selected_sequences = selected_sequences
        self.initUI()

    def initUI(self):

        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle("Boltz-1 Modeling")
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
        self.main_button.clicked.connect(lambda a=None: self.protocol.run_boltz_modeling())
        self.main_vbox.addWidget(self.main_button)
        self.main_button.setFixedWidth(self.main_button.sizeHint().width())


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


    def build_middle_widgets(self):

        self.middle_scroll_layout = QtWidgets.QVBoxLayout()
        self.middle_scroll_area.setLayout(self.middle_scroll_layout)

        # Ligand section 
        self.ligand_label = QtWidgets.QLabel("Ligand")

        # Get hetresidues if present 
        list_of_hetres = self.selected_sequences[0].get_heteroresidues()
        self.list_of_cb = []

        if list_of_hetres:
            for hetres in list_of_hetres:
                cb = QtWidgets.QCheckBox(hetres.three_letter_code)
                self.list_of_cb.append(cb)

        self.ligand_line_edit = QtWidgets.QLineEdit()
        # self.ligand_from_pymod.clicked.connect(self.show_widgets)
        # self.ligand_manually.clicked.connect(self.hide_widgets)

        self.middle_scroll_layout.addWidget(self.ligand_label)
        for cb in self.list_of_cb:
            self.middle_scroll_layout.addWidget(cb)
        self.middle_scroll_layout.addWidget(self.ligand_line_edit)

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
