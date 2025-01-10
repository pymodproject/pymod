# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

from pymol.Qt import QtWidgets, QtCore

"""
GUI for ESMFold in PyMod.
"""

class ESM_search_window_qt(QtWidgets.QMainWindow):
    """
    A class to represent the 'ESMFold Window' of PyMod.
    """

    is_pymod_window = True

    def __init__(self, parent, protocol):

        super(ESM_search_window_qt, self).__init__(parent)
        self.protocol = protocol
        self.initUI()

    def initUI(self):

        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle("ESMFold-DB search")

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
        self.main_button.clicked.connect(lambda a=None: self.protocol.launch_esm_search())
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
        self.sequence_rb = QtWidgets.QRadioButton("Sequence")

        self.middle_scroll_layout.addWidget(self.search_by_label)
        self.middle_scroll_layout.addWidget(self.sequence_id_rb)
        self.middle_scroll_layout.addWidget(self.sequence_rb)
