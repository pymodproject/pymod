# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.


from pymol.Qt import QtWidgets, QtCore
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtGui import QFont




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


    ##################MODIFICATION on 11/07/2025################
    #####Unique run prediction  button to start the prediction , in substtution to submit button##########
    ###NOT BOTTOM FRAME PRESENT #########
    def initUI(self):
        #########################
        # Configure the window. #
        #########################
        self.setWindowTitle("ESMFold")
        self.setGeometry(600, 300, 400, 200)  

        # Set central widget
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # Main vertical layout
        self.main_vbox = QtWidgets.QVBoxLayout()
        self.central_widget.setLayout(self.main_vbox)

        #################
        # Middle frame. #
        #################
        self.middle_widget = QtWidgets.QWidget()
        self.main_vbox.addWidget(self.middle_widget)

        self.middle_gridlayout = QtWidgets.QGridLayout()
        self.middle_widget.setLayout(self.middle_gridlayout)

        self.middle_scroll_area = QtWidgets.QScrollArea()
        self.middle_scroll_area.setWidgetResizable(True)
        self.middle_gridlayout.addWidget(self.middle_scroll_area)

        self.build_middle_widgets()

        # Optional: center alignment
        self.main_vbox.setAlignment(self.middle_widget, QtCore.Qt.AlignCenter)
    ###############END###############


    ##################MODIFICATION on 11/07/2025################
    #####Unique bar to run the prediction in the middle widgets##############
    def build_middle_widgets(self):
        # Create a large font object for the button
        big_font = QtGui.QFont()
        big_font.setPointSize(12)

        # Create a vertical layout for the scroll area and center-align its contents
        self.middle_scroll_layout = QtWidgets.QVBoxLayout()
        self.middle_scroll_layout.setAlignment(QtCore.Qt.AlignCenter)
        self.middle_scroll_area.setLayout(self.middle_scroll_layout)

        # ---------- "Run Prediction" Button ----------
        self.run_prediction_button = QtWidgets.QPushButton("Run Fold Sequence Prediction")
        self.run_prediction_button.setFont(big_font)
        self.run_prediction_button.setFixedHeight(50)
        self.run_prediction_button.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.run_prediction_button.setStyleSheet("font-weight: bold;")

        # Connect the button's click signal to the run_prediction method
        self.run_prediction_button.clicked.connect(self.run_prediction)

        # Add the button to the main vertical layout
        self.main_vbox.addWidget(self.run_prediction_button)

        # ---------- Progress Bar ----------
        self.prediction_progress_bar = QtWidgets.QProgressBar()
        self.prediction_progress_bar.setRange(0, 0)  # Set to indeterminate mode (infinite loop animation)
        self.prediction_progress_bar.setVisible(False)  # Initially hidden
        self.prediction_progress_bar.setFixedWidth(200)

        # ---------- Centered Layout Container ----------
        center_layout = QtWidgets.QVBoxLayout()
        center_layout.setAlignment(QtCore.Qt.AlignCenter)

        # Add the button and the progress bar to the center layout with spacing in between
        center_layout.addWidget(self.run_prediction_button)
        center_layout.addSpacing(15)
        center_layout.addWidget(self.prediction_progress_bar)

        # Create a container widget to hold the center layout
        center_widget = QtWidgets.QWidget()
        center_widget.setLayout(center_layout)

        # Set minimum widths to help layout adjust gracefully
        center_widget.setMinimumWidth(300)
        self.middle_scroll_area.setMinimumWidth(350)

        # Add vertical stretch before and after to vertically center the content
        self.middle_scroll_layout.addStretch()
        self.middle_scroll_layout.addWidget(center_widget)
        self.middle_scroll_layout.addStretch()


    def run_prediction(self):
        # Set progress bar to indeterminate mode and make it visible
        self.prediction_progress_bar.setRange(0, 0)
        self.prediction_progress_bar.setVisible(True)

        # Disable the button to prevent multiple clicks during processing
        self.run_prediction_button.setEnabled(False)

        # Start the actual prediction task after a short delay (1 second)
        # This allows the GUI to update and display the spinning bar first
        QtCore.QTimer.singleShot(1000, self._run_real_prediction)


    def _run_real_prediction(self):
        # Force the GUI to process any pending events (important for showing the progress bar)
        QtWidgets.QApplication.processEvents()

        # Run the real prediction function (this is blocking!)
        self.protocol.launch_esm_search()

        # After completion, stop the progress bar and re-enable the button
        self.prediction_complete()


    def prediction_complete(self):
        # Hide the progress bar
        self.prediction_progress_bar.setVisible(False)

        # Re-enable the prediction button so the user can run again
        self.run_prediction_button.setEnabled(True)
    ##############END################



    