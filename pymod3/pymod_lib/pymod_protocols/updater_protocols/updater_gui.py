# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import traceback
from socket import gaierror
import json

from pymol.Qt import QtWidgets, QtCore, QtGui


class Ui_UpdateDialog:

    n_cols = 5 # Number of columns in the updater window frame.
    component_col_idx = 0
    databases_col_idx = 1
    status_col_idx = 2
    source_col_idx = 3
    last_download_col_idx = 4

    def setupUi(self, UpdateDialog):
        UpdateDialog.setObjectName("UpdateDialog")
        UpdateDialog.resize(950, 470)
        self.verticalLayout = QtWidgets.QVBoxLayout(UpdateDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.select_comp_label = QtWidgets.QLabel(UpdateDialog)
        self.select_comp_label.setObjectName("select_comp_label")
        self.verticalLayout.addWidget(self.select_comp_label)
        self.components_tableWidget = QtWidgets.QTableWidget(UpdateDialog)
        default_font = QtGui.QFont()
        default_font.setPointSize(default_font.pointSize()-1)
        self.components_tableWidget.setFont(default_font)
        self.components_tableWidget.setProperty("showDropIndicator", False)
        self.components_tableWidget.setDragDropOverwriteMode(False)
        self.components_tableWidget.setAlternatingRowColors(True)
        self.components_tableWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.components_tableWidget.setGridStyle(QtCore.Qt.NoPen)
        self.components_tableWidget.setColumnCount(self.n_cols)
        self.components_tableWidget.setObjectName("components_tableWidget")
        self.components_tableWidget.setRowCount(1)
        vertical_header_item = QtWidgets.QTableWidgetItem()
        self.components_tableWidget.setVerticalHeaderItem(0, vertical_header_item)
        for i in range(self.n_cols):
            item = QtWidgets.QTableWidgetItem()
            self.components_tableWidget.setHorizontalHeaderItem(i, item)

        item = QtWidgets.QTableWidgetItem()
        self.components_tableWidget.setItem(0, self.component_col_idx, item)
        item = QtWidgets.QTableWidgetItem()
        self.components_tableWidget.setItem(0, self.databases_col_idx, item)
        self.components_tableWidget.horizontalHeader().setVisible(True)
        self.components_tableWidget.horizontalHeader().setCascadingSectionResizes(False)
        self.components_tableWidget.setColumnWidth(self.component_col_idx, 210)
        self.components_tableWidget.setColumnWidth(self.databases_col_idx, 180)
        self.components_tableWidget.setColumnWidth(self.status_col_idx, 190)
        # self.components_tableWidget.setColumnWidth(self.source_col_idx, 390)
        # self.components_tableWidget.setColumnWidth(self.last_download_col_idx, 250)

        self.components_tableWidget.horizontalHeader().setStretchLastSection(True)
        self.components_tableWidget.verticalHeader().setVisible(False)
        self.components_tableWidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.verticalLayout.addWidget(self.components_tableWidget)
        self.statusLabel = QtWidgets.QLabel(UpdateDialog)
        self.statusLabel.setObjectName("statusLabel")
        self.verticalLayout.addWidget(self.statusLabel)
        self.installation_progressBar = QtWidgets.QProgressBar(UpdateDialog)
        self.installation_progressBar.setEnabled(True)
        self.installation_progressBar.setProperty("value", 10)
        self.installation_progressBar.setObjectName("installation_progressBar")
        self.verticalLayout.addWidget(self.installation_progressBar)
        self.buttonsHorizontalLayout = QtWidgets.QHBoxLayout()
        self.buttonsHorizontalLayout.setObjectName("buttonsHorizontalLayout")
        self.installSel_button = QtWidgets.QPushButton(UpdateDialog)
        self.installSel_button.setObjectName("installSel_button")
        self.buttonsHorizontalLayout.addWidget(self.installSel_button)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.buttonsHorizontalLayout.addItem(spacerItem)
        self.cancel_button = QtWidgets.QPushButton(UpdateDialog)
        self.cancel_button.setObjectName("cancel_button")
        self.buttonsHorizontalLayout.addWidget(self.cancel_button)
        self.verticalLayout.addLayout(self.buttonsHorizontalLayout)

        UpdateDialog.setWindowTitle("Install and update databases")
        self.select_comp_label.setText("Select components to install or update")
        self.components_tableWidget.setSortingEnabled(True)
        self.components_tableWidget.horizontalHeaderItem(self.component_col_idx).setText("Component")
        self.components_tableWidget.horizontalHeaderItem(self.databases_col_idx).setText("Databases Names")
        self.components_tableWidget.horizontalHeaderItem(self.status_col_idx).setText("Status")
        self.components_tableWidget.horizontalHeaderItem(self.source_col_idx).setText("Source")
        self.components_tableWidget.horizontalHeaderItem(self.last_download_col_idx).setText("Last Downloaded")
        __sortingEnabled = self.components_tableWidget.isSortingEnabled()
        self.components_tableWidget.setSortingEnabled(False)
        self.components_tableWidget.setSortingEnabled(__sortingEnabled)
        self.statusLabel.setText("")
        self.cancel_button.setText("Cancel")
        self.installSel_button.setText("Install Selected")

        QtCore.QMetaObject.connectSlotsByName(UpdateDialog)


class InstallerUpdaterWindow(QtWidgets.QDialog, Ui_UpdateDialog):

    is_pymod_window = True

    def __init__(self, parent=None, installer_protocol=None):

        super(InstallerUpdaterWindow, self).__init__(parent)
        self.setupUi(self)
        self.installer_protocol = installer_protocol
        self._bind_thread_objects()                # sets the thread objects as children of this window

        self.view = self.components_tableWidget
        self.items_status_dict = {}
        self.items_last_download_dict = {}

        #puts the 'Install selected' button as default
        # self.installSel_button.setFocus()

        # shows the progress bar
        self.installation_progressBar.setValue(0)
        self.installation_progressBar.setVisible(True)

        # connections
        self._connections()

        # shows in the ListView the list of PyMod components
        self._fill_list()

        # shows the window itself
        self.show()

        # ping to server
        self.start_internet_check()


    def _connections(self):
        """Connecting signals with slots."""
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        self.installSel_button.clicked.connect(self.install_selected_slot)

        for comp in self.installer_protocol.components_list:
            # Used to show the size of the download before starting the donwload.
            comp.installer.retrieved_size.connect(self.set_starting_status_slot)
            # Updates the text in the "Status" column.
            comp.installer.set_update_status.connect(self.set_update_status_slot)
            # Terminates the installation.
            comp.installer.terminated_installation.connect(self.terminated_installation_slot)
            # Called when the installation of a component fails.
            comp.installer.critical_error.connect(self.critical_error_slot)
            # Shows a information messagebox.
            comp.installer.info_message.connect(self.info_slot)

    def on_cancel_button_click(self):
        self._terminate_threads()
        self.close()

    def _terminate_threads(self):
        for comp in self.installer_protocol.selected_comps:
            if comp.installer.isRunning():
                comp.installer.installation_status = "stopped"

    def closeEvent(self, evnt):
        self._terminate_threads()

    def show_message_in_label(self, message):
        self.statusLabel.setText(message)

    def activate_progressbar(self):
        self.installation_progressBar.setMinimum(0)
        self.installation_progressBar.setMaximum(0)
        # self.show_message_in_label("Updating %s components" % self.installer_protocol.pending_downloads)

    def terminated_installation_slot(self, component):
        # Updates the "Last Downloaded" column.
        self.items_last_download_dict[component.name].setText(component.last_downloaded)
        self.updated_gui_on_thread_finish("installation_terminated")

    def critical_error_slot(self, error_message, comp):
        """Aborts installation if an error occurred"""
        comp.installer.installation_status = "failed"
        comp.installer.exit(1)
        self.set_update_status_slot(comp, "Updated Failed: %s" % error_message, "red")
        self.updated_gui_on_thread_finish("failure")


    def updated_gui_on_thread_finish(self, status):

        if not status in ("download_terminated", "installation_terminated", "failure"):
            raise KeyError(status)

        # downloaded_comps = [comp for comp in self.installer_protocol.selected_comps if comp.installer.installation_status == "downloaded"]
        failed_comps = [comp for comp in self.installer_protocol.selected_comps if comp.installer.installation_status == "failed"]
        installed_comps = [comp for comp in self.installer_protocol.selected_comps if comp.installer.installation_status == "success"]
        running_comps = [comp for comp in self.installer_protocol.selected_comps if comp.installer.installation_status != "failed"]

        if status == "failure":
            if len(failed_comps) == len(self.installer_protocol.selected_comps):
                # self.show_message_in_label("Update failed.")
                self.installation_progressBar.setMaximum(100)
                self.installation_progressBar.setValue(0)

        elif status == "installation_terminated":
            # All the components have been installed.
            if len(installed_comps) == len(running_comps):
                # self.show_message_in_label("Update completed.")
                self.installation_progressBar.setMaximum(100)
                self.installation_progressBar.setValue(100)
                self.installer_protocol.finish_installation()


    def info_slot(self, message, comp):
        """shows an info message, does not abort the installation"""
        self.installer_protocol.pymod.main_window.show_info_message("Warning", message)


    def set_starting_status_slot(self, component, size_download):
        """It responds to the established_connection signal emitted by the thread.
        The signal carries the file size, that will be displayed in the window."""
        # gets the graphic object that displays the status, from a dictionary previously created
        status_item_object = self.items_status_dict[component.name]# [self.status_col_idx]
        status_item_object.setText('Available: '+str(int(size_download/1048576))+' MB')
        status_item_object.setForeground(QtGui.QBrush(QtGui.QColor(204, 255, 204)))


    def set_update_status_slot(self, component, text="new status", color="light_green"):
        status_item_object = self.items_status_dict[component.name]
        status_item_object.setText(text)
        if color == "light_green":
            status_item_object.setForeground(QtGui.QBrush(QtGui.QColor(204, 255, 204)))
        elif color == "green":
            status_item_object.setForeground(QtGui.QBrush(QtGui.QColor(10, 255, 10)))
        elif color == "gray":
            status_item_object.setForeground(QtGui.QBrush(QtGui.QColor(200, 200, 200)))
        elif color == "red":
            status_item_object.setForeground(QtGui.QBrush(QtGui.QColor(255, 0, 0)))
        else:
            raise KeyError(color)


    def install_selected_slot(self):
        """This is the action executed after the 'Install Selected' button is pressed. It calls
        the 'install_selected' method of the Updater Protocol class."""

        # Check if there are any selected databases.
        self.installer_protocol.selected_comps = []
        for comp in self.items_dict.keys():
            # Check the status of each checkbutton and get the ones having been selected by the user.
            if self.items_dict[comp].checkState() == QtCore.Qt.Checked:
                self.installer_protocol.selected_comps.append(comp)

        if not self.installer_protocol.selected_comps:
            self.info_slot("Please select at least one database to download.", None)
            return None

        # Check if the selected databases can be downloaded.
        for cmp in self.installer_protocol.selected_comps:

            # Check if the server has been reached while pinging before.
            if not cmp.can_be_downloaded:
                error_mess = "The %s database can not be currently downloaded. Please uncheck it." % (cmp.full_name)
                self.info_slot(error_mess, None)
                return None

            # Trying to figure out if there are hmmer executables. If not, blocks the installation,
            # because it can't index the databases without hmmpress.
            if cmp.name == 'hmmscan_databases':
                if not os.path.isfile(cmp.installer.hmmpress_exe_filepath):
                    error_mess = ("Cannot install PFAM databases since a 'hmmscan' executable was"
                                  " not found in the HMMER executables directory specified in the"
                                  " PyMod options window ('%s'). In order to update the PFAM database,"
                                  " please provide an 'hmmscan' executable, it is necessary for the"
                                  " indicization of the compressed database" % os.path.dirname(cmp.installer.hmmpress_exe_filepath))
                    self.info_slot(error_mess, None)
                    return None

            # Sets the target path for the components.
            self.installer_protocol.set_component_target_path(cmp)

            # Checks the target path.
            if cmp.target_installation_path == None:
                error_mess = ("Cannot install %s databases, because its databases directory is not"
                              " defined in the PyMod options. Please uncheck it or define it in the"
                              " PyMod options window." % (cmp.full_name))
                self.info_slot(error_mess, None)
                return None
            if not os.path.isdir(cmp.target_installation_path):
                error_mess = ("Cannot install %s databases, because its databases directory defined in"
                              " the PyMod options (%s) does not exists. Please uncheck it or define an"
                              " existing directory in the PyMod options window." % (cmp.full_name, cmp.target_installation_path))
                self.info_slot(error_mess, None)
                return None

        self.installSel_button.setEnabled(False) # Freezing the button.
        self.installer_protocol.install_selected()


    def _bind_thread_objects(self):
        """this method sets the window as the Thread parent, in order to work on the threads as GUI children """
        for c in self.installer_protocol.components_list:
            c.installer.setParent(self)


    def start_internet_check(self):
        """Ping method that also retrieves the file size"""
        for comp in self.installer_protocol.components_list:
            # this flag tells the installer not to download everything but only to ping the server
            comp.installer.install_mode = "ping"
            comp.installer.start()


    def _fill_list(self):
        """
        Create the list's data.
        """

        self.items_dict = {}

        # Checks the database log in order to obtain the date when each database was last downloaded.
        download_log_dict = {}
        if os.path.isfile(self.installer_protocol.download_log_filepath):
            with open(self.installer_protocol.download_log_filepath, "r") as l_fh:
                download_log_dict = json.loads(l_fh.read())

        # Configure the list of items is the 'all_components_list' from the PyMod Installer class.
        self.view.setRowCount(len(self.installer_protocol.components_list))
        for row_counter, component in enumerate(self.installer_protocol.components_list):
            # Create an item and set the component name as text.
            item = QtWidgets.QTableWidgetItem(component.full_name)
            # add a checkbox to the name
            item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            item.setCheckState(QtCore.Qt.Unchecked)
            # place the items in columns
            self.view.setItem(row_counter, self.component_col_idx, item)
            self.items_dict.update({component: item})


            # Set the names of the databases.
            self.view.setItem(row_counter, self.databases_col_idx, QtWidgets.QTableWidgetItem(component.databases_string))

            # Create another item displaying the status.
            graphic_status = 'Wait...'
            status = QtWidgets.QTableWidgetItem(graphic_status)
            status.setForeground(QtGui.QBrush(QtGui.QColor(191, 191, 191)))
            self.view.setItem(row_counter, self.status_col_idx, status)

            # Set the source URL.
            self.view.setItem(row_counter, self.source_col_idx, QtWidgets.QTableWidgetItem(component.remote_source))

            # Fill in the last downloaded column.
            if component.name in download_log_dict:
                last_downloaded_str = download_log_dict[component.name]
            else:
                last_downloaded_str = "Never"
            last_downloaded_item = QtWidgets.QTableWidgetItem(last_downloaded_str)
            self.view.setItem(row_counter, self.last_download_col_idx, last_downloaded_item)

            self.items_status_dict[component.name] = status
            self.items_last_download_dict[component.name] = last_downloaded_item


    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            pass
        else:
            QtWidgets.QDialog.keyPressEvent(self, event)
