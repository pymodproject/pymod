# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for a PyQt window with plots built with the pyqtgraph library (http://www.pyqtgraph.org/).
NumPy and pyqtgraph are required.
The PyMod plugin package integrates a modified version of pyqtgraph (which can be
imported through 'import pymod_lib.pymod_plot.pyqtgraph'). This version includes
a series of modifications to make pyqtgraph compatible with the current PyMOL
versions (see the comment lines starting with "PYMOD FIX" in the package).
"""

import io
import math

import numpy as np

from pymol.Qt import QtWidgets, QtCore, QtGui

from pymod_lib.pymod_plot import pyqtgraph
from pymod_lib.pymod_gui.shared_gui_components_qt import asksaveasfile_qt, small_font_style
from pymod_lib.pymod_os_specific import check_valid_pyqt


###############################################################################
# Plotting window in Pyqt for PyMod, using pyqtgraph.                         #
###############################################################################

class PyMod_plot_window_qt(QtWidgets.QMainWindow):
    """
    Class for a PyQt window with a plot built using the pyqtgraph library and
    a series of widgets to control the plot and interact with it.
    """

    is_pymod_window = True

    #---------------------------------
    # Configure the plotting window. -
    #---------------------------------

    control_colors = "black"
    plot_colors = ["#0000ff", "#00ff00", "#ff0000", "#00ffff",
                   "#ff00ff", "#ffff00", "#000000", "#cccccc"]

    def __init__(self, *args, **kwargs):

        super(PyMod_plot_window_qt, self).__init__(*args, **kwargs)

        # Central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        self.central_widget_layout = QtWidgets.QGridLayout()
        self.central_widget.setLayout(self.central_widget_layout)


        #------------------------------------------------
        # Upper frame (contains the plot and controls). -
        #------------------------------------------------

        expanding_size_policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                                      QtWidgets.QSizePolicy.Expanding)

        preferred_size_policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred,
                                                      QtWidgets.QSizePolicy.Preferred)

        # The upper frame contains three frames: info, plot and controls frames.
        # The infor and controls frames will be displayed only if the 'use_controls'
        # argument is set to 'True' when calling the 'build_plotting_area' method.
        self.upper_frame = QtWidgets.QFrame()
        self.upper_frame.setStyleSheet("background-color: #646464")
        self.upper_frame_layout = QtWidgets.QGridLayout()
        self.upper_frame.setLayout(self.upper_frame_layout)
        self.upper_frame.setSizePolicy(expanding_size_policy)
        self.central_widget_layout.addWidget(self.upper_frame, 0, 0)

        # Info frame, it contains the messagebar of the plot.
        self.info_frame = QtWidgets.QFrame()
        self.info_frame_layout = QtWidgets.QHBoxLayout()
        self.info_frame.setLayout(self.info_frame_layout)
        self.info_frame.setSizePolicy(preferred_size_policy)

        self.info_label = QtWidgets.QLabel("")
        self.info_frame_layout.addWidget(self.info_label)

        # Plot frame.
        self.plot_frame = QtWidgets.QFrame()
        # self.plot_frame.setStyleSheet("background-color: red")
        self.plot_frame_layout = QtWidgets.QGridLayout()
        self.plot_frame.setLayout(self.plot_frame_layout)
        self.plot_frame.setSizePolicy(expanding_size_policy)
        self.build_plot_widget()


        # Controls frame.
        self.controls_frame = QtWidgets.QWidget()
        self.controls_frame.setStyleSheet("background-color: #747474")
        self.controls_frame_layout = QtWidgets.QGridLayout()
        self.controls_frame.setLayout(self.controls_frame_layout)
        self.controls_frame_layout.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)

        self.controls_scrollarea = QtWidgets.QScrollArea()
        self.controls_scrollarea.setWidgetResizable(True)
        self.controls_scrollarea.setWidget(self.controls_frame)

        self.labels_title = QtWidgets.QLabel("Plots list")


        # Middle splitter.
        self.middle_splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.middle_splitter.setSizePolicy(expanding_size_policy)


        #---------------------------------------
        # Lower frame (contains some options). -
        #---------------------------------------

        self.lower_frame = QtWidgets.QFrame()
        self.lower_frame_layout = QtWidgets.QGridLayout()
        self.lower_frame.setLayout(self.lower_frame_layout)
        self.central_widget_layout.addWidget(self.lower_frame, 1, 0)

        # View buttons.
        self.view_label = QtWidgets.QLabel("View:")
        self.lower_frame_layout.addWidget(self.view_label, 0, 0)

        self.home_view_button = QtWidgets.QPushButton("Fit to data")
        self.home_view_button.clicked.connect(self.on_home_button_click)
        self.lower_frame_layout.addWidget(self.home_view_button, 0, 1)

        # On click behaviour. The buttons will be shown later, in the
        # 'build_plotting_area' metohd.
        self.interact_buttons_group = QtWidgets.QButtonGroup()

        self.on_click_label = QtWidgets.QLabel("Interact on click:")

        self.interact_button = QtWidgets.QPushButton("Yes")
        self.interact_button.setCheckable(True)
        self.interact_buttons_group.addButton(self.interact_button)
        self.interact_button.clicked.connect(self.on_interact_button_click)

        self.no_interaction_button = QtWidgets.QPushButton("No")
        self.no_interaction_button.setCheckable(True)
        self.no_interaction_button.setChecked(True)
        self.interact_buttons_group.addButton(self.no_interaction_button)
        self.no_interaction_button.clicked.connect(self.on_no_interaction_button_click)


        # Show/hide all buttons. They will be shown later, in the 'build_plotting_area'
        # method.
        self.show_label = QtWidgets.QLabel("Show:")

        self.show_all_button = QtWidgets.QPushButton("All")
        self.show_all_button.clicked.connect(self.show_all_command)

        self.hide_all_button = QtWidgets.QPushButton("None")
        self.hide_all_button.clicked.connect(self.hide_all_command)

        self.lower_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


        #---------------------
        # Build a main menu. -
        #---------------------

        self.save_to_csv_action = QtWidgets.QAction('Save to CSV', self)
        self.save_to_csv_action.triggered.connect(lambda a=None: self.save_to_csv_event())
        self.save_to_png_action = QtWidgets.QAction('Save to PNG', self)
        self.save_to_png_action.triggered.connect(lambda a=None: self.save_to_png_event())

        self.main_menubar = self.menuBar()
        self.file_menu = self.main_menubar.addMenu('File')


    def build_plot_widget(self):
        """
        Build and configure the PlotWidget of the plotting window.
        """

        self.plot_count = 0
        self.plot_items = {}
        self.plot_color_index = 0

        self.interact_mode = False
        self.selected_plot_item = None
        self.previous_highlight_dot = None
        self.previous_highlight_plot_id = None

        self.all_points = []
        self.all_points_info = []


        # PlotWidget from pyqtgraph.
        self.graphWidget = pyqtgraph.PlotWidget()
        # Configure the PlotWidget.
        self.graphWidget.setBackground('w')
        # self.graphWidget.setForeground(self.control_colors)
        self.graphWidget.hideButtons()

        self.plot_frame_layout.addWidget(self.graphWidget, 0, 0)


    def initialize(self, pymod, title):
        self.pymod = pymod
        self.setWindowTitle(title)


    #----------------------------------------------------------------
    # Configure the plotting area, add plots and iteract with them. -
    #----------------------------------------------------------------

    def build_plotting_area(self,
                            use_controls=True,
                            use_all_controls_buttons=True,
                            messagebar_initial_text=None,
                            update_messagebar=False,
                            messagebar_text_on_update="",
                            messagebar_vars_on_update=(),
                            on_click_action=None,
                            x_label_text="x",
                            y_label_text="y",
                            label_size=None,
                            use_save_to_csv=True,
                            hide_x_ticks=False,
                            hide_y_ticks=False,
                            # hide_x_label=False,
                            # hide_y_label=False,
                            highlight_points=True,
                            ):
        """
        Configures the plotting area of the window.

        # Arguments
            use_controls: adds to the plotting window a right column with checkbuttons
                to show/hide the lines plots in the drawing area.
            use_all_controls_buttons: adds to the plotting control column 'Show All'
                and 'Hide All' buttons.
            messagebar_initial_text: initial text to be displayed in the messagebar
                of the plotting window.
            update_messagebar: if 'True', the messagebar will be updated when clicking
                some point of the scatterplots.
            messagebar_text_on_update: text to be shown when the messagebar is update.
                If it contains the "__plot_name__", "__x__", "__y__" string, they
                will be substituted with the name of the plot, the x value and y value
                of the point respectively. If some 'additional_data' is provided
                for a plot, its data will also be used to update the messagebar.
                Each element of an 'additional_data' list has to correspond to a
                data point of a plot and must be a dictionary in which the keys are
                strings representing the names of the additional data series. If,
                for example, these dictionaries have a key named 'info' and if the
                'messagebar_text_on_update' contains the string "__info__", it will
                be substituted with the value associated to the 'info' key of that
                'additional_data' point.
            on_click_action: function to be called when a point of a scatterplot
                is clicked. This function must have the following argument:
                    point_data: a dictionary with additional data for the point
                        or a 'None' value.
            x_label_text: text for the x-axis label.
            y_label_text: text for the y-axis label.
            label_size: font size for the axis labels.
            use_save_to_csv: if 'True' add to the plotting window menu a command
                to save data in the csv format.
            hide_x_ticks: if 'True', hide the ticks and numbers of the x-axis.
            hide_y_ticks: if 'True', hide the ticks and numbers of the y-axis.
            highlight_points: if 'True', the scatterplot points will be highlighted
                when they are clicked with the mouse.

        """

        # Stores the arguments.
        self.use_controls = use_controls
        if self.use_controls:
            self.upper_frame_layout.addWidget(self.info_frame, 0, 0)

            self.middle_splitter.addWidget(self.plot_frame)
            self.middle_splitter.addWidget(self.controls_scrollarea)
            self.upper_frame_layout.addWidget(self.middle_splitter, 1, 0)

            self.controls_scrollarea.resize(230, self.controls_scrollarea.sizeHint().height())
            self.controls_frame_layout.addWidget(self.labels_title, 0, 0, 1, 2)

            self.lower_frame_layout.addWidget(self.on_click_label, 0, 2)
            self.lower_frame_layout.addWidget(self.interact_button, 0, 3)
            self.lower_frame_layout.addWidget(self.no_interaction_button, 0, 4)
        else:
            self.upper_frame_layout.addWidget(self.plot_frame, 0, 0)


        self.use_all_controls_buttons = use_all_controls_buttons
        if self.use_controls and self.use_all_controls_buttons:
            self.lower_frame_layout.addWidget(self.show_label, 0, 5)
            self.lower_frame_layout.addWidget(self.show_all_button, 0, 6)
            self.lower_frame_layout.addWidget(self.hide_all_button, 0, 7)

        self.messagebar_initial_text = messagebar_initial_text
        if self.messagebar_initial_text is not None:
            self.info_label.setText(self.messagebar_initial_text)

        self.update_messagebar = update_messagebar
        self.on_click_action = on_click_action
        if self.on_click_action is not None:
            if not hasattr(self.on_click_action, "__call__"):
                raise TypeError("'on_click_action' must be a function.")
        self.messagebar_text_on_update = messagebar_text_on_update
        self.messagebar_vars_on_update = messagebar_vars_on_update

        self.use_save_to_csv = use_save_to_csv
        if self.use_save_to_csv:
            self.file_menu.addAction(self.save_to_csv_action)
        self.file_menu.addAction(self.save_to_png_action)

        self.highlight_points = highlight_points


        # Configure the PlotWidget.
        self.graphWidget.setMenuEnabled(False)

        if label_size is not None:
            kwargs = {"font-size": "%spx" % label_size}
            font = QtGui.QFont()
            font.setPixelSize(label_size)
            # curve_pen = pyqtgraph.mkPen(width=2, color="r") # Color and width of the curve.
        else:
            kwargs = {}

        self.graphWidget.setLabel("left", y_label_text, **kwargs)
        self.graphWidget.setLabel("bottom", x_label_text, **kwargs)

        if label_size is not None:
            self.graphWidget.getAxis("left").tickFont = font
            self.graphWidget.getAxis("bottom").tickFont = font
            # self.graphWidget.getAxis("left").setPen(curve_pen)

        if hide_x_ticks:
            self.graphWidget.getAxis("bottom").setStyle(showValues=False, tickLength=0)
        if hide_y_ticks:
            self.graphWidget.getAxis("left").setStyle(showValues=False, tickLength=0)

        if self.on_click_action is not None:
            self.graphWidget.getViewBox().scene().sigMouseClicked.connect(self.on_scene_click)


    def on_scene_click(self, event):
        """
        Called whenever clicking on some region of the plot.
        """

        if not self.interact_mode:
            return None

        plot_point = self.graphWidget.getViewBox().mapSceneToView(event.scenePos())

        if hasattr(self.selected_plot_item, "_pymod_id"):
            self.on_point_click(plot_point)

        self.selected_plot_item = None


    def on_curve_click(self, item):
        """
        Sets the currently highlighted curve.
        """
        self.selected_plot_item = item


    def on_point_click(self, plot_point):
        """
        Called when a scatterplot point is clicked on the graph. If 'update_messagebar'
        if set to 'True', the messagebar will be updated. If an 'on_click_action'
        was provided, this function will also be executed.
        """

        if not self.interact_mode:
            return None

        plot_point_xy = (plot_point.x(), plot_point.y())

        min_dist = get_point_dist(self.all_points[0], plot_point_xy)
        min_id = 0

        for i, pi in enumerate(self.all_points):
            di = get_point_dist(pi, plot_point_xy)
            if di < min_dist:
                min_dist = di
                min_id = i
        plot_id, point_id = self.all_points_info[min_id]


        # Gets the data coordinates of the point.
        x_data = self.plot_items[plot_id]["x_data"][point_id]
        y_data = self.plot_items[plot_id]["y_data"][point_id]


        # Shows a circle around the selected point on the curve.
        if self.previous_highlight_dot is not None:
            self.graphWidget.removeItem(self.previous_highlight_dot)

        self.previous_highlight_dot = self.graphWidget.plot([x_data], [y_data],
                                                             symbol="o",
                                                             symbolSize=12,
                                                             symbolBrush="y")
        self.previous_highlight_plot_id = plot_id

        # Updates the message bar in the upper part of the plotting window.
        if self.update_messagebar:

            # Label of the plot.
            plot_label = self.plot_items[plot_id]["label"]

            # Build the new message.
            updated_text = self.messagebar_text_on_update.replace("__plot_name__", plot_label)
            updated_text = updated_text.replace("__x__", str(x_data))
            updated_text = updated_text.replace("__y__", str(y_data))
            # Additional data.
            plot_additional_data = self.plot_items[plot_id]["additional_data"]
            if plot_additional_data is not None:
                point_additional_data = plot_additional_data[point_id]
                for data_key in point_additional_data:
                    updated_text = updated_text.replace("__%s__" % data_key, str(point_additional_data[data_key]))
            else:
                point_additional_data = None

            # Set the new message.
            self.info_label.setText(updated_text)


        # Calls the custom function.
        if self.on_click_action is not None:
            self.on_click_action(point_additional_data)


    def add_plot(self, x_data, y_data,
                 label=None,
                 additional_data=None):
        """
        Adds a plot to the pyqtgraph PlotWidget.
        """

        # Check and prepare the data.
        if len(x_data) != len(y_data):
            raise ValueError(("The x series and the y series do not have the same"
                              " number of elements (%s and %s respectively)" % (len(x_data), len(y_data))))
        if additional_data:
            if not len(x_data) == len(additional_data):
                raise ValueError(("The 'additional_data' series does not have the"
                                  " same number of elements of the data to plot"
                                  " (%s and %s respectively)" % (len(x_data), len(additional_data))))

        _x_data, _y_data = self.remove_none(x_data, y_data)


        # Add the plot to the PlotWidget.
        plot_color = self.plot_colors[self.plot_color_index]
        self.change_plot_color_index()

        curve_pen = pyqtgraph.mkPen(width=2, color=plot_color) # Color and width of the curve.

        plot_item = self.graphWidget.plot(_x_data, _y_data,
                                          name=label,
                                          connect="finite",
                                          pen=curve_pen,
                                          clickable=True)
        plot_item._pymod_id = self.plot_count
        plot_item.curve._pymod_id = self.plot_count

        plot_item.curve.setClickable(True)
        # plot_item.curve.sigClicked.connect(self.on_curve_click)
        plot_item.sigClicked.connect(self.on_curve_click)


        # Store information about the plot.
        self.plot_items[self.plot_count] = {}
        # The 'plot' key will store the pyqtgraph object for the plot.
        self.plot_items[self.plot_count]["plot"] = plot_item
        # The 'state' will be 1 if the plot is shown, and 0 it is hidden.
        self.plot_items[self.plot_count]["state"] = 1
        # Add a label.
        if label is None:
            label = "Plot %s" % self.plot_count
        self.plot_items[self.plot_count]["label"] = label
        # Data.
        self.plot_items[self.plot_count]["x_data"] = x_data
        self.plot_items[self.plot_count]["y_data"] = y_data
        # Additional data.
        self.plot_items[self.plot_count]["additional_data"] = additional_data


        # Stores all the data in a single list.
        for idx, (xi, yi) in enumerate(zip(x_data, y_data)):
            if yi is not None:
                self.all_points.append((xi, yi))
                self.all_points_info.append((self.plot_count, idx))


        # Add a checkbox in the controls frame. Used for showing/hiding the plot.
        if self.use_controls:
            plot_checkbox = QtWidgets.QCheckBox(label)
            plot_checkbox.setStyleSheet(small_font_style)
            plot_checkbox.setChecked(True)
            plot_checkbox.clicked.connect(lambda a=None, i=self.plot_count: self.toggle_plot(i))
            plot_color_label = QtWidgets.QLabel("---") # \u2796") # "\u25CF"
            plot_color_label.setStyleSheet("color: %s; font-weight: bold" % plot_color)
            self.controls_frame_layout.addWidget(plot_color_label, self.plot_count+1, 0, 1, 1)
            self.controls_frame_layout.addWidget(plot_checkbox, self.plot_count+1, 1, 1, 1)

            self.plot_items[self.plot_count]["checkbox"] = plot_checkbox

        # Increase the plot counter.
        self.plot_count += 1


    def change_plot_color_index(self):
        if self.plot_color_index == len(self.plot_colors) - 1:
            self.plot_color_index = 0
        else:
            self.plot_color_index += 1


    def remove_none(self, x_values, y_values):
        """
        Used to remove 'None' values for a data series ('None' values are not
        supported by pyqtgraph, which supports numpy nan values instead).
        """

        _x_values = []
        _y_values = []
        has_valid_pyqt = check_valid_pyqt()

        # nan values are supported by PyQt.
        if has_valid_pyqt:

            for xi, yi in zip(x_values, y_values):
                if yi is None:
                    _y_values.append(np.nan)
                else:
                    _y_values.append(yi)
                _x_values.append(xi)

        # nan values are not supported by PyQt.
        else:

            for xi, yi in zip(x_values, y_values):
                if yi is None:
                    pass
                else:
                    _x_values.append(xi)
                    _y_values.append(yi)

        return _x_values, _y_values


    def toggle_plot(self, plot_idx):
        """
        Called when clicking on some checkbox on the control frame on the right
        of the plotting window.
        """
        plot_info = self.plot_items[plot_idx]
        if plot_info["state"] == 1:
            self.graphWidget.removeItem(plot_info["plot"])
            plot_info["state"] = 0
            if self.previous_highlight_plot_id == plot_idx:
                self.remove_highlight_dot()
        else:
            self.graphWidget.addItem(plot_info["plot"])
            plot_info["state"] = 1


    def remove_highlight_dot(self):
        self.graphWidget.removeItem(self.previous_highlight_dot)
        self.previous_highlight_dot = None
        self.previous_highlight_plot_id = None


    #----------------------------------------------------------------
    # Events called when pressing the control buttons of the window.-
    #----------------------------------------------------------------

    def on_home_button_click(self):
        self.graphWidget.autoBtnClicked()


    def on_interact_button_click(self):
        self.interact_mode = True

    def on_no_interaction_button_click(self):
        if self.previous_highlight_plot_id is not None:
            self.remove_highlight_dot()
        self.interact_mode = False


    def show_all_command(self):
        for plot_idx in self.plot_items:
            plot_info = self.plot_items[plot_idx]
            if plot_info["state"] == 0:
                self.graphWidget.addItem(plot_info["plot"])
                plot_info["state"] = 1
                plot_info["checkbox"].setChecked(True)

    def hide_all_command(self):
        for plot_idx in self.plot_items:
            plot_info = self.plot_items[plot_idx]
            if plot_info["state"] == 1:
                self.graphWidget.removeItem(plot_info["plot"])
                plot_info["state"] = 0
                plot_info["checkbox"].setChecked(False)
                if self.previous_highlight_plot_id == plot_idx:
                    self.remove_highlight_dot()


    #-----------------------
    # Save to file events. -
    #-----------------------

    def save_to_csv_event(self):
        output_filepath = asksaveasfile_qt("Save to CSV file",
                                           parent=self.pymod.get_qt_parent(),
                                           name_filter="*.csv")
        if not output_filepath:
            return None
        try:
            self._save_to_csv_file(output_filepath)
        except Exception as e:
            print("- WARNING: could not write a csv file: %s" % str(e))

    def _save_to_csv_file(self, output_filepath):

        output = io.StringIO()

        # Writes the header string.
        header_string = []
        for plot_idx in range(self.plot_count):
            header_string.append(self.plot_items[plot_idx]["label"])
            header_string.append(self.plot_items[plot_idx]["label"] + " info")
        header_string = ",".join(header_string)
        print(header_string, file=output)

        # Write each data point information.
        max_points_count = max([len(self.plot_items[idx]["x_data"]) for idx in range(self.plot_count)])
        for point_i in range(0, max_points_count):
            line_string = []
            for plot_idx in range(self.plot_count):
                try:
                    # Get the y value.
                    point_y = self.plot_items[plot_idx]["y_data"][point_i]
                    point_val = str(point_y) if point_y is not None else ""
                    # Get the additional data for the point.
                    adata = self.plot_items[plot_idx]["additional_data"]
                    point_additional_data = ""
                    if adata is not None:
                        point_adata = self.plot_items[plot_idx]["additional_data"][point_i]
                        if "export_label" in point_adata:
                            if point_adata["export_label"] != None:
                                point_additional_data = str(point_adata["export_label"])
                    line_string.extend([point_val, point_additional_data])

                except IndexError:
                    line_string.extend(["", ""])

            line_string = ",".join(line_string)
            print(line_string, file=output)

        contents = output.getvalue()
        output.close()

        output_file_handler = open(output_filepath, "w")
        print(contents, file=output_file_handler)
        output_file_handler.close()


    def save_to_png_event(self):
        output_filepath = asksaveasfile_qt("Save to PNG file",
                                           parent=self.pymod.get_qt_parent(),
                                           name_filter="*.png")
        if not output_filepath:
            return None

        try:
            from pymod_lib.pymod_plot.pyqtgraph.exporters import ImageExporter
            # Create an exporter instance, as an argument give it the whole plot.
            exporter = ImageExporter(self.graphWidget.plotItem)
            # Save to file.
            exporter.export(output_filepath)

        except Exception as e:
            print("- WARNING: could not write a pgn file: %s" % str(e))


def get_point_dist(point_i, point_j):
    return math.sqrt((point_i[0]-point_j[0])**2 + (point_i[1]-point_j[1])**2)


###############################################################################
# Builds a plot showing a distance tree.                                      #
###############################################################################

# The following code has been adapted from the 'draw' method of the 'Phylo' module
# of Biopython.
# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
'''
                 Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
'''
def draw_tree(tree, pymod, label_func=str, do_show=True, show_confidence=True,
              tree_lines_width=1, tree_color="#000000",
              labels_font_size="12px", labels_color="#404040"):
    """
    Draws a distance tree contained in a 'Tree' class object of Bio.Phylo using PyMod plotting
    engine. It mainly uses the algorithm implemented in the 'draw' method of 'Phylo'.
    """

    #############################################################################
    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Layout
    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = dict((tip, maxheight - i)
                       for i, tip in enumerate(reversed(tree.get_terminals())))

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (heights[clade.clades[0]] +
                              heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights
    #############################################################################

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)

    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    xlim = (-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    ylim = (0.2, max(y_posns.values()) + 0.8) # (max(y_posns.values()) + 0.8, 0.2)

    # Aesthetics
    plot_title = "Distance tree"
    if hasattr(tree, 'name') and tree.name:
        plot_title = tree.name

    # Builds the plot window and initializes the plotting area.
    cp = PyMod_plot_window_qt(pymod.main_window)
    cp.initialize(pymod=pymod, title=plot_title)
    cp.build_plotting_area(use_controls=False,
                           x_label_text="Branch length", y_label_text="Taxum",
                           hide_y_ticks=True,
                           use_all_controls_buttons=False,
                           use_save_to_csv=False)
    cp.show()

    #############################################################################
    def draw_clade_lines(use_linecollection=True, orientation='horizontal',
                         y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0,
                         color='black', lw='.1'):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == 'horizontal':
            horizontal_linecollections.append([x_start, y_here, x_here, y_here])

        elif orientation == 'vertical':
            vertical_linecollections.append([x_here, y_bot, x_here, y_top])


    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, 'color') and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, 'width') and clade.width is not None:
            lw = 1 # GX: clade.width * plt.rcParams['lines.linewidth']
        # Draw a horizontal line from start to here
        draw_clade_lines(use_linecollection=True, orientation='horizontal',
                         y_here=y_here, x_start=x_start, x_here=x_here, color=color, lw=lw)

        #----------------------------------------------------
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            html = "<span style='color: %s; font-size: %s'> %s</span>" % (labels_color,
                                                                          labels_font_size,
                                                                          label)
            text = pyqtgraph.TextItem(text=" %s" % label, anchor=(0, 0.5), html=html)
            cp.graphWidget.addItem(text)
            # text.setFont(QtGui.QFont(None, 8))
            text.setPos(x_here, y_here)

        # Add label above the branch (optional)
        show_confidence = False
        def conf2str(conf):
            if int(conf) == conf:
                return str(int(conf))
            return str(conf)
        if show_confidence:
            def format_branch_label(clade):
                if hasattr(clade, 'confidences'):
                    # phyloXML supports multiple confidences
                    return '/'.join(conf2str(cnf.value)
                                    for cnf in clade.confidences)
                if clade.confidence:
                    return conf2str(clade.confidence)
                return None
        else:
            def format_branch_label(clade):
                return None
        conf_label = format_branch_label(clade)
        if conf_label:
            # cp.draw_label(conf_label, [0.5 * (x_start + x_here), y_here]) # TODO: update this too.
            pass
        #----------------------------------------------------

        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(use_linecollection=True, orientation='vertical',
                             x_here=x_here, y_bot=y_bot, y_top=y_top, color=color, lw=lw)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)
    #############################################################################

    draw_clade(tree.root, 0, 'k', 1)

    def get_plot_data_from_coords(coords_list):
        x_list = []
        y_list = []
        for coords in coords_list:
            x_list.extend([coords[0], coords[2]])
            y_list.extend([coords[1], coords[3]])
        return x_list, y_list

    # If line collections were used to create clade lines, here they are added
    # to the pyqtgraph plot.
    tree_pen = pyqtgraph.mkPen(width=tree_lines_width, color=tree_color) # Color and width of the tree.

    x_list, y_list = get_plot_data_from_coords(horizontal_linecollections)
    # for coords in horizontal_linecollections:
    #     cp.draw_line(coords)
    cp.graphWidget.plot(x_list, y_list, connect="pairs", pen=tree_pen)

    x_list, y_list = get_plot_data_from_coords(vertical_linecollections)
    cp.graphWidget.plot(x_list, y_list, connect="pairs", pen=tree_pen)


###############################################################################
# Builds a plot showing a dendrogram from MODELLER.                           #
###############################################################################

def draw_modeller_dendrogram(dendrogram_filepath, pymod,
                             tree_lines_width=1, tree_color="#000000",
                             labels_font_size="12px", labels_color="#404040"):
    """
    Draw 'dendrogram_file' for SALIGN tree file using PyMod plot engine based on
    pyqtgraph.
    """

    with open(dendrogram_filepath, 'r') as fp:
        content = fp.read()

    width = max([len(sline) for sline in content.splitlines()])
    height = len(content.splitlines())

    # Builds the plot window and initializes the plotting area.
    cp = PyMod_plot_window_qt(pymod.main_window)
    cp.initialize(pymod=pymod, title="MODELLER Dendrogram")
    cp.build_plotting_area(use_controls=False,
                           x_label_text="Branch length", y_label_text="Taxum",
                           hide_y_ticks=True, hide_x_ticks=True,
                           use_all_controls_buttons=False,
                           use_save_to_csv=False)
    cp.show()


    # Functions used to draw the dendrogram.
    x_list = []
    y_list = []

    def add_single_line_data(coords):
        x_list.extend([coords[0], coords[2]])
        y_list.extend([coords[1], coords[3]])

    def draw_single_label(label, coords):
        html = "<span style='color: %s; font-size: %s'> %s</span>" % (labels_color,
                                                                      labels_font_size,
                                                                      label)
        text = pyqtgraph.TextItem(text=" %s" % label, anchor=(0, 0.5), html=html)
        cp.graphWidget.addItem(text)
        text.setPos(coords[0], coords[1])

    def draw_lines_recursively(sline, y=0, offset=0):
        if not sline.strip():
            return

        if sline.lstrip().startswith('.-'): # leaf
            x1=sline.find('.-')
            x2=sline.find('- ')+1
            if not x2 or x1>x2:
                x2=x1
            add_single_line_data([offset+x1,y,offset+x2,y])
            draw_single_label(sline[x2:].strip(),[offset+x2+1,y])
        elif sline.lstrip().startswith('+-') and sline[-2:]=='-+':
            x1=sline.find('+-')
            x2=len(sline)-1
            add_single_line_data([offset+x1,y,offset+x2,y])
            for j,c in enumerate(sline):
                if c=='+':
                    x=j
                    y1=height-i-0.1
                    y2=height-i+0.1
                    add_single_line_data([offset+x,y1,offset+x,y2])
        elif sline.lstrip()[0] in '0123456789':
            x1=0
            for j,c in enumerate(sline):
                if j and sline[j-1]==' ':
                    x1=j*1
                elif j+1==len(sline):
                    draw_single_label(sline[x1:],[offset+x1,y])
                elif j+1<len(sline) and sline[j+1]==' ':
                    draw_single_label(sline[x1:j],[offset+x1,y])
        elif sline.lstrip().startswith('|'): # branch
            x=sline.find('|')
            y1=height-i-1
            y2=height-i+1
            add_single_line_data([offset+x,y1,offset+x,y2])
            if x+1<len(sline):
                draw_lines_recursively(sline[x+1:],y,offset=x+1+offset)


    # Actually draw the dendrogram.
    for i, sline in enumerate(content.splitlines()):
        draw_lines_recursively(sline, y=height-i)

    tree_pen = pyqtgraph.mkPen(width=tree_lines_width, color=tree_color) # Color and width of the tree.
    cp.graphWidget.plot(x_list, y_list, connect="pairs", pen=tree_pen)


###############################################################################
# Minimal example for pyqtgraph.                                              #
###############################################################################

if __name__ == "__main__":

    from pymol.Qt import QtWidgets
    from pymod_lib.pymod_plot.pyqtgraph import PlotWidget, plot
    import pymod_lib.pymod_plot.pyqtgraph as pg

    class Plot_window(QtWidgets.QMainWindow):

        def __init__(self, *args, **kwargs):
            super(Plot_window, self).__init__(*args, **kwargs)

            self.graphWidget = pg.PlotWidget()
            self.setCentralWidget(self.graphWidget)

            hour = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            temperature = [30, 32, 34, 32, 33, 31, 29, 32, 35, 45]

            # plot data: x, y values
            self.graphWidget.plot(hour, temperature)

    w = Plot_window(self.main_window)
    w.show()
