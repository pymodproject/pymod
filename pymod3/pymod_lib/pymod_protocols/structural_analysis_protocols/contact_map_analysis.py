# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing a protocol for building and analyzing contact and distance maps in PyMod.
"""

import os
import math
import warnings

import numpy as np
from pymol import cmd

from pymol.Qt import QtWidgets, QtGui, QtCore

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_combobox_qt)
from pymod_lib.pymod_protocols.structural_analysis_protocols._base_structural_analysis import get_distance
from pymod_lib.pymod_protocols.structural_analysis_protocols._base_structural_analysis import Structural_analysis_mixin
from pymod_lib.pymod_vars import viridis_colors
from pymod_lib.pymod_threading import Protocol_exec_dialog

viridis_colors_rev = list(reversed(viridis_colors))


###################################################################################################
# Main class.                                                                                     #
###################################################################################################

class Contact_map_analysis(PyMod_protocol, Structural_analysis_mixin):
    """
    Class for implementing in PyMod contact and distance map analyses.
    """

    protocol_name = "contact_map_analysis"

    def __init__(self, pymod):
        PyMod_protocol.__init__(self, pymod)
        self.target_sequences = self.get_pymod_elements()


    def launch_from_gui(self):
        """
        Checks the condition for launching a contact map analysis.
        """

        error_title = "Selection Error"

        # Check that the selection is not empty.
        if len(self.target_sequences) == 0:
            self.pymod.main_window.show_error_message(error_title,
                "Please select at least one structure to display its contact map.")
            return None

        # Check that no nucleic acid structure has been selected.
        if not self.check_protein_selection(error_title):
            return None

        # Check the right conditions for the protocol.
        if len(self.target_sequences) == 1: # Single structure.

            if not self.target_sequences[0].has_structure():
                self.pymod.main_window.show_error_message(error_title,
                    "Please select only elements with a 3D structure loaded in PyMOL to display their contact maps.")
                return None
            self.target_sequence = self.target_sequences[0]

        else: # Multiple structures.

            # Check for multiple aligned structures in the same cluster.
            if not self.check_multiple_selection(error_title):
                return None

        # Opens the options window.
        self.contact_map_options_window = Contact_map_options_window_qt(self.pymod.main_window,
            protocol=self,
            title="Contact Map Options",
            upper_frame_title="Options for Contact Map Analysis",
            submit_command=self.contact_map_analysis_state)
        self.contact_map_options_window.show()


    def validate_input(self):

        try:
            # Feature type.
            self.feature_type, self.feature_type_full_name = self.contact_map_options_window.get_feature_type()
            # Interaction center.
            self.interaction_center = self.contact_map_options_window.get_interaction_center()
            # Distance threshold.
            self.dist_threshold = self.contact_map_options_window.dist_threshold_enf.getvalue(validate=True)
            # Pixel size.
            self.pixel_size = 5
            # Reference structure.
            if len(self.target_sequences) > 1:
                self.reference_element = self.target_sequences[self.contact_map_options_window.get_reference_id()]
            else:
                self.reference_element = self.target_sequence
        except (ValueError, KeyError) as e:
            return str(e)

        return None


    def contact_map_analysis_state(self):

        #------------------------------------
        # Gets the parameters from the GUI. -
        #------------------------------------

        input_error = self.validate_input()
        if input_error is not None:
            self.pymod.main_window.show_error_message("Input Error", input_error)
            return None

        self.contact_map_options_window.destroy()


        #-------------------------------------
        # Actually computes the contact map. -
        #-------------------------------------


        # Assign function to fill the target map.
        if self.feature_type == "contact":
            self.get_map_pixel = self._get_contact_val

        elif self.feature_type == "distance":
            self.get_map_pixel = self._get_distance_val

        elif self.feature_type == "distances_difference":
            self.get_map_pixel = self._get_distances_difference_val

        elif self.feature_type == "distances_mean":
            self.get_map_pixel = self._get_distances_mean

        elif self.feature_type == "distances_std":
            self.get_map_pixel = self._get_distances_std

        else:
            raise KeyError(self.feature_type)


        if not self.pymod.use_protocol_threads:
            self.compute_target_map()
        else:
            label_text = ("Computing %s map on %s. Please wait for the process to"
                          " complete..." % (self.feature_type_full_name, self.get_seq_text(self.target_sequences)))
            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.compute_target_map,
                                            args=(), lock=True,
                                            wait_start=0.25, wait_end=0.25,
                                            title="Computing %s map" % self.feature_type_full_name,
                                            label_text=label_text)
            p_dialog.exec_()


        #------------------------
        # Show the contact map. -
        #------------------------

        # Gets the title of the window.
        if self.feature_type in ("contact", "distance"):
            title = "%s %s Map." % (self.target_sequence.my_header, self.feature_type_full_name[0:-1])
        elif self.feature_type == "distances_difference":
            title = "Distances Difference Map."
        elif self.feature_type == "distances_mean":
            title = "Distances Mean Map."
        elif self.feature_type == "distances_std":
            title = "Distances Std Map."

        # Initializes the contact/distance map window.
        cp = Contact_map_analysis_window_qt(self.pymod.main_window)
        cp.initialize_map(pymod=self.pymod,
                          data_array=self.target_map,
                          pymod_elements=self.target_sequences,
                          ref_residues=self.ref_residues,
                          ref_selectors=self.ref_selectors,
                          title=title,
                          pixel_size=self.pixel_size,
                          feature_type=self.feature_type,
                          threshold=self.dist_threshold,
                          interaction_center=self.interaction_center)
        cp.show()


    def compute_target_map(self):

        # Reference residues and selectors lists. Used to interact with the 3D structures loaded in
        # PyMOL.
        self.target_map = []
        self.ref_residues = []
        self.ref_selectors = []

        #----------------------------------------------------------------------------------
        # Builds the array with the contact/distance maps to show for a single structure. -
        #----------------------------------------------------------------------------------

        if self.feature_type in ("contact", "distance"):

            # Get the coordinates of the residues.
            residues, coords_array, selectors = self.get_coords_array(self.target_sequence, self.interaction_center)

            self.target_map = np.zeros((len(residues), len(residues)))
            for res_idx_i, res_i in enumerate(residues):
                for res_idx_j, res_j in enumerate(residues):
                    if res_idx_i >= res_idx_j:
                        continue
                    dist = get_distance(coords_array[res_idx_i], coords_array[res_idx_j])
                    self.target_map[res_idx_i][res_idx_j] = self.get_map_pixel(dist)
                    self.target_map[res_idx_j][res_idx_i] = self.target_map[res_idx_i][res_idx_j]

            self.ref_residues = residues
            self.ref_selectors = selectors

        #---------------------------------------------------------------------------
        # Builds the array with the distance maps to show for multiple structures. -
        #---------------------------------------------------------------------------

        elif self.feature_type in ("distances_difference", "distances_mean", "distances_std"):

            # Get the coordinates of the residues for all the structures.
            ali_seq_dict_list = []
            res_dict_list = []

            for element in self.target_sequences:
                ali_seq_dict, res_dict = self._get_sda_data(element, self.interaction_center)
                ali_seq_dict_list.append(ali_seq_dict)
                res_dict_list.append(res_dict)
                if element is self.reference_element:
                    ref_ali_seq_dict = ali_seq_dict
                    ref_res_dict = res_dict


            # Iter through all the positions of the alignment.
            ali_len = len(self.target_sequences[0].my_sequence)

            self.target_map = np.zeros((ali_len, ali_len))
            self.target_map[:] = -1 # Diagonal elements and alignment posisions with gaps will have -1 as a value.

            for ali_pos_i in range(0, ali_len):

                # Get the reference residues.
                ref_res = self._get_res(ali_pos_i, ref_ali_seq_dict)
                ref_sel = self._get_selector(ref_res, ref_res_dict)
                if ref_sel is None:
                    self.ref_residues.append(None)
                    self.ref_selectors.append(None)
                else:
                    self.ref_residues.append(ref_res)
                    self.ref_selectors.append(ref_sel)

                # Computes the distances for each pair of alignment columns.
                for ali_pos_j in range(0, ali_len):

                    if ali_pos_i >= ali_pos_j:
                        continue

                    # For each alignment column get the corresponding residues and coordinates.
                    dist_list = []
                    for ali_seq_dict, res_dict in zip(ali_seq_dict_list, res_dict_list):

                        coords_i = self._get_coords(self._get_res(ali_pos_i, ali_seq_dict), res_dict)
                        if coords_i is None:
                            continue
                        coords_j = self._get_coords(self._get_res(ali_pos_j, ali_seq_dict), res_dict)
                        if coords_j is None:
                            continue
                        dist_list.append(get_distance(coords_i, coords_j))

                    if len(dist_list) > 1:
                        self.target_map[ali_pos_i][ali_pos_j] = self.get_map_pixel(dist_list)
                        self.target_map[ali_pos_j][ali_pos_i] = self.target_map[ali_pos_i][ali_pos_j]


    def _get_contact_val(self, v):
        """
        Binary values.
        """
        return int(v <= self.dist_threshold)

    def _get_distance_val(self, v):
        """
        Continuous values.
        """
        return v

    def _get_distances_difference_val(self, v):
        """
        Absolute difference bewteen two distances.
        """
        if len(v) < 2:
            return -1
        else:
            return abs(v[0]-v[1])

    def _get_distances_mean(self, v):
        """
        Mean between a list of distances.
        """
        return get_mean(v)

    def _get_distances_std(self, v):
        """
        Standard deviation between a list of distances.
        """
        return get_std(v)

def get_mean(v):
    return sum(v)/len(v)

def get_std(v):
    n = len(v)
    if n == 1:
        return 0
    m = get_mean(v)
    return math.sqrt(sum([i**2 for i in v])/n - m**2)


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

###############################################################################
# Options window.                                                             #
###############################################################################

class Contact_map_options_window_qt(PyMod_protocol_window_qt):
    """
    Window with options for the contact map analysis.
    """

    feature_types_list = ["Contacts", "Distances", "Distance Difference", "Distance Mean", "Distance Std"]
    interaction_centers_list = ["Carbon alpha", "Carbon beta"]
    default_contact_threshold = 8.0
    default_distance_threshold = 18.0
    default_distance_diff_threshold = 3.5
    default_distance_mean_threshold = default_distance_threshold
    default_distance_std_threshold = default_distance_diff_threshold

    def build_protocol_middle_frame(self):

        # Feature type selection.
        if len(self.protocol.target_sequences) == 1: # Contact and distance.
            sel_option_idx = 0
            all_options_idx = (0, 1)
        elif len(self.protocol.target_sequences) == 2: # Distance difference, mean and std.
            sel_option_idx = 2
            all_options_idx = (2, 3, 4)
        elif len(self.protocol.target_sequences) > 2: # Distance mean and std.
            sel_option_idx = 3
            all_options_idx = (3, 4)

        features_buttons = [self.feature_types_list[i] for i in all_options_idx]
        self.feature_type_rds = PyMod_radioselect_qt(label_text="Feature Type",
                                                     buttons=features_buttons)
        for button, option_idx in zip(self.feature_type_rds.get_buttons(), all_options_idx):
            button.clicked.connect(lambda e=None, o=option_idx: self.feature_type_state(o))
        self.feature_type_rds.setvalue(self.feature_types_list[sel_option_idx])
        self.middle_formlayout.add_widget_to_align(self.feature_type_rds)


        # Distance threshold for contacts.
        if len(self.protocol.target_sequences) == 1:
            threshold = self.default_contact_threshold
        elif len(self.protocol.target_sequences) == 2:
             threshold = self.default_distance_diff_threshold
        else:
            threshold = self.default_distance_threshold

        self.dist_threshold_enf = PyMod_entryfield_qt(label_text="Contact Threshold (%s)" % ("\u212B"),
                                                      value=str(threshold),
                                                      validate={'validator': 'real',
                                                                'min': 1.0, 'max': 100.0})
        self.middle_formlayout.add_widget_to_align(self.dist_threshold_enf)


        # Interaction center selection.
        self.interaction_center_rds = PyMod_radioselect_qt(label_text="Interaction Center",
                                                           buttons=self.interaction_centers_list)
        self.interaction_center_rds.setvalue(self.interaction_centers_list[0])
        self.middle_formlayout.add_widget_to_align(self.interaction_center_rds)


        # Reference structure combobox.
        if len(self.protocol.target_sequences) > 1:
            structures_list = [element.my_header for element in self.protocol.target_sequences]
            self.reference_combobox = PyMod_combobox_qt(label_text="Reference Structure",
                                                        items=structures_list)
            self.reference_combobox.combobox.setCurrentIndex(0)
            self.middle_formlayout.add_widget_to_align(self.reference_combobox)


        self.middle_formlayout.set_input_widgets_width("auto", padding=10)


    def feature_type_state(self, state):
        """
        Launched when the user changes the "Feature Type" in the Options window.
        """
        # Change the status of the Tkinter checkbutton.
        if state == 0:
            self.feature_type_rds.setvalue(self.feature_types_list[0])
            self.dist_threshold_enf.setvalue(str(self.default_contact_threshold))
        elif state == 1:
            self.feature_type_rds.setvalue(self.feature_types_list[1])
            self.dist_threshold_enf.setvalue(str(self.default_distance_threshold))
        elif state == 2:
            self.feature_type_rds.setvalue(self.feature_types_list[2])
            self.dist_threshold_enf.setvalue(str(self.default_distance_diff_threshold))
        elif state == 3:
            self.feature_type_rds.setvalue(self.feature_types_list[3])
            self.dist_threshold_enf.setvalue(str(self.default_distance_mean_threshold))
        elif state == 4:
            self.feature_type_rds.setvalue(self.feature_types_list[4])
            self.dist_threshold_enf.setvalue(str(self.default_distance_std_threshold))
        else:
            KeyError(state)


    def get_feature_type(self):
        feature_type = self.feature_type_rds.getvalue()
        if feature_type == self.feature_types_list[0]:
            return "contact", feature_type
        elif feature_type == self.feature_types_list[1]:
            return "distance", feature_type
        elif feature_type == self.feature_types_list[2]:
            return "distances_difference", feature_type
        elif feature_type == self.feature_types_list[3]:
            return "distances_mean", feature_type
        elif feature_type == self.feature_types_list[4]:
            return "distances_std", feature_type
        else:
            raise KeyError(feature_type)


    def get_interaction_center(self):
        int_center_string = self.interaction_center_rds.getvalue()
        if int_center_string == self.interaction_centers_list[0]:
            return "ca"
        elif int_center_string == self.interaction_centers_list[1]:
            return "cb"
        else:
            raise KeyError(int_center_string)


    def get_reference_id(self):
        try:
            return self.reference_combobox.get_index()
        except:
            print("- WARNING: could not obtain the reference structure id.")
            return 0


###############################################################################
# Results window.                                                             #
###############################################################################

class Contact_map_graphics_view(QtWidgets.QGraphicsView):
    pass


class Contact_map_pixel(QtWidgets.QGraphicsRectItem):
    """
    Custom class for drawing on a graphic scene of PyQt rectangles corresponding
    to pixels of a contact/distance map.
    """

    def __init__(self, x, y, w, h, contact_map_window, pen, brush, i, j):

        super(Contact_map_pixel, self).__init__(x, y, w, h)
        self.contact_map_window = contact_map_window
        self.setPen(pen)
        self.original_brush = brush
        self.setBrush(brush)
        self.setAcceptHoverEvents(True)
        self.i_idx = i
        self.j_idx = j


    def mousePressEvent(self, e):
        self.contact_map_window.canvas_plot_left_click_event(self.i_idx, self.j_idx)

    def hoverEnterEvent(self, e):
        self.setBrush(self.contact_map_window.highlight_brush)
        self.contact_map_window.canvas_plot_move_event(self.i_idx, self.j_idx)

    def hoverLeaveEvent(self, e):
        self.setBrush(self.original_brush)


class Contact_map_analysis_window_qt(QtWidgets.QMainWindow):
    """
    Class for a window containing a PyQt widget in which a contact/distance map
    will be drawn. Minimal PyQt implementation of a graphical contact map. NumPy
    is required.
    Note: for large protein, building every rectangle in the canvas takes a lot
    of time and there is room for optimization.
    """

    is_pymod_window = True

    distance_count = 0


    def initialize_map(self, pymod, data_array,
                       pymod_elements,
                       ref_residues,
                       ref_selectors,
                       title=None,
                       pixel_size=5,
                       feature_type="contact",
                       threshold=8.0,
                       interaction_center="ca"):

        # Sets the attributes.
        self.data_array = data_array
        self.pixel_size = pixel_size
        self.feature_type = feature_type
        self.threshold = threshold
        self.interaction_center = interaction_center
        if self.feature_type in ("contact", "distance"):
            self.pymod_elements = pymod_elements
            self.pymod_element = self.pymod_elements[0]
        else:
            self.pymod_elements = pymod_elements
        # Get the PyMod residues for each residue having an interaction center and the PyMOL selectors
        # for each residue.
        self.ref_residues = ref_residues
        self.ref_selectors = ref_selectors


        # Assign the methods to get the labels.
        if self.feature_type == "contact":
            self.get_value_label = self._get_value_label_contact
        elif self.feature_type == "distance":
            self.get_value_label = self._get_value_label_distance
        elif self.feature_type == "distances_difference":
            self.get_value_label = self._get_value_label_distance_diff
        elif self.feature_type == "distances_mean":
            self.get_value_label = self._get_value_label_distance_mean
        elif self.feature_type == "distances_std":
            self.get_value_label = self._get_value_label_distance_std
        else:
            raise KeyError(self.feature_type)

        # Set the canvas size.
        min_size = 150
        h = self.pixel_size*len(self.data_array)
        win_size = min((910, h))
        win_size = max((min_size, win_size))

        if title:
            self.setWindowTitle(title)


        # Set some appearance parameters.
        self.controls_padding = 4
        if self.feature_type in ("contact", "distance"):
            self.controls_font = "helvetica 11 bold"
        else:
            self.controls_font = "helvetica 10 bold"
        self.controls_config = {"fg": "black", "font": self.controls_font,
                                "padx": self.controls_padding,
                                "pady": self.controls_padding}
        self.labels_pack_config = {"side": "left", "pady": (0, 5), "padx": (5, 0)}
        self.buttons_pack_config = {"side": "left", "pady": (0, 5), "padx": (1, 0)}


        # Frame of the window containing a row for some control buttons, a row for
        # the plot and a row for a messagebar.
        self.plot_frame = QtWidgets.QWidget()
        self.plot_frame_layout = QtWidgets.QGridLayout()
        self.plot_frame.setLayout(self.plot_frame_layout)
        self.setCentralWidget(self.plot_frame)


        # Control frame.
        self.controls_frame = QtWidgets.QWidget()
        self.controls_frame_layout = QtWidgets.QGridLayout()
        self.controls_frame.setLayout(self.controls_frame_layout)
        self.plot_frame_layout.addWidget(self.controls_frame)

        self.delete_distances_button = QtWidgets.QPushButton("Delete all distances in PyMOL")
        self.delete_distances_button.setEnabled(False)
        self.delete_distances_button.clicked.connect(lambda a=None: self.clear_plot())
        self.controls_frame_layout.addWidget(self.delete_distances_button, 0, 0)

        self.scale_factor = 0
        self.scale_down_button = QtWidgets.QPushButton("Zoom out")
        try:
            self.scale_down_button.setIcon(QtGui.QIcon.fromTheme("go-down"))
        except:
            pass
        self.scale_down_button.clicked.connect(lambda a=None: self.scale_plot_down())
        self.controls_frame_layout.addWidget(self.scale_down_button, 0, 1)

        self.scale_up_button = QtWidgets.QPushButton("Zoom in")
        try:
            self.scale_up_button.setIcon(QtGui.QIcon.fromTheme("go-up"))
        except:
            pass
        self.scale_up_button.clicked.connect(lambda a=None: self.scale_plot_up())
        self.controls_frame_layout.addWidget(self.scale_up_button, 0, 2)

        self.controls_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


        # Frame containing the plot (with a scrollbar).
        self.canvas_plot_frame = QtWidgets.QWidget()
        self.canvas_plot_frame.setStyleSheet("background-color: white")
        self.canvas_plot_frame_layout = QtWidgets.QGridLayout()
        self.canvas_plot_frame.setLayout(self.canvas_plot_frame_layout)

        self.canvas_plot_scrollarea = QtWidgets.QScrollArea()
        self.canvas_plot_scrollarea.setWidgetResizable(True)
        self.canvas_plot_scrollarea.setWidget(self.canvas_plot_frame)
        self.plot_frame_layout.addWidget(self.canvas_plot_scrollarea)


        # Builds the scene where to draw the contact map.
        self.canvas_plot_scene = QtWidgets.QGraphicsScene()
        # Builds the graphics view containing the scene above.
        self.canvas_plot_view = Contact_map_graphics_view(self.canvas_plot_scene)
        self.canvas_plot_frame_layout.addWidget(self.canvas_plot_view)


        # A bottom frame fo the window, containing some buttons to interact with the graph.
        self.message_frame = QtWidgets.QFrame()
        self.message_frame_layout = QtWidgets.QHBoxLayout()
        self.message_frame.setLayout(self.message_frame_layout)
        self.plot_frame_layout.addWidget(self.message_frame)

        # Label to show which residue/position pair is currently being hovered by the mouse pointer.
        if self.feature_type in ("contact", "distance"):
            view_label_text = "Couple:"
        else:
            view_label_text = "Alignment positions:"
        self.view_label = QtWidgets.QLabel(view_label_text)
        # self.view_label.setStyleSheet(self.controls_config)
        self.message_frame_layout.addWidget(self.view_label)


        # Actually draws the contact map.
        self.draw_map()
        # self._test_plot()


    def draw_map(self):
        """
        Methods that actually draw the contact/distance map on a canvas widget.
        """
        w = self.pixel_size

        # Prepare the brushes.
        self.viridis_brushes = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for c in viridis_colors_rev:
                brush = QtGui.QBrush(QtGui.QColor(round(c[0]*255), round(c[1]*255), round(c[2]*255)))
                self.viridis_brushes.append(brush)
        self.default_brush = QtGui.QBrush(QtGui.QColor(242, 242, 242))
        self.highlight_brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))

        # Function to get the color of pixels.
        if self.feature_type == "contact":
            self.get_color = self._get_color_contact

        elif self.feature_type in ("distance", "distances_mean"):
            self._dist_bins = np.linspace(2.5, self.threshold, 63)
            self.get_color = self._get_color_distance

        elif self.feature_type in ("distances_difference", "distances_std"):
            self._dist_bins = np.linspace(0.0, self.threshold, 63)
            self.get_color = self._get_color_distance

        else:
            raise KeyError(self.feature_type)


        #-------------------------------
        # Draws the map on the canvas. -
        #-------------------------------

        # Use a transparent pen for the border pf the pixels.
        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 0), 0)

        for i in range(0, len(self.data_array)):

            for j in range(0, len(self.data_array)):

                if i <= j:

                    # Gets the color brush.
                    color_brush = self.get_color(self.data_array[i][j])

                    # Builds rectangles for both the upper and lower part of the
                    # matrix and adds them to the graphics scene.
                    pid = Contact_map_pixel(j*w, i*w, w, w, i=i, j=j,
                                            contact_map_window=self,
                                            pen=pen, brush=color_brush)
                    self.canvas_plot_scene.addItem(pid)

                    pid = Contact_map_pixel(i*w, j*w, w, w, i=j, j=i,
                                            contact_map_window=self,
                                            pen=pen, brush=color_brush)
                    self.canvas_plot_scene.addItem(pid)


        # Add events to the canvas.
        if self.feature_type in ("contact", "distance"):
            # Draws the map of a single structure on the canvas.
            self.canvas_plot_move_event = self.move_on_plot
            self.canvas_plot_left_click_event = self.click_on_plot

        else:
            # Draws the map of multiple structures on the canvas.
            self.canvas_plot_move_event = self.move_on_plot_ali
            self.canvas_plot_left_click_event = self.click_on_plot_ali


    # Click on the plot.
    def click_on_plot(self, i, j):
        """
        When clicking on a square, draw a distance in the PyMOL viewer between the corresponding
        residues. Used when analyzing a single structure.
        """
        res_1 = self.ref_residues[i]
        res_2 = self.ref_residues[j]

        if res_1 is res_2:
            return None

        sel_1 = self.ref_selectors[i]
        sel_2 = self.ref_selectors[j]
        cmd.distance("pymod_dist_%s" % self.distance_count, sel_1, sel_2)
        # cmd.center("pymod_distance_%s" % self.distance_count)
        self.distance_count += 1
        self.delete_distances_button.setEnabled(True)


    def click_on_plot_ali(self, i, j):
        """
        Used when analyzing multiple structures. The residues of the reference structures will be
        used in PyMOL.
        """
        if i == j:
            return None

        sel_1 = self.ref_selectors[i]
        sel_2 = self.ref_selectors[j]
        if sel_1 is None or sel_2 is None:
            return None
        cmd.distance("pymod_dist_%s" % self.distance_count, sel_1, sel_2)
        self.distance_count += 1
        self.delete_distances_button.setEnabled(True)


    # Move the mouse on the plot.
    def move_on_plot(self, i, j):
        """
        Used when showing the contact/distance map of a single structure.
        """
        val = self.data_array[i][j]
        res_1 = self.ref_residues[i]
        res_2 = self.ref_residues[j]
        self.view_label.setText("Couple: %s %s - %s %s (%s)." % (res_1.three_letter_code, res_1.db_index,
                                                                 res_2.three_letter_code, res_2.db_index,
                                                                 self.get_value_label(val)))

    def move_on_plot_ali(self, i, j):
        """
        Used when showing the distance map of multiple structures.
        """
        val = self.data_array[i][j]
        res_1 = self.ref_residues[i]
        res_2 = self.ref_residues[j]

        label = "Alignment positions: %s - %s (%s)." % (i+1, j+1, self.get_value_label(val))
        if res_1 is not None and res_2 is not None:
            label += " Reference: %s %s - %s %s" % (res_1.three_letter_code, res_1.db_index,
                                                    res_2.three_letter_code, res_2.db_index)
        self.view_label.setText(label)


    # Get the labels to show on the bottom of the window.
    def _get_value_label_contact(self, value):
        if value == 1:
            return "contact"
        else:
            return "non contact"

    def _get_value_label_distance(self, value):
        return str(round(value, 2)) + " \u212B"

    def _get_value_label_distance_diff(self, value):
        return self._get_value_label_distance_ali(value, "diff")

    def _get_value_label_distance_mean(self, value):
        return self._get_value_label_distance_ali(value, "mean")

    def _get_value_label_distance_std(self, value):
        return self._get_value_label_distance_ali(value, "std")

    def _get_value_label_distance_ali(self, value, label):
        if value == -1:
            return "not aligned"
        else:
            return "%s = %s %s" % (label, round(value, 2), " \u212B")


    # Get the colors for pixels in the map.
    def _get_color_contact(self, v):
        if v == 0:
            return self.viridis_brushes[-1] # "white"
        else:
            return self.viridis_brushes[0] # "gray"

    def _get_color_distance(self, v):
        """
        Assignes the bin for the distance value and returns the corresponding color.
        """
        if v == -1:
            return self.default_brush
        return self.viridis_brushes[np.digitize(v, self._dist_bins)]


    # Events influecing the whole plot.
    def clear_plot(self):
        """
        Remove all distances which have been drawn.
        """
        cmd.delete("pymod_dist_*")
        self.delete_distances_button.setEnabled(False)

    def scale_plot_up(self):
        if self.scale_factor > 10:
            return None
        self.canvas_plot_view.scale(1.25, 1.25)
        self.scale_factor += 1
        if self.scale_factor > 10:
            self.scale_up_button.setEnabled(False)
        self.scale_down_button.setEnabled(True)

    def scale_plot_down(self):
        if self.scale_factor < -10:
            return None
        self.canvas_plot_view.scale(0.8, 0.8)
        self.scale_factor -=1
        if self.scale_factor < -10:
            self.scale_down_button.setEnabled(False)
        self.scale_up_button.setEnabled(True)
