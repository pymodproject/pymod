# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import math
import time
import warnings

import pymod_lib.pymod_seq.seq_manipulation as pmsm
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol

from pymod_lib.pymod_vars import prot_standard_one_letter, prot_standard_three_letters


class Ramachandran_plot(PyMod_protocol):

    protocol_name = "ramachandran_plot"

    def __init__(self, pymod, target_sequences=None):
        PyMod_protocol.__init__(self, pymod)
        self.target_sequences = self.get_pymod_elements(target_sequences)


    def launch_from_gui(self):
        """
        PROCHEK style Ramachandran Plot.
        """
        if not len(self.target_sequences) == 1 or not self.target_sequences[0].has_structure():
            self.pymod.main_window.show_error_message("Selection Error", "Please select one structure to display its Ramachandran Plot.")
            return None

        if True in [e.polymer_type == "nucleic_acid" for e in self.target_sequences]:
            self.pymod.main_window.show_error_message("Selection Error", "Can not build a Ramachandran plot for nucleic acids structures.")
            return None

        if not len(str(self.target_sequences[0].my_sequence).replace('-','')):
            self.pymod.main_window.show_error_message("Selection Error", "No residue for Ramachandran Plot generation")
            return None

        self.target_sequence = self.target_sequences[0]

        self.title = self.target_sequence.my_header
        self.PDB_file = [self.target_sequence.get_structure_file(basename_only=False)]

        # Show an options window to choose kind of aa to plot.
        self.options_window = Ramachandran_plot_options_window_qt(
            parent=self.pymod.main_window,
            protocol=self,
            title="Ramachandran Plot Options",
            upper_frame_title="Options for Ramachandran Plot",
            submit_command=self.options_window_state)
        self.options_window.show()


    def options_window_state(self):

        #----------------
        # Checks input. -
        #----------------

        aa_list = None

        if self.options_window.get_aa_selection_mode() == "single":
            aa_list = ''
            for aa in prot_standard_one_letter:
                if self.options_window.aa_checkbutton[aa].isChecked():
                    aa_list += aa
            if aa_list:
                self.title += " (Amino Acids: " + aa_list + ")"
            else:
                self.pymod.main_window.show_error_message("Selection Error",
                    "No amino acids for Ramachandran Plot generation. Please select some amino acids.")
                return

        self.options_window.destroy()

        time.sleep(0.25)


        #---------------------------------------------------------------------
        # Extract the data from the 3D structures. This is performed here in -
        # order to quit the protocol if something goes wrong.                -
        #---------------------------------------------------------------------

        self.input_filepath = self.target_sequence.get_structure_file(basename_only=False)
        PDB_file = self.input_filepath

        if isinstance(PDB_file, str):
            PDB_file = [PDB_file]

        # Parses the PDB structure and gets the phi and psi angles.
        residues_tags = []

        try:
            for structure_idx, filepath in enumerate(PDB_file):

                structure = PDB.PDBParser(QUIET=True).get_structure("structure_%s" % structure_idx, filepath)

                for model in structure:
                    for chain in model:

                        polypeptides = PDB.CaPPBuilder().build_peptides(chain)
                        for poly_idx,poly in enumerate(polypeptides):

                            phi_psi = poly.get_phi_psi_list()
                            for res_idx, residue in enumerate(poly):

                                residue_tags = {}

                                if aa_list != None and not code_standard[residue.resname] in aa_list:
                                    continue
                                if not residue.resname in code_standard:
                                    continue

                                phi, psi = phi_psi[res_idx]

                                if not (phi and psi):
                                    residue_tags["region"] = "end_res"
                                    residues_tags.append(residue_tags)
                                    continue

                                het, resseq, icode = residue.id
                                phi_str = "%.2f" % (phi/math.pi*180)
                                psi_str = "%.2f" % (psi/math.pi*180)

                                # key = residue.resname + str(resseq)
                                residue_tags.update({"resname": residue.resname,
                                                     "position": str(resseq),
                                                     "phi": phi, "psi": psi,
                                                     "phi_str": phi_str, "psi_str": psi_str})

                                # Glycines.
                                if residue.resname == "GLY":
                                    residue_tags["region"] = "gly"

                                # Prolines.
                                elif residue.resname == "PRO":
                                    residue_tags["region"] = "pro"

                                # Other residues.
                                else:
                                    region = procheck[int(18-psi/math.pi*18)][int(phi/math.pi*18+18)]
                                    residue_tags["region"] = region

                                residues_tags.append(residue_tags)

        except Exception as e:
            title = "Error"
            message = ("The Ramachandran plot could not be computed beacause of the"
                       " following error: %s" % str(e))
            self.pymod.main_window.show_error_message(title, message)
            return None


        #------------------------
        # Show the plot window. -
        #------------------------

        ramachandra_plot_window = Ramachandran_plot_window_qt(self.pymod.main_window)
        ramachandra_plot_window.initialize_plot(pymod=self.pymod,
                                                target_element=self.target_sequence,
                                                residues_tags=residues_tags,
                                                plot_title=self.title,
                                                aa_list=aa_list)
        ramachandra_plot_window.show()


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

from pymol.Qt import QtWidgets, QtCore, QtGui

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_combobox_qt,
                                                          small_font_style)


class Ramachandran_plot_options_window_qt(PyMod_protocol_window_qt):

    def build_protocol_middle_frame(self):
        """
        Allow to choose between plotting all amino acids types or only a subset
        of them.
        """

        # Radioselect.
        aa_select_choices = ("Use all amino acids", "Select amino acids types")
        aa_select_choices_values = ["all", "single"]
        self.aa_select_choices_dict = dict([(k, v) for (k, v) in zip(aa_select_choices, aa_select_choices_values)])
        self.aa_select_rds = PyMod_radioselect_qt(label_text="Select Amino Acids",
                                                  buttons=aa_select_choices)
        self.aa_select_rds.setvalue(aa_select_choices[0])
        self.aa_select_rds.buttons_dict[aa_select_choices[0]].clicked.connect(self.hide_select_single_aa_frame)
        self.aa_select_rds.buttons_dict[aa_select_choices[1]].clicked.connect(self.show_select_single_aa_frame)
        self.middle_formlayout.add_widget_to_align(self.aa_select_rds)

        # Checkboxes for selecting single amino acids.
        self.aa_select_grid_layout = QtWidgets.QGridLayout()
        self.aa_select_rds.input.addLayout(self.aa_select_grid_layout)

        self.aa_checkbutton = {}
        self.aa_freq_dict = {}
        self.aa_checkbutton_list = []

        for i, aa in enumerate(prot_standard_one_letter):

            # Get the number of aa in the sequence.
            aa_freq = str(self.protocol.target_sequence.my_sequence).count(aa)
            self.aa_freq_dict[aa] = aa_freq

            # Build a checkbox for the aa.
            checkbox = QtWidgets.QCheckBox(pmsm.one2three(aa) + " (" + str(aa_freq) + ")")
            checkbox.setEnabled(False)
            self.aa_select_grid_layout.addWidget(checkbox, int(i%10), int(i/10))

            self.aa_checkbutton[aa] = checkbox
            self.aa_checkbutton_list.append(checkbox)

        self.aa_select_grid_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.middle_formlayout.set_input_widgets_width("auto", padding=40)


    def show_select_single_aa_frame(self):
        for aa in prot_standard_one_letter:
            # Only enable selection of aa present in protein sequence.
            if self.aa_freq_dict[aa] != 0:
                self.aa_checkbutton[aa].setEnabled(True)

    def hide_select_single_aa_frame(self):
        for checkbox in self.aa_checkbutton_list:
            checkbox.setEnabled(False)

    def get_aa_selection_mode(self):
        return self.aa_select_choices_dict[self.aa_select_rds.getvalue()]


###################################################################################################
# INTERACTIVE RAMACHANDRAN PLOT.                                                                  #
###################################################################################################

from Bio import PDB

from pymol import cmd

# Global variable for Ramachandran Plots
procheck=[ # 4 regions defined by PROCHECK
    "AFFFFFFFFFFAAAGGDDDDDGGGGGDDDDDGGGGA", # F - Favored
    "AFFFFFFFFFFFAAAGGDDDDDDDDDDDDDDGGGAA", # A - Additional allowed
    "AFFFFFFFFFFFAAAGGGGDDDDDDDDDDDDGGAAA", # G - Generously allowed
    "AAFFFFFFFFFFFAAAAGGDDDDDDDDDDDDGGGAA", # D - Disallowed
    "AAFFFFFFFFFFFFAAGGGDDDDDDDDDDDDGGGGA",
    "AAFFFFFFFFFFFAAAGGGDDDDDDDDDDDDDGGGA",
    "AAAFFFFFFFFFAAAAGGGDDDGGGGGGDDDDDGGA",
    "AAAAFFFFFFFAAAAAAGGGGGGGGGGGDDDDDGGA",
    "AAAAAFAAFAAAAAAAAGGGGGGGAAGGDDDDDGGA",
    "AAAAAAAAAAAAAAGGGGGGGAAAAGGGDDDDDGGG",
    "AAAAAAAAAAAAGGGGGGGGGAAAAGGGDDDDDGGG",
    "GAAAAAAAAAAAAGGDDDDGGGAAAGGGGDDDDGGG",
    "GGAAAAAAAAAAAGGDDDDGGAAAAAAGGDDDDGGG",
    "GGAAAAAAAAAAGGGDDDDGGAAFAAGGGDDDDGGG",
    "GAAAAAAAAAAAGGGDDDDGGGAFFAGGGDDDDGGG",
    "GAAAAAAFAAAAAGGGDDDGGGAAAAAGGDDDDGGG",
    "GAAAAAFFFFAAAGGGGDDDGGGAAAGGGDDDDGGG",
    "GAAAAAFFFFFAAAGGGDDDGGGGAAAGGDDDDGGG",
    "GAAAAFFFFFFFAAAGGGDDGGGAGAAGGDDDDGGG",
    "GAAAAAFFFFFFFAAGGGGDGGGGGGGGGDDDDGGG",
    "GGAAAAFFFFFFFFAAGGGDGGGGGGGGGDDDDGGG",
    "GGAAAAAFFFFFFFAAAGGGDDDDDDDDDDDDDGGG",
    "GGGAAAAAFFFFFFFAAGGGDDDDDDDDDDDDDGGG",
    "GAAAAAAAAFFFFFFAAAGGGDDDDDDDDDDDDGGG",
    "AAGAAAAAAAAFFFFAAAGGGDDDDDDDDDDDDGGG",
    "GGGAAAAAAAAAAAAAAAAGGDDDDDDDDDDDDGGG",
    "GGGGGAAAAAAAAAAAAAGGGDDDDDDDDDDDDGGG",
    "DGGGGAAAAAAAAAAGGGGGGDDDDDDDDDDDDDDD",
    "DDDGGGGGGGGAGGGGGGGGDDDDDDDDDDDDDDDD",
    "DDDGGGGAAGGGGGGGGDDDDDDDDDDDDDDDDDDD",
    "GGGGGAAGGGGGGGDDDDDDDGGGGGDDDDDDDDDD",
    "GGGGGGAAAAGGGGDDDDDDDGGGGGDDDDDDDDDD",
    "GAAAGAAAAAGGGGGDDDDDDGGAGGGDDDDDDDDD",
    "GAAAAAAAAAAGGGGDDDDDDGGAGGGDDDDDDGGG",
    "GAAAAAAAAAAAAGGDDDDDDGGAAGGDDDDDDGGG",
    "AAAAAAAAAAAAAGGDDDDDDGGGGGGDDDDDDGGA"]


code_standard = { # dict for 20 standard amino acids
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'}


class Ramachandran_plot_items_mixin:

    def common_initialization(self, i, j, tags, parent_window, pen, brush):
        self.parent_window = parent_window
        self.plot_tags = tags
        self.set_paiting_tools(pen, brush)
        self.setAcceptHoverEvents(True)
        self.i_val = i
        self.j_val = j

    def set_paiting_tools(self, pen, brush):
        if pen is None:
            pen = self.parent_window.default_pen
        self.setPen(pen)
        if brush is None:
            brush = self.parent_window.default_brush
        self.original_brush = brush
        self.setBrush(brush)


    def hoverEnterEvent(self, e):
        self.parent_window.canvas_plot_move_event(self.plot_tags)
        self.setBrush(self.parent_window.highlight_brush)

    def hoverLeaveEvent(self, e):
        self.setBrush(self.original_brush)
        self.parent_window.canvas_plot_title_reset(self.plot_tags, "resname")

    def highlight_region(self):
        self.setBrush(self.parent_window.highlight_region_brush)

    def inactivate_region(self):
        self.setBrush(self.original_brush)


    def mousePressEvent(self, e):
        if e.button() == QtCore.Qt.LeftButton:
            self.parent_window.click_residue_with_left_button(self.plot_tags)
        elif e.button() == QtCore.Qt.MiddleButton:
            self.parent_window.click_residue_with_middle_button(self.plot_tags)
        elif e.button() == QtCore.Qt.RightButton:
            self.parent_window.click_residue_with_right_button(self.plot_tags)

    def mouseReleaseEvent(self, e):
        if e.button() == QtCore.Qt.LeftButton:
            self.parent_window.release_residue_with_left_button(self.plot_tags)


class Ramachandran_plot_triangle(QtWidgets.QGraphicsPolygonItem, Ramachandran_plot_items_mixin):

    def __init__(self, i, j, tags, mark_size, parent_window, pen=None, brush=None):

        # Build the polygon object representing a triangle.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            polygon = QtGui.QPolygonF([QtCore.QPointF(i, j-2*mark_size),
                                       QtCore.QPointF(i-1.7*mark_size, j+mark_size),
                                       QtCore.QPointF(i+1.7*mark_size, j+mark_size)])
        # Initialize.
        super(Ramachandran_plot_triangle, self).__init__(polygon)
        self.common_initialization(i, j, tags, parent_window, pen, brush)


class Ramachandran_plot_rectangle(QtWidgets.QGraphicsRectItem, Ramachandran_plot_items_mixin):

    def __init__(self, i, j, tags, mark_size, parent_window, pen=None, brush=None):
        super(Ramachandran_plot_rectangle, self).__init__(i-mark_size, j-mark_size, mark_size*2, mark_size*2)
        self.common_initialization(i, j, tags, parent_window, pen, brush)


class Ramachandran_plot_circle(QtWidgets.QGraphicsEllipseItem, Ramachandran_plot_items_mixin):

    def __init__(self, i, j, tags, mark_size, parent_window, pen=None, brush=None):
        super(Ramachandran_plot_circle, self).__init__(i-mark_size, j-mark_size, mark_size*2, mark_size*2)
        self.common_initialization(i, j, tags, parent_window, pen, brush)


class Ramachandran_plot_info_labels(QtWidgets.QLabel):

    def __init__(self, text, message, region, active, parent_window):
        super(Ramachandran_plot_info_labels, self).__init__(text)
        self.message_to_show = message
        self.plot_region = region
        self.activate_plot = active
        self.parent_window = parent_window

    def enterEvent(self, e):
        if not self.activate_plot:
            return None
        self.parent_window.info_label_move_event(self.plot_region, self.message_to_show)

    def leaveEvent(self, e):
        if not self.activate_plot:
            return None
        self.parent_window.info_label_title_reset(self.plot_region)


class Ramachandran_plot_window_qt(QtWidgets.QMainWindow):
    """
    PyQt class for a window containing a PROCHECK style Ramachandran plot.
    Biopython is required.
    """

    is_pymod_window = True

    def initialize_plot(self, pymod, target_element, residues_tags, plot_title, aa_list):

        self.pymod = pymod
        self.target_element = target_element
        self.residues_tags = residues_tags
        self.plot_title = plot_title
        self.setWindowTitle("%s Ramachandran Plot" % self.target_element.my_header)
        self.aa_list = aa_list


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
        self.plot_frame_layout.addWidget(self.controls_frame, 0, 0)

        self.scale_factor = 0
        self.scale_down_button = QtWidgets.QPushButton("Zoom out")
        try:
            self.scale_down_button.setIcon(QtGui.QIcon.fromTheme("go-down"))
        except:
            pass
        self.scale_down_button.clicked.connect(lambda a=None: self.scale_plot_down())
        self.controls_frame_layout.addWidget(self.scale_down_button, 0, 0)

        self.scale_up_button = QtWidgets.QPushButton("Zoom in")
        try:
            self.scale_up_button.setIcon(QtGui.QIcon.fromTheme("go-up"))
        except:
            pass
        self.scale_up_button.clicked.connect(lambda a=None: self.scale_plot_up())
        self.controls_frame_layout.addWidget(self.scale_up_button, 0, 1)

        self.controls_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


        # Frame containing the plot (with a scrollbar).
        self.canvas_plot_frame = QtWidgets.QWidget()
        self.canvas_plot_frame.setStyleSheet("background-color: white")
        self.canvas_plot_frame_layout = QtWidgets.QGridLayout()
        self.canvas_plot_frame.setLayout(self.canvas_plot_frame_layout)

        self.canvas_plot_scrollarea = QtWidgets.QScrollArea()
        self.canvas_plot_scrollarea.setWidgetResizable(True)
        self.canvas_plot_scrollarea.setWidget(self.canvas_plot_frame)
        self.plot_frame_layout.addWidget(self.canvas_plot_scrollarea, 1, 0)

        self.default_pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255), 1)
        self.default_brush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 230))
        self.highlight_brush = QtGui.QBrush(QtGui.QColor(255, 0, 255))
        self.highlight_region_brush = QtGui.QBrush(QtGui.QColor(0, 255, 255))


        # Builds the scene where to draw the Ramachandran plot.
        self.canvas_plot_scene = QtWidgets.QGraphicsScene()
        # Builds the graphics view containing the scene above.
        self.canvas_plot_view = QtWidgets.QGraphicsView(self.canvas_plot_scene)
        self.canvas_plot_frame_layout.addWidget(self.canvas_plot_view)


        # A bottom frame fo the window, containing some buttons to interact with the graph.
        self.message_frame = QtWidgets.QFrame()
        self.message_frame_layout = QtWidgets.QFormLayout()
        self.message_frame.setLayout(self.message_frame_layout)
        self.plot_frame_layout.addWidget(self.message_frame, 1, 1)

        # Label to show which residue/position pair is currently being hovered by the mouse pointer.
        self.view_label_prefix = "Showing: "
        self.default_message = "Hover dots to color residues of the same type" # Hover over title.
        self.message_label = QtWidgets.QLabel(self.default_message)
        self.message_label.setStyleSheet(small_font_style)
        self.message_frame_layout.addRow(self.message_label)


        # Actually draws the plot.
        self.draw_plot()


        # Shows some data about the type of residues and their dihedral angles.
        self.regions_labels_dict = {}
        tot_regular_res = float(self.residues_count["T"])

        label_params = [("F", "Residues in the most favoured regions", "residues in most favoured regions", True, True),
                        ("A", "Residues in additional allowed regions", "residues in additional allowed regions", True, True),
                        ("G", "Residues in generously allowed regions", "residues in generously allowed regions", True, True),
                        ("D", "Residues in disallowed regions", "residues in disallowed regions", True, True),
                        ("T", "Non-gly and non-pro residues (circles)", "non-glycine and non-proline residues", True, True),
                        ("end_res", "End-residues", "end residues", False, False),
                        ("gly", "Gly residues (triangles)", "glycine residues", True, False),
                        ("pro", "Pro residues (squares)", "proline residues", True, False),
                        ("total", "Total number of residues", "all residues", True, False)]
        for region, label, message, active, use_ratio in label_params:
            region_label = Ramachandran_plot_info_labels(label, message, region, active, self)
            region_label.setStyleSheet(small_font_style)
            if use_ratio:
                text = "%s (%s%%)" % (self.residues_count[region], round(self.residues_count[region]/tot_regular_res*100, 1))
            else:
                text = str(self.residues_count[region])
            region_label_count = QtWidgets.QLabel(text)
            region_label_count.setStyleSheet(small_font_style)
            self.message_frame_layout.addRow(region_label, region_label_count)
            self.regions_labels_dict[region] = {"info": region_label, "data": region_label_count}


    def draw_plot(self):

        self.pymol_selection = self.target_element.get_pymol_selector()

        imin = 60
        jmin = 40
        imax = imin+360
        jmax = jmin+360

        mark_size = 2

        xticks_num = [-180,-135,-90,-45,0,45,90,135,180]
        yticks_num = [-180,-135,-90,-45,0,45,90,135,180]
        height = 10
        width = 10


        #--------------------------------------------------------------
        # Draws the background and the axes of the Ramachandran plot. -
        #--------------------------------------------------------------

        # Draws the background of the Ramachandran plot.
        line_pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255), 1)

        for ii in range(0,36):

            for jj in range(0,36):

                region = procheck[ii][jj]
                color = "#ffffff" #[ 1.0, 1.0, 1.0]
                if region == 'F':
                    color = "#f20000" #[.949, 0.0, 0.0]
                elif region == 'A':
                    color = "#f2f200" #[.949,.949, 0.0]
                elif region == 'G':
                    color = "#f2f2a9" #[.949,.949,.663]

                qcolor = QtGui.QColor(0, 0, 0)
                qcolor.setNamedColor(color)
                brush = QtGui.QBrush(qcolor)

                # edgecolor = color
                left = imin + jj*width
                top = jmin + ii*height

                # Actually draws the backgound.
                self.canvas_plot_scene.addRect(left, top, width, height, qcolor, brush)

                # Draw the countours of the various background regions.
                if ii: # top border
                    region_prev = procheck[ii-1][jj]
                    if ((region_prev=='F' and region!='F'   ) or
                        (region_prev=='A' and region in "GD") or
                        (region_prev=='G' and region=='D'   )):
                        self.canvas_plot_scene.addLine(left-.5, top, left+width+.5, top, line_pen)

                if jj: # left border
                    region_prev=procheck[ii][jj-1]
                    if ((region_prev=='F' and region!='F'   ) or
                        (region_prev=='A' and region in "GD") or
                        (region_prev=='G' and region=='D'   )):
                        self.canvas_plot_scene.addLine(left, top-.5, left, top+height+.5, line_pen)


        # Create axis lines.
        def add_text_to_canvas(x, y, text, anchor=None, font_color="black", font_size=6, html_font_size=10):
            _text = str(text)
            text_item = self.canvas_plot_scene.addText(_text)
            w = text_item.boundingRect().width()
            h = text_item.boundingRect().height()
            text_item.setPos(int(x-w/2.5), int(y-h/2.5))
            text_item.setFont(QtGui.QFont(text_item.font().family(), font_size))

        imid = (imin+imax)/2 # midpoint of X-axis
        jmid = (jmin+jmax)/2 # midpoint of Y-axis

        self.canvas_plot_scene.addLine(imin, jmax, imax, jmax)
        self.canvas_plot_scene.addLine(imin, jmin, imin, jmax)
        self.canvas_plot_scene.addLine(imin, jmin, imax, jmin)
        self.canvas_plot_scene.addLine(imax, jmin, imax, jmax)
        self.canvas_plot_scene.addLine(imid, jmin, imid, jmax)
        self.canvas_plot_scene.addLine(imin, jmid, imax, jmid)

        # Create tick marks and labels
        x_x_offset = 0 # 20
        x_y_offset = 14 # 5
        tic = imin
        for label in xticks_num:
            self.canvas_plot_scene.addLine(tic, jmax+5, tic, jmax)
            add_text_to_canvas(tic+x_x_offset, jmax+x_y_offset, label)
            if len(xticks_num)!=1:
                tic+=(imax-imin)/(len(xticks_num)-1)

        tic = jmax
        y_x_offset = 20 # 40
        y_y_offset = 0 # 10
        for label in yticks_num:
            self.canvas_plot_scene.addLine(imin , tic, imin-5, tic)
            add_text_to_canvas(imin-y_x_offset, tic-y_y_offset, text=label)
            if len(yticks_num)!=1:
                tic-=(jmax-jmin)/(len(yticks_num)-1)

        # Phi label.
        add_text_to_canvas((imin+imax)/2+5, jmax+35, text="\u03D5 (degrees)")
        # Psi label.
        add_text_to_canvas(imin/2-10, (jmin+jmax)/2, text="\u03A8")


        #---------------------------------------------------------
        # Actually plots the data for each residue on the scene. -
        #---------------------------------------------------------

        # Parses the PDB structure and gets the phi and psi angles.
        self.residues_count = {"F": 0,       # Favored
                               "A": 0,       # Additional allowed
                               "G": 0,       # Generously allowed
                               "D": 0,       # Disallowed
                               "gly": 0,     # Glycines
                               "pro": 0,     # Prolines
                               "end_res": 0, # end-residues
                               "total": 0,   # total number of residues
                              }

        self.regions_items_dict = {k: [] for k in self.residues_count}
        self.residues_items_dict = {k: [] for k in prot_standard_three_letters}

        for res_idx, res_tags in enumerate(self.residues_tags):

            self.residues_count["total"] += 1

            if res_tags["region"] == "end_res":
                self.residues_count["end_res"] +=1
                continue

            i = imin+(180+res_tags["phi"]/math.pi*180)
            j = jmin+(180-res_tags["psi"]/math.pi*180)


            # Glycines.
            if res_tags["resname"] == "GLY":
                item = Ramachandran_plot_triangle(i=i, j=j, tags=res_tags,
                                                  mark_size=mark_size,
                                                  parent_window=self)
                self.canvas_plot_scene.addItem(item)
                self.regions_items_dict["gly"].append(item)
                self.residues_count["gly"] += 1

            # Prolines.
            elif res_tags["resname"] == "PRO":
                item = Ramachandran_plot_rectangle(i=i, j=j, tags=res_tags,
                                                   mark_size=mark_size,
                                                   parent_window=self)
                self.canvas_plot_scene.addItem(item)
                self.regions_items_dict["pro"].append(item)
                self.residues_count["pro"] += 1

            # Other residues.
            else:
                item = Ramachandran_plot_circle(i=i, j=j, tags=res_tags,
                                                mark_size=mark_size,
                                                parent_window=self)
                self.canvas_plot_scene.addItem(item)
                if res_tags["region"] in self.residues_count:
                    self.regions_items_dict[res_tags["region"]].append(item)
                    self.residues_count[res_tags["region"]] += 1

            if res_tags["resname"] in self.residues_items_dict:
                self.residues_items_dict[res_tags["resname"]].append(item)


        self.residues_count["T"] = sum([self.residues_count[k] for k in ("F", "A", "G", "D")])
        self.regions_items_dict["T"] = [i for k in ("F", "A", "G", "D") for i in self.regions_items_dict[k]]


    def canvas_plot_move_event(self, residue_tags):
        residue_message = "%s %s (phi: %s, psi: %s)" % (residue_tags["resname"], residue_tags["position"],
                                                        residue_tags["phi_str"], residue_tags["psi_str"])
        self.message_label.setText(self.view_label_prefix + residue_message)
        self._activate_residues(residue_tags["resname"])


    def canvas_plot_title_reset(self, residue_tags, group): # mouse leaves title
        self.message_label.setText(self.default_message)
        if group == "resname":
            self._inactivate_residues(residue_tags["resname"])

        elif group == "region":
            self._inactivate_region(residue_tags["region"])

        else:
            raise KeyError("Unknown 'group': %s" % group)


    def info_label_move_event(self, region, message):
        self.message_label.setText(self.view_label_prefix + message)
        self._activate_region(region)

    def info_label_title_reset(self, region):
        self.message_label.setText(self.default_message)
        self._inactivate_region(region)


    def _activate_region(self, region):
        if region in self.regions_items_dict:
            for item in self.regions_items_dict[region]:
                item.highlight_region()

    def _inactivate_region(self, region):
        if region in self.regions_items_dict:
            for item in self.regions_items_dict[region]:
                item.inactivate_region()

    def _activate_residues(self, resname):
        if resname in self.residues_items_dict:
            for item in self.residues_items_dict[resname]:
                item.highlight_region()

    def _inactivate_residues(self, resname):
        if resname in self.residues_items_dict:
            for item in self.residues_items_dict[resname]:
                item.inactivate_region()

    def click_residue_with_left_button(self, residue_tags): # show residue in sticks
        sel = self.pymol_selection + " and resn " + residue_tags["resname"] + " and resi " + residue_tags["position"]
        try:
            cmd.show("sticks", sel)
        except:
            pass

    def release_residue_with_left_button(self, residue_tags): # back to cartoon
        sel = self.pymol_selection + " and resn " + residue_tags["resname"] + " and resi " + residue_tags["position"]
        try:
            cmd.hide("sticks", sel)
            cmd.hide("lines", sel)
            cmd.delete("pymod_selection")
        except:
            pass

    def click_residue_with_middle_button(self, residue_tags): # center & select residue
        sel = self.pymol_selection + " and resn " + residue_tags["resname"] + " and resi " + residue_tags["position"]
        try:
            cmd.center(sel)
            cmd.select("pymod_selection", sel)
            cmd.refresh() # causes the scene to be refresh as soon as it is safe to do so.
        except:
            pass

    def click_residue_with_right_button(self, residue_tags): # center residue
        sel = self.pymol_selection + " and resn " + residue_tags["resname"] + " and resi " + residue_tags["position"]
        try:
            cmd.center(sel, animate=1)
            cmd.delete("pymod_selection")
            cmd.refresh() # causes the scene to be refresh as soon as it is safe to do so.
        except:
            pass


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
