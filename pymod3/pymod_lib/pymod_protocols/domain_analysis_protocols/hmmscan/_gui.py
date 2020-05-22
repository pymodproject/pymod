# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
GUI for performing hmmscan searches in PyMod.
"""

from pymol.Qt import QtWidgets, QtCore, QtGui
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          active_entry_style,
                                                          inactive_entry_style,
                                                          small_font_style,
                                                          highlight_color)
from pymod_lib.pymod_vars import domain_colors_ordered, convert_rgb_to_hex


###################################################################################################
# Hmmscan options window.                                                                         #
###################################################################################################

class Hmmscan_options_window_qt(PyMod_protocol_window_qt):

    def build_protocol_middle_frame(self):

        # Add the buttons to choose the database in which to search for domain profiles.
        if self.protocol.father_protocol.domain_search_mode == 'remote':

            self.hmmer_database_rds = PyMod_radioselect_qt(label_text="Database Selection",
                                                           buttons=("PFAM", "Gene3D"))
            for button in self.hmmer_database_rds.get_buttons():
                button.clicked.connect(self.database_opt_cmd)

        elif self.protocol.father_protocol.domain_search_mode == 'local':

            # Build the list of database names.
            self.protocol.hmmscan_db_dict = {} # This dictionary associates a database code (displayed in the GUI) to its filename.
            db_list = []
            for db_filename in self.protocol.hmmscan_db_list:
                db_name = "".join(db_filename.split(".")[0:-2])
                db_list.append(db_name)
                self.protocol.hmmscan_db_dict[db_name] = db_filename

            # Packs the PHMMER database selection widget.
            self.hmmer_database_rds = PyMod_radioselect_qt(label_text="Database Selection",
                                                           buttons=db_list)

        self.middle_formlayout.add_widget_to_align(self.hmmer_database_rds)


        # E-value selection.
        self.e_value_threshold_enf = PyMod_entryfield_qt(label_text="E-value Threshold",
                                                         value="1.0",
                                                         validate={'validator': 'real',
                                                                   'min': 0.0, 'max': 1000.0})
        self.middle_formlayout.add_widget_to_align(self.e_value_threshold_enf)


        # Note about the Gene3D and Evalues.
        if self.protocol.father_protocol.domain_search_mode == 'remote':
            info_note = ('Note: The Gene3D online database will\n'
                         'ignore custom cut-off parameters since\n'
                         'they use a post processing step that\n'
                         'involves preset thresholds.')
            self.notelabel = QtWidgets.QLabel(info_note)
            self.middle_formlayout.addRow(self.notelabel)

        self.middle_formlayout.set_input_widgets_width(140)


    def database_opt_cmd(self):
        if self.hmmer_database_rds.getvalue() == 'Gene3D':
            self.e_value_threshold_enf.entry.setStyleSheet(inactive_entry_style)
            self.e_value_threshold_enf.entry.setEnabled(False)
        else:
            self.e_value_threshold_enf.entry.setStyleSheet(active_entry_style)
            self.e_value_threshold_enf.entry.setEnabled(True)


###################################################################################################
# Hmmscan results window.                                                                         #
###################################################################################################

# results_header_options = {'background': 'black', 'fg': 'red', 'height': 1, 'padx': 10, 'pady': 10, 'font': "comic 12"}
# results_row_options = {'background': 'black', 'fg': 'white', 'height': 1, 'highlightbackground': 'black', 'font': "comic 11"}

class Hmmscan_results_window_qt(QtWidgets.QMainWindow):
    """
    Window for showing similarity searches results.
    """

    is_pymod_window = True

    def __init__(self, parent, protocol):

        super(Hmmscan_results_window_qt, self).__init__(parent)
        self.protocol = protocol

        self.query_len = len(self.protocol.query_element.my_sequence.replace('-', ''))


        #########################
        # Configure the window. #
        #########################

        self.setWindowTitle("HMMSCAN Results")

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()

        # Parameters used to draw the 'QGraphicsView' widgets for showing domains.
        self.preferred_size_policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred,
                                                           QtWidgets.QSizePolicy.Preferred)
        self.view_bg_color = "transparent"
        self.full_seq_pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 0), 2)
        self.full_seq_color = "#7f7f7f"
        qcolor = QtGui.QColor(0, 0, 0)
        qcolor.setNamedColor(self.full_seq_color)
        self.full_seq_brush = QtGui.QBrush(qcolor)
        self.font_qcolor = QtGui.QColor(220, 220, 220, 255)
        self.font_size = 7


        ################
        # Upper frame. #
        ################

        self.upper_frame = QtWidgets.QFrame()
        self.upper_frame_layout = QtWidgets.QGridLayout()
        self.upper_frame.setLayout(self.upper_frame_layout)
        self.main_vbox.addWidget(self.upper_frame)

        if 'query_descr' in self.protocol.parsed_res[0] and self.protocol.parsed_res[0]['query_descr']:
            labelseq = self.protocol.query_element.my_header # + '\n' + querydescr
        else:
            try:
                if len(self.protocol.query_element.description) > 79:
                    labelseq = self.protocol.query_element.description[:78] + '...'
                else:
                    labelseq = self.protocol.query_element.description
            except TypeError:
                labelseq = self.protocol.query_element.my_header

        self.upper_frame_title = QtWidgets.QLabel("HMMSCAN search results for " + labelseq)
        self.upper_frame_layout.addWidget(self.upper_frame_title)

        #-------------------------
        # Domain graphics frame. -
        #-------------------------

        # Builds the scene where to draw the domain representations.
        self.canvas_plot_scene = QtWidgets.QGraphicsScene()
        self.canvas_plot_view = QtWidgets.QGraphicsView(self.canvas_plot_scene)
        self.canvas_plot_view.setFixedHeight(120)
        self.canvas_plot_view.setSizePolicy(self.preferred_size_policy)
        self.canvas_plot_view.setStyleSheet("background: %s" % self.view_bg_color)
        self.upper_frame_layout.addWidget(self.canvas_plot_view)

        # Draw a rectangle with the full sequence.
        self.x_init = 10
        y_init = 95 # 95
        self.domain_y_init = y_init-7
        self.full_seq_rect_w = 800
        full_seq_rect_h = 10
        self.canvas_plot_scene.addRect(self.x_init, y_init, self.full_seq_rect_w, full_seq_rect_h, self.full_seq_pen, self.full_seq_brush)

        # Draw the labels for the N- and C-terminal residues.
        text_offset_y = 15
        text_offset_x = 10
        text_n = self.canvas_plot_scene.addText("1")
        text_n.setPos(self.x_init-text_offset_x,
                      y_init+text_offset_y)
        text_n.setDefaultTextColor(self.font_qcolor)
        text_n.setFont(QtGui.QFont(text_n.font().family(), self.font_size))

        c_label = str(self.query_len)
        text_c = self.canvas_plot_scene.addText(c_label)
        text_offset_x_add = 5
        if len(c_label) > 2:
            text_offset_x_add = 10
        text_c.setPos(self.x_init+self.full_seq_rect_w-text_offset_x-text_offset_x_add,
                      y_init+text_offset_y)
        text_c.setDefaultTextColor(self.font_qcolor)
        text_c.setFont(QtGui.QFont(text_c.font().family(), self.font_size))


        #################
        # Middle frame. #
        #################

        # Scroll area which contains the widgets, set as the centralWidget.
        self.middle_scroll = QtWidgets.QScrollArea()
        self.main_vbox.addWidget(self.middle_scroll)
        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # Scroll area properties.
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        # QFormLayout in the middle frame.
        self.middle_formlayout = QtWidgets.QFormLayout()
        self.middle_widget.setLayout(self.middle_formlayout)


        #-----------------
        # Results frame. -
        #-----------------

        # Set the frame and its layout.
        self.results_frame = QtWidgets.QFrame()
        self.middle_formlayout.addRow(self.results_frame)
        self.results_grid = QtWidgets.QGridLayout()
        self.results_frame.setLayout(self.results_grid)

        # Calls a method which actually displays the similarity searches results.
        self.display_hmmscan_hits()

        # Align the gridded widgets to the left.
        self.results_grid.setAlignment(QtCore.Qt.AlignLeft)
        self.results_grid.setHorizontalSpacing(30)


        #################
        # Bottom frame. #
        #################

        self.main_button = QtWidgets.QPushButton("Submit")
        self.main_button.clicked.connect(lambda a=None: self.protocol.hmmer_results_state())
        self.main_vbox.addWidget(self.main_button)
        self.main_button.setFixedWidth(self.main_button.sizeHint().width())


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


    def display_hmmscan_hits(self):
        """
        This is used to display in the HMMSCAN results window information for
        each hit and a checkbutton to select it for importing it inside PyMod.
        """

        #------------------------------------
        # Shows the headers of the columns. -
        #------------------------------------

        headers_font_style = "%s; color: %s" % (small_font_style, highlight_color)
        headers_font_style = "%s; font-weight: bold" % (small_font_style)

        self.hmmscan_seq_label = QtWidgets.QLabel("Name")
        self.hmmscan_seq_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.hmmscan_seq_label, 0, 0)

        self.description_label = QtWidgets.QLabel("Description")
        self.description_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.description_label, 0, 1)

        self.hmmscan_e_val_label = QtWidgets.QLabel("E-Value")
        self.hmmscan_e_val_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.hmmscan_e_val_label, 0, 2)

        self.query_span_label = QtWidgets.QLabel("Query span")
        self.query_span_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.query_span_label, 0, 3)

        self.hmm_span_label = QtWidgets.QLabel("HMM span")
        self.hmm_span_label.setStyleSheet(headers_font_style)
        self.results_grid.addWidget(self.hmm_span_label, 0, 4)


        #----------------------
        # Prepare the colors. -
        #----------------------

        self.color_palette_dec = domain_colors_ordered
        # self.color_palette = [[j*255 for j in i[1]] for i in self.color_palette_dec]
        self.color_palette_hex = [convert_rgb_to_hex(i[1]) for i in self.color_palette_dec]

        color_idx = 0
        for hit in self.protocol.parsed_res:
            hit.update({'dom_color': self.color_palette_dec[color_idx],
                        'dom_color_hex':self.color_palette_hex[color_idx]})
            color_idx += 1
            if color_idx == len(self.color_palette_dec):
                color_idx = 0


        #---------------------------------------
        # Show the hits in the results window. -
        #---------------------------------------

        # Keep track of the rows in the grid.
        hmmscan_output_row = 1

        # This is going to contain the list of values of each checkbutton.
        self.color_square_lst = [] # all in a row
        self.domain_check_states = [] # clustered
        self.domain_widgets_dict = {}
        domain_counter = 0

        # Process all the hits. Eeach hit is a different HMM profile, that is, a
        # different domain type.
        for hit in self.protocol.parsed_res:

            # Hit name.
            hit_name_label = QtWidgets.QLabel(hit['id'])
            hit_name_label.setStyleSheet(small_font_style)
            self.results_grid.addWidget(hit_name_label, hmmscan_output_row, 0)

            # Hit description.
            if 'desc' in hit:
                descr_text = hit['desc'][:40] + '...' if len(hit['desc']) > 41 else hit['desc']
            else:
                descr_text = '-'
            hit_descr_label = QtWidgets.QLabel(descr_text)
            hit_descr_label.setStyleSheet(small_font_style)
            self.results_grid.addWidget(hit_descr_label, hmmscan_output_row, 1)

            # E-value info.
            evalue_label = QtWidgets.QLabel(str(hit['evalue']))
            evalue_label.setStyleSheet("%s; color: %s" % (small_font_style, hit['dom_color_hex'])) # results_row_options
            self.results_grid.addWidget(evalue_label, hmmscan_output_row, 2)

            hmmscan_output_row += 1

            self.domain_check_states.append([])

            for domain_hit in hit['location']:

                domain_widgets = {"checkbox": None, "rect": None}

                # Domain rectangle shown in the upper frame of the window.
                domain_rect = self.add_domain_representation(hit, domain_hit)
                domain_widgets["rect"] = domain_rect

                # Checkbox for selection.
                color_square = QtWidgets.QCheckBox(" ")
                self.color_square_lst.append(color_square)
                self.domain_check_states[-1].append(color_square)
                color_square.clicked.connect(lambda a=None, x=domain_counter: self.toggle_domain(x))
                self.results_grid.addWidget(color_square, hmmscan_output_row, 0)
                domain_widgets["checkbox"] = color_square

                # Grahical representation of the domain in the query sequence.
                graphics_view = self.create_little_canvas(hit, domain_hit)
                self.results_grid.addWidget(graphics_view, hmmscan_output_row, 1)

                # Individual E-value info.
                hsp_ievalue_label = QtWidgets.QLabel(str(domain_hit['evalue']))
                hsp_ievalue_label.setStyleSheet(small_font_style)
                self.results_grid.addWidget(hsp_ievalue_label, hmmscan_output_row, 2)

                # Query span info.
                span_info_label = QtWidgets.QLabel(str(domain_hit['start']) + ' - ' + str(domain_hit['end']))
                span_info_label.setStyleSheet(small_font_style)
                self.results_grid.addWidget(span_info_label, hmmscan_output_row, 3)

                # HMM span info.
                hspan_info_text = str(domain_hit['hmm_start']) + ' - ' + str(domain_hit['hmm_end'])
                hspan_info_label = QtWidgets.QLabel(hspan_info_text)
                hspan_info_label.setStyleSheet(small_font_style)
                self.results_grid.addWidget(hspan_info_label, hmmscan_output_row, 4)

                hmmscan_output_row += 1

                self.domain_widgets_dict[domain_counter] = domain_widgets
                domain_counter += 1


    def add_domain_representation(self, hit, domain_hit):

        queryspan_start = int(domain_hit['start'])
        queryspan_end = int(domain_hit['end'])
        domain_x = self.x_init + int(queryspan_start/float(self.query_len)*self.full_seq_rect_w)
        domain_y = self.domain_y_init
        domain_w = int((queryspan_end-queryspan_start)/float(self.query_len)*self.full_seq_rect_w)
        domain_h = 25

        domain_pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255), 1)
        qcolor = QtGui.QColor(0, 0, 0)
        qcolor.setNamedColor(hit['dom_color_hex'])
        domain_brush = QtGui.QBrush(qcolor)

        domain_rect = self.canvas_plot_scene.addRect(domain_x, domain_y, domain_w, domain_h, domain_pen, domain_brush)
        domain_rect.setVisible(False)

        return domain_rect


    def create_little_canvas(self, hit, domain_hit, default_width=300, default_height=11):

        # Builds the graphics view and scene.
        canvas_plot_scene = QtWidgets.QGraphicsScene()
        canvas_plot_view = QtWidgets.QGraphicsView(canvas_plot_scene)
        canvas_plot_view.setFixedHeight(default_height)
        canvas_plot_view.setFixedWidth(default_width)
        canvas_plot_view.setSizePolicy(self.preferred_size_policy)
        canvas_plot_view.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        canvas_plot_view.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        canvas_plot_view.setFrameShape(QtWidgets.QFrame.NoFrame)
        canvas_plot_view.setStyleSheet("border: 0px; background: %s" % self.view_bg_color)

        # Get the coordinates of the various graphics elements.
        one_res_span = default_width/float(self.query_len) # proporzione tra la lunghezza della seq e lo spazio grafico
        queryspan_start = int(domain_hit['start'])
        queryspan_end = int(domain_hit['end'])
        queryspan_start_graphic = int(queryspan_start*one_res_span)
        queryspan_end_graphic = int(queryspan_end*one_res_span)
        # canvas_true_width = int(queryspan_end_graphic-queryspan_start_graphic)
        # space_at_end = int(int(default_width)-(canvas_true_width+queryspan_start_graphic))

        # Draws a gray rectangle representing the full protein sequence.
        canvas_plot_scene.addRect(0, 3, default_width, 5, self.full_seq_pen, self.full_seq_brush)

        # Draws a colored rectangle representing the domain.
        line_pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255), 1)
        qcolor = QtGui.QColor(0, 0, 0)
        qcolor.setNamedColor(hit['dom_color_hex'])
        brush = QtGui.QBrush(qcolor)
        canvas_plot_scene.addRect(queryspan_start_graphic, 0, queryspan_end_graphic-queryspan_start_graphic, default_height, line_pen, brush)

        return canvas_plot_view


    def toggle_domain(self, domain_idx):
        domain_widgets = self.domain_widgets_dict[domain_idx]
        if domain_widgets["rect"].isVisible():
            domain_widgets["rect"].setVisible(False)
        else:
            domain_widgets["rect"].setVisible(True)
