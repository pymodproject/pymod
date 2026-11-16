# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import re
import shutil

import pymod_lib.pymod_os_specific as pmos
from pymol.Qt import QtCore, QtGui, QtWidgets
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_combobox_qt,
                                                          PyMod_hbox_option_qt)

from ._evolutionary_analysis_base import (Evolutionary_analysis_protocol,
                                          is_a_nucleic_acid_alignment)
from ._web_services_common import Web_services_common


class WebLogo_analysis(Evolutionary_analysis_protocol, Web_services_common):
    """
    Class implementing methods for accessing the WebLogo web service.
    """

    #Units list.
    units_list = ['Bits', 'Probability']
    #Color scheme list.
    colorscheme_list = ['Auto', '(AA) Charge', '(AA) Chemistry', '(AA default) Hydrophobicity', '(NA) Classic', '(NA default) Base pairing']

    # The scheme the window starts on, which is the one each half of the list calls its
    # default. The window used to start on base pairing whatever the alignment was made
    # of, so a protein logo was coloured by a scheme which names four nucleotides: the
    # four residues whose letters are also A, T, G and C were picked out, for a reason
    # which has nothing to do with them, and the rest came out in the default colour.
    default_protein_colorscheme = '(AA default) Hydrophobicity'
    default_nucleic_acid_colorscheme = '(NA default) Base pairing'

    def get_default_colorscheme(self):
        """
        The colour scheme to start on, from what the alignment is made of.
        """
        if is_a_nucleic_acid_alignment(self.input_cluster_element):
            return self.default_nucleic_acid_colorscheme
        return self.default_protein_colorscheme
    # Logo Format. 'PNG image' is the first (and default) choice, because Qt can
    # display a PNG logo inside PyMod, while a PDF one can only be opened with an
    # external viewer (see the 'show_logo_image' method).
    format_list = ['PNG image', 'PDF']

    # Magic numbers and content types of the file formats which the WebLogo server can
    # return. They are used to check that a response really holds the requested image:
    # when a request is not accepted, the server answers with the HTML page of its web
    # form instead of an image file.
    logo_magic_numbers = {'png_print': b'\x89PNG\r\n\x1a\n',
                          'pdf': b'%PDF'}
    logo_content_types = {'png_print': 'image/png',
                          'pdf': 'application/pdf'}

    def additional_initialization(self):
        self.ali_length = len(self.input_cluster_element.my_sequence)

    def launch_from_gui(self):
        self.build_logo_options_window()


    def build_logo_options_window(self):
        """
        Displayes a window with a series of widgets through which users can define WebLogo
        parameters.
        """

        self.logo_window = WebLogo_options_window_qt(self.pymod.main_window,
            protocol=self,
            title="WebLogo 3 web-application Options",
            upper_frame_title="Here you can modify options for WebLogo 3",
            submit_command=self.logo_state)
        self.logo_window.show()


    def check_logo_correct_parameters(self):
        '''
        Checks if the values that were insert in the LOGO window are correct.
        '''
        correct_input = True # This variable defines the status.
        try:
            # Checks if entries are integer numbers.
            start = int(self.logo_window.logo_start.value())
            end = int(self.logo_window.logo_end.value())
            # Get advanced options.
            if self.logo_window.showing_advanced_widgets:
                val = int(self.logo_window.logo_stacks_enf.getvalue())
            # Check on the logic of choosing extremities.
            if start >= end:
                correct_input = False
                errortitle = "Input Error"
                errormessage = "Start value cannot be greater than the end value. Please correct."
                self.pymod.main_window.show_error_message(errortitle, errormessage)
            elif start > self.ali_length or end > self.ali_length or start < 0 or end < 0:
                correct_input = False
                errortitle = "Input Error"
                errormessage = "Values cannot be greater than the sequence length and both must be greater then 0. Please correct."
                self.pymod.main_window.show_error_message(errortitle, errormessage)
        except Exception as e:
            correct_input=False
            errortitle = "Input Error"
            errormessage = "Non valid numeric input. Please correct."
            self.pymod.main_window.show_error_message(errortitle, errormessage)
        return correct_input


    def logo_state(self):
        """
        This method is called when the 'Submit' button on the LOGO window is pressed. It runs a
        check on the entries, if they are correct it calls the getLogo() function
        """
        if not self.check_logo_correct_parameters():
            return False
        self.getLogo()


    def getLogo(self):
        '''
        Generates a LOGO of the alignment, by using WebLogo 3 site.
        Requires active Internet connection.
        '''
        self.verbose = True

        #Units dictionary
        UNITS = {'Bits':'bits', 'Probability':'probability'}
        #Color scheme dictionary
        COLOR_SCHEME = {
            'Auto':'color_auto',
            '(NA default) Base pairing':'color_base_pairing',
            '(NA) Classic':'color_classic',
            '(AA default) Hydrophobicity':'color_hydrophobicity',
            '(AA) Chemistry':'color_chemistry',
            '(AA) Charge':'color_charge'
            }
        #Format dictionary
        FORMATS =  {'PNG image' : 'png_print',    'PDF' : 'pdf'}
        #switch format-extension
        extensions =  {'png_print': 'png',    'pdf' : 'pdf'}
        logo_yesno = {"Yes": "true", "No": "false"}

        #Options defined in the window
        LOGO_UNIT            = UNITS[self.logo_window.unit_combobox.get()]
        LOGO_COLOR           = COLOR_SCHEME[self.logo_window.color_combobox.get()]
        LOGO_RANGE_START     = self.logo_window.logo_start.value()
        LOGO_RANGE_END       = self.logo_window.logo_end.value()
        #Options defined in advanced options sub-window, not always visible. Here they are initialised.
        # Same default as the format combobox of the options window, so that the logo is
        # shown inside PyMod when the advanced options are not used.
        LOGO_FORMAT          = 'png_print'
        LOGO_TITLE           = ''
        LOGO_STACKS_PER_LINE = '80'
        LOGO_SCALE_STACKS    = 'false'
        LOGO_SHOW_ERRORBARS  = 'false'

        if self.logo_window.showing_advanced_widgets:
            LOGO_FORMAT          = FORMATS[self.logo_window.format_combobox.get()]
            LOGO_TITLE           = self.logo_window.logo_title_enf.getvalue()
            LOGO_STACKS_PER_LINE = self.logo_window.logo_stacks_enf.getvalue()
            LOGO_SCALE_STACKS    = logo_yesno[self.logo_window.scale_width_rds.getvalue()]
            LOGO_SHOW_ERRORBARS  = logo_yesno[self.logo_window.show_error_rds.getvalue()]
        self.logo_window.destroy()

        if self.verbose:
            print('- Running GetLogo...')

        # weblogo3 URL.
        # NOTE: the server now redirects http -> https. urllib re-issues a POST as a
        # bodyless GET on that redirect, so the server would answer with the HTML page
        # of its "create" form instead of the logo. Use https directly.
        weblogourl = 'https://weblogo.threeplusone.com/create.cgi'

        #Sets fields and arguments collecting values from the LOGO options window
        values = {'unit_name': LOGO_UNIT, 'color_scheme': LOGO_COLOR,
                  'logo_start': LOGO_RANGE_START, 'logo_end'  : LOGO_RANGE_END,
                  'format': LOGO_FORMAT, 'logo_title': LOGO_TITLE,
                  'stacks_per_line': LOGO_STACKS_PER_LINE,
                  'show_xaxis': 'true', 'show_yaxis': 'true',
                  'show_ends': 'true', 'show_fineprint': 'true', }
        values_update_scale = {'scale_width': LOGO_SCALE_STACKS}
        values_update_errorbars = {'show_errorbars': LOGO_SHOW_ERRORBARS}

        if LOGO_SCALE_STACKS != 'false':
            values.update(values_update_scale)
        if LOGO_SHOW_ERRORBARS != 'false':
            values.update(values_update_errorbars)

        # Builds an url with the multiple alingment and WebLogo parameters and sends a request to
        # the WebLogo server.
        form_upload_file_name = "sequences" # 'sequences_file'
        upload_response, response_content_type = self.upload_alignment(self.input_cluster_element, weblogourl, form_upload_file_name, other_values=values, show_error=False, return_content_type=True)

        #Check if valid response is given
        if upload_response:
            # Checks that the server really returned the requested image before writing
            # it on disk. Without this check, an unexpected response (such as the HTML
            # page of the WebLogo form) would be saved with an image extension and the
            # error would show up only later, when opening the file in an image viewer.
            if not self.check_logo_response(upload_response, response_content_type, LOGO_FORMAT):
                return None
            #Writes output content in a file with extension given by LOGO_FORMAT
            logofile = os.path.join(self.pymod.images_dirpath, 'logo_' + str(self.pymod.logo_image_counter) + '.' + extensions[LOGO_FORMAT])
            lf = open(logofile, 'wb')
            if self.verbose:
                print('- Creating file...')
            lf.write(upload_response)
            lf.close()
            self.pymod.logo_image_counter += 1
            print("- Sequence logo saved at: %s" % logofile)
            self.show_logo_image(logofile, LOGO_FORMAT)
            if self.verbose:
                print('- Done!')
        else:
            if self.verbose:
                print('- No response. Aborted.')


    def show_logo_image(self, logo_filepath, logo_format):
        """
        Displays the sequence logo which was just built. Raster images are shown in a
        window of PyMod, while other formats (PDF) can not be rendered by Qt and are
        opened with the default viewer of the operating system.
        """

        if logo_format == 'png_print':
            pixmap = QtGui.QPixmap(logo_filepath)
            # An unreadable image would give an empty window: in that case fall back
            # on the external viewer, as for the formats which Qt can not render.
            if not pixmap.isNull():
                self.logo_image_window = WebLogo_image_window_qt(self.pymod.main_window,
                    pixmap=pixmap, logo_filepath=logo_filepath,
                    title="Sequence logo: %s" % self.input_cluster_element.my_header)
                self.logo_image_window.show()
                return None
            print("- The sequence logo image could not be loaded in PyMod. Opening it"
                  " with the default viewer of the system.")

        pmos.open_document_with_default_viewer(logo_filepath)


    def check_logo_response(self, response_content, response_content_type, logo_format):
        """
        Checks that the WebLogo server returned the image which was requested. When a
        request is not accepted, the server answers with the HTML page of its web form
        and the response has to be discarded. Returns 'True' when the response holds a
        file of the expected format.
        """

        expected_content_type = self.logo_content_types[logo_format]

        # The magic number is the actual evidence: a declared content type may be
        # missing or wrong, while these formats always begin with a fixed signature.
        if response_content.startswith(self.logo_magic_numbers[logo_format]):
            return True

        # Reports on the PyMOL console what was received instead of the image.
        print("- The WebLogo server did not return a '%s' file." % expected_content_type)
        if response_content_type is not None:
            print("- Content type of the response: '%s'." % response_content_type)
        logo_error_message = self._get_logo_response_error(response_content)
        if logo_error_message is not None:
            print("- The response is the web page: '%s'." % logo_error_message)

        title = "WebLogo Error"
        message = ("The WebLogo server did not return the requested sequence logo."
                   " It answered with a '%s' response, while a '%s' file was expected."
                   "\n\nThe request was probably rejected by the server. Please check"
                   " your Internet connection and the options of the logo, then try"
                   " again." % (response_content_type or "unknown", expected_content_type))
        if logo_error_message is not None:
            message += "\n\nThe server returned the web page: '%s'." % logo_error_message
        self.pymod.main_window.show_error_message(title, message)

        return False


    def _get_logo_response_error(self, response_content):
        """
        Attempts to get the title of an HTML page returned by the WebLogo server, so
        that the page which was received instead of the logo can be identified. Returns
        'None' when the response does not look like an HTML page.
        """

        try:
            response_head = response_content[:4096].decode("utf-8", "replace")
        except Exception:
            return None

        title_match = re.search("<title>(.*?)</title>", response_head,
                                re.IGNORECASE | re.DOTALL)
        if title_match is None:
            return None

        return " ".join(title_match.group(1).split())


class WebLogo_options_window_qt(PyMod_protocol_window_qt):

    def add_middle_frame_widgets(self):

        # Units combobox.
        self.unit_combobox = PyMod_combobox_qt(label_text="Unit Selection",
                                               items=self.protocol.units_list)
        self.unit_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.unit_combobox)

        # Color combobox.
        self.color_combobox = PyMod_combobox_qt(label_text="Color Scheme Selection",
                                                items=self.protocol.colorscheme_list)
        self.color_combobox.combobox.setCurrentText(
            self.protocol.get_default_colorscheme())
        self.middle_formlayout.add_widget_to_align(self.color_combobox)

        # Sub-frame created to display entries for Logo Range option.
        self.range_subframe = PyMod_hbox_option_qt(label_text="Logo Range")

        # Logo start position spinbox.
        self.logo_start = QtWidgets.QSpinBox()
        self.logo_start.setRange(1, self.protocol.ali_length)
        self.range_subframe.hbox.addWidget(self.logo_start)
        # Separator dash.
        self.logo_range_dash = QtWidgets.QLabel(" - ")
        self.range_subframe.hbox.addWidget(self.logo_range_dash)
        # Logo end position spinbox.
        self.logo_end = QtWidgets.QSpinBox()
        self.logo_end.setRange(1, self.protocol.ali_length)
        self.logo_end.setValue(self.protocol.ali_length)
        self.range_subframe.hbox.addWidget(self.logo_end)
        self.range_subframe.set_auto_input_widget_width()
        self.middle_formlayout.add_widget_to_align(self.range_subframe)


        # ADVANCED OPTIONS.
        self.show_advanced_button()

        # Logo format combobox.
        self.format_combobox = PyMod_combobox_qt(label_text='Logo Format',
                                                 items=self.protocol.format_list)
        self.format_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.format_combobox, advanced_option=True)

        # LOGO title entry.
        self.logo_title_enf = PyMod_entryfield_qt(label_text="Logo Title",
                                                  value='')
        self.middle_formlayout.add_widget_to_align(self.logo_title_enf, advanced_option=True)

        # Stacks per line entry.
        self.logo_stacks_enf = PyMod_entryfield_qt(label_text="Stacks per line",
                                                   value='80',
                                                   validate={'validator': 'integer',
                                                             'min': 0, 'max': 100})
        self.middle_formlayout.add_widget_to_align(self.logo_stacks_enf, advanced_option=True)

        # Scale stacks width.
        self.scale_width_rds = PyMod_radioselect_qt(label_text="Scale stacks width",
                                                    buttons=('Yes', 'No'))
        self.scale_width_rds.setvalue("No")
        self.middle_formlayout.add_widget_to_align(self.scale_width_rds, advanced_option=True)

        # Show error bars.
        self.show_error_rds = PyMod_radioselect_qt(label_text="Show error bars",
                                                    buttons=('Yes', 'No'))
        self.show_error_rds.setvalue("No")
        self.middle_formlayout.add_widget_to_align(self.show_error_rds, advanced_option=True)

        self.middle_formlayout.set_input_widgets_width(175)


class WebLogo_image_window_qt(QtWidgets.QMainWindow):
    """
    A window displaying a sequence logo in a QGraphicsView, with zoom controls and
    the possibility of saving a copy of the image file. Built following the same
    scheme of the 'Alignment_figure_window_qt' class of the 'alignment_figure' module.
    """

    is_pymod_window = True

    def __init__(self, parent, pixmap, logo_filepath, title="Sequence logo", *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle(title)
        self.logo_filepath = logo_filepath

        self.scene = QtWidgets.QGraphicsScene()
        self.scene.addPixmap(pixmap)
        self.scene.setSceneRect(QtCore.QRectF(pixmap.rect()))

        self.view = QtWidgets.QGraphicsView(self.scene)
        self.view.setRenderHint(QtGui.QPainter.Antialiasing, True)
        self.view.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
        self.view.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("white")))

        central = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(central)

        button_bar = QtWidgets.QHBoxLayout()
        zoom_in = QtWidgets.QPushButton("Zoom +")
        zoom_in.clicked.connect(lambda: self.view.scale(1.2, 1.2))
        zoom_out = QtWidgets.QPushButton("Zoom -")
        zoom_out.clicked.connect(lambda: self.view.scale(1 / 1.2, 1 / 1.2))
        save_button = QtWidgets.QPushButton("Save image...")
        save_button.clicked.connect(self.save_image)
        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        button_bar.addWidget(zoom_in)
        button_bar.addWidget(zoom_out)
        button_bar.addStretch()
        button_bar.addWidget(save_button)
        button_bar.addWidget(close_button)

        layout.addLayout(button_bar)
        layout.addWidget(self.view)
        self.setCentralWidget(central)

        self.resize(min(1100, pixmap.width() + 60),
                    min(800, pixmap.height() + 90))


    def save_image(self):
        """
        Saves a copy of the logo image built by PyMod. The original file is copied, so
        that the image returned by the WebLogo server is preserved as it is.
        """

        extension = os.path.splitext(self.logo_filepath)[1] or ".png"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save sequence logo",
            os.path.expanduser("~/sequence_logo%s" % extension),
            "Image file (*%s);;All files (*)" % extension)
        if not path:
            return None
        try:
            shutil.copy(self.logo_filepath, path)
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self, "Save Error", "Could not save the sequence logo:\n%s" % str(e))
