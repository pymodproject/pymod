# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os

import pymod_lib.pymod_os_specific as pmos
from pymol.Qt import QtWidgets
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_combobox_qt,
                                                          PyMod_hbox_option_qt)

from ._evolutionary_analysis_base import Evolutionary_analysis_protocol
from ._web_services_common import Web_services_common


class WebLogo_analysis(Evolutionary_analysis_protocol, Web_services_common):
    """
    Class implementing methods for accessing the WebLogo web service.
    """

    #Units list.
    units_list = ['Bits', 'Probability']
    #Color scheme list.
    colorscheme_list = ['Auto', '(AA) Charge', '(AA) Chemistry', '(AA default) Hydrophobicity', '(NA) Classic', '(NA default) Base pairing']
    #Logo Format
    format_list = ['PDF', 'PNG image']

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
        LOGO_FORMAT          = 'pdf'
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

        #weblogo3 URL
        weblogourl = 'http://weblogo.threeplusone.com/create.cgi'

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
        upload_response = self.upload_alignment(self.input_cluster_element, weblogourl, form_upload_file_name, other_values=values, show_error=False)

        #Check if valid response is given
        if upload_response:
            #Writes output content in a file with extension given by LOGO_FORMAT
            logofile = os.path.join(self.pymod.images_dirpath, 'logo_' + str(self.pymod.logo_image_counter) + '.' + extensions[LOGO_FORMAT])
            lf = open(logofile, 'wb')
            if self.verbose:
                print('- Creating file...')
            lf.write(upload_response)
            lf.close()
            self.pymod.logo_image_counter += 1
            pmos.open_document_with_default_viewer(logofile)
            if self.verbose:
                print('- Done!')
        else:
            if self.verbose:
                print('- No response. Aborted.')


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
        self.color_combobox.combobox.setCurrentIndex(5)
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
