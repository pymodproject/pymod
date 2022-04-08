# Copyright 2020 by Dario Marzella, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import math

from pymod_lib import pymod_vars
from ._evolutionary_analysis_base import Evolutionary_analysis_protocol
from pymod_lib.pymod_exceptions import PyModMissingStructure

from pymol.Qt import QtWidgets
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          PyMod_scalebar_qt)

from pymol import cmd


###################################################################################################
# SCR_FIND protocol.                                                                              #
###################################################################################################

class SCR_FIND_analysis(Evolutionary_analysis_protocol):

    def additional_initialization(self):
        self.verbose = False
        self.matrices_initialized = False
        self.hide_non_scrs = True


    def launch_from_gui(self):

        # Checks the elements of the alignment.
        child_elements = [e for e in self.input_cluster_element.get_children() if e.has_structure()]
        if len(child_elements) < 2:
            self.pymod.main_window.show_error_message("Selection Error",
                ("A SCR_FIND analysis can only be performed on alignments in which"
                 " at least two elements have a 3D structure loaded in PyMOL."))
            return None

        if any([e.polymer_type == "nucleic_acid" for e in child_elements]):
            self.pymod.main_window.show_error_message("Selection Error",
                "Can not perform the analysis for nucleic acids structures.")
            return None

        if not self.input_cluster_element.algorithm in pymod_vars.structural_alignment_tools:
            message = pymod_vars.structural_alignment_warning % "SCR_FIND"
            self.pymod.main_window.show_warning_message("Alignment Warning", message)

        # Shows the option window.
        self.build_scr_find_window()


    def build_scr_find_window(self):
        """
        Builds a window with options for the SCR_FIND algorithm.
        """
        self.scr_find_window = SCR_FIND_window_qt(self.pymod.main_window,
            protocol=self,
            title="SCR_FIND algorithm options",
            upper_frame_title="Here you can modify options for SCR_FIND",
            submit_command=self.scr_window_submit)
        self.scr_find_window.show()


    def scr_window_submit(self, value = None):
        if self.verbose:
            print("- Computing centroids. This will take a bit, but just on first SCR_FIND Submit")
        self.hide_non_scrs = pymod_vars.yesno_dict[self.scr_find_window.scr_find_hide_non_scrs.getvalue()]
        try:
            self.GP = self.scr_find_window.scr_find_gap_penalty_enf.getvalue(validate=True)
            self.score_limit = float(self.scr_find_window.sc_scale.get())
            self.window_size = int(self.scr_find_window.sw_scale.get())
        except Exception as e:
            self.pymod.main_window.show_error_message("Input Error", str(e))
            return None
        self.scr_find_state()


    def scr_find_state(self):
        """
        Called when the "SUBMIT" button is pressed on the SCR_FIND window. Contains the code to compute
        SCR_FIND scores using the 'SCR_FIND' class.
        """

        if not self.matrices_initialized:

            if len(set([len(e.my_sequence) for e in self.input_cluster_element.get_children()])) > 1:
                raise Exception("The aligned sequences don't have the same length.")

            self.selected_elements = [e for e in self.input_cluster_element.get_children() if e.has_structure()]

            ###################################
            #Generating a list containing ordered residues coordinates
            ###################################

            self.matrix = []
            self.alignment_length = len(self.selected_elements[0].my_sequence)

            # First get all the coordinates of the alignment elements.
            all_coords_dict = {}

            for pymod_element in self.selected_elements:
                residues, coords = self.get_coords_array(pymod_element=pymod_element,
                                                         interaction_center="ca",
                                                         get_selectors=False)
                all_coords_dict[pymod_element] = dict([(r.db_index, c) for (r, c) in zip(residues, coords)])

            # Then build the matrix used to compute SCR scores.
            for ali_id in range(0, self.alignment_length):

                matrix_row = []

                for pymod_element in self.selected_elements:
                    pos = pymod_element.my_sequence[ali_id]
                    if not pos in ("-", "X"):
                        try:
                            residue = pymod_element.get_residue_by_index(ali_id, aligned_sequence_index=True)
                            coords = list(all_coords_dict[pymod_element][residue.db_index])
                            matrix_row.append((coords, residue))
                        except Exception as e:
                            matrix_row.append((["bugged_residue"], None))
                            message = ("- WARNING: Some problems occurred with"
                                       " structure %s, aminoacid %s, alignment"
                                       " position %s: %s" % (pymod_element.my_header, pos, ali_id, str(e)))
                            print(message)
                    else:
                         matrix_row.append((['-', '-', '-'], None))

                self.matrix.append(matrix_row)


            ###################################
            #Generating a list containing ordered centroids coordinates
            ###################################

            self.centroid_list = []
            for i in range(len(self.matrix)):
                x_list=[]
                y_list=[]
                z_list=[]
                for s in range(len(self.matrix[i])):
                    if "bugged_residue" not in self.matrix[i][s][0]:
                        if self.matrix[i][s][0][0] != '-':
                            x_list.append(self.matrix[i][s][0][0])
                        if self.matrix[i][s][0][1] != '-':
                            y_list.append(self.matrix[i][s][0][1])
                        if self.matrix[i][s][0][2] != '-':
                            z_list.append(self.matrix[i][s][0][2])
                if x_list == []:
                    datax= '-'
                    datay= '-'
                    dataz= '-'
                else:
                    datax= (sum(x_list))/(len(x_list))
                    datay= (sum(y_list))/(len(y_list))
                    dataz= (sum(z_list))/(len(z_list))
                self.centroid_list.append([datax, datay, dataz])

            self.matrices_initialized = True


        ###################################
        #Generating a SC score list
        ###################################

        # i position id, s structure, c coordinate (0 = x, 1 = y, 2 = z), N number of residues in current SCR
        self.score_list = []
        for i in range(len(self.matrix)):
            dc_list = []
            N=0
            gaps=0
            for s in range(len(self.matrix[i])):
                if '-' not in self.matrix[i][s][0] and "bugged_residue" not in self.matrix[i][s][0]:
                    for c in range(0,3):
                        dc = (self.matrix[i][s][0][c] - self.centroid_list[i][c])**2
                        dc_list.append(dc)
                    N+=1
                elif "-" in self.matrix[i][s][0]:
                    gaps+=1
                elif "bugged_residue" in self.matrix[i][s][0]:
                    pass
            if N == 0:
                SC = 1000 + (gaps*(float(self.GP)))
            else:
                SC = ((math.sqrt(sum(dc_list)/(N)))+(gaps*(float(self.GP))))
            for pos_in_str in self.matrix[i]:
                if pos_in_str[1] != None:
                    pos_in_str[1].scr_score = {"score": SC, "interval": None}
            self.score_list.append(SC)


        ################################
        #Finding SCRs with a sliding widow of lenght choosen by the user
        ################################

        #s defines the starting position of the sliding window, e defines its end.
        self.SCR_list=[]
        s = 0
        stn_dev = None
        while s in range((len(self.score_list))-self.window_size):
            e = s + self.window_size
            if e > (len(self.score_list)):
                break
            else:
                mean = (sum(self.score_list[s:e]))/(e-s)
                if mean <= self.score_limit:
                    while mean <= self.score_limit and e <= (len(self.score_list)) and (stn_dev == None or stn_dev <= 4):
                        e+=1
                        mean = (sum(self.score_list[s:e]))/(e-s)
                        devs = []
                        for score in self.score_list[s:e]:
                            dev = (score - mean)**2
                            devs.append(dev)
                        stn_dev = math.sqrt((sum(devs)/((e-s)-1)))
                    start = s+1
                    end = e-1
                    SCR = [start, end]
                    self.SCR_list.append(SCR)
                    s = e
                    stn_dev = None
                else:
                    s+=1

        if self.SCR_list == []:
            print('- WARNING: No SCRs found! try to change your parameters!')
        else:
            for element in self.selected_elements:
                for residue in element.get_polymer_residues():
                    residue.is_scr = False
                    ali_id = residue.get_id_in_aligned_sequence()
                    for SCR in self.SCR_list:
                        if SCR[0] <= ali_id+1 <= SCR[-1]:
                            residue.is_scr = True

        if self.verbose:
            print('- SCRs list:', self.SCR_list)


        #defines 10 color intervals, between the lowest and the highest value for SC score in any SCR
        min_list = []
        max_list = []
        for SCR in self.SCR_list:
            filtered_SCR = [ i for i in (self.score_list[SCR[0]:SCR[-1]])] # if i < self.GP #da 3 in su deviazioni standard
            partial_min = min(filtered_SCR)
            partial_max = max(filtered_SCR)
            min_list.append(partial_min)
            max_list.append(partial_max)
        if min_list == []:
            self.pymod.main_window.show_error_message("No SCR found",
                ("Your SC score limit is below the lowest SC score in the structure."
                " Please increase SC score limit or decrease minimum sliding window length"))
            return None
        else:
            glob_min = min(min_list)
            glob_max = max(max_list)
        click = (glob_max-glob_min)/10
        intervals = []
        for i in range(10):
            intervals.append(((glob_min+(i*click)), (glob_min+((i+1)*click))))

        #Assigns to each residue within ad SCR his color interval (scr_color_interval)

        for element in self.selected_elements:
            show_residue_list = []
            cmd.show("cartoon", element.get_pymol_selector())
            if self.hide_non_scrs:
                cmd.hide("everything", element.get_pymol_selector())
            for residue in element.get_polymer_residues():
                if residue.is_scr and residue.scr_score and residue.scr_score['score'] is not None:
                    color = False
                    i = 1
                    for interval in intervals:
                        if min(interval) <= residue.scr_score["score"] < max(interval):
                            residue.scr_score["interval"] = i
                            color = True
                            break
                        if residue.scr_score["score"] >= max(interval) and i == 10:
                            residue.scr_score["interval"] = 10
                            color = True
                            break
                        i += 1
                    if not color:
                        residue.scr_score["interval"] = 10
                    show_residue_list.append(str(residue.db_index))
                else:
                    residue.scr_score = {"score":None, "interval":None}
            cmd.show("cartoon", "%s and resi %s" % (element.get_pymol_selector(), self.pymod.main_window._join_residues_list(show_residue_list)))

        for element in self.selected_elements:
            self.pymod.main_window.color_element_by_scr_scores(element)


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class SCR_FIND_window_qt(PyMod_protocol_window_qt):

    def add_middle_frame_widgets(self):

        # SC-score scalebar.
        self.sc_scale = PyMod_scalebar_qt(label_text="SC-score Limit",
                                          slider_value=3.0,
                                          slider_from=0.5, slider_to=5.0,
                                          slider_resoution=0.25, slider_tickinterval=1.0,
                                          slider_use_float=True,
                                          slider_binding=self.protocol.scr_window_submit,
                                          slider_width=375)
        self.middle_formlayout.add_widget_to_align(self.sc_scale)

        # Sliding Window size scalebar.
        self.sw_scale = PyMod_scalebar_qt(label_text="Sliding Window Min. Size",
                                          slider_value=3,
                                          slider_from=2, slider_to=50,
                                          slider_resoution=1,  slider_tickinterval=5,
                                          slider_binding=self.protocol.scr_window_submit,
                                          slider_width=375)
        self.middle_formlayout.add_widget_to_align(self.sw_scale)


        # Gap penalty entry field.
        self.scr_find_gap_penalty_enf = PyMod_entryfield_qt(label_text="Gap Penalty",
                                                            value='100',
                                                            validate={'validator': 'integer',
                                                                      'min': 0, 'max': 1000},
                                                            enter_command=self.protocol.scr_window_submit)
        self.middle_formlayout.add_widget_to_align(self.scr_find_gap_penalty_enf)


        # Hide non-SCR residues or show them white.
        self.scr_find_hide_non_scrs = PyMod_radioselect_qt(label_text="Hide non SCRs",
                                                           buttons=('Yes', 'No'))
        self.scr_find_hide_non_scrs.setvalue("No")
        self.scr_find_hide_non_scrs.buttons_dict["No"].clicked.connect(self.protocol.scr_window_submit)
        self.scr_find_hide_non_scrs.buttons_dict["Yes"].clicked.connect(self.protocol.scr_window_submit)
        self.middle_formlayout.add_widget_to_align(self.scr_find_hide_non_scrs)


        self.middle_formlayout.set_input_widgets_width(110)
