# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Protocol to split a sequence in a series of subsequences, each representing a domain.
"""

import os

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_element_feature import Domain_feature

from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_entryfield_qt,
                                                          PyMod_radioselect_qt,
                                                          askyesno_qt)


class Split_into_domains_protocol(PyMod_protocol):

    protocol_name = "Domain Split"

    def __init__(self, pymod, pymod_element, output_directory=None):
        PyMod_protocol.__init__(self, pymod, output_directory)
        self.pymod_element = pymod_element


    def launch_from_gui(self):
        """
        Launches the domain splitting procedure from the PyMod GUI.
        """

        # If the sequence has already some derived domains, ask the users whether to repeat the
        # splitting operation.
        if self.pymod_element.derived_domains_list:
            confirmation = askyesno_qt("Confirm", "Do you want to overwrite the previous splitting operation?", parent=self.pymod.get_qt_parent())
            if not confirmation:
                return None

        # Show the options window.
        self.split_seq_offset_window = SplitSeqOffsetWindow_qt(self.pymod.main_window,
            protocol=self,
            title="Choose the Domain Offset",
            upper_frame_title="Choose how many flanking residues\nto append to each domain",
            submit_command=self.split_seq_offset_window_state)
        self.split_seq_offset_window.show()


    def split_seq_offset_window_state(self):
        """
        Actually split the query sequence in its domains and load the corresponing elements in PyMod.
        """

        # Gets the number of flanking residues from the GUI.
        try:
            n_c_term_offset = int(self.split_seq_offset_window.offset_1_enf.getvalue(validate=True))
        except Exception as e:
            self.pymod.main_window.show_error_message("Input Error", str(e))
            return None

        # Delete previous domains.
        if self.pymod_element.derived_domains_list:
            for el in self.pymod_element.derived_domains_list[:]:
                try:
                    el.delete()
                except ValueError: # If the user has already deleted it.
                    pass
            self.pymod_element.clear_derived_domains()
            self.pymod.main_window.gridder(update_clusters=True, update_menus=True)

        # Build the PyMod elements derived from the split domains.
        self.pymod_element.clear_derived_domains()

        for domain_idx, domain in enumerate(self.pymod_element.get_domains_features()):

            # Build a new PyMod element for each domain.
            ungapped = self.pymod_element.my_sequence.replace('-', '')
            new_startindex = max(0, domain.start-n_c_term_offset)
            new_endindex = min(len(ungapped), domain.end+n_c_term_offset)+1
            new_seq = ungapped[new_startindex:new_endindex]

            new_name = "domain_%s_%s_%s" % (domain_idx+1, self.pymod_element.domain_analysis_id+1, domain.full_name) # TODO.
            my_el = self.pymod.build_pymod_element_from_args(new_name, new_seq)

            # Take the new element representing the split domain and add to it a copy of its
            # corresponding domain feature.
            domain_info = domain.get_feature_dict()
            domain_info.update({"is_derived": True, "offset": (n_c_term_offset, n_c_term_offset)})
            domcopy = Domain_feature(**domain_info)

            my_el.add_domain_feature(domcopy)
            self.pymod_element.derived_domains_list.append(my_el)

            # Loads in PyMod.
            self.pymod.add_element_to_pymod(my_el)


        # Colors the new sequences (the 'gridder' method will be called in the method below).
        self.pymod.main_window.color_selection("multiple", self.pymod_element.derived_domains_list, "domains")

        # Finishes the process.
        self.evaluate_splitting()
        self.split_seq_offset_window.destroy()


    def evaluate_splitting(self):
        self.pymod.deselect_all_sequences()
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)


class SplitSeqOffsetWindow_qt(PyMod_protocol_window_qt):
    """
    Window for the header entry command 'Split Sequence Into Domains'.
    User select the N-term offset anc the C-term offset (in residues)
    to be left before the beginning and after the end of the domain itself.
    """

    default_offset = 10
    # self.geometry("400x250")

    def build_protocol_middle_frame(self):

        # Entryfield for offset selection.
        self.offset_1_enf = PyMod_entryfield_qt(label_text="N-Term and C-Term offset", # label_text = "N-Term offset",
                                                value=str(self.default_offset),
                                                validate={'validator': 'integer',
                                                          'min': 0, 'max': 1000})
        self.middle_formlayout.add_widget_to_align(self.offset_1_enf)

        # self.offset_2_enf = PyMod_entryfield_qt(label_text="C-Term offset", # label_text = "N-Term offset",
        #                                         value=str(self.default_offset),
        #                                         validate={'validator': 'integer',
        #                                                   'min': 0, 'max': 1000})
        # self.middle_formlayout.add_widget_to_align(self.offset_2_enf)

        self.middle_formlayout.set_input_widgets_width(140)
