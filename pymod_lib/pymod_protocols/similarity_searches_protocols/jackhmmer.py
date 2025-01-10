# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module for performing jackhmmer searches in PyMod. It mainly builds up on the code used to perform
BLAST searches and phmmer searches.
"""

import os

from pymod_lib.pymod_protocols.similarity_searches_protocols.phmmer import PHMMER_search, Phmmer_options_window_qt
from pymod_lib.pymod_os_specific import get_exe_file_name
from pymod_lib.pymod_gui.shared_gui_components_qt import PyMod_radioselect_qt, PyMod_entryfield_qt


################################################################################################
# JACKHMMER.                                                                                   #
################################################################################################

class Jackhmmer_search(PHMMER_search):

    blast_version = "jackhmmer"
    protocol_name = blast_version
    exe_filename = "jackhmmer"
    # JACKHMMER minimum inclusion E-value.
    min_inclusion_eval_default = 0.001


    def get_blast_window_class_qt(self):
        return Jackhmmer_options_window_qt


    def execute_hmmer_program(self, query_filepath, out_filepath, db_filepath, exe_filepath,
                              report_evalue=10e-6,
                              inclusion_evalue_threshold=10e-6,
                              n_iterations=3):
        """
        Execute the locally installed jackhmmer.
        """

        command_ls = [exe_filepath,
                      "-o", out_filepath, "--domE", str(report_evalue),
                      # "--incE", str(inclusion_evalue_threshold), # Full sequences.
                      "--incdomE", str(inclusion_evalue_threshold), # Domains.
                      "-N", str(n_iterations),
                      query_filepath, db_filepath]
        self.pymod.new_execute_subprocess(command_ls)


    def get_search_record(self):
        phmmer_query_result = parse_jackhmmer_output(self.result_filepath)
        return phmmer_query_result


    def get_hsp_query_seq(self, hsp_obj):
        return hsp_obj.query_seq

    def get_hsp_subj_seq(self, hsp_obj):
        return hsp_obj.hit_seq

    def get_hsp_evalue(self, hsp):
        return hsp.evalue, hsp.cond_evalue

    def get_subj_start(self, hsp):
        return hsp.hit_start

    def build_hsp_list(self):
        hsp_list = []
        query_seq = self.blast_query_element.my_sequence.replace("-", "")

        for hsp in self.search_record:
            hsp._title = hsp.hit_name

            # In jackhmmer the query sequence is changed (it becomes the profile's consensus), so it
            # has to be restored back to its original.
            hsp_query_seq = hsp.query_seq
            hsp_q_res_len = len(hsp_query_seq.replace("-", ""))
            hsp_q_start_index = hsp.query_start
            original_query_frag = query_seq[hsp_q_start_index:hsp_q_start_index+hsp_q_res_len]

            aligned_query_frag = []
            rc = 0
            for p in hsp_query_seq:
                if p != "-":
                    aligned_query_frag.append(original_query_frag[rc])
                    rc += 1
                else:
                    aligned_query_frag.append("-")
            aligned_query_frag = "".join(aligned_query_frag)
            hsp.query_seq = aligned_query_frag

            # Gets additional information on the hsp (such as its the id % and query span).
            hsp.pymod_info = self.get_hsp_info(hsp)
            hsp_list.append(hsp)

        # Sort by evalue.
        hsp_list = list(sorted(hsp_list, key=lambda h: h.evalue))

        return hsp_list


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Jackhmmer_options_window_qt(Phmmer_options_window_qt):
    """
    Window for JACKHMMER searches.
    """

    # geometry_string = "450x650"

    def build_additional_hmmer_widgets(self):

        # Jackhmmer iterations enf.
        self.jackhmmer_iterations_enf = PyMod_entryfield_qt(label_text="JACKHMMER Iterations",
                                                            value='3',
                                                            validate={'validator': 'integer',
                                                                      'min': 1, 'max': 5})
        self.middle_formlayout.add_widget_to_align(self.jackhmmer_iterations_enf, validate=True)


    def build_algorithm_advanced_options_widgets(self):
        # Inclusion E-value selection.
        self.inclusion_e_value_threshold_enf = PyMod_entryfield_qt(
            label_text="Max Inclusion c-Evalue",
            value=str(self.protocol.min_inclusion_eval_default),
            validate={'validator': 'real', 'min': 0.0, 'max': 1000.0})
        self.middle_formlayout.add_widget_to_align(self.inclusion_e_value_threshold_enf,
                                                   advanced_option=True, validate=True)


    def get_additional_hmmer_parameters(self):
        if not self.showing_advanced_widgets:
            inclusion_evalue_threshold = self.protocol.min_inclusion_eval_default
        else:
            inclusion_evalue_threshold = self.inclusion_e_value_threshold_enf.getvalue()
        return {"inclusion_evalue_threshold": inclusion_evalue_threshold,
                "n_iterations": self.jackhmmer_iterations_enf.getvalue()}


###################################################################################################
# Jackhmmer output parser.                                                                        #
###################################################################################################

import re


class Jackhmmer_hsp:

    def __init__(self, hit_name=None, evalue=None, cond_evalue=None):
        self.hit_name = hit_name
        self.cond_evalue = cond_evalue
        self.evalue = evalue
        self.query_seq = ""
        self.hit_seq = ""
        self.hit_start = None
        self.hit_end = None
        self.query_start = None
        self.query_end = None


def _get_int(field_string):
    try:
        return int(field_string)
    except ValueError:
        return None


def parse_jackhmmer_output(output_filepath):
    """
    Parses a jackhmmer output file and returns a list of 'Jackhmmer_hsp' object (one object for each
    domain hit).
    """

    last_round_line_id = None
    with open(output_filepath, "r") as r_fh:
        for i, l in enumerate(r_fh):
            if l.startswith("@@ Continuing to next round."):
                last_round_line_id = i
    if last_round_line_id is None:
        last_round_line_id = 0

    hits_names_list = []
    domains_list = []
    new_alignments_found = False
    hits_found = False
    query_line = None

    with open(output_filepath, "r") as r_fh:

        for i, l in enumerate(r_fh):

            if i < last_round_line_id:
                continue

            # A new sequence hit is found.
            if l.startswith(">>"):
                hits_names_list.append(l[3:].rstrip())
                new_alignments_found = False
                hits_found = True
                domains_i_value_list = []
                domain_count = 0

            elif l.startswith("  == domain"):
                # Add a new domain hit is found for the current sequence hit.
                hit_name = hits_names_list[-1].split(" ")[0]
                domains_list.append(Jackhmmer_hsp(hit_name=hit_name,
                                                  cond_evalue=float(l.split(":")[-1]),
                                                  evalue=domains_i_value_list[domain_count]))
                new_domain = domains_list[-1]
                new_alignments_found = True
                query_line = True
                domain_count += 1

            else:

                # Gets the aligned query and domain hit sequences.
                if new_alignments_found:
                    fields = l.rstrip().split()
                    if len(fields) == 0:
                        continue
                    if fields[-1] in ("RF", "PP"):
                        continue
                    if len(fields) < 3:
                        continue

                    _end = _get_int(fields[-1])
                    _start = _get_int(fields[-3])
                    if _start is not None and _end is not None:

                        seq_label = "Q" if query_line else "H"
                        if query_line:
                            if new_domain.query_start is None:
                                new_domain.query_start = _start - 1
                            new_domain.query_end = _end # + 1
                            new_domain.query_seq += fields[-2].upper().replace(".", "-")
                            query_line = False

                        else:
                            if new_domain.hit_start is None:
                                new_domain.hit_start = _start - 1
                            new_domain.hit_end = _end # + 1
                            new_domain.hit_seq += fields[-2].upper()
                            query_line = True

                else:

                    # Gets the i-evalues from the table summary table of each hit.
                    if hits_found:
                        fields = l.rstrip().split()
                        if len(fields) > 0 and re.match("[0-9]", fields[0]):
                            try:
                                domains_i_value_list.append(float(fields[5]))
                            except ValueError:
                                domains_i_value_list.append(1.0)

    return domains_list
