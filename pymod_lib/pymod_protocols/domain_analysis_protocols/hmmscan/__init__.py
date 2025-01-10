# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module implementing hmmscan (from the hmmer package) searches in PyMod. This is used to scan
profile-HMM libraries, such as the PFAM library.
"""

import os
import json

import urllib.request, urllib.error, urllib.parse

from Bio import SearchIO

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.domain_analysis_protocols.hmmscan._gui import Hmmscan_options_window_qt, Hmmscan_results_window_qt
from pymod_lib.pymod_element_feature import Domain_feature
from pymod_lib.pymod_os_specific import get_exe_file_name, check_network_connection

from pymod_lib.pymod_threading import Protocol_exec_dialog
from pymod_lib.pymod_exceptions import catch_protocol_exception


class Hmmscan_protocol(PyMod_protocol):

    save_domain_files = False
    protocol_name = "Domain Search"

    def __init__(self, pymod, father_protocol, output_directory=os.path.curdir):
        PyMod_protocol.__init__(self, pymod, pymod.domain_analysis_dirpath)
        self.father_protocol = father_protocol
        self.query_element = self.father_protocol.pymod_element


    def launch_from_gui(self):

        if not self.check_hmmscan_program():
            return None

        title = "HMMSCAN Options"
        if self.father_protocol.domain_search_mode == "remote":
            title = "EBI " + title

        self.hmmer_options_window = Hmmscan_options_window_qt(parent=self.pymod.main_window,
            protocol=self,
            submit_command=self.domain_search_state,
            title=title,
            upper_frame_title="Here you can modify the options for HMMSCAN",)
        self.hmmer_options_window.show()


    def check_hmmscan_program(self):

        if self.father_protocol.domain_search_mode == "local":

            exe_dirpath = self.pymod.hmmer_tool["exe_dir_path"].get_value()
            if not (exe_dirpath and os.path.isdir(exe_dirpath)):
                title = "Hmmer Suite Executable Error"
                message = ("The default Hmmer suite executables directory is missing. Please set one"
                           " in the 'Tools -> Options' menu.")
                self.pymod.main_window.show_error_message(title, message)
                return False

            exe_filename = get_exe_file_name("hmmscan")
            exe_filepath = os.path.join(exe_dirpath, exe_filename)
            if not os.path.isfile(exe_filepath):
                title = "Hmmer Suite Executable Error"
                message = ("A '%s' file is missing in the Hmmer executables directory. Please specify"
                           " in the 'Tools -> Options' menu a Hmmer executables directory where a"
                           " '%s' file is found." % (exe_filename, exe_filename))
                self.pymod.main_window.show_error_message(title, message)
                return False


            db_dirpath = self.pymod.hmmer_tool["hmmscan_db_dir_path"].get_value()
            if db_dirpath.replace(" ", "") == "":
                title = "No Database Directory"
                message = ("No databases directory is defined for HMMSCAN. Please define a database"
                           " directory in the PyMod options window in order to perform a HMMSCAN search.")
                self.pymod.main_window.show_error_message(title, message)
                return False

            if not os.path.isdir(db_dirpath):
                title = "Database Directory not Found"
                message = ("The specified database directory does not exists. Please specify an"
                           " existing one in the Tools -> Options menu.")
                self.pymod.main_window.show_error_message(title, message)
                return False

            self.hmmscan_db_list = [d for d in sorted(os.listdir(db_dirpath)) if d.endswith("hmm.h3m")]
            if len(self.hmmscan_db_list) == 0:
                title = "Hmmscan Database Error"
                message = ("No valid databases files were found in the hmmscan database directory."
                           " Please refer to the PyMod manual to know how to install them.")
                self.pymod.main_window.show_error_message(title, message)
                return False

            return True

        elif self.father_protocol.domain_search_mode == "remote":

            if not check_network_connection("https://google.com", timeout=3):
                title = "Connection Error"
                message = ("An internet connection is not available, can not connect to the EBI"
                           " server to run HMMSCAN.")
                self.pymod.main_window.show_error_message(title, message)
                return False

            return True


    def get_options_from_gui(self):

        # Gets the database on which to perform the search.
        self.hmmscan_db = self.hmmer_options_window.hmmer_database_rds.getvalue()
        if self.hmmscan_db == None:
            raise ValueError("Please select a HMM database perform a domain search.")

        # Gets the evalue cutoff.
        self.evalue_cutoff = self.hmmer_options_window.e_value_threshold_enf.getvalue(validate=True)


    @catch_protocol_exception
    def domain_search_state(self):
        """
        Launched when pressing the 'SUBMIT' button in the hmmscan option window.
        """

        #----------------------------
        # Get options from the GUI. -
        #----------------------------

        try:
            self.get_options_from_gui()
        except Exception as e:
            self.pymod.main_window.show_error_message("Input Error", str(e))
            return None

        self.hmmer_options_window.destroy()

        #----------------------------------------
        # Actually runs the searching protocol. -
        #----------------------------------------

        # Remote.
        if self.father_protocol.domain_search_mode == 'remote':
            self.search_protocol = Hmmscan_web_parsing_protocol(self.pymod, self.father_protocol)

        # Local.
        elif self.father_protocol.domain_search_mode == 'local':
            self.search_protocol = Hmmscan_local_parsing_protocol(self.pymod, self.father_protocol)
            self.hmmscan_db = os.path.join(self.pymod.hmmer_tool["hmmscan_db_dir_path"].get_value(),
                                           self.hmmscan_db_dict[self.hmmscan_db])


        if not self.pymod.use_protocol_threads:
            self.run_search_protocol_scan()
        else:
            if self.father_protocol.domain_search_mode == 'remote':
                title = "Running EBI HMMSCAN"
                label_text = "Connecting to the HMMSCAN EBI server. Please wait for the process to complete..."
            elif self.father_protocol.domain_search_mode == 'local':
                title = "Running HMMSCAN"
                label_text = "Running HMMSCAN. Please wait for the process to complete..."

            p_dialog = Protocol_exec_dialog(app=self.pymod.main_window, pymod=self.pymod,
                                            function=self.run_search_protocol_scan, args=(), wait_start=1,
                                            title=title, label_text=label_text)
            p_dialog.exec_()


        #----------------------------------------------------------
        # Parses the results and show them in the results window. -
        #----------------------------------------------------------

        self.parsed_res = self.run_search_protocol_parse()

        if not self.parsed_res:
            self.pymod.main_window.show_warning_message("Search completed", "No match found with enabled filters.")
            return None

        self.results_window = Hmmscan_results_window_qt(parent=self.pymod.main_window, protocol=self)
        self.results_window.show()


    def run_search_protocol_scan(self):
        self.domain_search_results = self.search_protocol.search_domains(self.query_element,
                                                                         evaluecutoff=self.evalue_cutoff,
                                                                         database=self.hmmscan_db)

    def run_search_protocol_parse(self):
        parser_generator = self.search_protocol.parse(self.domain_search_results)
        array = []
        for i in parser_generator:
            array.append(i.copy())
        sorted_array = sorted(array, key=lambda x: x['evalue'])
        return sorted_array


    #--------------------------------------------------------------------------
    # Import results in PyMod.                                                -
    #--------------------------------------------------------------------------

    def hmmer_results_state(self):
        """
        Called when the 'SUBMIT' button is pressed
        """
        # For each hsp takes the state of its tkinter checkbutton.
        my_domains_map = [[int(v.isChecked()) for v in hit_state] for hit_state in self.results_window.domain_check_states]
        selection_map = [int(v.isChecked()) for v in self.results_window.color_square_lst]

        # If the user selected at least one HSP.
        if 1 in selection_map:
            for hit_ix in range(len(my_domains_map)):
                # hsp_lst = self.results_window.pfam_data[hit_ix]['location']
                hsp_lst = self.parsed_res[hit_ix]['location']
                for loc_ix in range(len(hsp_lst)):
                    # coupling the HSP in self.pfam_data with the status of the checkbutton
                    hsp_lst[loc_ix].update({'selected': my_domains_map[hit_ix][loc_ix]})

            self.import_results_in_pymod()

        self.results_window.close()


    def import_results_in_pymod(self, color_sequence=True, color_structure=True):
        """
        Actually imports the domains in Pymod and assignes information about the selected domain to
        the query PyMod elements and its residues.
        """

        domain_count = 0

        for d_index in range(len(self.parsed_res)):

            d_item = self.parsed_res[d_index]
            d_hsp_lst = self.parsed_res[d_index]['location']

            for d in d_hsp_lst:
                if d['selected']:
                    startindex = int(d['start'])-1 # lo start originale conta da 1
                    endindex = int(d['end'])       # qui va bene, perche' l'ultimo indice e' esclusivo
                    new_domain = Domain_feature(id=d['hsp_number_id'],
                                                # name=d['hsp_number_id'],
                                                name=d['hsp_res_id'],
                                                start=startindex,
                                                end=endindex,
                                                evalue=d['evalue'],
                                                color=d_item['dom_color'],  #tupla
                                                description=d_item['desc'])

                    if self.save_domain_files:
                        o_fh = open(os.path.join(self.output_directory, "domain_%s.txt" % (domain_count+1)), "w")
                        o_fh.write(json.dumps(new_domain.get_feature_dict()))
                        o_fh.close()

                    self.query_element.add_domain_feature(new_domain)

                    domain_count += 1

        # Sorts the domain features list of the PyMod element.
        self.query_element.features["domains"].sort(key=lambda x: x.start)

        if color_sequence:
            self.pymod.main_window.color_selection("single", self.query_element, "domains")

        # Completes the process and stores the results in PyMod.
        self.father_protocol.evaluate_domain_search()


    #-------------------------------------------------------
    # Methods for storing the results of a hmmscan search. -
    #-------------------------------------------------------

    def _initializes_domains_search(self, query_element, **args):
        self.query_element_name = self.father_protocol.query_element_name
        self.query_element_seq = query_element.my_sequence.replace('-', '')
        self.query_filepath = self.father_protocol._pymod_element_seq_filepath


    def get_domain_name(self, parsed_output_item, locitem):
        return parsed_output_item['id']


###################################################################################################
# Remote hmmscan searches.                                                                        #
###################################################################################################

class Hmmscan_web_parsing_protocol(Hmmscan_protocol):

    hmm_url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'

    connection_error_title = "Connection Error"
    connection_error_message = "Can not connect to the EBI hmmscan server. Please check your Internet access."

    def search_domains(self, query, evaluecutoff, database='pfam'):

        self._initializes_domains_search(query)

        self.evaluecutoff = float(evaluecutoff)
        self.database = database

        # install a custom handler to prevent following of redirects automatically.
        class SmartRedirectHandler(urllib.request.HTTPRedirectHandler):
            def http_error_302(self, req, fp, code, msg, headers):
                return headers

        opener = urllib.request.build_opener(SmartRedirectHandler())
        urllib.request.install_opener(opener)

        parameters = {'hmmdb': database.lower(),
                      'seq': self.query_element_seq,
                      }

        enc_params = urllib.parse.urlencode(parameters).encode('utf-8')

        # post the search request to the server
        request = urllib.request.Request(self.hmm_url, enc_params)

        # Get the url where the search results can be fetched from.
        try:
            results_url = urllib.request.urlopen(request).get('location')
        except urllib.error.URLError as e:
            raise e

        # modify the range, format and presence of alignments in your results here
        res_params = {'output': 'json'}

        # add the parameters to your request for the results
        enc_res_params = urllib.parse.urlencode(res_params)
        modified_res_url = results_url + '?' + enc_res_params

        # send a GET request to the server
        results_request = urllib.request.Request(modified_res_url)
        data = urllib.request.urlopen(results_request)


        # Saves the results.
        response_content = data.read().decode('utf-8')
        results_file_name = self.query_element_name + '_web_' + database + '_output.' + res_params['output']
        results_file = os.path.join(self.output_directory, results_file_name)

        with open(results_file, 'w') as o_fh:
            o_fh.write(response_content)

        return results_file


    def parse(self, results_filepath):

        with open(results_filepath, "r") as r_fh:
            json_results = json.loads(r_fh.read())

        # Gets the list of "hits" (that is, the list of domains types identified in the query).
        matchlist = json_results["results"]["hits"]

        parsed_output_item = {}
        for match_idx, match in enumerate(matchlist):

            parsed_output_item.update(_readable_attrs(match))
            unique_id = parsed_output_item['id'] + '*' + str(match_idx).zfill(4)
            parsed_output_item.update({'unique_id': unique_id})
            loclist = match["domains"]

            locationattrs = []
            locitem = {}
            # Each "hit" may have different domains, for example, when a protein has multiple domains
            # of the same type (repeated domains).
            for loc_idx, loc in enumerate(loclist):
                locitem.update(_readable_attrs(loc))
                loc_id = parsed_output_item['id'] + '_hsp_' + str(loc_idx).zfill(3)
                loc_res = self.get_domain_name(parsed_output_item, locitem)
                locitem.update({'hsp_number_id':loc_id, 'hsp_res_id':loc_res})

                try:
                    if self.database.lower() != 'gene3d':
                        if float(locitem['evalue']) < self.evaluecutoff:
                            locationattrs.append(locitem.copy())  #
                    else:
                        locationattrs.append(locitem.copy())
                except AttributeError:
                    locationattrs.append(locitem.copy())

            if locationattrs:
                parsed_output_item.update({'location': locationattrs})  #
                yield parsed_output_item
            else:
                parsed_output_item = {}
                continue


def _readable_attrs(node):
    """
    Changes the keys of the json HMMSCAN ouptut to Biopython ones.
    """
    conversion_standard = {'name': 'id',
                           'aliL': 'length',
                           'ievalue': 'evalue',
                           'ienv': 'env_start',
                           'jenv': 'env_end',
                           'iali': 'ali_start',
                           'jali': 'ali_end',
                           'alisqfrom': 'start',
                           'alisqto': 'end',
                           'alihmmfrom': 'hmm_start',
                           'alihmmto': 'hmm_end',}
    new_node = {}
    for k in node:
        if k in conversion_standard:
            new_node[conversion_standard[k]] = node[k]
        else:
            new_node[k] = node[k]
    return new_node


###################################################################################################
# Local hmmscan searches.                                                                         #
###################################################################################################

class Hmmscan_local_parsing_protocol(Hmmscan_protocol):

    def search_domains(self, query_element, database, evaluecutoff):

        self._initializes_domains_search(query_element)
        self.evaluecutoff = evaluecutoff

        exe_filepath = os.path.join(self.pymod.hmmer_tool["exe_dir_path"].get_value(), get_exe_file_name("hmmscan"))
        out_filepath = os.path.join(self.output_directory, "hmmscan_out_" + self.query_element_name + ".txt")

        cline = [exe_filepath, "-o", out_filepath, "-E", str(evaluecutoff), database.replace('.h3m', ''),
                 self.query_filepath]

        # self.pymod.new_execute_subprocess(cline)
        self.pymod.execute_subprocess(cline, new_shell=False)

        return out_filepath


    def parse(self, file):

        # parsed_output = {}
        parsed_output_item = {}

        inputfile = open(file, 'r')
        for qr in SearchIO.parse(inputfile, 'hmmer3-text'):

            for hit in qr.hits:
                parsed_output_item.update({'id': hit.id, 'evalue': hit.evalue,
                                           'length': qr.seq_len, 'query_descr': qr.description,
                                           'desc': hit.description})
                unique_id = parsed_output_item['id'] + '*' + str(qr.hits.index(hit)).zfill(4)
                parsed_output_item.update({'unique_id': unique_id})
                hhits = hit.hsps

                locattrs = []
                locitem = {}

                for h in hhits:

                    # corresponding_hit = qr[h.hit_id]
                    locitem.update({'id': h.hit_id,
                                    'bitscore': h.bitscore,
                                    'evalue': h.evalue,
                                    'evalue_cond': h.evalue_cond,
                                    'env_start': h.env_start,
                                    'env_end': h.env_end,
                                    'start': int(h.query_start)+1,
                                    'end': h.query_end,
                                    'hmm_start': h.hit_start,
                                    'hmm_end': h.hit_end, })

                    loc_id = parsed_output_item['id'] + '_hsp_' + str(hit.hsps.index(h)).zfill(3)
                    loc_res = self.get_domain_name(parsed_output_item, locitem)
                    locitem.update({'hsp_number_id': loc_id, 'hsp_res_id': loc_res})

                    if locitem['evalue'] < self.evaluecutoff:
                        locattrs.append(locitem.copy())

                parsed_output_item.update({'location': locattrs})

                yield parsed_output_item

        inputfile.close()

        # return parsed_output
