# Copyright 2020 by Maria Giulia Prado, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import urllib.request, urllib.parse, urllib.error


###################################################################################################
# WEB SERVICES.                                                                                   #
###################################################################################################

class Web_services_common:

    verbose = False

    #################################################################
    # Common methods for interacting with web services.             #
    #################################################################

    def upload_alignment(self, alignment_element, url, form_upload_file_name, structure_element = None, other_values={}, show_error=False):
        '''
        This function creates a POST request to the URL 'url'. The 'form_upload_file_name' argument is the
        name of the form field that encodes the file to be uploaded. For instance: if in the upload form
        the field of the file is called "sequence_file", the 'form_upload_file_name' argument has to be set to
        'sequence_file'. It's equivalent to the 'name' variable of the UNIX command curl:
            curl --form name=@content
        The function saves the current alignment and sends it to the server. It may also send other data,
        encoded in 'other_values' dictionary (a dictionary containing the parameters normally sent by compiling
        a form in the HTML page). This argument is optional and by default is an empty dictionary.
        Returns the response given by the server as a string.
        '''
        response_content = ''

        #Saves alignment in FASTA format
        alignment_file_name = 'alignment_tmp'
        self.pymod.save_alignment_fasta_file(alignment_file_name, alignment_element.get_children(), first_element=structure_element)
        alignment_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.fasta')

        #Copy file content to a string
        al_file = open(alignment_file_path)
        alignment_string = al_file.read()
        al_file.close()
        os.remove(alignment_file_path)

        values = {form_upload_file_name: alignment_string}

        # Adds other values to the url.
        if other_values:
            values.update(other_values)
        # Uploads also a structure file.
        if structure_element != None:
            structure_file = open(os.path.join(self.pymod.structures_dirpath, structure_element.get_structure_file()))
            structure_file_string = structure_file.read()
            dbref_line = "DBREF %s" % (structure_element.my_header).ljust(80, " ")
            structure_file_string = dbref_line + "\n" + structure_file_string
            structure_file.close()
            values.update({"structure_file": structure_file_string})

        user_agent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_6_8)'
        headers = {'User-Agent': user_agent}

        try:
            #Creates a request
            data = urllib.parse.urlencode(values).encode("utf-8")
            # in python 3, params output from urlencode
            # must be encoded to bytes with 'encode' method before it is sent to urlopen as data
            req = urllib.request.Request(url, data, headers=headers)
            #Gets server response and reads it
            response = urllib.request.urlopen(req)
            response_content = response.read()
        except urllib.error.URLError:
            response_content = ''
            title = "Connection Error"
            message = "Can not access the server.\nPlease check your Internet access."
            self.pymod.main_window.show_error_message(title, message)

        return response_content
