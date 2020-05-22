# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Files input and output.
"""

import os

from pymod_lib.pymod_seq import seq_io
from pymod_lib import pymod_vars


class PyMod_files_managment:

    def build_sequence_file(self, elements, sequences_filename, new_directory=None, file_format="fasta", remove_indels=True, unique_indices_headers=False, use_structural_information=False, same_length=True, first_element=None, add_extension=True):
        """
        Wrapper for the 'build_sequence_file' in the 'seq_io' module.
        """

        alignment_extension = pymod_vars.alignment_extensions_dictionary[file_format]

        if new_directory == None:
            target_dirpath = self.alignments_dirpath
        else:
            target_dirpath = new_directory

        if add_extension:
            sequences_filepath = os.path.join(target_dirpath, "%s.%s" % (sequences_filename, alignment_extension))
        else:
            sequences_filepath = os.path.join(target_dirpath, sequences_filename)

        seq_io.build_sequence_file(elements=elements,
                                   sequences_filepath=sequences_filepath,
                                   file_format=file_format,
                                   remove_indels=remove_indels,
                                   unique_indices_headers=unique_indices_headers,
                                   use_structural_information=use_structural_information,
                                   same_length=same_length,
                                   first_element=first_element)


    def save_alignment_fasta_file(self, file_name, aligned_elements, first_element=None, custom_directory=None, unique_indices_headers=False):
        """
        Saves in the Alignments directory a .fasta alignment file containing the sequences of the
        "aligned_elements".
        """
        self.build_sequence_file(aligned_elements, file_name, file_format="fasta",
                                 remove_indels=False, first_element=first_element,
                                 new_directory=custom_directory,
                                 unique_indices_headers=unique_indices_headers)
