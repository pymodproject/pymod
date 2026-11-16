# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import shutil
import re

import pymod_lib.pymod_vars as pmdt


###################################################################################################
# SALIGN MIXINS.                                                                                  #
###################################################################################################

class SALIGN_alignment:
    """
    Mixin for all SALIGN alignments.
    """
    use_hetatm = False

    def alignment_program_exists(self):
        modeller_error = self.pymod.modeller.check_exception()
        if modeller_error is not None:
            message = "In order to use SALIGN, MODELLER must be installed and configured correctly. %s" % modeller_error
            self.pymod.main_window.show_error_message("MODELLER Error", message)
            return None
        return True

    def alignment_program_not_found(self):
        """
        This method does nothing, the error message is already showed in 'alignment_program_exists'
        in SALIGN protocols.
        """
        pass


class SALIGN_regular_alignment:
    """
    Mixin for SALIGN regular alignments (both sequence and structural).
    """
    def update_additional_information(self):
        """
        Sets the dendrogram file path once the alignment has been performed.
        """
        if len(self.elements_to_align) > 2 and self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
          # Builds a permanent copy of the original temporary .dnd file.
          temp_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, self.protocol_output_file_name + ".tree")
          new_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_%s_dendrogram.tree" % (self.pymod.alignments_files_names, self.alignment_element.unique_index))
          if os.path.isfile(temp_dnd_file_path):
              shutil.copy(temp_dnd_file_path, new_dnd_file_path)
          else:
              return None

          # Edit the new .dnd file to insert the actual names of the sequences.
          dnd_file_handler = open(new_dnd_file_path, "r")
          dnd_file_lines = dnd_file_handler.readlines()
          dnd_file_handler.close()
          new_dnd_file_lines = []
          for line in dnd_file_lines:
              for m in re.findall(pmdt.unique_index_header_regex, line):
                  line = line.replace(m, self.elements_to_align_dict[m].my_header)
              new_dnd_file_lines.append(line)
          dnd_file_handler = open(new_dnd_file_path, "w")
          for line in new_dnd_file_lines:
              dnd_file_handler.write(line)
          dnd_file_handler.close()

          self.alignment_element.set_tree_file_path(new_dnd_file_path)
