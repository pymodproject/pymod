# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol

class Evolutionary_analysis_protocol(PyMod_protocol):

    protocol_name = "evolutionary_analysis"

    def __init__(self, pymod, input_cluster_element, *args):
        self.input_cluster_element = input_cluster_element
        PyMod_protocol.__init__(self, pymod, *args)
