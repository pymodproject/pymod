# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Interactions with external tools.
"""

import os
import sys
import subprocess

from pymod_lib import pymod_os_specific as pmos


class PyMod_external:

    modeller_lock_title = "MODELLER warning"
    modeller_lock_message = ("Can not safely exit a MODELLER thread. Please wait"
                             " for it to complete (the only way to quit it is to"
                             " forcibly close PyMOL, but all PyMOL/PyMod data will"
                             " be lost).")


    def execute_subprocess(self, commandline,
                           new_stdout=subprocess.PIPE, new_stderr=subprocess.PIPE,
                           new_shell=(sys.platform!="win32"),
                           verbose=True,
                           executing_modeller=False):

        if verbose:
            print("- Executing the following command:", commandline)

        if not executing_modeller:
            subp = subprocess.Popen(commandline, stdout=new_stdout, stderr=new_stderr, shell=new_shell)
            out_std, err_std = subp.communicate()
            returncode = subp.returncode
            if verbose:
                print("- Stdout:", out_std)
            if returncode != 0:
                if verbose:
                    print("- Code:", returncode, ", Stderr:", err_std)
                raise Exception("Subprocess returned non-zero return code: %s (%s)" % (returncode, err_std))

        # Official PyMOL builds on Mac OS will crash if executing MODELLER through using the
        # 'subprocess' module. For this reason, the 'os' module will be used instead.
        else:
            os.system(commandline)


    def new_execute_subprocess(self, args, verbose=False):
        if verbose:
            print("- Executing the following command:", args)
        subprocess.check_call(args)
