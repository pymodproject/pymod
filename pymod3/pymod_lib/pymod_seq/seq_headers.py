# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Operations on the names and ids of the sequences loaded in PyMod.
"""

import re


headers_max_length = 50
max_compact_header_length = 35 # Old PyMod: 11, 10.

def get_header_string(header_str):
    header_str = header_str.split(" ")[0]
    header_str = header_str[0:headers_max_length]
    return header_str


def get_compact_header_string(header_str):
    """
    This is needed to build a compact version of the headers that will be used in various
    instances (for example when the identity matrix table for some alignment is displayed, or
    when assigning the name of ouputs files).
    """
    #----------------------------------
    # Uniprot (sp/tr) and gi entries. -
    #----------------------------------
    # From something like "sp|Q9H492|MLP3A_HUMAN" it will build "sp|Q92622".
    if header_str.startswith("sp|") or header_str.startswith("tr|") or header_str.startswith("gi|"):
        so = re.search(r'(\w{2})\|(\S{6,9})\|([\w\d\_\-|\.]{3,20})', header_str)
        # so = re.search(r'(\w{2})\|(\S{6,9})\|', header_str)
        if so:
            compact_header = so.group(1)+"|"+so.group(2) # +"|"
            # compact_header = compact_header+"|"+so.group(3) # +"|"
            # compact_header = compact_header.rstrip("|")
        else:
            compact_header = crop_header(header_str)
    #----------------------------------------------------------------------
    # Sequences imported from PDB files using the open_pdb_file() method. -
    #----------------------------------------------------------------------
    elif header_str[-8:-1] == "_chain_":
        if len(header_str) == 12: # If it is the name of sequence with a 4 letter PDB id.
            compact_header=header_str[0:4]+"_"+header_str[-1]
        else:
            compact_header=crop_header(header_str[0:-8])+"_"+header_str[-1]
    #-----------------------------------------------
    # Other kind of headers. Just crop the header. -
    #-----------------------------------------------
    else:
        compact_header = crop_header(header_str)

    return compact_header


def crop_header(h):
    return h[0:max_compact_header_length]


def is_uniprotkb_fasta_header(self, header):
    pass


def check_fetchable_pdb_from_header(header):

    # First check for the following format: "1UBI_A" or "1_1UBI_A"
    # match = re.match("([a-zA-Z0-9]{4})_([A-Z])$", header)
    match = re.match("([\d]+_)*([a-zA-Z0-9]{4})_([A-Z])$", header)
    if match:
        groups = match.groups()
        pdb_code = groups[1]
        pdb_chain = groups[2]

    # Then checks for the following format (found in BLAST databases files):
    #     gi|157830059|pdb|1ASH|A (at the time of the old PyMod version)
    #     pdb|1ASH|A (at the time of the new PyMod version)
    else:

        try:
            header_fields = header.split("|")

            if header_fields[0] == "pdb":
                pdb_code = header_fields[1]
                if header_fields[2] != "":
                    pdb_chain = header_fields[2][0]
                else:
                    return None

            elif header_fields[2] == "pdb":
                pdb_code = header_fields[3]
                if header_fields[4] != "":
                    pdb_chain = header_fields[4][0]
                else:
                    return None

            elif header_fields[4] == "pdb":
                pdb_code = header_fields[5]
                if header_fields[6][0] != "":
                    pdb_chain = header_fields[6][0]
                else:
                    return None

            else:
                return None

        except IndexError:
            return None

    return str(pdb_code), str(pdb_chain) # Unicode strings are not recognized by MODELLER.
