# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Sequences input and output.
"""

import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

from pymod_lib import pymod_vars
from pymod_lib.pymod_seq import seq_manipulation


def build_sequence_file(elements, sequences_filepath, file_format="fasta", remove_indels=True, unique_indices_headers=False, use_structural_information=False, same_length=True, first_element=None):
    """
    Builds a sequence file (the format is specified in the alignment_"format" argument) that will
    contain the sequences supplied in the "elements" which has to contain a list of
    "PyMod_element" class objects.
    """

    output_file_handler = open(sequences_filepath, 'w')

    if same_length:
        seq_manipulation.adjust_aligned_elements_length(elements)

    if first_element != None:
        elements.remove(first_element)
        elements.insert(0, first_element)

    if file_format == "fasta":
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            # Write an output in FASTA format to the output_file_handler given as argument.
            print(">"+header, file=output_file_handler)
            for i in range(0, len(sequence), 60):
                print(sequence[i:i+60], file=output_file_handler)
            print("", file=output_file_handler)

    elif file_format == "pir":
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            sequence += '*'
            structure=''
            if element.has_structure() and use_structural_information:
                structure=element.get_structure_file()
                chain=element.get_chain_id()
            if not structure: # sequence
                print(">P1;"+ header, file=output_file_handler)
                print("sequence:"+header+":::::::0.00:0.00", file=output_file_handler)
            else: # structure
                print(">P1;"+header, file=output_file_handler)
                print("structure:"+structure+":.:"+chain+":.:"+chain+":::-1.00:-1.00", file=output_file_handler)
            for ii in range(0,len(sequence),75):
                print(sequence[ii:ii+75].replace("X","."), file=output_file_handler)

    elif file_format in ("clustal", "stockholm"):
        records = []
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            records.append(SeqRecord(Seq(str(sequence)), id=header))
        SeqIO.write(records, output_file_handler, file_format)

    elif file_format == "pymod":
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            print(header, sequence, file=output_file_handler)

    else:
        output_file_handler.close()
        raise Exception("Unknown file format: %s" % file_format)

    output_file_handler.close()


def get_id_and_sequence_to_print(pymod_element, remove_indels=True, unique_indices_headers=False):
    sequence = pymod_element.my_sequence
    if remove_indels:
        sequence = sequence.replace("-","")
    if not unique_indices_headers:
        header = pymod_element.my_header
    else:
        header = pymod_element.get_unique_index_header()
        # child.my_header.replace(':','_')
    return header, sequence


def convert_sequence_file_format(input_filepath, input_format, output_format, output_filename=None):
    """
    Converts an sequence file specified in the 'input_format' argument in an alignment file
    in the format specified in the 'output_format'.
    """
    input_file_basename = os.path.basename(input_filepath)
    input_file_name = os.path.splitext(input_file_basename)[0]


    if not output_filename:
        output_file_basename = "%s.%s" % (input_file_name, pymod_vars.alignment_extensions_dictionary[output_format])
    else:
        output_file_basename = "%s.%s" % (output_filename, pymod_vars.alignment_extensions_dictionary[output_format])
    output_file_handler = open(os.path.join(os.path.dirname(input_filepath), output_file_basename), "w")


    if input_format == "pymod":
        input_file_handler = open(input_filepath, "r")
        records = [SeqRecord(Seq(l.split(" ")[1].rstrip("\n\r")), id=l.split(" ")[0]) for l in input_file_handler.readlines()]
    else:
        input_file_handler = open(input_filepath, "r")
        records = list(SeqIO.parse(input_file_handler, input_format, alphabet=SingleLetterAlphabet()))


    if output_format == "pymod":
        lines = []
        for i in [(rec.id, rec.seq) for rec in records]:
            lines.append(str(i[0])+'\n')
            lines.append(str(i[1])+'\n')
        output_file_handler.writelines(lines)
    else:
        SeqIO.write(records, output_file_handler, output_format)

    input_file_handler.close()
    output_file_handler.close()
