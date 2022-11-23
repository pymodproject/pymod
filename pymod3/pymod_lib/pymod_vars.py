# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module contaning variables and function used throughout the PyMod plugin.
"""

###################################################################################################
# Paths.                                                                                          #
###################################################################################################

# Name of the PyMod configuration directory.
pymod_cfg_dirname = ".pymod"

# Name of the Directory containing a subdir with the PyMod conda environment.
pymod_envs_dirname = "pymod_envs"
# Name of the PyMod conda environment.
pymod_env_name = "pymod_env"
# Name of the conda packages directory for PyMod.
pymod_pkgs_dirname = "pymod_pkgs"

# Tools directories names.
blast_databases_dirname = "blast_databases"
hmmer_databases_dirname = "hmmer_databases"
hmmscan_databases_dirname = "hmmscan_databases"

# Log filenames.
data_installer_log_filename = "download_log.txt"


###################################################################################################
# Define some variabiles used throughout the PyMod.                                               #
###################################################################################################

# Tuple containing ids of the algorithms used to perform sequence alignments.
sequence_alignment_tools = ("clustalw", "clustalo", "muscle","salign-seq")
# And structural alignments.
structural_alignment_tools = ("ce", "salign-str")

# Dictionaries for algorithms names.
algs_full_names_dict = {
    # Similarity searches.
    "blast": "NCBI BLAST",
    "blastp": "BLAST",
    "psi-blast": "PSI-BLAST",
    "phmmer": "PHMMER",
    "jackhmmer": "JACKHMMER",
    "hmmsearch": "HMMSEARCH",
    # Alignments.
    "clustalw": "ClustalW",
    "clustalo": "Clustal Omega",
    "muscle": "MUSCLE",
    "salign-seq": "SALIGN",
    "salign-str": "SALIGN",
    "ce": "CE-alignment",
    "imported": "Imported",
    "joined": "Joined"
}

can_show_rmsd_matrix = ("ce","salign-str")
can_show_guide_tree = ("clustalw","clustalo")
can_show_dendrogram = ("salign-seq","salign-str")
can_use_scr_find = ("ce","salign-str")


###################################################################################################
# PyMod elements information.                                                                     #
###################################################################################################

unique_index_header_formatted = "_%s_tpm"
unique_index_header_regex = r"_\d+_tpm"
structure_temp_name = "__%s_structure_full__.pdb" # "__%s_structure_temp__"
structure_chain_temp_name = "__%s_structure_chain_%s__.pdb" # "__%s_structure_temp_chain_%s__"
copied_chain_name = "cobj_%s"


###################################################################################################
# File formats supported in PyMod.                                                                #
###################################################################################################

# The keys are the name of the alignment file format and values are the extension for those
# alignment file.
alignment_extensions_dictionary = {
    "fasta":   "fasta",
    "pir":     "ali",
    "clustal": "aln",
    "stockholm": "sto",
    "pymod":   "txt"}


###################################################################################################
# GUI data.                                                                                       #
###################################################################################################

yesno_dict = {"Yes": True, "No": False}

structural_alignment_warning = ("It looks like that the input alignment you chose was not built through"
                                " a structural alignment algorithm. Please note that a %s analysis"
                                " is only suited for structural alignments where all the 3D structures are"
                                " superposed in PyMOL. It's results are only meaningful in such cases."
                                " If your input structures are already aligned in PyMOL, you can safely"
                                " ignore this warning.")


###################################################################################################
# Tool specific data.                                                                             #
###################################################################################################

# PSIPRED.
psipred_output_extensions = (".ss2", ".horiz")
psipred_element_dict = {"H": "alpha helix", "E": "beta sheet", "C": "aperiodic"}


###################################################################################################
# Substitution Matrices.                                                                             #
###################################################################################################

from Bio.Align import substitution_matrices

list_of_matrices = substitution_matrices.load()

# Create a dict of matrices. The keys are the matrices names #
# HOW TO GET MATRIX: e.g. dict_of_matrices["BLOSUM62"]
dict_of_matrices = {}

for matrix in list_of_matrices:
    items = substitution_matrices.load(matrix).items()
    dict_of_matrices[matrix] = {}
    for item in items:
         dict_of_matrices[matrix][item[0]] = item[1]

# PAM120 manually addedd, missing in Biopython.
dict_of_matrices["PAM120"] = {('W', 'F'): -1, ('L', 'R'): -4, ('S', 'P'): 1, ('V', 'T'): 0, ('Q', 'Q'): 6, ('N', 'A'): -1, ('Z', 'Y'): -5, ('W', 'R'): 1, ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 7, ('S', 'H'): -2, ('H', 'D'): 0, ('L', 'N'): -4, ('W', 'A'): -7, ('Y', 'M'): -4, ('G', 'R'): -4, ('Y', 'I'): -2, ('Y', 'E'): -5, ('B', 'Y'): -3, ('Y', 'A'): -4, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 8, ('G', 'N'): 0, ('E', 'C'): -7, ('Y', 'Q'): -5, ('Z', 'Z'): 4, ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -3, ('T', 'N'): 0, ('P', 'P'): 6, ('V', 'I'): 3, ('V', 'S'): -2, ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -4, ('V', 'Q'): -3, ('K', 'K'): 5, ('P', 'D'): -3, ('I', 'H'): -4, ('I', 'D'): -3, ('T', 'R'): -2, ('P', 'L'): -3, ('K', 'G'): -3, ('M', 'N'): -3, ('P', 'H'): -1, ('F', 'Q'): -6, ('Z', 'G'): -2, ('X', 'L'): -2, ('T', 'M'): -1, ('Z', 'C'): -7, ('X', 'H'): -2, ('D', 'R'): -3, ('B', 'W'): -6, ('X', 'D'): -2, ('Z', 'K'): -1, ('F', 'A'): -4, ('Z', 'W'): -7, ('F', 'E'): -7, ('D', 'N'): 2, ('B', 'K'): 0, ('X', 'X'): -2, ('F', 'I'): 0, ('B', 'G'): 0, ('X', 'T'): -1, ('F', 'M'): -1, ('B', 'C'): -6, ('Z', 'I'): -3, ('Z', 'V'): -3, ('S', 'S'): 3, ('L', 'Q'): -2, ('W', 'E'): -8, ('Q', 'R'): 1, ('N', 'N'): 4, ('W', 'M'): -6, ('Q', 'C'): -7, ('W', 'I'): -6, ('S', 'C'): 0, ('L', 'A'): -3, ('S', 'G'): 1, ('L', 'E'): -4, ('W', 'Q'): -6, ('H', 'G'): -4, ('S', 'K'): -1, ('Q', 'N'): 0, ('N', 'R'): -1, ('H', 'C'): -4, ('Y', 'N'): -2, ('G', 'Q'): -3, ('Y', 'F'): 4, ('C', 'A'): -3, ('V', 'L'): 1, ('G', 'E'): -1, ('G', 'A'): 1, ('K', 'R'): 2, ('E', 'D'): 3, ('Y', 'R'): -5, ('M', 'Q'): -1, ('T', 'I'): 0, ('C', 'D'): -7, ('V', 'F'): -3, ('T', 'A'): 1, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -2, ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -2, ('K', 'H'): -2, ('V', 'R'): -3, ('P', 'C'): -4, ('M', 'E'): -3, ('K', 'L'): -4, ('V', 'V'): 5, ('M', 'I'): 1, ('T', 'Q'): -2, ('I', 'G'): -4, ('P', 'K'): -2, ('M', 'M'): 8, ('K', 'D'): -1, ('I', 'C'): -3, ('Z', 'D'): 3, ('F', 'R'): -5, ('X', 'K'): -2, ('Q', 'D'): 1, ('X', 'G'): -2, ('Z', 'L'): -3, ('X', 'C'): -4, ('Z', 'H'): 1, ('B', 'L'): -4, ('B', 'H'): 1, ('F', 'F'): 8, ('X', 'W'): -5, ('B', 'D'): 4, ('D', 'A'): 0, ('S', 'L'): -4, ('X', 'S'): -1, ('F', 'N'): -4, ('S', 'R'): -1, ('W', 'D'): -8, ('V', 'Y'): -3, ('W', 'L'): -3, ('H', 'R'): 1, ('W', 'H'): -3, ('H', 'N'): 2, ('W', 'T'): -6, ('T', 'T'): 4, ('S', 'F'): -3, ('W', 'P'): -7, ('L', 'D'): -5, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1, ('B', 'T'): 0, ('L', 'L'): 5, ('Y', 'K'): -5, ('E', 'Q'): 2, ('Y', 'G'): -6, ('Z', 'S'): -1, ('Y', 'C'): -1, ('G', 'D'): 0, ('B', 'V'): -3, ('E', 'A'): 0, ('Y', 'W'): -2, ('E', 'E'): 5, ('Y', 'S'): -3, ('C', 'N'): -5, ('V', 'C'): -3, ('T', 'H'): -3, ('P', 'R'): -1, ('V', 'G'): -2, ('T', 'L'): -3, ('V', 'K'): -4, ('K', 'Q'): 0, ('R', 'A'): -3, ('I', 'R'): -2, ('T', 'D'): -1, ('P', 'F'): -5, ('I', 'N'): -2, ('K', 'I'): -3, ('M', 'D'): -4, ('V', 'W'): -8, ('W', 'W'): 12, ('M', 'H'): -4, ('P', 'N'): -2, ('K', 'A'): -2, ('M', 'L'): 3, ('K', 'E'): -1, ('Z', 'E'): 4, ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -2, ('X', 'F'): -3, ('K', 'C'): -7, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -4, ('F', 'C'): -6, ('Z', 'Q'): 4, ('X', 'Z'): -1, ('F', 'G'): -5, ('B', 'E'): 3, ('X', 'V'): -1, ('F', 'K'): -7, ('B', 'A'): 0, ('X', 'R'): -2, ('D', 'D'): 5, ('W', 'G'): -8, ('Z', 'F'): -6, ('S', 'Q'): -2, ('W', 'C'): -8, ('W', 'K'): -5, ('H', 'Q'): 3, ('L', 'C'): -7, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -5, ('W', 'S'): -2, ('S', 'E'): -1, ('H', 'E'): -1, ('S', 'I'): -2, ('H', 'A'): -3, ('S', 'M'): -2, ('Y', 'L'): -2, ('Y', 'H'): -1, ('Y', 'D'): -5, ('E', 'R'): -3, ('X', 'P'): -2, ('G', 'G'): 5, ('G', 'C'): -4, ('E', 'N'): 1, ('Y', 'T'): -3, ('Y', 'P'): -6, ('T', 'K'): -1, ('A', 'A'): 3, ('P', 'Q'): 0, ('T', 'C'): -3, ('V', 'H'): -3, ('T', 'G'): -1, ('I', 'Q'): -3, ('Z', 'T'): -2, ('C', 'R'): -4, ('V', 'P'): -2, ('P', 'E'): -2, ('M', 'C'): -6, ('K', 'N'): 1, ('I', 'I'): 6, ('P', 'A'): 1, ('M', 'G'): -4, ('T', 'S'): 2, ('I', 'E'): -3, ('P', 'M'): -3, ('M', 'K'): 0, ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 6, ('X', 'M'): -2, ('L', 'I'): 1, ('X', 'I'): -1, ('Z', 'B'): 2, ('X', 'E'): -1, ('Z', 'N'): 0, ('X', 'A'): -1, ('B', 'R'): -2, ('B', 'N'): 3, ('F', 'D'): -7, ('X', 'Y'): -3, ('Z', 'R'): -1, ('F', 'H'): -3, ('B', 'F'): -5, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4}


###################################################################################################
# Color dictionaries.                                                                             #
###################################################################################################

def convert_rgb_to_hex(rgb_tuple):
    """
    Converts an RGB color code (scaled from 0 to 1) to its corresponding HEX code.
    """
    rgb_tuple = [int(i*255) for i in rgb_tuple]
    hexa = [format(c, "02x") for c in rgb_tuple]
    return "#" + "".join(hexa)

def convert_hex_to_rgb(hex_string):
    hex_string = hex_string.lstrip("#")
    return tuple(int(hex_string[i:i+2], 16)/255.0 for i in (0, 2, 4))


#--------------------------
# Regular colors palette. -
#--------------------------

# PyMod color palette.
pymod_regular_colors_list = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',
                             'salmon', 'pink', 'magenta', 'orange', 'purple',
                             'firebrick', 'chocolate', 'gray', 'white']
pymod_regular_colors_dict = {'red': "#ff0000",
                             'green': "#00ff00",
                             'blue': "#0000ff",
                             'yellow': "#ffff00",
                             'violet': "#ee82ee",
                             'cyan': "#00ffff",
                             'salmon': "#fa8072",
                             'pink': "#ffc0cb",
                             'magenta': "#ff00ff",
                             'orange': "#ffa500",
                             'purple': "#a020f0",
                             'firebrick': "#b22222",
                             'chocolate': "#d2691e",
                             'gray': "#bebebe",
                             'white': "#ffffff",}

# PyMOL color cycle.
# pymol_full_color_cycle = [
#     "carbon",
#     "cyan",
#     "lightmagenta",
#     "yellow",
#     "salmon",
#     "hydrogen",
#     "slate",
#     "orange",
#     "lime",
#     "deepteal",
#     "hotpink",
#     "yelloworange",
#     "violetpurple",
#     "grey70",
#     "marine",
#     "olive",
#     "smudge",
#     "teal",
#     "dirtyviolet",
#     "wheat",
#     "deepsalmon",
#     "lightpink",
#     "aquamarine",
#     "paleyellow",
#     "limegreen",
#     "skyblue",
#     "warmpink",
#     "limon",
#     "violet",
#     "bluewhite",
#     "greencyan",
#     "sand",
#     "forest",
#     "lightteal",
#     "darksalmon",
#     "splitpea",
#     "raspberry",
#     "grey50",
#     "deepblue",
#     "brown"]

# PyMOL colors palette. Used to color structures loaded in PyMOL. Takes only a certain number of
# colors from the full cycle.
pymol_regular_colors_list = [
    "carbon",
    "cyan",
    "lightmagenta",
    "yellow",
    "salmon",
    "hydrogen",
    "slate",
    "orange",
    "lime",
    "deepteal",
    "hotpink",
    "yelloworange",
    "violetpurple",
    "grey70",
    "marine",
    "olive",
    "smudge",
    "teal",
    "dirtyviolet",
    "wheat"]

# Obtained from: https://pymolwiki.org/index.php/Color_Values
pymol_regular_colors_dict_rgb = {
    "carbon": (0.2, 1.0, 0.2),
    "cyan": (0.0, 1.0, 1.0),
    "lightmagenta": (1.0, 0.2, 0.8),
    "yellow": (1.0, 1.0, 0.0),
    "salmon": (1.0, 0.6, 0.6),
    "hydrogen": (0.9, 0.9, 0.9),
    "slate": (0.5, 0.5, 1.0),
    "orange": (1.0, 0.5, 0.0),
    "lime": (0.5, 1.0, 0.5),
    "deepteal": (0.1, 0.6, 0.6),
    "hotpink": (1.0, 0.0, 0.5),
    "yelloworange": (1.0, 0.87, 0.37),
    "violetpurple": (0.55, 0.25, 0.60),
    "grey70": (0.7, 0.7, 0.7),
    "marine": (0.0, 0.5, 1.0),
    "olive": (0.77, 0.70, 0.00),
    "smudge": (0.55, 0.70, 0.40),
    "teal": (0.00, 0.75, 0.75),
    "dirtyviolet": (0.70, 0.50, 0.50),
    "wheat": (0.99, 0.82, 0.65)}

# PyMOL light colors palette. Used to color multiple chains models.
pymol_light_colors_prefix = "l"
pymol_light_colors_list = [pymol_light_colors_prefix +"_"+c for c in pymol_regular_colors_list]
pymol_light_colors_dict_rgb = {
    pymol_light_colors_prefix + "_carbon": (0.80, 1.0, 0.80), # (0.94, 1.0, 0.94),
    pymol_light_colors_prefix + "_cyan": (0.80, 1.0, 1.0),
    pymol_light_colors_prefix + "_lightmagenta": (1.0, 0.75, 0.94),
    pymol_light_colors_prefix + "_yellow": (1.0, 1.0, 0.8),
    pymol_light_colors_prefix + "_salmon": (1.0, 0.86, 0.86),
    pymol_light_colors_prefix + "_hydrogen": (0.98, 0.98, 0.98),
    pymol_light_colors_prefix + "_slate": (0.85, 0.85, 1.0),
    pymol_light_colors_prefix + "_orange": (1.0, 0.9, 0.7),
    pymol_light_colors_prefix + "_lime": (0.85, 1.0, 0.85),
    pymol_light_colors_prefix + "_deepteal": (0.7, 0.95, 0.95),
    pymol_light_colors_prefix + "_hotpink": (1.0, 0.8, 0.98),
    pymol_light_colors_prefix + "_yelloworange": (1.0, 0.94, 0.70),
    pymol_light_colors_prefix + "_violetpurple": (0.58, 0.5, 0.60),
    pymol_light_colors_prefix + "_grey70": (0.87, 0.87, 0.87),
    pymol_light_colors_prefix + "_marine": (0.7, 0.98, 1.0),
    pymol_light_colors_prefix + "_olive": (0.77, 0.75, 0.63),
    pymol_light_colors_prefix + "_smudge": (0.65, 0.70, 0.62),
    pymol_light_colors_prefix + "_teal": (0.56, 0.75, 0.75),
    pymol_light_colors_prefix + "_dirtyviolet": (0.70, 0.62, 0.62),
    pymol_light_colors_prefix + "_wheat": (0.99, 0.92, 0.87)}

#-----------------------------
# Single amino acids colors. -
#-----------------------------

# Starts to define color dictionaries for amminoacids.
# residue_color_dict = {
#     "A": "blue",
#     "L": "blue",
#     "I": "blue",
#     "M": "blue",
#     "W": "blue",
#     "F": "blue",
#     "V": "blue",
#     "T": "green",
#     "N": "green",
#     "Q": "green",
#     "S": "green",
#     "P": "yellow",
#     "G": "orange",
#     "R": "red",
#     "K": "red",
#     "C": "pink",
#     "D": "magenta",
#     "E": "magenta",
#     "H": "cyan",
#     "Y": "cyan",
#     "X": "white",
#     "-": "white"}

#--------------------------------------------------------------------------
# Used to color residues according to their secondary structure. -
#--------------------------------------------------------------------------
# Observed secondary structure.
pymol_obs_sec_str_name = "pymod_oss"
sec_str_color_dict = {
          pymol_obs_sec_str_name + "_H": (1.0, 0.0, 0.0),    # PyMOL helix: red.
          pymol_obs_sec_str_name + "_S": (1.0, 1.0, 0.0), # PyMOL sheet: yellow.
          pymol_obs_sec_str_name + "_L": (0.0, 1.0, 0.0),  # PyMOL aperiodic: green.
          pymol_obs_sec_str_name + "_" : (1.0, 1.0, 1.0),
          pymol_obs_sec_str_name + "_None" : (1.0, 1.0, 1.0)}

# Predicted secondary structure.
pymol_psipred_color_name = "pymod_psipred"
# Generates names like: 'pymod_psipred_8_H' (the name of the color with which residues predicted in
# an helix with confidence score of 8 will be colored).
psipred_color_dict = {
      # Helices.
      pymol_psipred_color_name + "_9_H": (1.0, 0.0, 0.0),
      pymol_psipred_color_name + "_8_H": (1.0, 0.098039215686274508, 0.098039215686274508),
      pymol_psipred_color_name + "_7_H": (1.0, 0.20000000000000001, 0.20000000000000001),
      pymol_psipred_color_name + "_6_H": (1.0, 0.29803921568627451, 0.29803921568627451),
      pymol_psipred_color_name + "_5_H": (1.0, 0.40000000000000002, 0.40000000000000002),
      pymol_psipred_color_name + "_4_H": (1.0, 0.50196078431372548, 0.50196078431372548),
      pymol_psipred_color_name + "_3_H": (1.0, 0.59999999999999998, 0.59999999999999998),
      pymol_psipred_color_name + "_2_H": (1.0, 0.70196078431372544, 0.70196078431372544),
      pymol_psipred_color_name + "_1_H": (1.0, 0.80000000000000004, 0.80000000000000004),
      pymol_psipred_color_name + "_0_H": (1.0, 0.90196078431372551, 0.90196078431372551),
      # Sheets.
      pymol_psipred_color_name + "_9_E": (1.0, 1.0, 0.0),
      pymol_psipred_color_name + "_8_E": (1.0, 1.0, 0.098039215686274508),
      pymol_psipred_color_name + "_7_E": (1.0, 1.0, 0.20000000000000001),
      pymol_psipred_color_name + "_6_E": (1.0, 1.0, 0.29803921568627451),
      pymol_psipred_color_name + "_5_E": (1.0, 1.0, 0.40000000000000002),
      pymol_psipred_color_name + "_4_E": (1.0, 1.0, 0.50196078431372548),
      pymol_psipred_color_name + "_3_E": (1.0, 1.0, 0.59999999999999998),
      pymol_psipred_color_name + "_2_E": (1.0, 1.0, 0.70196078431372544),
      pymol_psipred_color_name + "_1_E": (1.0, 1.0, 0.80000000000000004),
      pymol_psipred_color_name + "_0_E": (1.0, 1.0, 0.90196078431372551),
      # Aperiodic.
      pymol_psipred_color_name + "_9_C": (0.0, 1.0, 0.0),
      pymol_psipred_color_name + "_8_C": (0.098039215686274508, 1.0, 0.098039215686274508),
      pymol_psipred_color_name + "_7_C": (0.20000000000000001, 1.0, 0.20000000000000001),
      pymol_psipred_color_name + "_6_C": (0.29803921568627451, 1.0, 0.29803921568627451),
      pymol_psipred_color_name + "_5_C": (0.40000000000000002, 1.0, 0.40000000000000002),
      pymol_psipred_color_name + "_4_C": (0.50196078431372548, 1.0, 0.50196078431372548),
      pymol_psipred_color_name + "_3_C": (0.59999999999999998, 1.0, 0.59999999999999998),
      pymol_psipred_color_name + "_2_C": (0.70196078431372544, 1.0, 0.70196078431372544),
      pymol_psipred_color_name + "_1_C": (0.80000000000000004, 1.0, 0.80000000000000004),
      pymol_psipred_color_name + "_0_C": (0.90196078431372551, 1.0, 0.90196078431372551)
      }

#---------------------------------------------------
# A dictionary containing colors for CAMPO scores. -
#---------------------------------------------------

pymol_campo_color_name = "pymod_campo"
campo_color_dict = {
    pymol_campo_color_name + "_None": (1,1,1), # (1,1,1),
    pymol_campo_color_name + "_1": (0.0, 0.1588235294117647, 1.0), # (0.0, 0.0, 0.5),
    pymol_campo_color_name + "_2": (0.0, 0.50392156862745097, 1.0), # (0.0, 0.0, 0.94563279857397498),
    pymol_campo_color_name + "_3": (0.0, 0.83333333333333337, 1.0), # (0.0, 0.29999999999999999, 1.0),
    pymol_campo_color_name + "_4": (0.21189120809614148, 1.0, 0.75585072738772952), # (0.0, 0.69215686274509802, 1.0),
    pymol_campo_color_name + "_5": (0.49019607843137247, 1.0, 0.47754585705249841), # (0.16129032258064513, 1.0, 0.80645161290322587),
    pymol_campo_color_name + "_6": (0.75585072738772918, 1.0, 0.2118912080961417), # (0.49019607843137247, 1.0, 0.47754585705249841),
    pymol_campo_color_name + "_7": (1.0, 0.9012345679012348, 0.0), # (0.80645161290322565, 1.0, 0.16129032258064513),
    pymol_campo_color_name + "_8": (1.0, 0.58169934640522891, 0.0), # (1.0, 0.7705156136528688, 0.0),
    pymol_campo_color_name + "_9": (1.0, 0.27668845315904156, 0.0), # (1.0, 0.40740740740740755, 0.0),
    pymol_campo_color_name + "_10": (0.8743315508021392, 0.0, 0.0) # (0.94563279857397531, 0.029774872912127992, 0.0)
}

#--------------------------
# SCR_FIND values colors. -
#--------------------------

pymol_scr_color_name = "pymod_scr"
scr_color_dict = {
    pymol_scr_color_name + "_None": (1,1,1),
    pymol_scr_color_name + "_1": (0.8743315508021392, 0.0, 0.0),
    pymol_scr_color_name + "_2": (1.0, 0.27668845315904156, 0.0),
    pymol_scr_color_name + "_3": (1.0, 0.58169934640522891, 0.0),
    pymol_scr_color_name + "_4": (1.0, 0.9012345679012348, 0.0),
    pymol_scr_color_name + "_5": (0.75585072738772918, 1.0, 0.2118912080961417),
    pymol_scr_color_name + "_6": (0.49019607843137247, 1.0, 0.47754585705249841),
    pymol_scr_color_name + "_7": (0.21189120809614148, 1.0, 0.75585072738772952),
    pymol_scr_color_name + "_8": (0.0, 0.83333333333333337, 1.0),
    pymol_scr_color_name + "_9": (0.0, 0.50392156862745097, 1.0),
    pymol_scr_color_name + "_10": (0.0, 0.1588235294117647, 1.0)}


#-----------------------------------------------------
# A dictionary containing colors for entropy scores. -
#-----------------------------------------------------

pymol_entropy_color_name = "pymod_entropy"
entropy_color_dict = {
    pymol_entropy_color_name + "_None": (1.0, 1.0, 1.0),
    pymol_entropy_color_name + "_9": (0.5490196078431373, 0.0392156862745098, 0.6431372549019608),
    pymol_entropy_color_name + "_8": (0.6627450980392157, 0.13725490196078433, 0.5843137254901961),
    pymol_entropy_color_name + "_7": (0.7607843137254902, 0.23921568627450981, 0.5019607843137255),
    pymol_entropy_color_name + "_6": (0.8431372549019608, 0.3411764705882353, 0.4196078431372549),
    pymol_entropy_color_name + "_5": (0.9137254901960784, 0.4470588235294118, 0.3411764705882353),
    pymol_entropy_color_name + "_4": (0.9686274509803922, 0.5215686274509804, 0.054901960784313725),
    pymol_entropy_color_name + "_3": (0.9921568627450981, 0.6862745098039216, 0.19215686274509805),
    pymol_entropy_color_name + "_2": (0.984313725490196, 0.8274509803921568, 0.1411764705882353),
    pymol_entropy_color_name + "_1": (0.9372549019607843, 0.9725490196078431, 0.12941176470588237)
}

# # Consurf color scheme.
# def get_rgb_tuple(t):
#     return tuple([i/255.0 for i in t])
#
# entropy_color_dict = {
#     pymol_entropy_color_name + "_None": (0.7, 0.7, 0.7),
#     pymol_entropy_color_name + "_9": get_rgb_tuple((15, 199, 207)),
#     pymol_entropy_color_name + "_8": get_rgb_tuple((143, 255, 255)),
#     pymol_entropy_color_name + "_7": get_rgb_tuple((208, 255, 255)),
#     pymol_entropy_color_name + "_6": get_rgb_tuple((224, 255, 255)),
#     pymol_entropy_color_name + "_5": get_rgb_tuple((255, 255, 255)),
#     pymol_entropy_color_name + "_4": get_rgb_tuple((255, 232, 240)),
#     pymol_entropy_color_name + "_3": get_rgb_tuple((240, 199, 223)),
#     pymol_entropy_color_name + "_2": get_rgb_tuple((235, 118, 157)),
#     pymol_entropy_color_name + "_1": get_rgb_tuple((159, 32, 95)),
# }

#----------------------------
# Pair conservation colors. -
#----------------------------

pymol_pc_color_name = "pymod_pc"
pc_color_dict = {
    pymol_pc_color_name + "_0": (1.0, 1.0, 1.0),
    pymol_pc_color_name + "_1": (1.0, 0.75, 0.75),
    pymol_pc_color_name + "_2": (1.0, 0.30, 0.30), # (1.0, 0.25, 0.25),
}

#----------------------
# DOPE values colors. -
#----------------------

pymol_dope_color_name = "pymod_dope"
dope_color_dict = {
    pymol_dope_color_name + "_None": (1,1,1),
    pymol_dope_color_name + "_1": (0.0, 0.1588235294117647, 1.0),
    pymol_dope_color_name + "_2": (0.0, 0.50392156862745097, 1.0),
    pymol_dope_color_name + "_3": (0.0, 0.83333333333333337, 1.0),
    pymol_dope_color_name + "_4": (0.21189120809614148, 1.0, 0.75585072738772952),
    pymol_dope_color_name + "_5": (0.49019607843137247, 1.0, 0.47754585705249841),
    pymol_dope_color_name + "_6": (0.75585072738772918, 1.0, 0.2118912080961417),
    pymol_dope_color_name + "_7": (1.0, 0.9012345679012348, 0.0),
    pymol_dope_color_name + "_8": (1.0, 0.58169934640522891, 0.0),
    pymol_dope_color_name + "_9": (1.0, 0.27668845315904156, 0.0),
    pymol_dope_color_name + "_10": (0.8743315508021392, 0.0, 0.0)}

#------------------
# Domains colors. -
#------------------

domain_colors_dict = {'01_limegreen': (0.5, 1.0, 0.5),
                      '02_yellow': (1.0, 1.0, 0.0),
                      '03_red': (1.0, 0.0, 0.0),
                      '04_marineblue': (0.0, 0.5, 1.0),
                      '05_magenta': (1.0, 0.2, 0.8),
                      '06_yelloworange': (1.0, 0.87, 0.37),
                      '07_lightblue': (0.0, 0.8333333333333334, 1.0),
                      '08_lightgreen': (0.7558507273877292, 1.0, 0.2118912080961417),
                      '09_lightyellow': (1.0, 1.0, 0.9019607843137255),
                      '10_orange': (1.0, 0.5816993464052289, 0.0),
                      '11_violetpurple': (0.55, 0.25, 0.6),
                      '12_slate': (0.5, 0.5, 1.0),
                      '13_smudge': (0.55, 0.7, 0.4)}
domain_color_orderedkeys = list(domain_colors_dict.keys())
domain_color_orderedkeys.sort()
domain_colors_ordered = [(k, domain_colors_dict[k]) for k in domain_color_orderedkeys]

#-------------------------------
# Hydrophobicity scale colors. -
#-------------------------------

# Initiliazes the hydrophobicity scale colors in PyMOL.
pymol_polarity_color_name = "pymod_h"

# Kyte and Doolittle: J Mol Biol. 1982 May 5;157(1):105-32.
kyte_doolittle_h_dictionary = {
    pymol_polarity_color_name + '_A': (1.0, 0.59607843137254901, 0.59607843137254901),
    pymol_polarity_color_name + '_C': (1.0, 0.4392156862745098, 0.4392156862745098),
    pymol_polarity_color_name + '_E': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_D': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_G': (0.90980392156862744, 0.90980392156862744, 1.0),
    pymol_polarity_color_name + '_F': (1.0, 0.37647058823529411, 0.37647058823529411),
    pymol_polarity_color_name + '_I': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_H': (0.28235294117647058, 0.28235294117647058, 1.0),
    pymol_polarity_color_name + '_K': (0.13333333333333333, 0.13333333333333333, 1.0),
    pymol_polarity_color_name + '_M': (1.0, 0.5725490196078431, 0.5725490196078431),
    pymol_polarity_color_name + '_L': (1.0, 0.14901960784313728, 0.14901960784313728),
    pymol_polarity_color_name + '_N': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_Q': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_P': (0.64313725490196072, 0.64313725490196072, 1.0),
    pymol_polarity_color_name + '_S': (0.82352941176470584, 0.82352941176470584, 1.0),
    pymol_polarity_color_name + '_R': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_T': (0.84705882352941175, 0.84705882352941175, 1.0),
    pymol_polarity_color_name + '_W': (0.80000000000000004, 0.80000000000000004, 1.0),
    pymol_polarity_color_name + '_V': (1.0, 0.062745098039215685, 0.062745098039215685),
    pymol_polarity_color_name + '_Y': (0.71372549019607845, 0.71372549019607845, 1.0),
    pymol_polarity_color_name + '_X': (1.0, 1.0, 1.0)}

# Fauchere and Pliska: Eur. J. Med. Chem. 18:369-375(1983).
fauchere_pliska_h_scale = {
    pymol_polarity_color_name + '_A': (1.0, 0.27450980392156865, 0.27450980392156865),
    pymol_polarity_color_name + '_C': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_E': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_D': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_G': (0.36862745098039218, 0.36862745098039218, 1.0),
    pymol_polarity_color_name + '_F': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_I': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_H': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_K': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_M': (1.0, 0.21176470588235319, 0.21176470588235319),
    pymol_polarity_color_name + '_L': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_N': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_Q': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_P': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_S': (0.12549019607843137, 0.12549019607843137, 1.0),
    pymol_polarity_color_name + '_R': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_T': (0.18823529411764706, 0.18823529411764706, 1.0),
    pymol_polarity_color_name + '_W': (0.062745098039215685, 0.062745098039215685, 1.0),
    pymol_polarity_color_name + '_V': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_Y': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_X': (1.0, 1.0, 1.0)}

# This is the color dictionary actually used in PyMod.
polarity_color_dict = kyte_doolittle_h_dictionary


#-----------------------------
# Residue type color scheme. -
#-----------------------------

_residue_type_color_dict = {"HY": (1.0, 1.0, 1.0), # "#ffffff", Hydrophobic.
                            "BA": (0.12156862745098039, 0.4666666666666667, 0.7058823529411765), # "#1f77b4". Basic.
                            "AC": (0.8392156862745098, 0.15294117647058825, 0.1568627450980392), # "#d62728". Acid.
                            "PO": (1.0, 0.4980392156862745, 0.054901960784313725), # "#ff7f0e". Polar.
                            "AR": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), # "#2ca02c". Aromatic.
                            # "FL": (0.5803921568627451, 0.403921568627451, 0.7411764705882353), # "#9467bd". Flexible.
                            "FL": (0.7372549019607844, 0.7411764705882353, 0.13333333333333333), # "#bcbd22". Flexible.
                            "X": (0.4980392156862745, 0.4980392156862745, 0.4980392156862745), # "#7f7f7f". Other residues.
                            }

pymol_residue_type_color_name = "pymod_restype"

residue_type_color_dict = {
    pymol_residue_type_color_name + '_A': _residue_type_color_dict["HY"],
    pymol_residue_type_color_name + '_C': _residue_type_color_dict["HY"],
    pymol_residue_type_color_name + '_E': _residue_type_color_dict["AC"],
    pymol_residue_type_color_name + '_D': _residue_type_color_dict["AC"],
    pymol_residue_type_color_name + '_G': _residue_type_color_dict["FL"],
    pymol_residue_type_color_name + '_F': _residue_type_color_dict["AR"],
    pymol_residue_type_color_name + '_I': _residue_type_color_dict["HY"],
    pymol_residue_type_color_name + '_H': _residue_type_color_dict["AR"],
    pymol_residue_type_color_name + '_K': _residue_type_color_dict["BA"],
    pymol_residue_type_color_name + '_M': _residue_type_color_dict["HY"],
    pymol_residue_type_color_name + '_L': _residue_type_color_dict["HY"],
    pymol_residue_type_color_name + '_N': _residue_type_color_dict["PO"],
    pymol_residue_type_color_name + '_Q': _residue_type_color_dict["PO"],
    pymol_residue_type_color_name + '_P': _residue_type_color_dict["FL"],
    pymol_residue_type_color_name + '_S': _residue_type_color_dict["PO"],
    pymol_residue_type_color_name + '_R': _residue_type_color_dict["BA"],
    pymol_residue_type_color_name + '_T': _residue_type_color_dict["PO"],
    pymol_residue_type_color_name + '_W': _residue_type_color_dict["AR"],
    pymol_residue_type_color_name + '_V': _residue_type_color_dict["HY"],
    pymol_residue_type_color_name + '_Y': _residue_type_color_dict["AR"],
    pymol_residue_type_color_name + '_X': _residue_type_color_dict["X"]}


#---------------------------------------------
# Viridis colormap, adapted from matplotlib. -
#---------------------------------------------

viridis_colors = [
    [0.267004, 0.004874, 0.329415],
    [0.272594, 0.025563, 0.353093],
    [0.277018, 0.050344, 0.375715],
    [0.280267, 0.073417, 0.397163],
    [0.282327, 0.094955, 0.417331],
    [0.283197, 0.11568, 0.436115],
    [0.282884, 0.13592, 0.453427],
    [0.281412, 0.155834, 0.469201],
    [0.278826, 0.17549, 0.483397],
    [0.275191, 0.194905, 0.496005],
    [0.270595, 0.214069, 0.507052],
    [0.265145, 0.232956, 0.516599],
    [0.258965, 0.251537, 0.524736],
    [0.252194, 0.269783, 0.531579],
    [0.244972, 0.287675, 0.53726],
    [0.237441, 0.305202, 0.541921],
    [0.229739, 0.322361, 0.545706],
    [0.221989, 0.339161, 0.548752],
    [0.214298, 0.355619, 0.551184],
    [0.206756, 0.371758, 0.553117],
    [0.19943, 0.387607, 0.554642],
    [0.192357, 0.403199, 0.555836],
    [0.185556, 0.41857, 0.556753],
    [0.179019, 0.433756, 0.55743],
    [0.172719, 0.448791, 0.557885],
    [0.166617, 0.463708, 0.558119],
    [0.160665, 0.47854, 0.558115],
    [0.154815, 0.493313, 0.55784],
    [0.149039, 0.508051, 0.55725],
    [0.143343, 0.522773, 0.556295],
    [0.13777, 0.537492, 0.554906],
    [0.132444, 0.552216, 0.553018],
    [0.127568, 0.566949, 0.550556],
    [0.123463, 0.581687, 0.547445],
    [0.120565, 0.596422, 0.543611],
    [0.119423, 0.611141, 0.538982],
    [0.120638, 0.625828, 0.533488],
    [0.12478, 0.640461, 0.527068],
    [0.132268, 0.655014, 0.519661],
    [0.143303, 0.669459, 0.511215],
    [0.157851, 0.683765, 0.501686],
    [0.175707, 0.6979, 0.491033],
    [0.196571, 0.711827, 0.479221],
    [0.220124, 0.725509, 0.466226],
    [0.24607, 0.73891, 0.452024],
    [0.274149, 0.751988, 0.436601],
    [0.304148, 0.764704, 0.419943],
    [0.335885, 0.777018, 0.402049],
    [0.369214, 0.788888, 0.382914],
    [0.404001, 0.800275, 0.362552],
    [0.440137, 0.811138, 0.340967],
    [0.477504, 0.821444, 0.318195],
    [0.515992, 0.831158, 0.294279],
    [0.555484, 0.840254, 0.269281],
    [0.595839, 0.848717, 0.243329],
    [0.636902, 0.856542, 0.21662],
    [0.678489, 0.863742, 0.189503],
    [0.720391, 0.87035, 0.162603],
    [0.762373, 0.876424, 0.137064],
    [0.804182, 0.882046, 0.114965],
    [0.845561, 0.887322, 0.099702],
    [0.886271, 0.892374, 0.095374],
    [0.926106, 0.89733, 0.104071],
    [0.964894, 0.902323, 0.123941]
]

viridis_colors_hex = [convert_rgb_to_hex(c) for c in viridis_colors]


###################################################################################################
# PyMod elements features.                                                                        #
###################################################################################################

feature_types = ["domains", "sequence"]


###################################################################################################
# Standard bioinformatics dictionaries and lists.                                                 #
###################################################################################################

# Containers of the letters representing amminoacids in sequences.
prot_one_to_three_code = {
    "G": "GLY",
    "A": "ALA",
    "L": "LEU",
    "I" : "ILE",
    "R" : "ARG",
    "K" : "LYS",
    "M" : "MET",
    "C" : "CYS",
    "Y" : "TYR",
    "T" : "THR",
    "P" : "PRO",
    "S" : "SER",
    "W" : "TRP",
    "D" : "ASP",
    "E" : "GLU",
    "N" : "ASN",
    "Q" : "GLN",
    "F" : "PHE",
    "H" : "HIS",
    "V" : "VAL"}

# Code for the 20 standard amminoacids.
code_standard = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'
    }

# Code for standard nucleic acids residues.
nucleic_acids_dictionary = {
    # Deoxyribonucleotides as defined in PDB records.
    # ' DA':'a',' DT':'t',' DG':'g',' DC':'c',
    ' DA':'A',' DT':'T',' DG':'G',' DC':'C',
    # Ribonucleotides.
    # '  A':'a','  U':'u','  G':'g','  C':'c'
    '  A':'A','  U':'U','  G':'G','  C':'C'
    }
code_standard.update(nucleic_acids_dictionary)


###################################################################################################
# Updated data dictionaries for the new version.                                                  #
###################################################################################################

#------------
# Proteins. -
#------------

# One letter.
prot_standard_one_letter = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
prot_standard_one_letter_set = set(prot_standard_one_letter)
prot_standard_and_x_one_letter = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X")
prot_standard_and_x_one_letter_set = set(prot_standard_and_x_one_letter)
non_prot_standard_regex = "[^ACDEFGHIKLMNPQRSTVWY-]"
non_prot_standard_and_x_regex = "[^ACDEFGHIKLMNPQRSTVWYX-]"

# Three letters.
prot_standard_three_letters = ("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")
prot_standard_three_letters_set = set(prot_standard_three_letters)
prot_standard_and_x_three_letters = ("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","XXX")
prot_standard_and_x_three_letters_set = set(prot_standard_and_x_three_letters)

# Dictionaries.
prot_standard_and_x_one_to_three_dict = {
    "G": "GLY", "A": "ALA", "L": "LEU", "I": "ILE", "R": "ARG", "K": "LYS", "M": "MET", "C": "CYS",
    "Y": "TYR", "T": "THR", "P": "PRO", "S": "SER", "W": "TRP", "D": "ASP", "E": "GLU", "N": "ASN",
    "Q": "GLN", "F": "PHE", "H": "HIS", "V": "VAL", "X": "XXX"
}
prot_standard_and_x_three_to_one_dict = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'ILE': 'I', 'LEU': 'L', 'ASP': 'D',
    'GLU': 'E', 'LYS': 'K', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H', 'CYS': 'C',
    'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G', 'XXX': 'X'
}

def get_prot_one_to_three(one_letter_code):
    if one_letter_code in prot_standard_and_x_one_letter_set:
        return prot_standard_and_x_one_to_three_dict.get(one_letter_code)
    else:
        return "XXX"

def get_prot_three_to_one(three_letter_code):
    if three_letter_code in prot_standard_and_x_three_letters_set:
        return prot_standard_and_x_three_to_one_dict.get(three_letter_code)
    else:
        return "X"


#--------------------
# Heteroatoms data. -
#--------------------

std_amino_acid_backbone_atoms = set(("N", "CA", "C"))
mod_amino_acid_backbone_atoms = set(("N2", "C1", "C2"))

modified_residue_one_letter = "X"
ligand_one_letter = "x"
water_one_letter = "w"

# Heteroatoms to be shown as spheres in PyMOL.
sphere_hetres_names = set((" BA", " BR", " CA", " CD", " CO", " CU", " CL", "IOD",
                           "  K", " MG", " MN", " NA", " NI", " ZN",
                           "OXY"))
