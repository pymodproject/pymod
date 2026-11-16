# Copyright 2016 by Chengxin Zhang, Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

import os
import subprocess

from pymol import cmd, stored

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol


class Secondary_structure_assignment(PyMod_protocol):
    """
    Observed secondary structure assignment protocol.
    """

    protocol_name = "secondary_structure_assignment"

    def __init__(self, pymod, pymod_element):
        PyMod_protocol.__init__(self, pymod)
        self.pymod_element = pymod_element


    def assign_secondary_structure(self, algorithm="ksdssp"):

        if algorithm == "ksdssp":
            # If ksdssp is present, use it by default.
            if self.pymod.ksdssp.tool_file_exists():
                self.assign_with_ksdssp()
            # Otherwise use PyMOL built-in dss algorithm.
            else:
                self.assign_with_pymol_dss()
        elif algorithm == "dss":
            self.assign_with_pymol_dss()
        else:
            raise Exception("Unknown secondary structure assignment algorithm.")

        self.pymod_element.assigned_secondary_structure = True


    ###############################################################################################
    # Ksdssp.                                                                                     #
    ###############################################################################################

    def assign_with_ksdssp(self):
        # Runs ksdssp.
        ksdssp_fp = self.pymod.ksdssp["exe_file_path"].get_value()
        input_filepath = os.path.join(self.pymod.structures_dirpath, self.pymod_element.get_structure_file())
        dssptext = self.runKSDSSP(input_filepath, ksdssp_exe=ksdssp_fp)
        # Parses ksdssp's output, that is, an series of pdb format 'HELIX' and 'SHEET' record lines.
        dsspout = dssptext.split("\n")
        helices = set() # A set to store the sequence numbers of the residues in helical conformation.
        sheets = set() # A set to store the sequence numbers of the residues in sheet conformation.
        for line in dsspout:
            if line.startswith("HELIX"):
                new_residues_set = set(range(int(line[21:25]), int(line[33:37])+1))
                helices.update(new_residues_set)
            elif line.startswith("SHEET"):
                new_residues_set = set(range(int(line[22:26]), int(line[33:37])+1))
                sheets.update(new_residues_set)
        # Assigns to the PyMod element the observed secondaey structure observed using ksdssp.
        for residue in self.pymod_element.get_polymer_residues():
            if residue.db_index in helices:
                self.assign_sec_str_to_residue(residue, "H")
                rsel = residue.get_pymol_selector()
                cmd.alter(rsel,"ss='H'") # Set the residue new conformation in PyMOL.
            elif residue.db_index in sheets:
                self.assign_sec_str_to_residue(residue, "S")
                rsel = residue.get_pymol_selector()
                cmd.alter(rsel,"ss='S'") # Set the residue new conformation in PyMOL.
            else:
                self.assign_sec_str_to_residue(residue, "L")
                rsel = residue.get_pymol_selector()
                cmd.alter(rsel,"ss='L'") # Set the residue new conformation in PyMOL.
        # Update in PyMOL.
        cmd.rebuild()


    def runKSDSSP(self, PDB_file, output_file=None, energy_cutoff=-0.5,
                  minimum_helix_length=3, minimum_strand_length=3,
                  summary_file=None, ksdssp_exe="ksdssp"):
        """Command line wrapper for ksdssp,
        an implementation of the algorithm described in Wolfgang Kabsch and
        Christian Sander, "Dictionary of Protein Secondary Structure: Pattern
        Recognition of Hydrogen-Bonded and Geometrical Features," Biopolymers,
        22, 2577-2637 (1983).

        http://www.cgl.ucsf.edu/Overview/software.html#ksdssp


        Example:

        >>> PDB_file = "1TSR.pdb"
        >>> output_file = "1TSR.ksdssp"
        >>> ksdssp_cline = runKSDSSP(PDB_file, output_file)

        Arguments:
        PDB_file
            The input Protein Data  Bank  (PDB)  file  may  contain  any
            legal  PDB records.   Only ATOM records will be used.  All others
            are silently discarded.

        output_file (default: None)
            The  output  of  ksdssp  is a set of PDB HELIX and SHEET records.
            If no output_file argument is given, the records will be returned
            by runKSDSSP

        energy_cutoff (default -0.5)
            The default energy cutoff for defining hydrogen bonds as
            recommended  by Kabsch  and  Sander  is  -0.5  kcal/mol.

        minimum_helix_length (default 3)
            Normally,  HELIX records for helices of length three residues or
            greater are generated.  This option allows the user to change the
            minimum  helix length.

        minimum_strand_length (default 3)
            Normally,  SHEET records for strands of length three residues or
            greater are generated.  This option allows the user to change the
            minimum strand length.   Reducing the minimum strand length to 1 is
            not recommended, as there are bridges in many structures  that
            confuse  the  algorithm  for defining sheets.

        summary_file (default None)
            Normally,  ksdssp silently discards all the hydrogen-bonding
            information after generating the HELIX and SHEET records.  This
            option makes  ksdssp print  the  information to a file.

        ksdssp_exe (default 'ksdssp')
            location of KSDSSP executable
        """
        PDB_file_isfile=True
        if os.path.isfile(PDB_file)==False:
            # Assume PDB_file is the content of a PDB file
            PDB_file_isfile=False
            fp=open(".runKSDSSP.PDB_file.tmp",'w')
            print(PDB_file, file=fp)
            fp.close()
            PDB_file=".runKSDSSP.PDB_file.tmp"

        if not os.path.isfile(ksdssp_exe):
            print("Warning! cannot find KSDSSP executable!")
            print("Specify ksdssp_exe parameter to point to KSDSSP.")

        cline=[ksdssp_exe,
               "-c", str(energy_cutoff),
               "-h", str(minimum_helix_length),
               "-s", str(minimum_strand_length)]


        if summary_file != None:
            cline.append("-S")
            cline.append(summary_file)

        cline.append(PDB_file)

        if output_file != None:
            cline.append(output_file)
            try:
                return_code = subprocess.call(cline)
            except:
                return_code = ''
        else:
            fp=open(".runKSDSSP.std.tmp",'w')
            try:
                return_code = subprocess.call(cline, stdout=fp)
            except:
                return_code = ''
            fp.close()
            fp=open(".runKSDSSP.std.tmp",'r')
            return_code=fp.read()
            fp.close()

            try:
                os.remove(".runKSDSSP.std.tmp")
            except:
                print("Fail to remove temporary file .runKSDSSP.std.tmp")

        if PDB_file_isfile==False:
            try:
                os.remove(".runKSDSSP.PDB_file.tmp")
            except:
                print("Fail to remove temporary file .runKSDSSP.PDB_file.tmp'")
        return return_code


    ###############################################################################################
    # PyMOL dss.                                                                                  #
    ###############################################################################################

    def assign_with_pymol_dss(self):
        """
        Uses PyMOL's DSS algorithm to assign the secondary structure to a sequence according to atom
        coordinates of its PDB file.
        """
        selection = "object %s and n. CA" % self.pymod_element.get_pymol_selector()
        stored.resi_set = set()
        stored.temp_sec_str = []
        stored.pymol_info = []
        stored.pymod_resi_set = set([res.db_index for res in self.pymod_element.get_polymer_residues()])
        def include_sec_str_val(ca_tuple):
            if not ca_tuple[1] in stored.resi_set and ca_tuple[1] in stored.pymod_resi_set:
                stored.temp_sec_str.append(ca_tuple[0])
                stored.resi_set.add(ca_tuple[1])
                stored.pymol_info.append(ca_tuple)
        stored.include_val = include_sec_str_val
        cmd.iterate(selection, "stored.include_val((ss, resv))")
        sec_str_results = list(stored.temp_sec_str)

        if not (len(sec_str_results) == len(self.pymod_element.get_polymer_residues())):
            # raise Exception("Error in secondary structure assignment by PyMOL dss.")
            pass # TODO.
        list(map(lambda t: self.assign_sec_str_to_residue(t[0], t[1]), list(zip(self.pymod_element.get_polymer_residues(), sec_str_results))))

    def assign_sec_str_to_residue(self, res, ssr):
        res.secondary_structure = ssr
