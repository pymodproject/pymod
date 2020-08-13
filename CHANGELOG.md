# PyMod 3 Changelog
This file contains descriptions of changes of new PyMod 3 releases.

## Version 3.0.2 (17/8/2020)
- Added the possibility to import results from HH-suite output files (HHR and A3M).
- Added compatibility for pyside2 (thanks to Thomas Holder).

## Version 3.0.1 (15/6/2020)
- Fixed the PDB fetching functionality (updated to HTTPS urls on RCSB).

## Version 3.0.0 (22/5/2020)
- Updated the PyMod GUI using PyQt5.
- Homology modeling using MODELLER. Implemented various changes, some of them are:
    - Customize the optimization schedule of MODELLER.
    - Customize the MODELLER objective function terms.
    - Added the GA341 scoring function for model evaluation.
    - When building oligomeric models, added the SOAP-PP scoring function for interfaces evaluation.
    - Added sortable tables for model evaluation.
    - Added the possibility to use parallel jobs when building multiple models.
    - Solved a series of bugs which did not allow to use some PDB structures as templates.
- Added the possibility to perform loop modeling (with MODELLER) for elements that have a 3D structure loaded in PyMod/PyMOL.
- Added the SCR_FIND protocol (for analyzing the conservation in multiple structural alignments).
- Added the possibility to use the following HMMER programs:
    - phmmer, jackhmmer, hmmsearch (for searching protein sequence databases)
    - hmmscan (for searching profile HMM databases). hmmscan can be used as a tool to assign Pfam domain families.
- Added an easy-to-use dialog to install all PyMod external tools.
- Added an easy-to-use dialog to update the sequence databases for the external tools of PyMod.
- Added a contact/distance map viewer feature.
- Added a series of other small functionalities and fixes (please refer to the updated PyMod manual for a full list of them).

## Version 2.0.8 (2/11/2017)
- Introduced a series of fixes to add compatibility to PyMOL 2.0.

## Version 2.0.7 (2/12/2016)
- Fixed a bug that changed the atom order of structures when performing CE alignments with PyMOL.

## Version 2.0.6 (21/11/2016)
- When coloring structures by residue properties, all atom types of a residue are now colored with the same color, allowing better visualization of surfaces.

## Version 2.0.5 (13/11/2016)
- Added back the possibility to use the NCBI QBLAST server to perform a BLAST search after the NCBI switched its services to the HTTPS protocol.

## Version 2.0.4 (7/11/2016)
- Added the possibility to save and load PyMod projects.
