# PyMod 3
PyMod 3 is an open source [PyMOL](https://github.com/schrodinger/pymol-open-source "PyMOL GitHub repository") plugin, designed to act as an interface between PyMOL and several bioinformatics tools (for example: BLAST+, HMMER, Clustal Omega, MUSCLE, PSIPRED and MODELLER). PyMod development has proceeded over the years with the aim of creating a simple, yet versatile tool for sequence/structure analysis and homology modeling within PyMOL. The current PyMod release, PyMod 3, is compatible with the most recent PyMOL versions (2.3 or higher) and has been extended with new functionalities. Starting from the amino acid sequence of a target protein, users may take advantage of PyMod 3 to carry out the three steps of the homology modeling process (that is, template searching, target-template sequence alignment and model building) in order to build a 3D atomic model of a target protein (or protein complex). Additionally, PyMod 3 may also be used outside the homology modeling context, in order to extend PyMOL with numerous types of  functionalities. Sequence similarity searches, multiple sequence/structure alignment building, evolutionary conservation analyses, domain parsing and loop modeling can be performed in the PyMod 3/PyMOL environment.

## Licence
PyMod 3 is licensed under the LGPL.

## Requirements
PyMod 3 runs on Windows, macOS and Linux versions of PyMOL. It is compatible with PyMOL version >= 2.3 and is supported by both incentive PyMOL builds (distributed by [Schrodinger](https://pymol.org/2/ "Schrodinger web-site")) and open source builds.

## Download
Download the latest PyMod 3 plugin file from [here](https://github.com/pymodproject/pymod/releases/download/v3.0/pymod3.zip "PyMod 3 plugin file for PyMOL")

## Installation
Please refer to the PyMod 3.0 User's Guide (see the link below) to learn how to install PyMod 3 on your system. Long story short: PyMod 3 can be installed as any other PyMOL plugin. Download the plugin ZIP file from the link above, and use the PyMOL plugin manager (**_Plugin -> Plugin Manager_** from the menu of its main window) to install it. The external tools of PyMod 3 can be obtained and configured through an easy-to-use installer dialog which can be launched by the plugin (**_Help -> Install PyMod Components_** from the main menu of PyMod 3). The way to configure MODELLER in PyMod may vary according to your PyMOL version (in the User's Guide we cover different scenarios).

## How to use PyMod
A series of introductory tutorials and in-depth information about all the PyMod 3 functionalities can be found in its User's Guide.

User's guide: download from [here](https://github.com/pymodproject/pymod/releases/download/v3.0/pymod_users_guide.pdf "PyMod 3 User's guide")

Home page (to be updated soon): http://schubert.bio.uniroma1.it/pymod/index.html

Old PyMod 2 packages (obsolete, they only work on PyMOL 1.X): https://github.com/pymodproject/pymod/releases/tag/v2.0

Contact: giacomo.janson@uniroma1.it
