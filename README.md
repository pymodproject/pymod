# PyMod 4
PyMod 4 is an open source [PyMOL](https://github.com/schrodinger/pymol-open-source "PyMOL GitHub repository") plugin that provides an integrated environment in which sequences, alignments and structures can be maintained and analyzed together. It acts as an interface between PyMOL and several bioinformatics tools and services (for example: BLAST+, HMMER, Clustal Omega, MUSCLE, PSIPRED, MODELLER, ESMFold2 and Boltz-2). PyMod 4 extends the PyMod 3 homology modeling framework (template searching, target-template alignment and model building) to current structure-generation, complex-modeling and protein-design workflows.

New in PyMod 4:
* **Nucleic acids**: DNA and RNA are supported as molecular objects, including comparative modeling of protein–nucleic acid complexes with MODELLER.
* **Structure retrieval and prediction**: precomputed models from the AlphaFold Protein Structure Database and ESM Atlas; template-free prediction with ESMFold2 and Boltz-2 (single proteins and protein–ligand complexes, with custom MSAs as input).
* **Structure and interface refinement**: clash detection combined with local optimization through Open Babel.
* **Phylogenetics**: maximum-likelihood trees with RAxML and IQ-TREE 2, visualized with the ETE toolkit.
* **Protein design**: structure-conditioned sequence design with ESM3 and binder design with BoltzProt-1.

Services accessed through an API (for example Boltz and ESM) require user-provided credentials.

## Requirements
PyMod 4 runs on Windows, macOS and Linux versions of PyMOL, with both incentive PyMOL builds distributed by [Schrodinger](https://pymol.org/ "Schrodinger website") and open source builds. It has been tested with the following PyMOL versions:

  | PyMOL version | Operating system |
  | --- | --- |
  | 3.1.9 | Ubuntu 24.04.4 LTS |
  | 3.1.9 | macOS Tahoe 26 |
  | 3.1.8 | Ubuntu 24.04.4 LTS |
  | 3.1.8 | macOS Tahoe 26 |
  | 3.1.4.1 | Ubuntu 24.04.4 LTS |
  | 3.0 | Ubuntu 24.04.4 LTS |
  | 2.5.8 | Ubuntu 24.04.4 LTS |
  | 2.5 | Ubuntu 24.04.4 LTS |


## Download
PyMod 4 plugin file (last update: 16/11/2026): [download from here](https://github.com/pymodproject/pymod/releases/download/v4.0/pymod4.zip "PyMod 4 plugin file for PyMOL")

## Installation
Please refer to the PyMod 4 [User's Guide](https://github.com/pymodproject/pymod/releases/download/v4.0/pymod_users_guide.pdf "PyMod 4 User's guide") to learn how to install PyMod 4 on your system. Long story short: PyMod 4 can be installed as any other PyMOL plugin. Download the PyMod 4 plugin file (a ZIP file named `pymod4.zip`, see the link above) and use the PyMOL plugin manager (**_Plugin -> Plugin Manager -> Install New Plugin_** from the menu of its main window) to install it. The external tools of PyMod 4 (BLAST+, HMMER, Clustal, MUSCLE, PSIPRED, ksdssp, IQ-TREE 2, RAxML), MODELLER and Open Babel can be obtained and configured through the installer dialog of the plugin (**_Help -> Install PyMod Components_** from its main menu). Installing MODELLER requires a MODELLER licence key.

## How to use PyMod
A series of introductory tutorials and in-depth information about all the PyMod 4 functionalities can be found in its User's Guide.

User's Guide: [download from here](https://github.com/pymodproject/pymod/releases/download/v4.0/pymod_users_guide.pdf "PyMod 4 User's guide")

PyMod 3 (previous version): [available from here](https://github.com/pymodproject/pymod/releases/tag/v3.0 "PyMod 3 release")

Old PyMod 2 packages (obsolete, they only work on PyMOL 1.X): if necessary [available from here](https://github.com/pymodproject/pymod/releases/tag/v2.0 "Legacy PyMod 2.0 installers")

Contact: serena.rosignoli@unimore.it, alessandro.paiardini@uniroma1.it

## Licence
PyMod 4 is licensed under the LGPL.
