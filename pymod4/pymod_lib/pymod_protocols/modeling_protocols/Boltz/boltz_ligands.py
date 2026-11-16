# Copyright 2025 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
The ligands PyMod already has, found and described so that a Boltz-2 request can
carry them without anything being typed.

A ligand reaches the service in one of two notations: a CCD code, which the
service looks up in the chemical component dictionary of the PDB, or a SMILES
string, which describes the molecule itself. A ligand which is on screen in PyMod
can be given in either, and which one is right is not a matter of taste:

  * a residue whose name is a CCD code -- ATP, MG, 0G6, CFL -- is sent as that
    code, and the service builds the molecule from the dictionary, with its
    hydrogens, its bond orders and its charges as the PDB records them;

  * a residue whose name is not in the dictionary -- a ligand a program invented,
    LIG, UNL, a name from a docking run -- has no entry to look up, so the
    molecule is read from the coordinates in PyMOL and written as SMILES.

The order matters: the code is tried first, because it carries chemistry which
coordinates alone do not. The dictionary is asked over the network, once per code
and cached; without a network the code is offered unverified, which is what PyMod
did with every ligand before this module existed.

Nothing here decides anything for the user: the ligands are described, the window
lists them, and the user ticks the ones the prediction should include.
"""

import os
import re
import tempfile

from pymol import cmd

try:
    import requests
    has_requests = True
except ImportError:
    has_requests = False

try:
    from openbabel import openbabel as ob
    has_openbabel = True
except ImportError:
    has_openbabel = False


# The chemical component dictionary of the PDB, one component at a time. 200 is an
# entry, 404 is a name which is not in the dictionary.
CCD_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/%s"

# Water is not a ligand of a prediction, and a nucleic acid is a chain of its own
# in a Boltz request, not a ligand.
SKIPPED_CODES = frozenset(("HOH", "WAT", "DOD", "H2O"))
SKIPPED_HETRES_TYPES = frozenset(("water", "nucleic_acid"))

# What a CCD code can look like: 1 to 5 letters or digits. A name which cannot be
# a code is not worth asking the dictionary about.
CCD_PATTERN = re.compile(r"^[A-Z0-9]{1,5}$")

# Names which say "nobody named this". The dictionary has three of them -- UNK is
# "UNKNOWN", UNL is "Unknown ligand", UNX is "UNKNOWN ATOM OR ION" -- so asking
# the dictionary about them succeeds, and sending one as a code would ask the
# service to model an unknown residue instead of the molecule which is on the
# screen. A molecule read from an SDF is called UNK or UNL by default, which is
# how this was found. These are always read from their coordinates instead.
PLACEHOLDER_CODES = frozenset(("UNK", "UNL", "UNX", "XXX", "LIG", "DRG", "MOL"))


def is_a_real_code(code):
    """Whether a residue name is worth asking the chemical dictionary about."""

    code = (code or "").strip().upper()
    return bool(code) and code not in PLACEHOLDER_CODES and bool(CCD_PATTERN.match(code))

# How a ligand is given to the service.
KIND_CCD = "ccd"
KIND_CCD_UNVERIFIED = "ccd?"
KIND_SMILES = "smiles"
KIND_UNKNOWN = "unknown"

# What was asked of the dictionary already, so a window which is opened twice, or
# a structure which carries a ligand ten times, costs one request.
_ccd_cache = {}


def clear_ccd_cache():
    """Forget what the dictionary answered. For the tests, and for a retry."""
    _ccd_cache.clear()


def ccd_entry(code, timeout=4.0):
    """
    What the dictionary says about this code: a dict with its name and formula if
    it is an entry, False if the dictionary does not have it, and None if the
    dictionary could not be asked at all.
    """

    code = (code or "").strip().upper()
    if not code or not CCD_PATTERN.match(code):
        return False
    if code in _ccd_cache:
        return _ccd_cache[code]
    if not has_requests:
        return None

    try:
        answer = requests.get(CCD_URL % code, timeout=timeout)
    except Exception:
        # No network, a proxy in the way, a service which is down: unknown, and
        # not remembered, so that a later attempt can still find out.
        return None

    if answer.status_code == 404:
        _ccd_cache[code] = False
        return False
    if answer.status_code != 200:
        return None

    try:
        component = answer.json().get("chem_comp", {})
    except ValueError:
        return None

    entry = {"name": (component.get("name") or "").strip(),
             "formula": (component.get("formula") or "").strip()}
    _ccd_cache[code] = entry
    return entry


def smiles_of_selection(selection):
    """
    The canonical SMILES of what is in this PyMOL selection, or None.

    PyMOL has no SMILES writer -- there is no cmd.get_smiles, and it refuses to
    save a .smi -- so the molecule goes through Open Babel, which PyMod already
    asks for elsewhere (the optimizers) and installs from its own dialog.

    It goes through a MOL file rather than a PDB one. A PDB carries no bond
    orders, so Open Babel has to perceive the whole molecule from the geometry,
    and it gets it wrong often enough to matter: the inhibitor of 3UP7 comes out
    of the PDB route as

        OC(=O)c1ccc(cc1)N[C]1=NC=CC(=[N]=1)Nc1ccccc1C(=O)O

    with two atoms of broken valence, and out of the MOL route as

        OC(=O)c1ccc(cc1)Nc1nccc(n1)Nc1ccccc1C(=O)O

    which is the molecule. A MOL file carries the bonds and their orders as PyMOL
    knows them -- exactly, for a molecule built in its own builder -- and Open
    Babel then has nothing to guess. The PDB route is kept as a fallback for
    anything PyMOL will not write as a MOL.

    (.mol2 would do as well, but PyMOL writes formal charges into it -- the
    carboxylates of that same inhibitor come out deprotonated -- and that is a
    decision about the molecule which belongs to the user, not to PyMod.)
    """

    if not has_openbabel:
        return None
    try:
        if not cmd.count_atoms(selection):
            return None
    except Exception:
        return None

    for file_format in ("mol", "pdb"):
        smiles = _smiles_through(selection, file_format)
        if smiles:
            return smiles
    return None


def _quiet_openbabel():
    """Silences the warnings of Open Babel, and says what to put back."""

    try:
        previous = ob.obErrorLog.GetOutputLevel()
        ob.obErrorLog.SetOutputLevel(0)
        return previous
    except Exception:
        return None


def _restore_openbabel(previous):
    if previous is None:
        return
    try:
        ob.obErrorLog.SetOutputLevel(previous)
    except Exception:
        pass


def _quiet_pymol(quiet):
    """Turns the warnings of the molecule reader of PyMOL off, and on again."""

    try:
        cmd.feedback("disable" if quiet else "enable",
                     "objectmolecule", "warnings")
    except Exception:
        pass


def _smiles_through(selection, file_format):
    """The SMILES of a selection, written out as 'file_format' on the way."""

    handle, path = tempfile.mkstemp(suffix="." + file_format,
                                    prefix="pymod_ligand_")
    os.close(handle)
    try:
        try:
            # PyMOL writes a warning to the shell for a residue it cannot
            # kekulise. The scan reads a lot of residues and the warnings are
            # about none of the user's doing, so they are kept quiet here.
            _quiet_pymol(True)
            try:
                cmd.save(path, selection, quiet=1)
            finally:
                _quiet_pymol(False)
        except Exception:
            return None
        if not os.path.getsize(path):
            return None

        try:
            # Open Babel says what it thinks of a molecule on the way through --
            # "Failed to kekulize aromatic bonds", and the like -- and the scan
            # reads every residue of every structure open. The reading is
            # answered by whether a SMILES comes out of it, not by the grumbling,
            # so the grumbling is kept out of the user's shell.
            previous = _quiet_openbabel()
            try:
                conversion = ob.OBConversion()
                conversion.SetInAndOutFormats(file_format, "can")
                molecule = ob.OBMol()
                if not conversion.ReadFile(molecule, path):
                    return None
                if not molecule.NumAtoms():
                    return None
                written = conversion.WriteString(molecule).strip()
            finally:
                _restore_openbabel(previous)
        except Exception:
            return None
    finally:
        try:
            os.remove(path)
        except OSError:
            pass

    # 'can' writes the SMILES and, after a tab, the title of the molecule.
    smiles = written.split("\t")[0].split()[0] if written else ""
    return smiles or None


class DetectedLigand:
    """One ligand PyMod has, and how a request should carry it."""

    SOURCE_PYMOD = "pymod"      # a ligand of a structure of the project.
    SOURCE_PYMOL = "pymol"      # a molecule which is open in PyMOL and nowhere else.

    def __init__(self, code, kind, value, element_name, chain, resi,
                 occurrences=1, name="", formula="", note="",
                 source=SOURCE_PYMOD, selector=""):
        # Where it was found, and what PyMod would send before the user edits it.
        self.source = source
        self.selector = selector
        self.detected_value = value
        self.code = code
        self.kind = kind                  # KIND_CCD, KIND_CCD_UNVERIFIED, ...
        self.value = value                # the code, or the SMILES string
        self.element_name = element_name
        self.chain = chain
        self.resi = resi
        self.occurrences = occurrences
        self.name = name                  # what the dictionary calls it
        self.formula = formula
        self.note = note

    @property
    def usable(self):
        return self.kind != KIND_UNKNOWN and bool(self.value)

    def notation(self):
        return {KIND_CCD: "CCD code",
                KIND_CCD_UNVERIFIED: "CCD code, not verified",
                KIND_SMILES: "SMILES"}.get(self.kind, "not usable")

    def label(self):
        """
        What to call it in the window: the name it goes by, and how many there
        are of it. Short on purpose -- where it came from and what it is are
        columns of their own, and the SMILES is in the field beside it, so
        writing any of them here only said everything twice.
        """

        # A molecule of PyMOL called UNK or UNL -- what a file calls a molecule
        # nobody named -- is known by the object it is in. One which carries a
        # real name keeps it: a structure opened in PyMOL can hold a dozen
        # ligands, and "PEG", "ACY", "MG" tell them apart where the name of the
        # object they share does not.
        if self.source == self.SOURCE_PYMOL and not is_a_real_code(self.code):
            text = self.element_name or "?"
            if text == "the PyMOL selection":
                text = "selection"
        else:
            text = self.code or "?"
        if self.occurrences > 1:
            text += " ×%d" % self.occurrences
        return text

    def where(self):
        """The second column: where this molecule was found."""

        if self.source == self.SOURCE_PYMOL:
            if not is_a_real_code(self.code):
                # The name column is already the object it is in.
                return "open in PyMOL"
            where = self.element_name or "PyMOL"
            if self.chain:
                where += " %s" % self.chain
                if self.resi:
                    where += "/%s" % self.resi
            elif self.resi:
                where += " %s" % self.resi
            return where

        where = self.element_name or "?"
        if self.chain:
            where += " %s" % self.chain
            if self.resi:
                where += "/%s" % self.resi
        elif self.resi:
            where += " %s" % self.resi
        return where

    def describes(self):
        """
        The third column, when there is no SMILES field in it: what the
        dictionary calls a code, or why a ligand cannot be sent at all.
        """

        if self.kind == KIND_UNKNOWN:
            return "cannot be sent"
        name = self.name or ""
        return name if len(name) <= 52 else name[:49] + "..."

    def edited(self, text):
        """
        The SMILES as the user wrote it. A molecule read from coordinates is a
        reading, not a fact -- a bond order or a stereocentre can come out wrong
        -- so the window lets it be corrected, and what is sent is what is shown.
        """

        self.value = (text or "").strip()

    def was_edited(self):
        return self.value != self.detected_value

    def tooltip(self):
        lines = ["%s, sent as %s" % (self.code, self.notation())]
        if self.source == self.SOURCE_PYMOL:
            lines.append("Open in PyMOL, and not part of the project: PyMod read"
                         " the molecule from its coordinates.")
        if self.name:
            lines.append(self.name)
        if self.formula:
            lines.append("Formula: %s" % self.formula)
        if self.kind == KIND_SMILES:
            lines.append("SMILES built from the coordinates in PyMOL:")
            lines.append(self.value)
        if self.note:
            lines.append(self.note)
        return "\n".join(lines)

    def __repr__(self):
        return "<DetectedLigand %s %s>" % (self.code, self.kind)


def _structure_elements(pymod, first=()):
    """
    The elements of the session which have a structure, the ones given first at
    the front. get_pymod_elements_list() returns clusters and alignments as well,
    and those have no structure to look into.
    """

    # A PyMod which cannot list its elements -- a stub, a caller which passed
    # something else -- still gets the elements it named.
    listing = getattr(pymod, "get_pymod_elements_list", None)
    everything = []
    if callable(listing):
        try:
            everything = list(listing() or [])
        except Exception:
            everything = []

    ordered, seen = [], set()
    for element in list(first) + everything:
        if id(element) in seen:
            continue
        seen.add(id(element))
        try:
            if element.has_structure():
                ordered.append(element)
        except Exception:
            continue
    return ordered


def _heteroresidues_of(element):
    """The heteroresidues of an element, waters and nucleic acids left out."""

    try:
        residues = element.get_heteroresidues() or []
    except Exception:
        return []

    kept = []
    for residue in residues:
        code = (getattr(residue, "three_letter_code", "") or "").strip().upper()
        if not code or code in SKIPPED_CODES:
            continue
        if getattr(residue, "hetres_type", None) in SKIPPED_HETRES_TYPES:
            continue
        kept.append((code, residue))
    return kept


def _objects_of_the_project(pymod):
    """The names of the PyMOL objects which belong to a structure of the project."""

    owned = set()
    for element in _structure_elements(pymod):
        try:
            selector = element.get_pymol_selector()
        except Exception:
            continue
        if selector:
            owned.add(str(selector).split(" and ")[0].strip("() "))
    return owned


def _residues_of(selection):
    """(chain, resi, resn) of everything in this selection, in PyMOL's order."""

    found = []
    seen = set()
    try:
        cmd.iterate(selection, "found.append((chain, resi, resn))",
                    space={"found": found})
    except Exception:
        return []

    ordered = []
    for residue in found:
        if residue in seen:
            continue
        seen.add(residue)
        ordered.append(residue)
    return ordered


def _molecules_of(selection):
    """
    The molecules of a selection, one selector each.

    A molecule is what is bonded together, which PyMOL knows and answers with
    'bymolecule'. Grouping by residue instead is wrong twice over: two molecules
    of one file often carry the same residue number, which made each of them
    select both and be offered as the other, and a ligand bound covalently to a
    second residue would be cut in half.
    """

    atoms = []
    try:
        cmd.iterate(selection, "atoms.append((model, index))",
                    space={"atoms": atoms})
    except Exception:
        return []

    molecules, assigned = [], set()
    for model, index in atoms:
        if (model, index) in assigned:
            continue
        group = []
        try:
            cmd.iterate("bymolecule (%s and model %s and index %s)"
                        % (selection, model, index),
                        "group.append((model, index))", space={"group": group})
        except Exception:
            group = [(model, index)]
        if not group:
            group = [(model, index)]

        assigned.update(group)
        indices = sorted(set(atom_index for _, atom_index in group))
        molecules.append("(%s) and model %s and index %s"
                         % (selection, model,
                            "+".join(str(number) for number in indices)))
    return molecules


def _describe_molecule(selection, where, ask_dictionary=True, timeout=4.0):
    """
    One molecule of PyMOL, described as a request should carry it.

    The rule is the one the ligands of the project follow: a name which the
    chemical component dictionary has is sent as that code, because a code
    carries chemistry which coordinates do not -- a magnesium read from its
    coordinates comes out as '[MgH2]', with hydrogens Open Babel had to invent,
    while 'MG' is the ion itself. Anything else is read from the coordinates.
    """

    residues = _residues_of(selection)
    if not residues:
        return None
    chain, resi, _ = residues[0]
    names = set((resn or "").strip().upper() for _, _, resn in residues)
    code = names.pop() if len(names) == 1 else ""
    if code in SKIPPED_CODES:
        return None

    if is_a_real_code(code) and ask_dictionary:
        answer = ccd_entry(code, timeout=timeout)
        if isinstance(answer, dict) and _atoms_agree_with(selection,
                                                          answer.get("formula", "")):
            return DetectedLigand(
                code, KIND_CCD, code, where, chain, resi,
                source=DetectedLigand.SOURCE_PYMOL, selector=selection,
                name=answer.get("name", ""), formula=answer.get("formula", ""),
                note="Open in PyMOL, and a name the chemical component"
                     " dictionary has, with the atoms to match: it is sent as"
                     " that code.")

    smiles = smiles_of_selection(selection)
    if not smiles:
        return None
    return DetectedLigand(code or "UNL", KIND_SMILES, smiles, where, chain, resi,
                          source=DetectedLigand.SOURCE_PYMOL, selector=selection,
                          note=_note_for_molecule(where, smiles))


def pymol_molecules(pymod, selection_name="sele", ask_dictionary=True,
                    timeout=4.0):
    """
    The molecules which are open in PyMOL and are not part of the project.

    PyMod knows the ligands of the structures it holds; it knows nothing of what
    is drawn in the builder of PyMOL, opened there by hand, or left behind by a
    docking run. Those are molecules a prediction may well want, and there is no
    reason to make anyone type their SMILES by hand when PyMOL has the atoms.

    Everything which is not polymer and not water is offered, one residue at a
    time, read through Open Babel. Nothing is decided here: they are described,
    and the window lists them unticked -- a session holds a good deal which was
    never meant for a prediction.
    """

    if not has_openbabel:
        return []

    try:
        objects = list(cmd.get_object_list() or [])
    except Exception:
        return []

    owned = _objects_of_the_project(pymod)
    molecules = []

    for name in objects:
        if name in owned:
            continue
        # A structure open in PyMOL alone still has its ligands taken from it,
        # and its protein left alone: a chain is a chain of the request, not a
        # ligand of it.
        selection = "(%s) and not polymer and not solvent" % name
        try:
            if not cmd.count_atoms(selection):
                continue
        except Exception:
            continue

        for molecule_selection in _molecules_of(selection):
            ligand = _describe_molecule(molecule_selection, name,
                                        ask_dictionary=ask_dictionary,
                                        timeout=timeout)
            if ligand is None:
                continue
            # Every copy is offered: two ions of one structure are two ligands
            # of the request, and which of them is wanted is not for PyMod to
            # decide.
            molecules.append(ligand)

    # What is selected is often one of the molecules just listed -- it was
    # clicked on to look at it -- and offering it twice helps nobody.
    # What is selected is usually one of the molecules just listed -- it was
    # clicked on to look at it -- and offering it again helps nobody.
    already = set(molecule.selector for molecule in molecules)
    for molecule in _molecule_of_selection(pymod, selection_name, owned,
                                           ask_dictionary=ask_dictionary,
                                           timeout=timeout):
        if molecule.selector in already:
            continue
        if any(molecule.value == existing.value for existing in molecules):
            continue
        molecules.append(molecule)

    return molecules


def _atoms_agree_with(selection, formula):
    """
    Whether the molecule in this selection is the one the dictionary describes,
    counted by its heavy atoms.

    The name of a molecule in PyMOL is worth much less than the name of a ligand
    in a structure file: the fragments of the builder carry names of their own --
    'butane' arrives as MET, which is methionine in the dictionary -- and sending
    that code would send a different molecule altogether. The name is believed
    only when the atoms say the same thing. Hydrogens are left out of it: a
    structure often has none, and the dictionary always does.
    """

    # The case of the symbol is what tells 'Mg' from 'M' and 'G', so the formula
    # is read as it is written: "C10 H16 N5 O13 P3", "Mg 1".
    wanted = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)",
                                     (formula or "").replace(" ", "")):
        symbol = element.upper()
        if symbol == "H":
            continue
        wanted[symbol] = wanted.get(symbol, 0) + int(count or 1)
    if not wanted:
        return False

    counted = {}
    try:
        cmd.iterate(selection, "counted[elem.upper()] = counted.get(elem.upper(), 0) + 1",
                    space={"counted": counted})
    except Exception:
        return False
    counted.pop("H", None)
    return counted == wanted


def _note_for_molecule(name, smiles):
    note = ("'%s' is open in PyMOL and is not part of the project, so the"
            " molecule was read from its coordinates." % name)
    if "." in smiles:
        note += (" It reads as more than one molecule: check it, and send one"
                 " ligand at a time.")
    return note


def _molecule_of_selection(pymod, selection_name, owned, ask_dictionary=True,
                           timeout=4.0):
    """
    What is selected in PyMOL right now, if anything is: the way a molecule just
    drawn in the builder is offered without having to be saved or named first.
    """

    if not selection_name:
        return []
    try:
        if selection_name not in (cmd.get_names("selections") or []):
            return []
        if not cmd.count_atoms(selection_name):
            return []
    except Exception:
        return []

    selection = "(%s) and not polymer and not solvent" % selection_name
    try:
        if not cmd.count_atoms(selection):
            return []
    except Exception:
        return []

    described = []
    for molecule_selection in _molecules_of(selection):
        ligand = _describe_molecule(molecule_selection, "the PyMOL selection",
                                    ask_dictionary=ask_dictionary,
                                    timeout=timeout)
        if ligand is not None:
            ligand.note = ("What is selected in PyMOL, read from its"
                           " coordinates.") if ligand.kind == KIND_SMILES \
                else ligand.note
            described.append(ligand)
    return described


def detect_ligands(pymod, selected_elements=(), timeout=4.0, ask_dictionary=True,
                   include_pymol_molecules=True):
    """
    The ligands of every structure open in PyMod, described for a Boltz request.

    Each ligand is described once, however many copies of it there are and however
    many structures carry it, and the first structure which had it is named. The
    ones of the selected sequences come first, because they are the ones the user
    is predicting around.

    A code the dictionary confirms is sent as a code; a code the dictionary does
    not have is written as SMILES from its coordinates; a code which could not be
    checked -- no network -- is offered as a code with that said. When the
    dictionary cannot be reached at all, it is not asked again for this scan: a
    session without a network would otherwise wait for every ligand in turn.
    """

    # One entry per residue, not one per name: two copies of a cofactor are two
    # ligands of the request, in two places, and each is ticked on its own.
    entries = []
    for element in _structure_elements(pymod, selected_elements):
        element_name = getattr(element, "my_header", None) or getattr(element, "name", "?")
        for code, residue in _heteroresidues_of(element):
            entries.append({"code": code, "element": element, "residue": residue,
                            "element_name": element_name, "occurrences": 1})

    dictionary_reachable = ask_dictionary
    ligands = []

    for entry in entries:
        code = entry["code"]
        residue = entry["residue"]
        chain = getattr(residue, "chain_id", "") or ""
        resi = getattr(residue, "db_index", "") or ""

        if not is_a_real_code(code):
            # A name which names nothing: read the molecule instead of looking
            # up a code which would mean "unknown" to the service.
            answer = False
        elif dictionary_reachable:
            answer = ccd_entry(code, timeout=timeout)
            if answer is None:
                # The dictionary could not be reached. Stop asking, and offer
                # every remaining code of this scan as it stands.
                dictionary_reachable = False
        else:
            answer = None

        if isinstance(answer, dict):
            ligands.append(DetectedLigand(
                code, KIND_CCD, code, entry["element_name"], chain, resi,
                occurrences=entry["occurrences"],
                name=answer.get("name", ""), formula=answer.get("formula", "")))
            continue

        if answer is None:
            note = ("The chemical component dictionary could not be reached, so"
                    " this code was not verified.")
            if not is_a_real_code(code):
                note = ("'%s' cannot be a CCD code, and no SMILES could be built"
                        " for it." % code)
                ligands.append(DetectedLigand(code, KIND_UNKNOWN, "",
                                              entry["element_name"], chain, resi,
                                              occurrences=entry["occurrences"],
                                              note=note))
                continue
            ligands.append(DetectedLigand(code, KIND_CCD_UNVERIFIED, code,
                                          entry["element_name"], chain, resi,
                                          occurrences=entry["occurrences"],
                                          note=note))
            continue

        # The dictionary does not have this name: read the molecule instead.
        selector = None
        try:
            selector = residue.get_pymol_selector()
        except Exception:
            selector = None
        smiles = smiles_of_selection(selector) if selector else None

        if smiles:
            if is_a_real_code(code):
                why = ("'%s' is not in the chemical component dictionary, so the"
                       " molecule was read from its coordinates." % code)
            else:
                why = ("'%s' names nothing in particular, so the molecule was"
                       " read from its coordinates instead of being looked up."
                       % code)
            ligands.append(DetectedLigand(
                code, KIND_SMILES, smiles, entry["element_name"], chain, resi,
                occurrences=entry["occurrences"], note=why))
        else:
            if not is_a_real_code(code):
                note = ("'%s' cannot be sent as a code and its molecule could"
                        " not be read from PyMOL." % code)
            elif not has_openbabel:
                note = ("'%s' is not in the chemical component dictionary and"
                        " Open Babel is not installed, so no SMILES could be"
                        " built for it." % code)
            else:
                note = ("'%s' is not in the chemical component dictionary and its"
                        " molecule could not be read from PyMOL." % code)
            ligands.append(DetectedLigand(code, KIND_UNKNOWN, "",
                                          entry["element_name"], chain, resi,
                                          occurrences=entry["occurrences"],
                                          note=note))

    # And what is open in PyMOL without being part of the project, which PyMod
    # would otherwise be blind to. After the ligands of the project: those are
    # the ones the user is predicting around.
    if include_pymol_molecules:
        ligands.extend(pymol_molecules(pymod, ask_dictionary=dictionary_reachable,
                                       timeout=timeout))

    return ligands
