# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Reading what a Foldseek search returns, and naming it the way PyMod names things.

One hit carries everything PyMod needs to make an element of it:

    target      1crn-assembly1.cif.gz_A WATER STRUCTURE OF A HYDROPHOBIC ...
    tSeq        TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN
    qAln/dbAln  the alignment, already made
    tCa         16.967,12.784,4.338,13.856,...   the alpha carbons, in order
    eval prob seqId score alnLength qStartPos qEndPos dbStartPos dbEndPos
    taxId taxName

**The name is what matters most.** A PyMod element decides from its own header
whether a structure can be fetched for it: a header like '1crn_A' makes PyMod
read the code and the chain and fetch that entry from the PDB, and a header
carrying a UniProt accession makes it fetch the model of that protein. A hit
named the way PyMod names things therefore behaves like every other element of
the project -- the structure association which follows is the ordinary one, not
something written for this search.

So each hit is given:

    pymod_header    1crn_A, or AF-P01542-F1-model_v6, or the target as it came
    pdb_id/chain    for a hit of the PDB
    uniprot         for a hit of an AlphaFold database

and the coordinates are kept as they arrived, so a hit can be drawn without
anything being downloaded.
"""

import re


# What the databases of the server call their entries.
#
#   pdb100        1crn-assembly1.cif.gz_A <description>
#   afdb*         AF-P01542-F1-model_v6 <description>
#   cath50        1crnA00 / 3.30.100.10 ...
#   the others    an identifier of their own
PDB_TARGET = re.compile(r"^(?P<code>[0-9][A-Za-z0-9]{3})"
                        r"(?:-assembly\d+)?(?:\.cif|\.pdb)?(?:\.gz)?"
                        r"[._-](?P<chain>[A-Za-z0-9]{1,4})\b")
AFDB_TARGET = re.compile(r"^AF-(?P<accession>[A-Z0-9]+)-F\d+(?:-model)?(?:_v\d+)?",
                         re.IGNORECASE)
UNIPROT = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")


class FoldseekResultError(ValueError):
    """Raised when the answer of the server can not be read."""
    pass


class FoldseekHit:
    """One structure which looks like the query."""

    def __init__(self, database, record):
        self.database = database
        self.record = record or {}

        target = str(self.record.get("target") or "").strip()
        self.target = target
        # The server puts the description after the identifier, separated by a
        # space; both are worth keeping, for the name and for the window.
        pieces = target.split(None, 1)
        self.identifier = pieces[0] if pieces else target
        self.description = pieces[1] if len(pieces) > 1 else ""

        self.sequence = str(self.record.get("tSeq") or "")
        self.query_alignment = str(self.record.get("qAln") or "")
        self.hit_alignment = str(self.record.get("dbAln") or "")

        self.evalue = _number(self.record.get("eval"))
        self.probability = _number(self.record.get("prob"))
        self.sequence_identity = _number(self.record.get("seqId"))
        self.score = _number(self.record.get("score"))
        self.alignment_length = _number(self.record.get("alnLength"))
        self.query_start = _number(self.record.get("qStartPos"))
        self.query_end = _number(self.record.get("qEndPos"))
        self.hit_start = _number(self.record.get("dbStartPos"))
        self.hit_end = _number(self.record.get("dbEndPos"))
        self.query_length = _number(self.record.get("qLen"))
        self.hit_length = _number(self.record.get("dbLen"))
        self.taxon_name = str(self.record.get("taxName") or "")
        self.taxon_id = self.record.get("taxId")

        self.pdb_id, self.pdb_chain = _pdb_of(self.identifier)
        # Only a name of the form 'AF-<accession>-F1...' is a model which can be
        # fetched. The accession pattern on its own matches identifiers of other
        # databases as well -- an MGnify entry such as MGYP000936678158 matches
        # it -- and PyMod would then try to fetch a model which does not exist.
        self.uniprot = _uniprot_of(self.identifier)
        self.is_alphafold_model = bool(AFDB_TARGET.match(self.identifier or ""))

    # ------------------------------------------------------------------
    # What PyMod makes of it.
    # ------------------------------------------------------------------

    @property
    def pymod_header(self):
        """
        The name to give the element: the one PyMod itself can act on.

        '1crn_A' is read by PyMod as the entry 1crn, chain A, and a name holding
        a UniProt accession is read as a model which can be fetched. Anything
        else keeps the identifier the database used.
        """

        if self.pdb_id and self.pdb_chain:
            return "%s_%s" % (self.pdb_id.lower(), self.pdb_chain.upper())
        if self.uniprot:
            return self.identifier
        return self.identifier

    @property
    def structure_source(self):
        """Where a structure for this hit would come from: 'pdb', 'afdb' or None."""

        if self.pdb_id and self.pdb_chain:
            return "pdb"
        if self.is_alphafold_model:
            return "afdb"
        return None

    def alpha_carbons(self):
        """
        [(x, y, z)] of the hit, as the search returned them, or [].

        The server sends them as plain numbers, three per residue, in the order
        of the sequence. They are enough to draw the backbone of a hit without
        fetching anything.
        """

        text = str(self.record.get("tCa") or "").strip()
        if not text:
            return []
        try:
            numbers = [float(value) for value in text.split(",") if value.strip()]
        except ValueError:
            return []
        if len(numbers) < 3:
            return []
        return list(zip(numbers[0::3], numbers[1::3], numbers[2::3]))

    def identity_percent(self):
        """
        The identity of the alignment, as a percentage.

        The server reports it as a percentage for a structure search (100 for a
        hit against itself), while other searches of PyMod carry a fraction; the
        two are told apart by size rather than by trust.
        """

        value = self.sequence_identity
        if value is None:
            return None
        return value * 100 if value <= 1 else value

    def coverage_percent(self):
        """
        How much of the query the alignment covers, as a percentage.

        The one which says whether a hit answers the whole question or only a
        part of it: a strong hit over a third of the query is a hit on one
        domain, and its E-value does not say so.
        """

        if not self.query_length or self.query_start is None or self.query_end is None:
            return None
        covered = abs(self.query_end - self.query_start) + 1
        return 100.0 * covered / self.query_length

    def hit_coverage_percent(self):
        """
        The same for the hit: how much of the structure which was found takes
        part in the alignment. Low here and high in the query means the hit is a
        larger protein of which the query is one piece.
        """

        if not self.hit_length or self.hit_start is None or self.hit_end is None:
            return None
        covered = abs(self.hit_end - self.hit_start) + 1
        return 100.0 * covered / self.hit_length

    def __repr__(self):
        return "<FoldseekHit %s e=%s>" % (self.identifier, self.evalue)


def _number(value):
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if number == int(number) and abs(number) < 1e15 else number


def _pdb_of(identifier):
    """('1crn', 'A') for a hit of the PDB, or (None, None)."""

    match = PDB_TARGET.match(identifier or "")
    if not match:
        return None, None
    return match.group("code").lower(), match.group("chain").upper()


def _uniprot_of(identifier):
    """The UniProt accession of an AlphaFold hit, or None."""

    match = AFDB_TARGET.match(identifier or "")
    if match:
        return match.group("accession").upper()
    # Some databases name an entry by the accession alone.
    match = UNIPROT.search(identifier or "")
    return match.group(0) if match else None


def read_results(payload):
    """
    The hits of a search, as {database: [FoldseekHit]}, in the order the server
    returned them -- which is the order of the score.
    """

    if not isinstance(payload, dict):
        raise FoldseekResultError("The search returned something which is not a"
                                  " set of results.")

    entries = payload.get("results")
    if not isinstance(entries, list):
        raise FoldseekResultError("The search returned no results at all.")

    found = {}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        database = entry.get("db") or "?"
        alignments = entry.get("alignments") or []
        # One query gives a list of lists; the inner list is what belongs to it.
        if alignments and isinstance(alignments[0], list):
            flat = []
            for group in alignments:
                flat.extend(group or [])
        else:
            flat = alignments
        found[database] = [FoldseekHit(database, record) for record in flat
                           if isinstance(record, dict)]
    return found


def query_of(payload):
    """What the server understood the query to be: (name, sequence)."""

    queries = (payload or {}).get("queries") or []
    if not queries:
        return "", ""
    first = queries[0] or {}
    return str(first.get("header") or ""), str(first.get("sequence") or "")


def hit_count(results):
    return sum(len(hits) for hits in results.values())


def as_pdb_string(hit, chain="A"):
    """
    The backbone of a hit as a PDB file of alpha carbons, to be shown in PyMOL
    without fetching anything.

    The residue names come from the sequence the search returned, so the trace
    can be superposed on the query by sequence like any other structure.
    """

    coordinates = hit.alpha_carbons()
    if not coordinates:
        return ""

    lines = []
    for index, (x, y, z) in enumerate(coordinates, start=1):
        letter = hit.sequence[index - 1] if index <= len(hit.sequence) else "X"
        lines.append("ATOM  %5d  CA  %-3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                     "           C" % (index, THREE_LETTER.get(letter, "UNK"),
                                        chain, index, x, y, z))
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


THREE_LETTER = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
                "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
                "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
                "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
                "X": "UNK"}
