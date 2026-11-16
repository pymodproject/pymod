# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Client for the ESM Atlas, the database of what ESMFold2 has already folded.

The companion of forge_client in this package, and its opposite: forge_client asks the
Forge API to predict a structure, with a key, while the Atlas holds the structures of
billions of proteins computed already and is free to read. Nothing is folded here, and a
sequence which is not in the Atlas has nothing to return.

It addresses a protein by the MD5 of its sequence, which has two consequences worth
knowing. The good one is that a hit covers exactly the sequence asked for -- there is no
alignment to make, nothing to trim and no numbering to reconcile. The other is that the
match is all or nothing: one residue of difference is a different hash and therefore a
miss, even for the same protein.

Reading the Atlas is all this module does. What PyMod does with what comes back -- the
element, the alignment, the structure, the cluster -- is atlas_protocols.
"""

import hashlib
import json
import re

import requests


# The alpha API of the Biohub platform, where the ESM Atlas moved. No key is needed to
# read it. It is documented at https://biohub.ai/esm/protein/atlas/api-docs/ and is
# declared experimental, so it is reached defensively.
ATLAS_BASE_URL = "https://biohub.ai/esm/protein/api/v1alpha1"

# Reading a record is a database lookup, not a computation, so it either answers quickly
# or something is wrong.
ATLAS_TIMEOUT = 60


def clean_sequence(sequence):
    """
    The sequence as the Atlas hashes it: the residues alone, in upper case, without the
    gaps of an alignment or any other decoration.
    """

    return re.sub(r'[^A-Za-z]', '', sequence or '').upper()


def atlas_hash(sequence):
    """The identifier of a protein in the Atlas: the MD5 of its sequence."""

    return hashlib.md5(clean_sequence(sequence).encode("utf-8")).hexdigest()


def atlas_url(sequence):
    return "%s/proteins/%s" % (ATLAS_BASE_URL, atlas_hash(sequence))


def record_url(record):
    """
    Where a record which has already been read lives.

    Built from the hash the record carries rather than from a sequence, because that is
    the only key the Atlas accepts: asking it for an accession answers 422, whichever
    accession it is -- even one the Atlas itself returned.
    """

    protein_hash = (record or {}).get("protein_hash")
    if not protein_hash:
        return None
    return "%s/proteins/%s" % (ATLAS_BASE_URL, protein_hash)


class Atlas_miss(Exception):
    """The sequence is not in the Atlas. Not an error: an answer."""


class Atlas_unavailable(Exception):
    """The Atlas could not be read, which says nothing about the sequence."""


def fetch_record(sequence, timeout=ATLAS_TIMEOUT, session=None):
    """
    The record of this sequence in the ESM Atlas.

    Raises Atlas_miss when the Atlas does not have it, and Atlas_unavailable when the
    service could not be read at all.
    """

    getter = (session or requests).get
    try:
        response = getter(atlas_url(sequence),
                          headers={"Accept": "application/json"},
                          timeout=timeout)
    except Exception as e:
        raise Atlas_unavailable("could not reach the ESM Atlas: %s" % e)

    if response.status_code == 404:
        raise Atlas_miss("the ESM Atlas has no protein with this exact sequence")
    if response.status_code != 200:
        raise Atlas_unavailable("the ESM Atlas answered %s" % response.status_code)

    try:
        record = response.json()
    except Exception:
        raise Atlas_unavailable("the ESM Atlas answered something which is not JSON")

    if not isinstance(record, dict) or not record.get("pdb"):
        raise Atlas_unavailable("the record of this protein carries no structure")
    return record


def pdb_text_of(record):
    """
    The coordinates of a record, with a HEADER line in front of them.

    The Atlas returns a file which begins with its licence, and the rest of PyMod expects
    a PDB file to begin with HEADER, so one is written here rather than relaxing the check
    everywhere else.
    """

    pdb = rescale_plddt(record["pdb"])
    if pdb.startswith("HEADER"):
        return pdb
    accession = record.get("accession") or record.get("protein_hash", "")[:12]
    header = "HEADER    PREDICTED STRUCTURE FROM THE ESM ATLAS    %s\n" % accession
    return header + pdb


def summary_of(record):
    """One line about a record, for the PyMOL console."""

    bits = []
    if record.get("accession"):
        bits.append("%s (%s)" % (record["accession"], record.get("source", "atlas")))
    if record.get("sequence_length"):
        bits.append("%s residues" % record["sequence_length"])
    if record.get("mean_plddt") is not None:
        # The Atlas reports the mean pLDDT from 0 to 1, while a PDB carries it from 0 to
        # 100 in the B-factor column.
        bits.append("mean pLDDT %.1f" % (float(record["mean_plddt"]) * 100))
    if record.get("ptm") is not None:
        bits.append("pTM %.2f" % float(record["ptm"]))
    if record.get("folded_on_demand"):
        bits.append("folded on demand")
    url = record_url(record)
    if url:
        bits.append(url)
    return ", ".join(bits)


# --------------------------------------------------------------------------------
# When the exact sequence is not in the Atlas.
# --------------------------------------------------------------------------------

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/%s.fasta"

# A UniProt accession, as it appears in the header of a sequence fetched by identifier.
# The boundaries are written by hand rather than with \b, because an underscore is a
# word character and headers in PyMod look like "Q9NZQ7_chain_A".
ACCESSION = re.compile(
    r'(?<![A-Za-z0-9])([OPQ][0-9][A-Z0-9]{3}[0-9]'
    r'|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?![A-Za-z0-9])')


def accession_in(text):
    """The UniProt accession named in a header, or None."""

    match = ACCESSION.search(text or "")
    return match.group(1) if match else None


def canonical_sequence(accession, timeout=ATLAS_TIMEOUT, session=None):
    """
    The canonical sequence of a UniProt entry, or None if it cannot be read.

    Used when the sequence in PyMod is not in the Atlas but is part of an entry which
    might be: a domain which was cut out of it, a mature chain, or a sequence which was
    edited. The Atlas is read by exact sequence, so the only way to ask it about a
    protein is to have that protein's exact sequence.
    """

    getter = (session or requests).get
    try:
        response = getter(UNIPROT_URL % accession, timeout=timeout)
        if response.status_code != 200:
            return None
        lines = response.text.splitlines()
    except Exception:
        return None
    return clean_sequence("".join(line for line in lines if not line.startswith(">")))


def locate_fragment(fragment, whole):
    """
    Where a fragment sits in a longer sequence, as (first, last) residue numbers counting
    from one, or None when it is not a contiguous piece of it.

    A domain cut out of a protein is a contiguous piece of it, which is what makes it
    possible to take the structure of the whole entry and keep only the part asked for.
    """

    fragment, whole = clean_sequence(fragment), clean_sequence(whole)
    if not fragment or not whole or fragment == whole:
        return None
    start = whole.find(fragment)
    if start < 0:
        return None
    return start + 1, start + len(fragment)


def trim_pdb(pdb_text, first, last):
    """
    The coordinates of residues 'first' to 'last' only, numbered as they were.

    The Atlas numbers the residues of a model from one, along the sequence of the entry,
    so a range of the sequence is the same range in the file.
    """

    kept = []
    for line in pdb_text.splitlines(True):
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 27:
            try:
                number = int(line[22:26])
            except ValueError:
                continue
            if not first <= number <= last:
                continue
        elif line.startswith(("TER", "CONECT", "MASTER")):
            continue
        kept.append(line)
    if not any(line.startswith("END") for line in kept):
        kept.append("END\n")
    return "".join(kept)


def rescale_plddt(pdb_text):
    """
    The B-factor column of an Atlas model, on the scale the rest of PyMod expects.

    The Atlas writes the per-residue confidence there, as ESMFold2 reports it, from 0 to 1.
    Everything which reads a pLDDT in PyMod -- the colouring of a model, the per-residue
    plot, the bins which say what is confident and what is not -- reads it from 0 to 100,
    the way a model of the AlphaFold database carries it. Left as it is, a model retrieved
    from the Atlas would be coloured as though none of it could be trusted.
    """

    lines = pdb_text.splitlines(True)
    values = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                values.append(float(line[60:66]))
            except ValueError:
                pass
    if not values or max(values) > 1.5:
        # Already on the 0 to 100 scale, or nothing to rescale.
        return pdb_text

    rescaled = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                value = float(line[60:66]) * 100
            except ValueError:
                rescaled.append(line)
                continue
            line = "%s%6.2f%s" % (line[:60], min(value, 999.99), line[66:])
        rescaled.append(line)
    return "".join(rescaled)
