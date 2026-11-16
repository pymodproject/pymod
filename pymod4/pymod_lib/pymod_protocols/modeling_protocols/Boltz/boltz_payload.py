# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Builds the request body of a Boltz "structure and binding" prediction.

This module is the ONLY place which knows the shape of the request, so that a
change in the API costs one file. It holds no GUI and no PyMod code: it works on
plain strings and dictionaries, which makes it testable on its own and lets the
same functions serve the protocol and its tests.

The schema below was verified against the live API (with the free
'estimate-cost' endpoint, which validates a request without running it):

    {"model": "boltz-2.1",
     "input": {
        "entities": [
            {"type": "protein" | "dna" | "rna" | "ligand_ccd" | "ligand_smiles",
             "chain_ids": ["A", ...],       # a list: one entity can fill several chains
             "value": "<sequence | CCD code | SMILES>",
             "msa": <see below>}            # proteins only
        ],
        "num_samples": 1..10,
        "binding": {"binder_chain_id": "L", "type": "ligand_protein_binding"}}}

Note that the sequence field is called 'value', that a ligand is not of type
'ligand' but of 'ligand_ccd' or 'ligand_smiles', and that the whole request is
wrapped in 'input'.

The MSA of a protein entity is either

    omitted                 -> the service generates one (the default),
    {"type": "empty"}       -> single sequence prediction,
    {"type": "custom", "format": "a3m", "source": {...}}

where a custom MSA may be sent inline, with no need to publish it anywhere:

    {"type": "base64", "data": "<base64>", "media_type": "text/x-a3m"}

The media type is mandatory and exact; 'text/plain' is refused by the service.

Two rules constrain the MSAs of a request, both enforced here so that a job is
not submitted only to be rejected:

  1. automatic MSA generation can not be mixed with custom MSAs: as soon as one
     protein carries a custom MSA, every other protein must carry a custom or an
     empty one;
  2. proteins with the same sequence must share the same MSA setting.
"""

import base64
import re


###################################################################################################
# Constants of the API.                                                                           #
###################################################################################################

# Model served by the structure and binding resource.
MODEL = "boltz-2.1"

TYPE_PROTEIN = "protein"
TYPE_DNA = "dna"
TYPE_RNA = "rna"
TYPE_LIGAND_CCD = "ligand_ccd"
TYPE_LIGAND_SMILES = "ligand_smiles"

POLYMER_TYPES = (TYPE_PROTEIN, TYPE_DNA, TYPE_RNA)
LIGAND_TYPES = (TYPE_LIGAND_CCD, TYPE_LIGAND_SMILES)

# The only kind of binding this model accepts. Protein-protein binding is
# documented but rejected by 'boltz-2.1', so it is deliberately not offered.
BINDING_LIGAND_PROTEIN = "ligand_protein_binding"

# Number of predicted samples the service allows.
MIN_SAMPLES = 1
MAX_SAMPLES = 10

# Media types the service requires for an inline custom MSA.
MSA_MEDIA_TYPES = {"a3m": "text/x-a3m",
                   "csv": "text/csv"}

# How the MSA of a protein is decided. 'automatic' leaves the key out.
MSA_AUTOMATIC = "automatic"
MSA_EMPTY = "empty"
MSA_CUSTOM = "custom"

# Residues accepted in a nucleic acid sequence, as documented by the service.
DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")

# Characters which are gaps in a PyMod alignment.
GAP_CHARS = "-."


class BoltzPayloadError(ValueError):
    """
    Raised when a request can not be built. The message is meant to be shown to
    the user as it is, so it always says which entity or chain is at fault.
    """
    pass


###################################################################################################
# Chain identifiers.                                                                              #
###################################################################################################

def chain_id_sequence(count):
    """
    Returns 'count' chain identifiers: A..Z, then AA, AB and so on.

    PyMod assigns these itself rather than reading them from a loaded structure,
    because a sequence being modelled usually has no structure at all (which is
    the whole point of predicting one). Taking them from 'element.structure'
    raised an AttributeError for exactly the common case.
    """

    if count < 0:
        raise BoltzPayloadError("A negative number of chains was requested.")

    ids = []
    letters = [chr(c) for c in range(ord("A"), ord("Z") + 1)]
    i = 0
    while len(ids) < count:
        name = ""
        n = i
        while True:
            name = letters[n % 26] + name
            n = n // 26 - 1
            if n < 0:
                break
        ids.append(name)
        i += 1
    return ids


###################################################################################################
# Sequences and ligand values.                                                                     #
###################################################################################################

def clean_sequence(sequence):
    """
    Removes the gaps and the whitespace a PyMod sequence may carry, and upper
    cases it. Ligand values are NOT passed through this: a SMILES string is case
    sensitive.
    """

    if sequence is None:
        return ""
    text = "".join(sequence.split())
    for gap in GAP_CHARS:
        text = text.replace(gap, "")
    return text.upper()


def guess_ligand_type(value):
    """
    Guesses whether a ligand was given as a CCD code or as a SMILES string, and
    returns the corresponding entity type.

    The two are NOT always distinguishable: 'CCO' is both a plausible CCD code
    and the SMILES of ethanol. A short alphanumeric string is therefore read as a
    CCD code, whatever its case, since a user typing 'atp' means the ligand ATP,
    and anything else as SMILES. Callers which know better should pass the type
    explicitly; this function only exists to pre-select the right option in the
    GUI, where the choice remains visible and can be corrected.
    """

    text = (value or "").strip().upper()
    if re.match(r"^[A-Z0-9]{1,5}$", text):
        return TYPE_LIGAND_CCD
    return TYPE_LIGAND_SMILES


###################################################################################################
# MSA descriptors.                                                                                 #
###################################################################################################

def msa_automatic():
    """The service builds the MSA. Encoded by leaving the key out entirely."""
    return None


def msa_empty():
    """Single sequence prediction, with no MSA at all."""
    return {"type": "empty"}


def msa_custom(text, msa_format="a3m"):
    """
    A custom MSA sent inline, base64 encoded, so that an alignment computed
    locally by PyMod can be used without publishing it anywhere.
    """

    if msa_format not in MSA_MEDIA_TYPES:
        raise BoltzPayloadError("Unsupported MSA format %r. Use one of: %s."
                                % (msa_format, ", ".join(sorted(MSA_MEDIA_TYPES))))
    if not text or not text.strip():
        raise BoltzPayloadError("The custom MSA is empty.")

    data = base64.b64encode(text.encode("utf-8")).decode("ascii")
    return {"type": "custom",
            "format": msa_format,
            "source": {"type": "base64",
                       "data": data,
                       "media_type": MSA_MEDIA_TYPES[msa_format]}}


###################################################################################################
# Entities.                                                                                        #
###################################################################################################

def polymer_entity(entity_type, chain_ids, sequence, msa=None):
    """
    Builds a protein, DNA or RNA entity. 'msa' is only meaningful for a protein
    and must be None for the others, since the service does not accept the key
    on a nucleic acid.
    """

    if entity_type not in POLYMER_TYPES:
        raise BoltzPayloadError("Unknown polymer type %r." % entity_type)

    chain_ids = _check_chain_ids(chain_ids)
    value = clean_sequence(sequence)
    if not value:
        raise BoltzPayloadError("Chain %s has an empty sequence."
                                % ", ".join(chain_ids))

    if entity_type == TYPE_DNA:
        _check_alphabet(value, DNA_ALPHABET, "DNA", chain_ids)
    elif entity_type == TYPE_RNA:
        _check_alphabet(value, RNA_ALPHABET, "RNA", chain_ids)

    entity = {"type": entity_type, "chain_ids": chain_ids, "value": value}

    if msa is not None:
        if entity_type != TYPE_PROTEIN:
            raise BoltzPayloadError(
                "An MSA can only be given for a protein, but one was given for"
                " the %s of chain %s." % (entity_type.upper(), ", ".join(chain_ids)))
        entity["msa"] = msa

    return entity


def ligand_entity(chain_ids, value, ligand_type=None):
    """
    Builds a ligand entity. 'ligand_type' is TYPE_LIGAND_CCD or
    TYPE_LIGAND_SMILES; when it is None the value is guessed, see
    'guess_ligand_type'.
    """

    chain_ids = _check_chain_ids(chain_ids)
    text = (value or "").strip()
    if not text:
        raise BoltzPayloadError("Chain %s has no ligand: give a CCD code (such as"
                                " ATP) or a SMILES string."
                                % ", ".join(chain_ids))

    if ligand_type is None:
        ligand_type = guess_ligand_type(text)
    if ligand_type not in LIGAND_TYPES:
        raise BoltzPayloadError("Unknown ligand type %r." % ligand_type)

    # A CCD code is a short alphanumeric identifier; the service upper cases it
    # anyway, and doing it here keeps a lower case entry working. A SMILES string
    # is left exactly as it was typed, because its case carries meaning.
    if ligand_type == TYPE_LIGAND_CCD:
        text = text.upper()
        if not re.match(r"^[A-Z0-9]{1,5}$", text):
            raise BoltzPayloadError(
                "%r is not a valid CCD code for chain %s. A CCD code is 1 to 5"
                " letters or digits, such as ATP or SAH. Choose SMILES instead if"
                " this is a molecule string." % (value, ", ".join(chain_ids)))

    return {"type": ligand_type, "chain_ids": chain_ids, "value": text}


def _check_chain_ids(chain_ids):
    if isinstance(chain_ids, str):
        chain_ids = [chain_ids]
    chain_ids = [str(c).strip() for c in (chain_ids or []) if str(c).strip()]
    if not chain_ids:
        raise BoltzPayloadError("An entity was built without any chain identifier.")
    return chain_ids


def _check_alphabet(value, alphabet, label, chain_ids):
    bad = sorted(set(value) - alphabet)
    if bad:
        raise BoltzPayloadError(
            "The %s sequence of chain %s contains residues which are not allowed"
            " (%s). Only %s are accepted."
            % (label, ", ".join(chain_ids), ", ".join(bad),
               ", ".join(sorted(alphabet))))


###################################################################################################
# Binding.                                                                                         #
###################################################################################################

def ligand_protein_binding(binder_chain_id):
    """
    Asks for the binding metrics of a ligand against the proteins of the complex.
    'binder_chain_id' must be a chain of a ligand entity.
    """

    binder = (binder_chain_id or "").strip()
    if not binder:
        raise BoltzPayloadError("No binder chain was given for the affinity"
                                " prediction.")
    return {"binder_chain_id": binder, "type": BINDING_LIGAND_PROTEIN}


###################################################################################################
# The request.                                                                                     #
###################################################################################################

###################################################################################################
# Constraints: where a chain should sit.                                                          #
###################################################################################################

# What the service expects of a pocket, and what is sensible to ask of it.
MIN_POCKET_DISTANCE = 1.0
MAX_POCKET_DISTANCE = 20.0
DEFAULT_POCKET_DISTANCE = 6.0


def pocket_constraint(binder_chain_id, contact_residues, entities,
                      max_distance=DEFAULT_POCKET_DISTANCE, force=False):
    """
    A pocket: keep 'binder_chain_id' near these residues of the receptor.

        {"type": "pocket", "binder_chain_id": "B",
         "contact_residues": {"A": [10, 11, 12]},
         "max_distance_angstrom": 6.0, "force": false}

    'contact_residues' is {chain_id: [0-based indices]}, as everything about
    residues in this API is; turning the numbering a user reads into those indices
    is the protocol's job. 'force' false biases the prediction towards the pocket,
    true requires it.
    """

    chains = {}
    for entity in entities:
        for chain_id in entity.get("chain_ids", []):
            chains[chain_id] = entity

    if binder_chain_id not in chains:
        raise BoltzPayloadError(
            "The pocket names chain %s as the one to place, but the complex has"
            " no such chain." % binder_chain_id)

    if not contact_residues:
        raise BoltzPayloadError(
            "A pocket needs the residues which line it: none were given.")

    cleaned = {}
    for chain_id, indices in contact_residues.items():
        entity = chains.get(chain_id)
        if entity is None:
            raise BoltzPayloadError(
                "The pocket names chain %s, which the complex does not have."
                % chain_id)
        if chain_id == binder_chain_id:
            raise BoltzPayloadError(
                "Chain %s cannot line the pocket it is placed in." % binder_chain_id)

        length = len(entity.get("value", "")) if entity.get("type") in POLYMER_TYPES else 0
        residues = []
        for index in indices:
            if not isinstance(index, int) or isinstance(index, bool):
                raise BoltzPayloadError("The pocket of chain %s holds %r, which is"
                                        " not a residue index." % (chain_id, index))
            if index < 0 or (length and index >= length):
                raise BoltzPayloadError(
                    "Chain %s has %s residues, so it has no residue %s."
                    % (chain_id, length, index + 1))
            if index not in residues:
                residues.append(index)
        if residues:
            cleaned[chain_id] = sorted(residues)

    if not cleaned:
        raise BoltzPayloadError(
            "A pocket needs the residues which line it: none were given.")

    try:
        distance = float(max_distance)
    except (TypeError, ValueError):
        raise BoltzPayloadError("The distance of a pocket must be a number, but"
                                " %r was given." % (max_distance,))
    if not MIN_POCKET_DISTANCE <= distance <= MAX_POCKET_DISTANCE:
        raise BoltzPayloadError(
            "The distance of a pocket must be between %s and %s A, but %s was"
            " given." % (MIN_POCKET_DISTANCE, MAX_POCKET_DISTANCE, distance))

    return {"type": "pocket",
            "binder_chain_id": binder_chain_id,
            "contact_residues": cleaned,
            "max_distance_angstrom": distance,
            "force": bool(force)}


def build_payload(entities, num_samples=None, binding=None, templates=None,
                  constraints=None, model=MODEL):
    """
    Assembles and validates the whole request body.

    Raises BoltzPayloadError, with a message meant for the user, rather than
    letting the service reject the job after it has been submitted.
    """

    if not entities:
        raise BoltzPayloadError("Select at least one sequence to predict a"
                                " structure.")

    _check_unique_chain_ids(entities)
    check_msa_consistency(entities)

    payload = {"model": model, "input": {"entities": list(entities)}}

    if num_samples is not None:
        if not isinstance(num_samples, int) or isinstance(num_samples, bool):
            raise BoltzPayloadError("The number of samples must be a whole number.")
        if not MIN_SAMPLES <= num_samples <= MAX_SAMPLES:
            raise BoltzPayloadError(
                "The number of samples must be between %s and %s, but %s was"
                " requested." % (MIN_SAMPLES, MAX_SAMPLES, num_samples))
        payload["input"]["num_samples"] = num_samples

    if binding is not None:
        _check_binding(binding, entities)
        payload["input"]["binding"] = binding

    if templates:
        payload["input"]["templates"] = list(templates)
    if constraints:
        payload["input"]["constraints"] = list(constraints)

    return payload


def _check_unique_chain_ids(entities):
    seen = {}
    for index, entity in enumerate(entities):
        for chain_id in entity.get("chain_ids", []):
            if chain_id in seen:
                raise BoltzPayloadError(
                    "Chain %s is used twice, by entity %s and entity %s."
                    % (chain_id, seen[chain_id] + 1, index + 1))
            seen[chain_id] = index


def _check_binding(binding, entities):
    binder = binding.get("binder_chain_id")
    ligand_chains = set()
    protein_chains = set()
    for entity in entities:
        if entity.get("type") in LIGAND_TYPES:
            ligand_chains.update(entity.get("chain_ids", []))
        elif entity.get("type") == TYPE_PROTEIN:
            protein_chains.update(entity.get("chain_ids", []))

    if binder not in ligand_chains:
        raise BoltzPayloadError(
            "The affinity prediction needs a ligand as its binder, but chain %s"
            " is not a ligand. Add a ligand to the complex, or turn the affinity"
            " prediction off." % binder)
    if not protein_chains:
        raise BoltzPayloadError(
            "The affinity prediction needs at least one protein in the complex"
            " for the ligand to bind to.")


def check_msa_consistency(entities):
    """
    Enforces the two rules the service applies to the MSAs of a request:

      1. as soon as one protein uses a custom MSA, no protein may use the
         automatic one;
      2. proteins with the same sequence must share the same MSA setting.

    Both are checked here so that the user is told before a job is submitted.
    """

    proteins = [e for e in entities if e.get("type") == TYPE_PROTEIN]
    if not proteins:
        return

    def setting_of(entity):
        msa = entity.get("msa")
        if msa is None:
            return MSA_AUTOMATIC
        return msa.get("type", MSA_AUTOMATIC)

    settings = [setting_of(e) for e in proteins]

    if MSA_CUSTOM in settings and MSA_AUTOMATIC in settings:
        automatic = [", ".join(e["chain_ids"]) for e, s in zip(proteins, settings)
                     if s == MSA_AUTOMATIC]
        raise BoltzPayloadError(
            "The service does not allow automatic and custom alignments in the"
            " same prediction. Chain %s would use an automatically generated"
            " alignment while another chain uses a custom one: give an alignment"
            " for every protein chain, or for none of them."
            % "; ".join(automatic))

    # Rule 2: group the proteins by sequence and compare their settings.
    by_sequence = {}
    for entity, setting in zip(proteins, settings):
        by_sequence.setdefault(entity.get("value"), []).append((entity, setting))
    for sequence, group in by_sequence.items():
        distinct = set(s for _, s in group)
        if len(distinct) > 1:
            chains = "; ".join(", ".join(e["chain_ids"]) for e, _ in group)
            raise BoltzPayloadError(
                "Chains %s have the same sequence, so they must use the same kind"
                " of alignment, but %s were requested for them."
                % (chains, " and ".join(sorted(distinct))))


###################################################################################################
# Summary for the GUI.                                                                             #
###################################################################################################

def describe_payload(payload):
    """
    A short, readable description of a request, to show in a confirmation dialog
    next to the cost estimate. Never includes the base64 of an MSA.
    """

    entities = payload.get("input", {}).get("entities", [])
    lines = []
    for entity in entities:
        chains = ", ".join(entity.get("chain_ids", []))
        entity_type = entity.get("type", "?")
        value = entity.get("value", "")
        if entity_type in LIGAND_TYPES:
            kind = "CCD code" if entity_type == TYPE_LIGAND_CCD else "SMILES"
            detail = "%s %s" % (kind, value if len(value) <= 40 else value[:37] + "...")
        else:
            msa = entity.get("msa")
            if msa is None:
                msa_text = "automatic alignment"
            elif msa.get("type") == MSA_EMPTY:
                msa_text = "no alignment"
            else:
                msa_text = "custom alignment"
            detail = "%s residues, %s" % (len(value), msa_text)
        lines.append("  chain %s: %s (%s)" % (chains, entity_type, detail))

    for constraint in payload.get("input", {}).get("constraints", []):
        if constraint.get("type") != "pocket":
            continue
        # Counted from zero in the request, read from one by a user.
        where = []
        for chain_id in sorted(constraint.get("contact_residues", {})):
            numbers = [index + 1
                       for index in constraint["contact_residues"][chain_id]]
            shown = ", ".join(str(number) for number in numbers[:10])
            if len(numbers) > 10:
                shown += ", ... (%s residues)" % len(numbers)
            where.append("chain %s residue %s" % (chain_id, shown))
        lines.append("  chain %s is kept within %s A of %s%s"
                     % (constraint.get("binder_chain_id", "?"),
                        constraint.get("max_distance_angstrom", "?"),
                        "; ".join(where),
                        " (required)" if constraint.get("force") else ""))

    num_samples = payload.get("input", {}).get("num_samples")
    if num_samples:
        lines.append("  %s samples" % num_samples)
    binding = payload.get("input", {}).get("binding")
    if binding:
        lines.append("  affinity of the ligand in chain %s"
                     % binding.get("binder_chain_id"))

    return "\n".join(lines)
