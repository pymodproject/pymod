# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Builds the request of a BoltzProt-1 protein design run.

As for 'boltz_payload', this is the only module which knows the shape of the
request, so that a change in the API costs one file, and it holds no GUI and no
PyMod code.

BoltzProt-1 designs binders against a target. The request verified against the
live service (with the free 'estimate-cost' endpoint) is:

    {"num_proteins": 10..1000000,
     "target": {"type": "no_template",
                "entities": [{"type": "protein", "chain_ids": ["T"],
                              "value": "<sequence>"}, ...]},
     "binder_specification": {
         "type": "no_template",
         "modality": "peptide"|"antibody"|"nanobody"|"custom_protein",
         "entities": [{"type": "designed_protein", "chain_ids": ["B"],
                       "value": "<binder specification>"}]}}

The surface to bind can be named, and then the request carries it as well:

    "target": {"type": "no_template",
               "entities": [...],
               "epitope_residues": {"T": [44, 45, 46]}}

**The indices are 0-based positions in the sequence which was sent**, not the
residue numbers of a PDB file: the API counts residues of a chain from zero
("All residue indices in constraints are 0-indexed"), and the sequence sent for a
target is its gapless sequence. Turning what a user reads in PyMod -- the
numbering of the structure -- into these indices is the protocol's job, not this
module's; what arrives here is already counted from zero and is checked against
the length of the chain it belongs to.

**Only a part of the API is covered here.** The full schema has 55 types: a binder
may also start from an uploaded structure ('structure_template') with per-chain
residue cropping and design motifs, and both target and binder accept covalent
bonds and pocket constraints. What is built here is the template-free case: design
a binder against sequences already in PyMod, optionally against a named epitope.
The rest can be added without touching anything else, since the request is
assembled in one place.

**The binder specification is a small language**, not a sequence. A number is a
run of designed residues, 'min..max' a run of variable length, and literal
one-letter codes are residues kept fixed:

    "30"                    30 designed residues
    "10..20"                between 10 and 20 designed residues
    "MKTAYI5..10VKSHFSRQ"   fixed MKTAYI, 5 to 10 designed, fixed VKSHFSRQ

**Costs are not those of a structure prediction.** The service refuses fewer than
10 designs and charges about $0.025 for each, so the smallest run costs around
$0.25 and a thousand designs around $25. 'estimate_cost' on the client gives the
real figure before anything is submitted.
"""

import re


# Endpoint resources of the design pipelines.
DESIGN_RESOURCE = "protein/design"

# What kind of binder is being designed.
MODALITY_PEPTIDE = "peptide"
MODALITY_ANTIBODY = "antibody"
MODALITY_NANOBODY = "nanobody"
MODALITY_CUSTOM = "custom_protein"
MODALITIES = (MODALITY_PEPTIDE, MODALITY_ANTIBODY, MODALITY_NANOBODY,
              MODALITY_CUSTOM)

# Binders Boltz maintains itself. Chosen instead of describing a binder by hand:
# the service picks from template libraries it keeps up to date, so no length is
# given and no modality either, since a curated binder is a different kind of
# request altogether.
CURATED_NANOBODY = "boltz_nanobody"
CURATED_ANTIBODY = "boltz_antibody"
CURATED_BINDERS = (CURATED_NANOBODY, CURATED_ANTIBODY)

# Bounds the service enforces on the size of a run.
MIN_PROTEINS = 10
MAX_PROTEINS = 1000000

# Entity types of the target, the same ones the structure prediction uses.
TYPE_PROTEIN = "protein"
TYPE_DNA = "dna"
TYPE_RNA = "rna"
TYPE_LIGAND_CCD = "ligand_ccd"
TYPE_LIGAND_SMILES = "ligand_smiles"

# Chain identifiers used by default: the target chains are named from T, the
# binder is B, so a design request reads clearly in the results.
BINDER_CHAIN_ID = "B"
TARGET_CHAIN_IDS = "TUVWXYZ"


class BoltzDesignError(ValueError):
    """
    Raised when a design request can not be built. The message is meant to be
    shown to the user as it is.
    """
    pass


###################################################################################################
# The binder specification language.                                                              #
###################################################################################################

# A specification is a sequence of blocks: a run of fixed residues, a fixed
# length, or a length range.
_BLOCK = re.compile(r"([A-Z]+)|(\d+)\.\.(\d+)|(\d+)")


def parse_binder_specification(text):
    """
    Checks a binder specification and returns (min_length, max_length).

    Raises BoltzDesignError, naming what is wrong, rather than letting the service
    reject the run after it has been submitted.
    """

    spec = "".join((text or "").split()).upper()
    if not spec:
        raise BoltzDesignError(
            "Describe the binder to design: a number of residues (for example 30),"
            " a range (10..20), or a mix of fixed residues and designed ones"
            " (MKTAYI5..10VKSHFSRQ).")

    position = 0
    minimum = maximum = 0
    designed = 0
    while position < len(spec):
        match = _BLOCK.match(spec, position)
        if match is None or match.start() != position:
            raise BoltzDesignError(
                "The binder specification %r can not be read, at %r. Use one-letter"
                " residue codes for the fixed parts, a number for a designed run,"
                " and min..max for a run of variable length."
                % (text, spec[position:position + 8]))
        residues, range_min, range_max, fixed_length = match.groups()
        if residues:
            minimum += len(residues)
            maximum += len(residues)
        elif fixed_length is not None:
            length = int(fixed_length)
            if length < 1:
                raise BoltzDesignError("A designed run must be at least 1 residue"
                                       " long, but %r was given." % text)
            minimum += length
            maximum += length
            designed += length
        else:
            low, high = int(range_min), int(range_max)
            if low < 1:
                raise BoltzDesignError("A designed run must be at least 1 residue"
                                       " long, but %r was given." % text)
            if low > high:
                raise BoltzDesignError(
                    "In the binder specification %r the range %s..%s runs backwards."
                    % (text, low, high))
            minimum += low
            maximum += high
            designed += high
        position = match.end()

    if designed == 0:
        raise BoltzDesignError(
            "The binder specification %r has no designed residue: it is entirely"
            " fixed, so there would be nothing to design. Add a number, such as 30,"
            " for the part to be designed." % text)

    return minimum, maximum


###################################################################################################
# Entities.                                                                                       #
###################################################################################################

def target_protein(chain_ids, sequence):
    return _target_polymer(TYPE_PROTEIN, chain_ids, sequence)


def target_dna(chain_ids, sequence):
    return _target_polymer(TYPE_DNA, chain_ids, sequence)


def target_rna(chain_ids, sequence):
    return _target_polymer(TYPE_RNA, chain_ids, sequence)


def _target_polymer(entity_type, chain_ids, sequence):
    chain_ids = _check_chain_ids(chain_ids)
    value = "".join((sequence or "").split()).replace("-", "").replace(".", "").upper()
    if not value:
        raise BoltzDesignError("The target chain %s has an empty sequence."
                               % ", ".join(chain_ids))
    return {"type": entity_type, "chain_ids": chain_ids, "value": value}


def target_ligand(chain_ids, value, ligand_type=TYPE_LIGAND_CCD):
    chain_ids = _check_chain_ids(chain_ids)
    text = (value or "").strip()
    if not text:
        raise BoltzDesignError("The target ligand in chain %s is empty."
                               % ", ".join(chain_ids))
    if ligand_type not in (TYPE_LIGAND_CCD, TYPE_LIGAND_SMILES):
        raise BoltzDesignError("Unknown ligand type %r." % ligand_type)
    if ligand_type == TYPE_LIGAND_CCD:
        text = text.upper()
    return {"type": ligand_type, "chain_ids": chain_ids, "value": text}


def curated_binder_specification(binder):
    """
    A binder taken from the libraries Boltz maintains, for a nanobody or an
    antibody. Neither a length nor a modality is given: the service chooses from
    its own templates, which is what makes this the right way to ask for those two
    classes rather than describing a single chain by hand.
    """

    if binder not in CURATED_BINDERS:
        raise BoltzDesignError("Unknown curated binder %r. Use one of: %s."
                               % (binder, ", ".join(CURATED_BINDERS)))
    return {"type": "boltz_curated", "binder": binder}


def designed_binder_specification(specification, modality=MODALITY_PEPTIDE,
                                  chain_id=BINDER_CHAIN_ID):
    """A binder described by hand: its composition and what class it belongs to."""

    if modality not in MODALITIES:
        raise BoltzDesignError("Unknown kind of binder %r. Use one of: %s."
                               % (modality, ", ".join(MODALITIES)))
    return {"type": "no_template",
            "modality": modality,
            "entities": [designed_binder(specification, chain_id)]}


def designed_binder(specification, chain_id=BINDER_CHAIN_ID):
    """The entity describing what has to be designed."""
    parse_binder_specification(specification)   # validates, raising on nonsense
    return {"type": "designed_protein",
            "chain_ids": [chain_id],
            "value": "".join((specification or "").split()).upper()}


def _check_chain_ids(chain_ids):
    if isinstance(chain_ids, str):
        chain_ids = [chain_ids]
    chain_ids = [str(c).strip() for c in (chain_ids or []) if str(c).strip()]
    if not chain_ids:
        raise BoltzDesignError("An entity was built without any chain identifier.")
    return chain_ids


def target_chain_ids(count, preferred=()):
    """
    Chain identifiers for the target, kept away from the binder's.

    A chain the target already has in PyMOL is kept, so that a structure loaded as
    chain A stays chain A in the design and in the complex which comes back. A
    preferred identifier is dropped when it is the binder's, when it is already
    taken, or when it is not a plain chain name; those chains fall back to T, U, V
    and so on.
    """

    preferred = list(preferred) + [None] * max(0, count - len(preferred))
    used = set()
    chosen = []
    for index in range(count):
        wanted = preferred[index]
        if isinstance(wanted, str):
            wanted = wanted.strip()
        if (wanted and wanted != BINDER_CHAIN_ID and wanted not in used
                and len(wanted) <= 4 and wanted.isalnum()):
            chosen.append(wanted)
            used.add(wanted)
        else:
            chosen.append(None)

    # The chains which kept nothing get the fallback names, skipping any already in
    # use so that no two chains collide.
    fallback = list(TARGET_CHAIN_IDS) + ["T%d" % (i + 1) for i in range(count)]
    available = (name for name in fallback
                 if name not in used and name != BINDER_CHAIN_ID)
    for index, name in enumerate(chosen):
        if name is None:
            chosen[index] = next(available)
            used.add(chosen[index])
    return chosen


###################################################################################################
# The request.                                                                                    #
###################################################################################################

###################################################################################################
# The epitope: which residues of the target the binder should bind.                               #
###################################################################################################

def parse_residue_numbers(text, chain_ids=()):
    """
    Reads what a user typed into residue numbers, per chain.

    Accepted: "45-52, 60" for a single target, and "T:45-52, U:88" when there is
    more than one, which is also the form the button writes from a PyMOL
    selection. A range is inclusive. Returns {chain_id: [numbers]}, the numbers as
    they were typed -- these are still the numbers of PyMod, not indices.
    """

    chain_ids = list(chain_ids)
    text = (text or "").strip()
    if not text:
        return {}

    found = {}
    current = chain_ids[0] if len(chain_ids) == 1 else None

    for piece in re.split(r"[,;]+", text.replace("\n", ",")):
        piece = piece.strip()
        if not piece:
            continue

        if ":" in piece:
            chain, piece = piece.split(":", 1)
            chain, piece = chain.strip(), piece.strip()
            if chain_ids and chain not in chain_ids:
                raise BoltzDesignError(
                    "There is no target chain %s. The target has %s."
                    % (chain, ", ".join(chain_ids)))
            current = chain
            if not piece:
                continue

        if current is None:
            raise BoltzDesignError(
                "The target has more than one chain (%s), so each residue has to"
                " say which one it belongs to, as in '%s:%s'."
                % (", ".join(chain_ids), chain_ids[0] if chain_ids else "T", piece))

        for number in _residue_range(piece):
            if number not in found.setdefault(current, []):
                found[current].append(number)

    return {chain: sorted(numbers) for chain, numbers in found.items() if numbers}


def _residue_range(piece):
    """'45' is one residue, '45-52' is eight."""

    match = re.match(r"^(-?\d+)\s*(?:-|\.\.)\s*(-?\d+)$", piece)
    if match:
        first, last = int(match.group(1)), int(match.group(2))
        if last < first:
            raise BoltzDesignError("The range %s ends before it starts." % piece)
        if last - first > 5000:
            raise BoltzDesignError("The range %s is too long to be an epitope."
                                   % piece)
        return list(range(first, last + 1))

    if re.match(r"^-?\d+$", piece):
        return [int(piece)]

    raise BoltzDesignError(
        "'%s' is not a residue number or a range of them. Give them as '45',"
        " '45-52', or 'T:45-52' when the target has more than one chain." % piece)


def epitope_annotation(indices_per_chain, target_entities, what="epitope"):
    """
    The 'epitope_residues' or 'non_binding_residues' of a request: {chain_id:
    [0-based indices]}, checked against the chains of the target and against the
    length of each of them. 'what' only names it in the messages.
    """

    if not indices_per_chain:
        return {}

    lengths = {}
    for entity in target_entities:
        if entity.get("type") in (TYPE_LIGAND_CCD, TYPE_LIGAND_SMILES):
            continue
        for chain_id in entity.get("chain_ids", []):
            lengths[chain_id] = len(entity.get("value", ""))

    annotation = {}
    for chain_id, indices in indices_per_chain.items():
        if chain_id not in lengths:
            raise BoltzDesignError(
                "There is no target chain %s to put an %s on. The target"
                " has %s." % (chain_id, what,
                              ", ".join(sorted(lengths)) or "no chain"))

        cleaned = []
        for index in indices:
            if not isinstance(index, int) or isinstance(index, bool):
                raise BoltzDesignError("The %s of chain %s holds %r, which is"
                                       " not a residue index."
                                       % (what, chain_id, index))
            if index < 0 or index >= lengths[chain_id]:
                raise BoltzDesignError(
                    "Chain %s has %s residues, so it has no residue %s."
                    % (chain_id, lengths[chain_id], index + 1))
            if index not in cleaned:
                cleaned.append(index)

        if cleaned:
            annotation[chain_id] = sorted(cleaned)

    return annotation


def build_design_payload(target_entities, binder_specification=None, num_proteins=None,
                         modality=MODALITY_PEPTIDE, binder_chain_id=BINDER_CHAIN_ID,
                         curated_binder=None, epitope_residues=None,
                         non_binding_residues=None):
    """
    Assembles and validates a design request.

    The binder is either described by hand, with 'binder_specification' and
    'modality', or taken from the curated libraries of the service by passing
    'curated_binder'. In the second case no length and no modality apply, and both
    are ignored.

    'epitope_residues' is {chain_id: [0-based indices]}: the residues of the target
    the binder should bind. Left out, the service is free to bind anywhere.
    'non_binding_residues' is the opposite: residues it should keep away from --
    a membrane face, a crystal contact, the site of another partner.
    """

    if not target_entities:
        raise BoltzDesignError("Select at least one sequence to design a binder"
                               " against.")

    if curated_binder is None and modality not in MODALITIES:
        raise BoltzDesignError("Unknown kind of binder %r. Use one of: %s."
                               % (modality, ", ".join(MODALITIES)))

    if not isinstance(num_proteins, int) or isinstance(num_proteins, bool):
        raise BoltzDesignError("The number of designs must be a whole number.")
    if num_proteins < MIN_PROTEINS:
        raise BoltzDesignError(
            "The service does not design fewer than %s binders at a time, but %s"
            " were asked for." % (MIN_PROTEINS, num_proteins))
    if num_proteins > MAX_PROTEINS:
        raise BoltzDesignError("At most %s designs can be asked for, but %s were."
                               % (MAX_PROTEINS, num_proteins))

    if curated_binder is None:
        # Only a binder of our own takes a chain, so only then can it clash.
        _check_unique_chain_ids(target_entities, binder_chain_id)
        specification = designed_binder_specification(binder_specification, modality,
                                                      binder_chain_id)
    else:
        specification = curated_binder_specification(curated_binder)

    target = {"type": "no_template", "entities": list(target_entities)}
    epitope = epitope_annotation(epitope_residues, target_entities)
    if epitope:
        target["epitope_residues"] = epitope
    avoided = epitope_annotation(non_binding_residues, target_entities,
                                 what="avoided surface")
    if avoided:
        target["non_binding_residues"] = avoided

    both = set()
    for chain_id, indices in epitope.items():
        both.update((chain_id, index) for index in indices
                    if index in avoided.get(chain_id, ()))
    if both:
        chain_id, index = sorted(both)[0]
        raise BoltzDesignError(
            "Residue %s of chain %s is named as part of the epitope and as a"
            " residue to avoid at the same time." % (index + 1, chain_id))

    return {
        "num_proteins": num_proteins,
        "target": target,
        "binder_specification": specification,
    }


def _check_unique_chain_ids(target_entities, binder_chain_id):
    seen = set()
    for index, entity in enumerate(target_entities):
        for chain_id in entity.get("chain_ids", []):
            if chain_id in seen:
                raise BoltzDesignError("The target uses chain %s twice." % chain_id)
            seen.add(chain_id)
    if binder_chain_id in seen:
        raise BoltzDesignError(
            "The binder would use chain %s, which the target already uses."
            % binder_chain_id)


def describe_design_payload(payload):
    """A readable summary, to show next to the cost before a run is started."""

    lines = []
    for entity in payload.get("target", {}).get("entities", []):
        chains = ", ".join(entity.get("chain_ids", []))
        value = entity.get("value", "")
        if entity.get("type") in (TYPE_LIGAND_CCD, TYPE_LIGAND_SMILES):
            detail = value if len(value) <= 40 else value[:37] + "..."
        else:
            detail = "%s residues" % len(value)
        lines.append("  target chain %s: %s (%s)"
                     % (chains, entity.get("type", "?"), detail))

    binder_spec = payload.get("binder_specification", {})
    if binder_spec.get("type") == "boltz_curated":
        curated = binder_spec.get("binder", "?")
        lines.append("  binder: %s, from the libraries maintained by Boltz"
                     % curated.replace("boltz_", ""))
    else:
        entities = binder_spec.get("entities", [])
        if entities:
            specification = entities[0].get("value", "")
            minimum, maximum = parse_binder_specification(specification)
            length = ("%s residues" % minimum if minimum == maximum
                      else "%s to %s residues" % (minimum, maximum))
            lines.append("  binder: %s, %s (%s)"
                         % (binder_spec.get("modality", "?"), length, specification))

    epitope = payload.get("target", {}).get("epitope_residues") or {}
    if epitope:
        # Counted from zero in the request, read from one by a user.
        for chain_id in sorted(epitope):
            numbers = [index + 1 for index in epitope[chain_id]]
            shown = ", ".join(str(n) for n in numbers[:12])
            if len(numbers) > 12:
                shown += ", ... (%s residues)" % len(numbers)
            lines.append("  epitope on chain %s: residue %s of the sequence sent"
                         % (chain_id, shown))
    else:
        lines.append("  no epitope given: the binder may bind anywhere")

    avoided = payload.get("target", {}).get("non_binding_residues") or {}
    for chain_id in sorted(avoided):
        numbers = [index + 1 for index in avoided[chain_id]]
        shown = ", ".join(str(n) for n in numbers[:12])
        if len(numbers) > 12:
            shown += ", ... (%s residues)" % len(numbers)
        lines.append("  kept away from chain %s: residue %s of the sequence sent"
                     % (chain_id, shown))

    lines.append("  %s designs" % payload.get("num_proteins"))
    return "\n".join(lines)
