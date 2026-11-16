# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Reads the results of a BoltzProt-1 design run.

The designs are NOT part of the job: a finished run reports them through a resource
of its own, one entry per design, fetched with 'BoltzClient.list_results'
(GET .../protein/design/{id}/results). Each entry carries

    id          the identifier of that design
    entities    every chain of the designed complex, the designed binder and the
                fixed chains of the target alike. The binder is the one whose
                chains are those the request asked to design.
    metrics     binding_confidence, structure_confidence, iptm,
                min_interaction_pae, helix_fraction, sheet_fraction, loop_fraction
    artifacts   structure and archive, as presigned URLs
    warnings    optional, each with a code and a message

**Unlike the structure prediction, this shape was NOT measured against a finished
run**: doing so means paying for a design run, which costs $0.25 at the very least.
It is taken from the type definitions of the official 'boltz-api' SDK, and the
endpoint was confirmed to exist. Everything here is therefore written defensively:
a field which is absent or of an unexpected type must not break the import.
"""


# The figures a design reports, in the order they are worth reading, with what they
# mean. 'binding_confidence' first: the documentation calls it the primary metric
# for hit discovery.
DESIGN_METRICS = [
    ("binding_confidence", "Binding confidence",
     "Confidence that the designed binder binds the target. The main figure for"
     " picking hits."),
    ("structure_confidence", "Structure confidence",
     "Confidence in the predicted structure of the design."),
    ("iptm", "ipTM",
     "Confidence in the interface between the binder and the target."),
    ("min_interaction_pae", "Interface PAE",
     "Smallest predicted aligned error at the interface, in Angstrom. Lower is"
     " better."),
    ("helix_fraction", "Helix fraction",
     "Fraction of the designed sequence forming alpha helices."),
    ("sheet_fraction", "Sheet fraction",
     "Fraction of the designed sequence forming beta sheets."),
    ("loop_fraction", "Loop fraction",
     "Fraction of the designed sequence in loops."),
]

# Figures where a small number is the good one.
LOWER_IS_BETTER = ("min_interaction_pae",)

# The figure designs are ranked by.
RANKING_METRIC = "binding_confidence"


class BoltzDesignResultError(ValueError):
    """
    Raised when the results of a design run can not be read. The message is meant
    to be shown to the user.
    """
    pass


def binder_chain_ids(job):
    """
    The chains the run was asked to design, read from the request the job carries.
    They identify which chain of a result is the designed binder rather than a copy
    of the target.
    """

    specification = ((job or {}).get("input") or {}).get("binder_specification") or {}
    chains = []
    for entity in specification.get("entities") or []:
        if (entity or {}).get("type") == "designed_protein":
            for chain_id in (entity.get("chain_ids") or []):
                if chain_id not in chains:
                    chains.append(chain_id)
    return chains


def results_of(response):
    """
    The list of designs out of what 'list_results' returned. The service pages its
    lists under 'data'; a bare list is accepted as well.
    """

    if isinstance(response, list):
        return response
    if isinstance(response, dict):
        data = response.get("data")
        if isinstance(data, list):
            return data
    raise BoltzDesignResultError(
        "The design run finished but its results could not be read.")


def designed_sequence_of(result, binder_chains=()):
    """
    The sequence of the designed binder of one result, or None when it can not be
    told apart from the rest of the complex.

    The binder is the entity holding the chains the request asked to design.

    When those chains are known, ONLY they are accepted: if no entity carries them,
    None is returned rather than falling back to another chain. Falling back would
    mean handing back a chain of the target as though it had been designed, which
    would be silently wrong. The fallback to "the only protein entity" is used just
    when the request could not be read at all, and never when there is more than
    one candidate.
    """

    entities = [e for e in ((result or {}).get("entities") or [])
                if isinstance(e, dict)]
    proteins = [e for e in entities if e.get("type") == "protein" and e.get("value")]

    if binder_chains:
        wanted = set(binder_chains)
        for entity in proteins:
            if set(entity.get("chain_ids") or []) & wanted:
                return entity.get("value")
        return None

    if len(proteins) == 1:
        return proteins[0].get("value")
    return None


def metrics_of(result):
    metrics = (result or {}).get("metrics")
    return metrics if isinstance(metrics, dict) else {}


def structure_url_of(result):
    """The URL of the structure of one design, or None when it carries none."""
    artifacts = (result or {}).get("artifacts") or {}
    structure = artifacts.get("structure") or {}
    return structure.get("url") or None


def warnings_of(result):
    """The warnings of one design, as readable lines."""
    lines = []
    for warning in ((result or {}).get("warnings") or []):
        if isinstance(warning, dict):
            message = warning.get("message") or warning.get("code")
            if message:
                lines.append(str(message))
        elif warning:
            lines.append(str(warning))
    return lines


def rank_results(response, job=None):
    """
    The designs, best first, each described by what PyMod needs to show and import.

    Returns a list of dictionaries with 'rank', 'id', 'sequence', 'metrics',
    'structure_url' and 'warnings'. Designs whose sequence can not be read are kept,
    with 'sequence' None, so that nothing disappears silently.
    """

    binder_chains = binder_chain_ids(job) if job else ()
    designs = []
    for result in results_of(response):
        designs.append({
            "id": (result or {}).get("id"),
            "sequence": designed_sequence_of(result, binder_chains),
            "metrics": metrics_of(result),
            "structure_url": structure_url_of(result),
            "warnings": warnings_of(result),
        })

    def sort_key(design):
        value = design["metrics"].get(RANKING_METRIC)
        try:
            return -float(value)
        except (TypeError, ValueError):
            # Designs with no figure to rank by go last, not first.
            return float("inf")

    designs.sort(key=sort_key)
    for position, design in enumerate(designs):
        design["rank"] = position + 1
    return designs


def format_value(name, value):
    if value is None:
        return "not reported"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if name in LOWER_IS_BETTER:
        return "%.2f A" % number
    return "%.3f" % number


def format_report(designs, run_name=None, top=5):
    """A readable summary of a design run, for the console."""

    lines = []
    if run_name:
        lines.append("BoltzProt-1 design run: %s" % run_name)
    usable = [d for d in designs if d["sequence"]]
    lines.append("%s designs were produced%s."
                 % (len(designs),
                    "" if len(usable) == len(designs)
                    else " (%s with a readable sequence)" % len(usable)))
    if not designs:
        return "\n".join(lines)

    lines.append("")
    lines.append("Best %s by %s:" % (min(top, len(designs)),
                                     dict((n, l) for n, l, _ in DESIGN_METRICS)[RANKING_METRIC].lower()))
    for design in designs[:top]:
        metrics = design["metrics"]
        sequence = design["sequence"] or ""
        lines.append("  %2d. %s  binding %s  structure %s  ipTM %s  %s residues"
                     % (design["rank"],
                        (design["id"] or "")[:22].ljust(22),
                        format_value("binding_confidence", metrics.get("binding_confidence")),
                        format_value("structure_confidence", metrics.get("structure_confidence")),
                        format_value("iptm", metrics.get("iptm")),
                        len(sequence) if sequence else "?"))
    warned = [d for d in designs if d["warnings"]]
    if warned:
        lines.append("")
        lines.append("%s designs came with a warning from the service." % len(warned))
    return "\n".join(lines)
