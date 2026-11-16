# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Reads the results of a finished Boltz prediction.

As 'boltz_payload' is the only module which knows the shape of a request, this is
the only one which knows the shape of a response, so that a change in the API costs
one file.

The shape below was measured by running a real prediction, since the free
'estimate-cost' endpoint returns no result:

    output.best_sample.metrics          the nine confidence figures
    output.best_sample.structure.url    the model, as mmCIF
    output.best_sample.ligand_structure the ligand alone, as SDF
    output.all_sample_results[]         the same, one entry per requested sample
    output.binding_metrics              binding_confidence, optimization_score
    output.archive.url                  a tar.gz with the metrics, the mmCIF and
                                        the predicted aligned error as .npz

Two things worth knowing about the files:

  - **the mmCIF already contains the ligand**, as HETATM records in its own chain,
    so the model is one complete file and the separate SDF is not needed to load
    the complex;
  - **the per-residue pLDDT is in the B-factor column** of that mmCIF (on the 0 to
    100 scale), which is where PyMod reads per-residue confidence from. The JSON
    carries no per-residue array: 'complex_plddt' is a single number for the whole
    complex.

Every result URL is presigned and **expires about thirty minutes after the job
completes**, while the data itself lives for a week. Results must therefore be
downloaded promptly; retrieving the job again yields fresh links.
"""

import io
import contextlib
import os


# The figures a prediction reports, with a readable name and what they mean, in
# the order they are worth reading. Keeping them here means the results window has
# no knowledge of the API.
METRICS = [
    ("structure_confidence", "Structure confidence",
     "Overall confidence in the predicted structure."),
    ("complex_plddt", "Complex pLDDT",
     "Per-residue confidence, averaged over the whole complex."),
    ("complex_iplddt", "Interface pLDDT",
     "The same, restricted to the residues at the interfaces between chains."),
    ("ptm", "pTM",
     "Predicted TM-score: how correct the overall fold is expected to be."),
    ("iptm", "ipTM",
     "The same for the interfaces between chains, so how correct their relative"
     " placement is expected to be."),
    ("protein_iptm", "Protein ipTM",
     "ipTM restricted to protein-protein interfaces."),
    ("ligand_iptm", "Ligand ipTM",
     "ipTM restricted to the interface with the ligand."),
    ("complex_pde", "Complex PDE",
     "Predicted distance error, in Angstrom. Lower is better."),
    ("complex_ipde", "Interface PDE",
     "The same at the interfaces. Lower is better."),
]

# The figures of a binding prediction, present only when one was asked for.
BINDING_METRICS = [
    ("binding_confidence", "Binding confidence",
     "Probability that the ligand binds at all."),
    ("optimization_score", "Optimization score",
     "Relative ranking of the strength of the binding, for lead optimisation."),
]

# Metrics where a small number is the good one.
LOWER_IS_BETTER = ("complex_pde", "complex_ipde")

# Metrics which only mean something when the complex has more than one chain: they
# all measure how correctly the chains are placed with respect to each other. For a
# single chain the service returns 0 for them, which reads as a terrible score when
# it actually means "there is no interface". They are reported as not applicable
# instead. 'complex_iplddt' is included because it comes back equal to
# 'complex_plddt' in that case, which is just as misleading.
INTERFACE_METRICS = ("iptm", "protein_iptm", "ligand_iptm", "complex_ipde",
                     "complex_iplddt")

# Shown in place of the value of a metric which does not apply.
NOT_APPLICABLE = "-"


class BoltzResultError(ValueError):
    """
    Raised when a finished job does not carry what is needed to import it. The
    message is meant to be shown to the user.
    """
    pass


###################################################################################################
# Reading the response.                                                                           #
###################################################################################################

def _output_of(job):
    output = (job or {}).get("output")
    if not isinstance(output, dict):
        raise BoltzResultError(
            "The prediction finished but returned no result. It may have been"
            " deleted on the server: predictions are kept for a limited time.")
    return output


def sample_results(job):
    """
    Every predicted sample, best one first. Each is a dictionary with 'metrics'
    and the files of that sample.
    """

    output = _output_of(job)
    samples = output.get("all_sample_results")
    if isinstance(samples, list) and samples:
        return samples
    best = output.get("best_sample")
    if isinstance(best, dict):
        return [best]
    raise BoltzResultError("The prediction returned no model.")


def best_sample(job):
    output = _output_of(job)
    best = output.get("best_sample")
    if isinstance(best, dict):
        return best
    return sample_results(job)[0]


def metrics_of(sample):
    """The confidence figures of one sample, as a plain dictionary."""
    metrics = (sample or {}).get("metrics")
    return metrics if isinstance(metrics, dict) else {}


def binding_metrics_of(job):
    """
    The figures of the binding prediction, or an empty dictionary when no affinity
    was asked for.
    """
    output = _output_of(job)
    metrics = output.get("binding_metrics")
    return metrics if isinstance(metrics, dict) else {}


def structure_url_of(sample):
    """The URL of the mmCIF of one sample. It contains the ligand as well."""
    structure = (sample or {}).get("structure")
    url = (structure or {}).get("url")
    if not url:
        raise BoltzResultError(
            "The prediction carries no structure file. If it finished a while ago,"
            " its download links may have expired: run it again, or retrieve it"
            " from the Boltz console.")
    return url


def archive_url_of(job):
    """
    The URL of the archive, which also holds the predicted aligned error. None
    when the job carries none.
    """
    output = _output_of(job)
    return ((output.get("archive") or {}).get("url")) or None


def url_expiry_of(sample):
    """When the links of a sample stop working, as the service reported it."""
    return ((sample or {}).get("structure") or {}).get("url_expires_at")


###################################################################################################
# Reporting.                                                                                      #
###################################################################################################

def chain_count(job):
    """
    How many chains the prediction was asked for, read back from the request the
    job carries. Returns 0 when it can not be told.
    """

    entities = ((job or {}).get("input") or {}).get("entities")
    if not isinstance(entities, list):
        return 0
    total = 0
    for entity in entities:
        chain_ids = (entity or {}).get("chain_ids")
        total += len(chain_ids) if isinstance(chain_ids, list) else 0
    return total


def has_interfaces(job):
    """
    True when the complex has more than one chain, and the interface figures
    therefore mean something. When the number of chains can not be told, the
    figures are shown rather than hidden.
    """
    count = chain_count(job)
    return count == 0 or count > 1


def applies(name, job):
    """Whether a metric means anything for this prediction."""
    if name in INTERFACE_METRICS and not has_interfaces(job):
        return False
    return True


def format_value(name, value, job=None):
    """
    Formats one figure. When 'job' is given, a figure which does not apply to that
    prediction is reported as such instead of as the zero the service returns.
    """

    if job is not None and not applies(name, job):
        return NOT_APPLICABLE
    if value is None:
        return "not reported"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if name in LOWER_IS_BETTER:
        return "%.2f A" % number
    return "%.3f" % number


def format_report(job, model_name=None, sample=None):
    """
    A readable report of a finished prediction, for the console and the results
    window. 'sample' selects which model is reported; the best one by default.
    """

    lines = []
    if model_name:
        lines.append("Boltz-2 prediction: %s" % model_name)
        lines.append("")

    samples = sample_results(job)
    metrics = metrics_of(sample if sample is not None else best_sample(job))
    lines.append("Confidence of the best model:")
    for name, label, _ in METRICS:
        if name in metrics:
            lines.append("  %-22s %s" % (label, format_value(name, metrics[name], job)))

    binding = binding_metrics_of(job)
    if binding:
        lines.append("")
        lines.append("Binding:")
        for name, label, _ in BINDING_METRICS:
            if name in binding:
                lines.append("  %-22s %s"
                             % (label, format_value(name, binding[name], job)))

    # Say why some figures are dashes, rather than leaving them to be read as a
    # failure.
    if not has_interfaces(job):
        lines.append("")
        lines.append("The figures marked '%s' measure how correctly the chains of a"
                     " complex are placed with respect to each other. This"
                     " prediction has a single chain, so they do not apply."
                     % NOT_APPLICABLE)

    if len(samples) > 1:
        lines.append("")
        lines.append("%s models were predicted, all of them loaded. The figures"
                     " above are those of one of them; each has its own, and its own"
                     " predicted aligned error, under Models -> Boltz-2."
                     % len(samples))

    return "\n".join(lines)


###################################################################################################
# Turning the model into something PyMod can load.                                                #
###################################################################################################

def pae_matrix_from_archive(archive_path, sample=0):
    """
    Reads the predicted aligned error out of the archive of a prediction, or
    returns None when it holds none.

    The archive is a tar.gz with 'prediction/sample_N_pae.npz' in it. The matrix is
    what the Models menu needs to draw a PAE plot, and it is not in the JSON of the
    job: only this file carries it.
    """

    import tarfile

    if not archive_path or not os.path.isfile(archive_path):
        return None

    try:
        import numpy
    except ImportError:
        return None

    try:
        with tarfile.open(archive_path) as archive:
            every = [n for n in archive.getnames() if n.endswith("_pae.npz")]
            names = [n for n in every if ("sample_%s_" % sample) in n]
            if not names:
                # Fall back to the only matrix in the archive ONLY when there is
                # exactly one and the first sample is the one being asked for. Any
                # looser fallback would hand the matrix of one model to another,
                # which is worse than having no plot: the axes would be labelled
                # with the residues of a different model.
                if len(every) == 1 and sample == 0:
                    names = every
                else:
                    return None
            handle = archive.extractfile(sorted(names)[0])
            if handle is None:
                return None
            with numpy.load(io.BytesIO(handle.read())) as data:
                for key in ("pae", "predicted_aligned_error", "arr_0"):
                    if key in data:
                        return data[key]
                keys = list(data.keys())
                if keys:
                    return data[keys[0]]
    except Exception as e:
        # A plot which can not be drawn is not worth failing an import for.
        print("- The predicted aligned error could not be read from the archive: %s."
              % e)
    return None


def load_model_in_pymod(pymod, pdb_path, model_name):
    """
    Loads a model as PyMod elements, one per polymer chain.

    Shared by the structure prediction and the protein design, which both end up
    with a predicted complex in a PDB file. The per-residue confidence the service
    put in the B-factor column comes along with the file, so PyMod can colour by it.
    A ligand does not become an element of its own: it stays in the PyMOL object as
    heteroatoms, as it does for any structure PyMod loads.
    """

    import pymod_lib.pymod_structure as pmstr
    import pymod_lib.pymod_vars as pmdt

    parsed = pmstr.Parsed_model_pdb_file(
        pymod, pdb_path,
        output_directory=pymod.structures_dirpath,
        new_file_name=model_name, model_root_name=model_name)

    colors = pmdt.pymol_light_colors_list
    elements = []
    for chain_number, element in enumerate(parsed.get_pymod_elements()):
        pymod.add_element_to_pymod(element,
                                   color=colors[chain_number % len(colors)])
        elements.append(element)
    return elements


@contextlib.contextmanager
def without_component_downloads():
    """
    Stops PyMOL fetching a definition of every residue while it reads an mmCIF.

    The mmCIF of a prediction names its components without describing them, so
    PyMOL asks the PDB for each one -- twenty-one files for a small protein,
    written into whatever directory PyMOL happens to have been started in, and
    printed as twenty-one lines of "Downloading .../SER.cif" for a model made of
    ordinary amino acids.

    They are fetched for the bonds, and the bonds are of no use here: the object
    exists to be saved as a PDB and deleted. Measured on a prediction of two
    chains, the conversion takes 7.18 seconds with them and 0.02 without, and the
    PDB which comes out is the same file to the byte.

    The downloader is PyMOL's own and is documented as internal, so it is put
    back as it was whatever happens.
    """

    try:
        import pymol as pymol_module
        from pymol import cmd as pymol_cmd
        real = pymol_module.internal.download_chem_comp
    except Exception:
        yield
        return

    def no_download(resn, quiet=1, _self=None):
        return ""

    pymol_module.internal.download_chem_comp = no_download
    try:
        pymol_cmd.download_chem_comp = no_download
    except Exception:
        pass
    try:
        yield
    finally:
        pymol_module.internal.download_chem_comp = real
        try:
            pymol_cmd.download_chem_comp = real
        except Exception:
            pass


def cif_to_pdb(cif_path, pdb_path):
    """
    Converts the mmCIF of a prediction to PDB, which is the format the structure
    parser of PyMod reads.

    PyMOL does the conversion, so nothing has to be installed. It keeps the
    B-factor column, which carries the per-residue pLDDT, and the heteroatoms of
    the ligand.

    Returns 'pdb_path'.
    """

    from pymol import cmd

    if not os.path.isfile(cif_path):
        raise BoltzResultError("The downloaded model file is missing: %s" % cif_path)

    # A name of its own, so an object already loaded in PyMOL is not disturbed.
    temporary_object = "pymod_boltz_tmp_model"
    try:
        cmd.delete(temporary_object)
    except Exception:
        pass

    try:
        with without_component_downloads():
            cmd.load(cif_path, temporary_object)
        if cmd.count_atoms(temporary_object) == 0:
            raise BoltzResultError("The model returned by the service has no atoms.")
        cmd.save(pdb_path, temporary_object)
    finally:
        try:
            cmd.delete(temporary_object)
        except Exception:
            pass

    if not os.path.isfile(pdb_path):
        raise BoltzResultError("The model could not be converted to the PDB format.")
    return pdb_path
