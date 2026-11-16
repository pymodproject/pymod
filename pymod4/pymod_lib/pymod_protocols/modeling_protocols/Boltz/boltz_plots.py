# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Plots of the quality of a Boltz-2 model, for the "Models" menu.

They use the same windows as the AlphaFold and the ESMFold ones, so they look and
behave the way a PyMod user already expects: the per-residue confidence in the plot
window of PyMod, which can highlight residues in PyMOL when clicked, and the
predicted aligned error in the contact map window.

A Boltz session, registered by the prediction protocol, carries what is needed:
'model_name', 'full_output_path' (the PDB of the model, whose B-factor column holds
the per-residue pLDDT), 'new_model_elements' (its chains, in order) and 'pae_matrix'
when the predicted aligned error could be read out of the results archive.

**One difference from the other services.** The predicted aligned error of Boltz
covers *tokens*, not residues: a ligand contributes one per atom, so a complex of 39
residues with an ATP comes back as a 70 x 70 matrix. Labelling those extra rows with
residues would be wrong, so the matrix is cropped to the residues of the model and
what was dropped is reported.
"""

import numpy as np

from pymod_lib.pymod_plot.pymod_plot_qt import PyMod_plot_window_qt
from pymod_lib.pymod_protocols.structural_analysis_protocols.contact_map_analysis import (
    Contact_map_analysis_window_qt)


def _click_handler():
    """
    The handler which selects a residue in PyMOL when a point of the profile is
    clicked, reused from the AlphaFold plots rather than copied, so that the plots of
    the three services behave identically.

    Imported here and not at the top of the module: the AlphaFold module reaches
    'pymod_main', which reaches the menu commands, which import the AlphaFold module
    again. Importing it while this module is being read therefore fails with a
    circular import, depending only on which module happened to be read first.
    """

    try:
        from pymod_lib.pymod_protocols.modeling_protocols.AlphaFold import (
            _highlight_in_pymol_from_plddt_plot_qt)
        return _highlight_in_pymol_from_plddt_plot_qt
    except ImportError as e:
        # Reached only if this module is used before the plugin is fully loaded,
        # which does not happen from the menu. The plot is still worth drawing
        # without the interaction, so it is not allowed to fail.
        print("- The plot will not be interactive: %s." % e)
        return None


def model_elements_of(session):
    """The chains of a model, in order, however the session recorded them."""

    elements = getattr(session, "new_model_elements", None)
    if not elements:
        single = getattr(session, "new_model_element", None)
        elements = [single] if single is not None else []
    return [element for element in elements if element is not None]


def chain_profiles(session):
    """
    The per-residue confidence of every chain, read from the B-factor column of the
    model in PyMOL, which is where the service put it.

    Returns a list of (element, residues, values) in chain order, skipping chains
    whose values could not be read.
    """

    from pymol import cmd

    profiles = []
    for element in model_elements_of(session):
        selector = element.get_pymol_selector()
        values = []
        try:
            cmd.iterate("%s and name CA" % selector, "values.append(b)",
                        space={"values": values})
        except Exception:
            values = []
        if not values:
            continue
        residues = element.get_polymer_residues()[:len(values)]
        profiles.append((element, residues, values[:len(residues)]))
    return profiles


def show_boltz_plots(pymod, session, plot_type):
    """
    Shows the per-residue confidence or the predicted aligned error of one Boltz
    session. Called from the Models menu.
    """

    if plot_type == "plddt":
        return _show_plddt_plot(pymod, session)
    if plot_type == "pae":
        return _show_pae_plot(pymod, session)
    raise ValueError("Unknown kind of plot: %r." % plot_type)


def _show_plddt_plot(pymod, session):
    """
    The per-residue confidence of the model, one curve per chain, on the 0 to 100
    scale of the file. The confidence is also stored on the residues, so the model
    can be coloured by it as an AlphaFold one can.
    """

    profiles = chain_profiles(session)
    if not profiles:
        pymod.main_window.show_error_message(
            "pLDDT not available",
            "The per-residue confidence of this model could not be read. It comes"
            " from the B-factor column of the model in PyMOL, so the model has to"
            " still be loaded.")
        return None

    model_name = getattr(session, "model_name", "the model")

    window = PyMod_plot_window_qt(pymod.main_window)
    window.initialize(pymod=pymod,
                      title="Predicted Local Distance Difference Test (pLDDT) of %s"
                            % model_name)
    window.build_plotting_area(
        messagebar_initial_text="Set the 'Interact on click' option to 'Yes' to"
                                " click on the plot and to highlight residues in"
                                " PyMOL.",
        update_messagebar=True,
        messagebar_text_on_update="Selected: __residue_name__"
                                  " __residue_pdb_position__ of __plot_name__,"
                                  " pLDDT value: __y__",
        # Without this the profile is drawn but does nothing when clicked: this is
        # what connects it to PyMOL, selecting and centring the residue of the point
        # which was clicked.
        on_click_action=_click_handler(),
        x_label_text="Positions",
        y_label_text="pLDDT")

    for element, residues, values in profiles:
        selector = element.get_pymol_selector()
        positions = [residue.db_index for residue in residues]
        additional_data = [
            {"residue_name": residue.three_letter_code,
             "residue_pdb_position": residue.db_index,
             "pymol_selector": "%s and resi %s" % (selector, residue.db_index),
             "export_label": "%s %s" % (residue.three_letter_code,
                                        residue.db_index)}
            for residue in residues]

        window.add_plot(positions, values,
                        label=element.compact_header,
                        additional_data=additional_data)

        # Store the confidence on the residues, as the AlphaFold and ESMFold
        # sessions do, so the model can also be coloured by it.
        _store_plddt_on_residues(element, residues, values)

    # The window has an "Interact on click" switch which starts on "No", so a fresh
    # plot does nothing when clicked until it is turned on. Since the point of this
    # plot is to find a doubtful residue and look at it in PyMOL, it starts on "Yes"
    # instead; the switch is still there to turn it off.
    _enable_interaction(window)

    window.show()
    return window


def _enable_interaction(window):
    """
    Turns the "Interact on click" switch of a plot window on, so clicking a point
    selects the residue in PyMOL straight away.
    """

    try:
        window.interact_mode = True
        if hasattr(window, "interact_button"):
            window.interact_button.setChecked(True)
    except Exception:
        # An interactive plot is the nicety here, not the plot itself.
        pass


def _store_plddt_on_residues(element, residues, values):
    """
    Puts the confidence of each residue on the residue itself, in the shape PyMod
    colours by: a score and the interval it falls in.
    """

    try:
        array = np.array(values, dtype=float)
        bins = np.linspace(array.min(), array.max(), num=10)
        intervals = np.digitize(array, bins)
        for residue, score, interval in zip(residues, array, intervals):
            residue.plddt_score = {"score": float(score), "interval": int(interval)}
        element._has_plddt_scores = True
    except Exception:
        # Colouring by confidence is a convenience: it must not stop the plot.
        pass


def _show_pae_plot(pymod, session):
    """
    The predicted aligned error of the model, in the contact map window, as the
    ESMFold sessions show theirs.
    """

    matrix = getattr(session, "pae_matrix", None)
    if matrix is None:
        pymod.main_window.show_error_message(
            "Predicted aligned error not available",
            "This model carries no predicted aligned error. It is read from the"
            " results archive of the prediction, so the archive has to have been"
            " downloaded with the model.")
        return None

    elements = model_elements_of(session)
    residues = []
    selectors = []
    for element in elements:
        selector = element.get_pymol_selector()
        for residue in element.get_polymer_residues():
            residues.append(residue)
            selectors.append("object %s and n. CA and i. %s"
                             % (selector, residue.db_index))

    matrix = [list(row) for row in matrix]
    if len(matrix) > len(residues) and residues:
        print("- The predicted aligned error of %s covers %s tokens while the model"
              " has %s residues: the extra ones are the atoms of its ligands, and"
              " the matrix is cropped to the residues."
              % (getattr(session, "model_name", "the model"), len(matrix),
                 len(residues)))

    # Keep the matrix and the reference lists the same length: the grid is sized
    # from the matrix, so a mismatch would mislabel it.
    size = min(len(matrix), len(residues)) if residues else len(matrix)
    if not size:
        pymod.main_window.show_error_message(
            "Predicted aligned error not available",
            "The residues of this model could not be matched to the predicted"
            " aligned error.")
        return None
    matrix = [row[:size] for row in matrix[:size]]
    residues = residues[:size]
    selectors = selectors[:size]

    window = Contact_map_analysis_window_qt(pymod.main_window)
    window.initialize_map(pymod=pymod,
                          data_array=matrix,
                          pymod_elements=elements,
                          ref_residues=residues,
                          ref_selectors=selectors,
                          title="Predicted Aligned Error (PAE) for %s"
                                % getattr(session, "model_name", "the model"),
                          pixel_size=5,
                          feature_type="pae",
                          threshold=35,
                          interaction_center="ca")
    window.show()
    return window
