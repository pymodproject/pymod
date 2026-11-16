# Copyright 2026 by the PyMod developers. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Local, in-app ESPript-style rendering of a sequence alignment.

This is an offline alternative to the web-based ESPript protocol (espript.py).
The alignment is drawn entirely with PyQt5 (the same toolkit PyMod's GUI already
uses, so NO new dependency is introduced) on a QGraphicsScene, with identity /
similarity shading and an optional secondary-structure topology row taken from a
selected 3D structure. The figure can be exported to PNG or PDF via QPainter.
"""

import os

from pymol.Qt import QtWidgets, QtCore, QtGui

from pymod_lib import pymod_vars
from pymod_lib.pymod_gui.shared_gui_components_qt import (PyMod_protocol_window_qt,
                                                          PyMod_combobox_qt,
                                                          PyMod_entryfield_qt)

from ._evolutionary_analysis_base import Evolutionary_analysis_protocol


# Per-residue conservation analyses whose results can be mapped onto the figure.
# Each entry: (menu label, element "_has_..." flag, residue score attribute,
# PyMod color dict, PyMod color-name prefix). The color dicts/prefixes are reused
# from pymod_vars so the figure matches the main-window coloring exactly.
SCORE_ANALYSES = [
    ("CAMPO conservation scores", "_has_campo_scores", "campo_score",
     pymod_vars.campo_color_dict, pymod_vars.pymol_campo_color_name),
    ("Entropy scores", "_has_entropy_scores", "entropy_score",
     pymod_vars.entropy_color_dict, pymod_vars.pymol_entropy_color_name),
    ("Pair conservation scores", "_has_pc_scores", "pc_score",
     pymod_vars.pc_color_dict, pymod_vars.pymol_pc_color_name),
]


###################################################################################################
# Conservation / color logic (pure, no GUI).                                                     #
###################################################################################################

# Physicochemical similarity groups (ESPript/Clustal-like). A column counts as
# "similar" when every non-gap residue in it belongs to a single group.
SIMILARITY_GROUPS = [
    set("AVLIM"),   # aliphatic / non-polar
    set("FYW"),     # aromatic
    set("ST"),      # small hydroxyl
    set("DE"),      # acidic
    set("NQ"),      # amide
    set("KRH"),     # basic
    set("G"),
    set("P"),
    set("C"),
]

# Colors for the residue background boxes and SS glyphs.
IDENTITY_COLOR = "#E03030"    # red box, white bold letter (fully conserved)
SIMILARITY_COLOR = "#FBE38A"  # pale yellow box (same-group conservation)
HELIX_COLOR = "#E03030"
STRAND_COLOR = "#3F8F3F"
RULER_COLOR = "#555555"

GAP_CHARS = ("-", ".")


def column_conservation(sequences):
    """
    Returns two lists (one entry per alignment column): 'identity' is True when
    every non-gap residue in the column is the same letter; 'similarity' is True
    when every non-gap residue falls in a single SIMILARITY_GROUPS set. Columns
    with a single non-gap residue are not counted as conserved.
    """
    length = len(sequences[0]) if sequences else 0
    identity = [False] * length
    similarity = [False] * length
    for col in range(length):
        residues = [seq[col].upper() for seq in sequences
                    if seq[col] not in GAP_CHARS]
        if len(residues) < 2:
            continue
        residue_set = set(residues)
        if len(residue_set) == 1:
            identity[col] = True
            similarity[col] = True
            continue
        for group in SIMILARITY_GROUPS:
            if residue_set <= group:
                similarity[col] = True
                break
    return identity, similarity


def secondary_structure_by_column(structure_element):
    """
    Maps a reference element's per-residue secondary structure onto the aligned
    columns. Returns a list (one entry per column) of normalized SS codes:
    'H' (helix), 'E' (strand) or '' (loop/coil/gap).
    """
    residues = structure_element.get_polymer_residues()
    col_ss = []
    res_index = 0
    for char in structure_element.my_sequence:
        if char in GAP_CHARS:
            col_ss.append("")
            continue
        ss = ""
        if res_index < len(residues):
            ss = (residues[res_index].secondary_structure or "")
        res_index += 1
        ss = ss.upper()
        if ss in ("H", "G", "I"):
            col_ss.append("H")
        elif ss in ("E", "S", "B"):
            col_ss.append("E")
        else:
            col_ss.append("")
    return col_ss


def contiguous_runs(ss_list, start, end):
    """
    Yields (run_type, run_start, run_end_inclusive) for the H/E stretches inside
    the [start, end) column slice, used to draw helix/strand glyphs spanning a
    whole element rather than one box per residue.
    """
    run_type = None
    run_start = None
    for col in range(start, end):
        kind = ss_list[col] if col < len(ss_list) else ""
        if kind != run_type:
            if run_type in ("H", "E"):
                yield run_type, run_start, col - 1
            run_type = kind
            run_start = col
    if run_type in ("H", "E"):
        yield run_type, run_start, end - 1


def scores_by_column(children, attr):
    """
    Reconstructs a per-column score array from the per-residue scores stored by a
    conservation analysis (CAMPO/Entropy/Pair Conservation). These analyses assign
    the same per-column value to every non-gap residue in a column, so reading the
    first non-gap residue that carries the score for each column is sufficient.
    Returns a list (one entry per column) of score values, or None where no score
    is available (gap-only / tossed columns).
    """
    length = len(children[0].my_sequence) if children else 0
    col_scores = [None] * length
    for child in children:
        residues = child.get_polymer_residues()
        res_index = 0
        for col, char in enumerate(child.my_sequence):
            if char in GAP_CHARS:
                continue
            if res_index < len(residues) and col_scores[col] is None:
                value = getattr(residues[res_index], attr, None)
                if value is not None:
                    col_scores[col] = value
            res_index += 1
    return col_scores


def _score_bin(value):
    """
    Extracts the integer bin index from a stored analysis score. CAMPO and Entropy
    store a dict like {"...-score": float, "interval": int} (and the color key is
    the 'interval'); Pair Conservation stores the bin int directly. Returns None
    when no bin is available (gap-only / tossed columns).
    """
    if isinstance(value, dict):
        value = value.get("interval")
    return value


def score_to_qcolor(value, color_dict, prefix):
    """
    Maps a stored analysis score to a QColor using a PyMod color dict (whose keys
    look like '<prefix>_<n>' and '<prefix>_None', with 0..1 RGB tuples).
    """
    bin_index = _score_bin(value)
    rgb = None
    if bin_index is not None:
        rgb = color_dict.get("%s_%d" % (prefix, int(bin_index)))
    if rgb is None:
        rgb = color_dict.get("%s_None" % prefix)
    if rgb is None:
        return None
    return QtGui.QColor.fromRgbF(*rgb)


def make_column_colors(children, attr, color_dict, prefix):
    """
    Builds a per-column list of QColor (or None) from a conservation analysis,
    ready to feed to build_alignment_scene as 'column_colors'.
    """
    return [score_to_qcolor(v, color_dict, prefix)
            for v in scores_by_column(children, attr)]


def _text_color_for(qcolor):
    """Black or white text, whichever is more legible on the given background."""
    luminance = (0.299 * qcolor.redF() + 0.587 * qcolor.greenF()
                 + 0.114 * qcolor.blueF())
    return QtGui.QColor("black") if luminance >= 0.5 else QtGui.QColor("white")


###################################################################################################
# Qt scene drawing.                                                                              #
###################################################################################################

def build_alignment_scene(headers, sequences, ss_list=None,
                          color_scheme="similarity", block_length=60,
                          column_colors=None):
    """
    Builds and returns a QGraphicsScene with an ESPript-style alignment, wrapped
    into horizontal blocks of 'block_length' columns. 'ss_list' is an optional
    per-column SS string (from secondary_structure_by_column).

    'column_colors', when given, is a per-column list of QColor (or None) used to
    shade the residue boxes (e.g. CAMPO/Entropy/Pair-Conservation scores). When it
    is None, the boxes are shaded by identity/similarity instead.
    """
    scene = QtWidgets.QGraphicsScene()
    n_seq = len(sequences)
    length = len(sequences[0]) if sequences else 0
    if length == 0:
        return scene

    identity, similarity = column_conservation(sequences)
    has_ss = bool(ss_list) and any(ss_list)

    # Fonts and cell geometry, derived from the monospaced font metrics so boxes
    # line up exactly under the letters regardless of the platform font.
    seq_font = QtGui.QFont("Courier New")
    seq_font.setStyleHint(QtGui.QFont.Monospace)
    seq_font.setPointSize(12)
    bold_font = QtGui.QFont(seq_font)
    bold_font.setBold(True)
    label_font = QtGui.QFont(seq_font)
    ruler_font = QtGui.QFont(seq_font)
    ruler_font.setPointSize(8)

    fm = QtGui.QFontMetricsF(seq_font)
    cell_w = fm.horizontalAdvance("W") + 4.0
    cell_h = fm.height() + 2.0

    label_fm = QtGui.QFontMetricsF(label_font)
    label_width = max((label_fm.horizontalAdvance(h) for h in headers), default=0.0)
    left_margin = label_width + 12.0

    ss_rows = 1 if has_ss else 0
    rows_per_block = ss_rows + 1 + n_seq + 1  # ss + ruler + sequences + spacer
    n_blocks = (length + block_length - 1) // block_length

    black = QtGui.QColor("black")
    white = QtGui.QColor("white")
    no_pen = QtGui.QPen(QtCore.Qt.NoPen)

    def add_centered_text(text, cx, cy, font, color):
        item = scene.addSimpleText(text, font)
        item.setBrush(QtGui.QBrush(color))
        br = item.boundingRect()
        item.setPos(cx - br.width() / 2.0, cy - br.height() / 2.0)
        return item

    for block in range(n_blocks):
        col_start = block * block_length
        col_end = min(col_start + block_length, length)
        block_top = block * rows_per_block * cell_h

        # --- Secondary-structure topology row. ---
        if has_ss:
            ss_y = block_top
            for kind, run_start, run_end in contiguous_runs(ss_list, col_start, col_end):
                x = left_margin + (run_start - col_start) * cell_w
                run_w = (run_end - run_start + 1) * cell_w
                mid_y = ss_y + cell_h / 2.0
                if kind == "H":
                    path = QtGui.QPainterPath()
                    path.addRoundedRect(QtCore.QRectF(x, ss_y + cell_h * 0.18,
                                                      run_w - 2.0, cell_h * 0.64),
                                        cell_h * 0.32, cell_h * 0.32)
                    scene.addPath(path, QtGui.QPen(black, 0.8),
                                  QtGui.QBrush(QtGui.QColor(HELIX_COLOR)))
                else:  # strand -> arrow
                    head = min(cell_w * 1.4, run_w * 0.5)
                    bx1 = x + run_w - head - 2.0
                    body = cell_h * 0.16
                    headh = cell_h * 0.34
                    poly = QtGui.QPolygonF([
                        QtCore.QPointF(x, mid_y - body),
                        QtCore.QPointF(bx1, mid_y - body),
                        QtCore.QPointF(bx1, mid_y - headh),
                        QtCore.QPointF(x + run_w - 2.0, mid_y),
                        QtCore.QPointF(bx1, mid_y + headh),
                        QtCore.QPointF(bx1, mid_y + body),
                        QtCore.QPointF(x, mid_y + body),
                    ])
                    scene.addPolygon(poly, QtGui.QPen(black, 0.6),
                                     QtGui.QBrush(QtGui.QColor(STRAND_COLOR)))

        # --- Column-number ruler (every 10 columns). ---
        ruler_y = block_top + ss_rows * cell_h
        for col in range(col_start, col_end):
            human_col = col + 1
            if human_col % 10 == 0:
                cx = left_margin + (col - col_start + 0.5) * cell_w
                add_centered_text(str(human_col), cx, ruler_y + cell_h / 2.0,
                                  ruler_font, QtGui.QColor(RULER_COLOR))

        # --- Sequence rows. ---
        for row, (header, sequence) in enumerate(zip(headers, sequences)):
            y = block_top + (ss_rows + 1 + row) * cell_h
            cy = y + cell_h / 2.0
            # Header, right-aligned against the residue grid.
            label = scene.addSimpleText(header, label_font)
            label.setBrush(QtGui.QBrush(black))
            lbr = label.boundingRect()
            label.setPos(left_margin - 6.0 - lbr.width(), cy - lbr.height() / 2.0)

            for col in range(col_start, col_end):
                x = left_margin + (col - col_start) * cell_w
                char = sequence[col]
                fill = None
                txt_color = black
                font = seq_font
                if char not in GAP_CHARS:
                    if column_colors is not None:
                        # Score-based shading (CAMPO / Entropy / Pair Conservation).
                        fill = column_colors[col]
                        if fill is not None:
                            txt_color = _text_color_for(fill)
                    elif identity[col]:
                        fill = QtGui.QColor(IDENTITY_COLOR)
                        txt_color = white
                        font = bold_font
                    elif color_scheme == "similarity" and similarity[col]:
                        fill = QtGui.QColor(SIMILARITY_COLOR)
                if fill is not None:
                    rect = scene.addRect(QtCore.QRectF(x, y, cell_w, cell_h),
                                         no_pen, QtGui.QBrush(fill))
                    rect.setZValue(-1)
                add_centered_text(char, x + cell_w / 2.0, cy, font, txt_color)

    scene.setSceneRect(scene.itemsBoundingRect().adjusted(-10, -10, 10, 10))
    return scene


def render_scene_to_file(scene, path):
    """
    Exports a QGraphicsScene to PNG or PDF (chosen by file extension), using
    QPainter only (no extra dependencies). Returns True on success.
    """
    rect = scene.sceneRect()
    ext = os.path.splitext(path)[1].lower()

    if ext == ".pdf":
        writer = QtGui.QPdfWriter(path)
        writer.setResolution(300)
        # Page size in points (1/72 inch); scene units are pixels at ~96 dpi.
        pt_w = rect.width() * 72.0 / 96.0
        pt_h = rect.height() * 72.0 / 96.0
        writer.setPageSizeMM(QtCore.QSizeF(pt_w * 25.4 / 72.0, pt_h * 25.4 / 72.0))
        painter = QtGui.QPainter(writer)
        scene.render(painter)
        painter.end()
        return True

    # Default: raster image (PNG, also JPG/BMP if requested).
    scale = 2.0  # supersample for crisp output
    image = QtGui.QImage(int(rect.width() * scale), int(rect.height() * scale),
                         QtGui.QImage.Format_ARGB32)
    image.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(image)
    painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
    painter.setRenderHint(QtGui.QPainter.TextAntialiasing, True)
    scene.render(painter, QtCore.QRectF(image.rect()), rect)
    painter.end()
    return image.save(path)


###################################################################################################
# Protocol.                                                                                       #
###################################################################################################

class Alignment_figure_analysis(Evolutionary_analysis_protocol):
    """
    Renders the alignment locally (PyQt5, in a Qt window) as an offline
    alternative to the web-based ESPript service.
    """

    def launch_from_gui(self):
        children = self.input_cluster_element.get_children()

        # Structures whose secondary structure can be drawn on top of the alignment.
        self.ss_structures_list = ["None"]
        self.ss_structure_elements_list = [None]
        for element in children:
            if element.has_structure() and element.polymer_type == "protein":
                self.ss_structures_list.append(element.my_header)
                self.ss_structure_elements_list.append(element)

        # Available "Color by" modes. The first two are always available; the
        # conservation-analysis modes are offered only when their scores have
        # already been computed (via Alignments -> Evolutionary conservation).
        # Each entry: (label, kind, attr, color_dict, prefix).
        self.color_modes = [
            ("Identity only", "identity", None, None, None),
            ("Identity + similarity", "similarity", None, None, None),
        ]
        for label, flag, attr, color_dict, prefix in SCORE_ANALYSES:
            if any(getattr(child, flag, False) for child in children):
                self.color_modes.append((label, "score", attr, color_dict, prefix))
        self.color_mode_labels = [mode[0] for mode in self.color_modes]

        self.options_window = Alignment_figure_options_window_qt(
            self.pymod.main_window,
            protocol=self,
            title="Render Alignment (local)",
            upper_frame_title="Options for the local alignment figure",
            submit_command=self.render_state)
        self.options_window.show()


    def render_state(self):
        # Read the options.
        ss_index = self.options_window.ss_combobox.get_index()
        ss_element = self.ss_structure_elements_list[ss_index]
        color_mode = self.color_modes[self.options_window.colorby_combobox.get_index()]
        try:
            block_length = int(self.options_window.block_entry.getvalue(validate=False))
        except (ValueError, TypeError):
            block_length = 60
        block_length = max(10, min(200, block_length))

        self.options_window.destroy()

        # Extract the alignment data.
        children = list(self.input_cluster_element.get_children())
        if not children:
            self.pymod.main_window.show_error_message(
                "Empty Alignment", "There are no sequences to render.")
            return
        headers = [child.compact_header for child in children]
        sequences = [child.my_sequence for child in children]

        # Resolve the chosen coloring: either identity/similarity, or per-column
        # shading taken from a conservation analysis already run on the alignment.
        _label, kind, attr, color_dict, prefix = color_mode
        color_scheme = "similarity"
        column_colors = None
        if kind in ("identity", "similarity"):
            color_scheme = kind
        else:
            column_colors = make_column_colors(children, attr, color_dict, prefix)

        ss_list = None
        if ss_element is not None:
            # Secondary structure is assigned lazily in PyMod (e.g. the first time
            # an element is colored by observed SS). A freshly loaded structure has
            # residue.secondary_structure == None, which would draw no SS row at all,
            # so assign it on demand here, exactly as color_element_by_obs_sec_str does.
            if not ss_element.has_assigned_secondary_structure():
                self.pymod.assign_secondary_structure(ss_element)
            ss_list = secondary_structure_by_column(ss_element)

        scene = build_alignment_scene(headers, sequences, ss_list=ss_list,
                                      color_scheme=color_scheme,
                                      block_length=block_length,
                                      column_colors=column_colors)

        self.figure_window = Alignment_figure_window_qt(
            self.pymod.main_window, scene=scene,
            title="Alignment figure: %s" % self.input_cluster_element.my_header)
        self.figure_window.show()


###################################################################################################
# GUI windows.                                                                                    #
###################################################################################################

class Alignment_figure_options_window_qt(PyMod_protocol_window_qt):

    def add_middle_frame_widgets(self):
        self.colorby_combobox = PyMod_combobox_qt(
            label_text="Color by",
            items=self.protocol.color_mode_labels)
        self.colorby_combobox.combobox.setCurrentIndex(1)  # Identity + similarity
        self.middle_formlayout.add_widget_to_align(self.colorby_combobox)

        self.ss_combobox = PyMod_combobox_qt(
            label_text="Secondary structure of",
            items=self.protocol.ss_structures_list)
        self.ss_combobox.combobox.setCurrentIndex(0)
        self.middle_formlayout.add_widget_to_align(self.ss_combobox)

        self.block_entry = PyMod_entryfield_qt(
            label_text="Residues per line", value="60")
        self.middle_formlayout.add_widget_to_align(self.block_entry)

        self.middle_formlayout.set_input_widgets_width(200)


class Alignment_figure_window_qt(QtWidgets.QMainWindow):
    """
    A window embedding the QGraphicsScene alignment figure in a QGraphicsView,
    with zoom controls and PNG/PDF export.
    """

    is_pymod_window = True

    def __init__(self, parent, scene, title="Alignment figure", *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle(title)
        self.scene = scene

        self.view = QtWidgets.QGraphicsView(scene)
        self.view.setRenderHint(QtGui.QPainter.Antialiasing, True)
        self.view.setRenderHint(QtGui.QPainter.TextAntialiasing, True)
        self.view.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
        self.view.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("white")))

        central = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(central)

        button_bar = QtWidgets.QHBoxLayout()
        zoom_in = QtWidgets.QPushButton("Zoom +")
        zoom_in.clicked.connect(lambda: self.view.scale(1.2, 1.2))
        zoom_out = QtWidgets.QPushButton("Zoom -")
        zoom_out.clicked.connect(lambda: self.view.scale(1 / 1.2, 1 / 1.2))
        save_button = QtWidgets.QPushButton("Save figure...")
        save_button.clicked.connect(self.save_figure)
        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        button_bar.addWidget(zoom_in)
        button_bar.addWidget(zoom_out)
        button_bar.addStretch()
        button_bar.addWidget(save_button)
        button_bar.addWidget(close_button)

        layout.addLayout(button_bar)
        layout.addWidget(self.view)
        self.setCentralWidget(central)

        rect = scene.sceneRect()
        self.resize(min(1100, int(rect.width()) + 60),
                    min(800, int(rect.height()) + 90))

    def save_figure(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save alignment figure", os.path.expanduser("~/alignment.png"),
            "PNG image (*.png);;PDF (*.pdf)")
        if not path:
            return
        try:
            ok = render_scene_to_file(self.scene, path)
        except Exception as e:
            ok = False
            QtWidgets.QMessageBox.critical(
                self, "Save Error", "Could not save the figure:\n%s" % str(e))
            return
        if not ok:
            QtWidgets.QMessageBox.critical(
                self, "Save Error", "Could not save the figure to:\n%s" % path)
