# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.


try:
    from ete3 import faces, COLOR_SCHEMES, AttrFace, CircleFace, RectFace, TextFace, SeqMotifFace, TreeStyle, PhyloTree, NodeStyle
except:
    pass

from collections import defaultdict
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QGraphicsRectItem, QGraphicsSimpleTextItem, QGraphicsEllipseItem, QFileDialog
from PyQt5.QtGui import QColor, QPen,  QBrush

import colorsys
import random
import copy

'''
 Method for common ancestor visualization
'''
def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize= 5)
        faces.add_face_to_node(N, node, 0, position="aligned")


'''
 Method for Bubble/support values tree visualization
'''

def bubble_ly(self):
    if self.is_leaf():
        # Add node name to laef nodes
        N = AttrFace("name", fsize=10, fgcolor="black")
        faces.add_face_to_node(N, self, 0)
    if "weight" in self.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        C = CircleFace(radius=self.weight, color="RoyalBlue", style="sphere")
        # Let's make the sphere transparent
        C.opacity = 0.2
        # And place as a float face over the tree
        faces.add_face_to_node(C, self, 0, position="float")


'''
 Class for Mutation Events visualization
'''
# Given a tree and its associated alignment, plot column changes for each internal branch.
class ColumnMutation():

    def mutation_columns(self,sequences):
        col2diffs = defaultdict(set)
        self.aln_length = len(sequences[0])

        for col in range(self.aln_length):
            for seq in sequences:
                col2diffs[col].add(seq[col])
            col2diffs[col].discard('-')
        subseqs = set()
        relevant_columns = []
        for col in range(self.aln_length):
            if len(col2diffs[col]) > 1:
                relevant_columns.append(col)
        for seq in sequences:
            subseqs.add(''.join([seq[col] for col in relevant_columns]))

        return subseqs, relevant_columns

    def get_mutation_event (self,tree, aln, single_event):
        ts = TreeStyle()
        # disable default PhyloTree Layout
        ts.layout_fn = lambda x: True
        tree.link_to_alignment(aln)
        node2content = tree.get_cached_content()
        for node in tree.traverse():
            node.img_style["size"] = 0

            if not node.is_leaf():
                leaves = node2content[node]
                # get columns with different aa
                subseqs, relevant_columns  = self.mutation_columns([lf.sequence for lf in leaves])
                if not single_event:
                    for seq in subseqs:
                        try:
                            f = SeqMotifFace(seq, seq_format="seq", width=10, height=8)
                            f.margin_top = 2
                            f.margin_right = 4
                            node.add_face(f, column=0, position="branch-bottom")
                        except ValueError:
                            pass

                    for j, col in enumerate(relevant_columns):
                        col_f = RectFace(10, 10, fgcolor=None, bgcolor=None,
                                         label={"text":str(col), "fonttype":"Courier", "color":"black", "fontsize":4})
                        node.add_face(col_f, column=j, position="branch-top")
                        col_f.margin_bottom = 2
                else:
                    for seq in subseqs:
                        if len(seq) == 1:
                            try:
                                f = SeqMotifFace(seq, seq_format="seq", width=10, height=8)
                                f.margin_top = 2
                                f.margin_right = 4
                                node.add_face(f, column=0, position="branch-bottom")
                            except ValueError:
                                pass
            else:
                f = SeqMotifFace(node.sequence, seq_format="seq", width=1)
                node.add_face(f, column=0, position="aligned")

        ts.draw_aligned_faces_as_table = False
        for colnum in range(self.aln_length):
            col_f = RectFace(10, 10, fgcolor=None, bgcolor=None,
                             label={"text":str(colnum), "fonttype":"Courier", "color":"black", "fontsize":4})
            ts.aligned_header.add_face(col_f, column=colnum)

        return tree, ts
###########################################################################################################################
