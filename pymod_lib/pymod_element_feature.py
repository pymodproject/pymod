# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

class Element_feature:

    def __init__(self, id, name, start, end, description='', feature_type='None', color=None):
        self.start = start
        self.end = end
        self.id = id
        self.name = name
        self.feature_type = feature_type
        self.description = description
        self.color = color
        self.show = False # Show in PyMod main window.

    def __repr__(self):
        return str(self.__dict__)

    def get_feature_dict(self):
        return {"start": self.start,
                "end": self.end,
                "id": self.id,
                "name": self.name,
                "feature_type": self.feature_type,
                "description": self.description,}


class Domain_feature(Element_feature):

    def __init__(self, id, name, start, end, evalue, color=None, description='', is_derived=False, offset=(0, 0)):

        Element_feature.__init__(self, id, name, start, end, description, feature_type='domain')
        self.domain_color = color
        self.evalue = evalue
        self.element = None

        self.full_name = "%s_%s_%s" % (self.name, self.start, self.end)

        # Initially is always set to 'False'. It will be set to 'True' when the copy of a
        # 'Domain_feature' object is assigned to a split domain element derived from an original
        # query sequence.
        self.is_derived = is_derived
        self.offset = offset # symmetric, nterm offset and cterm offset if the seq is split into domains

        if self.is_derived:
            # Store the original start/end values from the parent.
            self.parent_start = self.start
            self.parent_end = self.end
            # Updates the values to suite the new cropped domain sequence.
            new_startindex = max(0, self.start-self.offset[0])
            self.start -= new_startindex
            self.end -= new_startindex
        else:
            self.parent_start = None
            self.parent_end = None

        self.show = True

    def get_feature_dict(self):
        return {"id": self.id,
                "name": self.name,
                "start": self.start,
                "end": self.end,
                "evalue": self.evalue,
                "color": self.domain_color,
                "description": self.description,
                "is_derived": self.is_derived,
                "offset": self.offset,
                }
