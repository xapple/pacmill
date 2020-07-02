#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
Contact at www.sinclair.bio
MIT Licensed.
"""

# Built-in modules #

# Third party modules #

# First party modules #

# Internal modules #

###############################################################################
class Sample:
    """
    A Sample object consists first and foremost of either two FASTA files
    or two FASTQ files.

    These files contain paired sequences all coming from the same particular
    real-life lab sample. Might or might not correspond to an Illumina MID.

    In addition, a Sample has many other metadata associated with it.

        * self.short_name: the name of this sample.
        * self.parent:     a reference to a Project instance.
    """

    def __repr__(self):
        return '%s object code "%s"' % (self.__class__, self.short_name)

    def __init__(self, **kwargs):
        """
        A Sample object takes a dictionary of metadata as input.
        Each named parameter will be set as an attribute of this instance
        with the same name. Hence, the attributes of a Sample correspond
        directly to the column headers of the excel metadata file that
        was parsed by the Project object.
        """
        # Set the attributes of this instance with the given kwargs #
        for name, value in kwargs.items(): setattr(self, name, value)
        # We need to validate some attributes and also transform some #
        self.transform_attrs()
        self.validate_attrs()

    def transform_attrs(self):
        """
        This is where we add and change some attributes based on the
        information in the excel metadata file.
        """
        # The sample short name can just be called short_name #
        self.short_name = self.sample_short_name
        # Same for the number and long name #
        self.num = int(self.sample_num)
        self.long_name = self.sample_long_name

    def validate_attrs(self):
        """
        This is where we check that the information in the excel
        metadata file is consistent.
        """
        # Check the short name contains only alphanumerics and underscore #
        assert self.short_name.isidentifier()