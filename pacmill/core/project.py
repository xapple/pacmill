#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Third party modules #
import pandas

# First party modules #
from autopaths import Path
from plumbing.cache import property_cached

# Internal modules #
from pacmill.core.sample import Sample

###############################################################################
class Project:
    """A Project object regroups several Sample objects together."""

    def __repr__(self):
        return '%s object code "%s"' % (self.__class__, self.short_name)

    def __iter__(self):
        return iter(self.samples)

    def __len__(self):
        return len(self.samples)

    def __init__(self, short_name, *all_xlsx):
        """
        A Project object takes:

         * A project `short_name` attribute that will be used to pick only the
           samples that are associated to that particular project.
         * One or more <xlsx> files that describe samples locations and their
           metadata. See example excel file.

        Once the project is created you can access its attributes such as:

            >>> proj = Project('demo_project', '~/test/metadata.xlsx')
            >>> for sample in proj: print(sample.read_count)

        Attributes are the following:

            * self.short_name: the name of this project.
            * self.long_name:  the more lengthy description of this project.
            * self.all_xlsx:   a list of file paths (that are excel files).
            * self.metadata:   a pandas dataframe with all metadata
                               for this project combined.
            * self.samples:    a list of Sample objects.

        Other properties are described in their respective docstrings.
        """
        # The name of the project #
        self.short_name = short_name.lower()
        # Check it contains only alphanumerics and underscore #
        assert self.short_name.isidentifier()
        # You need at least one excel file #
        assert len(all_xlsx) > 0
        # Save all excel file paths as AutoPaths objects #
        self.all_xlsx = list(map(Path, all_xlsx))
        # Check all the excel paths actually exist #
        assert all(xlsx.exists for xlsx in self.all_xlsx)

    #----------------------------- Properties --------------------------------#
    @property_cached
    def metadata(self):
        """
        Return a pandas.DataFrame object describing the metadata of all
        samples contained in this project.
        """
        # Function to read one excel file #
        read_excel = lambda path: pandas.read_excel(str(path), header=1)
        # Read them all as data frames #
        all_dfs = [read_excel(path) for path in self.all_xlsx]
        # If there are several excel files, merge them together #
        metadata = pandas.concat(all_dfs, sort=False)
        # Filter and take only samples that match this project's short_name #
        query = f'project_short_name == "{self.short_name}"'
        metadata = metadata.query(query).copy()
        # Remove samples that are not marked as "yes" for "used" #
        metadata = metadata.query('used == "yes"').copy()
        # Return #
        return metadata

    @property_cached
    def samples(self):
        """Create all the Sample objects."""
        # Iterate over all rows of the input #
        rows = self.metadata.iterrows()
        # Make one Sample object per row #
        samples = [Sample(self, **dict(row)) for i, row in rows]
        # Add a reference the current project #
        for sample in samples: sample.parent = self
        # Check we have at least one sample #
        assert len(samples) > 0
        # Return #
        return samples

    @property_cached
    def long_name(self):
        """
        Get the project's long name which is a string describing the
        project title in a longer fashion than the short name.
        """
        # It's actually contained in the Sample objects #
        all_names = set(s.project_long_name for s in self)
        # Check that it doesn't diverge between samples #
        if not len(all_names) == 1:
            msg = "The project long name is not uniform across samples."
            raise ValueError(msg)
        # Return #
        return all_names.pop()
