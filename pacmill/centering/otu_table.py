#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
import pacmill.centering.otu_graphs

# First party modules #
from plumbing.cache import property_cached
from autopaths.file_path import FilePath

# Third party modules #
import pandas

# Constants #
class Dummy: pass

###############################################################################
class OtuTable:

    def __init__(self, tsv_path):
        # Record the tsv path #
        self.tsv_path = FilePath(tsv_path)
        # Pick a directory for storing the graphs #
        self.graphs_dir = self.tsv_path.directory + 'graphs/'

    @property_cached
    def df(self):
        """
        Parse the TSV file and return a pandas.DataFrame object.
        Returns a table with OTU as rows and samples as columns which tracks
        how many sequences where found from each sample in each OTU.
        """
        # Parse #
        df = pandas.read_csv(str(self.tsv_path), '\t', index_col=0)
        # Return #
        return df

    @property_cached
    def graphs(self):
        """
        Minor use of black magic here. The result is an object whose attributes
        are all the graphs found in `otu_graphs.py` initialized with this
        instance as only argument.
        """
        # Make a dummy object #
        result = Dummy()
        # Loop over graphs #
        for graph_name in pacmill.centering.otu_graphs.__all__:
            graph_cls = getattr(pacmill.centering.otu_graphs, graph_name)
            setattr(result, graph_cls.short_name, graph_cls(self))
        # Return #
        return result
