#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
from collections import defaultdict

# Internal modules #
from pacmill.taxonomy.taxa_graphs import TaxaBarstack, TaxaLegend

# First party modules #
from plumbing.cache import property_cached

# Third party modules #
import pandas

# Constants #
class Dummy: pass

###############################################################################
class TaxaTable:
    """
    Takes the OTU table along with the taxonomic assignment results to
    generate taxa tables at different ranks.
    """

    # Combine too many taxa together in the graphs #
    max_taxa_displayed = 18

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.taxonomy)

    def __init__(self, otu_table, taxonomy, base_dir):
        # Attributes #
        self.otu_table = otu_table
        self.taxonomy  = taxonomy
        self.base_dir  = base_dir
        # Short cuts #
        self.otu_df = self.otu_table.df
        # The rank names #
        if hasattr(self.taxonomy.database, 'rank_names'):
            self.rank_names = self.taxonomy.database.rank_names
        else:
            self.rank_names = list(map(str, range(1, 8)))

    @property
    def assignments(self):
        # Shortcuts #
        return self.taxonomy.results.assignments

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /taxa_table_domain.tsv
                /taxa_table_kingdom.tsv
                /taxa_table_phylum.tsv
                /taxa_table_class.tsv
                /taxa_table_order.tsv
                /taxa_table_family.tsv
                /taxa_table_tribe.tsv
                /taxa_table_genus.tsv
                /taxa_table_species.tsv
                /graphs/
                """

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.base_dir, self.all_paths)

    #------------------------------ Running ----------------------------------#
    def __call__(self, verbose=False):
        # Message #
        if verbose: print("Making all taxa tables in '%s'" % self.base_dir)
        # Make directory #
        self.base_dir.create_if_not_exists()
        # Loop over ranks #
        for i, rank_name in enumerate(self.rank_names):
            table = self.taxa_table_at_rank(i)
            tsv   = self.name_to_path(rank_name)
            table.to_csv(tsv.path, sep='\t', encoding='utf-8')
        # Return #
        return self.base_dir

    def name_to_path(self, rank_name):
        """Given a rank's name, return the path to the tsv file."""
        return self.base_dir + 'taxa_table_' + rank_name.lower() + '.tsv'

    def taxa_table_at_rank(self, rank):
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        # Loop over samples, then over OTUs #
        for sample_name, column in self.otu_df.T.iterrows():
            for otu_name, count in column.iteritems():
                # Because mothur renames OTUs #
                otu_name   = otu_name.replace(':', '_')
                # Retrieve a tuple unless it was discarded by post-processing #
                assignment = self.assignments.get(otu_name)
                # Get the assignment at this specific rank #
                if assignment is None:        taxa_term = "Unassigned"
                elif rank >= len(assignment): taxa_term = "Unassigned"
                elif assignment[rank] == '':  taxa_term = "Unassigned"
                else:                         taxa_term = assignment[rank]
                # Add the count we had #
                result[taxa_term][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Sort the table by sum #
        sums = result.sum()
        sums = sums.sort_values(ascending=False)
        result = result.reindex(sums.keys(), axis=1)
        # Return #
        return result

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the taxa tables were created already and the results
        are stored on the filesystem. Return False if it was not yet run.
        """
        return self.autopaths.phylum.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from taxa tables " \
                  "before generating them."
            raise Exception(msg)
        # Return the results #
        return TaxaTableResults(self)

###############################################################################
class TaxaTableResults:

    def __init__(self, parent):
        self.parent      = parent
        self.taxa_tables = parent

    def load_table(self, path):
        """Shortcut function to `pandas.read_csv`."""
        return pandas.read_csv(path, sep='\t', index_col=0, encoding='utf-8')

    @property_cached
    def taxa_tables_by_rank(self):
        """Return DataFrames in a list, one for each rank."""
        return [self.load_table(self.parent.name_to_path(n))
                for n in self.parent.rank_names]

    @property_cached
    def graphs(self):
        """
        The result is an object whose attributes are all the TaxaBarstack
        graphs (one at each rank) initialized with this instance as only
        argument. The graphs are also accessible in a list attribute named
        `by_rank`.

        Each graph also possesses an associated legend in a separate plot.
        These are available in the `legends` attribute.

        So to retrieve the taxonomic barstack graph and legend for the
        phylum rank, you could do the following:

            >>> print(proj.taxa_tables.results.graphs.taxa_barstack_phylum)
            >>> print(proj.taxa_tables.results.graphs.taxa_legend_phylum)
        """
        # Make a dummy object #
        result = Dummy()
        # Create a list attribute to hold each rank #
        result.by_rank = []
        # Create a list attribute to hold each legend #
        result.legends = []
        # Loop over ranks #
        for i, rank_name in enumerate(self.parent.rank_names):
            # The attributes of the graph we will create #
            attrs = dict(base_rank  = i,
                         short_name = 'taxa_barstack_' + rank_name.lower())
            # Create a graph type class for this specific rank #
            name = "TaxaBarstack" + rank_name.capitalize()
            clss = type(name, (TaxaBarstack,), attrs)
            # Instantiate the graph #
            graph = clss(self, base_dir=self.parent.autopaths.graphs_dir)
            # The attributes of the legend we will create #
            attrs = dict(base_rank  = i,
                         label      = rank_name,
                         short_name = 'taxa_legend_' + rank_name.lower())
            # Create a legend type class for this specific rank #
            name = "TaxaLegend" + rank_name.capitalize()
            clss = type(name, (TaxaLegend,), attrs)
            # Instantiate the legend #
            legend = clss(self, base_dir=self.parent.autopaths.graphs_dir)
            # Add them as an attribute of our result #
            setattr(result, graph.short_name,  graph)
            setattr(result, legend.short_name, legend)
            # Add them also to the lists #
            result.by_rank.append(graph)
            result.legends.append(legend)
        # Return #
        return result