#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# First party modules #
from plumbing.cache import property_cached
from plumbing.graphs import Graph
from plumbing.graphs.cool_colors import colors as cool_colors
from plumbing.graphs.solo_legend import SoloLegend

# Third party modules #
import numpy
from matplotlib import pyplot

################################################################################
class TaxaBarstack(Graph):
    """
    Distribution of named taxa by sample (at a given rank).
    The Y axis represents fractions of abundance in percent.
    The X axis contains sample short names.
    The legend contains the taxonomic titles themselves.
    """

    base_rank     = -1
    short_name    = 'taxa_barstack'
    width         = 12.0
    height        = 7.5
    y_label       = 'Relative abundances in percent'

    def plot(self, **kwargs):
        # Retrieve the corresponding legend to this graph #
        legend = [leg for leg in self.parent.graphs.legends
                  if leg.base_rank == self.base_rank][0]
        # Transpose, each row is now a taxonomic title #
        df = legend.df.T
        # Number of bars and numbers of categories (sub-bars) within bars #
        num_taxa, num_samples = df.shape
        # Record the bottoms for each successive new bar #
        cum_size = numpy.zeros(num_samples)
        # Space each bar one unit apart from the next #
        x_locations = list(range(num_samples))
        # Loop #
        for taxa_name, row in df.iterrows():
            pyplot.bar(x_locations,
                       row.values,
                       linewidth = 0.5,
                       edgecolor = 'k',
                       bottom = cum_size,
                       color  = legend.label_to_color[row.name])
            # Increase the plotting location for the next bar #
            cum_size += row.values
        # Retrieve current figure and axes #
        fig  = pyplot.gcf()
        axes = pyplot.gca()
        # Compute the plot title #
        title     = 'Taxonomic relative abundances per sample at rank %i (%s).'
        rank_name = self.parent.parent.rank_names[self.base_rank]
        title     = title % (self.base_rank, rank_name)
        # Set the plot title #
        axes.set_title(title)
        # Y limits should always be the same #
        axes.set_ylim([0, 100])
        # Add the sample names #
        axes.set_xticks(x_locations)
        axes.set_xticklabels(list(df.columns))
        # Change font of sample names #
        axes.set_xticklabels(axes.get_xticklabels(),
                             fontweight = 'bold',
                             size       = 'large',
                             fontfamily = 'monospace')
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        # Close #
        pyplot.close(fig)

################################################################################
class TaxaLegend(SoloLegend):

    # Base attributes #
    base_rank  = -1
    short_name = 'taxa_legend'

    # Number of columns #
    n_col = 4

    # Trim long names #
    max_name_len = 42

    @property_cached
    def df(self):
        """The taxa table, normalized, sorted, trimmed and renamed."""
        # Retrieve the taxa table for the current rank #
        taxa_table = self.parent.taxa_tables_by_rank[self.base_rank]
        # Normalize it #
        df = taxa_table.apply(lambda x: 100 * x / x.sum(), axis=1)
        # Sort columns with biggest first #
        df = df.reindex(df.mean().sort_values(ascending=False).index, axis=1)
        # Combine all columns above a certain number to 'Others' #
        max_taxa_count = self.parent.parent.max_taxa_displayed
        if len(df.columns) > max_taxa_count:
            cols_to_drop = df.columns[max_taxa_count-1:]
            cols_summed  = df[cols_to_drop].sum(axis=1)
            df = df.drop(columns=cols_to_drop)
            df['(All Others Combined)'] = cols_summed
        # Trim taxonomic names to a certain number of characters #
        new_columns = {col: col if len(col) <= self.max_name_len
        else col[:self.max_name_len-4] + ' ...'
                       for col in df.columns}
        df = df.rename(new_columns)
        # Return #
        return df

    @property_cached
    def label_to_color(self):
        """Mapping of taxonomic titles to colors as a dictionary."""
        # Get all the taxonomic titles #
        taxa = list(self.df.columns)
        # Modify the colors so that 'others' are always black #
        max_taxa = self.parent.parent.max_taxa_displayed - 1
        colors = cool_colors[:max_taxa] + ['#000000'] + cool_colors[max_taxa:]
        # Assign colors #
        from itertools import cycle
        result = dict(zip(taxa, cycle(colors)))
        # Return #
        return result