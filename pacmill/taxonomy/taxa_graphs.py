#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph
from plumbing.graphs.cool_colors import colors as cool_colors

# Third party modules #
from matplotlib import pyplot

################################################################################
class TaxaBarstack(Graph):
    """
    Distribution of named taxa by sample (at different ranks).
    The Y axis represents fractions.
    The X axis should be sample short names.
    The legend contains taxonomic titles.
    """

    base_rank     = -1
    short_name    = 'taxa_barstack'
    width         = 12.0
    height        = 7.5
    bottom        = 0.22
    top           = 0.94
    left          = 0.08
    right         = 0.96
    legend_anchor = -0.15
    y_label       = 'Relative abundances in percent'

    def plot(self, **kwargs):
        # Retrieve the taxa table for the current rank #
        taxa_table = self.parent.taxa_tables_by_rank[self.base_rank]
        # Normalize it #
        df = taxa_table.apply(lambda x: 100 * x / x.sum(), axis=1)
        # Special case where there is only one taxa e.g. only 'Bacteria' #
        if len(df.columns) < 2 : colors = 'gray'
        else:                    colors = cool_colors
        # Plot #
        axes = df.plot(kind='bar', stacked=True, color=colors)
        fig  = pyplot.gcf()
        # Other #
        title     = 'Taxonomic relative abundances per sample at rank %i (%s).'
        rank_name = self.parent.parent.rank_names[self.base_rank]
        title     = title % (self.base_rank, rank_name)
        axes.set_title(title)
        # Y limits should always be the same #
        axes.set_ylim([0, 100])
        # Put a legend below the current axis #
        axes.legend(loc            = 'upper center',
                    bbox_to_anchor = (0.5, self.legend_anchor),
                    fancybox       = True,
                    shadow         = True,
                    ncol           = 5)
        # Change font of sample names #
        axes.set_xticklabels(axes.get_xticklabels(),
                             fontweight = 'bold',
                             size       = 'large',
                             fontfamily = 'monospace')
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        # Close #
        pyplot.close(fig)