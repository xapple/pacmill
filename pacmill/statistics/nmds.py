#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# First party modules #
from autopaths.tmp_path import new_temp_file
from plumbing.graphs import Graph

# Third party modules #
import numpy, pandas
from matplotlib import pyplot

################################################################################
class GraphNMDS(Graph):
    """
    Stands for Non-metric multidimensional scaling.
    Using the information in the OTU table along with a distance metric
    such as the one developed by Horn 1966 (adapted from Morisita 1959),
    we can place every sample on a two-dimensional ordination plot.

    Pycogent used to do it, but the new Python 3 version dropped that feature:

    http://lira.no-ip.org:8080/doc/python-cogent-doc/html/examples/perform_nmds.html
    https://github.com/pycogent/pycogent/blob/master/cogent/cluster/nmds.py

    Qiime used to do it too, but they were just calling Pycogent:

    http://qiime.org/scripts/nmds.html
    https://github.com/biocore/qiime/blob/master/qiime/nmds.py
    https://github.com/qiime2/qiime2/search?q=nmds&unscoped_q=nmds

    So now scikit learn has an implementation:

    https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html

    And ecopy provides something that calls scikit learn:

    https://ecopy.readthedocs.io/en/latest/ordination.html

    Another possible implementation:

    https://github.com/jianshu93/NMDS

    Otherwise we can use rpy2 and the R vegan package.
    """

    short_name = 'nmds_horn'

    def plot(self, **kwargs):
        # We want to pass a transposed version to R #
        transposed = new_temp_file(suffix='.tsv')
        self.parent.df.T.to_csv(str(transposed), sep='\t')
        transposed.prepend('X')
        # Run via R #
        coord = self.run_via_r(transposed.path)
        # Data #
        x      = coord['NMDS1'].values
        y      = coord['NMDS2'].values
        names  = coord['NMDS1'].keys()
        # Make scatter #
        fig  = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_title('Non-Metric Multidimensional scaling')
        axes.set_xlabel('Dimension 1')
        axes.set_ylabel('Dimension 2')
        # Add annotations #
        for i in range(len(names)):
            bbox = {'boxstyle': 'round,pad=0.2', 'fc': 'yellow', 'alpha': 0.3}
            pyplot.annotate(names[i],
                            size       = 9,
                            xy         = (x[i], y[i]),
                            xytext     = (10, 0),
                            textcoords = 'offset points',
                            ha         = 'left',
                            va         = 'center',
                            bbox       = bbox)
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        # Close #
        pyplot.close(fig)

    def run_via_r(self, otu_table_path):
        # Module on demand #
        from rpy2 import robjects
        # Load dataframe #
        robjects.r("library(vegan)")
        # The command to read the table #
        cmd = "table = read.table('%s', sep='\t', header=TRUE, row.names='X')"
        robjects.r(cmd % otu_table_path)
        # Run computation #
        robjects.r("nmds = metaMDS(table, distance='horn', trymax=200)")
        # Extract result #
        robjects.r("coord = scores(nmds)")
        robjects.r("loadings = nmds$species")
        # Retrieve values #
        coord = robjects.r.coord
        # Convert #
        coord = self.r_matrix_to_dataframe(coord)
        # Return #
        return coord

    def r_matrix_to_dataframe(self, matrix):
        cols = list(matrix.colnames)
        rows = list(matrix.rownames)
        return pandas.DataFrame(numpy.array(matrix),
                                index=rows,
                                columns=cols)
