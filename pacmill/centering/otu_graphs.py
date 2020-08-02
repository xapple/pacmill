#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# First party modules #
from plumbing.graphs import Graph
from plumbing.common import split_thousands as thousands

# Third party modules #
import matplotlib
from matplotlib import pyplot

# Constants #
__all__ = ['OtuSizesDist', 'OtuSumsPerSample', 'SampleSumsPerOtu',
           'CumulativePresence']

################################################################################
class OtuSizesDist(Graph):
    """Distribution of OTU cluster sizes in log-log."""

    short_name   = 'otu_sizes_dist'
    y_grid       = True
    x_scale      = 'symlog'
    y_scale      = 'log'
    x_label      = 'Number of sequences in an OTU'
    y_label      = 'Number of OTUs with that many sequences in them'
    y_grid       = True
    x_grid       = True
    width        = 10
    height       = 6

    def plot(self, **kwargs):
        # Sum by row and count frequencies #
        distrib = self.parent.df.sum(axis=1).value_counts()
        # Get x and y values #
        x = distrib.keys()
        y = distrib.values
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_title('Distribution of sizes for %s OTUs' % thousands(sum(y)))
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        # Close #
        pyplot.close(fig)

################################################################################
class OtuSumsPerSample(Graph):

    short_name   = 'otu_sums_graph'
    title        = 'Histogram of OTU appearance sums per sample'
    x_label      = 'Number of OTUs present (i.e. non-null) in a sample'
    y_label      = 'Number of samples with that many OTUs in them'
    y_grid       = True
    x_grid       = False
    width        = 10
    height       = 6
    remove_frame = True

    def plot(self, **kwargs):
        # Sum by column #
        self.frame = self.parent.df.astype(bool).sum()
        # Make histogram #
        fig = pyplot.figure()
        axes = self.frame.hist(color='gray', bins=40)
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        # Close #
        pyplot.close(fig)

################################################################################
class SampleSumsPerOtu(Graph):

    short_name   = 'sample_sums_graph'
    title        = 'Histogram of OTU appearance sums per OTU'
    y_label      = 'Number of OTUs that appear in these many samples'
    width        = 10
    height       = 6
    remove_frame = True

    def plot(self, **kwargs):
        # Sum by column and count frequencies #
        df = self.parent.df.astype(bool).sum(axis=1)
        df = df.value_counts().sort_index()
        # Make histogram #
        fig = pyplot.figure()
        axes = df.plot(kind='bar', color='gray')
        # Set X label #
        msg = 'Number of samples an OTU appears in (max. %i)'
        axes.set_xlabel(msg % self.parent.df.T.shape[0])
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        # Close #
        pyplot.close(fig)

################################################################################
class CumulativePresence(Graph):
    """
    Cumulative graph of cluster presence in samples.

    This graph can be interpreted in the following fashion:

        "If I take N percent of all the samples, how many OTUs
         are non-null in every single sample selected (at maximum)."

    This often translates to something such as:

    - 0%  of all OTUs are present in 100% of the samples,
    - 10% of all OTUs are present in 90% of the samples,
    - 90% of all OTUs are present in 1% of the samples.
    - etc.
    """

    short_name   = 'cumulative_presence'
    y_label      = 'Number of OTUs that appear in that fraction of samples or more'
    y_grid       = True
    x_grid       = False
    width        = 10
    height       = 6
    remove_frame = True

    def plot(self, **kwargs):
        # Number of samples #
        num_of_otus    = self.parent.df.shape[0]
        num_of_samples = self.parent.df.shape[1]
        # Create an index #
        samples_index  = list(reversed(range(1, num_of_samples+1)))
        # Get value frequencies #
        counts = self.parent.df.astype(bool).sum(axis=1).value_counts()
        # Add missing values #
        for n in samples_index:
            if n not in counts:
                counts.at[n] = 0
        # Sort it #
        counts = counts.sort_index(ascending=False)
        # Cumulative sum #
        self.y = list(counts.cumsum())
        # Percentage of samples #
        self.x = [n/num_of_samples for n in samples_index]
        # Make a step plot #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.step(self.x, self.y, fillstyle='bottom')
        # Set titles #
        axes.set_title('Cumulative graph of OTU presence in samples for %s OTUs' % num_of_otus)
        axes.set_xlabel('Fraction of samples (100%% equates %i samples)' % num_of_samples)
        axes.invert_xaxis()
        # Fine tuning #
        axes.set_xticks([min(self.x) + (max(self.x)-min(self.x)) * n / 9 for n in range(10)])
        axes.set_yscale('log')
        # Set percentage #
        percentage = lambda x, pos: '%1.0f%%' % (x*100.0)
        axes.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(percentage))
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        # Close #
        pyplot.close(fig)
