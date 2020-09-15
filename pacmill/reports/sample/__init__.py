#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from pacmill.reports.base_template import ReportTemplate
from pacmill.reports.template      import Header, Footer

# First party modules #
from plumbing.cache    import property_cached, property_pickled
from plumbing.common   import split_thousands as thousands
from pymarktex         import Document
from pymarktex.figures import DualFigure, ScaledFigure
from autopaths         import Path

# Third party modules #
import pandas

###############################################################################
class SampleReport(Document):
    """A full report generated in PDF for every Sample object."""

    # Custom LaTeX headers and footers for the report #
    header_template = Header
    footer_template = Footer

    # Specific title for the sample report #
    params = {'title': 'Auto-generated sample report'}

    def __init__(self, sample, output_path):
        # Reference to parent objects #
        self.parent  = sample
        self.sample  = sample
        self.project = sample.project
        # The output location #
        self.output_path = Path(output_path)

    @property_cached
    def template(self):
        return SampleTemplate(self)

    def load_markdown(self):
        self.markdown = str(self.template)

###############################################################################
class SampleTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template."""

    delimiters = (u'{{', u'}}')

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, parent):
        # References to parent object #
        self.parent = parent
        self.report = parent
        # Reference to sample and project objects #
        self.sample  = parent.sample
        self.project = parent.project
        # Where to pickle properties that are costly to compute #
        self.cache_dir = self.sample.autopaths.report_cache_dir

    #------------------------- General information ---------------------------#
    def sample_short_name(self):     return self.sample.short_name
    def sample_long_name(self):      return self.sample.long_name
    def project_short_name(self):    return self.project.short_name
    def project_long_name(self):     return self.project.long_name
    def project_other_samples(self): return len(self.project) - 1

    #--------------------------- Metadata list -------------------------------#
    def all_metadata(self):
        """
        Summarize all the metadata given in the excel sheet for this
        current sample.
        """
        # All the keys provided #
        keys = self.sample.metadata_keys
        # Let's skip some that were already mentioned above #
        skip_keys = ['sample_num', 'used', 'sample_short_name',
                     'project_short_name', 'sample_long_name',
                     'project_long_name']
        # Filter them out #
        keys = [k for k in keys if k not in skip_keys]
        # Let's also remove those that have no value #
        keys = [k for k in keys if getattr(self.sample, k) is not None]
        # Format them as a bullet list #
        value  = lambda k: str(getattr(self.sample, k))
        bullet = lambda k: '* **`%s`**: `%s`' % (k, value(k))
        result = '\n'.join(bullet(k) for k in keys)
        # Return #
        return result

    #------------------------------ Quality ----------------------------------#
    def fastq_size(self): return str(self.sample.fastq.size)

    @property_pickled
    def fastq_count(self): return thousands(self.sample.fastq.count)

    @property_pickled
    def fastq_qual(self):  return "%.2f" % self.sample.fastq.avg_quality

    def fastqc_graphs(self):
        # Two paths #
        params = [self.sample.fastq.fastqc.results.per_base_qual,
                  self.sample.fastq.fastqc.results.per_seq_qual]
        # Caption #
        params += ["Per base quality", "Per sequence quality"]
        # Labels #
        params += ["per_base_qual",    "per_seq_qual"]
        # Main caption and main label #
        params += ["Quality graphs from FastQC", "fastqc_graphs"]
        # Return #
        return str(DualFigure(*params))

    #------------------------------ Lengths ----------------------------------#
    @property_pickled
    def shortest_seq(self):
        return thousands(min(self.sample.fastq.lengths))

    @property_pickled
    def longest_seq(self):
        return thousands(max(self.sample.fastq.lengths))

    def raw_len_dist(self):
        caption = "Distribution of sequence lengths in original reads."
        path    = self.sample.fastq.graphs.length_hist()
        label   = "raw_len_dist"
        return str(ScaledFigure(path, caption, label))

    #----------------------------- Filtering ---------------------------------#
    def filtering(self):
        return bool(self.sample.filter)

    def primer_max_dist(self):
        return self.sample.filter.primer_max_dist

    def mismatches_allowed(self):
        return self.sample.filter.primer_mismatches

    @property_pickled
    def primer_discard(self):
        before = self.sample.fastq.count
        after  = self.sample.filter.results.primers_fastq.count
        return thousands(before - after)

    @property_pickled
    def primer_left(self):
        return thousands(self.sample.filter.results.primers_fastq.count)

    @property_pickled
    def n_base_discard(self):
        before = self.sample.filter.results.primers_fastq.count
        after  = self.sample.filter.results.n_base_fastq.count
        return thousands(before - after)

    @property_pickled
    def n_base_left(self):
        return thousands(self.sample.filter.results.n_base_fastq.count)

    def min_read_length(self):
        return thousands(self.sample.filter.min_read_len)

    def max_read_length(self):
        return thousands(self.sample.filter.max_read_len)

    @property_pickled
    def length_discard(self):
        before = self.sample.filter.results.n_base_fastq.count
        after  = self.sample.filter.results.length_fastq.count
        return thousands(before - after)

    @property_pickled
    def length_left(self):
        return thousands(self.sample.filter.results.length_fastq.count)

    def phred_window_size(self):
        return thousands(self.sample.phred_window_size)

    def phred_threshold(self):
        return thousands(self.sample.phred_threshold)

    @property_pickled
    def score_discard(self):
        before = self.sample.filter.results.length_fastq.count
        after  = self.sample.filter.results.score_fastq.count
        return thousands(before - after)

    @property_pickled
    def score_left(self):
        return thousands(self.sample.filter.results.score_fastq.count)

    @property_pickled
    def percent_remaining(self):
        before = self.sample.fastq.count
        after  = self.sample.filter.results.clean.count
        percent = 100 * (after / before)
        return "%.1f%%" % percent

    #------------------------------ Chimeras ---------------------------------#
    def chimeras(self):
        return bool(self.sample.chimeras)

    @property_pickled
    def chimeras_discard(self):
        before = self.sample.filter.results.clean.count
        after  = self.sample.chimeras.results.count
        return thousands(before - after)

    @property_pickled
    def chimeras_left(self):
        return thousands(self.sample.chimeras.results.count)

    #------------------------------ Barrnap ----------------------------------#
    def barrnap(self):
        return bool(self.sample.barrnap)

    def barrnap_mode(self):
        return self.sample.barrnap_mode

    @property_pickled
    def barrnap_discard(self):
        before = self.sample.chimeras.results.count
        after  = self.sample.barrnap.results.count
        return thousands(before - after)

    @property_pickled
    def barrnap_left(self):
        return thousands(self.sample.barrnap.results.count)

    #------------------------------ Taxonomy ---------------------------------#
    def taxonomy(self):
        return bool(self.project.taxonomy.tables.silva)

    @property_pickled
    def taxa_table(self):
        # Pick the rank #
        rank_name = "Genus"
        # Get the taxa tables of choice (silva) #
        tables = self.project.taxonomy.tables.silva
        # Get the rank number #
        rank = tables.rank_names.index(rank_name)
        # Get the row for this sample #
        table = tables.results.taxa_tables_by_rank[rank]
        row   = table.loc[self.sample.short_name]
        row   = row.sort_values(ascending=False)
        # Make an empty dataframe #
        df = pandas.DataFrame(index=range(len(row)))
        # Populate dataframe #
        df['#']      = range(1, len(row)+1)
        df['Genera'] = row.index
        df['Reads']  = [thousands(r) for r in row.values]
        # Show only most abundant #
        df = df[0:20]
        df = dict(df)
        # Make it as text #
        from tabulate import tabulate
        table = tabulate(df, headers="keys", numalign="right", tablefmt="pipe")
        # Add caption #
        caption = "The 20 most abundant predicted genera in this sample" \
                  " predicted by silva."
        return table + "\n\n   : %s" % caption
