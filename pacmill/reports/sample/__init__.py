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
from pymarktex.figures import DualFigure
from autopaths         import Path

# Third party modules #
import pandas

###############################################################################
class SampleReport(Document):
    """A full report generated in PDF for every Sample object."""

    # Custom LaTeX headers and footers for the report #
    header_template = Header
    footer_template = Footer

    # Specific title for sample report #
    params = {'title': 'Auto-generated sample report'}

    def __init__(self, sample, output_path):
        # Reference to parent objects #
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
        # Reference to parent objects #
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
        # Let's also remove those that have no value (nan) #
        keys = [k for k in keys if not pandas.isna(getattr(self.sample, k))]
        # Format them as a bullet list #
        bullet = lambda k: '* **' + k + '**: ' + str(getattr(self.sample, k))
        result = '\n'.join(bullet(k) for k in keys)
        # Return #
        return result

    #----------------------------- Raw data ----------------------------------#
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

    #----------------------------- Filtering ---------------------------------#
    def filtering(self):
        if not self.sample.filter: return False

    def primer_max_dist(self):
        return self.sample.filter.primer_max_dist

    def mismatches_allowed(self):
        return self.sample.filter.primer_mismatches

    @property_pickled
    def primer_discard(self):
        before = self.sample.count
        after  = self.sample.filter.results.primers_fasta
        return thousands(len(before) - len(after))

    @property_pickled
    def primer_left(self):
        return thousands(len(self.sample.filter.results.primers_fasta))

    @property_pickled
    def n_base_discard(self):
        before = self.sample.filter.results.primers_fasta
        after  = self.sample.filter.results.n_base_fasta
        return thousands(len(before) - len(after))

    @property_pickled
    def n_base_left(self):
        return thousands(len(self.sample.filter.results.n_base_fasta))

    def min_read_length(self):
        return self.sample.filter.min_read_length

    def max_read_length(self):
        return self.sample.filter.max_read_length

    @property_pickled
    def length_discard(self):
        before = self.sample.filter.results.n_base_fasta
        after  = self.sample.filter.results.length_fasta
        return thousands(len(before) - len(after))

    @property_pickled
    def length_left(self):
        return thousands(len(self.sample.filter.results.length_fasta))

    @property_pickled
    def percent_remaining(self):
        percent = len(self.sample.filter.results.clean) / len(self.sample)
        percent *= 100
        return "%.1f%%" % percent
