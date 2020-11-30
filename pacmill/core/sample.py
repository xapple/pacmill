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
from plumbing.cache import property_cached
from autopaths      import Path

# Internal modules #

###############################################################################
class Sample:
    """
    A Sample object consists first and foremost of one FASTQ file.

    For a pipeline that uses paired reads (and hence two FASTQ files) such as
    the ones produced by Illumina sequencers please have a look at a sister
    project: www.github.com/xapple/sifes

    In addition, a Sample has many other metadata associated with it.
    Every column of the original excel file is now the name of an attribute.
    For instance (non exhaustive list):

        * self.num: the number of this sample in the project.
        * self.short_name: the name of this sample.
        * self.long_name: the name of this sample.
        * self.grouping: optional grouping information.

    In addition, the Sample object has the following:

        * self.project: a reference to a Project instance that owns this
                        sample.

    Other properties are described in their respective docstrings.
    """

    def __repr__(self):
        return '%s object code "%s"' % (self.__class__, self.short_name)

    def __init__(self, parent, **kwargs):
        """
        A Sample object takes a dictionary of metadata as input and a
        reference to the Project object that creates it.
        Each named parameter will be set as an attribute of this instance
        with the same name. Hence, the attributes of a Sample correspond
        directly to the column headers of the excel metadata file that
        was parsed by the Project object.
        """
        # Keep a reference to the parent Project object #
        self.parent  = parent
        self.project = parent
        # Record which metadata keys were passed to this sample #
        self.metadata_keys = [k for k in kwargs if not k.startswith("Unnamed")]
        # Set the attributes of this instance with the given kwargs #
        for key in self.metadata_keys: setattr(self, key, kwargs[key])
        # We need to transform some attributes #
        self.transform_attrs()
        # We need to validate some attributes #
        self.validate_attrs()
        # We need to set some attributes to defaults when absent #
        self.set_default_attrs()

    #------------------------------- Methods ---------------------------------#
    dna_keys = ['fwd_index_seq', 'rev_index_seq', 'fwd_primer_seq',
                'rev_primer_seq']

    int_keys = ['sample_num', 'fwd_read_count', 'rev_read_count', 'fwd_read_len', 'rev_read_len', 'primer_mismatches', 'primer_max_dist', 'min_read_len', 'max_read_len', 'phred_window_size', 'phred_threshold', 'otu_min_size', 'max_taxa']

    def transform_attrs(self):
        """
        This is where we add and change some attributes based on the
        information in the excel metadata file.
        """
        # The sample's short name can just be called short_name #
        self.short_name = self.sample_short_name
        # The sample's long name can just be called long_name #
        self.long_name = self.sample_long_name
        # The sample number can just be called num #
        self.num = int(self.sample_num)
        # Missing values should be None and not pandas.nan #
        for key in self.metadata_keys:
            if pandas.isna(getattr(self, key)): setattr(self, key, None)
        # If there are spaces in any DNA sequence, remove them #
        for key in self.dna_keys:
            seq = getattr(self, key)
            if seq: setattr(self, key, seq.replace(' ', ''))
        # Convert all values that are supposed to be integers
        for key in self.int_keys:
            value = getattr(self, key, None)
            if value is not None: setattr(self, key, int(value))

    def validate_attrs(self):
        """
        This is where we check that the information in the excel
        metadata file is consistent and usable.
        """
        # Check the short name contains only alphanumerics and underscore #
        assert self.short_name.isidentifier()
        # Check that all directories always end with a dash (/) #
        for key in self.metadata_keys:
            if not key.endswith('_dir'): continue
            value = getattr(self, key)
            if not isinstance(value, str) or not value.endswith("/"):
                msg = "The `%s` entry of <%s> must end with a slash. " \
                      "It currently does not: '%s'"
                msg = msg % (key, self.description, value)
                raise ValueError(msg)
        # Check that the FASTQ file exists #
        if not self.path.exists:
            msg = "The FASTQ file path of <%s> cannot be found. " \
                  "It should be located at: '%s'"
            msg = msg % (self.description, self.path)
            raise FileNotFoundError(msg)
        # Check the DNA sequences are properly formatted #
        for key in self.dna_keys:
            seq = getattr(self, key)
            if seq and '\n' in seq:
                msg = "The `%s` entry of <%s> contains illegal " \
                      "characters, please check: '%s'"
                msg = msg % (key, self.description, seq)
                raise ValueError(msg)
        # Check that the barrnap mode is a valid option #
        if hasattr(self, 'barrnap_mode'):
            assert self.barrnap_mode in ['off', 'filter', 'concat', 'trim']

    # Declare the default values #
    defaults = {
        'barrnap_mode':   'off',
        'otu_threshold':   0.97,
        'otu_min_size':    1,
        'run_silva':       True,
        'run_greengenes':  False,
        'run_rdp':         False,
        'run_crest':       False,
        'run_ncbi_blast':  False,
    }

    def set_default_attrs(self):
        """
        This is where we set default values for optional columns in the excel
        file. If the column is present, no change is made. If it is absent,
        we will add that particular attribute to the current instance.
        """
        # Check everyone of them #
        for key, value in self.defaults.items():
            if not hasattr(self, key): setattr(self, key, value)

    #----------------------------- Properties --------------------------------#
    @property
    def description(self):
        """A string describing the current sample."""
        desc = "Sample object %i in project '%s' with name '%s'"
        return desc % (self.num, self.project.short_name, self.short_name)

    @property_cached
    def path(self):
        """The path to the raw FASTQ reads file."""
        # Join the three components together #
        return Path(self.input_dir + self.suffix_dir + self.fwd_file_name)

    @property_cached
    def percent_lost(self):
        """Once all filtering is done, what fraction of reads were lost."""
        return 100 - 100 * self.chimeras.results.count / self.fastq.count

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /fastqc/
                /graphs/raw_len_dist.pdf
                /graphs/raw_len_hist.pdf
                /filtered/
                /chimeras/cleaned.fasta
                /chimeras/rejects.fasta
                /barrnap/results.gff
                /barrnap/only_genes.fasta
                /report/cache/
                /report/sample.pdf
                """

    @property_cached
    def base_dir(self):
        """
        The path to the directory where all results will be stored
        for this sample. We build it by joining three components together.
        See https://github.com/xapple/autopaths#directorypath-object
        """
        return Path(self.output_dir + 'samples/' + self.short_name + '/')

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.base_dir, self.all_paths)

    #---------------------------- Compositions -------------------------------#
    @property_cached
    def fastq(self):
        """
        The raw reads FASTQ object with convenience methods.
        See https://github.com/xapple/fasta#usage
        """
        # This class is taken from the `fasta` python package #
        from fasta import FASTQ
        fastq = FASTQ(self.path)
        # Change the location of the first FastQC, as we don't want to touch
        # the directory where the original reads are stored on the file system.
        # We might simply not have permission to write there.
        from fasta.fastqc import FastQC
        fastq.fastqc = FastQC(fastq, self.autopaths.fastqc_dir)
        # Change the location of the length distribution graphs too #
        fastq.graphs.length_dist.path = self.autopaths.len_dist_pdf
        fastq.graphs.length_hist.path = self.autopaths.len_hist_pdf
        # Return #
        return fastq

    @property_cached
    def primers(self):
        """
        Will return an object that holds the two primers of a
        given sample and has many convenience methods to parse and
        find the location of a primer inside all sequences.
        """
        from fasta.primers import TwoPrimers
        return TwoPrimers(self.fwd_primer_seq, self.rev_primer_seq)

    @property_cached
    def filter(self):
        """Takes care of filtering out unwanted sequences."""
        # Create filter object #
        from pacmill.filtering.seq_filter import SeqFilter
        seq_filter = SeqFilter(self)
        # Set parameters #
        seq_filter.primer_mismatches = self.primer_mismatches
        seq_filter.primer_max_dist   = self.primer_max_dist
        seq_filter.min_read_len      = self.min_read_len
        seq_filter.max_read_len      = self.max_read_len
        seq_filter.phred_window_size = self.phred_window_size
        seq_filter.phred_threshold   = self.phred_threshold
        # Return #
        return seq_filter

    @property_cached
    def chimeras(self):
        """Takes care of removing chimeric reads."""
        # Get file paths #
        source   = self.filter.results.clean
        cleaned  = self.autopaths.chimeras_cleaned
        rejects  = self.autopaths.chimeras_rejects
        # Create chimeras object #
        from pacmill.filtering.chimeras import Chimeras
        chimeras = Chimeras(source, cleaned, rejects)
        # Return #
        return chimeras

    @property_cached
    def barrnap(self):
        """Takes care of running the Barrnap program."""
        # In case the barrnap mode was turned off #
        if self.barrnap_mode == 'off': return None
        # Get file paths #
        source   = self.chimeras.results
        dest     = self.autopaths.barrnap_gff
        filtered = self.autopaths.barrnap_fasta
        # Create barrnap object depending on mode chosen #
        if self.barrnap_mode == 'filter':
            from pacmill.filtering.barrnap import BarrnapFilter
            barrnap = BarrnapFilter(source, dest, filtered)
        if self.barrnap_mode == 'concat':
            from pacmill.filtering.barrnap import RemoveITS
            barrnap = RemoveITS(source, dest, filtered)
        if self.barrnap_mode == 'trim':
            from pacmill.filtering.barrnap import BarrnapExtract
            barrnap = BarrnapExtract(source, dest, filtered)
        # Return #
        return barrnap

    @property_cached
    def final(self):
        """
        The final FASTQ for the sample before aggregation with other samples
        in a `Project` object. This is typically the output of the last step
        of quality control.
        """
        # If barrnap was the last step #
        if self.barrnap_mode != 'off': return self.barrnap.results
        # Otherwise chimeras was the last step #
        return self.chimeras.results

    @property_cached
    def report(self):
        """
        The SampleReport object is used to produce a PDF document containing
        all the information and graphs that concern this sample
        See https://github.com/xapple/pymarktex#usage
        """
        from pacmill.reports.sample import SampleReport
        return SampleReport(self, self.autopaths.report_pdf)
