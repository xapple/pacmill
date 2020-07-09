#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Third party modules #

# First party modules #
from plumbing.cache import property_cached
from autopaths      import Path

# Internal modules #

###############################################################################
class Sample:
    """
    A Sample object consists first and foremost of one FASTQ file.

    For a pipeline that uses paired reads (and hence two FASTQ files) please
    have a look at a sister project: www.github.com/xapple/sifes

    In addition, a Sample has many other metadata associated with it.
    Every column of the original excel file is now the name of an attribute.
    For instance (non exhaustive list):

        * self.num: the number of this sample in the project.
        * self.short_name: the name of this sample.
        * self.long_name: the name of this sample.
        * self.grouping: optional grouping information.

    In addition, the Sample object has the following:

        * self.proj: a reference to a Project instance.

    Other properties are described in their respective docstrings.
    """

    def __repr__(self):
        return '%s object code "%s"' % (self.__class__, self.short_name)

    def __init__(self, parent, **kwargs):
        """
        A Sample object takes a dictionary of metadata as input.
        Each named parameter will be set as an attribute of this instance
        with the same name. Hence, the attributes of a Sample correspond
        directly to the column headers of the excel metadata file that
        was parsed by the Project object.
        """
        # Keep a reference to the Project object #
        self.parent  = parent
        self.project = parent
        # Record which metadata keys were passed to this sample #
        self.metadata_keys = [k for k in kwargs if not k.startswith("Unnamed")]
        # Set the attributes of this instance with the given kwargs #
        for name, value in kwargs.items(): setattr(self, name, value)
        # We need to transform some attributes #
        self.transform_attrs()
        # We need to validate some attributes #
        self.validate_attrs()

    #------------------------------- Methods ---------------------------------#
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

    def validate_attrs(self):
        """
        This is where we check that the information in the excel
        metadata file is consistent and usable.
        """
        # Check the short name contains only alphanumerics and underscore #
        assert self.short_name.isidentifier()
        # Check that all directories always end with a dash (/) #
        for attribute in self.__dict__:
            if not attribute.endswith('_dir'): continue
            value = getattr(self, attribute)
            if not isinstance(value, str) or not value.endswith("/"):
                msg = "The `%s` entry of <%s> must end with a slash. " \
                          "It currently does not: '%s'"
                msg = msg % (attribute, self.description, value)
                raise ValueError(msg)
        # Check that the FASTQ file exists #
        if not self.path.exists:
            msg = "The FASTQ file path of <%s> is not found. " \
                  "It should be located at: '%s'"
            msg = msg % (self.description, self.path)
            raise FileNotFoundError(msg)
        # Check the DNA sequences #
        dna_keys = ['fwd_index_seq', 'rev_index_seq', 'fwd_primer_seq',
                    'rev_primer_seq']
        for key in dna_keys:
            seq = getattr(self, key)
            if '\n' in seq:
                msg = "The `%s` entry of <%s> contains illegal " \
                      "characters, please check: '%s'"
                msg = msg % (key, self.description, seq)
                raise ValueError(msg)

    #----------------------------- Properties --------------------------------#
    @property
    def description(self):
        """A string describing the current sample."""
        desc = "Sample object %i in project '%s' with name '%s'"
        return desc % (self.num, self.project.short_name, self.short_name)

    @property_cached
    def path(self):
        """The path to the FASTQ reads file."""
        # Join the three components together #
        return Path(self.input_dir + self.suffix_dir + self.fwd_file_name)

    #------------------------- Automatic paths -------------------------------#
    all_paths = """
                /fastqc/
                /graphs/raw_len_dist.pdf
                /graphs/raw_len_hist.pdf
                /filtered/
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
        The FASTQ object with convenience methods.
        See https://github.com/xapple/fasta#usage
        """
        # This class is taken from the `fasta` python package #
        from fasta import FASTQ
        fastq = FASTQ(self.path)
        # Change the location of the first FastQC, as we don't want to touch
        # the directory where the original reads are stored on the file system.
        from fasta.fastqc import FastQC
        fastq.fastqc = FastQC(fastq, self.autopaths.fastqc_dir)
        # Change the location of the length distribution graphs too.
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
        """Will filter out unwanted sequences."""
        from pacmill.filtering.seq_filter import SeqFilter
        return SeqFilter(self.fastq, self.autopaths.filtered_dir,
                         self.short_name, self.primers)

    @property_cached
    def report(self):
        """
        The SampleReport object is used to produce a PDF document containing
        all the informations and graphs that concern this sample
        See https://github.com/xapple/pymarktex#usage
        """
        from pacmill.reports.sample import SampleReport
        return SampleReport(self, self.autopaths.report_pdf)
