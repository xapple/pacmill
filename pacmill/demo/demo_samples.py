#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

This file is used to download samples that are used in the demonstration
project of the `pacmill` pipeline.

It describes 5 samples that were all part of the following study:

    ``Confident phylogenetic identification of uncultured prokaryotes through
    long read amplicon sequencing of the 16S-ITS-23S rRNA operon.''
    - Joran Martijn, [many others], Thijs Ettema.
    @Science for Life Laboratory, Uppsala University
    in Environmental microbiology, 2019
    https://doi.org/10.1111/1462-2920.14636

The publication is free to download and the reads are available through SRA
as a BioProject at the following address:

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA498591

The five samples used in the demo project are described as follows:

    * The mock community *
    Genomic DNA from 38 phylogenetically distinct and diverse bacteria
    and archaea.

    * the P19 sample *
    Sediment sample obtained from hot spring Radiata Pool, Ngatamariki,
    New Zealand.

    * the PM3 sample *
    Sediment sample taken from 1.25m below the sea floor using a gravity core
    at Aarhus Bay, Denmark.

    * the SALA sample *
    Black biofilm that was taken at 60m depth in an old silver mine near
    Sala, Sweden.

    * the TNS08 sample *
    Sediment sample taken from a shallow submarine hydrothermal vent field
    near Taketomi Island, Japan.
"""

# Built-in modules #

# First-party modules #
from autopaths.file_path import FilePath
from plumbing.cache import property_cached
from plumbing.scraping import download_from_url

# Internal modules #
from pacmill import repos_dir

###############################################################################
class DemoSample:
    """
    Typically you can use this class like such:

        >>> from pacmill.demo.demo_samples import samples
        >>> for sample in samples: print(sample.description)
        >>> for sample in samples: print(sample.path_fastq)
        >>> for sample in samples: print(sample.download())
        >>> for sample in samples: print(sample.extract())
        >>> for sample in samples: print(sample.short_name, sample.fastq.count)
    """

    # Default location for saving samples #
    base_dir = repos_dir + 'demo_project/demo_sequence_data/'

    def __repr__(self):
        return '%s object code "%s"' % (self.__class__, self.short_name)

    def __init__(self, short_name, num, url):
        """
        A DemoSample object takes an SRA download URL as input
        along with a sample short name and a number.
        """
        self.short_name = short_name
        self.url = url
        self.num = num

    #------------------------------- Methods ---------------------------------#
    def __call__(self):
        """Will download and extract this sample."""
        self.download()
        self.extract()

    def download(self):
        """Will download the SRA file from the internet."""
        return download_from_url(self.url, self.path_sra,
                                 stream=True, progress=True)

    def extract(self):
        """Will extract the FASTQ file from the SRA file."""
        return self.sra_dump()

    #----------------------------- Properties --------------------------------#
    @property
    def description(self):
        """A string describing the current sample."""
        desc = "Demo sample object %i with name '%s'"
        return desc % (self.num, self.short_name)

    @property_cached
    def path_fastq(self):
        """The path to the raw FASTQ reads file."""
        return FilePath(self.base_dir + self.short_name + '.fastq')

    @property_cached
    def path_sra(self):
        """The path to the raw SRA archive file."""
        return FilePath(self.base_dir + self.short_name + '.sra')

    #---------------------------- Compositions -------------------------------#
    @property_cached
    def sra_dump(self):
        """Takes care of running `fastq-dump` to extract sequences."""
        # Create filter object #
        from pacmill.demo.sra import DumpSRA
        sra_dump = DumpSRA(self.path_sra, self.path_fastq)
        # Return #
        return sra_dump

    @property_cached
    def fastq(self):
        """
        The raw reads FASTQ object with convenience methods.
        See https://github.com/xapple/fasta#usage
        """
        # This class is taken from the `fasta` python package #
        from fasta import FASTQ
        fastq = FASTQ(self.path_fastq)
        # Return #
        return fastq

###############################################################################
# Hardcode the download links #
samples = {
  'mock':   'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8113901/SRR8113901',
  'p_19':   'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8113902/SRR8113902',
  'pm_3':   'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8113898/SRR8113898',
  'sala':   'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8113900/SRR8113900',
  'tns_08': 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8113899/SRR8113899',
}

# Make DemoSample objects #
samples = [DemoSample(k, i+1, v) for i, (k, v) in enumerate(samples.items())]