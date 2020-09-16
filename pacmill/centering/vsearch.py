#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import multiprocessing

# First party modules #
from fasta import FASTA, FASTQ
from plumbing.check_cmd_found import check_cmd
from plumbing.cache           import property_cached
from autopaths.file_path      import FilePath

# Third party modules #
import sh

###############################################################################
class ClusterVsearch:
    """
    Takes care of clustering reads into OTUs by calling
    `vsearch --cluster_size`.

    From the man pages:

        Output cluster consensus sequences to filename. For each cluster,
        a multiple alignment is computed, and a consensus sequence is
        constructed by taking the majority symbol (nucleotide or gap) from
        each column of the alignment.

    Caveat: The Vsearch algorithm will name OTU differently in the FASTA
            file it outputs as compared to the TSV file it outputs

    In the `self.table` which is a TSV they are named:

        sample_1:1034      sample_1:1041       etc.

    In the `self.otus` which is a FASTA they are named:

       centroid=sample_1:1034;seqs=128      etc.

    After mothur classifies the OTU, they are named:

       centroid=sample_1_1034;seqs=128      etc.

    Which is almost the same except colons are underscores now.
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, threshold, min_size, otus=None, table=None):
        # Source is a FASTA to run the algorithm on #
        self.source = FASTQ(source)
        # Threshold is a float such as 0.97 #
        self.threshold = threshold
        # Minimum size is an integer such as 2 #
        self.min_size = min_size #TODO
        # self.otus is a FASTA file that contains the centroid sequences #
        if otus is None:
            otus = self.source.prefix_path + '.otus.fasta'
        self.otus = FASTA(otus)
        # self.table is a TSV file that contains the per sample counts #
        if table is None:
            table = self.source.prefix_path + '.otus.tsv'
        self.table = FilePath(table)

    #------------------------------ Running ----------------------------------#
    def __call__(self, cpus=None, verbose=True):
        # Message #
        if verbose:
            msg = "Running OTU creation on '%s'"
            print(msg % self.source)
        # Check it is installed #
        check_cmd('vsearch', True)
        # Check version #
        assert b"v2.14.1" in sh.vsearch('--version').stderr
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Pick the command parameters #
        command = ("--cluster_size", self.source,
                   "--consout",      self.otus,
                   "--id",           self.threshold,
                   "--otutabout",    self.table,
                   "--threads",      cpus)
        # Run the command on the input FASTA file #
        sh.vsearch(command)
        # Return #
        return self.otus

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the OTU creation software was run already and the
        results are stored on the filesystem. Return False if it was not yet
        run.
        """
        return self.otus.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from OTU creation " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return ClusteringResults(self.otus)

###############################################################################
class ClusteringResults(FASTA):
    """A file with the results."""
    pass