#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import multiprocessing, re, shutil

# First party modules #
from fasta import FASTA, FASTQ
from plumbing.check_cmd_found import check_cmd
from plumbing.cache           import property_cached
from autopaths.file_path      import FilePath

# Third party modules #
import sh, pandas

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
        self.min_size = min_size
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
        # Call vsearch #
        self.cluster(cpus, verbose)
        # Filter OTUs that are too small #
        if self.min_size > 1: self.dereplicate(verbose)
        # Return #
        return self.otus

    def cluster(self, cpus=None, verbose=True):
        """
        Will call the vsearch algorithm. Documentation for command line
        options are found at:
        https://manpages.debian.org/stretch/vsearch/vsearch.1.en.html
        """
        # Message #
        if verbose:
            msg = "Running OTU creation on '%s'"
            print(msg % self.source)
        # Check it is installed #
        check_cmd('vsearch', True)
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

    def dereplicate(self, verbose=True):
        """
        We will read the TSV table at `self.table` and remove OTUs that are
        below the self.min_size threshold. Following which, we will read the
        FASTA file at `self.otus` and remove the same sequences that are no
        longer required.
        """
        # Message #
        if verbose:
            msg = "Removing low abundance OTUs on '%s'"
            print(msg % self.otus)
        # Load table #
        df = pandas.read_csv(self.table, sep='\t', index_col=0)
        # Sum and sort #
        totals = df.sum(axis=1).sort_values(ascending=False)
        # Get the sequences IDs to drop and those to keep #
        self.keep_ids = list(totals[totals >= self.min_size].index)
        self.drop_ids = list(totals[totals <  self.min_size].index)
        # Filter the dataframe #
        df = df.loc[self.keep_ids]
        # Backup the old table and replace it #
        shutil.move(self.table, self.table + '.unfiltered')
        df.to_csv(self.table.path, sep='\t')
        # Function for filtering reads #
        def keep_reads_if(title):
            pattern  = r'\Acentroid=(.+);seqs=[0-9]+\Z'
            otu_name = re.findall(pattern, title)[0]
            if otu_name in self.keep_ids: return True
            return False
        # Filter the FASTA file #
        self.otus.copy(self.otus + '.unfiltered')
        self.otus.extract_sequences(keep_reads_if, in_place=True)
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