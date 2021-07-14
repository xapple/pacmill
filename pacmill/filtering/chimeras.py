#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, multiprocessing

# First party modules #
from fasta import FASTA, FASTQ
from plumbing.check_cmd_found import check_cmd
from plumbing.cache           import property_cached

# Third party modules #
import sh

###############################################################################
class Chimeras:
    """
    Takes care of detecting and removing chimeric reads, by calling
    `vsearch --uchime3_denovo`.
    Number of chunks cannot be set in the vsearch implementation.
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, cleaned=None, rejects=None):
        # Source is the FASTA or FASTQ file to run the algorithm on #
        self.source = FASTQ(source)
        # Cleaned is a FASTA file that contains the good reads #
        if cleaned is None:
            cleaned = self.source.prefix_path + '.no_chimeras.fasta'
        self.cleaned = FASTA(cleaned)
        # Rejects is a FASTA file that contains the bad chimeric reads #
        if rejects is None:
            rejects = self.source.prefix_path + '.yes_chimeras.fasta'
        self.rejects = FASTA(rejects)

    #------------------------------ Running ----------------------------------#
    def __call__(self, cpus=None, verbose=True):
        # Message #
        if verbose: print("Running chimeras detection on '%s'" % self.source)
        # Check it is installed #
        check_cmd('vsearch', True)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Pick the command parameters #
        command = ("--uchime3_denovo", self.source,
                   "-chimeras",        self.rejects,
                   "-nonchimeras",     self.cleaned,
                   "--threads",        cpus,
                   "-abskew",          1)
        # Run the command on the input FASTA file #
        sh.vsearch(command)
        # Sanity check the total #
        assert self.cleaned.count + self.rejects.count == self.source.count
        # Return #
        return self.cleaned

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the Chimeras software was run already and the results are
        stored on the filesystem. Return False if it was not yet run.
        """
        return self.cleaned.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from Chimeras " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return ChimerasResults(self.cleaned)

###############################################################################
class ChimerasResults(FASTA):
    """A file with the results from Chimeras."""
    pass