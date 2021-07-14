#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os

# First party modules #
from fasta import FASTQ
from plumbing.check_cmd_found import check_cmd
from plumbing.cache           import property_cached
from autopaths.tmp_path       import new_temp_dir, new_temp_path

# Third party modules #
import sh

###############################################################################
class MothurUchime:
    """
    Takes care of detecting and removing chimeric reads, by calling
    `mothur.uchime` as seen in:

    https://github.com/novigit/broCode/blob/master/pbamp/readCurationPipeline.sh#L103

    mothur "#chimera.uchime(fasta=roi.$sample.rhq.trim.fwdrev.pol.fasta,
                            reference=self, chunks=16, abskew=1, chimealns=T)"
    mothur "#remove.seqs(fasta=roi.$sample.rhq.trim.fwdrev.pol.fasta,
                         accnos=roi.$sample.rhq.trim.fwdrev.pol.denovo.uchime.accnos)"

    In case of missing output see:

    https://forum.mothur.org/t/chimera-uchime-does-not-produce-any-output/20621
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest=None):
        # Source is the FASTA or FASTQ file #
        self.source = FASTQ(source)
        # Destination is a FASTQ file that contains the results #
        if dest is None:
            dest = self.source.prefix_path + '.chimeras.fastq'
        self.dest = FASTQ(dest)

    #------------------------------ Running ----------------------------------#
    def __call__(self, verbose=True):
        # Message #
        if verbose: print("Running chimeras detection on '%s'" % self.source)
        # Check it is installed #
        check_cmd('mothur', True)
        # Make new temporary directory #
        tmp_dir = new_temp_dir()
        source = tmp_dir + 'reads.fasta'
        # If the input is a FASTQ we need to make a FASTA first #
        if self.source.endswith('fastq'):
            FASTQ(self.source).to_fasta(source)
        else:
            self.source.copy_to(source)
        # Make the long command as a string #
        command = "#chimera.uchime("    \
                  "fasta=%s,"           \
                  "reference=self,"     \
                  "chunks=16,"          \
                  "chimealns=T,"        \
                  "abskew=1)"
        # Mothur pollutes with so much files, have to cwd #
        current_dir = os.getcwd()
        os.chdir(tmp_dir)
        # Run the command on the input FASTA file #
        sh.mothur(command % source)
        # Restore current directory #
        os.chdir(current_dir)
        # Move files #
        pass #TODO
        # Return #
        return self.dest

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the chimeras software was run already and the results are
        stored on the filesystem. Return False if it was not yet run.
        """
        return self.dest.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from chimeras " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return ChimerasResults(self.dest)

###############################################################################
class ChimerasResults(FASTQ):
    """A file with the results from chimeras."""
    pass