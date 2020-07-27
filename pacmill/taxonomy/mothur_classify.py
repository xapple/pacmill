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
from fasta import FASTA
from plumbing.check_cmd_found import check_cmd
from plumbing.cache           import property_cached
from autopaths.file_path      import FilePath
from autopaths.dir_path       import DirectoryPath

# Third party modules #
import sh

###############################################################################
class MothurClassify:
    """
    Predicts taxonomic composition by using a reference database and calling
    `classify.seqs` in mothur as seen in:

    https://github.com/novigit/broCode/blob/master/pbamp/estimateTaxComposition.sh#L31

        mothur "#classify.seqs(fasta=rnammer.$sample.prok.16S.both.fasta,
        reference=silva.nr_v128.align, taxonomy=silva.nr_v128.tax,
        relabund=T, processors=30)"
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest_dir):
        # Source is the FASTA file containing OTU consensus sequences #
        self.source = FASTA(source)
        # Destination is a directory that contains all the results #
        self.dest_dir = DirectoryPath(dest_dir)

    #---------------------------- The database -------------------------------#
    @property
    def database(self):
        """A link to the `SilvaMothur` database object for convenience."""
        from seqsearch.databases.silva_mothur import silva_mothur
        return silva_mothur

    #----------------------------- Installing --------------------------------#
    def check_installed(self):
        """
        Try to determine if the Silva database is downloaded and
        accessible.
        """
        if not self.database:
            msg  = "The silva database does not seem to be accessible."
            msg += " More information follows.\n\n"
            msg += self.database.__doc__
            raise Exception(msg)

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /centers.fasta
                /assignments.txt
                /summary.txt
                /flipped.txt
                /stdout.txt
                /stderr.txt
                """

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.dest_dir, self.all_paths)

    #------------------------------ Running ----------------------------------#
    def __call__(self, cpus=None, verbose=True):
        # Message #
        if verbose:
            print("Running taxonomy classification on '%s'" % self.source)
        # Check the database is there #
        self.check_installed()
        # Check mothur is installed #
        check_cmd('mothur', True)
        # Remove the destination directory #
        self.dest_dir.remove()
        # Link our input OTU sequences to the destination directory #
        self.autopaths.centers.link_from(self.source)
        # Check version #
        assert "1.42.1" in sh.mothur('--version')
        # Make the long command as multiple strings #
        cmd = ("#classify.seqs(", # The command
               " fasta=%s,"     , # The input file
               " reference=%s," , # The database (was also called 'template')
               " taxonomy=%s,"  , # The taxonomy file
               " processors=%s,", # The number of threads
               " probs=F);")      # Disable the output of bootstrap values.
        # Turn it into a single string #
        cmd = ''.join(cmd)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Format the string to contain all parameters #
        cmd = cmd % (self.autopaths.centers,
                     self.database.alignment,
                     self.database.taxonomy,
                     cpus)
        # Mothur pollutes with so much files, have to cwd #
        current_dir = os.getcwd()
        os.chdir(self.dest_dir)
        # Run the command on the input FASTA file #
        sh.mothur(cmd,
                  _out = self.autopaths.stdout.path,
                  _err = self.autopaths.stderr.path)
        # Restore current directory #
        os.chdir(current_dir)
        # Check output #
        if "ERROR" in self.autopaths.stdout.contents:
            msg  = "The Mothur classify operation did not run correctly."
            msg += "Please check the log files."
            raise Exception(msg)
        # Outputs #
        asgnmts = self.dest_dir + "centers.%s.wang.taxonomy"
        summary = self.dest_dir + "centers.%s.wang.tax.summary"
        flipped = self.dest_dir + "centers.%s.wang.flip.accnos"
        # Make into FilePath objects #
        asgnmts = FilePath(asgnmts % self.database.nickname)
        summary = FilePath(summary % self.database.nickname)
        flipped = FilePath(flipped % self.database.nickname)
        # Rename files #
        asgnmts.move_to(self.autopaths.assignments)
        summary.move_to(self.autopaths.summary)
        if flipped: flipped.move_to(self.autopaths.flipped)
        # Return #
        return self.dest_dir

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the taxonomic software was run already and the results
        are stored on the filesystem. Return False if it was not yet run.
        """
        return self.autopaths.assignments.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from taxonomic classification " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return MothurClassifyResults(self.autopaths)

###############################################################################
class MothurClassifyResults:
    def __init__(self, autopaths):
        self.autopaths = autopaths
