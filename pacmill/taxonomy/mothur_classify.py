#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, multiprocessing, re

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

    This is the method that is implemented by the RDP and is described by
    Wang et al.
    """

    # Constants #
    short_name = 'mothur_classify'
    long_name  = 'Mothur Version 1.42.1'

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, database, dest_dir):
        # Source is the FASTA file containing OTU consensus sequences #
        self.source = FASTA(source)
        # Database is an object from the `seqsearch.databases` subpackage #
        self.database = database
        # Destination is a directory that contains all the results #
        self.dest_dir = DirectoryPath(dest_dir)

    #----------------------------- Installing --------------------------------#
    def check_installed(self):
        """
        Try to determine if the database chosen is downloaded and
        accessible.
        """
        if not self.database:
            msg  = "The taxonomic database does not seem to be accessible."
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
                /taxa_tables/
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
            message = "Running `%s` taxonomy classification on '%s'"
            print(message % (self.database.short_name, self.source))
        # Check the database is there #
        self.check_installed()
        # Check mothur is installed #
        check_cmd('mothur', True)
        # Remove the destination directory #
        self.dest_dir.remove()
        # Link our input OTU sequences to the destination directory #
        self.autopaths.centers.link_from(self.source)
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
            msg = "You can't access results from the mothur '%s' taxonomic" \
                  " classification before running the tool."
            raise Exception(msg % self.database.short_name)
        # Return the results #
        return MothurClassifyResults(self.autopaths)

###############################################################################
class MothurClassifyResults:

    def __init__(self, autopaths):
        self.autopaths = autopaths

    @property_cached
    def assignments(self):
        """
        Will parse the assignment file generated by mothur.
        Typically this attribute will return something like:

             'OTU-3576': ('Bacteria',
                          'Bacteria_unclassified',
                          'Bacteria_unclassified',
                          'Bacteria_unclassified',
                          'Bacteria_unclassified',
                          'Bacteria_unclassified',
                          '')
             'OTU-3577': ('Bacteria',
                          'Firmicutes',
                          'Clostridia',
                          'Clostridiales',
                          'Clostridiaceae_1',
                          'Clostridiaceae_1_unclassified',
                          '')
             'OTU-4086': ('unknown',
                          'unknown_unclassified',
                          'unknown_unclassified',
                          'unknown_unclassified',
                          'unknown_unclassified',
                          '')
        """
        # Initialize #
        result = {}
        with open(self.autopaths.assignments, 'r') as handle:
            for line in handle:
                # Get the OTU name and the assignment as a string #
                otu_name, species = line.split('\t')
                # Split the assignment into ranks #
                species = [i.strip('\n') for i in species.split(';')]
                # The OTU name matches what vsearch outputted in the FASTA #
                # And not what vsearch outputted in the TSV. We fix that #
                pattern  = '\Acentroid=(.+);seqs=[0-9]+\Z'
                otu_name = re.findall(pattern, otu_name)[0]
                # Assign that OTU to that tuple of ranks #
                result[otu_name] = tuple(species)
        # Return #
        return result

    @property_cached
    def count_unassigned(self):
        """
        Will count how many did not get a prediction at each level.
        NB: The greengenes results are one item longer than when using the
        silva or rdp databases, as they go until the species rank.
        """
        # Load the assignment values created by mothur #
        vals = list(self.assignments.values())
        # Calculate for each position in the tree of life #
        return [
            sum((1 for x in vals if x[0] == 'unknown')),      # Domain
            sum((1 for x in vals if 'unclassified' in x[1])), # Phylum
            sum((1 for x in vals if 'unclassified' in x[2])), # Class
            sum((1 for x in vals if 'unclassified' in x[3])), # Order
            sum((1 for x in vals if 'unclassified' in x[4])), # Family
            sum((1 for x in vals if 'unclassified' in x[5])), # Genus
            sum((1 for x in vals if 'unclassified' in x[6] or # Species
                                                      x[6] == ''))]

    @property_cached
    def count_assigned(self):
        """Will count how many did get a prediction at each level."""
        return [len(self.assignments) - x for x in self.count_unassigned]