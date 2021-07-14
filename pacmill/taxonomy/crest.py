#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, re

# First party modules #
from fasta import FASTA
from autopaths.dir_path import DirectoryPath
from plumbing.cache     import property_cached

# Third party modules #

###############################################################################
class CrestClassify:
    """
    Predicts taxonomy of DNA sequences by BLASTING against the 'SILVAMOD138'
    database and picking the lowest common ancestor amongst contradicting
    assignments in the tree of life.

    For instance if all hits for a given OTU only agree at the order level,
    the assignment stops at the order level.

    CREST stands for Classification Resources for Environmental Sequence Tags.

    The code is here:        https://github.com/xapple/crest4
    The publication is here: https://dx.plos.org/10.1371/journal.pone.0049334

    We have lowered the "drop-off from the highest bitscore" to half a percent
    by adding the "--score_drop 0.5" parameter. This should increase the
    resolution for long reads such as the ones we have.
    """

    # Constants #
    short_name = 'crest4'
    long_name  = 'CREST version 4.1.8'

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest_dir, otu_tsv):
        # Source is the FASTA file containing OTU consensus sequences #
        self.source = FASTA(source)
        # Destination is a directory that contains all the results #
        self.dest_dir = DirectoryPath(dest_dir)
        # Path to a TSV file containing the OTU table #
        # TODO The IDs in this file are like 'centroid=tns_08:1023;seqs=138'
        # TODO While the ones in the FASTA are like 'tns_08:1023'
        self.otu_tsv = otu_tsv

    #------------------------------ Running ----------------------------------#
    @property_cached
    def crest_obj(self):
        # Import #
        from crest4 import Classify
        # Create a new instance #
        return Classify(fasta       = self.source,
                        search_algo = 'blast',
                        num_threads = True,
                        output_dir  = self.dest_dir,
                        search_hits = self.dest_dir + "search.hits",
                        score_drop  = 0.5,
                        )

    @property
    def database(self):
        """
        The crest object has a database property that has all the
        information we need.
        """
        return self.crest_obj.database

    def __call__(self, run_search=True, verbose=True):
        # Message #
        if verbose:
            message = "Running `%s` taxonomy classification on '%s'"
            print(message % (self.short_name, self.source))
        # We want to regenerate the search results, unless specified #
        if run_search:
            self.crest_obj.search_hits.remove()
        # Run the similarity search and classification #
        self.crest_obj()
        # Check #
        self.crest_obj.search_hits.must_exist()
        if os.path.getsize(self.crest_obj.search_hits) == 0:
            msg = "Hits file empty. The BLAST process was probably killed."
            raise Exception(msg)

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the taxonomic software was run already and the results
        are stored on the filesystem. Return False if it was not yet run.
        """
        return bool(self.crest_obj.out_file)

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from the CREST taxonomic" \
                  " classification before running the tool."
            raise Exception(msg)
        # Return the results #
        return CrestResults(self.crest_obj)

###############################################################################
class CrestResults:

    def __init__(self, crest_obj):
        self.crest_obj = crest_obj

    @property_cached
    def assignments(self):
        """
        Typically an assignment from CREST looks like this:

            Main genome,Bacteria,Bacteria (superkingdom),Terrabacteria,
            Firmicutes,Bacilli,Bacillales,Bacillaceae,Bacillus

        We skip the first entry as it seems to be always "Main genome"
        and classified as rank "meta".
        """
        result = {}
        with open(self.crest_obj.out_file, 'r') as handle:
            for line in handle:
                # Split by tabs #
                code, species = line.strip('\n').split('\t')
                # Split by semi-colons #
                species = species.split(';')
                # Remove useless whitespace #
                species = tuple(map(lambda s: s.strip(), species))
                # Get the original OTU name #
                pattern  = r'centroid=(.+);seqs=[0-9]+\Z'
                otu_name = re.findall(pattern, code)[0]
                # Conform to the mothur standard of names #
                otu_name = otu_name.replace(':', '_')
                # Make a dictionary #
                result[otu_name] = species
        # Check size #
        assert len(result) == len(self.crest_obj.fasta)
        # Return #
        return result

    @property_cached
    def count_unassigned(self):
        """Will count how many did not get a prediction at each level."""
        # Load the assignment values created by CREST #
        taxas = list(self.assignments.values())
        # Iterate over every rank number #
        rank_numbers = list(range(len(self.crest_obj.database.rank_names)))
        # Calculate for each position in the tree of life #
        result = [sum(1 for x in taxas if len(x) <= i) for i in rank_numbers]
        # Return #
        return result

    @property_cached
    def count_assigned(self):
        """Will count how many did get a prediction at each level."""
        return [len(self.assignments) - x for x in self.count_unassigned]