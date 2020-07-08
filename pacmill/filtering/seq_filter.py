#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import re
from collections import Counter

# Internal modules #

# First party modules #
from fasta                import FASTA
from autopaths.auto_paths import AutoPaths
from plumbing.cache       import property_cached

###############################################################################
class SeqFilter(object):
    """
    This class takes care of filtering sequences from a FASTA or FASTQ file
    depending on certain criteria.

    - Filter primers:
       * Check that the primers are found where they should be found.
       * Check that the primers have the sequence they should have.

    - Filter based on N bases.

    - Filter based on maximum and minimum sequence lengths.

    You can adjust several parameters:

    * primer_mismatches: Number of mismatches allowed before discarding.
    * primer_max_dist:   Maximum distance for primer presence before discarding.
                         (Counted in base pairs from sequence start)
    * min_read_length:   Minimum sequence length
    * max_read_length:   Maximum sequence length
    """

    # Attributes #
    short_name = 'seq_filter'

    # Parameters #
    primer_mismatches = 2
    primer_max_dist   = 70
    min_read_length   = None
    max_read_length   = None

    # Paths #
    all_paths = """
                /primers.fasta
                /n_base.fasta
                /length.fasta
                /region.fasta
                /renamed.fasta
                """

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.input)

    def __init__(self, input, result_dir, sample_name, primers):
        # Save attributes #
        self.input       = input
        self.result_dir  = result_dir
        self.sample_name = sample_name
        self.primers     = primers
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # The different files #
        self.primers_fasta = FASTA(self.p.primers)
        self.n_base_fasta  = FASTA(self.p.n_base)
        self.length_fasta  = FASTA(self.p.length)
        self.renamed_fasta = FASTA(self.p.renamed)
        # The final result #
        self.clean = self.renamed_fasta

    def __call__(self, verbose=False):
        # Message #
        if verbose: print("Filtering sample '%s'" % self.sample_name)
        # Primers #
        self.primer_filter()
        # N bases #
        self.n_base_filter()
        # Length #
        self.len_filter()
        # Rename with a number #
        self.length_fasta.rename_with_num(self.sample_name + ':', self.renamed_fasta)
        # Check #
        if len(self.length_fasta) == 0:
            raise Exception("No results left after filtering the sample '%s'" % self.sample_name)
        # Return #
        return self.results.clean

    # Primers #
    def primer_filter(self):
        """Will take only reads that have both primers and will trim the primers."""
        def good_primer_iterator(reads):
            for r in reads:
                if r.fwd_start_pos is None or r.rev_start_pos is None: continue
                if r.fwd_start_pos > self.primer_max_dist:             continue
                if r.rev_start_pos < -self.primer_max_dist:            continue
                yield r.read[r.fwd_end_pos:r.rev_end_pos]
        # Useful function, returns a GenWithLength #
        parse_primers = lambda: self.input.parse_primers(self.primers, self.primer_mismatches, revcompl=True)
        # Do it #
        self.primers_fasta.write(good_primer_iterator(parse_primers()))

    # N base #
    def n_base_filter(self):
        def good_n_base_iterator(reads):
            for r in reads:
                if 'N' in r: continue
                yield r
        self.n_base_fasta.write(good_n_base_iterator(self.primers_fasta))

    # Length #
    def len_filter(self):
        def good_len_iterator(reads, verbose=False):
            for r in reads:
                if self.min_read_length > 0:
                    if len(r.seq) < self.min_read_length:
                        if verbose: print("Discard")
                        continue
                if self.max_read_length > 0:
                    if len(r.seq) > self.max_read_length:
                        if verbose: print("Discard")
                        continue
                if verbose: print("Keep")
                yield r
        self.length_fasta.write(good_len_iterator(self.n_base_fasta))

    #-------------------------------------------------------------------------#
    @property_cached
    def primer_positions(self):
        """Useful for diagnostics."""
        # Count positions #
        all_fwd_pos, all_rev_pos = Counter(), Counter()
        parse_primers = lambda: self.input.parse_primers(self.primers, self.primer_mismatches, revcompl=True)
        for r in parse_primers():
            if r.fwd_start_pos is not None: all_fwd_pos.update((r.fwd_start_pos,))
            if r.rev_start_pos is not None: all_rev_pos.update((r.rev_start_pos,))
        # Return results #
        return all_fwd_pos, all_rev_pos

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the filtering was run already and the results are
        stored on the filesystem. Return False if it was not yet run.
        """
        return bool(self.renamed_fasta)

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from the filtering procedure " \
                  "before running the algorithm."
            raise Exception(msg)
        # Return the results #
        return SeqFilterResults(self)

###############################################################################
class SeqFilterResults(object):

    def __init__(self, parent):
        # Reference to the SeqFilter instance #
        self.parent = parent
        # The intermediary steps #
        self.primers_fasta = parent.primers_fasta
        self.n_base_fasta  = parent.n_base_fasta
        self.length_fasta  = parent.length_fasta
        self.renamed_fasta = parent.renamed_fasta
        # The final result #
        self.clean = self.renamed_fasta