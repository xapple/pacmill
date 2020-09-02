#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
from collections import Counter

# Internal modules #

# First party modules #
from fasta import FASTQ
from plumbing.cache import property_cached

# Third party modules #
import pandas

###############################################################################
class SeqFilter:
    """
    This class takes care of filtering sequences from a FASTA or FASTQ file
    depending on certain criteria.

    - Filter primers:
       * Check that the primers are found where they should be found.
       * Check that the primers have the sequence they should have.

    - Filter based on presence of N bases.

    - Filter based on maximum and minimum sequence lengths.

    - Filter based on the minimum PHRED score of a sliding window.

    You can adjust several parameters:

    * primer_mismatches: Number of mismatches allowed before discarding.
    * primer_max_dist:   Maximum distance for primer presence before
                         discarding (counted in base pairs from sequence
                         start or end).
    * min_read_length:   Minimum sequence length.
    * max_read_length:   Maximum sequence length.
    * phred_window_size: The size of the rolling window (in base pairs) inside
                         which scores are averaged.
    * phred_threshold:   The quality score cutoff point for the average score
                         within any window.

    The final cleaned FASTQ is available at `self.results.clean`
    """

    # Attributes #
    short_name = 'seq_filter'

    # Default parameters #
    primer_mismatches = None
    primer_max_dist   = None
    min_read_len      = None
    max_read_len      = None
    phred_window_size = None
    phred_threshold   = None

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.sample)

    def __init__(self, sample):
        # Save the reference to a sample object #
        self.sample = sample
        # The different files #
        self.primers_fastq = FASTQ(self.autopaths.primers)
        self.n_base_fastq  = FASTQ(self.autopaths.n_base)
        self.length_fastq  = FASTQ(self.autopaths.length)
        self.score_fastq   = FASTQ(self.autopaths.score)
        self.renamed_fastq = FASTQ(self.autopaths.renamed)
        # The final result #
        self.clean = self.renamed_fastq

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /primers.fastq
                /n_base.fastq
                /length.fastq
                /score.fastq
                /renamed.fastq
                """

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.sample.autopaths.filtered_dir, self.all_paths)

    #-------------------------------- Running --------------------------------#
    def __call__(self, verbose=False):
        # Message #
        if verbose: print("Filtering sample '%s'" % self.sample.short_name)
        # Primers #
        self.primer_filter()
        # N bases #
        self.n_base_filter()
        # Length #
        self.len_filter()
        # Score #
        self.score_filter()
        # Rename with a number #
        self.score_fastq.rename_with_num(self.sample.short_name + ':',
                                         self.renamed_fastq)
        # Check #
        if len(self.score_fastq) == 0:
            msg = "No results left after filtering the sample '%s'."
            raise Exception(msg % self.sample.short_name)
        # Return #
        return self.results.clean

    #-------------------------------- Primers --------------------------------#
    def primer_gen(self, reads, verbose=False, debug=False):
        """
        We will uses regex patterns to search every read for both primers.
        We will record the start and end positions of primers when they are
        found. Both the forward and reverse primers are searched for.
        Both the original sequences and their reverse complements are
        searched for, in case the read is on the opposite strand.
        """
        # Select verbosity #
        import tqdm
        wrapper = tqdm.tqdm if verbose else lambda x: x
        # Shorter name for the distance allowed #
        dist = self.primer_max_dist
        # Loop #
        for r in wrapper(reads):
            # Use this for debugging purposes #
            if debug: print(r.pretty_visualization)
            # Did we find both primers? #
            fwd_found = r.fwd_srt is not None or r.fwd_rc_srt is not None
            rev_found = r.rev_srt is not None or r.rev_rc_srt is not None
            # Skip reads that don't pass this criteria #
            if not fwd_found or not rev_found: continue
            # These situations should not occur, could be long chimeras #
            if r.fwd_srt is not None    and r.rev_srt is not None:    continue
            if r.fwd_rc_srt is not None and r.rev_rc_srt is not None: continue
            # We are in a forward sequence situation #
            if r.fwd_srt is not None and r.rev_rc_srt is not None:
                if dist and r.fwd_srt                 > dist: continue
                if dist and len(r.seq) - r.rev_rc_end > dist: continue
                out = r.read[r.fwd_end:r.rev_rc_srt]
                if len(out) == 0: continue
            # We are in a reverse complemented situation #
            if r.fwd_rc_srt is not None and r.rev_srt is not None:
                if dist and r.rev_srt                 > dist: continue
                if dist and len(r.seq) - r.fwd_rc_end > dist: continue
                out = r.read[r.rev_end:r.fwd_rc_srt].reverse_complement()
                if len(out) == 0: continue
            # Return #
            yield out

    def primer_filter(self):
        """
        Will take only reads that have both primers in the correct position
        and with the correct sequence. Will also trim the primers once found.
        """
        # This will return a generator object with a length property #
        all_reads = self.sample.fastq.parse_primers(self.sample.primers,
                                                    self.primer_mismatches)
        # Write #
        self.primers_fastq.write(self.primer_gen(all_reads))
        # Close #
        self.sample.fastq.close()

    #------------------------------ N bases ----------------------------------#
    def n_base_gen(self, reads):
        for r in reads:
            if 'N' in r: continue
            yield r

    def n_base_filter(self):
        self.n_base_fastq.write(self.n_base_gen(self.primers_fastq))

    #------------------------------- Length ----------------------------------#
    def len_gen(self, reads, verbose=False):
        for r in reads:
            if self.min_read_len > 0:
                if len(r.seq) < self.min_read_len:
                    if verbose: print("Discard")
                    continue
            if self.max_read_len > 0:
                if len(r.seq) > self.max_read_len:
                    if verbose: print("Discard")
                    continue
            if verbose: print("Keep")
            yield r

    def len_filter(self):
        # Optionally bypass this step #
        if self.min_read_len is None and self.max_read_len is None:
            self.length_fastq.copy(self.n_base_fastq)
        # Perform the length filtering #
        self.length_fastq.write(self.len_gen(self.n_base_fastq))

    #-------------------------------- Score ----------------------------------#
    def score_gen(self, reads):
        # Parameters #
        window    = self.phred_window_size
        threshold = self.phred_threshold
        # Loop #
        for r in reads:
            quality = pandas.Series(r.letter_annotations['phred_quality'])
            rolling = quality.rolling(window).mean().dropna().tolist()
            if any(score < threshold for score in rolling): continue
            yield r

    def score_filter(self):
        # Optionally bypass this step #
        if self.phred_window_size is None or self.phred_threshold is None:
            self.score_fastq.copy(self.length_fastq)
        # Perform the score filtering #
        self.score_fastq.write(self.score_gen(self.length_fastq))

    #------------------------------ Debugging --------------------------------#
    @property_cached
    def primer_positions(self):
        """Useful for diagnostics. Returns the primer positions counts."""
        # Count positions #
        all_fwd_pos,    all_rev_pos    = Counter(), Counter()
        all_fwd_rc_pos, all_rev_rc_pos = Counter(), Counter()
        # Make a generator #
        reads = self.sample.fastq.parse_primers(self.sample.primers,
                                                self.primer_mismatches)
        # Select verbosity #
        verbose = True
        import tqdm
        wrapper = tqdm.tqdm if verbose else lambda x: x
        # Iterate over the generator #
        for r in wrapper(reads):
            # Did we find both primers? #
            fwd_found = r.fwd_srt is not None or r.fwd_rc_srt is not None
            rev_found = r.rev_srt is not None or r.rev_rc_srt is not None
            # Skip reads that don't pass this criteria #
            if not fwd_found or not rev_found: continue
            # Use this for debugging purposes #
            if False: print(r.pretty_visualization)
            # These situations should not occur, could be long chimeras #
            if r.fwd_srt is not None    and r.rev_srt is not None:    continue
            if r.fwd_rc_srt is not None and r.rev_rc_srt is not None: continue
            # Add counts #
            if r.fwd_srt is not None:
                all_fwd_pos.update((r.fwd_srt,))
            if r.rev_srt is not None:
                all_rev_pos.update((r.rev_srt,))
            if r.fwd_rc_srt is not None:
                all_fwd_rc_pos.update((r.fwd_rc_srt,))
            if r.rev_rc_srt is not None:
                all_rev_rc_pos.update((r.rev_rc_srt,))
        # Return results #
        return all_fwd_pos, all_rev_pos, all_fwd_rc_pos, all_rev_rc_pos

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the filtering was run already and the results are
        stored on the filesystem. Return False if it was not yet run.
        """
        return bool(self.renamed_fastq)

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
class SeqFilterResults:

    def __init__(self, parent):
        # Reference to the SeqFilter instance #
        self.parent = parent
        # The intermediary steps #
        self.primers_fastq = parent.primers_fastq
        self.n_base_fastq  = parent.n_base_fastq
        self.length_fastq  = parent.length_fastq
        self.score_fastq   = parent.score_fastq
        self.renamed_fastq = parent.renamed_fastq
        # The final result #
        self.clean = self.renamed_fastq