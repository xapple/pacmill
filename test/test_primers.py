#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Script to test the parsing and filtering of primers in sequence data.
"""

# Built-in modules #

# Internal modules #

# First party modules #
from fasta.primers import TwoPrimers
from fasta import FASTQ

# Third party modules #
import pytest, Bio

###############################################################################
cases = [{

    'fwd_primer': "ATTTA",                         # This is a standard case
    'rev_primer':                     "AGGGA",     # Only forward and reverse
    'seq':        "ATTTAGGGGGGGCGGGGGGGTCCCT",     # primers are found
    'expected':        "GGGGGGGCGGGGGGG",

}, {

    'fwd_primer': "AAA",                           # Here both are found but
    'rev_primer': "TTT",                           # they are found in the same
    'seq':        "AAAAAAAA",                      # spot causing a 0 length

}, {

    'fwd_primer':                     "AATAA",     # The sequence is in the
    'rev_primer': "AAAAA",                         # opposite direction
    'seq':        "AAAAAGGGGGGCGCGGGGGGTTATT",     # (reverse complement)
    'expected':        "CCCCCCGCGCCCCCC",

}, {

    'fwd_primer': "ATTTA",                         # Mismatch in the reverse
    'rev_primer':                     "AGGGA",     # primer
    'seq':        "ATTTAGGGGGGGCGGGGGGGTCTCT",

}, {

    'fwd_primer': "ATTTA",                        # Authorizing a mismatch here
    'rev_primer':                     "AGGGA",    # causes the reverse primer
    'seq':        "ATTTAGGGGGGGCGGGGGGGTCCCT",    # to be found twice at pos
    'mismatches': 1                               # 4 and at pos 20

}, {

    'fwd_primer': "ATTAATTA",
    'rev_primer':                    "AGAAGA",    # Authorizing a mismatch here
    'seq':        "TTTAATTAGGGGCGGGGGGTCTTTT",    # causes both primers to be
    'expected':           "GGGGCGGGGGG",          # found correctly
    'mismatches': 1

}, {

    'fwd_primer':    "ATTAATTA",
    'rev_primer':                       "AGAAGA",     # Primers are further
    'seq':        "GGGTTTAATTAGGGGCTGGGGGTCTTTTGGG",  # inside the read but
    'expected':              "GGGGCTGGGGG",           # are found correctly

    'mismatches': 1
}]

input_names = ('fwd_primer', 'rev_primer', 'seq', 'expected', 'mismatches')
params = [[case.get(name) for name in input_names] for case in cases]

###############################################################################
@pytest.mark.parametrize(input_names, params)
def test_primers(seq_filter, fwd_primer, rev_primer, seq, expected,
                 mismatches):
    # Create primers #
    primers = TwoPrimers(fwd_primer, rev_primer)
    # Create the sequences for the  fake FASTQ #
    reads = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id='pytest_read')]
    # Make a fake FASTQ #
    fastq = FASTQ('/dummy_path.fastq')
    # Set some attributes (monkey-patch) #
    fastq.parse = lambda: reads
    fastq.count = 1
    # Call the primer parsing routines #
    reads_w_primers = fastq.parse_primers(primers, mismatches=mismatches)
    # Get the result of the function #
    result = seq_filter.primer_gen(reads_w_primers)
    # Case where the sequence should be eliminated #
    if expected is None:
        with pytest.raises(StopIteration): next(result)
        return
    # Get the result #
    output = str(next(result).seq)
    # Assert #
    assert output == expected
