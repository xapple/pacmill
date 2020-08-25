#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Script to test the parsing and filtering of sequence after
running them through the barrnap algorithm.
"""

# Built-in modules #

# Internal modules #
from pacmill.filtering.barrnap import RemoveITS

# First party modules #
from fasta import FASTA

# Third party modules #

###############################################################################
def test_remove_its(this_script_dir):
    # Get file paths #
    test_fasta = this_script_dir + 'test_sequence.fasta'
    # Length of original sequence #
    orig_len = len(FASTA(test_fasta).first)
    # Create object #
    barrnap = RemoveITS(test_fasta)
    # Run it #
    barrnap()
    # Parse result #
    result_len = len(barrnap.filtered.first)
    # The non-gene regions to be removed #
    its_len = 222
    # Assert #
    assert orig_len - result_len == its_len
