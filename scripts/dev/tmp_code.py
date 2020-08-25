#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Development script to test some of the methods in `pacmill`
and try out different things. This script can safely be ignored
and is meant simply as a sandbox.

Typically you would run this file from a command line like this:

     ipython3 -i -- ~/deploy/pacmill/scripts/dev/tmp_code.py
"""

# Built-in modules #
import os, inspect

# Internal modules #
from pacmill.core.project import Project

# Third party modules #
from tqdm import tqdm

# Constants #
proj_name = os.environ.get("PACMILL_PROJ_NAME", "No project has been set.")
proj_xls = os.environ.get("PACMILL_PROJ_XLS", "No excel path has been set.")

###############################################################################
# Import #
from pacmill.filtering.barrnap import RemoveITS
# This directory #
from autopaths import Path
this_file = Path((inspect.stack()[0])[1])
this_dir  = this_file.directory
# Get file paths #
test_fastq = this_dir + '../../test/barrnap/test_sequence.fasta'
# Create object #
barrnap = RemoveITS(test_fastq)
# Run it #
barrnap()
# Parse result #
seq_length = len(barrnap.filtered.first)
# Assert #
assert seq_length == 1556
