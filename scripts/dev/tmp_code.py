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
import os

# Internal modules #
from pacmill.core.project import Project

# Third party modules #
from tqdm import tqdm

# Constants #
proj_name = os.environ.get("PACMILL_PROJ_NAME", "No project has been set.")
proj_xls = os.environ.get("PACMILL_PROJ_XLS", "No excel path has been set.")

###############################################################################
# Create project #
proj = Project(proj_name, proj_xls)

# Show samples #
for sample in proj: print(sample.short_name)

# Test PHRED window size #
for sample in proj:
    import pandas
    quality = pandas.Series(range(100))
    rolling = quality.rolling(sample.phred_window_size).mean().dropna().tolist()
    print(rolling)