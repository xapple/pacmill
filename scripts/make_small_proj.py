#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

This script takes a pacmill project along with all the information describing
the samples, and creates a new 'small_test' project that is identical to the
original project except for the number of sequences of the raw sample FASTQ
which is randomly sub-sampled to 5'000.

Such a test project is useful for development, as all the steps of
the pipeline can be run in a much shorter time frame.

Typically you would run this file from a command line like this:

     ipython3 -i -- ~/deploy/pacmill/scripts/subsampling/make_small_proj.py
"""

# Built-in modules #
import os

# Internal modules #
from pacmill.core.project import Project

# Constants #
proj_name = os.environ.get("PACMILL_PROJ_NAME", "No base project has been set")
proj_xls = os.environ.get("PACMILL_PROJ_XLS", "No base project has been set")

###############################################################################
# Create project #
proj = Project(proj_name, proj_xls)

# Loop #
for i, sample in enumerate(proj):
    # Get location for the new input files #
    new_path = '/data/small_test/raw/small_sample_%i.fastq' % (i+1)
    print(new_path)
    # Do it #
    sample.fastq.subsample(down_to=5000, new_path=new_path)
