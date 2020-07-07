#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Script to run the `pacmill` pipeline on the first test project
called `aj_skin`.

Typically you would run this file from a command line like this:

     ipython3 -i -- ~/deploy/pacmill/scripts/running/run_aj_skin.py
"""

# Built-in modules #

# Internal modules #
from pacmill.core.project import Project

# Constants #
proj_xls = "/home/sinclair/deploy/collab_sinclair/pacmill_projects/aj_skin/metadata_aj_skin.xlsx"

###############################################################################
# Create project #
proj = Project('aj_skin', proj_xls)

# Validate the format of the FASTQs #
for sample in proj: print(sample.fastq.validator())

# Run FastQC on the samples individually #
for sample in proj:
    print(sample.fastq.fastqc())
