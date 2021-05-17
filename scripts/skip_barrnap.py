#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Development script to skip some of the steps in the `pacmill` pipeline
and try out different things. This script can safely be ignored.

Typically you would run this file from a command line like this:

     ipython3 -i -- ~/deploy/pacmill/scripts/dev/skip_barrnap.py
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

# Copy file to new location directly #
for sample in proj:
    sample.chimeras.results.copy(sample.barrnap.results)
