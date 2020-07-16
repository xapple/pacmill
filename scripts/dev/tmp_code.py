#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Development script to test some of the methods in `pacmill`
and try out different things. This script can safely be ignored.

Typically you would run this file from a command line like this:

     ipython3 -i -- ~/deploy/pacmill/scripts/dev/tmp_code.py
"""

# Built-in modules #

# Internal modules #
from pacmill.core.project import Project

# Constants #
proj_xls = "/home/sinclair/deploy/collab_sinclair/pacmill_projects/aj_skin/metadata_aj_skin.xlsx"

###############################################################################
# Create project #
proj = Project('aj_skin', proj_xls)

for sample in proj:
   print(sample.filter(True))

for sample in proj:
   print(sample.report())
