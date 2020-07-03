#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
Contact at www.sinclair.bio
MIT Licensed.

Development script to test some of the methods in `pacmill`

Typically you would run this file from a command line like this:

     ipython3 -i -- ~/deploy/pacmill/scripts/test_dev.py
"""

# Built-in modules #

# Internal modules #
from pacmill.core.project import Project

# Constants #
proj_xls = "/home/sinclair/deploy/collab_sinclair/pacmill_projects/aj_skin/metadata_aj_skin.xlsx"

###############################################################################
proj = Project('aj_skin', proj_xls)
for sample in proj: print(sample.path.size)
