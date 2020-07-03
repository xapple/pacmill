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
proj_xls = "~/repos/collab_sinclair/pacmill_projects/aj_skin/200702_metadata_sinclair_lm.xlsx"

###############################################################################
proj = Project('aj_skin', '~/test/metadata.xlsx')
for sample in proj: print(sample.read_count)
