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

# Constants #
proj_name = os.environ.get("PACMILL_PROJ_NAME", "No base project has been set")
proj_xls = os.environ.get("PACMILL_PROJ_XLS", "No base project has been set")

###############################################################################
# Create project #
proj = Project(proj_name, proj_xls)

#for sample in proj:
#    print(sample.filter(True))
#    break

#for sample in proj:
#    result = sample.chimeras()
#    print('\n'.join(result.directory.flat_contents))
#    break

#for sample in proj:
#   print(sample.report())

#proj.otus()

#proj.barrnap.extract()

#proj.taxonomy()

#print(proj.otu_table.graphs.otu_sums_graph(rerun=True))
#print(proj.otu_table.graphs.sample_sums_graph(rerun=True))
#print(proj.otu_table.graphs.cumulative_presence(rerun=True))

#print(proj.taxa_tables())

#for g in proj.taxa_tables.results.graphs.by_rank: print(g(rerun=True))

proj.nmds_graph()

#proj.report()

#proj.bundle()
#print(proj.bundle.results.rsync)