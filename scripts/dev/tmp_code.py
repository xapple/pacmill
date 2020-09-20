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

proj.taxonomy.tables.silvamod128()

#for sample in proj:
#    print(sample.fastq.fastqc())
#    break

#for report in proj.taxonomy.reports.all:
#    print(report())

#proj.otus()
#proj.taxonomy()

#for tables in proj.taxonomy.tables.all:
#    for graph in tables.results.graphs.by_rank:
#        print(graph(rerun=True))
#    for legend in tables.results.graphs.legends:
#        print(legend(rerun=True))

#for i, rank in tqdm(enumerate(proj.taxa_tables.rank_names)):
#    proj.taxa_tables.results.graphs.by_rank[i](rerun=True)
#    proj.taxa_tables.results.graphs.legends[i](rerun=True)

#print(proj.otu_table.graphs.otu_sums_graph(rerun=True))
#print(proj.otu_table.graphs.sample_sums_graph(rerun=True))
#print(proj.otu_table.graphs.cumulative_presence(rerun=True))

#for i, rank in tqdm(enumerate(proj.taxa_tables.rank_names)):
#    proj.taxa_tables.results.graphs.by_rank[i](rerun=True)
#    proj.taxa_tables.results.graphs.legends[i](rerun=True)

#for sample in proj:
#    print(sample.report())
#    break

#proj.nmds_graph()

#proj.report()

#proj.bundle()
#print(proj.bundle.results.rsync)
