#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Script to run the `pacmill` pipeline.
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

# Validate the format of the FASTQs #
for sample in proj: print(sample.fastq.validator())

# Run FastQC on the samples individually #
for sample in proj:
    print(sample.fastq.fastqc())

# Filter reads in every sample based on several criteria #
for sample in proj:
    print(sample.filter(verbose=True))

# Detect presence of rRNA genes #
for sample in proj:
    print(sample.barrnap())
    print(sample.barrnap.filter())

# Remove chimeric reads #
for sample in proj:
    print(sample.chimeras())

# Concatenate reads from all samples to one file #
print(proj.combine_reads())

# Pick OTUS #
print(proj.otus())

# Extract 16S location within OTUs #
print(proj.barrnap())
print(proj.barrnap.extract())

# Assign taxonomy #
print(proj.taxonomy())

# Make all taxa tables #
print(proj.taxa_tables())

# Regenerate the graphs for samples #
for sample in proj:
    print(sample.fastq.graphs.length_hist(rerun=True))

# Regenerate the graphs for the project #
print(proj.otu_table.graphs.otu_sums_graph(rerun=True))
print(proj.otu_table.graphs.sample_sums_graph(rerun=True))
print(proj.otu_table.graphs.cumulative_presence(rerun=True))

# Regenerate the graphs for taxa bar-stacks #
for g in proj.taxa_tables.results.graphs.by_rank: print(g(rerun=True))

# Clear the cache #
for sample in proj:
    print(sample.report.cache_dir.remove())

# Create the PDF reports for each sample #
for sample in proj:
    print(sample.report())

# Create the PDF report for the project #
print(proj.report())