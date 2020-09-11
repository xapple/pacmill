#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Script to run the `pacmill` pipeline.

Typically you would run this file from a command line like this:

    ipython3 -i -- ~/deploy/pacmill/scripts/running/run_pacmill.py \
                   small_test small_test/metadata_small_test.xlsx
"""

# Built-in modules #
import os, argparse

# Internal modules #
from pacmill.core.project import Project

# First party modules #
from plumbing.timer import Timer
from plumbing.processes import prll_map

# Constants #
env_name = os.environ.get("PACMILL_PROJ_NAME")
env_xls  = os.environ.get("PACMILL_PROJ_XLS")

###############################################################################
if __name__ == "__main__":
    # Some strings #
    proj_name_help = "Short name of the project to run."
    proj_xls_help  = "The path to the excel metadata file."
    description    = "Script to run a given project through the whole" \
                     " pacmill pipeline."

    # Parse the shell arguments #
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("proj_name", help=proj_name_help,
                        type=str, nargs='?', default=None)
    parser.add_argument("proj_xls",  help=proj_xls_help,
                        type=str, nargs='?', default=None)
    args = parser.parse_args()

    # Get defaults for project name #
    proj_name = args.proj_name or env_name
    if proj_name is None:
        raise Exception("No short project name has been given.")

    # Get defaults for project excel file #
    proj_xls = args.proj_xls or env_xls
    if proj_xls is None:
        raise Exception("No excel file path has been given.")

    # Create project #
    proj = Project(proj_name, proj_xls)

    # Message #
    print("------------------------------------------")
    print("Running pacmill on project '%s'." % proj.short_name)
    print("------------------------------------------")

    # Start the timer #
    timer = Timer()
    timer.print_start()

    print("# Validate the format of the FASTQs #")
    prll_map(lambda s: s.fastq.validator(), proj)
    timer.print_elapsed()

    print("# Run FastQC on the samples individually #")
    for sample in proj:
        print(sample.fastq.fastqc())
    timer.print_elapsed()

    print("# Regenerate the graphs for samples #")
    for sample in proj:
        print(sample.fastq.graphs.length_hist(rerun=True))
    timer.print_elapsed()

    print("# Filter reads in every sample based on several criteria #")
    prll_map(lambda s: s.filter(), proj)
    timer.print_elapsed()

    print("# Remove chimeric reads #")
    prll_map(lambda s: s.chimeras(cpus=1), proj)
    timer.print_elapsed()

    print("# Detect presence of rRNA genes (optional) #")
    for sample in proj:
        if sample.barrnap_mode == 'on':
            print(sample.barrnap())
    timer.print_elapsed()

    print("# Concatenate reads from all samples into one file #")
    print(proj.combine_reads())
    timer.print_elapsed()

    print("# Pick OTUS #")
    print(proj.otus())
    timer.print_elapsed()

    print("# Assign taxonomy and make all taxa tables #")
    print(proj.taxonomy())
    timer.print_elapsed()

    print("# Regenerate the graphs for the project #")
    print(proj.otu_table.graphs.otu_sums_graph(rerun=True))
    print(proj.otu_table.graphs.sample_sums_graph(rerun=True))
    print(proj.otu_table.graphs.cumulative_presence(rerun=True))
    print(proj.nmds_graph(rerun=True))
    timer.print_elapsed()

    print("# Regenerate the graphs and legends for taxa bar-stacks #")
    for tables in proj.taxonomy.tables.all:
        for graph in tables.results.graphs.by_rank:
            print(graph(rerun=True))
        for legend in tables.results.graphs.legends:
            print(legend(rerun=True))
    timer.print_elapsed()

    # Clear the cache #
    for sample in proj:
        print(sample.report.template.cache_dir.remove())

    print("# Create the PDF reports for each sample #")
    for sample in proj:
        print(sample.report())

    print("# Create the PDF reports for each taxonomic classification #")
    for report in proj.taxonomy.reports.all:
        print(report())

    print("# Create the PDF report for the project #")
    print(proj.report())
    timer.print_elapsed()

    print("# Create a bundle for distribution #")
    proj.bundle()
    timer.print_elapsed()

    # Success message #
    print("The pipeline was successfully run on project '%s'." % proj_name)
    print("------------------------------------------")

    # Timer message #
    timer.print_end()
    timer.print_total_elapsed()

    # Bundle message #
    print("\n To quickly download the sample and project reports to your local"
          " computer, you can use the following command:"
          " \n\n%s\n" % proj.bundle.results.rsync)
