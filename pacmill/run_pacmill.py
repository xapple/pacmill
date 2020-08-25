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
    proj_name = env_name or args.proj_name
    if proj_name is None:
        raise Exception("No short project name has been given.")

    # Get defaults for project excel file #
    proj_xls = args.proj_xls or env_xls
    if proj_xls is None:
        raise Exception("No excel file path has been given.")

    # Create project #
    proj = Project(args.proj_name, args.proj_xls)

    # Validate the format of the FASTQs #
    for sample in proj: print(sample.fastq.validator())

    # Run FastQC on the samples individually #
    for sample in proj:
        print(sample.fastq.fastqc())

    # Regenerate the graphs for samples #
    for sample in proj:
        print(sample.fastq.graphs.length_hist(rerun=True))

    # Filter reads in every sample based on several criteria #
    for sample in proj:
        print(sample.filter(verbose=True))

    # Remove chimeric reads #
    for sample in proj:
        print(sample.chimeras())

    # Detect presence of rRNA genes #
    for sample in proj:
        print(sample.barrnap())

    # Concatenate reads from all samples into one file #
    print(proj.combine_reads())

    # Pick OTUS #
    print(proj.otus())

    # Assign taxonomy and make all taxa tables #
    print(proj.taxonomy())

    # Regenerate the graphs for the project #
    print(proj.otu_table.graphs.otu_sums_graph(rerun=True))
    print(proj.otu_table.graphs.sample_sums_graph(rerun=True))
    print(proj.otu_table.graphs.cumulative_presence(rerun=True))
    print(proj.nmds_graph(rerun=True))

    # Regenerate the graphs and legends for taxa bar-stacks #
    for tables in proj.taxonomy.tables.all:
        for graph in tables.results.graphs.by_rank:
            print(graph(rerun=True))
        for legend in tables.results.graphs.legends:
            print(legend(rerun=True))

    # Clear the cache #
    for sample in proj:
        print(sample.report.template.cache_dir.remove())

    # Create the PDF reports for each sample #
    for sample in proj:
        print(sample.report())

    # Create the PDF reports for each taxonomic classification #
    for report in proj.taxonomy.reports.all:
        print(report())

    # Create the PDF report for the project #
    print(proj.report())

    # Create a bundle for distribution #
    proj.bundle()

    # Success message #
    print("The pipeline was successfully run on project '%s'." % proj_name)

    # Bundle message #
    print("To quickly download the sample and project reports to your local"
          " computer, you can use the following command:"
          " \n\n%s\n" % proj.bundle.results.rsync)
