#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Script to run the `pacmill` pipeline.

Typically you would run this file from a command line like this:

    ipython3 -i -- ~/deploy/pacmill/pacmill/run_pacmill.py \
                   small_test small_test/metadata_small_test.xlsx
"""

# Matplotlib #
import matplotlib
matplotlib.use('Agg')

# Built-in modules #
import os, sys, argparse

# Internal modules #
from pacmill.core.project import Project

# First party modules #
from plumbing.timer import Timer

# Third party modules #
from p_tqdm import p_umap

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

    # Duplicate all output to a text log created in the project directory #
    class UnbufferedLoggerDuplicator:
        def __init__(self, path):
            self.orig = sys.stdout
            self.path = path
            self.log  = open(self.path, "w")

        def write(self, message):
            self.orig.write(message)
            self.log.write(message)
            self.orig.flush()
            self.log.flush()

        def flush(self):
            self.orig.flush()
            self.log.flush()

    sys.stdout = UnbufferedLoggerDuplicator(proj.autopaths.log)

    # Message #
    print("------------------------------------------")
    print("Running pacmill on project '%s'." % proj.short_name)
    print("------------------------------------------")

    # Start the timer #
    timer = Timer()
    timer.print_start()

    # Start the pipeline #
    print("# Validate the format of the FASTQs #")
    p_umap(lambda s: s.fastq.validator(), proj)
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
    p_umap(lambda s: s.filter(), proj)
    timer.print_elapsed()

    print("# Remove chimeric reads #")
    p_umap(lambda s: s.chimeras(cpus=1), proj)
    timer.print_elapsed()

    print("# Detect presence of rRNA genes (optional) #")
    for sample in proj:
        if sample.barrnap_mode != 'off':
            print(sample.barrnap())
    timer.print_elapsed()

    print("# Concatenate reads from all samples into one file #")
    print(proj.combine_reads())
    timer.print_elapsed()

    print("# Pick OTUS #")
    print(proj.otus())
    timer.print_elapsed()

    if proj.check_homogeneous("run_ncbi_blast"):
        print("# BLAST sequences against the NCBI 16S database #")
        proj.ncbi_blast()
        timer.print_elapsed()

    print("# Assign taxonomy and make all taxa tables #")
    proj.taxonomy()
    timer.print_elapsed()

    print("# Regenerate the graphs for the project #")
    print(proj.otu_table.graphs.otu_sums_graph(rerun=True))
    print(proj.otu_table.graphs.sample_sums_graph(rerun=True))
    print(proj.otu_table.graphs.cumulative_presence(rerun=True))
    print(proj.nmds_graph(rerun=True))
    timer.print_elapsed()

    print("# Regenerate the graphs and legends for taxa bar-stacks #")
    for tables in proj.taxonomy.tables.all:
        if not tables.taxonomy.should_run: continue
        for graph in tables.results.graphs.by_rank:
            print(graph(rerun=True))
        for legend in tables.results.graphs.legends:
            print(legend(rerun=True))
    timer.print_elapsed()

    # Clear the cache #
    for sample in proj:
        sample.report.template.cache_dir.remove()

    print("# Create the PDF reports for each sample #")
    for sample in proj:
        print(sample.report())
    timer.print_elapsed()

    print("# Create the PDF reports for each taxonomic classification #")
    for report in proj.taxonomy.reports.all:
        if report.tax.should_run:
            print(report())
    timer.print_elapsed()

    if proj.check_homogeneous("run_ncbi_blast"):
        print("# Create the PDF report for the NCBI BLAST #")
        proj.ncbi_blast.results.report()
        timer.print_elapsed()

    print("# Create the PDF report for the project #")
    print(proj.report())
    timer.print_elapsed()

    print("# Create a bundle for distribution #")
    print(proj.bundle())
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
