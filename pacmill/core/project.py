#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Third party modules #
import pandas

# First party modules #
from autopaths import Path
from autopaths.dir_path import DirectoryPath
from plumbing.cache import property_cached

# Internal modules #

###############################################################################
class Project:
    """A Project object regroups several Sample objects together."""

    def __repr__(self):
        return '%s object code "%s"' % (self.__class__, self.short_name)

    def __iter__(self):
        return iter(self.samples)

    def __len__(self):
        return len(self.samples)

    def __init__(self, short_name, *all_xlsx):
        """
        A Project object takes:

         * A project `short_name` attribute that will be used to pick only the
           samples that are associated to that particular project.
         * One or more <xlsx> files that describe samples locations and their
           metadata. See example excel file for more information.

        Once the project is created you can access its attributes such as:

            >>> proj = Project('demo_project', '~/test/metadata.xlsx')
            >>> for sample in proj: print(sample.read_count)

        Attributes are the following:

            * self.short_name: the name of this project.
            * self.long_name:  the more lengthy description of this project.
            * self.all_xlsx:   a list of file paths (that are excel files).
            * self.metadata:   a pandas dataframe with all metadata
                               for this project combined.
            * self.samples:    a list of Sample objects.

        Other properties are described in their respective docstrings.
        """
        # The name of the project #
        self.short_name = short_name.lower()
        # Check it contains only alphanumerics and underscore #
        if not self.short_name.isidentifier():
            msg = "The short name of a project can only contain " \
                  "alphanumerical characters and underscores."    \
                  "Currently it is: \n\n   `%s`\n"
            raise ValueError(msg % short_name)
        # You need at least one excel file #
        assert len(all_xlsx) > 0
        # Save all excel file paths as Paths objects #
        self.all_xlsx = list(map(Path, all_xlsx))
        # Check all the excel paths actually exist #
        assert all(xlsx.exists for xlsx in self.all_xlsx)

    #----------------------------- Properties --------------------------------#
    @property_cached
    def metadata(self):
        """
        Return a pandas.DataFrame object describing the metadata of all
        samples contained in this project.
        """
        # Function to read one excel file #
        read_excel = lambda path: pandas.read_excel(str(path), header=1)
        # Read all excel files as data frames #
        all_dfs = [read_excel(path) for path in self.all_xlsx]
        # If there are several excel files, merge them together #
        metadata = pandas.concat(all_dfs, sort=False)
        # Filter and take only samples that match this project's short_name #
        query = f'project_short_name == "{self.short_name}"'
        metadata = metadata.query(query).copy()
        # Return #
        return metadata

    @property_cached
    def samples(self):
        """Create all the Sample objects."""
        # Remove samples that are not marked as "yes" for "used" #
        metadata = self.metadata.query('used == "yes"').copy()
        # Iterate over all rows of the input #
        rows = metadata.iterrows()
        # Import object #
        from pacmill.core.sample import Sample
        # Make one Sample object per row #
        samples = [Sample(self, **dict(row)) for i, row in rows]
        # Add a reference to the current project #
        for sample in samples: sample.parent = self
        # Check we have at least one sample #
        assert len(samples) > 0
        # Return #
        return samples

    @property_cached
    def long_name(self):
        """
        Get the project's long name which is a string describing the
        project title in a longer fashion than the short name.
        """
        # It's actually contained in the Sample objects #
        all_names = set(s.project_long_name for s in self)
        # Check that it doesn't diverge between samples #
        if not len(all_names) == 1:
            msg = "The project long name is not uniform across samples."
            raise ValueError(msg)
        # Return #
        return all_names.pop()

    @property_cached
    def output_dir(self):
        """Get the project's output directory."""
        # It's actually contained in the Sample objects #
        all_output_dirs = set(s.output_dir for s in self)
        # Check that it doesn't diverge between samples #
        if not len(all_output_dirs) == 1:
            msg = "The project output directory is not uniform across samples."
            raise ValueError(msg)
        # Return #
        return DirectoryPath(all_output_dirs.pop())

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /reads/all_reads.fasta
                /otus/consensus.fasta
                /otus/consensus.tsv
                /barrnap/results.gff
                /barrnap/only_16s.fasta
                /taxonomy/
                /taxa_tables/
                /graphs/
                /report/cache/
                /report/project.pdf
                /bundle/
                """

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.output_dir, self.all_paths)

    #------------------------------- Methods ---------------------------------#
    def combine_reads(self, verbose=True, check=False):
        """
        This method will concatenate all the cleaned reads from every sample
        within this project into a singe FASTA file for further processing.
        """
        # Print a message #
        if verbose:
            msg = "Combining all reads for project '%s' (%i samples)"
            msg = msg % (self.short_name, len(self.samples))
            print(msg)
        # Get all input paths #
        inputs = [sample.chimeras.results for sample in self]
        # Make the command to be run #
        cmd = 'cat %s > %s' % (' '.join(inputs), self.fasta)
        # We don't want python to be buffering the text for speed #
        from shell_command import shell_output
        shell_output(cmd)
        # Sanity check the total #
        if check:
            before = sum(len(s.chimeras.results) for s in self)
            after  = self.fasta.count
            assert before == after
        # Return #
        return self.fasta

    #---------------------------- Compositions -------------------------------#
    @property_cached
    def fasta(self):
        """
        The reads from all samples combined into the same file accessed through
        a FASTA object with convenience methods.
        See https://github.com/xapple/fasta#usage
        """
        from fasta import FASTA
        return FASTA(self.autopaths.all_reads)

    @property_cached
    def otus(self):
        """
        Takes care of running a single-pass, greedy centroid-based clustering
        algorithm on all sequences to determine consensus sequences of OTUs.
        """
        from pacmill.centering.vsearch import ClusterVsearch
        return ClusterVsearch(self.fasta,
                              self.autopaths.otus_fasta,
                              self.autopaths.otus_tsv)

    @property_cached
    def otu_table(self):
        """
        Takes care of storing a table with OTU as rows and samples as columns
        and tracks how many sequences where found from each sample in each OTU.
        Also responsible for drawing graphs of OTU distribution.
        """
        from pacmill.centering.otu_table import OtuTable
        return OtuTable(self.otus.table)

    @property_cached
    def barrnap(self):
        """
        Takes care of extracting the 16S rRNA portion from each OTU sequence
        in preparation for a similarity search against a taxonomic database.
        """
        from pacmill.filtering.barrnap import Barrnap
        return Barrnap(self.otus.results,
                       self.autopaths.barrnap_gff,
                       self.autopaths.barrnap_16s)

    @property_cached
    def taxonomy(self):
        """
        Will compare all OTU sequences against a database of curated 16S genes
        to associate a taxonomic assignment where possible.
        """
        from pacmill.taxonomy.mothur_classify import MothurClassify
        return MothurClassify(self.barrnap.results,
                              self.autopaths.taxonomy_dir)

    @property_cached
    def taxa_tables(self):
        """
        By using the OTU table along with the taxonomic assignment results,
        we can generate taxa tables at different ranks.
        """
        from pacmill.taxonomy.taxa_tables import TaxaTable
        return TaxaTable(self.otu_table,
                         self.taxonomy,
                         self.autopaths.taxa_tables_dir)

    @property_cached
    def nmds_graph(self):
        """
        Stands for Non-metric multidimensional scaling.
        Using the information in the OTU table along with a distance metric
        such as the one developed by Horn 1966 (adapted from Morisita 1959),
        we can place every sample on a two-dimensional ordination plot.
        """
        from pacmill.statistics.nmds import GraphNMDS
        return GraphNMDS(self.otu_table, self.autopaths.graphs_dir)

    @property_cached
    def report(self):
        """
        The ProjectReport object is used to produce a PDF document containing
        all the information and graphs that concern this project.
        """
        from pacmill.reports.project import ProjectReport
        return ProjectReport(self, self.autopaths.report_pdf)

    @property_cached
    def bundle(self):
        """
        The Bundle object is used to regroup all PDF reports together in
        the same directory/zip-file for easy distribution to stakeholders.
        """
        from pacmill.distribute.bundle import Bundle
        return Bundle(self, self.autopaths.bundle_dir)
