#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from pacmill.taxonomy.mothur_classify import MothurClassify

# First party modules #
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class MultiTaxDatabases:
    """
    This object should be attached to a Project object and will enable OTU
    sequences to be matched for taxonomic assignments on several databases.

    Currently we support:

        * Silva.
        * Greengenes.
        * RDP.
        * CREST.
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.proj)

    def __init__(self, proj, base_dir):
        self.proj     = proj
        self.base_dir = base_dir

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /silva/
                /greengenes/
                /rdp/
                /crest/
                /reports/
                """

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.base_dir, self.all_paths)

    #------------------------------ Running ----------------------------------#
    @property_cached
    def taxonomies(self):
        return [self.silva, self.greengenes, self.rdp]

    def __call__(self, verbose=True):
        # Loop over databases #
        for tax in self.taxonomies:
            tax(verbose=verbose)
        # Make all tables #
        for table in self.tables.all:
            table(verbose=verbose)

    #---------------------------- Compositions -------------------------------#
    @property_cached
    def silva(self):
        from seqsearch.databases.mothur.silva import silva_mothur
        return MothurClassify(self.proj.otus.results,
                              silva_mothur,
                              self.autopaths.silva_dir)

    @property_cached
    def greengenes(self):
        from seqsearch.databases.mothur.greengenes import gg_mothur
        return MothurClassify(self.proj.otus.results,
                              gg_mothur,
                              self.autopaths.greengenes_dir)

    @property_cached
    def rdp(self):
        from seqsearch.databases.mothur.rdp import rdp_mothur
        return MothurClassify(self.proj.otus.results,
                              rdp_mothur,
                              self.autopaths.rdp_dir)

    @property_cached
    def crest(self):
        from pacmill.taxonomy.crest import CrestClassify
        return CrestClassify(self.proj.otus.results,
                            self.autopaths.crest_dir)

    #--------------------------------- Tables --------------------------------#
    @property_cached
    def tables(self):
        """
        By using the OTU table along with the taxonomic assignment results,
        we can generate taxa tables at different ranks.

        Typically you can access these tables like this:

            >>> print(project.taxonomy.tables.silva)
            >>> print(project.taxonomy.tables.rdp)
        """
        # Import #
        from pacmill.taxonomy.taxa_tables import TaxaTable
        # Make a dummy object #
        result = type('TaxaTablesCollection', (), {})
        # Prepare an empty list #
        result.all = []
        # Loop over taxonomies #
        for tax in self.taxonomies:
            # Instantiate #
            table = TaxaTable(self.proj.otu_table, tax,
                              tax.autopaths.taxa_tables_dir)
            # Adjust the number of taxa displayed #
            max_taxa = getattr(self.proj.samples[0], 'max_taxa', None)
            if max_taxa: result.max_taxa_displayed = max_taxa
            # Set attribute #
            setattr(result, tax.database.tag, table)
            # Add it to list #
            result.all.append(table)
        # Return #
        return result

    #-------------------------------- Reports --------------------------------#
    @property_cached
    def reports(self):
        """
        Each different taxonomic classification method has its own PDF
        report automatically generated.
        """
        # Import #
        from pacmill.reports.taxonomy import TaxonomyReport
        # Make a dummy object #
        result = type('TaxaReportCollection', (), {})
        # Prepare an empty list #
        result.all = []
        # Loop over taxonomies #
        for tax in self.taxonomies:
            # Path to report #
            path =  self.autopaths.reports_dir + tax.database.tag + '.pdf'
            # Matching taxonomy tables #
            tables = getattr(self.tables, tax.database.tag)
            # Instantiate #
            report = TaxonomyReport(tax, tables, self.proj, path)
            # Set attribute #
            setattr(result, tax.short_name, report)
            # Add it to list #
            result.all.append(report)
        # Return #
        return result
