#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from pacmill.reports.taxonomy import TaxonomyReport, TaxonomyTemplate

# First party modules #
from plumbing.cache    import property_cached
from autopaths         import Path

# Third party modules #

###############################################################################
class NcbiBlastReport(TaxonomyReport):
    """
    A full report generated in PDF for the NCBI blast taxonomic classification.
    """

    def __init__(self, tax, output_path):
        # References to parent objects #
        self.tax    = tax
        self.proj   = tax.proj
        # For string representation we have to select a parent #
        self.parent = tax
        # The output location #
        self.output_path = Path(output_path)

    @property_cached
    def template(self):
        return NcbiBlastTemplate(self)

###############################################################################
class NcbiBlastTemplate(TaxonomyTemplate):
    """All the parameters to be rendered in the markdown template."""

    def __init__(self, parent):
        # References to parent object #
        self.report = parent
        # Reference to the taxonomy object #
        self.tax = self.report.tax
        # Reference to the project object #
        self.project = self.report.proj
        # Reference to the samples #
        self.samples = self.project.samples

    def min_e_value(self):
        return self.tax.min_e_value

    def min_perc_identity(self):
        return self.tax.min_perc_identity

    def max_target_seqs(self):
        return self.tax.max_target_seqs

    def blast_result_table(self):
        # Call the dataframe function #
        df = self.tax.all_otus_df(1).head(10)
        # Filter columns we want #
        df = df[['abundance', 'description', 'e_value', 'length', 'identity']]
        # Rename some columns #
        df = df.rename(columns={'abundance': 'size'})
        # Convert to text #
        table = df.to_markdown(index    = False,
                               numalign = "right",
                               tablefmt = "pipe")
        # Add the caption #
        return table + "\n\n   : Classification summary for OTUs."
