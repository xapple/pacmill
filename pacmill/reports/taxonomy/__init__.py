#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from pacmill.reports.base_template import ReportTemplate
from pacmill.reports.template      import Header, Footer

# First party modules #
from plumbing.cache    import property_cached
from plumbing.common   import split_thousands as thousands
from pymarktex         import Document
from pymarktex.figures import ScaledFigure, BareFigure
from autopaths         import Path

# Third party modules #

###############################################################################
class TaxonomyReport(Document):
    """
    A full report generated in PDF for every taxonomic database that
    we assign against.
    """

    # Custom LaTeX headers and footers for the report #
    header_template = Header
    footer_template = Footer

    # Specific title for the project report #
    params = {'title': 'Auto-generated taxonomy report'}

    def __init__(self, tax, tables, proj, output_path):
        # References to parent objects #
        self.tax    = tax
        self.tables = tables
        self.proj   = proj
        # For string representation we have to select a parent #
        self.parent = tax
        # The output location #
        self.output_path = Path(output_path)

    @property_cached
    def template(self):
        return TaxonomyTemplate(self)

    def load_markdown(self):
        self.markdown = str(self.template)

###############################################################################
class TaxonomyTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template."""

    delimiters = (u'{{', u'}}')

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.report)

    def __init__(self, parent):
        # References to parent object #
        self.report = parent
        # Reference to the taxonomy object #
        self.tax = self.report.tax
        # Reference to the taxonomy tables object #
        self.tables = self.report.tables
        # Reference to the project object #
        self.project = self.report.proj
        # Reference to the samples #
        self.samples = self.project.samples

    #------------------------- General information ---------------------------#
    def tax_tag(self):
        return self.tax.database.tag

    def proj_short_name(self):
        return self.project.short_name

    def proj_long_name(self):
        return self.project.long_name

    #------------------------------- Taxonomy --------------------------------#
    def taxonomy(self):
        return bool(self.tax) and bool(self.tables)

    def tax_long_method(self):
        return self.tax.long_name

    def tax_long_database(self):
        return self.tax.database.long_name

    def otus_count(self):
        return thousands(self.project.otus.results.count)

    def otu_classified_table(self):
        # The functions #
        def rank(i): return "**" + self.tax.database.rank_names[i] + "**"
        def classified(i): return self.tax.results.count_assigned[i]
        def unclassified(i):  return self.tax.results.count_unassigned[i]
        # The columns #
        info = {
            'Rank':         rank,
            'Classified':   classified,
            'Unclassified': unclassified,
        }
        # The table contents, row after row #
        table = [[i+1] + [f(i) for f in info.values()]
                 for i in range(len(self.tax.database.rank_names))]
        # The title row #
        headers = ['#'] + list(info.keys())
        # Make it as text #
        from tabulate import tabulate
        table = tabulate(table, headers, numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Classification summary for OTUs."

    #--------------------------- Taxa Table Graphs ---------------------------#
    def taxa_barstack_at_rank(self, rank, label=None):
        # Default label #
        if label is None: label = "taxa_barstack_%i" % rank
        # Get the graph itself #
        graphs = self.tables.results.graphs.by_rank
        graph  = [g for g in graphs if g.base_rank == rank][0]
        # Get the legend #
        legends = self.tables.results.graphs.legends
        legend  = [leg for leg in legends if leg.base_rank == rank][0]
        # Format the graph #
        graph  = str(BareFigure(graph()))
        # Format the legend #
        caption = "Relative abundances per sample on the '%s' level"
        legend  = str(ScaledFigure(legend(), caption % legend.label, label))
        # Combine both graphs #
        return graph + '\n\n' + legend

    def taxa_barstacks(self):
        # Get the parameter in the excel file #
        ranks = getattr(self.samples[0], 'taxa_barstacks')
        # If none are included skip this step #
        if ranks is None: return "\n"
        # Split on the comma and remove spaces #
        ranks = [r.strip() for r in ranks.split(',')]
        # Convert named ranks to numbered ranks #
        names = [n.lower() for n in self.tables.rank_names]
        ranks = [names.index(r) for r in ranks]
        # Produce as many graphs as needed #
        graphs = [self.taxa_barstack_at_rank(r) for r in ranks]
        # Separate graphs with newlines #
        graphs = '\n\n'.join(graphs)
        # Return #
        return graphs