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
class ProjectReport(Document):
    """A full report generated in PDF for every Project object."""

    # Custom LaTeX headers and footers for the report #
    header_template = Header
    footer_template = Footer

    # Specific title for the project report #
    params = {'title': 'Auto-generated project report'}

    def __init__(self, project, output_path):
        # Reference to parent objects #
        self.parent  = project
        self.project = project
        # The output location #
        self.output_path = Path(output_path)

    @property_cached
    def template(self):
        return ProjectTemplate(self)

    def load_markdown(self):
        self.markdown = str(self.template)

###############################################################################
class ProjectTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template."""

    delimiters = (u'{{', u'}}')

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, parent):
        # Reference to parent objects #
        self.parent = parent
        self.report = parent
        # Reference to the project object #
        self.project = parent.project
        # Reference to the samples #
        self.samples = self.project.samples
        # Where to pickle properties that are costly to compute #
        self.cache_dir = self.project.autopaths.report_cache_dir

    #------------------------- General information ---------------------------#
    def proj_short_name(self):
        return self.project.short_name

    def proj_long_name(self):
        return self.project.long_name

    def count_samples(self):
        return len(self.project.samples)

    def count_sequences(self):
        return thousands(self.project.fasta.count)

    def input_length_dist(self):
        caption = "Distribution of sequence lengths at input"
        path    = self.project.fasta.graphs.length_hist()
        return str(ScaledFigure(path, caption, "input_length_dist"))

    #----------------------------- Summary table -----------------------------#
    def sample_table(self):
        # The functions #
        name = lambda s: "**" + s.short_name + "**"
        desc = lambda s: s.long_name
        lost = lambda s: "%.1f%%" % s.percent_lost
        left = lambda s: thousands(len(s.filter.results.clean))
        # The columns #
        info = {
          'Name':        name,
          'Description': desc,
          'Reads lost':  lost,
          'Reads left':  left,
        }
        # The table contents, row after row #
        table = [[i+1] + [f(s) for f in info.values()]
                 for i, s in enumerate(self.samples)]
        # The title row #
        headers = ['#'] + list(info.keys())
        # Make it as text #
        from tabulate import tabulate
        table = tabulate(table, headers, numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Summary information for all samples."

    #------------------------------ OTU making -------------------------------#
    def otus(self):
        return bool(self.project.otus)

    def otus_threshold(self):
        return "%.1f%%" % (self.project.otus.threshold * 100)

    def otus_total(self):
        return thousands(self.project.otus.results.count)

    def otu_sums_graph(self):
        caption = "Distribution of OTU presence per OTU"
        path    = self.project.otu_table.graphs.otu_sums_graph()
        return str(ScaledFigure(path, caption, 'otu_sums_graph'))

    def sample_sums_graph(self):
        caption = "Distribution of OTU presence per sample"
        path    = self.project.otu_table.graphs.sample_sums_graph()
        return str(ScaledFigure(path, caption, 'sample_sums_graph'))

    def cumulative_presence(self):
        caption = "Cumulative number of reads by OTU presence"
        path    = self.project.otu_table.graphs.cumulative_presence()
        return str(ScaledFigure(path, caption, 'cumulative_presence'))

    #------------------------------- Taxonomy --------------------------------#
    def taxonomy(self):
        return bool(self.project.taxonomy) and bool(self.project.taxa_tables)

    def classify_citation(self):
        return "the '%s' method" % self.project.taxonomy.long_name

    def classify_database(self):
        return "'" + self.project.taxonomy.database.long_name + "'"

    def otu_classified_table(self):
        # The functions #
        def rank(i):
            return "**" + self.project.taxonomy.database.rank_names[i] + "**"
        def classified(i):
            return self.project.taxonomy.results.count_assigned[i]
        def unclassified(i):
            return self.project.taxonomy.results.count_unassigned[i]
        # The columns #
        info = {
            'Rank':         rank,
            'Classified':   classified,
            'Unclassified': unclassified,
        }
        # The table contents, row after row #
        table = [[i+1] + [f(i) for f in info.values()]
                 for i in range(len(self.project.taxonomy.database.rank_names))]
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
        graphs = self.project.taxa_tables.results.graphs.by_rank
        graph  = [g for g in graphs if g.base_rank == rank][0]
        # Get the legend #
        legends = self.project.taxa_tables.results.graphs.legends
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
        names = [n.lower() for n in self.project.taxa_tables.rank_names]
        ranks = [names.index(r) for r in ranks]
        # Produce as many graphs as needed #
        graphs = [self.taxa_barstack_at_rank(r) for r in ranks]
        # Separate graphs with newlines #
        graphs = '\n\n'.join(graphs)
        # Return
        return graphs

    #------------------------------ Comparison -------------------------------#
    def comparison(self):
        if len(self.samples) < 2: return False
        return bool(self.project.nmds_graph)

    def otu_nmds(self):
        caption = "NMDS using the OTU table for %i samples" % len(self.samples)
        graph   = self.project.nmds_graph()
        return str(ScaledFigure(graph, caption, 'otu_nmds'))
