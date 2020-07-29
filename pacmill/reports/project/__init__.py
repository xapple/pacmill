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
from pymarktex.figures import ScaledFigure
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
        label   = "input_length_dist"
        return str(ScaledFigure(path, caption, label))

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
        return str(ScaledFigure(path, caption))

    def sample_sums_graph(self):
        caption = "Distribution of OTU presence per sample"
        path    = self.project.otu_table.graphs.sample_sums_graph()
        return str(ScaledFigure(path, caption))

    def cumulative_presence(self):
        caption = "Cumulative number of reads by OTU presence"
        path    = self.project.otu_table.graphs.cumulative_presence()
        return str(ScaledFigure(path, caption))

    #------------------------------ Barrnap ----------------------------------#
    def barrnap(self):
        return bool(self.project.barrnap)

    def barrnap_discard(self):
        before = self.project.otus.results.count
        after  = self.project.barrnap.results.count
        return thousands(before - after)

    def barrnap_left(self):
        return thousands(self.project.barrnap.results.count)

    #------------------------------- Taxonomy --------------------------------#
    def taxonomy(self):
        return bool(self.project.taxonomy)

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
    def taxa_barstack_at_rank(self, rank, label):
        graphs  = self.project.taxa_tables.results.graphs.by_rank
        graph   = [g for g in graphs if g.base_rank == rank][0]
        caption = "Relative abundances per sample on the '%s' level"
        return str(ScaledFigure(graph(), caption % graph.label, label))

    def level_one_barstack(self):
        return self.taxa_barstack_at_rank(2, "level_one_barstack")

    def level_two_barstack(self):
        return self.taxa_barstack_at_rank(3, "level_two_barstack")

    def level_three_barstack(self):
        return self.taxa_barstack_at_rank(4, "level_three_barstack")

    #------------------------------ Comparison -------------------------------#
    def comparison(self):
        return False
