#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import re

# First party modules #
from fasta import FASTA
from autopaths.dir_path       import DirectoryPath
from plumbing.cache           import property_cached
from plumbing.check_cmd_found import check_cmd

# Third party modules #
import pandas

###############################################################################
class BlastClassify:
    """
    Use BLAST and the NCBI 16S database to identify close matches between
    OTU centroid sequences and cultured bacterial/archaeal 16 sequences.
    """

    # Constants #
    short_name = 'ncbi_blast'
    long_name  = 'BLAST against the NCBI 16S database'

    # The database #
    db_short_name   = "ncbi_16s_rna"
    db_long_name    = 'The NCBI 16S RNA database'

    # Parameters #
    min_e_value        = 1e-5
    min_perc_identity  = 97
    max_target_seqs    = 5

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, database, dest_dir, proj):
        # Source is the FASTA file containing OTU consensus sequences #
        self.source = FASTA(source)
        # Database is an object from the `seqsearch.databases` subpackage #
        self.database = database
        # Destination is a directory that contains all the results #
        self.dest_dir = DirectoryPath(dest_dir)
        # Keep a reference to the project #
        self.proj = proj

    #----------------------------- Installing --------------------------------#
    apt_packages = ['ncbi-blast+']

    def check_installed(self, exception=True):
        """
        Try to determine if the BLAST software is installed and accessible.
        Also try to determine if the database chosen is downloaded and
        accessible.
        """
        # The database #
        if not self.database:
            msg  = "The taxonomic database does not seem to be accessible."
            msg += " More information follows.\n\n"
            msg += self.database.__doc__
            raise Exception(msg)
        # The command #
        return check_cmd('blastn', exception)

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /blast/hits.xml
                /blast/stdout.txt
                /blast/stderr.txt
                /summary/single_hit.xlsx
                /summary/five_hits.xlsx
                /report.pdf
                """

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.dest_dir, self.all_paths)

    #------------------------------ Running ----------------------------------#
    @property_cached
    def seq_search(self):
        """
        The SeqSearch object is used for quickly launching BLAST search queries
        against a database of choice starting from a FASTA object.
        See https://github.com/xapple/seqsearch
        """
        # Import #
        from seqsearch.search import SeqSearch
        # Extra options for the BLAST search #
        params = {'-evalue':          self.min_e_value,
                  '-perc_identity':   self.min_perc_identity,
                  '-max_target_seqs': self.max_target_seqs}
        # Create the object #
        return SeqSearch(input_fasta = self.source,
                         database    = self.database,
                         seq_type    = 'nucl',
                         out_path    = self.autopaths.hits,
                         _out        = self.autopaths.stdout,
                         _err        = self.autopaths.stderr,
                         params      = params)

    def __call__(self, cpus=None, verbose=True):
        # Message #
        if verbose:
            message = "Running `%s` taxonomy classification on '%s'."
            print(message % (self.short_name, self.source))
        # Check blast is installed #
        self.check_installed()
        # Run #
        self.seq_search.run()
        # Message #
        if verbose:
            message = "Generating BLAST summary Excel files at '%s'."
            print(message % self.autopaths.summary_dir)
        # Generate summary lists #
        one  = self.all_otus_df(1)
        five = self.all_otus_df(5)
        # Create excel reports #
        one.to_excel(self.autopaths.xlsx_single_hit.path, index=False)
        five.to_excel(self.autopaths.xlsx_five_hits.path, index=False)

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the taxonomic software was run already and the results
        are stored on the filesystem. Return False if it was not yet run.
        """
        return self.autopaths.hits.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from the BLAST taxonomic" \
                  " classification before running the tool."
            raise Exception(msg)
        # Return the results #
        return BlastResults(self)

    #------------------------------ Dataframes -------------------------------#
    # Define columns needed #
    columns = ['abundance', 'otu_id', 'description', 'score', 'e_value',
               'accession', 'length', 'cover', 'identity']

    def all_otus_df(self, num_hits_per_otu=5):
        """
        Returns a dataframe object summarizing all the top hits for all OTUs.
        The results looks like:

        abundance       otu_id  hit_number  ...  accession  length  identity
             4233    small_1:0           1  ...  NR_113957    1475     99.9%
             4233    small_1:0           2  ...  NR_036904    1472    100.0%
             4233    small_1:0           3  ...  NR_113405    1488     99.5%
             4233    small_1:0           4  ...  NR_118997    1524     98.6%
        """
        # Parse the BLAST XML output file with BioPython #
        all_otus_hits = self.seq_search.results
        # Produce one mini-dataframe for every OTU #
        all_otus_dfs = [self.one_otu_df(record, num_hits_per_otu)
                        for record in all_otus_hits]
        # Case there are no hits in any OTUs at all #
        if all(df is None for df in all_otus_dfs):
            return pandas.DataFrame()
        # Combine them all into a big dataframe #
        df = pandas.concat(all_otus_dfs)
        # Keep the hit number as a seperate column #
        df = df.reset_index()
        col = df.pop("index")
        df.insert(2, 'hit_number', col)
        df['hit_number'] += 1
        # Sort by abundance but keep OTUs grouped of course #
        df = df.sort_values(by=['abundance', 'otu_id', 'hit_number'],
                            ascending=[False, False, True])
        # Return #
        return df

    def one_otu_df(self, record, max_hits=5):
        """
        Returns a dataframe object summarizing the top hits for a given OTU.

        The record object is of type: <Bio.Blast.Record.Blast>

        The title looks like: "gi|631252759|ref|NR_113957.1| Staphylococcus..."
        The query looks like: "centroid=small_1:0;seqs=4233"

        To calculate the percentage identity see:

            https://github.com/peterjc/galaxy_blast/blob/master/tools/
            ncbi_blast_plus/blastxml_to_tabular.py#L230
        """
        # Case there are no hits #
        if not record.descriptions: return
        # Get OTU name and its abundance #
        pattern = '\Acentroid=(.+);seqs=([0-9]+)\Z'
        otu, abund = re.findall(pattern, record.query)[0]
        # Cast to integer #
        abund = int(abund)
        # The description field contains lots of info we don't need #
        shorten = lambda s: s.split('|')[-1].strip()
        # Generate one line per hit #
        def one_row_per_hit():
            for align, desc in zip(record.alignments, record.descriptions):
                # Calculate all parameters #
                score       = desc.score
                e_value     = desc.e
                acc         = desc.accession
                title       = shorten(desc.title)
                hsps        = align.hsps[0]
                ident       = hsps.identities
                # The lengths available #
                qry_length  = record.query_length
                hit_length  = align.length
                algn_length = hsps.align_length
                # The ratios defined #
                p_cover = float(qry_length) / float(algn_length)
                p_ident = float(ident)      / float(algn_length)
                # Format percentages as strings #
                p_cover = "%0.1f%%" % (100 * p_cover)
                p_ident = "%0.1f%%" % (100 * p_ident)
                # Yield one row #
                yield abund, otu, title, score, e_value, acc,\
                      algn_length, p_cover, p_ident
        # Convert to dataframe #
        df = pandas.DataFrame(one_row_per_hit(), columns=self.columns)
        # Sort by E-value and also by score #
        df = df.sort_values(by=['e_value', 'score'], ascending=False)
        # Filter and take only the top hits #
        df = df.head(max_hits)
        # Return #
        return df

###############################################################################
class BlastResults:

    def __init__(self, tax):
        self.tax = tax

    @property_cached
    def report(self):
        """
        The NcbiBlastReport object is used to produce a PDF document containing
        all the information that concern this particular taxonomic
        classification strategy.
        """
        # Imports #
        from pacmill.reports.ncbi_blast import NcbiBlastReport
        # Output directory #
        path = self.tax.proj.autopaths.taxonomy_dir + 'reports/'
        path = path + self.tax.short_name + '.pdf'
        # Return #
        return NcbiBlastReport(self.tax, path)
