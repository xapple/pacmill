#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, multiprocessing, shutil

# First party modules #
from plumbing.cache import property_cached
from plumbing.check_cmd_found import check_cmd

# Third party modules #
import sh

###############################################################################
class CrestClassify:
    """
    Predicts taxonomy of DNA sequences by BLASTING against the 'SILVAMOD'
    database and picking the lowest common ancestor amongst contradicting
    assignments in the tree of life.

    For instance if all hits for a given OTU only agree at the order level,
    the assignment stops at the order level.

    CREST stands for Classification Resources for Environmental Sequence Tags.
    It is also known as "LCAClassifier".

    The code is here:        https://github.com/lanzen/CREST
    The publication is here: http://dx.plos.org/10.1371/journal.pone.0049334
    """

    short_name = 'crest'
    long_name  = 'LCAClassifier/CREST version 2.0.4'
    executable = 'classify'

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, centers, result_dir):
        # Parent #
        self.centers    = centers
        self.database   = database
        self.result_dir = result_dir
        # Find where the database is on the file system #
        self.database_path  = which('classify').physical_path.directory.directory
        self.database_path += 'parts/flatdb/%s/%s.fasta' % (database, database)

    #----------------------------- Installing --------------------------------#
    apt_packages = ['ncbi-blast+']
    zip_url = "https://github.com/lanzen/CREST/releases/download/3.2.0/" \
              "LCAClassifierV3.2.0.tar.gz"

    @classmethod
    def check_installed(cls, exception=True):
        """
        Try to determine if the CREST software is installed and
        accessible.
        """
        return check_cmd(cls.executable, exception, cls.install.__doc__)

    @classmethod
    def install(cls, prefix="~/programs/crest/"):
        """
        To automatically download and install the CREST classifier software on
        this computer and for the current user, type these commands in python:

            >>> from pacmill.taxonomy.crest import CrestClassify
            >>> CrestClassify.install()
        """
        # Start with required apt packages #
        get_apt_packages(cls.apt_packages, verbose=True)
        # Make a temporary directory #
        tmp_dir = new_temp_dir()
        # Download tarball #
        zip_loc = download_from_url(cls.zip_url, tmp_dir, stream=True,
                                    progress=True)
        # Uncompress #
        zip_loc.unzip_to(tmp_dir, single=False)
        src_dir = tmp_dir.sub_directory
        src_dir.move_to(prefix)
        # Set executable permissions #
        bin_loc = src_dir + 'bin/barrnap'
        bin_loc.permissions.make_executable()
        bin_loc = src_dir + 'binaries/linux/nhmmer'
        bin_loc.permissions.make_executable()
        # The directory that contains the executable #
        bin_dir = src_dir.with_tilda[:-1].replace('~', '$HOME')
        # Suggest adding to the $PATH #
        print("\nBarrnap was installed successfully. You should now "
              "add this line to your .bash_profile: \n\n    "
              "export PATH=%s/bin:$PATH\n" % bin_dir)

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /otu_table.csv
                /otu_table_norm.csv
                /centers.fasta
                /graphs/
                /stats/
                /comp_phyla/
                /comp_tips/
                /comp_order/
                /comp_class/
                /blast/db_hits.xml
                /blast/stdout.txt
                /blast/stderr.txt
                /crest/assignments.tsv
                /crest/composition.tsv
                /crest/tree.txt
                /crest/Relative_Abundance.tsv
                /crest/Richness.tsv
                /crest/stdout.txt
                /crest/stderr.txt
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
    def run(self, cpus=None):
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run #
        sh.blastn('-task',              'megablast',
                  '-num_threads',       cpus,
                  '-query',             self.centers,
                  '-db',                self.database_path,
                  '-out',               self.p.db_hits,
                  '-max_target_seqs',   '100',
                  '-outfmt',            '5',
                  _out=self.p.blast_stdout, _err=self.p.blast_stderr)
        # Check #
        if os.path.getsize(self.p.db_hits) == 0:
            raise Exception("Hits file empty. The MEGABLAST process was probably killed.")
        # CREST #
        self.p.crest_dir.remove()
        sh.classify('--verbose', '--rdp',
                    '-o', self.base_dir + 'crest/',
                    '-d', self.database,
                    self.p.db_hits,
                    out=self.p.crest_stdout.path, _err=self.p.crest_stderr.path)
        # Move #
        shutil.move(self.p.db_hits.prefix_path + '_Composition.tsv', self.p.crest_composition)
        shutil.move(self.p.db_hits.prefix_path + '_Tree.txt', self.p.crest_tree)
        shutil.move(self.p.db_hits.prefix_path + '_Assignments.tsv', self.p.crest_assignments)
        # Clean up #
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0: os.remove("error.log")

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the taxonomic software was run already and the results
        are stored on the filesystem. Return False if it was not yet run.
        """
        return self.autopaths.assignments.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from the CREST taxonomic" \
                  " classification before running the tool."
            raise Exception(msg)
        # Return the results #
        return CrestResults(self.autopaths)

###############################################################################
class CrestResults(object):

    def __init__(self, autopaths):
        self.autopaths = autopaths

    @property_cached
    def assignments(self):
        result = {}
        with open(self.autopaths.assignments, 'r') as handle:
            for line in handle:
                code, species = line.split('\t')
                result[code] = tuple(species.strip('\n').split(';'))[:8]
        return result

    @property
    def count_assigned(self):
        """How many got a position?"""
        return len([s for s in self.assignments.values() if s != ('No hits',)])