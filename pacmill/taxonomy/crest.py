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
from fasta import FASTA
from autopaths.dir_path       import DirectoryPath
from autopaths.tmp_path       import new_temp_dir
from plumbing.cache           import property_cached
from plumbing.check_cmd_found import check_cmd
from plumbing.apt_pkg         import get_apt_packages
from plumbing.scraping        import download_from_url

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

    # Constants #
    short_name = 'crest'
    long_name  = 'LCAClassifier/CREST version 3.2.0'
    executable = 'classify'

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest_dir):
        # Source is the FASTA file containing OTU consensus sequences #
        self.source = FASTA(source)
        # Destination is a directory that contains all the results #
        self.dest_dir = DirectoryPath(dest_dir)

    @property_cached
    def database(self):
        # Find where the database is on the file system #
        self.database_path  = which('classify').physical_path.directory.directory
        self.database_path += 'parts/flatdb/%s/%s.fasta' % (database, database)

    #----------------------------- Installing --------------------------------#
    apt_packages = ['ncbi-blast+', 'python2']
    tgz_url = "https://github.com/lanzen/CREST/releases/download/3.2.0/" \
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
        tgz_loc = download_from_url(cls.tgz_url, tmp_dir, stream=True,
                                    progress=True)
        # Uncompress #
        tgz_loc.untargz_to(tmp_dir)
        src_dir = tmp_dir.sub_directory
        # Move it to where the user wants it #
        src_dir.move_to(prefix)
        # This 'bootstrap' stuff can't detect paths it seems, have to cwd #
        current_dir = os.getcwd()
        os.chdir(src_dir.path)
        # CREST needs to be updated, until then the next commands won't work #
        print(1/0)
        # Call the install command - step one #
        bootstrap_script = src_dir + 'bootstrap.py'
        sh.python2(bootstrap_script)
        # Call the install command - step two #
        buildout_cmd = src_dir + 'bin/buildout'
        sh.python2(buildout_cmd)
        # Restore current directory #
        os.chdir(current_dir)
        # The directory that contains the executable #
        bin_dir = src_dir.with_tilda[:-1].replace('~', '$HOME')
        # Suggest adding to the $PATH #
        print("\nCREST was installed successfully. You should now "
              "add this line to your .bash_profile: \n\n    "
              "export PATH=%s/bin:$PATH\n" % bin_dir)

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /graphs/
                /stats/
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
    def __call__(self, cpus=None):
        # Check blast is installed #
        check_cmd('blastn', True)
        # Check crest is installed #
        check_cmd('classify', True)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run #
        sh.blastn('-task',             'megablast',
                  '-num_threads',      cpus,
                  '-query',            self.source,
                  '-db',               self.database_path,
                  '-out',              self.autopaths.db_hits,
                  '-max_target_seqs',  '100',
                  '-outfmt',           '5',
                  _out = self.autopaths.blast_stdout,
                  _err = self.autopaths.blast_stderr)
        # Check #
        if os.path.getsize(self.autopaths.db_hits) == 0:
            msg = "Hits file empty. The MEGABLAST process was probably killed."
            raise Exception(msg)
        # Remove directory #
        self.autopaths.crest_dir.remove()
        # Run algorithm #
        sh.classify('--verbose',
                    '--rdp',
                    '-o', self.base_dir + 'crest/',
                    '-d', self.database,
                    self.autopaths.db_hits,
                    _out = self.autopaths.crest_stdout.path,
                    _err = self.autopaths.crest_stderr.path)
        # Move #
        shutil.move(self.autopaths.db_hits.prefix_path + '_Composition.tsv', self.autopaths.crest_composition)
        shutil.move(self.autopaths.db_hits.prefix_path + '_Tree.txt', self.autopaths.crest_tree)
        shutil.move(self.autopaths.db_hits.prefix_path + '_Assignments.tsv', self.autopaths.crest_assignments)
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