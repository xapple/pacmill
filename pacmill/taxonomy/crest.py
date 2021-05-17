#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, sys, multiprocessing, shutil, re

# First party modules #
from fasta import FASTA
from autopaths.dir_path       import DirectoryPath
from autopaths.file_path      import FilePath
from autopaths.tmp_path       import new_temp_dir
from plumbing.cache           import property_cached
from plumbing.apt_pkg         import get_apt_packages
from plumbing.check_cmd_found import check_cmd
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

    We have lowered the "drop-off from the highest bitscore" to half a percent
    by adding the "-r 0.5" parameter. This should increase the resolution for
    long reads such as the ones we have.
    """

    # Constants #
    short_name = 'crest'
    long_name  = 'LCAClassifier/CREST version 3.2.0'
    executable = 'classify'
    hard_path  = '~/programs/crest/bin/classify'

    # The database #
    db_version_name = "silvamod128"
    db_short_name   = "silvamod"
    db_long_name    = 'Silva version 128 modified for CREST'

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest_dir):
        # Source is the FASTA file containing OTU consensus sequences #
        self.source = FASTA(source)
        # Destination is a directory that contains all the results #
        self.dest_dir = DirectoryPath(dest_dir)

    @property
    def rank_names(self):
        return ['Domain',         # 1
                'Superkingdom',   # 2
                'Kingdom',        # 3
                'Phylum',         # 4
                'Class',          # 5
                'Order',          # 6
                'Family',         # 7
                'Genus',          # 8
                'Species']        # 9

    @property_cached
    def database(self):
        # Find where the database is on the file system #
        path  = os.path.expanduser("~/programs/crest/parts/flatdb/")
        path += '%s/%s.fasta' % (self.db_short_name, self.db_version_name)
        # Make an object with a path and rank_names #
        attributes = {'path':       path,
                      'rank_names': self.rank_names,
                      'tag':        self.db_version_name,
                      'long_name':  self.db_long_name}
        return type('CrestDatabase', (), attributes)

    #----------------------------- Installing --------------------------------#
    apt_packages = ['ncbi-blast+', 'python2', 'python-dev']
    tgz_url = "https://github.com/lanzen/CREST/releases/download/3.2.2/" \
              "LCAClassifierV3.2.2.tar.gz"

    @classmethod
    def check_installed(cls, exception=True):
        """
        Try to determine if the CREST software is installed and
        accessible.
        """
        if not FilePath(cls.hard_path).exists and exception:
            msg = "The executable '%s' is required for this operation." \
                  " Unfortunately it cannot be found."
            msg = msg % cls.hard_path + '\n\n' + cls.install.__doc__
            raise Exception(msg)

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
        print(cls.tgz_url)
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
        # Modify the buildout configuration #
        buildout_config = src_dir + 'buildout.cfg'
        old_line = "#setuptools=44.0.0"
        new_line = "setuptools=44.0.0"
        buildout_config.replace_line(old_line, new_line)
        # Call the install command - step one #
        print("Calling bootstrap command")
        bootstrap_script = src_dir + 'bootstrap.py'
        sh.python2(bootstrap_script,
                   '--setuptools-version=44.0.0',
                   '--buildout-version=2.12.0',
                   _out=sys.stdout,
                   _err=sys.stderr)
        # Call the install command - step two #
        print("Calling buildout command")
        buildout_cmd = src_dir + 'bin/buildout'
        sh.python2(buildout_cmd,
                   _out=sys.stdout,
                   _err=sys.stderr)
        # Restore current directory #
        os.chdir(current_dir)
        # Success message #
        print("\nCREST was installed successfully.")

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /blast/db_hits.xml
                /blast/stdout.txt
                /blast/stderr.txt
                /output/
                /results/stdout.txt
                /results/stderr.txt
                /results/assignments.tsv
                /taxa_tables/
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
    def __call__(self, cpus=None, verbose=True):
        # Message #
        if verbose:
            message = "Running `%s` taxonomy classification on '%s'"
            print(message % (self.short_name, self.source))
        # Check blast is installed #
        check_cmd('blastn', True)
        # Check crest is installed #
        self.check_installed()
        # Check crest is at the right location #
        crest = sh.Command(self.hard_path)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run #
        sh.blastn('-task',             'megablast',
                  '-num_threads',      cpus,
                  '-query',            self.source,
                  '-db',               self.database.path,
                  '-out',              self.autopaths.db_hits,
                  '-max_target_seqs',  '100',
                  '-outfmt',           '5',
                  _out = self.autopaths.blast_stdout,
                  _err = self.autopaths.blast_stderr)
        # Check #
        if os.path.getsize(self.autopaths.db_hits) == 0:
            msg = "Hits file empty. The MEGABLAST process was probably killed."
            raise Exception(msg)
        # Remove directory otherwise crest complains #
        self.autopaths.output_dir.remove()
        # Run algorithm #
        crest('--verbose',
              '-o', self.dest_dir + 'output/',
              '-d', self.db_version_name,
              '-r', 0.5,
              self.autopaths.db_hits,
              _out = self.autopaths.results_stdout.path,
              _err = self.autopaths.results_stderr.path)
        # Get generated files #
        assignments = self.autopaths.output_dir + 'otus.csv'
        # Move files into place #
        shutil.move(assignments, self.autopaths.results_assignments)
        # Clean up #
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0:
            os.remove("error.log")

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
        return CrestResults(self.autopaths, self.database)

###############################################################################
class CrestResults:

    def __init__(self, autopaths, database):
        self.autopaths = autopaths
        self.database = database

    @property_cached
    def assignments(self):
        """
        Typically an assignment from CREST looks like this:

            Main genome;Bacteria;Bacteria (superkingdom);Terrabacteria;
            Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus

        We skip the first entry as it seems to be always "Main genome"
        and classified as rank "meta".
        """
        result = {}
        with open(self.autopaths.assignments, 'r') as handle:
            for line in handle:
                # Split by tabs #
                code, count, species = line.strip('\n').split('\t')
                # First line is just the column titles #
                if species == "classification": continue
                # Split by semi-colons #
                species = tuple(species.strip('\n').split(';'))
                # Get the original OTU name #
                pattern  = r':centroid=(.+);seqs=[0-9]+\Z'
                otu_name = re.findall(pattern, code)[0]
                # Conform to the mothur standard of names #
                otu_name = otu_name.replace(':', '_')
                # Make a dictionary #
                result[otu_name] = species[1:]
        # Return #
        return result

    @property_cached
    def count_unassigned(self):
        """Will count how many did not get a prediction at each level."""
        # Load the assignment values created by CREST #
        values = list(self.assignments.values())
        # Iterate over every rank number #
        rank_numbers = list(range(len(self.database.rank_names)))
        # Calculate for each position in the tree of life #
        result = [sum(1 for x in values if len(x) <= i) for i in rank_numbers]
        # Return #
        return result

    @property_cached
    def count_assigned(self):
        """Will count how many did get a prediction at each level."""
        return [len(self.assignments) - x for x in self.count_unassigned]