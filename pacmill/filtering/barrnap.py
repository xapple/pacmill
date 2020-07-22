#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, multiprocessing

# First party modules #
from fasta import FASTQ
from autopaths.file_path      import FilePath
from autopaths.tmp_path       import new_temp_dir, new_temp_path
from plumbing.cache           import property_cached
from plumbing.check_cmd_found import check_cmd
from plumbing.apt_pkg         import get_apt_packages
from plumbing.scraping        import download_from_url

# Third party modules #
import sh, tag

###############################################################################
class Barrnap:
    """
    Takes care of running the Barrnap program on a given FASTQ file.
    See https://github.com/tseemann/barrnap
    Barrnap predicts the presence and location of ribosomal RNA genes in
    genomes. Expects version 0.9.

    The output of the program is a GFF3 and looks something like:

        P.marinus  barrnap:0.9  rRNA  353314  354793  0       +  .  Name=16S_rRNA
        P.marinus  barrnap:0.9  rRNA  355464  358334  0       +  .  Name=23S_rRNA
        P.marinus  barrnap:0.9  rRNA  358433  358536  9.6e-07 +  .  Name=5S_rRNA
    """

    # Proportional length threshold to reject prediction #
    reject_threshold = 0.4

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest=None, filtered=None):
        # Source is the FASTA or FASTQ file on which barrnap will be run #
        self.source = FilePath(source)
        # Destination is a GFF file that contains the results #
        if dest is None:
            dest = self.source.prefix_path + '.barrnap.gff'
        self.dest = FilePath(dest)
        # Filtered is a FASTA file containing only raw reads with rRNA genes #
        if filtered is None:
            filtered = self.source.prefix_path + '.barrnap.fastq'
        self.filtered = FASTQ(filtered)

    #---------------------------- Installing ---------------------------------#
    apt_packages = ['bedtools']
    zip_url = "https://github.com/tseemann/barrnap/archive/0.9.zip"

    @classmethod
    def check_installed(cls, exception=True):
        """
        Try to determine if the Barrnap software is installed and
        accessible.
        """
        return check_cmd('barrnap', exception, cls.install.__doc__)

    @classmethod
    def install(cls, prefix="~/programs/barrnap/"):
        """
        To automatically download and install the Barrnap software on this
        computer and for the current user, type these commands in python:

            >>> from pacmill.filtering.barrnap import Barrnap
            >>> Barrnap.install()
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

    #------------------------------ Running ----------------------------------#
    def __call__(self, cpus=None, verbose=True):
        # Message #
        if verbose: print("Running barrnap on '%s'" % self.source)
        # Check it is installed #
        self.check_installed()
        # Check version #
        assert b"0.9" in sh.barrnap('--version').stderr
        # If the input is a FASTQ we need to make a FASTA first #
        if self.source.endswith('fastq'):
            source = new_temp_path(suffix='.fasta')
            FASTQ(self.source).to_fasta(source)
        else: source = self.source
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run it #
        sh.barrnap(source, '--threads', cpus, '--reject',
                   self.reject_threshold, _out=str(self.dest))
        # Return #
        return self.dest

    def filter(self, verbose=True):
        """
        Using the original reads file and the GFF output of barrnap,
        we will create a new FASTQ file containing only the original reads
        that had a hit for both rRNA genes.
        """
        # Message #
        if verbose:
            msg = "Extracting sequences with both rRNA genes from '%s'"
            print(msg % self.dest)
        # Get IDs that passed (only for 16S) #
        reader = tag.GFF3Reader(infilename=self.dest)
        reader = tag.select.features(reader, type='rRNA')
        ids_16s = [rec.seqid for rec in reader if '16S_rRNA' in rec.attributes]
        # Get IDs that passed (only for 23S) #
        reader = tag.GFF3Reader(infilename=self.dest)
        reader = tag.select.features(reader, type='rRNA')
        ids_23s = [rec.seqid for rec in reader if '23S_rRNA' in rec.attributes]
        # Get those that had both genes #
        ids_both = set(ids_16s) & set(ids_23s)
        # Extract those IDs #
        FASTQ(self.source).extract_sequences(ids_both, self.filtered, verbose)
        # Return #
        return self.filtered

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the Barrnap software was run already and the results are
        stored on the filesystem. Return False if it was not yet run.
        """
        return self.filtered.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from Barrnap " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return BarrnapResults(self.filtered)

###############################################################################
class BarrnapResults(FASTQ):
    """A file with the results from Barrnap."""
    pass