#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import multiprocessing

# First party modules #
from fasta import FASTA, FASTQ
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

        P.marinus  barrnap:0.9  rRNA  3314  4793  0       +  .  Name=16S_rRNA
        P.marinus  barrnap:0.9  rRNA  5464  8334  0       +  .  Name=23S_rRNA
        P.marinus  barrnap:0.9  rRNA  8433  8536  9.6e-07 +  .  Name=5S_rRNA

    You can either use barrnap in filter mode, or in extract mode.
    In addition, there is a "remove ITS" mode too.
    See the three child classes below.
    """

    # Proportional length threshold to reject prediction #
    reject_threshold = 0.4

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest=None, filtered=None):
        # Source is the FASTA or FASTQ file on which barrnap will be run #
        if source.endswith('fasta'): self.source = FASTA(source)
        else:                        self.source = FASTQ(source)
        # Destination is a GFF file that contains the results #
        if dest is None:
            dest = self.source.prefix_path + '.barrnap.gff'
        self.dest = FilePath(dest)
        # Filtered is a file containing only raw reads with rRNA genes #
        if filtered is None:
            extension = self.source.extension
            filtered  = self.source.prefix_path + '.barrnap' + extension
        # It can be either a FASTA or a FASTQ #
        if filtered.endswith('fasta'): self.filtered = FASTA(filtered)
        else:                          self.filtered = FASTQ(filtered)

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

        If you are on macOS you can just type:
            $ brew install brewsci/bio/barrnap
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
    def run(self, cpus=None, verbose=True):
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

    #-------------------------------- Parsing --------------------------------#
    def parse_hit_ids(self, attr_text):
        """
        Using the GFF output of barrnap and the `tag` parser, we retrieve
        all IDs that contained a given text in their attributes.
        """
        reader = tag.GFF3Reader(infilename=self.dest)
        reader = tag.select.features(reader, type='rRNA')
        return [rec.seqid for rec in reader if attr_text in rec.attributes]

    @property_cached
    def ids_16s(self):
        """Return read IDs for reads that contained a 16S gene."""
        return frozenset(self.parse_hit_ids('16S_rRNA'))

    @property_cached
    def ids_23s(self):
        """Return read IDs for reads that contained a 23S gene."""
        return frozenset(self.parse_hit_ids('23S_rRNA'))

    def parse_hit_location(self, attr_text):
        """
        Using the GFF output of barrnap and the `tag` parser, we retrieve
        the start and end position of rRNA genes.
        """
        # Get reads that had a rRNA hit #
        reader = tag.GFF3Reader(infilename=self.dest)
        reader = tag.select.features(reader, type='rRNA')
        # Make a lookup table for those with hits #
        return {rec.seqid: (rec.start, rec.end)
                for rec in reader if attr_text in rec.attributes}

    @property_cached
    def loc_16s(self):
        """Return start and end location of reads that contained a 16S gene."""
        return self.parse_hit_location('16S_rRNA')

    @property_cached
    def loc_23s(self):
        """Return start and end location of reads that contained a 23S gene."""
        return self.parse_hit_location('23S_rRNA')

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
        return self.filtered

###############################################################################
class BarrnapFilter(Barrnap):
    """
    With this subclass, we will filter reads after running barrnap.

    Using the original reads file and the GFF output of barrnap,
    we will create a new FASTQ file containing only the original reads
    that had a hit for both rRNA genes on the same read.

    The sequences themselves however are not modified.
    """

    def __call__(self, cpus=None, verbose=True):
        # Call the parent class method #
        self.run(cpus, verbose)
        # Extra processing #
        return self.filter()

    #--------------------------- Presence of gene ----------------------------#
    def filter(self, verbose=True):
        # Message #
        if verbose:
            msg = "Extracting sequences with both rRNA genes from '%s'"
            print(msg % self.dest)
        # Get those that had both genes #
        ids_both = self.ids_16s & self.ids_23s
        # Extract those IDs #
        FASTQ(self.source).extract_sequences(ids_both, self.filtered)
        # Return #
        return self.filtered

###############################################################################
class BarrnapExtract(Barrnap):
    """
    With this subclass, we will extract reads regions after running barrnap.

    Using the original reads file and the GFF output of barrnap, we will create
    a new FASTA file containing only the portion of the original reads that
    contains the 16S rRNA gene.
    """

    def __call__(self, cpus=None, verbose=True):
        # Call the parent class method #
        self.run(cpus, verbose)
        # Extra processing #
        return self.extract()

    #--------------------------- Location of gene ----------------------------#
    def extract(self, verbose=True):
        # Message #
        if verbose:
            msg = "Extracting the 16S rRNA portion of sequences from '%s'"
            print(msg % self.dest)
        # Function to yield only the good part of each read #
        def only_16s_portion(reads):
            for r in reads:
                start_and_end = self.loc_16s.get(r.id)
                if start_and_end is None: continue
                start, end = start_and_end
                yield r[start:end]
        # Write new FASTA file #
        self.filtered.write(only_16s_portion(self.source))
        # Return #
        return self.filtered

###############################################################################
class RemoveITS(Barrnap):
    """
    With this subclass, we will extract the 16S region and the 23S region
    from every read that had both, and concatenate them one next to
    each other, effectively removing the ITS regions.
    """

    def __call__(self, cpus=None, verbose=True):
        # Call the parent class method #
        self.run(cpus, verbose)
        # Extra processing #
        return self.remove_its()

    #------------------------ Location of both genes -------------------------#
    def remove_its(self, verbose=True):
        # Message #
        if verbose:
            msg = "Removing the ITS portion of sequences from '%s'"
            print(msg % self.dest)
        # Function to yield concatenated read #
        def concat_16s_23s(reads):
            for r in reads:
                # Retrieve positions #
                loc_16s = self.loc_16s.get(r.id)
                loc_23s = self.loc_23s.get(r.id)
                # Skip this read if no 16S found #
                if loc_16s is None: continue
                # Get only the 16S part in a new sequence #
                start, end = loc_16s
                seq = r[start:end]
                # Add the 23S part to sequence if it was found #
                if loc_23s is not None:
                    start, end = loc_23s
                    seq += r[start:end]
                # Return the newly created sequence #
                yield seq
        # Write new FASTA file #
        self.filtered.write(concat_16s_23s(self.source))
        # Return #
        return self.filtered
