#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import uuid

# First party modules #
from fasta import FASTQ
from autopaths.file_path import FilePath
from autopaths.tmp_path  import new_temp_dir
from plumbing.cache      import property_cached
from plumbing.check_cmd_found import check_cmd
from plumbing.scraping   import download_from_url

# Third party modules #
import sh

###############################################################################
class DumpSRA:
    """
    Takes care of running the `sra-toolkit` program to extract a FASTQ from
    an SRA file. See:

    https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

    Until Ubuntu 19, `sra-toolkit` was a package in the apt universe. This
    is not the case for Ubuntu 20. See:

    https://askubuntu.com/questions/1232028/sra-toolkit-for-ubuntu-20-04-tls

    If you are on macOS you can just type "brew install sratoolkit".
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest=None):
        # Source is the SRA on which the sra-toolkit will be run #
        self.source = FilePath(source)
        # Destination is a FASTQ file that will contain the results #
        if dest is None: dest = self.source.prefix_path + '.fastq'
        self.dest = FASTQ(dest)

    #---------------------------- Installing ---------------------------------#
    url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/" \
          "sratoolkit.2.10.8-ubuntu64.tar.gz"

    @classmethod
    def check_installed(cls, exception=True):
        """
        Try to determine if the `sra-toolkit` software is installed and
        accessible.
        """
        return check_cmd('fastq-dump', exception, cls.install.__doc__)

    @classmethod
    def install(cls, prefix="~/programs/sra-toolkit/"):
        """
        To automatically download and install the `sra-toolkit` software on
        this computer and for the current user, type these commands in python:

            >>> from pacmill.demo.sra import DumpSRA
            >>> DumpSRA.install()
        """
        # Make a temporary directory #
        tmp_dir = new_temp_dir()
        # Download tarball #
        zip_loc = download_from_url(cls.url, tmp_dir, stream=True,
                                    progress=True)
        # Uncompress #
        zip_loc.untargz_to(tmp_dir)
        src_dir = tmp_dir.sub_directory
        src_dir.move_to(prefix)
        # The directory that contains the executable #
        bin_dir = src_dir.with_tilda[:-1].replace('~', '$HOME')
        # Mandatory configuration #
        cls.vdb_config_workaround()
        # Suggest adding to the $PATH #
        print("\nThe sra-toolkit was installed successfully. You should now "
              "add this line to your .bash_profile: \n\n    "
              "export PATH=%s/bin:$PATH\n" % bin_dir)

    @classmethod
    def vdb_config_workaround(cls):
        """
        This method needs to be run to avoid an error when running any of
        the `sra-toolkit` tools. This is due to some very strange and
        nonsensical change introduced in version 2.10.3
        See: https://github.com/ncbi/sra-tools/issues/291
        """
        # The directory where the configuration file must be stored #
        ncbi_settings = FilePath("~/.ncbi/user-settings.mkfg")
        ncbi_settings.directory.create_if_not_exists()
        # Generate a new random ID #
        new_guid = uuid.uuid4()
        # Add this to the file #
        ncbi_settings.write('/LIBS/GUID = "%s"\n' % new_guid)

    #------------------------------ Running ----------------------------------#
    def __call__(self, cpus=None, verbose=True):
        # Message #
        if verbose: print("Running fastq-dump on '%s'" % self.source)
        # Check it is installed #
        self.check_installed()
        # Some questionable change they did in recent versions #
        self.vdb_config_workaround()
        # Get the command #
        dump = sh.Command("fastq-dump")
        # Make a temporary directory #
        tmp_dir = new_temp_dir()
        # Run it #
        dump('--outdir', tmp_dir, self.source)
        # Move files #
        result = tmp_dir + self.source.prefix + '.fastq'
        result.move_to(self.dest)
        # Return #
        return self.dest

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the `sra-toolkit` software was run already and the
        results are stored on the filesystem. Return False if it was not yet
        run.
        """
        return self.dest.exists

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from `sra-toolkit` " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return self.dest
