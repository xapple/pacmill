#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import socket

# Internal modules #

# First party modules #
from autopaths.dir_path import DirectoryPath
from plumbing.cache     import property_cached

# Third party modules #

###############################################################################
class Bundle:
    """
    A bundle object will regroup various result files and PDF reports from
    one project into a single zip file, useful for delivery
    and distribution.

    Once the bundle is ready, you can simply download it to your local
    computer with a simple rsync command.
    """

    def __init__(self, parent, base_dir, archive=None):
        # Keep a reference to the parent Project object #
        self.parent  = parent
        self.project = parent
        # Where the results will be aggregated #
        self.base_dir = DirectoryPath(base_dir)
        # Where the zip archive will be placed #
        if archive is None:
            archive = self.base_dir.path[:-1] + '.zip'
        self.archive = archive

    #-------------------------- Automatic paths ------------------------------#
    all_paths = """
                /projects/report.pdf
                /samples/
                /taxonomies/
                /metadata/samples.xlsx
                """

    @property_cached
    def autopaths(self):
        """
        The AutoPaths object is used for quickly assessing the filesystem paths
        of various file inputs/outputs and directories.
        See https://github.com/xapple/autopaths#autopaths-object
        """
        from autopaths.auto_paths import AutoPaths
        return AutoPaths(self.base_dir, self.all_paths)

    #------------------------------ Running ----------------------------------#
    def __call__(self):
        # Remove the directory if it was created previously #
        self.base_dir.remove()
        self.base_dir.create()
        # Loop every sample report #
        for s in self.project:
            destination = self.autopaths.samples_dir + s.short_name + '.pdf'
            s.report.output_path.copy(destination)
        # Check for any early exit #
        if not self.project.otus: return
        # Report for project #
        self.project.report.output_path.copy(self.autopaths.report)
        # Report for taxonomies #
        for r in self.project.taxonomy.reports.all:
            if not r.tax.should_run: continue
            r.output_path.copy(self.autopaths.taxonomies_dir)
        # Zip it #
        self.base_dir.zip_to(self.archive)
        # Remove the directory #
        self.base_dir.remove()
        # Return #
        return self.archive

    #------------------------------- Results ---------------------------------#
    @property_cached
    def results(self):
        results = BundleResults(self)
        message = "You can't access results from a bundle before" \
                  "creating it."
        if not results: raise Exception(message)
        return results

###############################################################################
class BundleResults:

    def __nonzero__(self): return bool(self.autopaths.x)

    def __init__(self, parent):
        self.parent    = parent
        self.base_dir  = parent.base_dir
        self.autopaths = parent.autopaths

    @property_cached
    def rsync(self):
        """
        Will return an rsync bash command as a string that can be later used
        to download the bundle.
        """
        # Get the hostname #
        host = socket.gethostname()
        # Get the path to the gzip archive #
        gzip = self.parent.archive
        # Get the destination name #
        name = self.parent.project.short_name + '.zip'
        # Make the command #
        cmd = 'rsync -avz %s:%s ~/Downloads/pacmill/%s'
        cmd = cmd % (host, gzip, name)
        # Return #
        return cmd