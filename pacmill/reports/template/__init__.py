#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os

# Internal modules #
from pymarktex.templates.sinclair_bio import HeaderTemplate, FooterTemplate

###############################################################################
class Header(HeaderTemplate):
    """All the parameters to be rendered in the LaTeX header template."""

    orig_header    = r"From the \textbf{pacmill} project"
    orig_subheader = r"Written by consultants at \url{www.sinclair.bio}"
    orig_link      = r"Hosted at \url{www.github.com/xapple/pacmill}"

    def header(self):
        """This text appears first in the header, on the left hand side."""
        return os.environ.get('PACMILL_HEADER', self.orig_header)

    def subheader(self):
        """This text appears second in the header, on the left hand side."""
        return os.environ.get('PACMILL_SUBHEADER', self.orig_subheader)

    def link(self):
        """This text appears third in the header, on the left hand side."""
        return os.environ.get('PACMILL_LINK', self.orig_link)

    def title(self):
        """This text appears first in the header, on the right hand side."""
        return self.options.get('title', 'Auto-generated report')

###############################################################################
class Footer(FooterTemplate):
    """All the parameters to be rendered in the LaTeX footer template."""
    pass