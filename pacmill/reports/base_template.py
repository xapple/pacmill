#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os
import socket

# Internal modules #
import pacmill

# First party modules #
from pymarktex.templates import Template
from plumbing.common     import pretty_now

###############################################################################
class ReportTemplate(Template):
    """Things that are common to most reports in pacmill."""

    def project_name(self):      return pacmill.project_name
    def project_url(self):       return pacmill.project_url
    def project_version(self):   return pacmill.__version__
    def now(self):               return pretty_now()

    def hostname(self):
        host = os.environ.get('PACMILL_HOSTNAME')
        if host is not None: return host
        return socket.gethostname()

    def git(self):
        if not pacmill.git_repo: return False
        return {'git_hash'  : pacmill.git_repo.hash}
