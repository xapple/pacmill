#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Special variables #
__version__ = '0.5.3'

# Built-in modules #
import os, sys

# First party modules #
from autopaths    import Path
from plumbing.git import GitRepo

# Constants #
project_name = 'pacmill'
project_url  = 'https://github.com/xapple/pacmill'

# Get paths to module #
self       = sys.modules[__name__]
module_dir = Path(os.path.dirname(self.__file__))

# The repository directory #
repos_dir = module_dir.directory

# The module is maybe in a git repository #
git_repo = GitRepo(repos_dir, empty=True)

# If we are in development mode, check dependencies #
setup_py = repos_dir + 'setup.py'
if setup_py.exists:
    from plumbing.dependencies import check_setup_py
    check_setup_py(setup_py)
