#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from pacmill.core.project         import Project
from pacmill.core.sample          import Sample
from pacmill.filtering.seq_filter import SeqFilter

# First party modules #
from autopaths.dir_path import DirectoryPath

# Third party modules #
import mock, pytest

###############################################################################
@pytest.fixture(scope="module")
def this_script_dir(request):
    """Return the directory of the currently running test script."""
    return DirectoryPath(request.fspath.dirname)

###############################################################################
@pytest.fixture
# Remove the constructor #
@mock.patch.object(Project, '__init__', lambda *args: None)
def project(tmp_path):
    """Returns a fake Project object ready to be used."""
    # Make a new instance #
    project = Project('pytest_project')
    # Set the name #
    project.short_name = 'pytest_project'
    # Return #
    return project

###############################################################################
@pytest.fixture
# Remove some methods #
@mock.patch.object(Sample, 'transform_attrs', lambda *args: None)
@mock.patch.object(Sample, 'validate_attrs',  lambda *args: None)
def sample(project, tmp_path):
    """Returns a fake Sample object ready to be used."""
    # Make a new instance #
    sample = Sample(project)
    # Set the name #
    sample.short_name = 'pytest_sample'
    # Set the output directory #
    sample.output_dir = str(tmp_path)
    # Return #
    return sample

###############################################################################
@pytest.fixture
def seq_filter(sample):
    """Returns a fake SeqFilter object ready to be used."""
    # Make a mock SeqFilter object #
    seq_filter = SeqFilter(sample)
    # Return #
    return seq_filter
