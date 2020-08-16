# `pacmill` version 0.2.8

The `pacmill` python package is a bioinformatics pipeline that is developed to process microbial 16S amplicon sequencing data. It is specialized in the analysis of long reads such as those provided by PacBio sequencers.

## Prerequisites

Since `pacmill` is written in python, it is compatible with all operating systems: Linux, macOS and Windows. The only prerequisite is `python3` (which is often installed by default) along with the `pip3` package manager.

To check if you have `python3` installed, type the following on your terminal:

    $ python3 -V

If you do not have `python3` installed, please refer to the section [obtaining python3](docs/installing_tips.md#obtaining-python3).

To check if you have `pip3` installed, type the following on your terminal:

    $ pip3 -V

If you do not have `pip3` installed, please refer to the section [obtaining pip3](docs/installing_tips.md#obtaining-pip3).

## Installing

To install the `pacmill` package, simply type the following commands on your terminal:

    $ pip3 install --user pacmill

Alternatively, if you want to install it for all users of the system:

    $ sudo pip3 install pacmill

These commands will also automatically install all the other python modules on which `pacmill` depends.

## Usage

### Metadata

The first thing to do when starting a new analysis is to fill in a metadata file that details all there is to know about the samples being processed.

A template for such a file is found under this repository at "./metadata/metadata_blank.xlsx". In addition, another file called "./metadata/metadata_example.xlsx" shows typical values that the fields are supposed to take.

### Loading the project

Bellow are some examples to illustrate the various ways there are to use this package.

    # This example is not completed yet. TODO.

### Customizing report headers

To change the text that appears inside the header of the PDF reports generated, you can adjust these three environment variables to your liking. Credit is appreciated where credit is due, but the software has a very permissive license that lets you decide what is best.

    $ export PACMILL_HEADER="From the \textbf{pacmill} project"
    $ export PACMILL_SUBHEADER="Written by consultants at \url{www.sinclair.bio}""
    $ export PACMILL_LINK="Hosted at \url{www.github.com/xapple/pacmill}"

## Demo project

To run the demo project, do the following:

    # This example is not completed yet. TODO.

## Flowchart

<p align="center">
<img src="docs/flowchart.png?raw=true">
</p>

## Extra documentation

More documentation is available at:

<http://xapple.github.io/pacmill/pacmill>

This documentation is simply generated from the source code with:

    $ pdoc --html --output-dir docs --force pacmill