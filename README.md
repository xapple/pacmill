[![PyPI version](https://badge.fury.io/py/pacmill.svg)](https://badge.fury.io/py/pacmill)

# `pacmill` version 0.3.5

The `pacmill` python package is a bioinformatics pipeline that is developed to process microbial 16S amplicon sequencing data. It is specialized in the analysis of long reads such as those provided by PacBio sequencers.

## Prerequisites

Since `pacmill` is written in python, it is compatible with all operating systems: Linux, macOS and Windows. The only prerequisite is `python3` (which is often installed by default) along with the `pip3` package manager.

To check if you have `python3` installed, type the following on your terminal:

    $ python3 -V

If you do not have `python3` installed, please refer to the section [obtaining python3](docs/markdown/installing_tips.md#obtaining-python3).

To check if you have `pip3` installed, type the following on your terminal:

    $ pip3 -V

If you do not have `pip3` installed, please refer to the section [obtaining pip3](docs/markdown/installing_tips.md#obtaining-pip3).

## Installing

To install the `pacmill` package, simply type the following commands on your terminal:

    $ pip3 install --user pacmill

Alternatively, if you want to install it for all users of the system:

    $ sudo pip3 install pacmill

These commands will also automatically install all the other python modules on which `pacmill` depends.

## External programs

The `pacmill` pipeline also depends on several shell commands being available. The following executables should be present in your `$PATH` environment variable.

* `fastQValidator`, `fastqc`, `barrnap`, `vsearch`, `mothur`, `xelatex`, `fastq-dump`

If any of these required external programs are missing, you will be prompted to install them and given easy instructions to do so.

## Usage

### Metadata

The first thing to do when starting a new analysis is to fill in a metadata file that details all there is to know about the samples being processed.

An empty template for such a file is found under this repository at "./metadata/metadata_blank.xlsx". You can make a copy of this file for every new project.
 
 In addition, another file called "./metadata/metadata_example.xlsx" shows typical values that the fields are supposed to take along with a short documentation for each entry. A excerpt of this file is shown below:
 
 <p align="center">
 <img src="docs/images/metadata_screenshot.png?raw=true">
 </p>

### Loading the project

Bellow are some examples to illustrate the various ways there are to use this package.

    # This example is not completed yet. TODO.

### Customizing report headers

To change the text that appears inside the header of the PDF reports generated, you can adjust these three environment variables to your liking. Credit is appreciated where credit is due, but the software has a very permissive license that lets you decide what is best.

    $ export PACMILL_HEADER="From the \textbf{pacmill} project"
    $ export PACMILL_SUBHEADER="Written by consultants at \url{www.sinclair.bio}""
    $ export PACMILL_LINK="Hosted at \url{www.github.com/xapple/pacmill}"

## Demo project

In order to test and evaluate the pipeline, we have provided a demonstration project ready to be processed. This enables the user to see what type of outputs are generated by `pacmill` without having to bring his own sequence data. Five samples are included and are taken from the following publication:

* "Confident phylogenetic identification of uncultured prokaryotes through long read amplicon sequencing of the 16S-ITS-23S rRNA operon.''
* Joran Martijn, [many others], Thijs Ettema.
* Science for Life Laboratory, Uppsala University
* https://doi.org/10.1111/1462-2920.14636

The samples are publicly accessible on SRA and are described as follows:

* `mock`: Genomic DNA from 38 phylogenetically distinct and diverse bacteria and archaea.
* `p19`: Sediment sample obtained from hot spring Radiata Pool, Ngatamariki, New Zealand.
* `pm3`: Sediment sample taken from 1.25m below the sea floor using a gravity core at Aarhus Bay, Denmark.
* `sala`: Black biofilm that was taken at 60m depth in an old silver mine near Sala, Sweden.
* `tns08`: Sediment sample taken from a shallow submarine hydrothermal vent field near Taketomi Island, Japan.

To run the demo project, do the following:

    # This example is not completed yet. TODO.

## Example graphs 

The `pacmill` pipeline produces a multitude of graphs and visualizations after having processed the sequence data. Below are two examples. Firstly a sequence length distribution of cleaned reads. Secondly, a barstack of taxonomic assignments for five different samples at the phylum level.

<p align="center">
<img src="docs/images/demo_len_hist.png?raw=true">
</p>

<p align="center">
<img src="docs/images/demo_barstack.png?raw=true">
</p>

<p align="center">
<img src="docs/images/demo_legend.png?raw=true">
</p>

## Example reports

After running the pipeline on a set of FASTQ files, several PDF reports are auto-generated. Examples of two three reports are given below. The first concerns an individual sample while the second details the results of a project containing several samples. The third focuses on taxonomic assignment results and visualizations.

<p align="center">
<a href="https://xapple.github.io/pacmill/demo_reports/project.pdf" class="image fit" target="_blank">
<img src="docs/images/pdf_icon.png" width="120em">
<p>Project report</p>
</a>

<a href="https://xapple.github.io/pacmill/demo_reports/sample.pdf" class="image fit" target="_blank">
<img src="docs/images/pdf_icon.png" width="120em">
<p>Sample report</p>
</a>

<a href="https://xapple.github.io/pacmill/demo_reports/taxonomy.pdf" class="image fit" target="_blank">
<img src="docs/images/pdf_icon.png" width="120em">
</a>
<p>Taxonomy report</p>
</p>

## Flowchart

Below is presented a flowchart detailing the multiple processing steps that occur in the `pacmill` pipeline in a chronological order.

<p align="center">
<img src="docs/images/flowchart.png?raw=true">
</p>


## Extra documentation

More documentation is available at:

<http://xapple.github.io/pacmill/pacmill>

This documentation is simply generated from the source code with:

    $ pdoc --html --output-dir docs --force pacmill