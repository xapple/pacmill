#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

This script will download the raw input data for the demonstration project
of the `pacmill` pipeline.

Typically you would run this file from a command line like this:

    ipython3 -i -- ~/deploy/pacmill/demo/download.py
"""

# Built-in modules #

# Internal modules #

# Constants #

###############################################################################
if __name__ == "__main__":
    # Import #
    from pacmill.demo.demo_samples import samples
    # For debugging purposes #
    if False: samples = samples[3:4]
    # Message #
    print("Downloading the 5 demo samples which are:\n")
    print('\n  * '.join(sample.description for sample in samples), '\n')
    # Download #
    for sample in samples: print(sample.download())
    # Extract #
    for sample in samples: print(sample.extract())
    # Compress #
    print("Compressing FASTQ files.")
    for sample in samples: print(sample.fastq.compress(remove_orig=True))
    # Remove the sra archive #
    for sample in samples: sample.path_sra.remove()
    # Success #
    print("Done.")
