#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

This script will download the raw input data for the demonstration project
of the `pacmill` pipeline.

Typically you would run this file from a command line like this:

    ipython3 -i -- ~/deploy/pacmill/pacmill/demo/download.py
"""

# Built-in modules #

# Internal modules #

# Constants #

###############################################################################
if __name__ == "__main__":
    # Import #
    from pacmill.demo.demo_samples import samples
    # Check we have the required program for uncompressing first #
    from pacmill.demo.sra import DumpSRA
    DumpSRA.check_installed()
    # Message #
    print("\n Downloading the 5 demo samples which are:")
    texts = (sample.description for sample in samples)
    print('\n  * ' + '\n  * '.join(texts), '\n')
    print("This will take a few minutes.")
    # Download #
    for sample in samples: print(sample.download())
    # Extract #
    for sample in samples: print(sample.extract())
    # Compress #
    print("Compressing all FASTQ files.")
    for sample in samples: print(sample.fastq.compress(remove_orig=True))
    # Remove the SRA archive #
    for sample in samples: sample.path_sra.remove()
    # Success #
    print("Done. Results are in '%s'." % samples[0].base_dir)
