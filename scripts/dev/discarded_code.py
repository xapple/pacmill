#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Placeholder script to contain pieces of code that are not used anymore
or discarded. This script can safely be ignored.
"""

# Built-in modules #

# Internal modules #

# Constants #

###############################################################################
def all_metadata(self):
   # Initialize result #
   result = "\n\\begin{itemize}\n"
   # Format them as a bullet list #
   value =  lambda k: str(getattr(self.sample, k))
   bullet = lambda k: r'\item \texttt{%s} \texttt{%s}' % (k, value(k))
   # Join every item #
   result += '\n'.join(bullet(k) for k in keys)
   # Close the result text #
   result += "\n\\end{itemize}\n"
   # Return #
   return result