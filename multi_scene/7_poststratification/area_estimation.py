#!/usr/bin/env python

"""Area calculation with standard errors

Usage:
    estimate_area.py <error matrix csv)
"""
from docopt import docopt

import pandas as pd
import os
import sys
import numpy as np

# Make stdout unbuffered
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

def read_csv(csvfile):
    """
    Reads CSV file with error matrix as numpy array
    """
    # Read in file
    with open(csvfile, 'rb') as f:
        csvreader = csv.reader(f)
        if header:
            csvreader.next()
        errormatrix = { int(row[0]) : int(row[1]) for row in csvreader }


return lut


