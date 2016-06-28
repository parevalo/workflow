#!/usr/bin/env python

"""Area calculation with standard errors. Assumes input table has row and 
   column labels.   

Usage:
    estimate_area.py <error matrix csv)
"""
#TODO:
#- Add back multiple options for specifying row columns/labels
#- Add options to determine what to export, e.g. area proportions, CI, etc

from docopt import docopt

import pandas as pd
import os
import sys
import numpy as np

# Make stdout unbuffered
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

def read_csv(csvfile):
    """
    Reads CSV file with error matrix as pandas dataframe
    """
    # Read in file
    df = pd.read_csv(csvfile, sep=',', header=0, index_col=0)

    return df


def calculate_areas(df):
    # Convert to numpy array to make it easier
    table = df.values

    # Get number of rows, assumed to be number of classes
    numclasses = table.shape[0]

    # Calculate terms of the equation. 
    S = np.zeros((numclasses))
    ni = table[:, 12]
    w = table[:, 15]
    p = w * table[:,0:12] / ni
   
    for i in range(numclasses):
        # Iterate over columns (Attepmt again using removing p line below and indexing directly)
        p = w * table[:, i] / ni #Delete if line below works
        S[i] = np.sqrt(np.sum((w * p[:, i] - p[:, i]**2) / (ni  - 1)))
		
    	print S

    return S
