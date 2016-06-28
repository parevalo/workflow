#!/usr/bin/env python

"""Area calculation with standard errors. Assumes input table has row and 
   column labels.   

Usage:
    estimate_area.py <csv>
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

    # Get number of rows, assumed to be number of classes, and column labels
    numclasses = table.shape[0]
    headers = df.dtypes.index[0:numclasses]
    
    # Calculate total area
    total_area = np.sum(table[:, 14]

    # Calculate terms of the equation. 
    S = np.zeros((numclasses))
    ni = table[:, 12]
    w = table[:, 15]
   
    for i in range(numclasses):
        # Iterate over columns 
        p = w * table[:, i] / ni 
        S[i] = np.sqrt(np.sum((w * p - p**2) / (ni  - 1)))
   
    print headers
    print S*1.96*100

    return S

def main():
    csvfile = arguments['<csv>']
    if not os.path.exists(csvfile):
            print 'Could not find input file {0}'.format(csvfile)
            sys.exit(1)

    # Read CSV
    df = read_csv(csvfile)

    # Calculate areas
    S = calculate_areas(df)

if __name__ == '__main__':
    arguments = docopt(__doc__)
    main()

