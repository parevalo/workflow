#!/usr/bin/env python

"""Area calculation with standard errors. For now it assumes that input table 
   has row and column labels, and that the last four columns are samples per 
   class, number of map pixel per class, area in ha and class weigth.

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
    """
    Calculates area proportions and their standard error, as well as estimated
    area and it's confidence interval
    """
    # Convert to numpy array to make it easier to index
    table = df.values

    # Get number of rows, assumed to be number of classes, and column labels
    numclasses = table.shape[0]
    headers = df.dtypes.index[0:numclasses]
    
    # Calculate total area
    total_area = np.sum(table[:, 14])

    # Initialize empty arrays for area proportions and standard error
    p = np.zeros((numclasses, numclasses))
    se = np.zeros((numclasses))

    # Calculate fixed terms of the equation. 
    ni = table[:, 12]
    w = table[:, 15]
   
    # Calculate area proportions and standard error of them
    for i in range(numclasses):
        # Iterate over columns (e.g. reference samples)
        p[:, i] = w * table[:, i] / ni 
        se[i] = np.sqrt(np.sum((w * p[:, i] - p[:, i]**2) / (ni  - 1)))
   
    # Calculate estimated area by class and its confidence interval

    class_area = p.sum(axis=0) * total_area
    class_ci = se * 1.96 * total_area


    return (se, p, class_area, class_ci)

def main():
    csvfile = arguments['<csv>']
    if not os.path.exists(csvfile):
            print 'Could not find input file {0}'.format(csvfile)
            sys.exit(1)

    # Read CSV
    df = read_csv(csvfile)

    # Calculate areas
    se, p, class_area, class_ci = calculate_areas(df)
    
    # Print/export values
    print "-------> Estimated area with lower and upper 95% CI"
    for i in range(se.shape[0]):
        print df.dtypes.index[i] + " = {:0.2f} {:0.2f} {:0.2f}".format(\
                                    class_area[i], \
                                    class_area[i] - class_ci[i], \
                                    class_area[i] + class_ci[i])

if __name__ == '__main__':
    arguments = docopt(__doc__)
    main()

