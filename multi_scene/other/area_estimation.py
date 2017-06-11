#!/usr/bin/env python

"""Area calculation with standard errors. For now it assumes that input table 
   has row and column labels, and that the last four columns are samples per 
   class, number of map pixel per class, area in ha and class weigth.

   This is an old script and hasn't been updated in a while. All of these steps
   are being done in the calculate_strate_per_year.R for now.


Usage:
    area_estimation.py <input_csv> <output_csv>
"""
#TODO:
#- Add back multiple options for specifying row columns/labels
#- Add options to determine what to export?
#- Remove the need to specify areas and weights

from docopt import docopt

import pandas as pd
import os
import sys
import numpy as np

# Make stdout unbuffered
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

def read_csv(csv_in):
    """
    Reads CSV file with error matrix as pandas dataframe
    """
    # Read in file
    df = pd.read_csv(csv_in, sep=',', header=0, index_col=0)

    return df


def calculate_areas(df):
    """
    Calculates area proportions and their standard error, as well as estimated
    area and its confidence interval
    """
    # Convert to numpy array to make it easier to index
    table = df.values

    # Get number of rows, assumed to be number of classes, and column labels
    nclass = table.shape[0]
    headers = df.dtypes.index[0:nclass]
    
    # Calculate total area
    total_area = np.sum(table[:, 14])

    # Initialize empty arrays for area proportions and standard error
    p = np.zeros((nclass, nclass))
    se = np.zeros((nclass))

    # Calculate fixed terms of the equation. 
    ni = table[:, 12]
    w = table[:, 15]
   
    # Calculate area proportions and standard error of them
    for i in range(nclass):
        # Iterate over columns (e.g. reference samples)
        p[:, i] = w * table[:, i] / ni 
        se[i] = np.sqrt(np.sum((w * p[:, i] - p[:, i]**2) / (ni  - 1)))
   
    # Calculate estimated area by class and its confidence interval

    class_area = p.sum(axis=0) * total_area
    class_ci = se * 1.96 * total_area


    return (se, p, class_area, class_ci, nclass)

def main():
    # Input csv
    csv_in = arguments['<input_csv>']
    if not os.path.exists(csv_in):
            print 'Could not find input file {}'.format(csv_in)
            sys.exit(1)

    # Output csv
    csv_out = arguments['<output_csv>'] 
    if os.path.dirname(csv_out) == '':
        csv_out = './' + csv_out
    
    # Read CSV
    df = read_csv(csv_in)

    # Calculate areas and CI
    se, p, class_area, class_ci, nclass = calculate_areas(df)
    upper_ci = class_area + class_ci
    lower_ci = class_area - class_ci

    # Print/export values (redundant, remove or change) 

    print "-------> Estimated area with lower and upper 95% CI"
    for i in range(nclass):
        print df.dtypes.index[i] + " = {:0.2f} {:0.2f} {:0.2f}".format(class_area[i], 
                                                                lower_ci[i], 
                                                                upper_ci[i])

    # Save output to csv
    rownames = ['Estimated_area','lower_ci', 'upper_ci']
    out_df = pd.DataFrame(data=[class_area,lower_ci, upper_ci], 
                                columns=df.dtypes.index[0:nclass], 
                                index=rownames)
    out_df.to_csv(csv_out)

if __name__ == '__main__':
    arguments = docopt(__doc__)
    main()

