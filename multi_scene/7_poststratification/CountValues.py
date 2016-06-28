#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Written by Eric Bullock
# Modified by Paulo Arevalo


"""Script to count the number of pixels per class in a given raster file. 
   Outputs a csv with the class numbers and pixel count.
   Current version only really checks the first band. 

Usage:
  CountValues.py <filename> <output>


"""

from docopt import docopt
import gdal, ogr, osr, os
import numpy as np
from osgeo import gdal
import sys
import math


if __name__ == '__main__':
    args = docopt(__doc__, version='0.1.0')


path = args['<filename>']
out_csv = args['<output>']

 
gdalData = gdal.Open(path)
if gdalData is None:
  sys.exit( "ERROR: can't open raster" )
 
# get width and heights of the raster
xsize = gdalData.RasterXSize
ysize = gdalData.RasterYSize
 
# get number of bands
bands = gdalData.RasterCount
 
# process the raster
for i in xrange(1, bands + 1):
    band_i = gdalData.GetRasterBand(i)
    raster = band_i.ReadAsArray()

    # count unique values for the given band (missing values get a count of 0)
     
    flatraster = raster.flatten()
    stats = np.bincount(flatraster)
    max = flatraster.max()
    min = flatraster.min()
    outarray = np.zeros((max+1, 2))
    outarray[:, 0] = np.arange(max+1)
    outarray[:, 1] = stats


# Print/export class and number of pixels
for i in (range(min, max+1)):
    print "Class {0}: {1}".format(i, stats[i])

np.savetxt(out_csv, outarray, delimiter=',', header="class,pixels", fmt='%01d')

