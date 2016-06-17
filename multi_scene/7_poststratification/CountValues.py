#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Written by Eric Bullock
# Modified by Paulo Arevalo


"""HTML Cleaner.

Usage:
  CountValues.py <filename>


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

    # count unique values for the given band
  
    flatraster = raster.flatten()
    stats = np.bincount(flatraster)
    max = flatraster.max()
    min = flatraster.min()

# Print class and number of pixels
for i in (range(min, max+1)):
    print "Class {0}: {1}".format(i, stats[i])
