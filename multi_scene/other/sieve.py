#!/usr/bin/env python

import os
import subprocess
import numpy
import rasterio
from rasterio.features import sieve, shapes

os.chdir('/projectnb/landsat/projects/Colombia/Mosaics/M2B/')

# Register GDAL and OGR drivers.
with rasterio.drivers():
    
    # Read a raster to be sieved.
    with rasterio.open('2001-01-01_seq.tif') as src:
        shade = src.read(1)
    
    # Sieve out features 13 pixels or smaller.
    sieved = sieve(shade, 4, out=numpy.zeros(src.shape, src.dtypes[0]))

    # Print the number of shapes in the sieved raster.
    print("Sieved (13) shapes: %d" % len(list(shapes(sieved))))

    # Write out the sieved raster.
    kwargs = src.meta
    kwargs['transform'] = kwargs.pop('affine')
    with rasterio.open('example-sieved.tif', 'w', **kwargs) as dst:
        dst.write(sieved, indexes=1)

# Dump out gdalinfo's report card and open (or "eog") the TIFF.
#print(subprocess.check_output(
#    ['gdalinfo', '-stats', 'example-sieved.tif']))
#subprocess.call(['open', 'example-sieved.tif'])

