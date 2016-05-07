#!/bin/bash -l

# Script to cut the classification maps to the WRS2 boundary, because there
# are extra areas with less data and lower classification quality
# caused by the bigger footprint of L8 images.

module load gdal/1.11.1
cd /projectnb/landsat/projects/Colombia/images/004058/Results/M1/Class
poly=/projectnb/landsat/projects/Colombia/vector/WRS2_amazon_selection.shp

gdalwarp -tr 30 30 -srcnodata 0 -cutline $poly -cl WRS2_amazon_selection \
 -cwhere 'PTRW=458' ClassM1B_2001.tif Cutexample.tif

