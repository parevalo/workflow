#!/bin/bash
# Basic version of mosaicking script for a single date. A more advanced version
# will be needed to specify the order of the rasters according the the scenes
# with the highest number of images
# MAGNA SIRGAS is EPSG:4686
# UTM18-19N are EPSG:32618-32619

cd /projectnb/landsat/projects/Colombia/Mosaics/M1A

#module load gdal/1.10.0, not working!

qsub -pe omp 4 -V -j y -b y gdalwarp -multi -co NBITS=4 -wt Byte -t_srs EPSG:4686 -srcnodata 0  \
 "/projectnb/landsat/projects/Colombia/images/*/Results/M1/Class/ClassM1A_2008.tif" 2008_MAGNA_SIRGAS.tif

