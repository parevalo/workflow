#!/bin/bash -l
# Basic version of mosaicking script for a single date. A more advanced version
# will be needed to specify the order of the rasters according the the scenes
# with the highest number of images
# MAGNA SIRGAS is EPSG:4686
# UTM18-19N are EPSG:32618-32619

cd /projectnb/landsat/projects/Colombia/Mosaics/M1B

module load gdal/1.10.0

# Quick mosaic, doesnt consider order of images. Gdalwarp is actually using 
# the 8 cores and it speeds up the mosaicking drastically!

for yr in $(seq -w 02 02); do
    qsub -pe omp 8 -V -N mosaic_$yr -j y -b y gdalwarp --config \
    GDAL_CACHEMAX 500 -wm 500 -multi -co NBITS=4 -wt Byte -t_srs EPSG:4686 \
    -srcnodata 0 \
    "/projectnb/landsat/projects/Colombia/images/*/Results/M1/Class/ClassM1B_2002.tif" \
    20$yr"_MAGNA.tif"
done
