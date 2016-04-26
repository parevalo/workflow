#!/bin/bash
# Basic version of mosaicking script for a single date


qsub -V -j y -b y gdalwarp -co BIGTIFF=YES -wt Byte -t_srs EPSG:32619 -srcnodata 0  \
 "/projectnb/landsat/projects/Colombia/images/*/Results/M1/Class/ClassM1A_2000.tif" 2000.tif
