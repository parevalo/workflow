#!/bin/bash -l

module load python/2.7.5_nopath
module load gdal/1.11.1

cd /projectnb/landsat/projects/Colombia/Mosaics/M2B

#gdal_proximity.py 2001-01-01_seq.tif Mask3_2001_2-3.tif -co NBITS=2 \
# -ot Byte -values 2,3 -maxdist 2 -nodata 3 -fixed-buf-val 1 -use_input_nodata YES

gdal_calc.py -A Mask3_2001_2-3.tif --outfile=Mask3B_2001_2-3.tif \
 --calc="logical_or(A == 0, A == 1)" --NoDataValue=3 --type=Byte --co="NBITS=2"

