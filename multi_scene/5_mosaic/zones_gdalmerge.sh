#!/bin/bash -l
# Merge west and east zones, both in UTM18

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal
#module load python/2.7.5_nopath
#module load gdal/1.11.1

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

qsub -V -N mergezones -j y -b y \
 gdal_merge.py -o finalmergeexample.tif -co COMPRESS=PACKBITS -co NBITS=4 \
 -ot Byte -n 0 -v eastUTM18tapnewextent.tif westUTM18_2016.tif


