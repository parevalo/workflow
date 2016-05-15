#!/bin/bash -l

# Module to sieve the mosaics. Will be changed to sieve individual images,
# given that the mosaic sieving is not working correctly (job never finishes/
# takes too long). 

module load python/2.7.5_nopath
module load gdal/1.11.3

cd /projectnb/landsat/projects/Colombia/images/007059/Results/M2/Class 
export GDAL_CACHEMAX=2048

for yr in $(seq -w 01 01); do
    qsub -V -N sieve_$yr -j y -b y \
    gdal_sieve.py -st 8 -8 ClassM2B_20$yr"-01-01_crop_class.tif" -mask 2001_2-3.tif  sievetest.tif
done

