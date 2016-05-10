#!/bin/bash -l

# Module to sieve the mosaics

module load python/2.7.5_nopath
module load gdal/1.11.1

cd /projectnb/landsat/projects/Colombia/Mosaics/M2B
export GDAL_CACHEMAX=2048

for yr in $(seq -w 01 01); do
    qsub -V -N sieve_$yr -j y -b y \
    gdal_sieve.py -st 8 -8  20$yr"-01-01_seq.tif" -mask Mask3B_2001_2-3.tif sievetest2.tif
done

