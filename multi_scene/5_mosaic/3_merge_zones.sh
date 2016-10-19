#!/bin/bash -l
# Merge west and east zones, both in UTM18

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

for yr in $(seq -w 01 16); do
    qsub -V -N mrgezones_$yr -j y -b y \
     gdal_merge.py -o 20$yr"_final.tif" -co COMPRESS=PACKBITS -co NBITS=4 \
     -ot Byte -n 0 -v eastUTM18_20$yr".tif" westUTM18_20$yr".tif"
done

