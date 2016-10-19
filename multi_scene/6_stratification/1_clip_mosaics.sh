#!/bin/bash -l

# Script to clip mosaic to Amazon boundary. 

cd /projectnb/landsat/projects/Colombia/Mosaics/M3
limit=/projectnb/landsat/projects/Colombia/vector/amazon_boundary.shp

#Do clipping

for yr in $(seq -w 01 16); do
    qsub -j y -V -N clip$yr -b y \
    gdalwarp -cutline $limit -cl amazon_boundary \
    -ot Byte -co COMPRESS=PACKBITS -overwrite 20$yr"_final.tif" \
    20$yr"_final_crop.tif"
done


