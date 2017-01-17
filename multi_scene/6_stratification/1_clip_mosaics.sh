#!/bin/bash -l

# Script to clip mosaic to Amazon boundary. 

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

# Create mask (this was done manually, but recreate script here)
# Created by reclassifying the cropped map of 2001 that had been created with
# the previous code that used a vector as cutline. Maybe recreate that code here:
# use cutline to cut 2001, reclassify and use that as the mask, for completeness
# sake. Mask has values of 15 for nodata and 1 in the valid region. 

# Do clipping with a mask in order to cut to the amazon boundary AND
# separate the true NoData from the pixels with a break and no label afterwards
# (unclassified), bc we will assume those are mostly pastures. This will be
# also reflected in the stratification. This is done this way in order to 
# be able to create the annual change maps instead of cummulative maps 
# (e.g. 2004-2005 instead of 2001-2005).

for yr in $(seq -w 01 16); do
    qsub -j y -V -N clip$yr -b y \
     gdal_calc.py -A 20$yr"_final.tif" -B "true_data_mask.tif" \
      --outfile=20$yr"_final_crop.tif" \
      --calc='"logical_and(A >= 0 , B == 0)*15 + logical_and(A >=0, B == 1)*A"' \
      --type=Byte --co="COMPRESS=PACKBITS" --overwrite

done




