#!/bin/bash -l

# Script to clip mosaic to Amazon boundary. 

cd /projectnb/landsat/projects/Colombia/Mosaics/M3
limit=/projectnb/landsat/projects/Colombia/vector/amazon_boundary.shp

# Use cutline to cut 2001 and reclassify to generate a mask of valid area that
# we will use for the rest of the dates. Mask has values of 0 for nodata 
# and 1 in the valid region. 

# Do clipping
clippedraster=2001_final_cutline_cropped.tif
maskname=true_data_mask.tif

qsub -j y -V -N cutlineclip -b y  gdalwarp -cutline $limit -cl amazon_boundary \
 -ot Byte -co COMPRESS=PACKBITS -overwrite "2001_final.tif" \
    $clippedraster

# Create mask

qsub -j y -V -N createmask -b y -hold_jid cutlineclip \
 gdal_calc.py -A $clippedraster --outfile=$maskname \
  --calc='"(A==0)*0 + (A>0)*1"' --type=Byte --co="COMPRESS=PACKBITS" --overwrite 

# Do clipping with a mask in order to cut to the amazon boundary AND
# separate the true NoData from the pixels with a break and no label afterwards
# (unclassified), bc we will assume those are mostly pastures. This will be
# also reflected in the stratification. This is done this way in order to 
# be able to create the annual change maps instead of cummulative maps 
# (e.g. 2004-2005 instead of 2001-2005).

for yr in $(seq -w 01 16); do
    qsub -j y -V -N clip$yr -b y -hold_jib createmask \
     gdal_calc.py -A 20$yr"_final.tif" -B $maskname \
      --outfile=20$yr"_final_crop.tif" \
      --calc='"logical_and(A >= 0 , B == 0)*15 + logical_and(A >=0, B == 1)*A"' \
      --type=Byte --co="COMPRESS=PACKBITS" --overwrite

done




