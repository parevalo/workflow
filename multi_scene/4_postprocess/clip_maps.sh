#!/bin/bash -l

# Script to cut the classification maps to the WRS2 boundary, because there
# are extra areas with less data and lower classification quality
# caused by the bigger footprint of L8 images.

# List of scenes to be processed

scn_list="007059" #003058 003059 004057 004058 004059 004061 004062 005057 005058 \
#          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
#          007060 007061 008058 008059 008060 009059 009060"

# General settings

module load gdal/1.11.1
poly=/projectnb/landsat/projects/Colombia/vector/WRS2_amazon_selection.shp
pref=mergedmaps_20
suf=_postsieved.tif

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # cd to the corresponding class folder
    cd /projectnb/landsat/projects/Colombia/images/$s/Results/M3/Class
    
    
    # Submit the clipping job
    for yr in $(seq -w 01 01); do
#        qsub -j y -V -N clip$pt$rw"-"$yr -b y \
         gdalwarp -tr 30 30 -srcnodata 0 -cutline $poly -cl WRS2_amazon_selection \
          -cwhere "'PTRW=$pt$rw'" -ot Byte -co NBITS=4 -overwrite $pref$yr$suf \
           mergedmaps_20$yr"_final.tif"
    done
done

