#!/bin/bash -l
#$ -V
#$ -b y
#$ -j y
#$ -N edit_chgmap

# This script UNSETS the NoData information from the changemap and mergedmaps
# layers, making the # -9999 and 0, respectively,  usable in gdalcalc, 
# and also avoiding the problem of NoData areas in the calculations. 
# REQUIRES GDAL >= 2.1

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings

rootdir=/projectnb/landsat/projects/Colombia/images
pre=numchange_2001-2016_

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    cd $rootdir/$s/Results/M3/Class
        
    # Unset nodata in changemap
    gdal_edit.py -unsetnodata $pre$pt$rw".tif"

    # Unset nodata in mergedmaps
    for yr in $(seq -w 01 16); do
        gdal_edit.py -unsetnodata mergedmaps_20$yr"-01-01.tif"
    done
done

