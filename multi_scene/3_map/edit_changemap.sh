#!/bin/bash -l
#$ -V
#$ -b y
#$ -j y
#$ -N edit_chgmap

# This script changes the nodata value of the changemap to 0, in order to
# be able to use the -9999 as an actual number to make a comparison in gdalcalc 

module load python/2.7.5_nopath
module load gdal/1.11.1

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings

rootdir=/projectnb/landsat/projects/Colombia/images
pre=numchange_2001-2015_

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    cd $rootdir/$s/Results/M3/Class
    
    # Modify nodata from changemap so that we can use it in gdalcalc
    gdal_edit.py -a_nodata 0 $pre$pt$rw".tif"
        
done

