#!/bin/bash -l

# Script to split the categorical and prediction probability bands into
# separate files.  

# List of scenes to be processed

scn_list="006059" #004058 004061 004062 005057 005058 005059 005060 \
          #005061 006058 006059 006060 006061 007058 007059 007060 \
          #008058 008059 008060 009059" #004057

# General settings

module load gdal/1.11.1
dt="-01-01"

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # cd to the corresponding class folder
    cd /projectnb/landsat/projects/Colombia/images/$s/Results/M2/Class
    
    # Submit the clipping job
    for yr in $(seq -w 00 01); do
        qsub -j y -V -N split$pt$rw"-"$yr -b y \
         gdal_translate -b 1 -ot Byte -co NBITS=4 ClassM2B_20$yr$dt"_crop.tif" \
         ClassM2B_20$yr$dt"_crop_class.tif"
    done
done

