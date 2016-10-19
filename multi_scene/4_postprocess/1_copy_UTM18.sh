#!/bin/bash -l

# Script to create copies of the maps in UTM18 zone to follow the naming
# convention of the projected files, in order to make mosaicking easier. 

# List of scenes to be processed (those already in UTM18)

scn_list="006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings

pref=mergedmaps_20
suf=_final.tif

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # cd to the corresponding class folder
    cd /projectnb/landsat/projects/Colombia/images/$s/Results/M3/Class
    
    # Create the copies
    for yr in $(seq -w 01 01); do
        cp -v $pref$yr$suf $pref$yr"_final_UTM18N.tif"
    done
done

