#!/bin/bash

# This script runs the MAPPING script for a given list of scenes.
# It also automates finding the example image, which is the first 
# Landsat 7 image found in the folder

# TODO
# Make this and the other files more verbose as a way to keep track of the
# workflow

# List of scenes to be processed

scn_list="004057 004058 004061 004062 005057 005058 005059 005060 \
          005061 006058 006059 006060 006061 007058 007059 007060 \
          008058 008059 008060 009059"

# General setting: path to template, root dir, etc

rootdir=/projectnb/landsat/projects/Colombia/images

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # Set scene, TS(images) and results paths for the map script 
    scn_path=$rootdir/$s
    ts_path=$scn_path/images
    res_path=$scn_path/Results/M1/TSR

    # Find example image
    img=$(find $ts_path -maxdepth 1 -type d -name "*LE7*" | head -1 )
    example_img=$(basename $img)
    
    # Set variables for map script
    img_path=$ts_path/$example_img/$example_img"_stack" 
    class_day=-01-01
 
    # CD to Class folder
    cd /projectnb/landsat/projects/Colombia/images/$s/Results/M1/Class

    # Run map script for multiple dates
    for yr in $(seq -w 0 15); do    
        qsub -j y -V -N map_$pt$rw"-"$yr -b y \
         yatsm -v map --root $ts_path --result $res_path --image $img_path\
          class 20$yr$class_day ClassM1A_20$yr$class_date".tif"
    done
done 
