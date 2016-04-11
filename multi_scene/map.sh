#!/bin/bash

# This script runs the MAPPING script for a given list of scenes.
# It also automates finding the example image, which is the first 
# Landsat 7 image found in the folder

# TODO
# Make this and the other files more verbose as a way to keep track of the
# workflow

# List of scenes to be processed

scn_list="005058 006058"

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
    res_path=$scn_path/Results/FIT1/TSR

    # Find example image
    img=$(find $ts_path -maxdepth 1 -type d -name "*LE7*" | head -1 )
    example_img=$(basename $img)
    
    # Set variables for map script
    img_path=$ts_path/$example_img/$example_img"_stack" 
    class_date=2001-01-31
 
    # CD to Class folder
    cd /projectnb/landsat/projects/Colombia/images/$s/Results/FIT1/Class

    # Run map script
    qsub -j y -V -N map_$pt$rw -b y \
     yatsm -v map --root $ts_path --result $res_path --image $img_path\
      class $class_date Class3_"$class_date".tif    
    
done 
