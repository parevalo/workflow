#!/bin/bash
# This script submits the jobs to generate the CSV list of the dates and full
# filepaths for all of the scenes

# List of scenes to be processed
scn_list="004057 004058 004061 004062 005057 005058 005059 005060 \
          005061 006058 006059 006060 006061 007058 007059 007060 \
          008058 008059 008060 009059" 

ROOTDIR=/projectnb/landsat/projects/Colombia/images
scriptpath=/projectnb/landsat/users/parevalo/yatsm/scripts

# Iterate over scenes
for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}  
    rw=${s:4:2}
    
    # Set variables for image folder and csv names
    OUTPUT=$ROOTDIR/$s/$pt$rw"_input.csv" 
    IMG=$ROOTDIR/$s/images
    
    # Run the script
    $scriptpath/gen_date_file.sh -o -v $IMG $OUTPUT

done

