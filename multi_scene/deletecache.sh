#!/bin/bash
#$ -V
#$ -j y

# This script is meant to delete the old cache files in order to create
# the new ones with the recently downloaded images. I tried to use the
# --update option but zsh threw and error because of the pattern and I
# couldn't figure it out, hence this script. 

# List of scenes to be processed
#scn_list="004061 004062 005057 005058 005059 005060 \
#          005061 006058 006059 006060 006061 007058 007059 007060 \
#          008058 008059 008060 009059" #004057 004058 
scn_list="006061 005060 008060 009059"

ROOTDIR=/projectnb/landsat/projects/Colombia/images

# Iterate over scenes
for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}  
    rw=${s:4:2}
    
    # Set cache folder and delete files
    cachefolder=$ROOTDIR/$s/images/.cache
    cd $cachefolder
    rm -f yatsm_r*

done

