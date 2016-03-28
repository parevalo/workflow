#!/bin/bash
# This script takes the specified folders in the DB list and submits
# the corresponding yatsm line jobs to be run. Change FIT folder if necessary.

# List of scenes to be processed
scn_list="006058" #004058 006060 007059

# General settings: path to template, root dir, num of jobs 

yconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/yatsm_config.yaml
export ROOTDIR=/projectnb/landsat/projects/Colombia/images
njob=400

# Iterate over scenes
for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}  
    rw=${s:4:2}
    
    # Export all relevant variables for the yaml file
    export INPUT=$ROOTDIR/$s/$pt$rw"_input.csv" 
    export RESULTS=$ROOTDIR/$s/Results/FIT1/TSR
    export CONFIG=$ROOTDIR/$s/Results/FIT1
    export IMG=$ROOTDIR/$s/images
    #export TRAINING=To be determined
    
    # CD to log folder for that path-row
    cd /projectnb/landsat/projects/Colombia/logs/$pt$rw
 
    # Run yatsm
    for job in $(seq 1 $njob); do
        qsub -j y -V -N y$pt$rw"_"$job -b y \
         yatsm -v line --check_cache --resume $yconfig $job $njob
    done
       
done
