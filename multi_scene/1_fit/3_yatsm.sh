#!/bin/bash
# This script takes the specified folders in the DB list and submits
# the corresponding yatsm line jobs to be run. Change FIT folder if necessary.

# List of scenes to be processed
scn_list="003058 003059 004059 007061 009060" # 004057 004058 004059 004061 004062 005057 005058 \
          #005059 005060 005061 006058 006059 006060 006061 007058 007060 \
          #007061 008058 008059 008060 009059 009060" #007059 006061

# General settings: path to template, root dir, num of jobs 

yconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/yatsm_config.yaml
export ROOTDIR=/projectnb/landsat/projects/Colombia/images
njob=512

# Iterate over scenes
for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}  
    rw=${s:4:2}
    
    # Export all relevant variables for the yaml file
    export INPUT=$ROOTDIR/$s/$pt$rw"_input.csv" 
    export RESULTS=$ROOTDIR/$s/Results/M2/TSR
    export IMG=$ROOTDIR/$s/images
    #export TRAINING=To be determined
    
    # CD to log folder for that path-row
    cd /projectnb/landsat/projects/Colombia/logs/$pt$rw/M2
 
    # Run yatsm
    for job in $(seq 1 $njob); do
        qsub -j y -V -N y$pt$rw"_"$job -b y \
         yatsm -v line --check_cache $yconfig $job $njob
    done
    
    # For debugging purposes only
    #qsub -j y -V -b y -N checkTS \
    # yatsm -v line --check_cache $yconfig 121 $njob
  
done

