#!/bin/bash
# This script submits the jobs to cache the data for all of the scenes

# List of scenes to be processed
#scn_list="004058 004061 004062 005057 005058 005059 005060 \
#          005061 006058 006059 006060 006061 007058 007059 007060 \
#          008058 008059 008060 009059" #004057

scn_list="006061 005060 008060 009059"
# General settings: path to template, root dir, num of jobs 

yconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/yatsm_config.yaml
export ROOTDIR=/projectnb/landsat/projects/Colombia/images
njob=100

# Iterate over scenes
for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}  
    rw=${s:4:2}
    
    # Export all relevant variables for the yaml file
    export INPUT=$ROOTDIR/$s/$pt$rw"_input.csv" 
    export IMG=$ROOTDIR/$s/images
    
    # CD to log folder for that path-row
    cd /projectnb/landsat/projects/Colombia/logs/$pt$rw
 
    # Run yatsm
    for job in $(seq 1 $njob); do
        qsub -j y -V -N cache$pt$rw"_"$job -b y \
         yatsm -v cache $yconfig $job $njob
    done
    
    # For debugging purposes only
    #qsub -j y -V -b y yatsm -v cache --update "yatsm_r/*" $yconfig 1 $njob
  
done
