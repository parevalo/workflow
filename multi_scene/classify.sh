#!/bin/bash

# This script runs the classification script for a given list of scenes.
# It is NOT designed to change the model parameters inside the YAML file,
# as they are not taken into account for this step

# List of scenes to be processed

scn_list="005058 006058"

# General setting: path to template, root dir, etc

yconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/yatsm_config.yaml
algopath=/projectnb/landsat/projects/Colombia/classifiers
export ROOTDIR=/projectnb/landsat/projects/Colombia/images
njob=200

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
    export TRAINING=$ROOTDIR/$s/images/Training1.tif
    
    # Test for scenes and change the start and end train date accordingly   
    # THIS SECTION IS PROBABLY UNNECESSARY 
    if [ $s = "005058" ]; then 
        export TRAINSTART="2000-078"
        export TRAINEND="2001-032"
    elif [ $s = "006058" ]; then 
        export TRAINSTART="2001-031"
        export TRAINEND="2003-013"
    else
        echo "Date in conditions don't match those in the list"
        exit
    fi

    # CD to classifiers folder
    cd /projectnb/landsat/projects/Colombia/logs/$pt$rw
    

    # Run classification, verify algorithm being used
    for job in $(seq 1 $njob); do
        qsub -j y -V -N class$pt$rw"_"$job -b y \
         yatsm -v classify $yconfig $algopath/mergedtrain_859-658-558.pkl $job $njob 
         # yatsm -v classify $yconfig $algopath/$pt$rw"_trainRF1.pkl" $job $njob 
        
    done
done 
