#!/bin/bash

# This script runs the training script for a given list of scenes in order to
# get the cache files, which then will be merged into a single file. This will
# be used to produce a single pickle used to classify all of the images.

# List of scenes to be processed, arrays of training dates (start and end)

scn_list="004057 005058 005060 006058 006059 006060 007058 007059 008059"
trn_str=(2001-033 2000-078 2000-350 2001-031 2000-005 2000-005 2001-062 2003-004 2001-005) 
trn_end=(2002-004 2001-032 2002-003 2003-013 2001-031 2001-031 2003-004 2004-015 2002-008)

# General setting: path to template, root dir, etc

yconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/yatsm_config.yaml
rfconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/RandomForest.yaml
train_dir=/projectnb/landsat/projects/Colombia/Training
export ROOTDIR=/projectnb/landsat/projects/Colombia/images

# Iterate over scenes with counter for the array
i=0

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # Export all relevant variables for the yaml file
    export INPUT=$ROOTDIR/$s/$pt$rw"_input.csv"
    export RESULTS=$ROOTDIR/$s/Results/M3/TSR
    export IMG=$ROOTDIR/$s/images
    export TRAINING=$train_dir/$s/Training1.tif
    export TRAINCACHE=$train_dir/train_cached/M3/cache$pt$rw".npz"
    
    export TRAINSTART=${trn_str[$i]}
    export TRAINEND=${trn_end[$i]}
    let i+=1

    # CD to classifiers folder
    cd /projectnb/landsat/projects/Colombia/classifiers/M3

    # Run training, verify number correspond to training raster number
    qsub -j y -V -N "trainM3_"$pt$rw -b y \
     yatsm -v train --diagnostics $yconfig $rfconfig M3_train$pt$rw".pkl" 

done 
