#!/bin/bash

# This script runs the training using the merged cache file. 
# It's been modified from the cache_training, so many things were left
# unchanged

# We can do this with one scene only 

scn=004057

# General setting: path to template, root dir, etc

yconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/yatsm_config.yaml
rfconfig=/projectnb/landsat/projects/Colombia/workflow/multi_scene/RandomForest.yaml
train_dir=/projectnb/landsat/projects/Colombia/Training
export ROOTDIR=/projectnb/landsat/projects/Colombia/images

# Get path and row in short version
pt=${s:2:1}
rw=${s:4:2}

# Export all relevant variables for the yaml file. Note that variables like
# TRAINING are needed even if they're not used

export INPUT=$ROOTDIR/$scn/$pt$rw"_input.csv"
export RESULTS=$ROOTDIR/$scn/Results/M3/TSR
export IMG=$ROOTDIR/$scn/images
export TRAINING=$train_dir/$scn/Training1.tif
export TRAINCACHE=$train_dir/train_cached/M3/M3_full_traincache.npz

# CD to classifiers folder
cd /projectnb/landsat/projects/Colombia/classifiers/M3

# Run training, verify number correspond to training raster number
qsub -j y -V -N trainM3_full -b y \
 yatsm -v train --diagnostics $yconfig $rfconfig M3_fulltrain_noweights.pkl

 
