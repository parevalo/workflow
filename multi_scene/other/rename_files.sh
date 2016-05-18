#!/bin/bash
#$ -V
#$ -j y
#$ -N rename_files

# Script to rename files (mostly classification results) if needed

# List of scenes to be processed

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings
pref=ClassM3_20
dt="-01-01"

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # cd to the corresponding class folder
    cd /projectnb/landsat/projects/Colombia/images/$s/Results/M3/Class
    
    # Rename the files, be verbose
    for yr in $(seq -w 01 15); do
        mv -v $pref$yr$dt".tif" $pref$yr$dt"_M3train.tif"
    done

done

