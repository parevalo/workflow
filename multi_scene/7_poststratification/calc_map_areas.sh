#!/bin/bash

# Calculate class areas for all the strata maps
spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene/7_poststratification 

# Find all strata maps
cd /projectnb/landsat/projects/Colombia/Mosaics/M3

for i in $(ls final_strata*.tif); do
    yr=${i:13:5}
    qsub -j y -b y -V -N str_area_$yr -l mem_per_core=8G\
    $spath/CountValues.py $i "strata_"$yr"_pixcount.csv"

done

