#!/bin/bash

# Calculate class areas for all the comparison maps

cd /projectnb/landsat/projects/Colombia/Mosaics/M3/IDEAM
spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene/7_poststratification 

period_list="03-05 05-07 07-09 09-11 11-13 13-14 14-15"

for p in $period_list; do
    qsub -j y -b y -V -N compare_area_$p -l mem_per_core=8G\
    $spath/CountValues.py comparison_$p".tif" "comparison_"$p"_pixcount.csv"

done

