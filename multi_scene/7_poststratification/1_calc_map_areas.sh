#!/bin/bash

# Calculate class areas for all the strata maps
spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene/7_poststratification 

# Find all strata maps with a given filename suffix
cd /projectnb/landsat/projects/Colombia/Mosaics/M3
matchstring='buffered3_final_strata_annual*UTM18N.tif'
outprefix="strata_buffered_"
outsuffix="_pixcount.csv"

# Iterate over strata maps and get the year pair to use in the output names
for i in $(ls $matchstring); do
   yr=$(ls $i | grep -Eo "[0-9]{2}_[0-9]{2}")
    qsub -j y -b y -V -N str_area_$yr -l mem_per_core=8G\
     $spath/count_pixels.py $i $outprefix$yr$outsuffix

done

