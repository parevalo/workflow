#!/bin/bash

# Script to create the samples for each of the biannual maps (including
# a buffer class (16) around the forest to pasture areas. Ran in interactive
# mode with stdout-stderr redirection bc it wasn't working with qsub.

# Set paths and file names
cd /projectnb/landsat/projects/Colombia/biannual_samples_june2017

inpath=/projectnb/landsat/projects/Colombia/Mosaics/M3
spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene/6_stratification/helper_scripts
seed_base=10000

# Set year and step info
first_yr=1
last_yr=16
step=2
fy=$(printf %02d $first_yr)
ly=$(printf %02d `expr $last_yr - $step`)


# Submit jobs for each biannual map
for yr in $(seq -w $fy $step $ly); do

    # Set filenames
    yr2=$(printf %02d `expr $yr + $step`)
    bname="buffered3_final_strata_annual_"$yr"_"$yr2"_UTM18N.tif"
    seedval=$(expr $seed_base + $yr2)

    # Run script 
    $spath/sample_map.py -v --size 1050 \
     --allocation "50 400 75 50 75 50 50 50 50 50 50 50 50" \
      --mask "15 255" --ndv 255 --vector sample_$yr"_"$yr2".shp" \
       --seed_val $seedval stratified $inpath/$bname &>> biannual_samples_out.txt
done

