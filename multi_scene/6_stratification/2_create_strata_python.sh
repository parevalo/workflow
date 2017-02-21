#!/bin/bash -l

# Script to create the strata from mosaics using python3/rasterio

# CD to folder and define suffix of input tif files

cd /projectnb/landsat/projects/Colombia/Mosaics/M3
suffix="_final_crop.tif"
spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene/

# Define the start and end year we want calculate the strata on

first_yr=7
last_yr=8

# Define if we want the strata to be annual or cumulative (e.g. 2001 to 2005)
# and set variables accordingly. Make sure they are zero padded

annual=true
fy=$(printf %02d $first_yr)

if [ $annual = true ]; then
    ly=$(printf %02d `expr $last_yr - 1`)      
else
    ly=$(printf %02d $last_yr)
fi

# Create the strata following those parameters. Here we will assume stable 
# unclassified (class 0) is stable pasture.

for yr in $(seq -w $fy $ly); do
    if [ $annual = true ]; then
        first=20$yr$suffix
        yr2=$(printf %02d `expr $yr + 1`)
        second=20$yr2$suffix
        outfile="final_strata_annual_"$yr"_"$yr2"_UTM18N.tif"
    else
        first=20$fy$suffix
        second=20$yr$suffix
        outfile="final_strata_"$fy"_"$yr"_UTM18N.tif"
    fi
    
    qsub -j y -b y -V -N strata_$yr \
     python $spath/6_stratification/create_strata.py $first $second \
     $spath/data/original_lut.csv strata_python_test.tif

done

