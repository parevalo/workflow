#!/bin/bash -l

# Script to create the strata from mosaics using python3/rasterio
# Overwrites by default

# CD to folder and define suffix of input tif files

cd /projectnb/landsat/projects/Colombia/Mosaics/M3
suffix="_final_crop.tif"
spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene
lut=for-nofor_lut_A.csv
# e.g _defor, or nothing for "original" stratification
lutsuf=_lutA

# Define the start and end year we want calculate the strata on, and the step
# (i.e. annual, biannual, etc)

first_yr=1
last_yr=16
step=1

# Define if we want the strata to be annual or cumulative (e.g. 2001 to 2005)
# and set variables accordingly. Make sure they are zero padded

annual=true
fy=$(printf %02d $first_yr)

if [ $annual = true ]; then
    ly=$(printf %02d `expr $last_yr - $step`)      
else
    ly=$(printf %02d $last_yr)
fi

# Create the strata following those parameters.
 
for yr in $(seq -w $fy $step $ly); do
    if [ $annual = true ]; then
        first=20$yr$suffix
        yr2=$(printf %02d `expr $yr + $step`)
        second=20$yr2$suffix
        outfile="final_strata_annual_"$yr"_"$yr2"_UTM18N"$lutsuf".tif"
    else
        first=20$fy$suffix
        second=20$yr$suffix
        outfile="final_strata_"$fy"_"$yr"_UTM18N"$lutsuf".tif"
    fi
    
    qsub -j y -b y -V -N strata_$yr -l mem_total=94G \
     python $spath/6_stratification/create_strata.py $first $second \
     $spath/data/$lut $outfile

done

# Create "original" strata between first and last year if we haven't done so
if [ $annual = true ]; then
    y1=$(printf %02d $first_yr)
    y2=$(printf %02d $last_yr)
    outfile="final_strata_"$y1"_"$y2"_UTM18N"$lutsuf".tif"

    qsub -j y -b y -V -N strata_$lutsuf -l mem_total=94G \
        python $spath/6_stratification/create_strata.py \
         20$y1$suffix 20$y2$suffix $spath/data/$lut $outfile

fi

