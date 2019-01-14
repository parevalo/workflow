#!/bin/bash

# Script to create biannual stratification maps that incorporate a buffer 
# INTO the stable forest class  (after sieving it to avoid a lot of spurious buffers)
# This is aimed at tackling the forest loss omission errors. 
# Specifying the LUT allows to run this over
# the original classif. or over the forest-noforest maps. 

export GDAL_CACHEMAX=2048

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

# Set distance in pixels and other parameters
npix=3
spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene
lut=poststrat_buffer_lut_original.csv
# define lutsuf suffix, leave empty for "original" stratification
lutsuf=

# Set year and step info
first_yr=1
last_yr=16
step=2
fy=$(printf %02d $first_yr)
ly=$(printf %02d `expr $last_yr - $step`)

# Iterate over biannual files and run the following processes
for yr in $(seq -w $fy $step $ly); do
    
    # Set filenames
    yr2=$(printf %02d `expr $yr + $step`)
    bname="final_strata_annual_"$yr"_"$yr2"_UTM18N"$lutsuf
    infile=$bname".tif"
    out_sieve=$bname"_sieved4buff.tif"
    out_buff=$bname"_buffer"$npix    
    out_strata="buffered"$npix"_"$bname".tif"

    # Sieve strata raster
    qsub -j y -b y -V -N  sieve_$yr gdal_sieve.py -st 6 -8 $infile $out_sieve

    # Buffer from classes of interest (e.g. forest to pastures) INTO the forest
    qsub -j y -b y -V -N buff_$yr -hold_jid sieve_$yr \
     gdal_proximity.py $out_sieve $out_buff".tif" \
      -values 8 -maxdist $npix \
       -fixed-buf-val 1 -nodata 0 -ot Byte -co "COMPRESS=PACKBITS"

    # Create buffered version of the stratification
    qsub -j y -b y -V -N poststrat -l mem_total=94G -hold_jid buff_$yr \
     $spath/6_stratification/helper_scripts/create_strata.py \
      $bname".tif" $out_buff".tif" $spath/data/$lut $out_strata
done

