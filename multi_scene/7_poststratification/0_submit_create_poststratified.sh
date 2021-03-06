#!/bin/bash

# Script to create a buffer INTO the stable forest class of the
# 2001-2016 stratification (after sieving it to avoid a lot of spurious buffers)
# in order to do a poststratification, mostly aimed at
# tackling the omission errors. Specifying the LUT allows to run this over
# the original classif. or over the forest-noforest maps. 

export GDAL_CACHEMAX=2048

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

# Set distance in pixels and file names
npix=3
bname=final_strata_01_16_UTM18N_lutA
in_name=$bname".tif"
out_sieve=$bname"_sieved_for_buffer.tif"
out_buff=$bname"_buffer"$npix"B"
out_strata="buffered"$npix"B_"$bname".tif"

spath=/projectnb/landsat/projects/Colombia/workflow/multi_scene
#lut=poststrat_buffer_lut_original.csv
lut=poststrat_buffer_lutA.csv

# Sieve strata raster
qsub -j y -b y -V -N  sieve_strata gdal_sieve.py -st 6 -8 $in_name $out_sieve

# Buffer from classes of interest (e.g. pastures)  INTO the forest.
# Original classes used were -values 4,5,8,9,11,13,14 \

qsub -j y -b y -V -N buff_gdal -hold_jid sieve_strata \
 gdal_proximity.py $out_sieve $out_buff".tif" \
  -values 3 -maxdist $npix \
   -fixed-buf-val 1 -nodata 0 -ot Byte -co "COMPRESS=PACKBITS"

# Create strata with the new buffer

qsub -j y -b y -V -N poststrat -l mem_total=94G -hold_jid buff_gdal \
 $spath/6_stratification/create_strata.py \
  $bname".tif" $out_buff".tif" $spath/data/$lut $out_strata

