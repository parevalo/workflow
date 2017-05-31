#!/bin/bash
#$ -b y
#$ -j y
#$ -N create_sample
#$ -V

# Submit sample creation. 

fpath=/projectnb/landsat/projects/Colombia/Mosaics/M3

#python sample_map.py -v --size 850 \
#--allocation "50 200 75 50 75 50 50 50 50 50 50 50 50" \
#--mask "15 255" --ndv 255 --vector $fpath/sample_may2016.shp --seed_val 10000 \
#stratified $fpath/strata_01_16_UTM18N.tif

# Additional forest sample

python sample_map.py -v --size 200 \
--allocation "0 200 0 0 0 0 0 0 0 0 0 0 0" \
--mask "15 255" --ndv 255 --vector $fpath/additional_forest_sample_UTM18N.shp \
--seed_val 10001 stratified $fpath/strata_01_16_UTM18N.tif

