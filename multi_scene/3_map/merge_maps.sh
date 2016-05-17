#!/bin/bash -l

# This script merges the output of models M2B and M3. This is required
# because M2B is better at estimating stable forest but less than ideal
# at estimating deforestation in very dynamic areas and viceversa. 
# TODO Test if this runs fast directly in the  mosaics 

module load python/2.7.5_nopath
module load gdal/1.11.1

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings

rootdir=/projectnb/landsat/projects/Colombia/images
pre=ClassM3_20
suf=-01-01.tif

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    cd $rootdir/$s/Results/M3/Class
    
    for yr in $(seq -w 01 01); do
 		# Create variable for model M2
		m2_path=$rootdir/$s/Results/M2/Class/ClassM2B_20$yr$suf

		# Merge the maps of the two models (A = M3, B = M2) TEST!

        qsub -j y -V -N merge_$pt$rw"_"$yr -b y \
         gdal_calc.py -A $pre$yr$suf -B $m2_path \
           --outfile=Finalclass_$yr$suf \
           --calc='"logical_and(A == 4, B == 1)*B +'\
          'logical_and(A != 4, B != 1)*A"'
    done
done

