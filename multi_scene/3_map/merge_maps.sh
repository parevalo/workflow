#!/bin/bash -l

# This script was originally intended to  merge the output of models M2B and M3.
# This was required because M2B is better at estimating stable forest but less 
# than ideal at estimating deforestation in very dynamic areas and viceversa. 
# However, model M2B won't be used anymore, and instead, two versions of M3, 
# one classified with classifier M1 and the other with M3.

module load python/2.7.5_nopath
module load gdal/1.11.1

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings

rootdir=/projectnb/landsat/projects/Colombia/images
pre=ClassM3_20
suf="-01-01.tif"

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    cd $rootdir/$s/Results/M3/Class
    export GDAL_CACHEMAX=2048

    for yr in $(seq -w 16 16); do

	    # Replace grassl. in M1 with whatever is in M3 as long as there is 
        # no change in the ENTIRE  period

        qsub -V -j y -b y -N mapmerge_$pt$rw \
        gdal_calc.py -A $pre$yr"-01-01_M1train.tif" -B $pre$yr"-01-01_M3train.tif" \
         -C numchange_2001-2016_$pt$rw".tif" --outfile=mergedmaps_20$yr$suf \
         --calc='"(C == -9999)*((A == 2)*B + (A != 2)*A) + ((C != -9999)*A)"' \
         --type=Byte --co NBITS=4 --overwrite --NoDataValue=0

    done
done

