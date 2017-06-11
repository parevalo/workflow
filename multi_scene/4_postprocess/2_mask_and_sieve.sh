#!/bin/bash -l

# This script run the mask over areas of grasslands, croplands and urban
# and avoid any modification over other land cover classes, or over pixels
# that experienced change during the entire time period in any of the three 
# classes of interest

# This script goes BEFORE the crop script and the input maps MUST NOT
# have set a NoData value


scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings

rootdir=/projectnb/landsat/projects/Colombia/images
pre=mergedmaps_20
suf=-01-01.tif

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    cd $rootdir/$s/Results/M3/Class
    
    # We're only doing this for the last year, then backpropagating changes 
    for yr in $(seq -w 16 16); do

        # Create mask to remove areas we don't want to be included in the 
        # sieving, as described above. 
       
#        qsub -j y -V -N m_$pt$rw"_"$yr -b y \
         gdal_calc.py -A $pre$yr$suf -B numchange_2001-2016_$pt$rw".tif" \
          --outfile=sieve_mask_$yr".tif" \
          --calc="logical_and(logical_or(logical_or(A == 2, A == 4), A == 3), (B == -9999))" \
          --overwrite 

        # Run sieve with that mask

        export GDAL_CACHEMAX=2048
#        qsub -j y -V -N sv_$pt$rw"_"$yr -hold_jid m_$pt$rw"_"$yr -b y \
         gdal_sieve.py -st 10 -8  $pre$yr$suf -mask sieve_mask_$yr".tif" \
          $pre$yr"_sieved".tif

    done
done

