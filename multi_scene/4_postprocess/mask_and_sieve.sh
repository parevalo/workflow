#!/bin/bash -l

# This script creates a mask for a given class, with a given distance in pixels
# around it. The area of the mask will be the valid area for gdal_sieve to
# operate upon. The masking process is split in two because gdal_proximity
# is not giving the expected output, so it needs to be reclassified.
# This script goes BEFORE the crop script.

module load python/2.7.5_nopath
module load gdal/1.11.1

scn_list="007059" #003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          #005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          #007060 007061 008058 008059 008060 009059 009060"

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
    
    for yr in $(seq -w 01 01); do

        # Calculate proximity in order to generate first mask

        qsub -j y -V -N mA_$pt$rw"_"$yr -b y \
         gdal_proximity.py $pre$yr$suf maskA_$yr"_2-3.tif" -co NBITS=2 \
          -ot Byte -values 2,3 -maxdist 2 -fixed-buf-val 1 -nodata 3 #-use_input_nodata YES

        # Reclassify proximity because it's not giving the output the way it's 
        # supposed to. Double quotes needed in the calc section to avoid 
        # problems when qsub'd

        qsub -j y -V -N mB_$pt$rw"_"$yr -hold_jid mA_$pt$rw"_"$yr -b y \
         gdal_calc.py -A maskA_$yr"_2-3.tif" --outfile=maskB_$yr"_2-3.tif" \
          --calc='"logical_or(A == 0, A == 1)"' --NoDataValue=3 --type=Byte \
           --co="NBITS=2" --overwrite

        # Create mask to remove areas we don't want to have in the sieving mask,
		# so that they cannot be modified e.g. all water areas and areas
        # were yatsm detected changes during the entire study period. CHECK
        
        qsub -j y -V -N mC_$pt$rw"_"$yr -hold_jid mB_$pt$rw"_"$yr -b y \
         gdal_calc.py -A $pre$yr$suf -B numchange_2001-2015_$pt$rw".tif" \
		  --outfile=maskC_$yr"_6.tif" --calc='"(A != 6) * (B == -9999)"' \
          --NoDataValue=3 --type=Byte --co="NBITS=2" --overwrite

        # Get final mask
                
        qsub -j y -V -N mD_$pt$rw"_"$yr -hold_jid mC_$pt$rw"_"$yr -b y \
         gdal_calc.py -A maskB_$yr"_2-3.tif" -B maskC_$yr"_6.tif" \
           --outfile=maskD_$yr".tif" --calc='"A*B"' --NoDataValue=3 \
           --type=Byte --co="NBITS=2" --overwrite

        # Run sieve with that mask

        export GDAL_CACHEMAX=2048
        qsub -j y -V -N sv_$pt$rw"_"$yr -hold_jid mD_$pt$rw"_"$yr -b y \
         gdal_sieve.py -st 10 -8  $pre$yr$suf -mask maskD_$yr".tif" \
          $pre$yr"_sieved".tif

        # Post sieving-post processing, in order to revert any changes
		# made to other land cover types (including pastures)

        qsub -j y -V -N postsv_$yr -hold_jid sv_$pt$rw"_"$yr -b y \
         gdal_calc.py -A $pre$yr$suf -B $pre$yr"_sieved.tif" \
           --outfile=$pre$yr"_postsieved.tif" \
           --calc='"logical_and(A != 2, A != 3)*A +'\
          '(A == 2)*B + (A == 3)*B"' --overwrite
    done
done

