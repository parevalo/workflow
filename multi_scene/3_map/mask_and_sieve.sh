#!/bin/bash -l

# This script creates a mask for a given class, with a given distance in pixels
# around it. The area of the mask will be the valid area for gdal_sieve to
# operate upon. The masking process is split in two because gdal_proximity
# is not giving the expected output, so it needs to be reclassified.

module load python/2.7.5_nopath
module load gdal/1.11.1

scn_list="007059" #004058 004061 004062 005057 005058 005059 005060 \ 
          #005061 006058 006059 006060 006061 007058 007059 007060 \  
          #008058 008059 008060 009059" #004057

# General settings

rootdir=/projectnb/landsat/projects/Colombia/images
pre=ClassM2B_20
suf=-01-01_crop_class.tif

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    cd $rootdir/$s/Results/M2/Class
    
    for yr in $(seq -w 01 01); do

        # Calculate proximity in order to generate first mask

        qsub -j y -V -N mA_$pt$rw"_"$yr -b y \
         gdal_proximity.py $pre$yr$suf maskA_$yr"_2-3.tif" -co NBITS=2 \
          -ot Byte -values 2,3 -maxdist 4 -fixed-buf-val 1 -nodata 3 #-use_input_nodata YES

        # Reclassify proximity because it's not giving the output the way it's 
        # supposed to. Double quotes needed in the calc section to avoid 
        # problems when qsub'd

        qsub -j y -V -N mB_$pt$rw"_"$yr -hold_jid mA_$pt$rw"_"$yr -b y \
         gdal_calc.py -A maskA_$yr"_2-3.tif" --outfile=maskB_$yr"_2-3.tif" \
          --calc='"logical_or(A == 0, A == 1)"' --NoDataValue=3 --type=Byte \
           --co="NBITS=2"

        # Create water mask to remove from the sieving mask to make the sieving
        # more targeted
        
        qsub -j y -V -N mC_$pt$rw"_"$yr -hold_jid mB_$pt$rw"_"$yr -b y \
         gdal_calc.py -A $pre$yr$suf --outfile=maskC_$yr"_6.tif" \
          --calc='"(A != 6)"' --NoDataValue=3 --type=Byte --co="NBITS=2"

        # Get final mask
                
        qsub -j y -V -N mD_$pt$rw"_"$yr -hold_jid mC_$pt$rw"_"$yr -b y \
         gdal_calc.py -A maskB_$yr"_2-3.tif" -B maskC_$yr"_6.tif" \
           --outfile=maskD_$yr".tif" \
            --calc='"A*B"' --NoDataValue=3 --type=Byte --co="NBITS=2"

        # Run sieve with that mask

        export GDAL_CACHEMAX=2048
        qsub -j y -V -N sv_$pt$rw"_"$yr -hold_jid mD_$pt$rw"_"$yr -b y \
         gdal_sieve.py -st 12 -8  $pre$yr$suf -mask maskD_$yr".tif" \
          ClassM2B_$yr"_sievedF".tif

        # Post sieving-post processing...
        qsub -j y -V -N postsv_$yr -hold_jid sv_$pt$rw"_"$yr -b y \
         gdal_calc.py -A $pre$yr$suf -B ClassM2B_$yr"_sievedF.tif" \
           --outfile=fixed_sieveF.tif \
           --calc='"(logical_and(logical_and(A != 2, A != 3), A != 4))*A +'\
          '(A == 2)*B + (A == 3)*B + (A == 4)*B"'
    done
done

