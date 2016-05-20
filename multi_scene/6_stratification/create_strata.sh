#!/bin/bash -l

# Script to create the strata from mosaics

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

initial=2001_final_crop.tif

# Get strata for each of the three end years we want to check 

for yr in $(seq -w 15 16); do
    qsub -j y -V -N strata_$yr -b y \
     gdal_calc.py -A $initial -B 20$yr"_final_crop.tif" \
      --outfile=strata_01_$yr".tif" \
      --calc='"logical_and(A == 1, B==1)*1 + logical_and(A == 2, B==2)*2 +' \
              'logical_and(A == 3, B==3)*3 + logical_and(A == 4, B==4)*4 +' \
              'logical_and(A == 5, B==5)*5 + logical_and(A == 6, B==6)*6 +' \
              'logical_and(A == 7, B==7)*7 + logical_and(A == 1, B==4)*8 +' \
              'logical_and(A == 1, B==5)*9 +' \
              'logical_and(A == 1, logical_or(B==2, B==3))*10 +' \
              'logical_and(A == 1, logical_or(B==6, B==7))*10 +' \
              'logical_and(A == 4, B==5)*11 + logical_and(A != 1, B==1)*12 +' \
              'logical_and(A != 0, B==0)*13 + logical_and(A == 0, B==0)*14"' \
      --type=Byte --co="COMPRESS=PACKBITS" --overwrite
done
