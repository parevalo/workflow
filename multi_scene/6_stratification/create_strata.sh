#!/bin/bash -l

# Script to create the strata from mosaics

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

initial=2001_final_crop.tif

# Get strata for each year we might want to check 

for yr in $(seq -w 16 16); do
    qsub -j y -V -N strata_$yr -b y \
     gdal_calc.py -A $initial -B 20$yr"_final_crop.tif" \
      --outfile=final_strata_01_$yr"_UTM18N.tif" \
      --calc='"logical_and(A == 1, B==1)*1 + logical_and(A == 2, B==2)*2 +' \
              'logical_and(A == 3, B==3)*3 + logical_and(A == 4, B==4)*4 +' \
              'logical_and(A == 5, B==5)*5 + logical_and(A == 5, B==1)*5 +' \
              'logical_and(A == 6, B==6)*6 +' \
              'logical_and(A == 7, B==7)*3 + logical_and(A == 1, B==4)*8 +' \
              'logical_and(A == 1, B==5)*9 +' \
              'logical_and(A == 1, logical_or(B==2, B==3))*0 +' \
              'logical_and(A == 1, logical_or(B==6, B==7))*0 +' \
              'logical_and(A == 1, B == 0)*8 +' \
              'logical_and(A == 4, B==5)*11 +' \
              'logical_and(logical_or(A == 2, A==3), B==5)*11 +' \
              'logical_and(logical_or(A == 6, A==7), B==5)*11 +' \
              'logical_and(logical_and(logical_and(A!=0,A!=1), A!=5), B==0)*13 +' \
              'logical_and(A == 5, logical_and(B != 5, B!= 1))*14 +' \
              'logical_and(A == 0, B==0)*15"' \
      --type=Byte --co="COMPRESS=PACKBITS" --overwrite
done
