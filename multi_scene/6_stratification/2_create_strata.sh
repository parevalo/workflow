#!/bin/bash -l

# Script to create the strata from mosaics

# CD to folder and define suffix of input tif files

cd /projectnb/landsat/projects/Colombia/Mosaics/M3
suffix="_final_crop.tif"


# Define the start and end year we want calculate the strata on

first_yr=15
last_yr=16

# Define if we want the strata to be annual or cumulative (e.g. 2001 to 2005)
# and set variables accordingly. Make sure they are zero padded

annual=true
fy=$(printf %02d $first_yr)

if [ $annual = true ]; then
    ly=$(printf %02d `expr $last_yr - 1`)      
else
    ly=$(printf %02d $last_yr)
fi

# Create the strata following those parameters. Here we will assume stable 
# unclassified (class 0) is stable pasture.

for yr in $(seq -w $fy $ly); do
    if [ $annual = true ]; then
        first=20$yr$suffix
        yr2=$(printf %02d `expr $yr + 1`)
        second=20$yr2$suffix
        outfile="final_strata_annual_"$yr"_"$yr2"_UTM18N.tif"
    else
        first=20$fy$suffix
        second=20$yr$suffix
        outfile="final_strata_"$fy"_"$yr"_UTM18N.tif"
    fi
    
    qsub -j y -V -N strata_$yr -b y \
     gdal_calc.py -A $first -B $second --outfile=$outfile \
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
              'logical_and(A == 15, B==15)*15 +' \
              'logical_and(A == 0, B==0)*4"' \
      --type=Byte --co="COMPRESS=PACKBITS" --overwrite
done

