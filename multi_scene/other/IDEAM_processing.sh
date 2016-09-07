#!/bin/bash -l

# Script to process the maps to compare them to those of IDEAM
# Has three main steps
# 1) Create custom strata with only forest/non forest classes
# (All "other" changes will still be labeled as zero)
# 2) Reproject to MAGNA-SIRGAS and match grids.
# 3) Sieve polygons of less than 1 ha (or around 11.11 pixels)

cd /projectnb/landsat/projects/Colombia/Mosaics/M3/IDEAM
export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal
export GDAL_CACHEMAX=2048

# Create a list of years we want to get strata from, omitting 00-02
# Periods have +1 year to make them more comparable to my maps (which correspond
# to the first day of the year
period_list="05-07 07-09 09-10 11-13 13-14 14-15" #03-05

# Calculate strata for each year we might want to check
# We can classify regrowing vegetation as forest or not
# How to label change from 1 to 2, 3, 6, 7?
# How to label all to unclassified? (13)

# Without regrowing vegetation as forest
for p in $period_list; do
    f="custom__strata_"$p"_UTM18N.tif"
    
    # 1) Create custom strata
    y1=${p:0:2}
    y2=${p:3:2}
    qsub -j y -V -N strata_$p -b y \
     gdal_calc.py -A ../20$y1"_final_crop.tif" -B ../20$y2"_final_crop.tif" \
      --outfile=$f \
      --calc='"logical_and(A == 1, B==1)*1 + logical_and(A == 2, B==2)*2+' \
              'logical_and(A == 3, B==3)*2 + logical_and(A == 4, B==4)*2 +' \
              'logical_and(A == 5, B==5)*2 + logical_and(A == 5, B==1)*1+' \
              'logical_and(A == 6, B==6)*2 +' \
              'logical_and(A == 7, B==7)*2 + logical_and(A == 1, B==4)*3 +' \
              'logical_and(A == 1, B==5)*3 +' \
              'logical_and(A == 1, logical_or(B==2, B==3))*3 +' \
              'logical_and(A == 1, logical_or(B==6, B==7))*3 +' \
              'logical_and(A == 1, B == 0)*3 +' \
              'logical_and(A == 4, B==5)*2 +' \
              'logical_and(logical_or(A == 2, A==3), B==5)*2 +' \
              'logical_and(logical_or(A == 6, A==7), B==5)*2 +' \
              'logical_and(logical_and(logical_and(A!=0,A!=1), A!=5), B==0)*2 +' \
              'logical_and(A == 5, logical_and(B != 5, B!= 1))*2 +' \
              'logical_and(A == 0, B==0)*15"' \
      --type=Byte --co="COMPRESS=PACKBITS" --overwrite

    # 2) Reproject to MAGNA-SIRGAS and grid
    # This section works fine except because the output file has a 2 meter 
    # offset in the lower left Y coordinate for an unknown reason. Also, the 
    # annual IDEAM files have a sligtly different origin than the others


    # Use 15 as input and output nodata, because the "other to other" class is 
    # still labeled as zero
    qsub -V -b y -j y -N reproj_IDEAM -hold_jid strata_$p \
     gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte \
      -te $(gdal_extent cambio_2002_2004_v6_amazonia.tif) -te_srs EPSG:3116 \
       -t_srs EPSG:3116 -tr 30.717 -30.2624 -srcnodata 15 -dstnodata 15 \
        -overwrite $f $(basename $f .tif)"_MAGNA.tif"
    

    # 3) Sieve "polygons" of less than 1 ha (around 11.11 pixels)
    qsub -V -b y -j y -N sieve_$p -hold_jid strata_$p  gdal_sieve.py -st 11 -8 \
     $(basename $f .tif)"_MAGNA.tif" $(basename $f .tif)"_MAGNA_sieved.tif"


done

