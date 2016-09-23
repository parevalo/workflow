#!/bin/bash -l

# Script to process the maps to compare them to those of IDEAM
# Has three main steps
# 1) Create custom strata with only forest/non forest classes
# (All "other" changes will still be labeled as zero)
# 2) Reproject to MAGNA-SIRGAS and match grids.
# 3) Sieve polygons of less than 1 ha (or around 11.11 pixels)
# 4) Use gdalcalc to compare the maps

cd /projectnb/landsat/projects/Colombia/Mosaics/M3/IDEAM
export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal
export GDAL_CACHEMAX=2048

# Create a list of years we want to get strata from, omitting 00-02
# Periods have +1 year to make them more comparable to my maps (which correspond
# to the first day of the year
period_list="03-05 05-07 07-09 09-10 11-13 13-14 14-15"

# Calculate strata for each year we might want to check

# 1 stable forest, 2 forest to others (forest loss), 3 forest to regrowth, 
# 4 regrowth to forest (regrowth), 5 Stable nonforest, 
# 6 Stable regrowth, 7 Anything to regrowth (regrowth), 8 Regrowth to others 
# 9 others to forest (regrowth)
for p in $period_list; do
    f="custom_strata__"$p"_UTM18N.tif"
    
    # 1) Create custom strata 
    y1=${p:0:2}
    y2=${p:3:2}
#    qsub -j y -V -N strata_$p -b y \
#     gdal_calc.py -A ../20$y1"_final_crop.tif" -B ../20$y2"_final_crop.tif" \
#      --outfile=$f \
#      --calc='"logical_and(A == 1, B==1)*1 + logical_and(A == 2, B==2)*5+' \
#              'logical_and(A == 3, B==3)*5 + logical_and(A == 4, B==4)*5 +' \
#              'logical_and(A == 5, B==5)*6 + logical_and(A == 5, B==1)*4+' \
#              'logical_and(A == 6, B==6)*5 +' \
#              'logical_and(A == 7, B==7)*5 + logical_and(A == 1, B==4)*2 +' \
#              'logical_and(A == 1, B==5)*3 +' \
#              'logical_and(A == 1, logical_or(B==2, B==3))*2 +' \
#              'logical_and(A == 1, logical_or(B==6, B==7))*2 +' \
#              'logical_and(A == 1, B == 0)*2 +' \
#              'logical_and(A == 4, B==5)*7 +' \
#              'logical_and(logical_or(A == 2, A==3), B==5)*7 +' \
#              'logical_and(logical_or(A == 6, A==7), B==5)*7 +' \
#              'logical_and(logical_and(logical_and(A!=0,A!=1), A!=5), B==0)*5 +' \
#              'logical_and(A == 5, logical_and(B != 5, B!= 1))*8 +' \
#              'logical_and(logical_or(A == 2, A==3), B==1)*9 +' \
#              'logical_and(logical_or(A == 4, A==6), B==1)*9 +' \
#              'logical_and(A == 7, B==1)*9 +' \
#              'logical_and(A == 0, B==0)*15"' \
#      --type=Byte --co="COMPRESS=PACKBITS" --overwrite
#
    # 2) Reproject to MAGNA-SIRGAS and grid
    # Use 15 as input and output nodata, because the "other to other" class is 
    # still labeled as zero. User -ts instead of -tr because that was creating
    # a displacement in the grid for some reason. IDEAM annual maps have a 
    # slightly different extent and number of rows

    if [ $p == 13-14 ] || [ $p == 14-15 ] ;then                                 
        height=61620
        refextent=cambio_2013_2014_v6_amazonia.tif                                                       
    else                                                                        
        height=61618    
        refextent=cambio_2002_2004_v6_amazonia.tif                                                  
    fi         

    qsub -V -b y -j y -N reproj_$p -hold_jid strata_$p \
     gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte \
      -te $(gdal_extent $refextent) -te_srs EPSG:3116 \
       -t_srs EPSG:3116 -ts 45206 $heigth -srcnodata 15 -dstnodata 15 \
        -overwrite $f $(basename $f .tif)"_MAGNA.tif"

    # 3) Sieve "polygons" of less than 1 ha (around 11.11 pixels)
    qsub -V -b y -j y -N sieve_$p -hold_jid reproj_$p  gdal_sieve.py -st 11 -8 \
     $(basename $f .tif)"_MAGNA.tif" $(basename $f .tif)"_MAGNA_sieved.tif"


    # 4) Use gdalcalc to compare the the maps
    # Revert to IDEAM's numbering
    y1p=`expr $y1 - 1`
    y2p=`expr $y2 - 1`
    printf -v yr1 '%02d' $y1p 
    printf -v yr2 '%02d' $y2p

    #IDEAM labels: 1 forest, 2 defor, 3 no info, 4 regrowth, 5 stable non-forest        
    qsub -j y -V -N compare_$p -hold_jid sieve_$p -b y \
     gdal_calc.py -A "cambio_20"$yr1"_20"$yr2"_v6_amazonia.tif" \
      -B $(basename $f .tif)"_MAGNA_sieved.tif" \
      --outfile=comparison_$p".tif" \
       --calc='"logical_and(A == 1, B==1)*1 + logical_and(A == 2, B==2)*2+' \
               'logical_and(A == 4, logical_or(B==3, B==4 ))*4 +' \
               'logical_and(A == 4, logical_or(B==6, B==7 ))*4 +' \
               'logical_and(A == 4, B==9 )*4 +' \
               'logical_and(A == 5, B==5 )*5 +' \
               'logical_and(A == 2, B==1 )*6 +' \
               'logical_and(A == 2, logical_or(B==3, B==8 ))*7 +' \
               'logical_and(A == 5, B==1 )*8 +' \
               'logical_and(A == 5, B==6 )*9 +' \
               'logical_or(A == 0 , B==15 )*14"' \
       --type=Byte --NoDataValue=15 --co="COMPRESS=PACKBITS" --co NBITS=4 --overwrite

    #1 forest agreem, 2 defor agreem, 4 regrowth agreem, 5 nonforest agreem, 
    #6 IDEAM defor - strata forest, 7 IDEAM-defor - strata- 3 and 8 ("loss" of forest)
    #8 IDEAM NF - strata forest, 9 IDEAM NF -strata stable regrowth.
    #15 total nodata areas from both layers
done

