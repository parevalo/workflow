#!/bin/bash -l

# Script to create area mosaics of change dates, using the
# new YATSM results for the more recent period available (e.g. not
# using the TS before 1997), without merging the segments after 
# running the Chow test, unlike the original results.
# The mosaics are created with the scenes in sequence, like the 
# original maps. VRTs cant be used for the entire script bc it is not posible
# to write VRTs from the python utilities.

# Steps:
# 1) mosaic individual zones with as vrt
# 2) reproject zone 19
# 3) scale to integer to reduce final size and processing time
# 4) mosaic whole area
# 5) clip to amazon boundary

cd /projectnb/landsat/projects/Colombia/Mosaics/M3_recent_period/sp_error_nocomm
imgf=/projectnb/landsat/projects/Colombia/images
prefix=firstchg_

for yr in $(seq -w 2003 2 2013); do
    yr2=$(expr $yr + 2)    

    ###### STEP 1: Zone mosaics
    echo "Creating zone mosaic for period $yr - $yr2"
    fld=/Results/M3/chg_maps/$prefix$yr-$yr2"_clipped.tif"
    fld2=/Results/M3_recent_period/chg_maps/$prefix$yr-$yr2"_clipped.tif"
    out_west=$prefix"westUTM18_"$yr-$yr2.vrt
    gdalbuildvrt $out_west \
    $imgf/006059/$fld2 \
    $imgf/006060/$fld2 $imgf/006061/$fld $imgf/007058/$fld \
    $imgf/007059/$fld $imgf/007060/$fld2 $imgf/007061/$fld2 \
    $imgf/008058/$fld $imgf/008059/$fld $imgf/008060/$fld \
    $imgf/009059/$fld $imgf/009060/$fld 

    out_east=$prefix"eastUTM19_"$yr-$yr2.vrt
    gdalbuildvrt $out_east \
    $imgf/003058/$fld $imgf/003059/$fld \
    $imgf/004057/$fld $imgf/004058/$fld2 $imgf/004059/$fld \
    $imgf/004061/$fld2 $imgf/004062/$fld $imgf/005057/$fld \
    $imgf/005058/$fld2 $imgf/005059/$fld $imgf/005060/$fld \
    $imgf/005061/$fld $imgf/006058/$fld 

    ###### STEP 2: Reproject east zone to UTM18
    echo "Reprojecting east zone for period $yr-$yr2"
    east18=$prefix"eastUTM18_"$yr-$yr2.vrt
    gdalwarp -te 725295 -429765 1506435 595185 -te_srs EPSG:32618 \
     -t_srs EPSG:32618 -tr 30 -30 -of VRT $out_east $east18

    ###### STEP 3: Create full mosaic.
    echo "Submit creation of full mosaic for period $yr-$yr2"
    out_mosaic=$prefix"full_"$yr-$yr2".tif"
    qsub -j y -V -N prefix_mosaic$yr-$yr2 -b y \
     gdal_merge.py -o $out_mosaic -v $east18 $out_west \
      -co "COMPRESS=PACKBITS" -co "BIGTIFF=YES" -a_nodata -9999

    ##### STEP 4: Unset NoData. Required to fix values of 0 showing as nodata. 
    # It is not required for analysis or sampling but makes visualization easier
    # Even though it runs very quickly, it needs to be submitted to use hold_jid
    echo "Submit unset NoData"
    qsub -j y -V -N unsetnd$yr-$yr2 -hold_jid prefix_mosaic$yr-$yr2 -b y \
     gdal_edit.py -unsetnodata $out_mosaic 

    ###### STEP 5: Clip to amazon boundary and set true nodata
    echo "Submit clipping to amazon boundary for period $yr-$yr2"
    mask=/projectnb/landsat/projects/Colombia/Mosaics/M3/true_data_mask_nodata_int.tif
    final_mosaic=$prefix"final_"$yr-$yr2".tif"
    qsub -j y -V -N clip$yr-$yr2 -hold_jid unsetnd$yr-$yr2 -b y \
     gdal_calc.py -A $out_mosaic -B $mask --outfile=$final_mosaic \
      --calc='"(A*0 + (A > 0)*A)*B"' --co="COMPRESS=PACKBITS" --co="BIGTIFF=YES" \
       --NoDataValue=-9999
done

