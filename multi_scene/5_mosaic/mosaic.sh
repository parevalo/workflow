#!/bin/bash -l

# Script to create mosaics for each of the UTM zones. The mosaics are created
# with the scenes in sequence because this improves the appearance, and because
# most of the training data and areas of interest are located towards the west

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal

cd /projectnb/landsat/projects/Colombia/Mosaics/M3
imgf=/projectnb/landsat/projects/Colombia/images
suf=_final.tif
fld=Results/M3/Class/mergedmaps_20

# Mosaics for WEST ZONE (UTM18N) 

for yr in $(seq -w 01 16); do
    qsub -V -N mrgUTM18_$yr -j y -b y \
     gdal_merge.py -o westUTM18_20$yr".tif" -co COMPRESS=PACKBITS -co NBITS=4 \
     -ot Byte -n 0 \
    $imgf/006059/$fld$yr$suf \
    $imgf/006060/$fld$yr$suf $imgf/006061/$fld$yr$suf $imgf/007058/$fld$yr$suf \
    $imgf/007059/$fld$yr$suf $imgf/007060/$fld$yr$suf $imgf/007061/$fld$yr$suf \
    $imgf/008058/$fld$yr$suf $imgf/008059/$fld$yr$suf $imgf/008060/$fld$yr$suf \
    $imgf/009059/$fld$yr$suf $imgf/009060/$fld$yr$suf 
done

# Mosaics for EAST ZONE (UTM19N)
for yr in $(seq -w 01 16); do
    qsub -V -N mrgUTM19_$yr -j y -b y \
     gdal_merge.py -o eastUTM19_20$yr".tif" -co COMPRESS=PACKBITS -co NBITS=4 \
     -ot Byte -n 0 \
    $imgf/003058/$fld$yr$suf $imgf/003059/$fld$yr$suf \
    $imgf/004057/$fld$yr$suf $imgf/004058/$fld$yr$suf $imgf/004059/$fld$yr$suf \
    $imgf/004061/$fld$yr$suf $imgf/004062/$fld$yr$suf $imgf/005057/$fld$yr$suf \
    $imgf/005058/$fld$yr$suf $imgf/005059/$fld$yr$suf $imgf/005060/$fld$yr$suf \
    $imgf/005061/$fld$yr$suf $imgf/006058/$fld$yr$suf 
done

