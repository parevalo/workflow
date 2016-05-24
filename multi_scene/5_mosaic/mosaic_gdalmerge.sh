#!/bin/bash -l
# Mosaic script with sequential or ordered scenes. (sequential seems to work
# better) 
# MAGNA SIRGAS is EPSG:4686. Could also use MS-Bogota zone: 3116
# UTM18-19N are EPSG:32618-32619

module load python/2.7.5_nopath
module load gdal/1.11.1

cd /projectnb/landsat/projects/Colombia/Mosaics/M3
imgf=/projectnb/landsat/projects/Colombia/images
suf=_final_UTM18N.tif
fld=Results/M3/Class/mergedmaps_20

# Quick mosaic, stacks images in sequential order based on the folder pathrow

for yr in $(seq -w 16 16); do
    qsub -V -N mosaicmrg_$yr -j y -b y \
     gdal_merge.py -o merge_20$yr$suf -co COMPRESS=PACKBITS -co NBITS=4 \
     -ot Byte $imgf/003058/$fld$yr$suf $imgf/003059/$fld$yr$suf \
    $imgf/004057/$fld$yr$suf $imgf/004058/$fld$yr$suf $imgf/004059/$fld$yr$suf \
    $imgf/004061/$fld$yr$suf $imgf/004062/$fld$yr$suf $imgf/005057/$fld$yr$suf \
    $imgf/005058/$fld$yr$suf $imgf/005059/$fld$yr$suf $imgf/005060/$fld$yr$suf \
    $imgf/005061/$fld$yr$suf $imgf/006058/$fld$yr$suf $imgf/006059/$fld$yr$suf \
    $imgf/006060/$fld$yr$suf $imgf/006061/$fld$yr$suf $imgf/007058/$fld$yr$suf \
    $imgf/007059/$fld$yr$suf $imgf/007060/$fld$yr$suf $imgf/007061/$fld$yr$suf \
    $imgf/008058/$fld$yr$suf $imgf/008059/$fld$yr$suf $imgf/008060/$fld$yr$suf \
    $imgf/009059/$fld$yr$suf $imgf/009060/$fld$yr$suf 
done

