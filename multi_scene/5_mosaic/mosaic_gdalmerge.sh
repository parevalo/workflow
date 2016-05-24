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

# Quick mosaic, stacks images in sequential order based on the folder pathrow

for yr in $(seq -w 16 16); do
    qsub -V -N mosaicmrg_$yr -j y -b y \
     gdalmerge.py -o 20$yr$suf -co COMPRESS=PACKBITS -co NBITS=4 \
     -ot Byte  "$imgf"/*/Results/M3/Class/mergedmaps_20"$yr$suf" \
done



