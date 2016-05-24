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
# Using 8 cores and more RAM to speed up the process.

for yr in $(seq -w 16 16); do
    qsub -pe omp 4 -V -N mosaic_$yr -j y -b y gdalwarp --config \
     GDAL_CACHEMAX 4000 -wm 4000 -multi -co COMPRESS=PACKBITS -co NBITS=4 \
     -wt Byte -wo NUM_THREADS=4 -tr 30 30 -srcnodata 0 \
     "$imgf"/*/Results/M3/Class/mergedmaps_20"$yr$suf" \
     20$yr$suf
done


# Mosaic with scenes ordered with the highest number of images at the end
# of the list (in order to be on top). Didn't provide good results but left
# in case it's needed.

#sp1=/projectnb/landsat/projects/Colombia/images
#sp2=Results/M1/Class/ClassM1B_20

#for yr in $(seq -w 01 01); do
#    qsub -pe omp 4 -V -N mosaic_$yr -j y -b y gdalwarp --config \
#    GDAL_CACHEMAX 500 -wm 500 -multi -co NBITS=4 -wt Byte -t_srs EPSG:4686 \
#    -srcnodata 0 \
#    $sp1/005061/$sp2$yr".tif" $sp1/006061/$sp2$yr".tif" $sp1/005060/$sp2$yr".tif" \
#    $sp1/008060/$sp2$yr".tif" $sp1/007059/$sp2$yr".tif" $sp1/008059/$sp2$yr".tif" \
#    $sp1/007058/$sp2$yr".tif" $sp1/004062/$sp2$yr".tif" $sp1/007060/$sp2$yr".tif" \
#    $sp1/009059/$sp2$yr".tif" $sp1/006058/$sp2$yr".tif" $sp1/008058/$sp2$yr".tif" \
#    $sp1/006060/$sp2$yr".tif" $sp1/004061/$sp2$yr".tif" $sp1/005057/$sp2$yr".tif" \
#    $sp1/005059/$sp2$yr".tif" $sp1/006059/$sp2$yr".tif" $sp1/004057/$sp2$yr".tif" \
#    $sp1/005058/$sp2$yr".tif" $sp1/004058/$sp2$yr".tif" \
#    20$yr"_new_MAGNA.tif"
#done

