#!/bin/bash -l
# Basic version of mosaicking script for a single date. A more advanced version
# will be needed to specify the order of the rasters according the the scenes
# with the highest number of images
# MAGNA SIRGAS is EPSG:4686. Could also use MS-Bogota zone: 3116
# UTM18-19N are EPSG:32618-32619


module load gdal/1.11.1 
cd /projectnb/landsat/projects/Colombia/Mosaics/M1B

# Quick mosaic, stacks images in sequential order based on the folder pathrow
# Using 8 cores and more RAM to speed up the process

for yr in $(seq -w 03 15); do
    qsub -pe omp 4 -V -N mosaic_$yr -j y -b y gdalwarp --config \
    GDAL_CACHEMAX 500 -wm 500 -multi -co NBITS=4 -wt Byte -t_srs EPSG:4686 \
    -srcnodata 0 \
    "/projectnb/landsat/projects/Colombia/images/*/Results/M1/Class/ClassM1B_20"$yr".tif" \
    20$yr"_seq.tif"
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

