#!/bin/bash 

# This script rasterizes the training data. It gets the extent of the 
# first L7 image it finds and uses those values for the raster output.
# Its also hardcoded to rasterize all the training shapefiles with the
# same name

# Scenes for which there is training data
scn_list="004057" #005058 005060 006058 006059 006060 007058 007059 008059

# Folder setting
rootdir=/projectnb/landsat/projects/Colombia/images
train_dir=/projectnb/landsat/projects/Colombia/Training

# Iterate over scenes
for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # Set scene TS folder and find example image
    ts_path=$rootdir/$s/images
    img=$(find $ts_path -maxdepth 1 -type d -name "*LE7*" | head -1)
    example_img=$(basename $img)
    img_path=$ts_path/$example_img/$example_img"_stack"

    # Get extent, ugly awk but does the job for now
    LL_X=$(gdalinfo $img_path | grep "Lower Left" | awk -F " " '{ print $4 }' | awk -F "," '{print $1}')
    LL_Y=$(gdalinfo $img_path | grep "Lower Left" | awk -F " " '{ print $5 }' | awk -F ")" '{print $1}')
    UR_X=$(gdalinfo $img_path | grep "Upper Right" | awk -F " " '{ print $4 }'| awk -F "," '{print $1}')
    UR_Y=$(gdalinfo $img_path | grep "Upper Right" | awk -F " " '{ print $5 }'| awk -F ")" '{print $1}')

    cd /projectnb/landsat/projects/Colombia/logs/$pt$rw
    
    # Submit rasterizing job
    qsub -j y -b y -V -N rstrze_$pt$rw \
     gdal_rasterize -a Class -init 0 -a_nodata 0 -te $LL_X $LL_Y $UR_X $UR_Y \
      -tr 30 30 -ot UInt16 $train_dir/$s/Training1.shp $train_dir/$s/Training1.tif

done
