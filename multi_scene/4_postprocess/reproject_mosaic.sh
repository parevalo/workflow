#!/bin/bash -l

#module load python/2.7.5_nopath
#module load gdal/1.11.1

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal
# Last step
# Script to reproject the maps from UTM19 to UTM18, where most of the changes
# are concentrated 

# cd to the corresponding class folder
cd /projectnb/landsat/projects/Colombia/Mosaics/M3
    
# Project to UTM18N 
#qsub -V -b y -j y -N reproj_UTM18N \
#gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte -overwrite \
#-te 725295 -429765 1506435 595185 -te_srs EPSG:32618 \
#-t_srs EPSG:32618 -tr 30 -30 eastUTM19_2016.tif eastUTM18tapnewextent.tif

# (Re)project to UTM19N
qsub -V -b y -j y -N reproj_mUTM19N \
gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte -overwrite \
-te $(gdal_extent eastUTM19_2016.tif) -te_srs EPSG:32619 \
-t_srs EPSG:32619 -tr 30 -30 eastUTM18tapnewextent.tif east_reprojUTM19.tif
