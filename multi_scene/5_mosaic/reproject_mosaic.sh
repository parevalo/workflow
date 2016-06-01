#!/bin/bash -l

# Reprojects rasters in east zone to UTM18, where most of the changes are
# concentrated. The output file is "snapped" to the grid of the west zone
# mosaic to simplify merging the two.

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal

# cd to the corresponding class folder
cd /projectnb/landsat/projects/Colombia/Mosaics/M3
    
# Project to UTM18N. Extent calculated manually in order to "snap" the new
# raster to the west zone grid. 

for yr in $(seq -w 01 16); do
    qsub -V -b y -j y -N repr18N_$yr \
    gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte -overwrite \
    -te 725295 -429765 1506435 595185 -te_srs EPSG:32618 \
    -t_srs EPSG:32618 -tr 30 -30 eastUTM19_20$yr".tif" eastUTM18_20$yr".tif"
done

# (Re)project to UTM19N (in case it's needed later)
#qsub -V -b y -j y -N reproj_mUTM19N \
#gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte -overwrite \
#-te $(gdal_extent eastUTM19_2016.tif) -te_srs EPSG:32619 \
#-t_srs EPSG:32619 -tr 30 -30 eastUTM18tapnewextent.tif east_reprojUTM19.tif
