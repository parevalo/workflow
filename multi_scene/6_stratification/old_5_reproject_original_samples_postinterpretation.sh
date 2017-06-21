#!/bin/bash -l

# Reprojects sample in east zone to UTM18 in order to create a single
# sample file after samples have been interpreted. The output file is 
# "snapped" to the  grid of the west zone in the same way it was done for 
# the mosaics so that they can be merged easily and match the strata maps.
# This is more of a legacy file, since this step was probably unnecesary
# to begin with, because the interpreted table could have been joined
# back to the original sample file.

# cd to the corresponding class folder
cd /projectnb/landsat/projects/Colombia/Mosaics/M3

# Rasterize east samples (original and additional)
qsub -j y -b y -V -N rast_samp1 \
    gdal_rasterize -a ID -a_nodata -9999 \
    -te $(gdal_extent eastUTM19_2016.tif) -tr 30 30 -ot Int16 \
    -co COMPRESS=PACKBITS sample_east_UTM19N_ID_PR_ZONE.shp \
    sample_east_UTM19N_ID_PR_ZONE.tif

qsub -j y -b y -V -N rast_samp2 \
    gdal_rasterize -a ID -a_nodata -9999 \
    -te $(gdal_extent eastUTM19_2016.tif) -tr 30 30 -ot Int16 \
    -co COMPRESS=PACKBITS additional_forest_sample_east_UTM19N_ID_PR_ZONE.shp \
    additional_forest_sample_east_UTM19N_ID_PR_ZONE.tif

# Reproject east samples to UTM18
qsub -j y -b y -V -N repr_sample1 -hold_jid rast_samp1 \
    gdalwarp -co COMPRESS=PACKBITS -wt Int16 -overwrite \
    -te 725295 -429765 1506435 595185 -te_srs EPSG:32618 \
    -t_srs EPSG:32618 -tr 30 -30 sample_east_UTM19N_ID_PR_ZONE.tif \
    sample_east_UTM18N_ID_PR_ZONE.tif

qsub -j y -b y -V -N repr_sample2 -hold_jid rast_samp2 \
    gdalwarp -co COMPRESS=PACKBITS -wt Int16 -overwrite \
    -te 725295 -429765 1506435 595185 -te_srs EPSG:32618 \
    -t_srs EPSG:32618 -tr 30 -30 additional_forest_sample_east_UTM19N_ID_PR_ZONE.tif \
    additional_forest_sample_east_UTM18N_ID_PR_ZONE.tif


# Poligonize those files
qsub -j y -b y -V -N pol_sample1 -hold_jid repr_sample1 \
    gdal_polygonize.py sample_east_UTM18N_ID_PR_ZONE.tif -f '"ESRI Shapefile"' \
    sample_east_UTM18N_ID_PR_ZONE.shp sample_east_UTM18N_ID_PR_ZONE ID

qsub -j y -b y -V -N pol_sample2 -hold_jid repr_sample2 \
    gdal_polygonize.py additional_forest_sample_east_UTM18N_ID_PR_ZONE.tif -f '"ESRI Shapefile"' \
    additional_forest_sample_east_UTM18N_ID_PR_ZONE.shp additional_forest_sample_east_UTM18N_ID_PR_ZONE ID


# Then manually join back the attribute table from the east zone for the two
# shapefiles and merge the shapefiles from the two zones...


