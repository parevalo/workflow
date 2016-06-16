#!/bin/bash -l

# Reprojects sample in east zone to UTM18. The output file is "snapped" to the 
# grid of the west zone in the same way it was done for the mosaics
#  so that they can be merged easily

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal

# cd to the corresponding class folder
cd /projectnb/landsat/projects/Colombia/Mosaics/M3

# Rasterize east sample, need to define column to raterize
qsub -j y -b y -V -N rast_samp1 \
    gdal_rasterize -a strata -init 0 -a_nodata 255 \
    -te $(gdal_extent eastUTM19_2016.tif) -tr 30 30 -ot Byte \
    -co COMPRESS=PACKBITS sample_east_UTM19N_ID_PR_ZONE.shp \
    sample_east_UTM19N_ID_PR_ZONE.tif

# Rasterize west sample, need to define column to rasterize
qsub -j y -b y -V -N rast_samp2 \
    gdal_rasterize -a strata -init 0 -a_nodata 255 \
    -te $(gdal_extent westUTM18_2016.tif) -tr 30 30 -ot Byte \
    -co COMPRESS=PACKBITS sample_west_UTM18N_ID_PR_ZONE.shp \
    sample_west_UTM18N_ID_PR_ZONE.tif

# Reproject east sample to UTM18
qsub -j y -b y -V -N repr_sample -hold_jid rast_samp2 \
    gdalwarp -co COMPRESS=PACKBITS -wt Byte -overwrite \
    -te 725295 -429765 1506435 595185 -te_srs EPSG:32618 \
    -t_srs EPSG:32618 -tr 30 -30 sample_east_UTM19N_ID_PR_ZONE.tif \
    sample_east_UTM18N_ID_PR_ZONE.tif

# Re-merge samples

qsub -j y -b y -V -N merge_sample -hold_jid repr_sample \
    gdal_merge.py -o interpreted_sample_UTM18N.tif -co COMPRESS=PACKBITS \
    -ot Byte -n o -v sample_west_UTM18N_ID_PR_ZONE.tif \
    sample_east_UTM18N_ID_PR_ZONE.tif


