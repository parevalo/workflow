#!/bin/bash -l

# Generic script with parameters to properly project east zone to
# UTM19. Useful to reproject samples and strata for interpretation.

# cd to the corresponding class folder
cd /projectnb/landsat/projects/Colombia/Mosaics/M3
    
# (Re)project to UTM19N (needed to reproject strata and samples)
qsub -V -b y -j y -N reproj_UTM19N \
gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte \
-te $(gdal_extent eastUTM19_2016.tif) -te_srs EPSG:32619 \
-t_srs EPSG:32619 -tr 30 -30 additional_forest_sample_UTM18N.tif additional_forest_sample_UTM19N.tif
#-t_srs EPSG:32619 -tr 30 -30 sample_may2016_UTM18N.tif sample_may2016_UTM19N.tif

