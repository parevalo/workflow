#!/bin/bash -l

# Submit crosstabulation script created by Chris Holden. Requires docopt

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal

cd /projectnb/landsat/projects/Colombia/Mosaics/M3

# Submit the job

qsub -j y -b y -V -N crosstab ~/misc/maps/crosstab.py --attribute=strata -v \
 sample_may2016_UTM18N.tif final_sample_merge_UTM18N.shp crosstab.csv

