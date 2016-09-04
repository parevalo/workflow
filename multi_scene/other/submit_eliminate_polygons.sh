#!/bin/bash -l
#$ -b y
#$ -V
#$ -j y
#$ -N elim_poly


set -x

module load python/2.7.5_nopath
module load qgis/2.6.1

export PYTHONPATH=/share/pkg/qgis/2.6.1/install/share/qgis/python

# Polygonize with GDAL before using the QGIS script

cd /projectnb/landsat/projects/Colombia/vector

# Use QGIS to process the vector file
python /projectnb/landsat/projects/Colombia/workflow/multi_scene/other/eliminate_small_polygons.py 

echo "Something"

# Use gdal to rasterize that polygon

# Clip and reproject to a commmon extent in order to be able to compare

