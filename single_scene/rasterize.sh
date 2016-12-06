#!/bin/bash
#$ -V
#$ -N rasterize
#$ -j y

vector=/projectnb/landsat/projects/Colombia/Training/007058/Training1.shp
raster=/projectnb/landsat/projects/Colombia/images/007058/images/Training1.tif

gdal_rasterize -a Class -init 0 -a_nodata 0 -te  564285 213585  818715 447015\
 -tr 30 30 -ot UInt16 $vector $raster
