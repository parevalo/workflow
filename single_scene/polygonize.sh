#!/bin/bash
#$ -V
#$ -N polygonize
#$ -m e

module load gdal/1.10.0

path=/projectnb/landsat/projects/Colombia/images/007058/Results/FIT3/Class

cd $path

gdal_polygonize.py Class_2004-01-15_Train5_updt2.tif -f "ESRI Shapefile" Class_2004-01-15_Train5Updt2.shp Class
