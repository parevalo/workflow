#!/bin/bash
#$ -V
#$ -N change_758
#$ -j y

cd /projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/Class

ts_path=/projectnb/landsat/projects/Colombia/images/007058
example_img=

yatsm -v changemap --root $ts_path/images --result $ts_path/Results/FIT1/TSR\
 --image $ts_path/images/$example_img --magnitude 2003-01-04 2008-12-27 Change03-08.tif
