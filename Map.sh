#!/bin/bash
#$ -V
#$ -N Map_758
#$ -j y

cd /projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/Class

ts_path=/projectnb/landsat/projects/Colombia/images/007058
example_img=

yatsm -v map --root $ts_path/images --result $ts_path/Results/FIT1/TSR\
 --image $ts_path/images/$example_img class 2003-01-04 Class2_2003-01-04.tif

