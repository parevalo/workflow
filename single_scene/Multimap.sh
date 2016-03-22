#!/bin/bash

cd /projectnb/landsat/projects/Colombia/images/008059/Results/FIT1/Class

ts_path=/projectnb/landsat/projects/Colombia/images/008059
example_img=

for yr in $(seq -w 0 12); do
    qsub -j y -V -N ymap_20$yr -b y  \
	  yatsm -v map --root $ts_path/images --result $ts_path/Results/FIT1/TSR\
	  --image $ts_path/images/$example_img --after \
	  class 20$yr-01-01 Class_20$yr-01-01_mult_aft.tif
done

