#!/bin/bash -l
# Script intended to test the polygonize and buffer python script. In the end
# it was replaced with a GDAL script that is much faster and easier to use, but
# left here for future reference in case it's needed


cd /projectnb/landsat/projects/Colombia/workflow/multi_scene/7_poststratification

qsub -pe omp 4 -j y -b y -N buffer -l mem_total=94G python polygonize_and_buffer.py
 
