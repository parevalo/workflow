#!/bin/bash -l

cd /projectnb/landsat/projects/Colombia/workflow/multi_scene/7_poststratification

qsub -pe omp 4 -j y -b y -N buffer -l mem_total=94G python polygonize_and_buffer.py
 
