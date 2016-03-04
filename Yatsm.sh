#!/bin/bash

cd /projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/logs
yconfig=/projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/758_FIT1.yaml

njob=400

for job in $(seq 1 $njob); do
    qsub -j y -V -N yatsm_$job -b y \
	 yatsm line --check_cache $yconfig $job $njob
  done

