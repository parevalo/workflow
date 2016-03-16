#!/bin/bash

cd /projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/logs

cfg_path=/projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/758_FIT1.yaml
algo_path=/projectnb/landsat/projects/Colombia/classifiers/859_trainRF1.pkl

njob=200

for job in $(seq 1 $njob); do
    qsub -j y -V -N yclass_$job -b y \
	  yatsm -v classify $cfg_path $algo_path $job $njob
done

