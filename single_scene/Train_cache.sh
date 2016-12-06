#!/bin/bash
#$ -V
#$ -l h_rt=48:00:00
#$ -N train_859
#$ -j y

cd /projectnb/landsat/projects/Colombia/images/008059/Results/FIT1/logs

yatsm_cfg=/projectnb/landsat/projects/Colombia/images/008059/Results/FIT1/859_FIT1.yaml
rf_cfg=/projectnb/landsat/projects/Colombia/images/008059/Results/FIT1/RandomForest.yaml

yatsm train $yatsm_cfg $rf_cfg trainRF2.pkl 


