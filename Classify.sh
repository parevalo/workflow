cd /projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/logs

cfg_path=/projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/758_FIT1.yaml
algo_path=/projectnb/landsat/projects/Colombia/images/007058/Results/FIT1/trainRF1.pkl

njob=400

for job in $(seq 1 $njob); do
    qsub -j y -V -l -N yclass_$job -b y \
	  yatsm -v classify $cfg_path $algo_path $job $njob
  done
