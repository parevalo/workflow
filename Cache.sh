cd /projectnb/landsat/projects/Colombia/images/006058/Results/FIT1/logs

cfg_path=/projectnb/landsat/projects/Colombia/images/006058/Results/FIT1/658_FIT1.yaml

njob=400

for job in $(seq 1 $njob); do
    qsub -j y -V -l -N cache_$job -b y yatsm -v cache $cfg_path $job $njob
done

