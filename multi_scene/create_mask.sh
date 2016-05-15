#!/bin/bash -l

module load python/2.7.5_nopath
module load gdal/1.11.1

cd /projectnb/landsat/projects/Colombia/images/007059/Results/M2/Class

#for yr in $(seq -w 01 01); do
#    qsub -V -j y -b y gdal_calc.py -A 20$yr"-01-01_seq.tif" \
#    --outfile=$yr"_2-3.tif" --calc="logical_or(A == 2, A == 3)" \
#    --NoDataValue=3 --type=Byte --co="NBITS=2"
#done

# Not working with qsub for some reason
gdal_calc.py -A ClassM2B_2001-01-01_crop_class.tif \
 --outfile=2001_2-3.tif --calc="logical_or(A == 2, A == 3)" \
 --NoDataValue=3 --type=Byte --co="NBITS=2"

