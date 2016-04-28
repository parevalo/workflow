#!/bin/zsh
#$ -b y
#$ -V
#$ -j y
#$ -N mosaic_compare

# Script to submit the python scripts that compare either two or multiple
# mosaic files

cd /projectnb/landsat/projects/Colombia/Mosaics/M1A

# Compare two mosaics
python ~/Scripts/Python/MosaicCompare_two.py

# Compare multiple mosaics (sequential dates)
#for y in $(seq -w 00 02 02); do
#	for r in $(seq -w 01 02 03); do
#		qsbin script mosaic_20$y.tif mosaic_20$r.tif
#	done
#done

