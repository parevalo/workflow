#!/bin/bash

# Script to count the total number of images per year

cd /projectnb/landsat/projects/Colombia/images 

# Look for the year preceded by a number to avoid erroneous count in 859
# bc of a different folder naming scheme

for yr in $(seq 2001 2015); do
    echo $yr, $(ls -d 00*/images/*[0-9]$yr*/ | wc -l) 

done
