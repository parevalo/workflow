#!/bin/bash

# Script to count the total number of images per scene

cd /projectnb/landsat/projects/Colombia/images 

for folder in $(ls -d */); do
    echo $folder, $(ls -d $folder"images/"L* | wc -l) 

done
