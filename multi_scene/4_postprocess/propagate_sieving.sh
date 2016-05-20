#!/bin/bash -l
#$ -V
#$ -j y

# This script back propagates the sieving changes made to the last map
# in the timeseries (e.g. 2015 or 2016) to all of the older maps.

module load python/2.7.5_nopath
module load gdal/1.11.1

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# General settings

rootdir=/projectnb/landsat/projects/Colombia/images

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}
    
    cd $rootdir/$s/Results/M3/Class
    
    # Create mask for sieved areas  (i.e where sieved and original differ)
    gdal_calc.py -A mergedmaps_2016-01-01.tif -B mergedmaps_2016_sieved.tif \
      --outfile=sieved_pixels_2016.tif --calc="A != B" \
      --NoDataValue=0 --type=Byte --co="NBITS=2" --overwrite

    # Propagate the changes
    #for yr in $(seq -w 01 15); do    
    #    gdal_calc.py -A mergedmaps_20$yr"-01-01.tif" -B mergedmaps_2016_sieved.tif \
    #     -C sieved_pixels_2016.tif --outfile=mergedmaps_20yr"_sieved.tif \
    #     --calc="(C == 0)*A + (C == 1)*B"
    #done

done

