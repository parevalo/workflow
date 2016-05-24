#!/bin/bash -l

module load python/2.7.5_nopath
module load gdal/1.11.1

# Script to reproject the maps from UTM19 to UTM18, where most of the changes
# are concentrated 

# List of scenes in UTM19

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058"

# General settings

pref=mergedmaps_20
suf=_final.tif

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}

    # cd to the corresponding class folder
    cd /projectnb/landsat/projects/Colombia/images/$s/Results/M3/Class
    
    # Submit the clipping job
    for yr in $(seq -w 16 16); do
        qsub -V -b y -j y -N reproj_$pt$rw"_"$yr \
        gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte \
        -t_srs EPSG:32618 -tr 30 30 $pref$yr$suf $pref$yr"_final_UTM18N.tif"
    done
done

