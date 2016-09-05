#!/bin/bash -l

# Script to reproject the mosaics to the IDEAM coordinate system and grid
# in order to allow the comparison of areas. Script is working fine, 
# except because the output file has a 2 meter offset in the lower left Y
# coordinate for an unknown reason. Also, the annual IDEAM files have a 
# sligtly different origin than the others 

export GDAL_DATA=/usr3/graduate/parevalo/miniconda2/envs/GDAL_ENV/share/gdal
cd /projectnb/landsat/projects/Colombia/Mosaics/M3/IDEAM 

# Submit the warping job with reprojection and resampling. Use 15 as input 
# and output nodata, because the "other to other" class is still labeled as zero

for f in $(ls custom__*.tif); do
    qsub -V -b y -j y -N reproj_IDEAM \
     gdalwarp -co COMPRESS=PACKBITS -co NBITS=4 -wt Byte \
      -te $(gdal_extent cambio_2002_2004_v6_amazonia.tif) -te_srs EPSG:3116 \
       -t_srs EPSG:3116 -tr 30.717 -30.2624 -srcnodata 15 -dstnodata 15 \
        -overwrite $f $(basename $f .tif)"_MAGNA.tif"
done    


