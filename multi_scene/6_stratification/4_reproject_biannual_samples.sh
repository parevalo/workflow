#!/bin/bash -l

# Script that performs the following operations:
# 1) Split samples in the east and west zones into separate shapefiles
# 2) Rasterize samples. Done bc it is easier to reproject rasters to make sure
# that the grids match
# 3) Reprojects sample in east zone to UTM19. The output file is "snapped" to the 
# grid of the east zone in the same way it was done for the mosaics
# so that it matches the grid individual scenes in that zone for TS interpretation
# 4) Poligonize those reprojected samples
# 5) TODO: Add fields to those samples (fiona/shapely script, iterate over scenes
# to add scene info"

# cd to the biannual samples folder
cd /projectnb/landsat/projects/Colombia/biannual_samples_june2017
zpath=/projectnb/landsat/projects/Colombia/Mosaics/M3
vpath=/projectnb/landsat/projects/Colombia/vector
scriptpath=/projectnb/landsat/projects/Colombia/workflow/multi_scene

# Set year and step info
first_yr=1
last_yr=4
step=2
fy=$(printf %02d $first_yr)
ly=$(printf %02d `expr $last_yr - $step`)

for yr in $(seq -w $fy $step $ly); do
    # Set vars and names
    yr2=$(printf %02d `expr $yr + $step`)    

    # Split samples
    sname=sample_$yr"_"$yr2
    qsub -j y -b y -V -N split_$yr \
        $scriptpath/6_stratification/split_samples.py $sname".shp" \
         $vpath/amazon_selection_EAST_UTM18_clipped.shp \
          $vpath/amazon_selection_WEST_UTM18.shp \
           $sname"_east.shp" $sname"_west.shp"

    # Rasterize east samples
    qsub -j y -b y -V -N rast_$yr -hold_jid split_$yr \
        gdal_rasterize -a ID -a_nodata -9999 \
         -te $(gdal_extent $zpath/eastUTM18_2016.tif) -tr 30 30 -ot Int16 \
          -co COMPRESS=PACKBITS $sname"_east.shp" \
           $sname"_east_UTM18N.tif"

    # Reproject east samples to UTM19
    qsub -j y -b y -V -N reproj_$yr -hold_jid rast_$yr \
        gdalwarp -co COMPRESS=PACKBITS -wt Int16 -overwrite \
         -te $(gdal_extent $zpath/eastUTM19_2016.tif) -te_srs EPSG:32619 \
          -t_srs EPSG:32619 -tr 30 -30 $sname"_east_UTM18N.tif" \
           $sname"_east_UTM19N.tif"

    # Poligonize those files. shapes are not overwritten with this tool
    qsub -j y -b y -V -N polig_$yr -hold_jid reproj_$yr \
        gdal_polygonize.py $sname"_east_UTM19N.tif" -f '"ESRI Shapefile"' \
         $sname"_east_UTM19N.shp" $sname"_east_UTM19N" ID

    # Assign path-row info to these and save as new. Submit one for east
    # and one for west. west doesn't need hold_jid
    qsub -j y -b y -V -N pr_W_$yr \
        $scriptpath/other/assign_path-row.py $sname"_west.shp" \
         $vpath/WRS2_amazon_selection_UTM18_no-overlap.shp \
          $sname_"west_PR.shp"  
    
    qsub -j y -b y -V -N pr_E_$yr -hold_jid polig_$yr \
        $scriptpath/other/assign_path-row.py $sname"_east_UTM19.shp" \
         $vpath/WRS2_amazon_selection_UTM19_no-overlap.shp \
          $sname_"east_PR.shp"  
done


