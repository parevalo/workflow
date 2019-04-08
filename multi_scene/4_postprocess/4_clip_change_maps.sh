#!/bin/bash -l

# Script to cut the change  maps to the WRS2 boundary, because there
# are extra areas with less data and lower classification quality
# caused by the bigger footprint of L8 images.

# List of scenes to be processed

scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# List of scenes for which recent subfolders were created, and then the rest
recent="004058 004061 005058 006059 006060 007060 007061"
rest="003058 003059 004057 004059 004062 005057 005059 \
      005060 005061 006058 006061 007058 007059 \
      008058 008059 008060 009059 009060"

# General settings
rootdir=/projectnb/landsat/projects/Colombia/images
poly=/projectnb/landsat/projects/Colombia/vector/WRS2_amazon_selection.shp

# Iterate over scenes

for s in $scn_list; do
    # Get path and row in short version
    pt=${s:2:1}
    rw=${s:4:2}
    scn_path=$rootdir/$s

    # Set scene, TS(images) and results paths for the map script
    case $rest in
    *"$s"*)
        ts_path=$scn_path/images
        outdir=$scn_path/Results/M3/chg_maps
    esac

    case $recent in
    *"$s"*)
        echo "$s will run with the 'recent' settings"
        ts_path=$scn_path/images/recent_period
        outdir=$scn_path/Results/M3_recent_period/chg_maps
    esac

    # cd to the corresponding class folder
    cd $outdir
    
    # Submit the clipping job
    i=1
    for f in $(find ./ -maxdepth 1 -name "*.tif"); do
        # Extract basename and construct out name
        bname=$(basename $f)
        outname=${bname%.*}"_clipped.tif"
        echo "Processing file $bname, output file $outname"
        qsub -j y -V -N clip$pt$rw"_"$i -b y \
         gdalwarp -tr 30 30 -srcnodata 0 -cutline $poly -cl WRS2_amazon_selection \
          -cwhere "'PTRW=$pt$rw'" $bname $outname
        let i+=1
    done
done

