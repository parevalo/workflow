#!/bin/bash -l

# This script runs the change map  script for a given list of scenes.
# It also automates finding the example image, which is the first 
# Landsat 7 image found in the folder

# List of scenes to be processed

## UPDATE THIS TO DO IT PER BIANNUAL PERIOD, PER SCENE. THEN CREATE A MODIFIED
## VERSION OF THE MOSAIC SCRIPT. WITH THESE NEW MAPS, I CAN EXTRACT THOSE VALUES
## FOR THE SAMPLES (FIRST AND SECOND) AND KNOW FOR SURE WHEN YATSM DETECTED THE
## CHANGE. 


scn_list="003058 003059 004057 004058 004059 004061 004062 005057 005058 \
          005059 005060 005061 006058 006059 006060 006061 007058 007059 \
          007060 007061 008058 008059 008060 009059 009060"

# List of scenes for which recent subfolders were created, and then the rest
recent="004058 004061 005058 006059 006060 007060 007061"
rest="003058 003059 004057 004059 004062 005057 005059 \
      005060 005061 006058 006061 007058 007059 \
      008058 008059 008060 009059 009060"


# General settings: path to base folder, step to do the change map.

rootdir=/projectnb/landsat/projects/Colombia/images
start_yr=2001
end_yr=2013
step=2

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
        res_path=$scn_path/Results/M3/TSR_nocomm
        outdir=$scn_path/Results/M3/chg_maps
    esac

    case $recent in
    *"$s"*)
        echo "$s will run with the 'recent' settings"
        ts_path=$scn_path/images/recent_period
        res_path=$scn_path/Results/M3_recent_period/TSR_nocomm
        outdir=$scn_path/Results/M3_recent_period/chg_maps
    esac

    # Check if folder exists, create it otherwise
    if [ ! -d $outdir ]; then
        echo "Creating output folder $outdir"
        mkdir -p $outdir
        cd $outdir
    else
        echo "Changing directory to the specified path"
        cd $outdir
    fi

    # Find example image
    img=$(find $ts_path -maxdepth 1 -type d -name "*LE7*" | head -1 )
    example_img=$(basename $img)
    
    # Set variables for map script
    img_path=$ts_path/$example_img/$example_img"_stack" 
 
    # Run change map script. End date is non inclusive. 
    for yr in $(seq -w $start_yr 2 $end_yr); do
        yr2=$(expr $yr + $step)
        qsub -j y -V -N cm_$pt$rw'_'$yr -b y \
         yatsm -v changemap --root $ts_path --result $res_path --image $img_path \
          --magnitude first $yr-01-01 $yr2-01-02 firstchg_$yr-$yr2.tif
    done
done
