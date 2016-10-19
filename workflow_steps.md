# Steps followed to obtain the maps and the unbiased area estimation

The steps in this workflow were performed in the Shared Computer Cluster at Boston University. The fitting of time series models, classification and map creation was performed using YATSM (https://github.com/ceholden/yatsm). Most of the scripts are loops that submit the same jobs for each of the landsat scenes used in the study.

## 1_fit 

1. Submit jobs to generate a CSV list of dates and filepaths of each image available for each Landsat scene. This assumes that all images are stored following the file location convention specified in the README.md file. The images are assumed to be stacks of surface reflectance with the same extent per scene. Bands are assumed tto be stacked in the following order: Blue, Green, Red, NIR, SWIR1, SWIR2, brightness temperature and FMask. Preprocessing steps can be found [here](https://github.com/parevalo/landsat_process) 
2. Cache images on disk to make them faster to read, using yatsm cache.
3.  Fit segments on time series per pixel, using yatsm line.

## 2_train-class

1. Rasterize the training data collected in shapefile format. 
2. Cache the training data in individual cache files.
3. Merge these individual cache files into a single file.
4. Use the single merged cache fil to train a RandomForest classifier.
5. User that single classifier to classify the fitted segments of all the scenes.

## 3_map

1. Create the land cover maps based on the classified results from previous step.
2. Merge the TWO resulting land cover maps (I used two classifiers, need to add documentation on that)
3. Create changemap between 2001 and 2016, required for the sieving operation (see below)
4. Modify the metadata to make the maps easier to use with gdalcalc.

## 4_postprocess

1. Copy files in UTM18 zone to follow naming convention and make it easier to process them.
2. Perform a masking and sieving operation on the files without modifying the pixels that experienced any change.
3. Propagate the sieving to all of the maps in the study period
4. Clip maps to their WRS2 footprint to make them easier to stack.

## 5_mosaic

1. Mosaic the scenes present in each of the two UTM zones in the area (18 and 19).
2. Reproject the east zone mosaic (in UTM19) to UTM18
3. Merge the two mosaics into a single one.

## 6_stratification

1. Clip mosaics to the amazon region boundary
2. Create strata based on the classes shown in the README file.
3. Submit the script to draw the stratified sample based on the parameters in the `notes.txt` file.

## 7_poststratification

TODO

