# Steps followed to obtain the maps and the unbiased area estimation

The steps in this workflow were performed in the Shared Computer Cluster at Boston University. The fitting of time series models, classification and map creation was performed using YATSM (https://github.com/ceholden/yatsm). Most of the scripts are loops that submit the same jobs for each of the landsat scenes used in the study. Most scripts were written to run under GDAL>=2.1. The environment file can be used with conda to recreate a working environment for all the Python dependent scripts, and works with yatsm v0.6.2.

## 1_fit 

1. Submit jobs to generate a CSV list of dates and filepaths of each image available for each Landsat scene. This assumes that all images are stored following the file location convention specified in the README.md file. The images are assumed to be stacks of surface reflectance with the same extent per scene. Bands are assumed tto be stacked in the following order: Blue, Green, Red, NIR, SWIR1, SWIR2, brightness temperature and FMask. Preprocessing steps can be found [here](https://github.com/parevalo/landsat_process) 
2. Cache images on disk to make them faster to read, using yatsm cache.
3. Fit segments on time series per pixel, using yatsm line.

## 2_train-class

1. Rasterize the training data collected in shapefile format. 
2. Cache the training data in individual cache files.
3. Merge these individual cache files into a single file.
4. Use the single merged cache fil to train a RandomForest classifier.
5. User that single classifier to classify the fitted segments of all the scenes.

## 3_map

1. Create the land cover maps based on the classified results from previous step.
2. Merge the TWO resulting land cover maps. Technically only the results of the model M3 were used to create 
 the maps, but those records were classified two times, using two different classifiers. A map was created for
 each  one of those classifications, then those two maps were merged. The first classifier was obtained from 
 model M1 (which has minor differences from model M3) and the second came from using the training data over
 the M3 records. The first was used because it produced better results for forest mapping, the second was used
 because it produced better results to distinguish between pastures and grasslands. Therefore, the map from the first approach was preserved, but areas mapped as grasslands were replaced by whatever label was assigned in
 the second approach, as long as there was no change in the entire period. This procedure is more of a stopgap  to distinguish pastures and grasslands properly, but in the future simpler methods should be favored.
 
3. Create changemap between 2001 and 2016, required for the sieving operation (see below)
4. Modify the metadata to make the maps easier to use with gdalcalc.

## 4_postprocess

1. Copy files in UTM18 zone to follow naming convention and make it easier to process them.
2. Perform a masking and sieving operation on the files without modifying the pixels that experienced any change during the entire study period.
3. Propagate the sieving to all of the maps in the study period
4. Clip maps to their WRS2 footprint to make them easier to stack.

## 5_mosaic

1. Mosaic the scenes present in each of the two UTM zones in the area (18 and 19).
2. Reproject the east zone mosaic (in UTM19) to UTM18.
3. Merge the two mosaics into a single one.

## 6_stratification

1. Clip mosaics to the amazon region boundary
2. Create strata based on the classes shown in the README file. This can be done either annually
 or biannually and including a buffer class to capture omission errors around areas where forest to
 pasture conversion is expected. Both scripts use the `create_strata.py` script.
3. Submit the script to draw the stratified random samples either for the stratification of the entire period 
 (original, 2001 - 2016) or the biannual stratification. Both scripts use the `sample_map.py` script. 
4. Pre-interpretation reprojection and processing of the samples that intersect east zone to UTM19 in order 
 to aid the interpretation using TSTools. This process was made mostly manually for the origina samples but
 was automated for the biannual samples.

 The biannual script uses the `split_samples.py` script to separate those in east and west  zones and then 
 assigns the correct path rows using the `assign_path-row.py` script and a customized version  of the scenes 
 footprints to deal properly with the overlap zones. The process of path-row assignment may result in samples 
 being located in the line boundary between two scenes, in which case the sample will be duplicated and a 
 manual decision will have to be made regarding which scene to use. Therefore, it is possible that after 
 this step the total number of samples is temporarily increased. 
5. Calculation of areas and accuracies for each biannual period using unbiased estimators. 

All python scripts are located in /helper_scripts.


## 7_poststratification

0. Generate a poststratified version of the current strata that creates a buffer into the forest around
 forest to pasture areas, in order to capture omissions errors.
1. Calculate areas for each class in each of the strata maps (annual, bianual, original). Uses the script
 `count_pixels.py`. This is required to calculate the unbiased area estimates and accuracies.
calculate_strata_per_year.R along with other things.
2. calculate_strata_per_year.R: This script can be run from start to end, specifying multiple input files
 with different configurations that refer to stratification schemes (e.g original, forest-noforest, original
 where class 8 and 9 are merged together).The main purpose is to calculate the unbiased area estimations, 
 standard errors and margins of error per period and create the plots. Other sections to create useful tables,
 as well as other calculations and plots need to be rewritten. 


- The strata_calc_config folder contains different configuration files that allowed calculating the
areas and accuracies for multiple scenarios (e.g. with and without buffer, using different buffer 
sizes, etc). Most of them were only used during the exploration phase. 

TODO: Organize scripts to replicate the figures in the Arevalo et al. 2019 paper in a separate, 
organized folder.

## data

This folder contains lookup tables used for map reclassification and assignment of vector labels.
It also contains the test data provided in Stehman et al. (2014), used to verify the scripts
that calculte areas, accuracies and their standard errors using indicator functions (ratio estimator).

## reference
This folder contains the logs from the creation of the original samples and the number of images per year.

## other
This folder contains miscelaneous scripts and files that were used at some point in the exploratory
phase of the analysis and are not needed for the calculation of areas and accuracies. 
