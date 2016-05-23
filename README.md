# workflow
Repository to include all the basic scripts needed to run all steps of yatsm for each of the scenes

**UPDATE WITH EXPLANATIONS OF ALL OF THE STEPS!**

### Stratification:

* 0 - Other classes to other classes
* 1 - Stable forest
* 2 - Stable grassland
* 3 - Stable urban
* 4 - Stable pasture/cropland
* 5 - Stable regrowth
* 6 - Stable water
* 7 - Stable other (sandbanks, rocks)
* 8 - Forest to pastures
* 9 - Forest to regrowth
* 10 - Forest to all others
* 11 - Pastures to regrowth
* 12 - All other to forest
* 13 - ALL CLASSES to unclassified
* 14 - True NoData

### File locations:

**Images**: `/projectnb/landsat/projects/Colombia/<path><row>/images`

**Image cache**: `/projectnb/landsat/projects/Colombia/<path><row>/images/.cache` 

**Time Series Results**: `/projectnb/landsat/projects/Colombia/<path><row>/M3/TSR`
