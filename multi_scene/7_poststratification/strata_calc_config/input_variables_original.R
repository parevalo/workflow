# File with path to necessary files to run the area estimation analysis with the original parameters

# Paths and suffixes/prefixes
lutpath = "/home/paulo/workflow/multi_scene/data/original_lut.csv"
lut_name = "original_lut"
orig_stratif = "final_strata_01_16_UTM18N"
rast_prefix = "final_strata_annual_"
rast_suffix = "_UTM18N"
pixcount_suffix = "_pixcount.csv"
pixcount_strata = paste0("strata_01_16", pixcount_suffix)
savepath = "/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/original_sampling/results/original/"

# Reference strata names (Ordered names of the classes that will get their areas estimated based on the reference samples)
strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                 "Stable pasture-cropland", "Stable secondary forest", "Stable water", "Forest to pasture", 
                 "Forest to secondary forest", "Gain of secondary forest", "Loss of secondary forest")

# List of original strata names (Ordered names of all the classes from the original stratification, including those for which
# areas won't be estimated)
orig_strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                      "Stable pasture-cropland", "Stable secondary forest", "Stable water", "Forest to pasture", 
                      "Forest to secondary forest", "Gain of secondary forest", "All to unclassified", "Loss of secondary forest")

cr = c(7, 10, 12, 15) # Classes to ignore from the loaded area count tables
step = 2   # Number of years to do the analysis over

# Run estimation with additional number of rows with perfect forest classification
add_samples = FALSE 
nfor = 0 # of perfect forest samples in stable forest class
nbuf = 0 # of perfect forest samples in buffer (if applicable)  
