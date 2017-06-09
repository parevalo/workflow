# File with path to necessary files to run the area estimation analysis with the original parameters

# Paths and suffixes/prefixes
savepath = "/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/results/for_nofor/"
lutpath = "/home/paulo/workflow/multi_scene/data/for-nofor_lut_A.csv"
lut_name = "lutA"
orig_stratif = paste0("final_strata_01_16_UTM18N_", lut_name)
rast_prefix = "final_strata_annual_"
rast_suffix = paste0("_UTM18N_", lut_name)
pixcount_suffix = paste0("_pixcount_", lut_name, ".csv")
pixcount_strata = paste0("strata_01_16", pixcount_suffix)

# Reference strata names and other variables
strata_names = c("Stable forest", "Stable non-forest", "Forest loss", "Forest gain")
cr = c(0,seq(5,15)) # Classes to ignore from the loaded area count tables
step = 1  # Number of years to do the analysis over

# Run estimation with additional number of rows with perfect forest classification
add_samples = FALSE 
nfor = 500 # of perfect forest samples in stable forest class
nbuf = 0 # of perfect forest samples in buffer (if applicable)  