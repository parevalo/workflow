### This script has three main sections:
### 1) Read reference strata shapefile and calculate labels per year
### 2) Read yearly strata raster and extract values to the SpatialPoints object (i.e. shapefile) 
### 3) Calculate reference class area proportions and variance of ref samples per strata
### 4) Calculate unbiased standard error for proportion of reference class areas 

require(rgdal)
require(raster)

## 1) READ SHAPEFILE AND CALCULATE REFERENCE LABELS PER YEAR
# This section takes the sample shapefile and assigns the proper strata for each year, depending
# on when the model break happened. For now it's written to filter the pixels with only ONE CHANGE (1020 pts) 

# Read shapefile with reference strata and change info
setwd("C:/OneDrive/Lab/sample_may2016/interpreted_w_strata_23062016")
full_samples <- readOGR(".", "final_extended_sample_merge_UTM18N_point")

# Subset to only use records with one or no change
samples <- full_samples[full_samples@data$NUMCHANGES <= 1, ]

# We need to convert the dates to characters bc they are factors right now
samples$CHGDATE <- as.character(samples$CHGDATE)

# Set up important loop variables
rows = nrow(samples)
start = 2001
end = 2016
years = seq(start, end)
field_names = paste("ref_", years, sep="")

# Initialize empty vectors to store the data
df = vector()
field = vector()

# Iterate over years ("columns")
for (y in years){
  # Iterate over rows
  for (row in 1:rows){
      if (is.na(samples$CHGDATE[row]) == TRUE) {
          field[row] = samples$strata[row]
      } else {
        # Find FIRST date and make the comparison. 
        chgdate = na.omit(as.numeric(unlist(strsplit(samples$CHGDATE[row], "[^0-9]"))))
        if (y <= chgdate[1]){
          field[row] = samples$CODE1[row]
        } else {
          field[row] = samples$strata[row]
        }
      }
  }
  # Find and replace values of 7 with 3, because that's how this strata was defined
  field[which(field == 7)] = 3
  # Paste vectors iteratively
  df <- as.data.frame(cbind(df, field))
}

names(df) = field_names

# Attach table to shapefile 
samples@data[,field_names] <- df
#writeOGR(samples, "sample_yearly_strata", "sample_yearly_strata", driver="ESRI Shapefile", overwrite_layer = T)

## 2) READ YEARLY STRATA RASTERS AND EXTRACT THEIR VALUES TO THE SHAPEFILE
# Create a list with the raster names
years_short = seq(02,16)
rast_names = paste0("final_strata_01_", sprintf("%02d",years_short), "_UTM18N")

# Iterate over names and extract to shapefile
for (r in rast_names){
  map = raster(paste0("C:/test/", r, ".tif"))
  samples = extract(map, samples, sp = TRUE) 
}

# Crosstab reference sample vs final strata (DOESN'T MATCH EXCEL FILE BC THE RECORDS WERE FILTERED ABOVE. 
# If run with original, confusion matrices are the same)
ct = table(samples$final_strata_01_16_UTM18N, samples$strata)

# LOAD the total strata sample size produced from CountValues.py, bc calculating it here with hist() takes forever...
ss = read.csv("C:/test/strata_count.csv", header=TRUE, col.names=c("stratum", "pixels"))
# Classes to be removed/ignored
cr = c(7, 10, 12, 15)
# Filter classes NOT in that list
ss = ss[!(ss$stratum %in% cr),] 

# Calculate total number of samples per ORIGINAL stratum
strata_pixels = aggregate(samples$final_strata_01_16_UTM18N, by=list(samples$final_strata_01_16_UTM18N), length)

## 3) CALCULATE REFERENCE CLASS AREA PROPORTIONS AND VARIANCE OF REFERENCE SAMPLES PER STRATA
# Returns a vector with length equal to the number of reference classes
# Reorganize to make more legible

calc_area_prop = function(strata, reference){

  # Generate numbered sequence from the unique values in the REFERENCE field
  ref_codes = sort(unique(reference))
  ref_ind=seq_along(ref_codes)
  
  # Obtain unique values in the STRATA field
  str_codes = sort(unique(strata))
  str_ind=seq_along(str_codes)
  
  # Check if values in strata are the same as in reference, or check i
  
  # Initialize empty df for proportions per strata, and for sample variance per strata
  ref_prop = data.frame()
  ref_var = data.frame()
  
  # Compare the fields
  for (s in str_ind){
    for (r in ref_ind){
      # Compared fields and get TRUE or FALSE on each row
      cond_bool = strata == str_codes[s] & reference == ref_codes[r]
      # Get row numbers that meet that condition
      ind = which(cond_bool)
      # Get location of rows for the current strata
      str_ref = which(strata == str_codes[s])
      # Get proportion on reference class present on that strata
      ref_prop[s, r] = length(ind)/strata_pixels$x[strata_pixels$Group.1 == str_codes[s]]
      # Calculate variance of current reference code in current strata (needed later)
      # THIS IS PRODUCING NA'S IN 01-02, REF-2002 IN STRATA 8 BC THERE IS ONLY ONE OBSERVATION, HOW TO DEAL
      # WITH THIS?
      ref_var[s, r] = var(cond_bool[str_ref])
    }
  }
  
  # Assign column and row names for easier reference
  rownames(ref_prop) = paste0("strat_",str_codes)
  colnames(ref_prop) = paste0("ref_", ref_codes)
  rownames(ref_var) = paste0("strat_",str_codes)
  colnames(ref_var) = paste0("ref_", ref_codes)
  
  # Calculate reference class proportions (i.e. by columns) using total, original strata areas.
  class_prop = vector()
  # Filter only total sample sizes that are present in the strata for that year
  fss = ss[ss$stratum %in% str_codes,]
  for (r in 1:ncol(ref_prop)){
    # LEAVE THE SUM OF THE ENTIRE STRATA?
    class_prop[r] = sum(fss$pixels * ref_prop[,r])/sum(ss$pixels)
  }
  return(list(class_prop, ref_var, fss, ref_prop)) 
}

# Call the function for every reference year we want and get area proportions and sample variance
# NOTE the double square brackets to allow for substitution
ca = calc_area_prop(samples$final_strata_01_02_UTM18N, samples$ref_2002)

area_prop = data.frame()
ref_var_list = list()
filtered_ss = list()
for (i in (1:length(rast_names))){
  # Compare year strata with year reference. Field names must start at 2002, hence i+1
  out = calc_area_prop(samples[[rast_names[i]]], samples[[field_names[i+1]]]) 
  ap = out[[1]]
  rv = out[[2]]
  fs = out[[3]]
  area_prop = rbind(area_prop, ap)
  ref_var_list[[i]] = rv
  filtered_ss[[i]] = fs
}

# Assign names to make it easier to interprete
rownames(area_prop) = years[2:length(years)]
refcodes = sort(unique(samples$ref_2016))
colnames(area_prop) = refcodes
#Not all the rows in area_prop add to 1, check! Maybe there is some sort of double counting in the labels??

## 4) CALCULATE UNBIASED STANDARD ERROR FOR PROPORTION OF REFERENCE CLASS AREAS
# add function to iterate over each class, then call the fnc for each year until 2015
# get filtered sample size from filtered_ss for each year (only needed if problem with 2002 and 2010 is fixed)

# Formula works correctly, tested with unfiltered sample (i.e 1048 records) and it gives roughly the same CI
se_prop = data.frame()
#Iterate over years
for (y in 1:length(ref_var_list)){
  # Iterate over classes
  for (c in 1:length(refcodes)){
  N=sum(ss$pixels)
  v = 1/N^2 * (sum(ss$pixels^2 * (1 - strata_pixels$x/ss$pixels) * (ref_var_list[[y]][,c] / strata_pixels$x)))
  se_prop[y,c] = sqrt(v)
  }
}

# Reassign names...
rownames(se_prop) = years[2:length(years)]
colnames(se_prop) = refcodes
## Years 2002 and 2010 have NA's because stratum 8 and 13 (respectively) only have 1 value and thus no variance

# Total area in ha
N_ha = N * 30^2 / 100^2
# Calculate area in ha from area proportions
area_ha = area_prop * N_ha
# Calculate confidence interval in ha
area_ci = se_prop * 1.96 * N_ha
#Upper and lower CI
area_ha + area_ci
area_ha - area_ci
area_ci

# Reference sample count per year
yearly_ref_class = data.frame()
for (f in 2:(length(field_names))){
  a = table(samples[[field_names[f]]])
  yearly_ref_class = rbind(yearly_ref_class, a)
}

yearly_ref_class = rbind(0, yearly_ref_class)
yearly_ref_class[1,2:11] = table(samples$ref_2001)
rownames(yearly_ref_class) = years
colnames(yearly_ref_class) = refcodes
