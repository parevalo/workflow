## This script takes the sample shapefile and assigns the proper strata for each year, depending
## on when the model break happened. The way it's written now it only works for the FIRST change, 
## so the input shapefile must be that with samples with only stable or 1 change pixels
## The code is not very efficient (too many loops and conditionals) but it does the job

require(rgdal)
require(raster)

## READ SHAPEFILE AND CALCULATE REFERENCE LABELS PER YEAR

setwd("C:/OneDrive/Lab/sample_may2016/interpreted_w_strata_23062016")

##  CALCULATE UNBIASED AREA PROPORTIONSfull_samples <- readOGR(".", "final_extended_sample_merge_UTM18N_point")

# Subset to only use records with one or no change
samples <- full_samples[full_samples@data$NUMCHANGES <= 1, ]


# Split with anything that is not a number, not used here bc is inside the loop for now
#a <- na.omit(as.numeric(unlist(strsplit(samples$CHGDATE, "[^0-9]"))))

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
        if (y < chgdate[1]){
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

## READ YEARLY STRATA RASTERS AND EXTRACT THEIR VALUES TO THE SHAPEFILE
# Create a list with the raster names
years_short = seq(02,16)
rast_names = paste0("final_strata_01_", sprintf("%02d",years_short), "_UTM18N.tif")

# Iterate over names and extract to shapefile
for (r in rast_names){
  map = raster(paste0("C:/test/", r))
  samples = extract(map, samples, sp = TRUE) 
}

# Crosstab reference sample vs final strata (DOESN'T MATCH EXCEL FILE BC THE RECORDS WERE FILTERED ABOVE. If run with original, confusion
# matrices are the same)
ct = table(samples$final_strata_01_16_UTM18N, samples$strata)

# LOAD the total strata sample size produced from CountValues.py, bc calculating it here with hist() takes forever...
ss = read.csv("C:/test/strata_count.csv", header=TRUE, col.names=c("stratum", "pixels"))
# Classes to be removed
cs = c(7, 10, 12, 15)
# Filter classes NOT in that list
ss = ss[!(ss$stratum %in% cs),] 

# Calculate total number of samples per ORIGINAL stratum
strata_pixels = aggregate(samples$final_strata_01_16_UTM18N, by=list(samples$final_strata_01_16_UTM18N), length)

## Calculate reference class area proportions and variance of ref samples per strata. 
## Returns a vector with length equal to the number of reference classes
## Reorganize to make more legible

calc_area_prop = function(strata, reference){

  # Generate numbered sequence from the unique values in the REFERENCE field
  ref_codes = sort(unique(reference))
  ref_ind=seq_along(ref_codes)
  
  # Obtain unique values in the STRATA field
  str_codes = sort(unique(strata))
  str_ind=seq_along(str_codes)
  
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
      ref_var[s, r] = var(cond_bool[str_ref])
    }
  
  }
  
  # Assign column and row names for easier reference
  rownames(ref_prop) = paste0("strat_",str_codes)
  colnames(ref_prop) = paste0("ref_", ref_codes)
  rownames(ref_var) = paste0("strat_",str_codes)
  colnames(ref_var) = paste0("ref_", ref_codes)
  
  # Calculate reference class proportions (i.e. by columns) using total, original strata areas. Vectorize to avoid the ugly loop
  class_prop = vector()
  for (r in 1:ncol(ref_prop)){
    class_prop[r] = sum(ss$pixels * ref_prop[,r])/sum(ss$pixels)
  }
  return(list(class_prop, ref_var)) 
}

# Call the function for every reference year we want and get area proportions and sample variance
# NOTE the double square brackets to allow for substitution
ca = calc_area_prop(samples$final_strata_01_16_UTM18N, samples$ref_2001)
cp = ca[[1]]
rvar = ca[[2]]

area_prop = data.frame()
ref_var_list = list()
for (i in (1:length(field_names))){
  out = calc_area_prop(samples$final_strata_01_16_UTM18N, samples[[field_names[i]]]) #ALSO ITERATE OVER MULTI YEAR STRATA!?!?
  ap = out[[1]]
  rv = out[[2]]
  area_prop = rbind(area_prop, ap)
  ref_var_list[[i]] = rv
}

# Assign names to make it easier to interprete
rownames(area_prop) = years
colnames(area_prop) = ref_codes
#No all the rows in area_prop add to 1, check!

## CALCULATE UNBIASED STANDARD ERROR FOR PROPORTION OF REFERENCE CLASS AREAS
function to iterate overeach class, then call the fnc for each year
N=sum(ss$pixels)
v = 1/N^2 * (sum(ss$pixels^2 * (1 - strata_pixels$x/ss$pixels) * (ref_var_list[[1]][,1] / strata_pixels$x)))
se = sqrt(v)

# To check which numbers on a list are part of others (x %in% y). To "merge" two lists of numbers, in case it's needed
# to complete the list of available classes union(x,y)
