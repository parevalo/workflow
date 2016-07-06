## This script takes the sample shapefile and assigns the proper strata for each year, depending
## on when the model break happened. The way it's written now it only works for the FIRST change, 
## so the input shapefile must be that with samples with only stable or 1 change pixels
## The code is not very efficient (too many loops and conditionals) but it does the job

require(rgdal)
require(raster)

## READ SHAPEFILE AND CALCULATE REFERENCE LABELS PER YEAR

setwd("C:/OneDrive/Lab/sample_may2016/interpreted_w_strata_23062016")
full_samples <- readOGR(".", "final_extended_sample_merge_UTM18N_point")

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

# Attach table to shapefile and write out, with overwrite (OPTIONAL)
samples@data[,field_names] <- df
#writeOGR(samples, "sample_yearly_strata", "sample_yearly_strata", driver="ESRI Shapefile", overwrite_layer = T)

## READ YEARLY STRATA RASTERS AND EXTRACT THEIR VALUES TO THE SHAPEFILE
# Create a list with the raster names
years_short = seq(01,16)
rast_names = paste("final_strata_01_", sprintf("%02d",years_short), "_UTM18N.tif", sep="")

# Iterate over names and extract to shapefile
for (r in rast_names){
  map = raster(r)
  extract(map, samples, sp = TRUE) 
}

##  CALCULATE UNBIASED AREA PROPORTIONS

# Load or calculate area of strata in pixels and total area here. field STRATA used in this block as example, but
# it needs to contain the actual ORIGINAL stratification (eg. final_strata_01_16)
strata_pixels = aggregate(samples$strata, by=list(samples$strata), length)

# Mock MAP class, REPLACE WITH ACTUAL FIELD NAME!
samples$map_2006 = samples$ref_2006

## Calculate sample size per stratum. Count how many times a single REFERENCE class is contained on EACH stratum
calc_area_prop = function(strata, reference){

  # Generate numbered sequence from the unique values in the REFERENCE field
  ref_codes = sort(unique(reference))
  ref_ind=seq_along(ref_codes)
  
  # Obtain unique values in the STRATA field
  str_codes = sort(unique(strata))
  str_ind=seq_along(str_codes)
  
  # Initialize empty lists for total proportions per strata/reference, and for references per strata
  ref_props = list()
  p = vector()
  
  # Compare the fields
  for (s in str_ind){
    for (r in ref_ind){
      # Compared fields and get TRUE or FALSE on each row
      cond_bool = strata == str_codes[s] & reference == ref_codes[r]
      # Get row numbers that meet that condition
      ind = which(cond_bool)
      # Get location of rows for the current strata
      str_ref = which(strata == str_codes[s])
      # Calculate variance of current reference code in current strata (needed later)
      var(cond_bool[str_ref])
      # Get proportion on reference class present on that strata
      p[r] = length(ind)/strata_pixels$x[strata_pixels$Group.1 == str_codes[s]]
      # Save all the reference proportions for this SINGLE strata
      #ref_props = rbind(ref_props, p)
      ref_props[[s]] = p
     
    }
  }
  
  # Calculate reference class proportions, get that for the # of values across the strata (e.g # of values in the p vector).
  # Most likely possible to vectorize but don't see how to make it easily. Put into DF and add class labels or something.
  
  class_prop = vector()
  for (r in ref_ind){
    class_prop[r] = sum(strata_pixels$x * unlist(lapply(ref_props, `[`, r)))/sum(strata_pixels$x)
  }
  return(class_prop)
}

# Call the function for every reference year we want. NOTE the double square brackets to allow for substitution
ca = calc_area_prop(samples$strata, samples$ref_2001)
area_prop = data.frame()
for (f in (field_names)){
  ap = calc_area_prop(samples$strata, samples[[f]])
  area_prop = rbind(area_prop, ap)
}

# Assign names to make it easier to interprete
rownames(area_prop) = field_names
colnames(area_prop) = ref_codes
#No all the rows in area_prop add to 1, check!

# Calculate reference sample variance per stratum. Done with presence-absence vector


# To check which numbers on a list are part of others (x %in% y). To "merge" two lists of numbers, in case it's needed
# to complete the list of available classes union(x,y)
