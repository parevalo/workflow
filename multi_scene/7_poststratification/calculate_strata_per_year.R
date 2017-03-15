### This script has seven main sections:
### 0) Set global variables
### 1) Read reference strata shapefile and calculate labels per year
### 2) Read yearly strata raster and extract values to the SpatialPoints object (i.e. shapefile) 
### 3) Calculate reference class area proportions and variance of ref samples per strata
### 4) Calculate unbiased standard error for proportion of reference class areas 
### 5) Create some useful tables
### 6) Create some useful plots

require(rgdal)
require(raster)
require(reshape2)
require(ggplot2)
require(gtable)
require(grid)
require(gridExtra)
require(xtable)
require(matrixcalc)

##############################################################################################################
#0) SET VARIABLES/FOLDERS

# Working directory and result (aux) files from the cluster. Move them to Onedrive and exclude from sync!
wd = "C:/OneDrive/Lab/area_calculation/final_sample"
auxpath = "C:/test/"

if( .Platform$OS.type == "unix" )
    wd = "/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/final_sample"
    auxpath = "/media/paulo/785044BD504483BA/test/"
    stratpath = "/home/paulo/workflow/multi_scene/7_poststratification/"
    
setwd(wd)
source(paste0(stratpath, "functions.R"))
source(paste0(stratpath, "input_variables_original.R")) # CHANGE THIS FILE TO RUN WITH OTHER INPUT PARAMETERS!


# Set up important global variables
start = 2001
end = 2016

years = seq(start, end, 1) 
ref_names = paste("ref_", years, sep="")

# Grid table theme, only used to display some ancillary tables
tt=ttheme_default(core=list(fg_params=list(font="Times", fontface="plain", fontsize=14)),
                  colhead=list(fg_params=list(font="Times", fontface="bold", fontsize=14, parse=TRUE)),
                  rowhead=list(fg_params=list(font="Times", fontface="plain", fontsize=14)))

#############################################################################################################
#1) READ SHAPEFILE AND CALCULATE REFERENCE LABELS PER YEAR
# This section takes the sample shapefile and assigns the proper strata for each year, depending
# on when the model break happened. In this version we use all of the points because when done 
# annually, it does not matter how many changes we found in the TS.  

# Read shapefile with reference strata and change info. 
samples <- readOGR(".", "final_extended_sample_merge_UTM18N_point")

# We need to convert the dates to characters bc they are factors right now
samples$CHGDATE <- as.character(samples$CHGDATE)

# Convert reference sample info into vectors containing proper labels per each of
# the years in our time period. This takes into account that we could be doing the
# analysis every two years or more. 

# Initialize empty vectors to store the data, current change and class code (start with the first!) 
# and name of the fields that store the class codes. Load LUT.
# df stores the labels per year with strata label when there is a change
# df2 stores only the labels per year, with no strata info.
# This whole massive loop could be rewritten using the break_calc function used
# for the break analysis.

rows = nrow(samples)
df = vector()
df2 = vector()
field = vector()
field2 = vector()
codelist = c("CODE1", "CODE2", "CODE3", "CODE4")
lut = read.table(lutpath, header = T, sep = ",")

# Iterate over rows (samples, easier this way)
for (row in 1:rows){
  current_change = 1
  current_code=1
  # Iterate over years
  for (i in 1:length(years)){
  
      # If there are no changes, just use the only class code for both years
      if (is.na(samples$CHGDATE[row]) == TRUE) {
        field[i] = calc_strata(samples$CODE1[row], samples$CODE1[row], lut)
        field2[i] = samples$CODE1[row]
      }
    
      # If there is a change, compare each year to the current change year and update that one accordingly
      else {
        # Get the current date of change (initialized as the first change available)
        chg_date = unlist(strsplit(as.vector(samples$CHGDATE[row]), ','))[current_change] #742
        # Extract the YEAR from the date of change
        chg_year = na.omit(as.numeric(unlist(strsplit(chg_date, "[^0-9]"))))[1]
        
        # If we haven't reached change year yet
        if (years[i] < chg_year) {
          field[i] = calc_strata(samples@data[codelist[current_code]][row,], samples@data[codelist[current_code]][row,], lut) 
          field2[i] = samples@data[codelist[current_code]][row,]
        } 
        
        # If we JUST reached a change year
        else if (years[i] == chg_year) {
          field[i] = calc_strata(samples@data[codelist[current_code]][row,], samples@data[codelist[current_code+1]][row,], lut) 
          field2[i] = samples@data[codelist[current_code+1]][row,]
          # Check if we haven't reached the max number of recorded changes
          if (current_change < samples$NUMCHANGES[row]) {
            current_change = current_change + 1 
          }
          current_code = current_code + 1
        }
        
        # If we went past the last change date, use the last code 
        else if (years[i] > chg_year) {
          field[i] = calc_strata(samples@data[codelist[samples$endcodecol[row]]][row,], samples@data[codelist[samples$endcodecol[row]]][row,], lut)
          field2[i] = samples@data[codelist[current_code]][row,]
        }
      } 
  }
  df <- as.data.frame(rbind(df, field))
  df2 <- as.data.frame(rbind(df2, field2))
  
}

names(df) = ref_names
names(df2) = ref_names

# Recalculate "original stratification", needed when we use a LUT different than the original
strata = vector()
for (row in 1:rows){
  strata[row] = calc_strata(samples@data["CODE1"][row,], samples@data[codelist[samples$endcodecol[row]]][row,], lut)
}

# Attach table to shapefile 
samples@data[,ref_names] <- df
#writeOGR(samples, "sample_yearly_strata", "sample_yearly_strata", driver="ESRI Shapefile", overwrite_layer = T)

#############################################################################################################
# 2) READ ORIGINAL AND ANNUAL STRATA RASTERS AND EXTRACT THEIR VALUES TO THE SHAPEFILE
# Also calculate original strata size, weights, area proportions, and accuracies

# Iterate over names and extract to shapefile. We need these to calculate 
# confusion matrices per year. We DON'T NEED THESE for the area estimation.
# If we want to create confusion matrices for aggregated years, then we NEED TO
# CREATE THOSE RASTERS FIRST e.g. (2001-2003 and so on)
map_names = character()
short_years = substr(years, 3,4) # Get years in two digit format
for (y in 1:(length(years)-1)){
  map_names[y] = paste0(rast_prefix, short_years[y], "_", short_years[y+1], rast_suffix)
  map = raster(paste0(auxpath, map_names[[y]], ".tif"))
  samples = extract(map, samples, sp = TRUE) 
}

# Read original stratification. Done last so that the reference and map fields are contiguous
samples = extract(raster(paste0(auxpath, orig_stratif, ".tif")), samples, sp=TRUE)

# Get unique ref codes and map codes for all the years. Also get unique codes from
# stratification layer (because it may have a buffer that the individual maps don't have)
# Then create a single set of unique codes to be used in the confusion matrices
ref_codes = sort(unique(unlist(samples@data[ref_names])))
map_codes = sort(unique(unlist(samples@data[map_names])))
strat_codes = sort(unique(unlist(samples@data[orig_stratif])))
class_codes = sort(union(ref_codes, map_codes))
class_codes = sort(union(class_codes, strat_codes))

# Crosstab final strata and reference strata (for the same period 01-16) 
# Returns a square matrix with all the reference and map codes in the samples

ct = calc_ct(samples[[orig_stratif]], strata, class_codes)

# Load mapped areas for each individual stratification map (total strata sample size) produced from CountValues.py, 
# REQUIRED for comparison between mapped and estimated areas.
mapped_areas_list = list()
filenames = dir(auxpath, pattern=(paste0("*", pixcount_suffix)))
for(i in 1:length(filenames)){
  mapped_areas_list[[i]] = read.csv(paste0(auxpath,filenames[i]), header=TRUE, col.names=c("stratum", "pixels"))
}

# Load mapped area of the ORIGINAL Stratification (e.g. 01-16)

ss=read.csv(paste0(auxpath, pixcount_strata), header=TRUE, col.names=c("stratum", "pixels"))

# Filter classes NOT in the list of classes from the pixel count files to be ignored. 
# TODO: This probably has to be done to all csv files if we decide to use them
ss = ss[!(ss$stratum %in% cr),] 

# Calculate total number of samples per ORIGINAL stratum. 
# ASSUMES that there is at least one sample per original stratum

strata_pixels = aggregate(samples[[orig_stratif]], by=list(samples[[orig_stratif]]), length)

# Calculate original strata weights and area proportions
tot_area_pix = sum(ss$pixels)
str_weight = ss$pixels / tot_area_pix

# Calculate area proportions for original strata (2001-2016). NEED to apply over MARGIN 2 (i.e. columns)
aprop = apply(ct, 2, function(x) x * str_weight / strata_pixels$x)

# Calculate accuracies
tot_acc = sum(diag(aprop)) * 100
usr_acc = diag(aprop) / rowSums(aprop) *100
prod_acc = diag(aprop) / colSums(aprop) *100

# Format and save tables
nsamp = rowSums(ct)
ct_save = cbind(ct, nsamp, str_weight)

suffix = paste0("_step", step, "_", lut_name, ".csv")
write.csv(ct_save, file=paste0(savepath, "confusion_matrix_counts_and_weights", suffix))
write.csv(aprop, file=paste0(savepath, "confusion_matrix_area_prop", suffix))
write.csv(cbind(usr_acc, prod_acc, tot_acc), file=paste0(savepath, "accuracies", suffix))


##############################################################################################################
# 3) CALCULATE REFERENCE CLASS AREA PROPORTIONS AND VARIANCE OF REFERENCE SAMPLES PER STRATA 
# 4) CALCULATE UNBIASED STANDARD ERROR FOR PROPORTION OF REFERENCE CLASS AREAS

# Initialize variables for area proportions, variances, etc (Step 1)
area_prop = matrix(0, nrow=length(years)-1, ncol=length(ref_codes), dimnames=list(years[2:length(years)], ref_codes))
ref_var_list = list()
filtered_ss = list()
ref_prop_list = list()

# Initialize empty matrix to store standard error proportions (Step 2)
se_prop = matrix(0, nrow=length(years)-1, ncol=length(ref_codes), dimnames=list(years[2:length(years)], ref_codes))

# Initialize matrices for area calculations (Step 3)
area_ha = matrix(0, nrow=length(years)-1, ncol=length(ref_codes), dimnames=list(years[2:length(years)], ref_codes))
area_ci = matrix(0, nrow=length(years)-1, ncol=length(ref_codes), dimnames=list(years[2:length(years)], ref_codes))
area_upper = matrix(0, nrow=length(years)-1, ncol=length(ref_codes), dimnames=list(years[2:length(years)], ref_codes))
area_lower = matrix(0, nrow=length(years)-1, ncol=length(ref_codes), dimnames=list(years[2:length(years)], ref_codes))
margin_error = matrix(0, nrow=length(years)-1, ncol=length(ref_codes), dimnames=list(years[2:length(years)], ref_codes))

#' Run all the calculations. Call the function for every reference year we want and get area proportions
#' and sample variance, then standard errors on proportions, then areas and their confidence intervals and
#' margin of errors. NOTE the double square brackets to allow for substitution. Only requires the original 
#' strata, so no need to iterate over yearly rasters, that's only needed for accuracies. 

for (y in (1:(length(years)-1))){
  # Compare year strata with year reference. Field names MUST start at 2002, hence i+1. Other indexed variables DO need to
  # start at 1 though.
  prop_out = calc_area_prop(samples[[orig_stratif]], samples[[ref_names[y+1]]], ss, strata_pixels, ref_codes) 
  
  # Assign outputs of Step 1
  area_prop[y,] = prop_out[[1]]
  ref_var_list[[y]] = prop_out[[2]]
  filtered_ss[[y]] = prop_out[[3]]
  ref_prop_list[[y]] = prop_out[[4]]
  tot_area_pix = prop_out[[5]] # Will be overwritten with the same value anyway...
  
  # Run step two
  se_prop[y,] = calc_se_prop(ss, strata_pixels, ref_var_list[[y]], ref_codes, tot_area_pix)
  
  # Run and assign outputs of Step 3
  areas_out = calc_unbiased_area(tot_area_pix, area_prop[y,], se_prop[y,]) 
  area_ha[y,] = areas_out[[1]]
  area_ci[y,] = areas_out[[2]]
  area_upper[y,] = areas_out[[3]]
  area_lower[y,] = areas_out[[4]]
  margin_error[y,] = areas_out[[5]]
  
}

# Calculate total area per stratum, and a filtered version with the reference
# codes only
stratum_areas= ss$stratum *30^2 / 100^2
filtered_stratum_areas= ss$stratum[ss$stratum %in% ref_codes] *30^2 / 100^2
map_bias = filtered_stratum_areas - area_ha['2016',] # doesn't work in all cases bc 
# 2016 is not always created if we aggregate the years.

# Write results to csv 
suffix = paste0("_step", step, "_", lut_name, ".csv")

write.csv(area_ha, file=paste0(savepath, "area_ha", suffix))
write.csv(area_lower, file=paste0(savepath, "area_lower", suffix))
write.csv(area_upper, file=paste0(savepath, "area_upper", suffix))

# Calculate areas and margin of errors for ORIGINAL STRATIFICATION (2001-2016) and save
# Using the same equations givee quivalent results to doing it manually via the
# confusion matrix.

prop_out_orig = calc_area_prop(samples[[orig_stratif]], strata, ss, strata_pixels, ref_codes) 
se_prop_orig = calc_se_prop(ss, strata_pixels, prop_out_orig[[2]], ref_codes, tot_area_pix)
areas_out_orig = calc_unbiased_area(tot_area_pix, prop_out_orig[[1]], se_prop_orig)

areas_orig = data.frame(t(sapply(areas_out_orig,c)))
colnames(areas_orig) = ref_codes
rownames(areas_orig) = c("area_ha", "area_ci", "area_upper", "area_lower", "margin_error")
write.csv(areas_orig,  file=paste0(savepath, "area_ha_orig_strata.csv"))


##############################################################################################################
# 5) CREATE SOME USEFUL TABLES

# Export table of  sample count, areas, percentages

stratum_percentages=round(ss$pixels / tot_area_pix * 100, digits=3) 
strata_table = as.data.frame(cbind(stratum_areas, stratum_percentages, strata_pixels$x))
strata_table$stratum_areas = format(strata_table$stratum_areas, scientific = FALSE, big.mark = ",")
rownames(strata_table) = orig_strata_names 
# Need to escape special characters, including backslash itself (e.g. $\\alpha$)
colnames(strata_table) = c("Area [ha]", "Area / $W_h$ [\\%]", "Sample size ($n_h$)") 
# Create table in Latex instead, and produce the pdf there, much easier than grid.table
print(xtable(strata_table, digits=c(0,2,2,0)),type = "latex",sanitize.text.function=function(x){x})

# Calculate map bias and create accuracy table with margin of error, only for strata 2001-2016
# Also print to latex to be converted to pdf

# TODO: Fix this block, not working at all
accuracies = as.data.frame(cbind(usr_acc, prod_acc))
accuracy_table=cbind(accuracies[-11,], t(map_bias), t(margin_error['2016',]*100))
colnames(accuracy_table) = c("User's accuracy", "Producer's accuracy", "Map bias", "Margin of error")
rownames(accuracy_table) = orig_strata_names[-11]
accuracy_table = format(accuracy_table, scientific = FALSE, big.mark = ",", digits=2)
print(xtable(accuracy_table, digits=c(0,2,2,0,1)),type = "latex",sanitize.text.function=function(x){x})

# Calculate reference sample count per year. Initialize zero matrix with year * class dims and proper row and column names
ref_sample_count = matrix(0, nrow=length(years), ncol=length(ref_codes), dimnames=list(years, ref_codes))

for (f in 1:(length(ref_names))){
  # Get table, then check if unique classes is in that table, then use the boolean to assign values
  a = table(samples[[ref_names[f]]])
  class_check = ref_codes %in% names(a)
  ref_sample_count[f,class_check] = a
}

# Format and show/save
grid.newpage()
tt =  ttheme_default(base_size=20)
grid.table(round(ref_sample_count), theme=tt) 
png(paste0(savepath, "numchange_strata_ref.png"), width=1000, height = 1000, units = "px"); grid.table(ref_sample_count, theme=tt); dev.off()

## BREAK ANALYSIS
## Calculate when break occurred in maps and compare to break in reference samples

# Function to return the indices of the locations where a there is a change in 
# class code of a vector. Useful for detecting years of change of reference and
# map labels. 

break_calc = function(values){
  out = vector()
  unq = unique(as.numeric(values))
  # When there are no changes
  if (length(unq) == 1) {
    out = 0
  # When there is one or more changes
  } else {
    out = c(1 + which(diff(as.numeric(values))!=0))
  }
  return(out)
}
    
        
# We need to load the original maps in order to avoid having to reclassify
# the strata rasters. 

alt_samples = samples
rast_names = character()
for (y in 1:length(years)){
 rast_names[y] = paste0(years[y],"_final_crop")
 map = raster(paste0(auxpath, rast_names[y], ".tif"))
 alt_samples = extract(map, alt_samples, sp = TRUE) 
}

# After reading, prepend an X to be able to call columns by name
rast_names = paste0("X", rast_names)

# Get break indices (column) for REFERENCE data in ORIGINAL class codes (e.g no strata codes)
bi_ref = apply(df2, MARGIN = 1, FUN = break_calc)
# Get the max number of elements in a subelement of the list
max_l <- max(sapply(bi_ref, length))
# Fill with NA's using the max number of elements
l <- lapply(bi_ref, function(v) { c(v, rep(NA, max_l - length(v)))})
# Bind into a dataframe
bi_ref = do.call(rbind, l)

# Get break indices (column) for the original MAPPED data. 

bi_map = apply(alt_samples@data[rast_names], MARGIN = 1, FUN = break_calc) 
# Get the max number of elements in a subelement of the list
max_l <- max(sapply(bi_map, length))
# Fill with NA's using the max number of elements
l <- lapply(bi_map, function(v) { c(v, rep(NA, max_l - length(v)))})
# Bind into a dataframe
bi_map = do.call(rbind, l)
# Substract 1 and replace -1 with zero. We need this beause the maps really represent the 
# dynamic of the year before them. This also allows us to index from STRATA maps directly
# (because the max == 15)
bi_map = bi_map - 1
bi_map[bi_map == -1] = 0

# Get rows where there is only one or no changes in the REFERENCE AND MAP labels
# to simplify the analysis.
# Use the NUMCHANGES column because it is a better source of break info

ref_chg = samples@data["NUMCHANGES"] <= 1 
numchg_map = apply(bi_map, MARGIN = 1, FUN = function(x) length(unique(na.omit(x))))
map_chg = numchg_map <= 1
subind = which(ref_chg & map_chg)
subind2 = which(!(ref_chg & map_chg)) # Get the complementary indices, we will need them

# Subset reference and map labels based on those indices
# Remove all the NA's (effectively making it a vector)
bi_ref_sub = bi_ref[subind, ]
bi_ref_sub = bi_ref_sub[!is.na(bi_ref_sub)]

bi_map_sub = bi_map[subind, ]
bi_map_sub = bi_map_sub[!is.na(bi_map_sub)]

# Function to get the strata codes using those indices, both for reference and maps.
# There must be a way to vectorize that but I don't see it
get_break_strata = function(rw, index){
  # When there is no change for reference (=0)
  if(index <= 0){
    code = rw[,1]
  } else {
    code = rw[,index]  
  }
  return(code)
}

# Initialize vectors and get reference and map labels using the column indices
break_ref = vector(length = length(subind))
break_map = vector(length = length(subind))

for(i in 1:length(subind)){
  break_ref[i] = get_break_strata(df[subind[i],], bi_ref_sub[i])  
  break_map[i] = get_break_strata(samples@data[subind[i], map_names], bi_map_sub[i])  
}

# Create dataframe with ref breaks map breaks and their corresponding labels.
# Create a larger df with these results and the row id so that we can go back to the 
# original samples easily. 
break_ref_df = cbind(bi_ref_sub, break_ref)
break_map_df = cbind(bi_map_sub, break_map)
break_compare = as.data.frame(cbind(subind, break_ref_df, break_map_df, samples@data[subind, "PTRW"]))
colnames(break_compare) = c("row_id", "ref_year", "ref_label", "map_year", "map_label", "ptrw")

# Convert the indices to actual years to make subsequent analysis easier
# (Just add 2000 and replace 2000 with zeros!)
break_compare$ref_year = break_compare$ref_year + 2000
break_compare$ref_year[break_compare$ref_year == 2000] = 0

break_compare$map_year = break_compare$map_year + 2000
break_compare$map_year[break_compare$map_year == 2000] = 0

# Create confusion matrix of breaks and labels between map and reference
cm_breaks = calc_ct(break_compare$map_year, break_compare$ref_year)
cm_labels = calc_ct(break_compare$map_label, break_compare$ref_label)

# Add row and column totals and display/save
cm_breaks = cbind(cm_breaks, total=rowSums(cm_breaks))
cm_breaks = rbind(cm_breaks, total=colSums(cm_breaks)) 

# Format table and display
tt= ttheme_default(base_size=18)
grid.newpage()
grid.table(round(cm_breaks, digits=2), theme=tt)
png("numchange_strata_ref.png", width=1000, height = 1000, units = "px"); grid.table(cm_breaks, theme=tt); dev.off()

# Find records for any given set of map year (from 1 to 16) and change year and get comparison of labels for it
# Eg show label distribution for BREAK OMISSION ERRORS. 
records = break_compare[break_compare$map_year == 0 & break_compare$ref_year != 0,]
calc_ct(records$map_label, records$ref_label)

# Eg show label distribution for BREAK COMMISSION ERRORS
records = break_compare[break_compare$map_year != 0 & break_compare$ref_year == 0,]
calc_ct(records$map_label, records$ref_label)

# Eg show label distribution for EXACT CHANGE MATCH BREAK DETECTION
records = break_compare[(break_compare$map_year == break_compare$ref_year) & (break_compare$map_year !=0 & break_compare$ref_year != 0),]
calc_ct(records$map_label, records$ref_label)

# Eg show label distribution for STABLE MATCH
records = break_compare[break_compare$map_year == 0 & break_compare$ref_year == 0,]
calc_ct(records$map_label, records$ref_label)

# Can also be done for labels, to find years of break
records = break_compare[break_compare$map_label == 8 & break_compare$ref_label == 8,]
calc_ct(records$map_year, records$ref_year)

# Analyze strata breaks vs ref breaks per path/row
break_compare$breakdif = break_compare$ref_year - break_compare$map_year

# Calculate type of break difference per sample: ommited change, commited change, 
# detected change with small ofset in timing (either positive or negative), 
# exact change match (first change detected in the same year) or stable match
# (no changes detected in reference or maps)

break_compare$type[break_compare$breakdif > 2000] = "omission"
break_compare$type[break_compare$breakdif < -2000] = "comission" 
break_compare$type[(break_compare$breakdif > 0) & (break_compare$breakdif < 2000)] = "pos_offset" 
break_compare$type[(break_compare$breakdif < 0) & (break_compare$breakdif > -2000)] = "neg_offset" 
# Assume rows == 0 are exact change matches, then overwrite those that are stable matches
break_compare$type[break_compare$breakdif == 0] = "exact_change_match" 
break_compare$type[(break_compare$ref_year == 0) & (break_compare$map_year == 0)] = "stable_match"

# Create new dataframe with the original number of rows and populate with the
# results. Fill in the rows with no results (those that we didnt analyze)
long_break_compare = data.frame(matrix(ncol = ncol(break_compare), nrow = nrow(df)))
colnames(long_break_compare) = colnames(break_compare)
long_break_compare[subind, ] = break_compare
long_break_compare[subind2, "type"] = "not_analyzed"
long_break_compare[subind2, "ptrw"] = samples@data[subind2,"PTRW"]

# Aggregate results and calculate total points per path-row. 
breakdif_count = aggregate(long_break_compare$type, by=list(long_break_compare$ptrw), FUN="summary")
total_ptrw_pts = apply(breakdif_count[,2], MARGIN = 1, FUN = "sum")

# # Calculate number of each type of change as ratios of the total
bd_ratios = t(apply(breakdif_count[,2], MARGIN = 1, function(x) x/sum(x))) 
rownames(bd_ratios) = breakdif_count$Group.1
bd_ratios_df = cbind(bd_ratios, total_ptrw_pts)

# Write aggregated results to CSV and row results to shapefile
write.csv(bd_ratios_df, paste0(savepath, "LC_change_ratios_pathrow.csv"))
samples@data[,names(long_break_compare)] <- long_break_compare
writeOGR(samples, paste0(savepath, "sample_yearly_strata_break_analysis"), "sample_yearly_strata_break_analysis", driver="ESRI Shapefile", overwrite_layer = T)


# Create confusion matrices per year between reference and map labels. USING FULL DATA HERE
cm_list = list()
for (y in 1:(length(years) - 1)){
  cm_list[[y]] = calc_ct(samples[[map_names[y]]], samples[[ref_names[y]]], class_codes)
}


# Compare total breaks detected per year
total_breaks = cbind(cm_breaks["total",], cm_breaks[,"total"])
colnames(total_breaks) = c("ref_breaks", "map_breaks")


##############################################################################################################
# 6) CREATE SOME USEFUL PLOTS
# Most of these plots were recreated in python (add name of the script) bc the lack of dual axis support in ggplot. This code is left for 
# reference.

## Plot area estimates with CI and margin of error in separate plot. Couldn't figure out how to do it with facets
# so I did it with grids. 

for(i in 1:length(ref_codes)){
  
  # Define data and variables to use
  tempdf <- as.data.frame(cbind(years[2:length(years)], area_ha[,i], area_lower[,i], area_upper[,i], margin_error[,i]))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Margin_error")
  
  # Plot areas with CI
  # Specify data, add "ribbon" with lower and higher CI and fill, then plot the estimated area with a line.
  # Other custom settings for number of breaks, label formatting and title.
  
  a <- ggplot(data=tempdf, aes(x=Years, y=Area_ha)) + 
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3) + geom_line() + 
    scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL) + 
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
    ylab("Area and 95% CI [ha]") + ggtitle(strata_names[[i]]) + expand_limits(y=0) +
    theme(plot.title = element_text(size=19), axis.title=element_text(size=16), axis.text=element_text(size=16))
  
  
  # Plot margin of error
  b <- ggplot(data=tempdf, aes(x=Years, y=Margin_error * 100)) + geom_line(size=1.1) + 
    scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL) + expand_limits(y=0) + ylab("Margin of error [%]") +
    theme(axis.title=element_text(size=16), axis.text=element_text(size=16))
  
  # To remove grid and background add this:
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()
  
  # Use gtable to stack plots together with matching extent and save  
  g1 <- ggplotGrob(a)
  g2 <- ggplotGrob(b)
  g <- rbind(g1, g2, size="first") # stack the two plots
  g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
  
  grid.newpage()
  grid.draw(g)
  filename = paste0(savepath, strata_names[[i]], "_areas_me_step", step, "_", lut_name, ".png")
  png(filename, width=1000, height = 1000, units = "px"); plot(g); dev.off()

}


## Plot net regrowth (net secondary forest). Yearly loss in 5 equals yearly gain in 14

net_rg = area_ha[,6] + area_ha[,10] + area_ha[,9] - area_ha[,11]

# Get only gain classes 
regr_area = area_ha[,c(6,9,10)]
names(regr_area) = c("Stable secondary forest", "Forest to secondary forest", "Other gains of secondary forest") #, "Loss of regrowth")

# Melt and plot
regr_area_melt = melt(as.matrix(regr_area))

regr_plot <- ggplot(regr_area_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
  geom_area(position="stack", alpha=0.8) + #geom_line(aes(x=Var1, y=732031.3), linetype=2) + 
  scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Area [ha]") + xlab("Years")+
  scale_fill_brewer(palette="GnBu",  breaks=levels(regr_area_melt$Var2), guide = guide_legend(reverse=T)) + 
  theme(axis.title=element_text(size=15), axis.text=element_text(size=13), legend.text=element_text(size=13),
  legend.title=element_blank()) #+ 
  #annotate("segment", x = 2016, xend=2016, y = 732031.3, yend = 416495.7, colour = "black", linetype=2) +
  #annotate("text", x = 2014,  y = 600000, label="Loss of regrowth")  

print(regr_plot)
#ggsave("net_regrowth_secondaryforest.png", plot=regr_plot, device="png") 


## Plot forest converstion to pasture and secondary forest 
# Forest loss = total area of classs 8 + 9 + whatever is left
for_loss = area_ha[,c(8,9)]
names(for_loss) = c("Forest to pasture", "Forest to secondary forest")

# Melt and plot. 
for_loss_melt = melt(as.matrix(for_loss))
forest_loss_plot <- ggplot(for_loss_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
  geom_area(position="stack", alpha=0.8) + 
  scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Loss of primary forest [ha]") + xlab("Years") +
  scale_fill_brewer(palette="GnBu", breaks=levels(for_loss_melt$Var2), guide = guide_legend(reverse=T)) + 
  theme(legend.title=element_blank()) +
  theme(axis.title=element_text(size=15), axis.text=element_text(size=13), legend.text=element_text(size=13))
  
print(forest_loss_plot)
#ggsave("forest_loss.png", plot=forest_loss_plot, device="png") 


# COMPARE MAPPED AREAS WITH ESTIMATES

# Create empty matrix to store mapped values and fill
mapped_areas = matrix(0, nrow=length(years[2:length(years)]), ncol=nrow(mapped_areas_list[[1]]), byrow=T, dimnames = list(years[2:length(years)], mapped_areas_list[[1]][,1]))

for (i in 1:length(mapped_areas_list)){
  mapped_areas[i,] = mapped_areas_list[[i]][,2]  
}


# Convert to ha
mapped_areas = mapped_areas * 30^2 / 100^2

# Calculate yearly area change and rate change (percentage of total area that is changing)
# Initialize zero matrix with year (03 to 16) * class (16) dimensions and proper row and column names
mapped_chg = matrix(0, nrow=length(years)-2, ncol=ncol(mapped_areas), dimnames=list(years[3:16], colnames(mapped_areas)))
mapped_rate = matrix(0, nrow=length(years)-2, ncol=ncol(mapped_areas), dimnames=list(years[3:16], colnames(mapped_areas)))

# Iterate over strata labels
for(i in 1:ncol(mapped_areas)){
  mapped_chg[,i] = diff(mapped_areas[,i])
  mapped_rate[,i] = mapped_chg[,i] / mapped_areas[1:14,i] * 100 #14 years
}


mapped_for_change = mapped_chg[,c(2, 9, 10)] # This matrix has extra cols

# Calculate how much of deforestation is not caused by pastures or regrowth and get absolute values
mfor_to_other = rowSums(mapped_for_change[,1:3])
mapped_for_change = cbind(mapped_for_change, mfor_to_other)
mapped_for_change = abs(mapped_for_change[,2:4])
colnames(mapped_for_change) = c("To pasture", "To regrowth", "To other")
mapped_for_chg = mapped_chg[,c(8,9)] 

## Melt and plot annual areas of change per year 
mfor_change_melt = melt(mapped_for_change)
mforest_plot <- ggplot(mfor_change_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
  geom_area(position="stack") + 
  scale_x_continuous(breaks=years[3:16], minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Annual forest conversion to other classes [ha]") + xlab("Years")+
  scale_fill_brewer(palette="GnBu", breaks=levels(mfor_change_melt$Var2), guide = guide_legend(reverse=T)) + 
  theme(legend.title=element_blank()) +
  theme(axis.title=element_text(size=15), axis.text=element_text(size=13), legend.text=element_text(size=13))

print(mforest_plot)
#ggsave("mapped_forest_change.png", plot=mforest_plot, device="png") 

## Plot MAPPED areas of forest to pasture, forest to regrowth and deforestation
mapped_plots = as.data.frame(mapped_areas[,c(9, 10)])
mapped_plots = cbind(years[2:length(years)], mapped_plots)
colnames(mapped_plots) = c("Years", "Forest_to_pasture",  "Forest_to_regrowth")

mfp <- ggplot(data=mapped_plots, aes(x=Years, y=Forest_to_pasture)) + geom_line(size=1.1) + 
  scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL) + expand_limits(y=0) + ylab("Area [ha]") +
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=16))

mfr <- ggplot(data=mapped_plots, aes(x=Years, y=Forest_to_regrowth)) + geom_line(size=1.1) + 
  scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL) + expand_limits(y=0) + ylab("Area [ha]") +
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=16))

mdefor <- ggplot(data=mapped_plots, aes(x=Years, y=Forest_to_pasture + Forest_to_regrowth)) + geom_line(size=1.1) + 
  scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL) + expand_limits(y=0) + ylab("Area [ha]") +
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=16))


print(mdefor)
#ggsave("mapped_forest_to_pasture.png", plot=mfp, device="png") 
#ggsave("mapped_forest_to_secondary_forest.png", plot=mfr, device="png") 
#ggsave("mapped_deforestation.png", plot=mdefor, device="png") 


## Plot estimated areas with CI along with mapped areas. The way it is right now is only taking the mapped area of forest to pasture,
# so the rest of the plots are not correct. I need to fix that.

for(i in 1:length(area_ha)){
  
  # Define data and variables to use
  tempdf <- as.data.frame(cbind(years[2:length(years)], area_ha[,i], area_lower[,i], area_upper[,i], mapped_areas[,'8']))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Mapped")
  
  # Specify data, add "ribbon" with lower and higher CI and fill, then plot the estimated area with a line.
  # Other custom settings for number of breaks, label formatting and title.
  a<-ggplot(data=tempdf, aes(x=Years, y=Area_ha))+ geom_line(aes(x=Years, y=Mapped),linetype=4 )+
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3) + geom_line() + 
    scale_x_continuous(breaks=years[2:length(years)], minor_breaks = NULL) + 
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
    ylab("Area in ha") + ggtitle(strata_names[[i]]) + expand_limits(y=0)
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=13))
  filename = paste0(strata_names[[i]], "_areasCI", ".png")
  
  print(a)
  #ggsave(filename, plot=a, device="png") 
  
}

##############################################################################################################
# MISCELANEOUS
# Plot number of forest to pasture reference samples over time.
fpc = vector()
for (f in ref_names){
  a = samples@data[,f] == 8
  fpc = cbind(fpc, sum(a))
}

dtf = as.data.frame(cbind(years, t(fpc)))
colnames(dtf) = c("Years", "Count")


