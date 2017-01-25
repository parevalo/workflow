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
wd = "C:/OneDrive/Lab/sample_may2016/interpreted_w_strata_23062016"
auxpath = "C:/test"

if( .Platform$OS.type == "unix" )
    wd = "/media/paulo/785044BD504483BA/OneDrive/Lab/sample_may2016/interpreted_w_strata_23062016"
    auxpath = "/media/paulo/785044BD504483BA/test/"
    funcpath = "/home/paulo/workflow/multi_scene/7_poststratification/strata_calculation.R"

setwd(wd)
source(funcpath)

# If TRUE, uses class 8 and 9 together as deforestation, labeled as class 17. Otherwise keeps them separate
deformode = FALSE

# List of original strata names
orig_strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                      "Stable pasture-cropland", "Stable secondary forest", "Stable water", "Forest to pasture", 
                      "Forest to secondary forest", "Gain of secondary forest", "All to unclassified", "Loss of secondary forest")


# Strata names depending on deformode. They DON'T have the "all to unclass." class bc it disappears when the ref. samples are collected.
if (deformode == FALSE){
  strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                   "Stable pasture-cropland", "Stable secondary forest", "Stable water", "Forest to pasture", 
                   "Forest to secondary forest", "Gain of secondary forest", "Loss of secondary forest")
}else{
  strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                   "Stable pasture-cropland", "Stable regrowth", "Stable water", "All others to regrowth", 
                   "Loss of regrowth", "Deforestation")
}

# Grid table theme, only used to display some ancillary tables
tt=ttheme_default(core=list(fg_params=list(font="Times", fontface="plain", fontsize=14)),
                  colhead=list(fg_params=list(font="Times", fontface="bold", fontsize=14, parse=TRUE)),
                  rowhead=list(fg_params=list(font="Times", fontface="plain", fontsize=14)))

#############################################################################################################
#1) READ SHAPEFILE AND CALCULATE REFERENCE LABELS PER YEAR
# This section takes the sample shapefile and assigns the proper strata for each year, depending
# on when the model break happened. In this version we use all of the points because when done 
# annually, it does not matter how many changes we found in the TS. However, given that in the shapefile
# the strata was calculated only betwen first and last label, we need to create as many strata as
# NUMCHANGES. For this reason a new R function will be introduced to do this. 

# Read shapefile with reference strata and change info. Test if working from ubuntu or windows
samples <- readOGR(".", "final_extended_sample_merge_UTM18N_point")

# We need to convert the dates to characters bc they are factors right now
samples$CHGDATE <- as.character(samples$CHGDATE)

# Set up important loop variables
rows = nrow(samples)
start = 2001
end = 2016
years = seq(start, end)
field_names = paste("ref_", years, sep="")

# Initialize empty vectors to store the data, current change and class code (start with the first!) 
# and name of the fields that store the class codes.

df = vector()
field = vector()
codelist = c("CODE1", "CODE2", "CODE3", "CODE4")
current_change = 1
current_code=1

# Iterate over rows (easier to do calculations by row)
for (row in 1:rows){
  current_change = 1
  current_code=1
  # Iterate over years
  for (i in 1:length(years)){
    
      # If there are no changes, just use the only class code for both years
      if (is.na(samples$CHGDATE[row]) == TRUE) {
        field[i] = calculate_strata(samples$CODE1[row], samples$CODE1[row])
      }
    
      # If there is a change, compare each year to the current change year and update that one accordingly
      else {
        # Get the current date of change (initialized as the first change available)
        chg_date = unlist(strsplit(as.vector(samples$CHGDATE[row]), ','))[current_change] #742
        # Extract the YEAR from the date of change
        chg_year = na.omit(as.numeric(unlist(strsplit(chg_date, "[^0-9]"))))[1]
        
        # If we haven't reached change year yet
        if (years[i] < chg_year) {
          field[i] = calculate_strata(samples@data[codelist[current_code]][row,], samples@data[codelist[current_code]][row,]) 
          #print(paste0(row, " cond1 ", years[i], " ", chg_year," ", field[i]))
        } 
        
        # If we JUST reached a change year
        else if (years[i] == chg_year) {
          field[i] = calculate_strata(samples@data[codelist[current_code]][row,], samples@data[codelist[current_code+1]][row,]) 
          # Check if we haven't reached the max number of recorded changes
          if (current_change < samples$NUMCHANGES[row]) {
            current_change = current_change + 1 
          }
          current_code = current_code + 1
          #print(paste0(row, " cond2 ",years[i], " ", chg_year," ", field[i]))
        }
        
        # If we went past the last change date, use the last code 
        else if (years[i] > chg_year) {
          field[i] = calculate_strata(samples@data[codelist[samples$endcodecol[row]]][row,], samples@data[codelist[samples$endcodecol[row]]][row,])
          #print(paste0(row, " cond3 ",years[i], " ", chg_year," ", field[i]))
        }
      } 
  }
  df <- as.data.frame(rbind(df, field))
}

names(df) = field_names

# Attach table to shapefile 
samples@data[,field_names] <- df
#writeOGR(samples, "sample_yearly_strata", "sample_yearly_strata", driver="ESRI Shapefile", overwrite_layer = T)

#############################################################################################################
# 2) READ ORIGINAL AND ANNUAL STRATA RASTERS AND EXTRACT THEIR VALUES TO THE SHAPEFILE
# Also calculate original strata size, weights, area proportions, and accuracies

# Create a list with the annual raster names
years_short = seq(02,16)

# Iterate over names and extract to shapefile. We only need these to calculate accuracies per year.
for (y in years_short){
  rast_name = paste0("final_strata_annual_", sprintf("%02d",y-1), "_", sprintf("%02d",y), "_UTM18N")
  map = raster(paste0(auxpath, rast_name, ".tif"))
  samples = extract(map, samples, sp = TRUE) 
}

# Read original stratification. Done last so that the reference and map fields are contiguous
samples = extract(raster(paste0(auxpath, "/final_strata_01_16_UTM18N.tif")), samples, sp=TRUE)


#################
# Artificially modify reference data. THIS SECTION IS EXPERIMENTAL. 
# Results follow the same behavior found in the cumulative version of this script (e.g. a much higher number of correct samples
# of stable forest and forest to pasture decrease the widht of CI and make it bigger than zero, reduces the spikes and the margin of
# error, but that one is still higher than one might expect). Higher number of stable forest reduces CI, higher number of correct
# change samples reduces margin of error. For example, 1000 more stable forest sample and only 48 more samples of correct forest to
# pastures make the CI get above zero and its total width around 200K. Increasing that forest sample to 2000 reduces the total CI
# to about 125K.

# Create 16 ref columns and 16 map columns, with a given number of rows. Determine proportion of samples that will be
# right or wrong, and with which classes.
#samples@data = rbind(samples@data, samples@data)

# create backup
#samples_backup = samples
samples = samples_backup

# Add many correct forest samples. EXTENDED bc we need to include original stratification!
df <- data.frame(matrix(1,ncol = 32, nrow = 2000))
colnames(df) = colnames(samples@data)[23:54]

# Create matrix with 1, 4 and diag = 8
repet = 10
mat = matrix(1,  ncol=16, nrow=16)
mat[upper.tri(mat)] = 4
diag(mat) = 8
mat = matrix(rep(t(mat), repet) , ncol=ncol(mat) , byrow=TRUE)
df5 = cbind(mat, mat[,2:16], rep(8, 16*repet)) # We need to fill the original strata column
colnames(df5) = colnames(samples@data)[23:54]

# Same but shifting map data 1 year. DOESNT HAVE ANY EFFECT SINCE INDIVIDUAL MAPS ARE NOT USED IN CALCULATIONS
# EXCEPT IF WE WANT TO CALCULATE ANNUAL ACCURACIES....
#mat_shift = shift.right(mat, 1, 1)
#df5 = cbind(mat, mat_shift[,1:15], rep(8, 16*repet)) # We need to fill the original strata column
#colnames(df5) = colnames(samples@data)[23:54]

# Create matrix with 1, 5 and diag = 9
mat = matrix(1,  ncol=16, nrow=16)
mat[upper.tri(mat)] = 5
diag(mat) = 9
mat = matrix(rep(t(mat), repet) , ncol=ncol(mat) , byrow=TRUE)
df6 = cbind(mat, mat[,2:16], rep(8, 16*repet)) # We need to fill the original strata column
colnames(df6) = colnames(samples@data)[23:54]

samples@data = rbind(samples@data[,23:54], df)
samples@data = rbind(samples@data[,23:54], df, df5, df6)

# Create artificial sample with perfect reference/map matching and same number of total samples as the real data
# Shows same behavior than in cummulative version of this script. Even with perfect sample, a class like secondary
# forest has a CI that goes below zero. Why?? Maybe too few samples?

sample_size = strata_pixels$x # Relies on original strata pixels, fix!
sample_size[2] = 400
sample_size[8] = 80
sample_size[9] = 50
ref_matrix = matrix(, ncol=17, nrow=0) # Initialize empty matrix
for (i in 1:dim(strata_pixels)[1]){
  code = strata_pixels$Group.1[i]
  repet = ceiling(sample_size[i]/16)
  if (code <= 6){
    mat = matrix(code, ncol=16, nrow=16)
  } else if (code == 8 ) {
    mat = matrix(1, ncol=16, nrow=16)
    mat[upper.tri(mat)] = 4
    diag(mat) = code
  } else if (code == 9 ) {
    mat = matrix(1, ncol=16, nrow=16)
    mat[upper.tri(mat)] = 5
    diag(mat) = code
  } else if (code == 11) {
    mat = matrix(4, ncol=16, nrow=16)
    mat[upper.tri(mat)] = 5
    diag(mat) = code
  } else if (code == 13) {
    next
  } else if (code == 14) {
    mat = matrix(5, ncol=16, nrow=16)
    mat[upper.tri(mat)] = 4 # Lets assume it all changes to pastures 
    diag(mat) = code
  }
  mat = cbind(mat, rep(code, 16))
  print(ncol(mat))
  mat2 = matrix(rep(t(mat), repet) , ncol=ncol(mat), byrow=TRUE)
  ref_matrix = rbind(ref_matrix, mat2)
}

full_matrix = cbind(ref_matrix[,1:16], ref_matrix[,2:16], ref_matrix[,17]) # If we dont want to shift the samples
colnames(full_matrix) = names(samples_backup[23:54])
samples@data = as.data.frame(full_matrix)

map_shifted = shift.right(ref_matrix, 1, 1) # If we want to shift the samples
full_matrix = cbind(ref_matrix, map_shifted[,1:15])

###############

# Crosstab final strata and reference strata (for the same period 01-16) 
ct = table(samples$final_strata_01_16_UTM18N, samples$strata)

# Load mapped areas (total strata sample size) produced from CountValues.py, bc calculating it here with hist() takes forever..
# TODO: CALCULATE NEW ONES FOR THE NEW MAPS AND CHANGE THE CODE ACCORDINGLY.
mapped_areas_list = list()
filenames = dir(auxpath, pattern="*_pixcount.csv")
for(i in 1:length(filenames)){
  mapped_areas_list[[i]] = read.csv(paste0(auxpath,filenames[i]), header=TRUE, col.names=c("stratum", "pixels"))
}

# Get only the total area for the original strata (i.e 01-16)
ss = mapped_areas_list[[15]]
# Classes to be removed/ignored
cr = c(7, 10, 12, 15)
# Filter classes NOT in that list
ss = ss[!(ss$stratum %in% cr),] 

# Calculate total number of samples per ORIGINAL stratum
strata_pixels = aggregate(samples$final_strata_01_16_UTM18N, by=list(samples$final_strata_01_16_UTM18N), length)

# Calculate original strata weights and area proportions
tot_area_pix = sum(ss$pixels)
str_weight = ss$pixels / tot_area_pix
# Add column for class 13 with zeroes to obtain a square matrix and make everything easier
cmatr = cbind(ct[,1:10], matrix(0, nrow=nrow(ct), ncol=1, byrow=T, dimnames = list(rownames(ct), "13")), ct[,11,drop=F])
# Calculate area proportions for original strata (2001-2016). NEED to apply over MARGIN 2 (e.g. columns)
aprop = apply(cmatr, 2, function(x) x * str_weight / strata_pixels$x)

# Calculate accuracies
tot_acc = sum(diag(aprop)) * 100
usr_acc = diag(aprop) / rowSums(aprop)
prod_acc = diag(aprop) / colSums(aprop)

##############################################################################################################
# 3) CALCULATE REFERENCE CLASS AREA PROPORTIONS AND VARIANCE OF REFERENCE SAMPLES PER STRATA 

# Takes vectors of sample strata, sample references, their totals and the reference codes (strata codes are calculated) 
# Returns a list of vectors with length equal to the number of reference classes. The list contains: 
# Area proportions per class, sample variance per sample strata, total strata size present in that year, 
# proportion of sample ref class present on each samp strata
calc_area_prop = function(samp_strata, samp_reference, strata_totals, sample_totals, rfcodes){

  # Obtain unique values in the samp_strata field
  str_codes = sort(unique(samp_strata))
  str_ind = seq_along(str_codes)
  
  # Get a sequence for the reference codes
  ref_ind = seq_along(rfcodes)
  
  # Initialize empty df for proportions per samp_strata, and for sample variance per samp_strata
  ref_prop = data.frame()
  ref_var = data.frame()
  
  # Compare the fields, iterate over "samp_strata" and "samp_reference" classes
  for (s in str_ind){
    for (r in ref_ind){
      # Compared fields and get TRUE or FALSE on each row
      cond_bool = samp_strata == str_codes[s] & samp_reference == rfcodes[r]
      # Get row numbers that meet that condition
      ind = which(cond_bool)
      # Get location of rows for the current samp_strata
      str_ref = which(samp_strata == str_codes[s])
      # Get proportion of samp_reference class present on that samp_strata
      ref_prop[s, r] = length(ind)/sample_totals$x[sample_totals$Group.1 == str_codes[s]]
      # Calculate SAMPLE variance of current samp_reference code in current samp_strata (needed later)
      # which is the same formula specified in the paper
      ref_var[s, r] = var(cond_bool[str_ref]) 
      
    }
  }
  
  # Assign column and row names for easier samp_reference
  rownames(ref_prop) = paste0("strat_",str_codes)
  colnames(ref_prop) = paste0("ref_", rfcodes)
  rownames(ref_var) = paste0("strat_",str_codes)
  colnames(ref_var) = paste0("ref_", rfcodes)
  
  # Calculate samp_reference class proportions (i.e. by columns) using total, original samp_strata areas.
  class_prop = vector()
  # Filter only total sample sizes that are present in the samp_strata for that year
  fss = strata_totals[strata_totals$stratum %in% str_codes,]
  for (r in 1:ncol(ref_prop)){
    # LEAVE THE SUM OF THE ENTIRE samp_strata
    class_prop[r] = sum(fss$pixels * ref_prop[,r])/tot_area_pix
  }
  return(list(class_prop, ref_var, fss, ref_prop)) 
}

# IF DEFORMODE = TRUE THESE STEPS ARE REQUIRED BC WE NEED ONE CLASS INSTEAD OF TWO
# a) calculate deforestation combined (e.g 8 + 9): Reclassify samples$final and samples$ref<year>, ss, strata_pixels.
defor_str = samples$final_strata_01_16_UTM18N
defor_str[defor_str == 8 | defor_str == 9] = 17

# b) Find reference samples = 8 or 9 for each year and change value to 17
samples_defor = samples
for (f in 1:(length(field_names))){
  defor_ind = samples_defor[[field_names[f]]] == 8 | samples_defor[[field_names[f]]] == 9
  samples_defor@data[defor_ind,field_names[f]] = 17
}

# c) Sum class 8 and 9 to produce "class 17", then delete 8 and 9. Do that for total pixels and sample counts
defor_str_totals = ss 
class17a = defor_str_totals[defor_str_totals$stratum == 8,] + defor_str_totals[defor_str_totals$stratum == 9,]
defor_str_totals = rbind(defor_str_totals, class17a)
defor_str_totals = defor_str_totals[!(defor_str_totals$stratum == 8 | defor_str_totals$stratum == 9),]

defor_samp_totals = strata_pixels 
class17b = defor_samp_totals[defor_samp_totals$Group.1 == 8,] + defor_samp_totals[defor_samp_totals$Group.1 == 9,]
defor_samp_totals = rbind(defor_samp_totals, class17b)
defor_samp_totals = defor_samp_totals[!(defor_samp_totals$Group.1 == 8 | defor_samp_totals$Group.1 == 9),]

# Deformode modifies output from HERE

# Get unique classes through all the reference years (This won't have class 13 for that reason) and get numbered sequence
if (deformode == TRUE){
  ref_codes = sort(unique(unlist(samples_defor@data[field_names])))
}else { 
  ref_codes = sort(unique(unlist(samples@data[field_names])))
}


# Call the function for every reference year we want and get area proportions and sample variance
# NOTE the double square brackets to allow for substitution. Only requires the original strata, so no need to 
# iterate over yearly rasters, that's only needed for accuracies. 

area_prop = data.frame()
ref_var_list = list()
filtered_ss = list()
ref_prop_list = list()
for (i in (1:length(years_short))){
  # Compare year strata with year reference. Field names MUST start at 2002, hence i+1.
  if (deformode == FALSE){
    out = calc_area_prop(samples$final_strata_01_16_UTM18N, samples[[field_names[i+1]]], ss, strata_pixels, ref_codes) 
  } else {
    out = calc_area_prop(defor_str, samples_defor[[field_names[i+1]]], defor_str_totals, defor_samp_totals, ref_codes) 
  }
  ap = out[[1]]
  rv = out[[2]]
  fs = out[[3]]
  rp = out[[4]]
  area_prop = rbind(area_prop, ap)
  ref_var_list[[i]] = rv
  filtered_ss[[i]] = fs
  ref_prop_list[[i]] = rp
}

# Assign names to make it easier to interprete
rownames(area_prop) = years[2:length(years)]
colnames(area_prop) = ref_codes

##############################################################################################################
# 4) CALCULATE UNBIASED STANDARD ERROR FOR PROPORTION OF REFERENCE CLASS AREAS

# Function to calculate standard error of area proportion of reference classes
calc_se_prop = function(strata_totals, sample_totals){

  # Formula works correctly, tested
  se_pr = data.frame()
  
  #Iterate over years
  for (y in 1:length(ref_var_list)){
    # Iterate over classes
    for (c in 1:length(ref_codes)){
    v = 1/tot_area_pix^2 * (sum(strata_totals$pixels^2 * (1 - sample_totals$x/strata_totals$pixels) * (ref_var_list[[y]][,c] / sample_totals$x)))
    se_pr[y,c] = sqrt(v)
    }
  }
  return(se_pr)
}

# Check deforestation mode and get calculation accordingly
if (deformode == FALSE){
  se_prop = calc_se_prop(ss, strata_pixels)
  }else{
  se_prop = calc_se_prop(defor_str_totals, defor_samp_totals)
}
  
  
# Reassign names...
rownames(se_prop) = years[2:length(years)]
colnames(se_prop) = ref_codes

# Total area in ha
N_ha = tot_area_pix * 30^2 / 100^2
# Calculate area in ha from area proportions
area_ha = area_prop * N_ha
# Calculate confidence interval in ha
area_ci = se_prop * 1.96 * N_ha
#Upper and lower CI
area_upper = area_ha + area_ci
area_lower = area_ha - area_ci

# Write results to csv 
if (deformode == TRUE){
   suffix = "annual_defor.csv"  
 } else { 
   suffix = "annual_regular.csv"
}

#write.csv(area_ha, file=paste0("area_ha", suffix))
#write.csv(area_lower, file=paste0("area_lower", suffix))
#write.csv(area_upper, file=paste0("area_upper", suffix))


##############################################################################################################
# 5) CREATE SOME USEFUL TABLES

# Export table of  sample count, areas, percentages
stratum_areas= ss$pixels *30^2 / 100^2
stratum_percentages=round(ss$pixels / tot_area_pix * 100, digits=3) 
strata_table = as.data.frame(cbind(stratum_areas, stratum_percentages, strata_pixels$x))
strata_table$stratum_areas = format(strata_table$stratum_areas, scientific = FALSE, big.mark = ",")

# Create table in Latex instead, and produce the pdf there, much easier than grid.table
# Need to escape special characters, including backslash itself (e.g. $\\alpha$)
rownames(strata_table) = orig_strata_names
colnames(strata_table) = c("Area [ha]", "Area / $W_h$ [\\%]", "Sample size ($n_h$)") 
print(xtable(strata_table, digits=c(0,2,2,0)),type = "latex",sanitize.text.function=function(x){x})

# Calculate map bias and create accuracy table with margin of error, only for strata 2001-2016
# Also print to latex to be converted to pdf
margin_error = area_ci / area_ha 
map_bias = stratum_areas - area_ha['2016',]
accuracies = as.data.frame(cbind(usr_acc, prod_acc))
accuracy_table=cbind(accuracies[-11,], t(map_bias), t(margin_error['2016',]*100))
colnames(accuracy_table) = c("User's accuracy", "Producer's accuracy", "Map bias", "Margin of error")
rownames(accuracy_table) = orig_strata_names[-11]
accuracy_table = format(accuracy_table, scientific = FALSE, big.mark = ",", digits=2)
print(xtable(accuracy_table, digits=c(0,2,2,0,1)),type = "latex",sanitize.text.function=function(x){x})

# Calculate reference sample count per year. Initialize zero matrix with year * class dims and proper row and column names
ref_sample_count = matrix(0, nrow=length(years), ncol=length(ref_codes), dimnames=list(years, ref_codes))

for (f in 1:(length(field_names))){
  # Get table, then check if unique classes is in that table, then use the boolean to assign values
  if (deformode == FALSE){
    a = table(samples[[field_names[f]]])
  }else{
    a = table(samples_defor[[field_names[f]]])
  }
  
  class_check = ref_codes %in% names(a)
  ref_sample_count[f,class_check] = a
}

# Format and show
grid.newpage()
tt =  ttheme_default(base_size=14)
grid.table(round(ref_sample_count), theme=tt) # TODO REFORMAT AND SHOW! THIS IS PROBABLY THE KEY TO WIDE CI'S

# Calculate yearly area change and rate change (percentage of total area that is changing)
# Initialize zero matrix with year (03 to 16) * class (11) dimensions and proper row and column names
chg_area = matrix(0, nrow=length(years)-2, ncol=length(ref_codes), dimnames=list(years[3:16], ref_codes))
chg_rate = matrix(0, nrow=length(years)-2, ncol=length(ref_codes), dimnames=list(years[3:16], ref_codes))

# Iterate over strata labels
for(i in 1:length(area_ha)){
  chg_area[,i] = diff(area_ha[,i])
  chg_rate[,i] = chg_area[,i] / area_ha[1:14,i] * 100 #14 years
}

#Format and show
grid.newpage()
tt =  ttheme_default(base_size=14)
grid.table(round(chg_rate, digits=2), theme=tt)


## Calculate when break occurred in maps and compare to break in reference samples

break_calc = function(dataset){
  unq = unique(as.numeric(dataset))
  if (length(unq) == 1){
    return(0)
  }  
  else{
    # Find last column with the first code
    break_col = tail(which(dataset == unq[1]), n=1)
    # Calculate year of change (e.g if break was 12, break year was 2013) and format accordingly
    break_year = as.numeric(paste0("20",sprintf("%02d",break_col+1)))
    #print(sprintf("Break happened in %s", break_year))
    return(break_year)
  }
}

# Apply functions over rows on strata columns only
strata_breaks = apply(samples@data[, 39:53], MARGIN = 1, FUN = break_calc)

# Split CHGDATE to get year, then replace NA with 0 and substract, or convert to date and compare.
# Get YEAR of change from reference samples
ref_breaks = strsplit(samples$CHGDATE, "-")
# Get the first element of each list (i.e the year), unlist the output, pass as numeric
ref_breaks = as.numeric(unlist(lapply(ref_breaks, '[[', 1)))
# Replace NA's with zeros
ref_breaks[is.na(ref_breaks)] = 0
# Create table with the two datasets and add totals
change_cm = table(strata_breaks, ref_breaks)
change_cm = cbind(change_cm, total=rowSums(change_cm))
change_cm = rbind(change_cm, total=colSums(change_cm)) 
# Compare total breaks detected per year
total_ref = change_cm["total",]
total_ref = total_ref[-2]
total_breaks = cbind(total_ref, change_cm[,"total"])
colnames(total_breaks) = c("ref_breaks", "strata_breaks")

# Format table and display
tt= ttheme_default(base_size=18)
grid.newpage()
grid.table(round(change_cm, digits=2), theme=tt)
#png("numchange_strata_ref.png", width=1000, height = 1000, units = "px"); grid.table(change_cm, theme=tt); dev.off()


##############################################################################################################
# 6) CREATE SOME USEFUL PLOTS
# Most of these plots were recreated in python (add name of the script) bc the lack of dual axis support in ggplot. This code is left for 
# reference.

area_ha = cum_area
area_ci = cum_ci
area_lower = cum_lower
area_upper = cum_upper
margin_error = cum_me

## Plot area estimates with CI and margin of error in separate plot. Couldn't figure out how to do it with facets
# so I did it with grids. 

for(i in 1:length(area_ha)){
  
  # Define data and variables to use
  tempdf <- as.data.frame(cbind(years[2:16], area_ha[,i], area_lower[,i], area_upper[,i], margin_error[,i]))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Margin_error")
  
  # Plot areas with CI
  # Specify data, add "ribbon" with lower and higher CI and fill, then plot the estimated area with a line.
  # Other custom settings for number of breaks, label formatting and title.
  
  a <- ggplot(data=tempdf, aes(x=Years, y=Area_ha)) + 
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3) + geom_line() + 
    scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + 
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
    ylab("Area and 95% CI [ha]") + ggtitle(strata_names[[i]]) + expand_limits(y=0) +
    theme(plot.title = element_text(size=19), axis.title=element_text(size=16), axis.text=element_text(size=16))
  
  
  # Plot margin of error
  b <- ggplot(data=tempdf, aes(x=Years, y=Margin_error * 100)) + geom_line(size=1.1) + 
    scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + expand_limits(y=0) + ylab("Margin of error [%]") +
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
  filename = paste0(strata_names[[i]], "_areas_me", ".png")
  #png(filename, width=1000, height = 1000, units = "px"); plot(g); dev.off()

}


## Plot forest change per year. Yearly loss in 1 mostly equals yearly gain in 8+9, with the exception of a couple years
# Get only classes we're interested in
for_change = chg_area[,c(2, 8, 9)]
# Calculate how much of deforestation is not caused by pastures or regrowth and get absolute values
# Positive values (and also small values) in the "to other" column are attributed to rounding errors
# While larger ones are attributed to changes from forest to other land cover type (class 0)
for_to_other = rowSums(for_change[,1:3])
for_change = cbind(for_change, for_to_other)
for_change = abs(for_change[,2:4])
colnames(for_change) = c("To pasture", "To secondary forest", "To other")

# Melt and plot. 
for_change_melt = melt(for_change)
forest_plot <- ggplot(for_change_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
  geom_area(position="stack") + 
  scale_x_continuous(breaks=years[3:16], minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Annual forest conversion to other classes [ha]") + xlab("Years")+
  scale_fill_brewer(palette="GnBu", breaks=levels(for_change_melt$Var2), guide = guide_legend(reverse=T)) + 
  theme(legend.title=element_blank()) +
  theme(axis.title=element_text(size=15), axis.text=element_text(size=13), legend.text=element_text(size=13))

print(forest_plot)
#ggsave("forest_change.png", plot=forest_plot, device="png") 


## Plot net regrowth (net secondary forest). Yearly loss in 5 equals yearly gain in 14

net_rg = area_ha[,6] + area_ha[,10] + area_ha[,9] - area_ha[,11]

# Get only gain classes 
regr_area = area_ha[,c(6,9,10)]
names(regr_area) = c("Stable secondary forest", "Forest to secondary forest", "Other gains of secondary forest") #, "Loss of regrowth")

# Melt and plot
regr_area_melt = melt(as.matrix(regr_area))

regr_plot <- ggplot(regr_area_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
  geom_area(position="stack", alpha=0.8) + #geom_line(aes(x=Var1, y=732031.3), linetype=2) + 
  scale_x_continuous(breaks=years[2:16], minor_breaks = NULL)  + 
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
  scale_x_continuous(breaks=years[2:16], minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Loss of primary forest [ha]") + xlab("Years") +
  scale_fill_brewer(palette="GnBu", breaks=levels(for_loss_melt$Var2), guide = guide_legend(reverse=T)) + 
  theme(legend.title=element_blank()) +
  theme(axis.title=element_text(size=15), axis.text=element_text(size=13), legend.text=element_text(size=13))
  
print(forest_loss_plot)
#ggsave("forest_loss.png", plot=forest_loss_plot, device="png") 


# COMPARE MAPPED AREAS WITH ESTIMATES

# Create empty matrix to store mapped values and fill
mapped_areas = matrix(0, nrow=length(years[2:16]), ncol=nrow(mapped_areas_list[[1]]), byrow=T, dimnames = list(years[2:16], mapped_areas_list[[1]][,1]))

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
mapped_plots = cbind(years[2:16], mapped_plots)
colnames(mapped_plots) = c("Years", "Forest_to_pasture",  "Forest_to_regrowth")

mfp <- ggplot(data=mapped_plots, aes(x=Years, y=Forest_to_pasture)) + geom_line(size=1.1) + 
  scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + expand_limits(y=0) + ylab("Area [ha]") +
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=16))

mfr <- ggplot(data=mapped_plots, aes(x=Years, y=Forest_to_regrowth)) + geom_line(size=1.1) + 
  scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + expand_limits(y=0) + ylab("Area [ha]") +
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=16))

mdefor <- ggplot(data=mapped_plots, aes(x=Years, y=Forest_to_pasture + Forest_to_regrowth)) + geom_line(size=1.1) + 
  scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + expand_limits(y=0) + ylab("Area [ha]") +
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
  tempdf <- as.data.frame(cbind(years[2:16], area_ha[,i], area_lower[,i], area_upper[,i], mapped_areas[,'8']))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Mapped")
  
  # Specify data, add "ribbon" with lower and higher CI and fill, then plot the estimated area with a line.
  # Other custom settings for number of breaks, label formatting and title.
  a<-ggplot(data=tempdf, aes(x=Years, y=Area_ha))+ geom_line(aes(x=Years, y=Mapped),linetype=4 )+
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3) + geom_line() + 
    scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + 
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
for (f in field_names){
  a = samples@data[,f] == 8
  fpc = cbind(fpc, sum(a))
}

dtf = as.data.frame(cbind(years, t(fpc)))
colnames(dtf) = c("Years", "Count")


# Analyze strata breaks vs ref breaks per path/row
break_compare = as.data.frame(cbind(ref_breaks, strata_breaks, samples@data$PTRW))
colnames(break_compare) = c("ref_breaks", "strata_breaks", "ptrw")
break_compare$breakdif = break_compare$ref_breaks - break_compare$strata_breaks

# Function to calculate frequency of break differences. Sum gives total of pts
calc_break_time <- function(freqlist){
  pos = sum(freqlist > 2000) # Map didn't detect break but reference did (omission)
  neg = sum(freqlist < 0) # Reference didn't detect break but map did (commission)
  sm = sum(freqlist > 0 & freqlist< 2000)
  zr = sum(freqlist == 0)
  return(rbind(pos, neg, sm, zr))
}

breakdif_count = by(break_compare$breakdif, break_compare$ptrw, calc_break_time, simplify = FALSE) 
total_ptrw_pts = unlist(lapply(breakdif_count, sum))
bd_ratios = lapply(breakdif_count, function(x) x/sum(x)) # Calculate as ratios of the total
bd_ratios_melt =melt(bd_ratios) # Get that out of the ugly list
bd_ratios_df = dcast(bd_ratios_melt, L1~Var1) # Reshape to get an easier to manage df
bd_ratios_df = cbind(bd_ratios_df, total_ptrw_pts)
write.csv(bd_ratios_df, "break_ratios.csv")
        
#TODO

# - Plot of primary forest loss disagregated by class (requires operating over original mosaics)
