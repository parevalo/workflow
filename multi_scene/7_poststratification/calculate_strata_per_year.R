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
auxpath = "C:/test/"

if( .Platform$OS.type == "unix" )
    wd = "/media/paulo/785044BD504483BA/OneDrive/Lab/sample_may2016/interpreted_w_strata_23062016"
    auxpath = "/media/paulo/785044BD504483BA/test/"
    stratpath = "/home/paulo/workflow/multi_scene/7_poststratification/"
    #lutpath = "/home/paulo/workflow/multi_scene/data/for-nofor_lut_A.csv"
    lutpath = "/home/paulo/workflow/multi_scene/data/original_lut.csv"
    
    
setwd(wd)
source(paste0(stratpath, "functions.R"))
#source(paste0(stratpath, "testsource.R"))

#Set input names and suffixes. Add/remove lutA or lutB at the end, except for the prefix
# lut_name = "lutA"
# orig_stratif = paste0("final_strata_01_16_UTM18N_", lut_name)
# rast_prefix = "final_strata_annual_"
# rast_suffix = paste0("_UTM18N_", lut_name)
# pixcount_suffix = paste0("_pixcount_", lut_name, ".csv")
# pixcount_strata = paste0("strata_01_16", pixcount_suffix)

# Set input names and suffixes
lut_name = "original_lut"
orig_stratif = "buffered10_final_strata_01_16_UTM18N"
#orig_stratif = "final_strata_01_16_UTM18N"
rast_prefix = "final_strata_annual_"
rast_suffix = "_UTM18N"
pixcount_suffix = "_pixcount.csv"
pixcount_strata = paste0("buffered10_strata_01_16", pixcount_suffix)
#pixcount_strata = paste0("strata_01_16", pixcount_suffix)

# Set up important global variables
start = 2001
end = 2016
step = 1  # Number of years to do the analysis over
years = seq(start, end, step) 
ref_names = paste("ref_", years, sep="")
#cr = c(0,seq(5,15)) # For lutA
cr = c(7, 10, 12, 15) # Classes to ignore from the loaded area count tables

# List of original strata names
orig_strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                      "Stable pasture-cropland", "Stable secondary forest", "Stable water", "Forest to pasture", 
                      "Forest to secondary forest", "Gain of secondary forest", "All to unclassified", "Loss of secondary forest")


# # Strata names List doesn't have the "all to unclass." class bc it disappears when the ref. samples are collected.
strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                   "Stable pasture-cropland", "Stable secondary forest", "Stable water", "Forest to pasture", 
                   "Forest to secondary forest", "Gain of secondary forest", "Loss of secondary forest")

#strata_names = c("Stable forest", "Stable non-forest", "Forest loss", "Forest gain")

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

# Convert reference sample info into vectors containing proper labels per each of
# the years in our time period.
# Initialize empty vectors to store the data, current change and class code (start with the first!) 
# and name of the fields that store the class codes. Load LUT.

rows = nrow(samples)
df = vector()
field = vector()
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
        field[i] = calculate_strata(samples$CODE1[row], samples$CODE1[row], lut)
      }
    
      # If there is a change, compare each year to the current change year and update that one accordingly
      else {
        # Get the current date of change (initialized as the first change available)
        chg_date = unlist(strsplit(as.vector(samples$CHGDATE[row]), ','))[current_change] #742
        # Extract the YEAR from the date of change
        chg_year = na.omit(as.numeric(unlist(strsplit(chg_date, "[^0-9]"))))[1]
        
        # If we haven't reached change year yet
        if (years[i] < chg_year) {
          field[i] = calculate_strata(samples@data[codelist[current_code]][row,], samples@data[codelist[current_code]][row,], lut) 
        } 
        
        # If we JUST reached a change year
        else if (years[i] == chg_year) {
          field[i] = calculate_strata(samples@data[codelist[current_code]][row,], samples@data[codelist[current_code+1]][row,], lut) 
          # Check if we haven't reached the max number of recorded changes
          if (current_change < samples$NUMCHANGES[row]) {
            current_change = current_change + 1 
          }
          current_code = current_code + 1
        }
        
        # If we went past the last change date, use the last code 
        else if (years[i] > chg_year) {
          field[i] = calculate_strata(samples@data[codelist[samples$endcodecol[row]]][row,], samples@data[codelist[samples$endcodecol[row]]][row,], lut)
        }
      } 
  }
  df <- as.data.frame(rbind(df, field))
}

names(df) = ref_names

# Recalculate "original stratification", needed when we use a LUT different than the original
strata = vector()
for (row in 1:rows){
  strata[row] = calculate_strata(samples@data["CODE1"][row,], samples@data[codelist[samples$endcodecol[row]]][row,], lut)
}

# Attach table to shapefile 
samples@data[,ref_names] <- df
#writeOGR(samples, "sample_yearly_strata", "sample_yearly_strata", driver="ESRI Shapefile", overwrite_layer = T)

#############################################################################################################
# 2) READ ORIGINAL AND ANNUAL STRATA RASTERS AND EXTRACT THEIR VALUES TO THE SHAPEFILE
# Also calculate original strata size, weights, area proportions, and accuracies

# Get unique ref codes and map codes for all the years. Then create a single
# set of unique codes to be used in the confusion matrices
ref_codes = sort(unique(unlist(samples@data[ref_names])))
map_codes = sort(unique(unlist(samples@data[map_names])))
class_codes = sort(union(ref_codes, map_codes))

# Iterate over names and extract to shapefile. We need these to calculate 
# confusion matrices per year
map_names = character()
for (y in 1:(length(years)-1)){
  map_names[y] = paste0(rast_prefix, sprintf("%02d",y), "_", sprintf("%02d",y+1), rast_suffix)
  map = raster(paste0(auxpath, map_names[[y]], ".tif"))
  samples = extract(map, samples, sp = TRUE) 
}

# Read original stratification. Done last so that the reference and map fields are contiguous
samples = extract(raster(paste0(auxpath, orig_stratif, ".tif")), samples, sp=TRUE)

# Crosstab final strata and reference strata (for the same period 01-16) 
# Returns a square matrix with all the reference and map codes in the samples

ct = calc_ct(samples[[orig_stratif]], strata, class_codes)

# Create confusion matrices per year between reference and map labels
cm_list = list()
for (y in 1:(length(years) - 1)){
  cm_list[[y]] = calc_ct(samples[[map_names[y]]], samples[[ref_names[y]]], class_codes)
}

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
usr_acc = diag(aprop) / rowSums(aprop)
prod_acc = diag(aprop) / colSums(aprop)

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

# Calculate total area per stratum, making sure we only filter the classes
# present in the reference code list.
stratum_areas= ss$stratum[ss$stratum %in% ref_codes] *30^2 / 100^2
# map_bias = stratum_areas - area_ha['2016',] # doesn't work in all cases bc FIX!!
# 2016 is not always created if we aggregate the years.

# Write results to csv 
suffix = paste0("_step", step, "_", lut_name, ".csv")

#write.csv(area_ha, file=paste0("area_ha", suffix))
#write.csv(area_lower, file=paste0("area_lower", suffix))
#write.csv(area_upper, file=paste0("area_upper", suffix))


##############################################################################################################
# 5) CREATE SOME USEFUL TABLES

# Export table of  sample count, areas, percentages

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
  if (deformode == FALSE){
    a = table(samples[[ref_names[f]]])
  }else{
    a = table(samples_defor[[ref_names[f]]])
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
# THIS APPROACH IS IGNORING THE FACT THAT SOME SAMPLES HAVE MULTIPLE TRANSITONS

break_calc = function(dataset){
  unq = unique(as.numeric(dataset))
  if (length(unq) == 1){
    return(0)
  }  
  else{
    # Find last column with the first code
    break_col = tail(which(dataset == unq[1]), n=1)
    # Calculate REAL year of change (e.g. if last column with the first code is 
    # 12, then change happened in 2013-2014. Given that maps are created at the start of the year
    # then the actual change date is 2013. This is temporary, a more robust approach
    # would be better
    break_year = years[break_col] + 1
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
temprow = rep(0, length(years)+1)
change_cm = rbind(change_cm[1,],temprow, change_cm[2:16,])
rownames(change_cm)[1:2] = c(0,2001)
change_cm = cbind(change_cm, total=rowSums(change_cm))
change_cm = rbind(change_cm, total=colSums(change_cm)) 

# Format table and display
tt= ttheme_default(base_size=18)
grid.newpage()
grid.table(round(change_cm, digits=2), theme=tt)
png("numchange_strata_ref.png", width=1000, height = 1000, units = "px"); grid.table(change_cm, theme=tt); dev.off()

# Compare total breaks detected per year
total_ref = change_cm["total",]
total_breaks = cbind(total_ref, change_cm[,"total"])
colnames(total_breaks) = c("ref_breaks", "strata_breaks")


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
  #filename = paste0(strata_names[[i]], "_areas_me_step", step, ".png")
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


# Analyze strata breaks vs ref breaks per path/row
break_compare = as.data.frame(cbind(ref_breaks, strata_breaks, samples@data$PTRW))
colnames(break_compare) = c("ref_breaks", "strata_breaks", "ptrw")
break_compare$breakdif = break_compare$ref_breaks - break_compare$strata_breaks

# Calculate type of break difference per row: ommited change, commited change, 
# detected change with small ofset in timing (either positive or negative), 
# exact change match (first change detected in the same year) or stable match
# (no changes detected in reference or maps)
# THIS ANALYSIS ONLY COMPARES THE FIRST BREAK FOR SIMPLICITY'S SAKE
# DOESNT TAKE LABELS INTO ACCOUNT, JUST TIMING OF BREAKS

break_compare$type[break_compare$breakdif > 2000] = "omission"
break_compare$type[break_compare$breakdif < -2000] = "comission" # small negative numbers ARE ALSO SMALL OFFSET (e.g 2002 - 2004)
break_compare$type[(break_compare$breakdif > 0) & (break_compare$breakdif < 2000)] = "pos_offset" 
break_compare$type[(break_compare$breakdif < 0) & (break_compare$breakdif > -2000)] = "neg_offset" 
# Assume rows == 0 are exact change matches, then overwrite those that are stable matches
break_compare$type[break_compare$breakdif == 0] = "exact_change_match" 
break_compare$type[(break_compare$ref_breaks == 0) & (break_compare$strata_breaks == 0)] = "stable_match"

# Aggregate results and calculate total points per path-row 
breakdif_count = aggregate(break_compare$type, by=list(break_compare$ptrw), FUN="summary")
total_ptrw_pts = apply(breakdif_count[,2], MARGIN = 1, FUN = "sum")

# # Calculate number of each type of change as ratios of the total
bd_ratios = t(apply(breakdif_count[,2], MARGIN = 1, function(x) x/sum(x))) 
rownames(bd_ratios) = breakdif_count$Group.1
bd_ratios_df = cbind(bd_ratios, total_ptrw_pts)

# Write aggregated results to CSV and row results to shapefile
write.csv(bd_ratios_df, "LC_change_ratios_pathrow.csv")
samples@data[,names(break_compare)] <- break_compare
writeOGR(samples, "sample_yearly_strata_break_analysis", "sample_yearly_strata_break_analysis", driver="ESRI Shapefile", overwrite_layer = T)
     
# Re check change_cm and this, just in case
   
#TODO

# - Plot of primary forest loss disagregated by class (requires operating over original mosaics)
