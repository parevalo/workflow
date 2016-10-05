### This script has three main sections:
### 1) Read reference strata shapefile and calculate labels per year
### 2) Read yearly strata raster and extract values to the SpatialPoints object (i.e. shapefile) 
### 3) Calculate reference class area proportions and variance of ref samples per strata
### 4) Calculate unbiased standard error for proportion of reference class areas 

require(rgdal)
require(raster)
require(reshape2)
require(ggplot2)
require(gtable)
require(grid)
require(gridExtra)
library(xtable)

## 1) READ SHAPEFILE AND CALCULATE REFERENCE LABELS PER YEAR
# This section takes the sample shapefile and assigns the proper strata for each year, depending
# on when the model break happened. For now it's written to filter the pixels with only ONE CHANGE (1020 pts) 

# Read shapefile with reference strata and change info
#setwd("/home/paulo/sample_may2016/interpreted_w_strata_23062016")
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
  #map = raster(paste0("/home/paulo/test/", r, ".tif"))
  samples = extract(map, samples, sp = TRUE) 
}

# Crosstab reference sample vs final strata (DOESN'T MATCH EXCEL FILE BC THE RECORDS WERE FILTERED ABOVE. 
# If run with original, confusion matrices are the same)
ct = table(samples$final_strata_01_16_UTM18N, samples$strata)

# LOAD the total strata sample size produced from CountValues.py, bc calculating it here with hist() takes forever...
#ss = read.csv("/home/paulo/test/strata_01_16_pixcount.csv", header=TRUE, col.names=c("stratum", "pixels"))
ss = read.csv("C:/test/strata_01_16_pixcount.csv", header=TRUE, col.names=c("stratum", "pixels"))
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


## 3) CALCULATE REFERENCE CLASS AREA PROPORTIONS AND VARIANCE OF REFERENCE SAMPLES PER STRATA
# Returns a vector with length equal to the number of reference classes
# Reorganize to make more legible

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
      # Get proportion on samp_reference class present on that samp_strata
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
    # LEAVE THE SUM OF THE ENTIRE samp_strata?
    class_prop[r] = sum(fss$pixels * ref_prop[,r])/tot_area_pix
  }
  return(list(class_prop, ref_var, fss, ref_prop)) 
}

# Call the function for every reference year we want and get area proportions and sample variance
# NOTE the double square brackets to allow for substitution
# THE SCRIPT REQUIRES THE STRATA TO BE THE ORIGINAL STRATA, SO NO NEED TO ITERATE OVER RASTERS!! 
# THOSE RASTERS ARE ONLY NEEDED IF WE WANT TO CALCULATE THE ACCURACIES

# To calculate deforestation combined (e.g 8 + 9): Reclassify samples$final and samples$ref<year>, ss, strata_pixels.
defor_str = samples$final_strata_01_16_UTM18N
defor_str[defor_str == 8 | defor_str == 9] = 17

# Find reference samples = 8 or 9 for each year and change value to 17
samples_defor = samples
for (f in 1:(length(field_names))){
  defor_ind = samples_defor[[field_names[f]]] == 8 | samples_defor[[field_names[f]]] == 9
  samples_defor@data[defor_ind,field_names[f]] = 17
}

# Sum class 8 and 9 to produce "class 17", then delete 8 and 9. Do that for total pixels and sample counts
defor_str_totals = ss 
class17a = defor_str_totals[defor_str_totals$stratum == 8,] + defor_str_totals[defor_str_totals$stratum == 9,]
defor_str_totals = rbind(defor_str_totals, class17)
defor_str_totals = defor_str_totals[!(defor_str_totals$stratum == 8 | defor_str_totals$stratum == 9),]

defor_samp_totals = strata_pixels 
class17b = defor_samp_totals[defor_samp_totals$Group.1 == 8,] + defor_samp_totals[defor_samp_totals$Group.1 == 9,]
defor_samp_totals = rbind(defor_samp_totals, class17b)
defor_samp_totals = defor_samp_totals[!(defor_samp_totals$Group.1 == 8 | defor_samp_totals$Group.1 == 9),]

# CONDITION THAT MODIFIES THE OUTPUT FROM HERE! Temporary fix to make it easier to run both cases if necessary
# If true, uses class 8 and 9 together as deforestation, labeled as class 17. Otherwise keeps them separate
deformode = FALSE

# Get unique classes through all the reference years (This won't have class 13 for that reason) and get numbered sequence
if (deformode == TRUE){
  ref_codes = sort(unique(unlist(samples_defor@data[field_names])))
}else { 
  ref_codes = sort(unique(unlist(samples@data[field_names])))
}

area_prop = data.frame()
ref_var_list = list()
filtered_ss = list()
ref_prop_list = list()
for (i in (1:length(rast_names))){
  # Compare year strata with year reference. Field names must start at 2002, hence i+1
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

### 4) CALCULATE UNBIASED STANDARD ERROR FOR PROPORTION OF REFERENCE CLASS AREAS

calc_se_prop = function(strata_totals, sample_totals){

  # Formula works correctly, tested with unfiltered sample (i.e 1048 records) and it gives roughly the same CI
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
  suffix = "_defor.csv"  
}else { 
  suffix = "_regular.csv"
}

write.csv(area_ha, file=paste0("area_ha", suffix))
write.csv(area_lower, file=paste0("area_lower", suffix))
write.csv(area_upper, file=paste0("area_upper", suffix))

# Export table of  sample count, areas, percentages
orig_strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                      "Stable pasture-cropland", "Stable secondary forest", "Stable water", "Forest to pasture", 
                      "Forest to secondary forest", "Gain of secondary forest", "All to unclassified", "Loss of secondary forest")
stratum_areas= ss$pixels *30^2 / 100^2
stratum_percentages=round(ss$pixels / tot_area_pix * 100, digits=3) 
strata_table = as.data.frame(cbind(stratum_areas, stratum_percentages, strata_pixels$x))
strata_table$stratum_areas = format(strata_table$stratum_areas, scientific = FALSE, big.mark = ",")
#Math notation in column name can only be displayed if it's on its own, and not with mixed text...
colnames(strata_table) = c("Area (ha)", "Area / W[h]", "Sample size (nh)") 
rownames(strata_table) = orig_strata_names
windowsFonts(Times=windowsFont("TT Times New Roman"))  #clearly, only required for windows machines
# Parse = TRUE required to display math notation
tt=ttheme_default(core=list(fg_params=list(font="Times", fontface="plain", fontsize=14)),
                  colhead=list(fg_params=list(font="Times", fontface="bold", fontsize=14, parse=TRUE)),
                  rowhead=list(fg_params=list(font="Times", fontface="plain", fontsize=14)))
grid.newpage()
grid.table(strata_table, theme=tt)  

# Create table in Latex instead, and produce the pdf there, much easier than grid.table
# Need to escape special characters, including backslash itself (e.g. $\\alpha$)
colnames(strata_table) = c("Area [ha]", "Area / $W_h$ [\\%]", "Sample size ($n_h$)") 
print(xtable(strata_table, digits=c(0,2,2,0)),type = "latex",sanitize.text.function=function(x){x})

#calculate map bias and create accuracy table with margin of error, only for strata 2001-2016
margin_error = area_ci / area_ha 
map_bias = stratum_areas - area_ha['2016',]
accuracies = as.data.frame(cbind(usr_acc, prod_acc))
accuracy_table=cbind(accuracies[-11,], t(map_bias), t(margin_error['2016',]*100))
colnames(accuracy_table) = c("User's accuracy", "Producer's accuracy", "Map bias", "Margin of error")
rownames(accuracy_table) = orig_strata_names[-11]
accuracy_table = format(accuracy_table, scientific = FALSE, big.mark = ",", digits=2)
grid.newpage()
grid.table(accuracy_table, theme=tt)  
print(xtable(accuracy_table, digits=c(0,2,2,0,1)),type = "latex",sanitize.text.function=function(x){x})

# Reference sample count per year. Use something like this above to deal with varying number of classes per year

# Initialize zero matrix with year * class dimensions and proper row and column names
m = matrix(0, nrow=length(years), ncol=length(ref_codes), dimnames=list(years, ref_codes))

# Calculate the sample count, not sure why I did this, not used yet!
for (f in 1:(length(field_names))){
  # Get table, then check if unique classes is in that table, then use the boolean to assign values
  if (deformode == FALSE){
    a = table(samples[[field_names[f]]])
  }else{
    a = table(samples_defor[[field_names[f]]])
  }
  
  class_check = ref_codes %in% names(a)
  m[f,class_check] = a
}

# Calculate yearly area change and rate change (percentage of total area that is changing)
# Initialize zero matrix with year (03 to 16) * class (11) dimensions and proper row and column names
chg_area = matrix(0, nrow=length(years)-2, ncol=length(ref_codes), dimnames=list(years[3:16], ref_codes))
chg_rate = matrix(0, nrow=length(years)-2, ncol=length(ref_codes), dimnames=list(years[3:16], ref_codes))

# Iterate over strata labels
for(i in 1:length(area_ha)){
  chg_area[,i] = diff(area_ha[,i])
  chg_rate[,i] = chg_area[,i] / area_ha[1:14,i] * 100 #14 years
}

#Format and save table?
tt =  ttheme_default(base_size=14)
grid.table(round(chg_rate, digits=2), theme=tt)


# Calculate when break occurred in maps and compare to break in reference samples
# (trying to get rid of loops now...)

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

# Format table and save
tt= ttheme_default(base_size=18)
grid.newpage()
png("numchange_strata_ref.png", width=1000, height = 1000, units = "px"); grid.table(change_cm, theme=tt); dev.off()


## 5) PLOTS

# Nice, good looking plots

# Create vector with name of classes, doesn't include class 13
if (deformode == FALSE){
  strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                 "Stable pasture-cropland", "Stable regrowth", "Stable water", "Forest to pasture", 
                 "Forest to regrowth", "All others to regrowth", "Loss of regrowth")
  }else{
  strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                   "Stable pasture-cropland", "Stable regrowth", "Stable water", "All others to regrowth", 
                   "Loss of regrowth", "Deforestation")
  
}

# Create each plot and save it to png, make y axis start at 0, break the axis for stable forest
# MODIFIED TO INCLUDE MAPPED DEFORESTATION, FIX TO INCLUDE EACH OF THE MAPPED CLASSES!

for(i in 1:length(area_ha)){
  
  # Define data and variables to use
  tempdf <- as.data.frame(cbind(years[2:16], area_ha[,i], area_lower[,i], area_upper[,i], mapped_areas[,'8'] + mapped_areas[,'9']))
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

# Same plot but with two panels, one including the margin error. Plots weren't created with ggplot facets bc
# I couldn't find/figure out how to modify properties of individual facets

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
  png(filename, width=1000, height = 1000, units = "px"); plot(g); dev.off()

}


# Forest change. Yearly loss in 1 mostly equals yearly gain in 8+9, with the exception of a couple years
# Get only classes we're interested in
for_change = chg_area[,c(2, 8, 9)]
# Calculate how much of deforestation is not caused by pastures or regrowth and get absolute values
# Positive values (and also small values) in the "to other" column are attributed to rounding errors
# While larger ones are attributed to changes from forest to other land cover type (class 0)
for_to_other = rowSums(for_change[,1:3])
for_change = cbind(for_change, for_to_other)
for_change = abs(for_change[,2:4])
colnames(for_change) = c("To pasture", "To regrowth", "To other")

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
ggsave("forest_change.png", plot=forest_plot, device="png") 

# Net regrowth. Yearly loss in 5 equals yearly gain in 14

net_rg = area_ha[,6] + area_ha[,10] + area_ha[,9] - area_ha[,11]

plot(years[2:16], area_ha[,6])
lines(years[2:16], net_rg)

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
ggsave("net_regrowth_secondaryforest.png", plot=regr_plot, device="png") 


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
  + geom.line
print(forest_loss_plot)
ggsave("forest_loss.png", plot=forest_loss_plot, device="png") 

# Load and plot MAPPED rates to compare with estimated rates

setwd("C:/test")
mapped_areas_list = list()
filenames = dir("C:/test/", pattern="*_pixcount.csv")
for(i in 1:length(filenames)){
  mapped_areas_list[[i]] = read.csv(filenames[i], header=TRUE, col.names=c("stratum", "pixels"))
  
}

# Create empty matrix to store values and fill

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

# Melt and plot. 
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
ggsave("mapped_forest_change.png", plot=mforest_plot, device="png") 

# Plot MAPPED areas of forest to pasture, forest to regrowth and deforestation
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
ggsave("mapped_forest_to_pasture.png", plot=mfp, device="png") 
ggsave("mapped_forest_to_regrowth.png", plot=mfr, device="png") 
ggsave("mapped_deforestation.png", plot=mdefor, device="png") 


print(mfp)
ggsave("mapped_forest_to_pasture.png", plot=mfp, device="png") 

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
