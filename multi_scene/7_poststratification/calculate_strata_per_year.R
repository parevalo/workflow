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

# Get unique classes through all the reference years (This won't have class 13 for that reason)
# and get numbered sequence
ref_codes = sort(unique(unlist(samples@data[field_names])))
ref_ind = seq_along(ref_codes)

calc_area_prop = function(strata, reference){

  # Obtain unique values in the STRATA field
  str_codes = sort(unique(strata))
  str_ind=seq_along(str_codes)
  
  # Initialize empty df for proportions per strata, and for sample variance per strata
  ref_prop = data.frame()
  ref_var = data.frame()
  
  # Compare the fields, iterate over "strata" and "reference" classes
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
# THE SCRIPT REQUIRES THE STRATA TO BE THE ORIGINAL STRATA, SO NO NEED TO ITERATE OVER RASTERS!! 
# THOSE RASTERS ARE ONLY NEEDED IF WE WANT TO CALCULATE THE ACCURACIES
ca = calc_area_prop(samples$final_strata_01_16_UTM18N, samples$ref_2003)

area_prop = data.frame()
ref_var_list = list()
filtered_ss = list()
for (i in (1:length(rast_names))){
  # Compare year strata with year reference. Field names must start at 2002, hence i+1
  out = calc_area_prop(samples$final_strata_01_16_UTM18N, samples[[field_names[i+1]]]) 
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

### 4) CALCULATE UNBIASED STANDARD ERROR FOR PROPORTION OF REFERENCE CLASS AREAS

# Formula works correctly, tested with unfiltered sample (i.e 1048 records) and it gives roughly the same CI
se_prop = data.frame()
#Iterate over years

for (y in 1:length(ref_var_list)){
  # Iterate over classes
  for (c in 1:length(refcodes)){
  N=sum(ss$pixels)
  # This is giving warnings of strata_pixels not being of same lenght of ref_var_list (refvar is shorter)!! CHECK!!
  # Probably means that class 14 is being incorrectly calculated?
  v = 1/N^2 * (sum(ss$pixels^2 * (1 - strata_pixels$x/ss$pixels) * (ref_var_list[[y]][,c] / strata_pixels$x)))
  se_prop[y,c] = sqrt(v)
  }
}

# Reassign names...
rownames(se_prop) = years[2:length(years)]
colnames(se_prop) = refcodes

# Total area in ha
N_ha = N * 30^2 / 100^2
# Calculate area in ha from area proportions
area_ha = area_prop * N_ha
# Calculate confidence interval in ha
area_ci = se_prop * 1.96 * N_ha
#Upper and lower CI
area_upper = area_ha + area_ci
area_lower = area_ha - area_ci

# Write results to csv
write.csv(area_ha, file="area_ha.csv")
write.csv(area_lower, file="area_lower.csv")
write.csv(area_upper, file="area_upper.csv")

# Reference sample count per year. Use something like this above to deal with varying number of classes per year

# Initialize zero matrix with year * class dimensions and proper row and column names
m = matrix(0, nrow=length(years), ncol=length(ref_codes), dimnames=list(years, ref_codes))

# Calculate the sample count
for (f in 1:(length(field_names))){
  # Get table, then check if unique classes is in that table, then use the boolean to assign values
  a = table(samples[[field_names[f]]])
  class_check = ref_codes %in% names(a)
  m[f,class_check] = a
  
}

# Calculate yearly area change and rate change (percentage of total area that is changing)
# Initialize zero matrix with year (03 to 16) * class (11) dimensions and proper row and column names
chg_area = matrix(0, nrow=length(years)-2, ncol=length(ref_codes), dimnames=list(years[3:16], ref_codes))
chg_rate = matrix(0, nrow=length(years)-2, ncol=length(ref_codes), dimnames=list(years[3:16], ref_codes))

for(i in 1:length(area_ha)){
  chg_area[,i] = diff(area_ha[,i])
  chg_rate[,i] = chg_area[,i] / area_ha[1:14,i] * 100
}


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

## 5) Plots

# Basic, ugly plot
plot(area_ha[,2], type="l", xlab="Years", ylab="Area in ha")
lines(area_upper[,2], col="red")
lines(area_lower[,2], col="red")

# Nice, good looking plots

# Create vector with name of classes, doesn't include class 13

strata_names = c("Other to other", "Stable forest", "Stable grassland", "Stable Urban + Stable other", 
                 "Stable pasture-cropland", "Stable regrowth", "Stable water", "Forest to pasture", 
                 "Forest to regrowth", "Pasture to all others", "Loss of regrowth")

# Calculate margin of error, plot along with the areas with CI (right below), use righ axis
margin_error = area_ci / area_ha

# Create each plot and save it to png, make y axis start at 0, break the axis for stable forest

for(i in 1:length(area_ha)){
  
  # Define data and variables to use
  tempdf <- as.data.frame(cbind(years[2:16], area_ha[,i], area_lower[,i], area_upper[,i]))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper")
  
  # Specify data, add "ribbon" with lower and higher CI and fill, then plot the estimated area with a line.
  # Other custom settings for number of breaks, label formatting and title.
  a<-ggplot(data=tempdf, aes(x=Years, y=Area_ha)) + 
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3) + geom_line() + 
    scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + 
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
    ylab("Area in ha") + ggtitle(strata_names[[i]]) + expand_limits(y=0)
    theme(plot.title = element_text(size=18), axis.title=element_text(size=13), axis.text=element_text(size=11))
  filename = paste0(strata_names[[i]], "_areasCI", ".png")
  print(a)
  ggsave(filename, plot=a, device="png") 
  
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
    ylab("Area in ha") + ggtitle(strata_names[[i]]) + expand_limits(y=0)
  theme(plot.title = element_text(size=18), axis.title=element_text(size=14), axis.text=element_text(size=12))
  
  
  # Plot margin of error
  b <- ggplot(data=tempdf, aes(x=Years, y=Margin_error)) + geom_line(size=1.1) + 
    scale_x_continuous(breaks=years[2:16], minor_breaks = NULL) + expand_limits(y=0)
  
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
for_to_other = rowSums(for_change[,1:3])
for_change = cbind(for_change, for_to_other)
for_change = abs(for_change[,2:4])
colnames(for_change) = c("To pasture", "To regrowth", "To other")

# Melt and plot. Fix palette
for_change_melt = melt(for_change)
forest_plot <- ggplot(for_change_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
  geom_area(position="stack", alpha=0.8) +
  scale_x_continuous(breaks=years[3:16], minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Annual forest conversion to other classes [ha]") + xlab("Years")+
  scale_fill_brewer(palette="Greens", breaks=levels(for_change_melt$Var2)) + theme(legend.title=element_blank())

print(forest_plot)

# Net regrowth. Yearly loss in 5 equals yearly gain in 14

net_rg = area_ha[,6] + area_ha[,10] + area_ha[,9] - area_ha[,11]

plot(years[2:16], area_ha[,6])
lines(years[2:16], net_rg)

# Get only classes we're interested in
regr_area = area_ha[,c(6,9,10,11)]
names(regr_area) = c("Stable regrowth", "Forest to regrowth", "Other to regrowth", "Loss of regrowth")
# Melt and plot
regr_area_melt = melt(as.matrix(regr_area))
# Substract loss of regrowth from the total, and then plot that but not stacked, maybe as a single line

regr_plot <- ggplot(regr_area_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
  geom_area(position="stack", alpha=0.8) + 
  scale_x_continuous(breaks=years[2:16], minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Total annual area in regrowth [ha]") + xlab("Years")+
  scale_fill_brewer(palette="Greens", breaks=levels(regr_area_melt$Var2)) + theme(legend.title=element_blank())

print(regr_plot)

#TODO
# - Find out WHERE the biggest omission and comission errors are happening, and their percentage with respect to the
# total sample size in that path-row
# - For the paper, remove plot titles and add to the y label (e.g. forest to pasture converstion [ha])
# - Make plot of net change in secondary forest
# - Plot of primary forest loss disagregated by class (requires operating over original mosaics)
