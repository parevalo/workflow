# Script to process the biannual samples and calculate area estimates and
# accuracies

require(tidyverse)
require(rgdal)
require(raster)
require(rgeos)
require(grid)
require(gridExtra)

# Set working directories and vars
auxpath = "/media/paulo/785044BD504483BA/test/"
setwd("/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/biannual_sampling/")
stratpath = "/home/paulo/workflow/multi_scene/7_poststratification/"
source(paste0(stratpath, "functions.R"))
source(paste0(stratpath, "input_variables_original_buffer3B.R"))
lut = read.table(lutpath, header = T, sep = ",")

start = 2001
end = 2016
step = 2
years = seq(start, end, step) 

# Load original sample shapefiles and csvs with interpretations
# CSV files should have only 1050 rows, i.e. filtered the samples that intersect scenes

samples_names = character()
csv_names = character()
shp_list = list()
csv_list = list()
short_years = substr(years, 3,4) # Get years in two digit format
periods = paste0(short_years[-8], "_", short_years[-1])
pixcount_list = list()

for (i in 1:(length(periods))){
  fname = paste0("sample_", periods[i])

  # Load shapes with map strata
  samples_names[i] = fname
  shp_list[[i]] = readOGR(paste0(auxpath, "biannual_samples/"), samples_names[i])

  # Load csvs with reference strata
  csv_list[[i]] <- read.csv(paste0(samples_names[i], ".csv"))
  
  # Load csvs with map pixel count
  fname = paste0("strata_buffered_", periods[i], "_pixcount.csv")
  pixcount_list[[i]] = read.csv(paste0(auxpath,"biannual_samples/", fname), header=TRUE, col.names=c("stratum", "pixels"))  
  pixcount_list[[i]] = pixcount_list[[i]][!(pixcount_list[[i]]$stratum %in% cr),]     

}


# Get number of unique ID's per file to verify they add to 1050.
samples_uniqueids = as.data.frame(do.call(rbind, lapply(csv_list, function(x) length(unique(x[,"ID"])))))
colnames(samples_uniqueids) = "count"
samples_uniqueids$period = samples_names
if (all(samples_uniqueids$count == 1050)){
  print("Unique ID count matches total number of samples")
  }else{
  errorcount = whichs(samples_uniqueids$count != 1050)
  print(paste0(samples_uniqueids[errorcount, 'period'], "does not have 1050 unique ID's"))
}

# TEMPORARY, USE TO CHANGE CODES AND SEE THE RESULTS IN AREAS
# chg_ids = get_condition_rows(1,14,8)$ID
# ids = which(csv_list[[1]]$ID %in% chg_ids)
# csv_list[[1]][ids, 'CODE1'] = 5


# Fnc to calculate strata for each year, given that there may or may not be
# land cover change in each period. Then do the calculation.

calc_strata_aux = function(df, lut){
  refstrata = vector(length = length(df$CODE1))
  ind_na = is.na(df$CODE2)
  # If second label is NA, just use the first to calculate stratum
  refstrata[ind_na] = calc_strata(df$CODE1, df$CODE1, lut)[ind_na]
  refstrata[!ind_na] = calc_strata(df$CODE1, df$CODE2, lut)[!ind_na]
  refstrata_id = as.data.frame(cbind(df$ID, refstrata))
  colnames(refstrata_id) = c("ID", "ref_strata")
  return(refstrata_id)
}

ref_strata = lapply(csv_list, calc_strata_aux, lut=lut)

# Join map strata and reference strata to calculate accuracies

join_ref_map_strata = function(map_shp, refstrata_id){
  map_shp@data = inner_join(map_shp@data, refstrata_id, by="ID", copy=T)
  return(map_shp)
}

# Use mapply to use each corresponding element of the two lists
shp_list_ref = mapply(join_ref_map_strata, shp_list, ref_strata)

# # TEMPORARY DEFORMODE. collapse ref and map labels, and pixcount 
# apply_deformod = function(df){
#   df$STRATUM[df$STRATUM == 9] = 8
#   df$ref_strata[df$ref_strata == 9] = 8
#   return(df)
# }
# 
# apply_deformode2 = function(df){
#   w8 = df[,1] == 8
#   w9 = df[,1] == 9
#   df[w8, 2] = df[w8, 2] + df[w9, 2]
#   df = df[-(which(w9)),]
#   return(df)
# }
# 
# # Overwrite files to avoid creating new calls to the functions below
# pixcount_list = lapply(pixcount_list, apply_deformode2)
# shp_list_ref = lapply(shp_list_ref, apply_deformod)


# Get lists of ref and map codes per period. This can probably be simplified
get_ref_codes = function(shp){
  refcodes = sort(unique(shp@data$ref_strata))
  return(refcodes)
}

get_map_codes = function(shp){
  mapcodes = sort(unique(shp@data$STRATUM))
  return(mapcodes)
}

ref_codes = lapply(shp_list_ref, get_ref_codes)
map_codes = lapply(shp_list_ref, get_map_codes)
ref_codes_all = sort(unique(unlist(ref_codes)))
map_codes_all = sort(unique(unlist(map_codes)))

# Crosstab map and ref labels for each period
create_cm = function(shape, refcodes, mapcodes){
  class_codes = sort(union(refcodes, mapcodes))
  cm = calc_ct(shape@data$STRATUM, shape@data$ref_strata, class_codes)
  return(cm)
}

cm_list = mapply(create_cm, shp_list_ref, ref_codes, map_codes, SIMPLIFY = F)

# Check that the sample allocation is correct
map_sample_count = as.data.frame(do.call(rbind, lapply(cm_list, rowSums)))

# Check all the maps have the same total area in pixels and use that value
map_total_pix = as.data.frame(do.call(rbind, lapply(pixcount_list, colSums)))
tot_area_pix = map_total_pix[1,2]

# Create single variable with sample allocation
strata_pixels = aggregate(shp_list_ref[[1]]$STRATUM, by=list(shp_list_ref[[1]]$STRATUM), length)

# Calculate strata weights in percentage as an aid to interpret the area plots.
strata_weights = as.data.frame(do.call(rbind, lapply(pixcount_list, function(x) (x[,2]/tot_area_pix)*100)))
colnames(strata_weights) = map_codes_all


# Calculate areas and accuracies, cant vectorize it the way the fncs are written
prop_out = list()
se_prop = list()
areas_out = list()
accuracies_out = list()

# ref_codes_all used to make sure output tables have same dimensions
for (p in 1:length(periods)){
  prop_out[[p]] = calc_props_and_vars(shp_list_ref[[p]]$STRATUM, shp_list_ref[[p]]$ref_strata, shp_list[[p]]$STRATUM, 
                                  pixcount_list[[p]], strata_pixels, ref_codes_all)
  
  se_prop[[p]] = calc_se_prop(pixcount_list[[p]], strata_pixels, prop_out[[p]][[2]], ref_codes_all, tot_area_pix)
  areas_out[[p]] = calc_unbiased_area(tot_area_pix, prop_out[[p]][[11]], se_prop[[p]]) 
  accuracies_out[[p]] = calc_accuracies(pixcount_list[[p]], strata_pixels, ref_codes_all, tot_area_pix,
                                    prop_out[[p]][[1]], prop_out[[p]][[2]], prop_out[[p]][[3]], prop_out[[p]][[4]],
                                    prop_out[[p]][[5]], prop_out[[p]][[6]], 
                                    prop_out[[p]][[7]], prop_out[[p]][[8]],
                                    prop_out[[p]][[9]], prop_out[[p]][[10]])

}

# Get some key results in a readable format

usr_acc = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 4)))
prod_acc = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 7)))

area_ha = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 1)))
ci_ha = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 2)))
area_upper = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 3)))
area_lower = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 4)))
margin_error = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 5)))

### PLOT AREAS

plot_areas = function(xlabels, areaha, lower, upper, me){
  #
  tempdf <- as.data.frame(cbind(seq(1,length(xlabels)), areaha, lower, upper, me))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Margin_error")
  
  ind_area = which(tempdf$Area_ha < 0)
  ind_lower = which(tempdf$Lower < 0)
  
  tempdf[union(ind_area, ind_lower), 'Area_ha'] = NA
  tempdf[ind_lower, 'Lower'] = NA
  tempdf[ind_lower, 'Upper'] = NA

  
  print(tempdf)
  # Plot areas with CI. Add dashed lines on top of ribbon for places where CI < 0
  # could also use geom_pointrange
  
  area_plot <- ggplot(data=tempdf, aes(x=Years, y=Area_ha)) +  
    geom_line(data=tempdf[!is.na(tempdf$Area_ha),], aes(x=Years, y=Lower), linetype=8) +
    geom_line(data=tempdf[!is.na(tempdf$Area_ha),], aes(x=Years, y=Upper), linetype=8) +
    geom_line(data=tempdf[!is.na(tempdf$Area_ha),], aes(x=Years, y=Area_ha), linetype=8) +
    #geom_line(data=tempdf, aes(x=Years, y=Upper), color="deepskyblue4") +
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3) + geom_line() + geom_point(shape=3, size=4) +
    scale_x_continuous(breaks=seq(1,length(periods)), labels=periods, minor_breaks = NULL) + 
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
    ylab("Area and 95% CI [ha]") + ggtitle(strata_names[[i]]) + expand_limits(y=0) +
    theme(plot.title = element_text(size=19), axis.title=element_text(size=16), axis.text=element_text(size=16))
  
  # Put plots together on a list and change some properties in order to save them together as a single plot
  area_plot2 = area_plot + theme(axis.title=element_blank(), axis.text.x=element_text(size=7), axis.text.y=element_text(size=10)) 
  
  # Plot margin of error
  me_plot <- ggplot(data=tempdf, aes(x=Years, y=Margin_error * 100)) + geom_line(size=1.1) + 
    scale_x_continuous(breaks=seq(1,length(periods)), labels=periods, minor_breaks = NULL) + expand_limits(y=0) + ylab("Margin of error [%]") +
    theme(axis.title=element_text(size=16), axis.text=element_text(size=16))
  
  # To remove grid and background add this:
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()
  
  # Use gtable to stack plots together with matching extent and save  
  # g1 <- ggplotGrob(area_plot)
  # g2 <- ggplotGrob(me_plot)
  # g <- rbind(g1, g2, size="first") # stack the two plots
  # g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
  # 
  grid.newpage()
  grid.draw(area_plot)
  
  filename = paste0("plots/", strata_names[[i]], "_areas_me_", lut_name, ".png")
  #png(filename, width=1000, height = 1000, units = "px"); plot(g); dev.off()
  
}


testplot = plot_areas(periods, area_ha[,11], area_lower[,11], area_upper[,11], margin_error[,11])


plot_mapped = F
allplots = list()

for(i in 1:length(ref_codes_all)){
  
  # Define data and variables to use
  tempdf <- as.data.frame(cbind(seq(1,length(periods)), area_ha[,i], area_lower[,i], area_upper[,i], margin_error[,i]))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Margin_error")
  
  # Plot areas with CI. Add dashed lines on top of ribbon for places where CI <
  
  a <- ggplot(data=tempdf, aes(x=Years, y=Area_ha)) + 
    geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3) + geom_line() +
    scale_x_continuous(breaks=seq(1,length(periods)), labels=periods, minor_breaks = NULL) + 
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) +
    ylab("Area and 95% CI [ha]") + ggtitle(strata_names[[i]]) + expand_limits(y=0) +
    theme(plot.title = element_text(size=19), axis.title=element_text(size=16), axis.text=element_text(size=16))
  
  # Plot MAPPEd areas
  if(plot_mapped == T){
    m = geom_line(mapping = aes(x=Years, y=mapped_areas[,i]), linetype=2)
    a = a+m
  }
  
  # Put plots together on a list and change some properties in order to save them together as a single plot
  a2 = a + theme(axis.title=element_blank(), axis.text.x=element_text(size=7), axis.text.y=element_text(size=10)) 
  allplots[[i]] = a2
  
  # Plot margin of error
  b <- ggplot(data=tempdf, aes(x=Years, y=Margin_error * 100)) + geom_line(size=1.1) + 
    scale_x_continuous(breaks=seq(1,length(periods)), labels=periods, minor_breaks = NULL) + expand_limits(y=0) + ylab("Margin of error [%]") +
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
  
  if(plot_mapped == T){
    filename = paste0("plots/", strata_names[[i]], "_areas_me_", lut_name, "mapped_areas.png")
  } else {
    filename = paste0("plots/", strata_names[[i]], "_areas_me_", lut_name, ".png")
  }
  #png(filename, width=1000, height = 1000, units = "px"); plot(g); dev.off()
  
}

# Arrange plots into a single one to make it easier to compare differences between classes when a change is made
# NOT YET WORKING WITH THE MAPPED AREAS FOR SOME REASON BUT INDIVIDUALLY PLOTS ARE CORRECT
multiplots = grid.arrange(grobs=allplots, ncol=4)
#ggsave(paste0("plots/", "ALL_", lut_name, ".png"), plot=multiplots,  width = 20, height = 10) 








# Fnc to extract and put together the same column name from multiples df on a list
get_dflist_columns = function(df, colname){
  subcols = as.data.frame(do.call(rbind, lapply(df, function(x) x[,colname])))
  return(t(subcols))
}

ref2 = get_dflist_columns(cm_list, '2')
ref3 = get_dflist_columns(cm_list, '3')
ref4 = get_dflist_columns(cm_list, '4')
ref5 = get_dflist_columns(cm_list, '5')
ref6 = get_dflist_columns(cm_list, '6')
ref8 = get_dflist_columns(cm_list, '8')
ref9 = get_dflist_columns(cm_list, '9')
ref11 = get_dflist_columns(cm_list, '11')
ref14 = get_dflist_columns(cm_list, '14')
  

# Function to find the rows that meet a condition of map and reference labels for a given period
# Useful to check them manually in TSTools and find errors.
get_condition_rows = function(period, mapcode, refcode){
  match_rows = which(shp_list_ref[[period]]@data[,"STRATUM"] == mapcode & shp_list_ref[[period]]@data[,"ref_strata"] == refcode)
  match_ids = shp_list_ref[[period]]@data[match_rows, "ID"]
  return(csv_list[[period]][csv_list[[period]]$ID %in% match_ids,])
}

for(i in 1:length(periods)){
  print(periods[i])
  print(get_condition_rows(i, 5, 14))
}

get_condition_rows(2,4,9)

# Determine what transitions are being calculated in the maps for class 0, and if most of them are forest to grassland
# then add them to class 8 and see what happens. 

# Change plots to be error bars with stars where there is an actual data point and nothing when the estimate is not different from zero. Connect
# dots and the error bars with lines. Change the labels to show the middle year instead of the period
  