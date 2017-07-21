# Script to join the tables of interpreted samples to the original shapefiles

require(tidyverse)
require(rgdal)
require(raster)
require(rgeos)
require(grid)
require(gridExtra)

# Set working directories and vars
auxpath = "/media/paulo/785044BD504483BA/test/"
setwd("/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/biannual_samples/")
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


# Get row count per sample file to verify
samples_nrows = as.data.frame(do.call(rbind, lapply(csv_list, nrow)))
colnames(samples_nrows) = "count"
samples_nrows$period = samples_names

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
strata_pixels = aggregate(shp_list[[1]]$STRATUM, by=list(shp_list[[1]]$STRATUM), length)

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
area_upper = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 3)))
area_lower = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 4)))
margin_error = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 5)))

### PLOT AREAS

plot_mapped = F
allplots = list()

for(i in 1:length(ref_codes_all)){
  
  # Define data and variables to use
  tempdf <- as.data.frame(cbind(seq(1,7), area_ha[,i], area_lower[,i], area_upper[,i], margin_error[,i]))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Margin_error")
  
  # Plot areas with CI
  # Specify data, add "ribbon" with lower and higher CI and fill, then plot the estimated area with a line.
  # Other custom settings for number of breaks, label formatting and title.
  
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
    filename = paste0(savepath, strata_names[[i]], "_areas_me_step", step, "_", lut_name, "mapped_areas.png")
  } else {
    filename = paste0(savepath, strata_names[[i]], "_areas_me_step", step, "_", lut_name, ".png")
  }
  #png(filename, width=1000, height = 1000, units = "px"); plot(g); dev.off()
  
}

# Arrange plots into a single one to make it easier to compare differences between classes when a change is made
# NOT YET WORKING WITH THE MAPPED AREAS FOR SOME REASON BUT INDIVIDUALLY PLOTS ARE CORRECT
multiplots = grid.arrange(grobs=allplots, ncol=4)
#ggsave(paste0(savepath, "ALL_step", step, "_", lut_name, add_samples_suffix, ".png"), plot=multiplots,  width = 20, height = 10) 

