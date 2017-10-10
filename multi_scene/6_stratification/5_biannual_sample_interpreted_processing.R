# Script to process the biannual samples and calculate area estimates and
# accuracies

require(tidyverse)
require(rgdal)
require(raster)
require(rgeos)
require(grid)
require(gridExtra)
require(reshape2)
require(xtable)

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
periods_long = paste0(years[-8], "-", years[-1])
pixcount_list = list()
mapped_areas = list()

# Classes to be removed from the pixcount list bc they are not estimated
cr_extra = c(13,16) 

for (i in 1:(length(periods))){
  fname = paste0("sample_", periods[i])

  # Load shapes with map strata
  samples_names[i] = fname
  shp_list[[i]] = readOGR(paste0(auxpath, "biannual_samples/"), samples_names[i])

  # Load csvs with reference strata
  csv_list[[i]] = read.csv(paste0("katelyn_revised/", samples_names[i], ".csv"))
  
  # Add period column to facilitate some operations
  csv_list[[i]]$period = periods[[i]]
  
  # Load csvs with map pixel count
  fname = paste0("strata_buffered_", periods[i], "_pixcount.csv")
  pixcount_list[[i]] = read.csv(paste0(auxpath,"biannual_samples/", fname), header=TRUE, col.names=c("stratum", "pixels"))  
  pixcount_list[[i]] = pixcount_list[[i]][!(pixcount_list[[i]]$stratum %in% cr),] 
  mapped_areas[[i]] = pixcount_list[[i]][!(pixcount_list[[i]]$stratum %in% cr_extra),] 

}

# Reformat and actually calculate mapped areas
mapped_areas = as.data.frame(do.call(rbind, lapply(mapped_areas, '[[', 2))) * 30^2 / 100^2


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
# land cover change in each period. Then do the calculation. ADD FOR-NOTFOR LUT HERE!

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

## TEMPORARY DEFORMODE. collapse ref and map labels, and pixcount 
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

## END of DEFORMODE


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
tot_area_ha = tot_area_pix * 30^2 / 100^2

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

overall_acc = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 1)))
overall_acc_lower = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 2)))
overall_acc_upper = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 3)))

usr_acc = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 4)))
usr_acc_lower = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 5)))
usr_acc_upper = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 6)))

prod_acc = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 7)))
prod_acc_lower = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 8)))
prod_acc_upper = as.data.frame(do.call(rbind, lapply(accuracies_out, '[[', 9)))

area_ha = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 1)))
ci_ha = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 2)))
area_upper = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 3)))
area_lower = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 4)))
margin_error = as.data.frame(do.call(rbind, lapply(areas_out, '[[', 5)))


############################# PLOT AREAS

plot_areas = function(totareaha, xlabels, areaha, lower, upper, mappedarea, me, miny, maxy, plot_title, plotmode, biglabels){
  # Need two copies bc of the complexity of the graph
  tempdf = as.data.frame(cbind(seq(1,length(xlabels)), areaha, lower, upper, mappedarea, me))
  names(tempdf) = c("Years", "Area_ha", "Lower", "Upper", "Mapped_area", "Margin_error")
  tempdf2 = tempdf
  
  # Find rows where the CI or the area go below 0 and assign NA's
  ind_area = which(tempdf$Area_ha < 0)
  ind_lower = which(tempdf$Lower < 0)
  
  tempdf[union(ind_area, ind_lower), 'Area_ha'] = NA
  tempdf[ind_lower, 'Lower'] = NA
  tempdf[ind_lower, 'Upper'] = NA
    
  # Define breaks manually so that the two y axes grids match
  bks = seq(miny, maxy, length.out = 6)
  bks2 = bks / totareaha * 100
  
  # Conditonals for pontus modes, yscale is shared for the two 
  yscale = scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}, breaks=bks, limits=c(miny,maxy), 
                              sec.axis = sec_axis(~./totareaha * 100, breaks=bks2, labels=function(n){format(n, digits=1)}))
  ribbon = geom_ribbon(data=tempdf, aes(x=Years, ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3)
  
  # Remove CI and area when it intersects with zero
  if (plotmode == 1){
    
    lowerline = geom_line(data=tempdf[!is.na(tempdf$Area_ha),], aes(x=Years, y=Lower), linetype=8)
    upperline = geom_line(data=tempdf[!is.na(tempdf$Area_ha),], aes(x=Years, y=Upper), linetype=8) 
    centerline = geom_line(data=tempdf[!is.na(tempdf$Area_ha),], aes(x=Years, y=Area_ha), linetype=8)
  
  # Or keep the original outlines but remove the fill
  } else if (plotmode == 2) {
    
    lowerline = geom_line(data=tempdf2, aes(x=Years, y=Lower), linetype=8)
    upperline = geom_line(data=tempdf2, aes(x=Years, y=Upper), linetype=8) 
    centerline = geom_line(data=tempdf2, aes(x=Years, y=Area_ha), linetype=8)
  
  # Or keep the simplified original with changes in axes or data
  } else {
    lowerline = geom_blank()
    upperline = geom_blank()
    centerline = geom_line(data=tempdf2, aes(x=Years, y=Area_ha))
    yscale = scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) 
    ribbon = geom_ribbon(data=tempdf2, aes(x=Years, ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3)
  }
    
  # Plot areas with CI. Add dashed lines on top of ribbon for places where CI < 0
  area_plot = ggplot() +  
    lowerline + upperline + centerline + ribbon +
    geom_line(data=tempdf2, aes(x=Years, y=Mapped_area), colour="red") + 
    geom_point(data=tempdf, aes(x=Years, y=Area_ha), shape=3, size=4, stroke=1) +
    scale_x_continuous(breaks=seq(1,length(xlabels)), labels=xlabels, minor_breaks = NULL) +
    yscale + ylab("Area and 95% CI [ha]") +
    ggtitle(plot_title)  + geom_hline(yintercept = 0, size=0.3) +
    theme(plot.title = element_text(size=19), axis.title=element_text(size=16), axis.text=element_text(size=16)) 
  
  # Put plots together on a list and change some properties in order to save them together as a single plot
  area_plot_small = area_plot + theme(axis.title=element_blank(), axis.text.x=element_text(size=7), axis.text.y=element_text(size=10)) 
  
  # Plot margin of error
  me_plot = ggplot(data=tempdf, aes(x=Years, y=Margin_error * 100)) + geom_line(size=1.1) + 
    scale_x_continuous(breaks=seq(1,length(xlabels)), labels=xlabels, minor_breaks = NULL) + ylim(0, 200) + 
    ylab("Margin of error [%]") + ggtitle(plot_title) +
    theme(axis.title=element_text(size=16), axis.text=element_text(size=16))
  
  if (plotmode == T){
      return_list = list(area_plot_small, me_plot)
  } else {
    if (biglabels == T){
      return_list = list(area_plot + expand_limits(y=0), me_plot)
    } else {
      return_list = list(area_plot_small + expand_limits(y=0), me_plot)
    }
  }

  return(return_list)

}

# Vector of max and min y axis values for pontus modes
# Selected to guarantee that one of the breaks (6 total) is zero
maxy_vect1 = c(12000, 45000000, 4500000, 300000, 4500000, 4500000, 4500000, 300000, 300000, 300000, 300000)
maxy_vect2 = c(12000, 45000000, 4500000, 400000, 4500000, 4500000, 4500000, 400000, 400000, 400000, 400000)
miny_vect2 = c(-3000, 0, 0, -100000, 0, 0, 0, -100000, -100000, -100000, -100000)


# Create each plot in the original order
plot_list1 = list()
plot_list2 = list()
plot_list3 = list()
gpl1 = list()
gpl2 = list()
gpl3 = list()
mep1 = list()
mep2 = list()
mep3 = list()
widths1 = list()
widths2 = list()
widths3 = list()
widths1me = list()
widths2me = list()
widths3me = list()
plot_periods = seq(2002,2014,2)

# Get AREA PLOTS  in the original order, for both plot modes plus regular
for(i in 1:length(strata_names)){
  plot_list1[[i]] = plot_areas(tot_area_ha, plot_periods, area_ha[,i], area_lower[,i], area_upper[,i], mapped_areas[,i],
                               margin_error[,i], 0, maxy_vect1[i], strata_names[i], plotmode=1, biglabels=F)  
  plot_list2[[i]] = plot_areas(tot_area_ha, plot_periods, area_ha[,i], area_lower[,i], area_upper[,i], mapped_areas[,i],
                               margin_error[,i], miny_vect2[i], maxy_vect2[i], strata_names[i], plotmode=2, biglabels=F)  
  plot_list3[[i]] = plot_areas(tot_area_ha, plot_periods, area_ha[,i], area_lower[,i], area_upper[,i], mapped_areas[,i],
                               margin_error[,i], miny_vect2[i], maxy_vect2[i], strata_names[i], plotmode=3, biglabels=F)  
  
  gpl1[[i]] = ggplotGrob(plot_list1[[i]][[1]])
  gpl2[[i]] = ggplotGrob(plot_list2[[i]][[1]])
  gpl3[[i]] = ggplotGrob(plot_list3[[i]][[1]])
  mep1[[i]] = ggplotGrob(plot_list1[[i]][[2]])
  mep2[[i]] = ggplotGrob(plot_list2[[i]][[2]])
  mep3[[i]] = ggplotGrob(plot_list3[[i]][[2]])
  widths1[[i]] = gpl1[[i]]$widths[2:5]
  widths2[[i]] = gpl2[[i]]$widths[2:5]
  widths3[[i]] = gpl3[[i]]$widths[2:5]
  widths1me[[i]] = mep1[[i]]$widths[2:5]
  widths2me[[i]] = mep2[[i]]$widths[2:5]
  widths3me[[i]] = mep3[[i]]$widths[2:5]
}

# Calculate max width among all the grobs for each case and use that value for all of them
# This ensures the plotted areas match despite different y axis widths.
maxwidth1 = do.call(grid::unit.pmax, widths1)
maxwidth2 = do.call(grid::unit.pmax, widths2)
maxwidth3 = do.call(grid::unit.pmax, widths3)
maxwidth1me = do.call(grid::unit.pmax, widths1me)
maxwidth2me = do.call(grid::unit.pmax, widths2me)
maxwidth3me = do.call(grid::unit.pmax, widths3me)

for (i in 1:length(gpl1)){
  gpl1[[i]]$widths[2:5] = as.list(maxwidth1)
  gpl2[[i]]$widths[2:5] = as.list(maxwidth2)
  gpl3[[i]]$widths[2:5] = as.list(maxwidth3)
  mep1[[i]]$widths[2:5] = as.list(maxwidth1me)
  mep2[[i]]$widths[2:5] = as.list(maxwidth2me)
  mep3[[i]]$widths[2:5] = as.list(maxwidth3me)
}

# Arrange AREA PLOTS in the NEW grouping order and save multiplots
pontus_multiplot1 = grid.arrange(textGrob(""), gpl1[[1]], gpl1[[2]], gpl1[[4]], 
                         gpl1[[3]], gpl1[[5]], gpl1[[6]], gpl1[[7]],
                         gpl1[[8]], gpl1[[9]], gpl1[[10]], gpl1[[11]],ncol=4, 
                         left="Area [ha]", right="Percentage of total area", bottom="Time")

ggsave(paste0("plots/post_katelyn/", "ALL_Pontus1_", lut_name, ".png"), plot=pontus_multiplot1,  width = 20, height = 10) 

pontus_multiplot2 = grid.arrange(textGrob(""), gpl2[[1]], gpl2[[2]], gpl2[[4]], 
                                 gpl2[[3]], gpl2[[5]], gpl2[[6]], gpl2[[7]],
                                 gpl2[[8]], gpl2[[9]], gpl2[[10]], gpl2[[11]],ncol=4, 
                                 left="Area [ha]", right="Percentage of total area", bottom="Time")

ggsave(paste0("plots/post_katelyn/", "ALL_Pontus2_", lut_name, ".png"), plot=pontus_multiplot2,  width = 20, height = 10) 

regular_multiplot = grid.arrange(textGrob(""), gpl3[[1]], gpl3[[2]], gpl3[[4]], 
                                 gpl3[[3]], gpl3[[5]], gpl3[[6]], gpl3[[7]],
                                 gpl3[[8]], gpl3[[9]], gpl3[[10]], gpl3[[11]],ncol=4, 
                                 left="Area [ha]", right="Percentage of total area", bottom="Time")

ggsave(paste0("plots/post_katelyn/", "ALL_regular_", lut_name, ".png"), plot=regular_multiplot,  width = 20, height = 10) 


# Arrange MARGIN OF ERROR PLOTS in the NEW grouping order and save multiplots
pontus_multiplotme1 = grid.arrange(textGrob(""), mep1[[1]], mep1[[2]], mep1[[4]], 
                                 mep1[[3]], mep1[[5]], mep1[[6]], mep1[[7]],
                                 mep1[[8]], mep1[[9]], mep1[[10]], mep1[[11]],ncol=4)

ggsave(paste0("plots/post_katelyn/", "ALL_Pontus1me_", lut_name, ".png"), plot=pontus_multiplotme1,  width = 20, height = 10) 

pontus_multiplotme2 = grid.arrange(textGrob(""), mep2[[1]], mep2[[2]], mep2[[4]], 
                                 mep2[[3]], mep2[[5]], mep2[[6]], mep2[[7]],
                                 mep2[[8]], mep2[[9]], mep2[[10]], mep2[[11]],ncol=4)

ggsave(paste0("plots/post_katelyn/", "ALL_Pontus2me_", lut_name, ".png"), plot=pontus_multiplotme2,  width = 20, height = 10) 

regular_multiplot_me = grid.arrange(textGrob(""), mep3[[1]], mep3[[2]], mep3[[4]], 
                                 mep3[[3]], mep3[[5]], mep3[[6]], mep3[[7]],
                                 mep3[[8]], mep3[[9]], mep3[[10]], mep3[[11]],ncol=4)

ggsave(paste0("plots/post_katelyn/", "ALL_regular_me_", lut_name, ".png"), plot=regular_multiplot_me,  width = 20, height = 10) 



# Individual regular sized figures for separate saving with margin of error
ap = list()
mep = list()

for(i in 1:length(strata_names)){
  ap[[i]] = ggplotGrob(plot_list3[[i]][[1]])
  mep[[i]] = ggplotGrob(plot_list3[[i]][[2]])
  g = rbind(ap[[i]], mep[[i]], size="first") 
  g$widths = unit.pmax(ap[[i]]$widths, mep[[i]]$widths)
  grid.newpage()
  grid.draw(g)
  
  filename = paste0("plots/post_katelyn/", strata_names[[i]], "_areas_me_", lut_name, ".png")
  png(filename, width=1000, height = 1000, units = "px"); plot(g); dev.off()
}


###### OTHER PLOTS

# Create a single big tidy df to facilitate creating the complex plots
list_vars = list(area_ha, ci_ha, area_upper, area_lower, margin_error)
list_var_names = c("area_ha", "ci_ha", "area_upper", "area_lower", "margin_error")
melted_list = list()
plot_periods = seq(2002,2014,2)

# Assign names to classes, create year column
for (i in 1:length(list_vars)){
  names(list_vars[[i]]) = strata_names
  list_vars[[i]]$year = plot_periods
  list_vars[[i]]$varname = list_var_names[[i]]
}

melted_vars = melt(list_vars, id.vars = c('year', 'varname'), variable.name = 'class', value.name = 'value')
melted_vars = dplyr::select(melted_vars,-L1) # Drop variable created by melt_list
plot_vars = spread(melted_vars, varname, value)

# Find rows where CI or area cross the zero line
plot_vars$zero = plot_vars$area_ha < 0 | plot_vars$area_lower < 0


##### Regrowth plot
# Get only regrowth related classes and column with values != 0
regr_area = subset(plot_vars, class %in% strata_names[c(6,10,9)])

# Get cumsum of variables so that we can plot lines or points on top of the area plot
cumsumvars = regr_area %>%
  group_by(., year) %>%
  mutate(., cumsum_area_ha = cumsum(area_ha))

cumsumvars$area_lower = cumsumvars$cumsum_area_ha - cumsumvars$ci_ha
cumsumvars$area_upper = cumsumvars$cumsum_area_ha + cumsumvars$ci_ha

# Calculate the loss of regrowth based on total cumulative area, in order to
# be able to plot it. Ignore the wrong values that are created for other classes
# since we don't need them
cumsumvars$templossregr = cumsumvars$cumsum_area_ha - plot_vars$area_ha[plot_vars$class == "Loss of secondary forest"]

# Can only show the full area, there is not way to represent the "gaps" in this type of plot
regr_plot <- ggplot(regr_area, aes(x=year,y=area_ha,group=class,fill=class)) + 
  geom_area(position=position_stack(reverse = T), alpha=0.8) +
  #geom_line(data=regr_area[regr_area$class %in% c("Stable secondary forest"),], aes(x=year, y=area_lower), linetype=8) + 
  #geom_line(data=regr_area[regr_area$class %in% c("Stable secondary forest"),], aes(x=year, y=area_upper), linetype=8) +
  #geom_line(data=cumsumvars[cumsumvars$class == "Gain of secondary forest",], aes(x=year, y=templossregr), colour='red') +
  geom_point(data=cumsumvars[cumsumvars$zero == FALSE,], aes(x=year, y=cumsum_area_ha )) +
  #geom_errorbar(data=cumsumvars[cumsumvars$zero == FALSE & cumsumvars$class != "Stable secondary forest" ,], aes(x=year, ymin=area_lower, ymax=area_upper)) +
  scale_x_continuous(breaks=regr_area$year, labels = regr_area$year, minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Area [ha]") + xlab("Years")+
  scale_fill_brewer(palette="GnBu",  breaks=levels(as.factor((regr_area$class))), guide = guide_legend(reverse=T)) + 
  theme(axis.title=element_text(size=15), axis.text=element_text(size=13), legend.text=element_text(size=13),
        legend.title=element_blank()) 


print(regr_plot)
filename = paste0("plots/post_katelyn/", "secondary_forest_dynamics", ".png")
ggsave(filename, plot=regr_plot, device="png")

##### DEFOR plot
# Get only regrowth related classes and column with values != 0
defor_area = subset(plot_vars, class %in% strata_names[c(8,9)])

# Get cumsum of variables so that we can plot lines or points on top of the area plot
cumsumvars = defor_area %>%
  group_by(., year) %>%
  mutate(., cumsum_area_ha = cumsum(area_ha))

cumsumvars$area_lower = cumsumvars$cumsum_area_ha - cumsumvars$ci_ha
cumsumvars$area_upper = cumsumvars$cumsum_area_ha + cumsumvars$ci_ha

# Can only show the full area, there is not way to represent the "gaps" in this type of plot
defor_plot <- ggplot(defor_area, aes(x=year,y=area_ha,group=class,fill=class)) + 
  geom_area(position=position_stack(reverse = T), alpha=0.8) +
  #geom_line(data=defor_area[defor_area$class == "Forest to pasture",], aes(x=year, y=area_lower), linetype=8) + 
  #geom_line(data=defor_area[defor_area$class == "Forest to pasture",], aes(x=year, y=area_upper), linetype=8) +
  geom_point(data=cumsumvars[cumsumvars$zero == FALSE,], aes(x=year, y=cumsum_area_ha ), shape=3, size=3, stroke=1) +
  scale_x_continuous(breaks=defor_area$year, labels = defor_area$year, minor_breaks = NULL)  + 
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")}) + 
  ylab("Area [ha]") + xlab("Years")+
  scale_fill_brewer(palette="GnBu",  breaks=levels(as.factor((defor_area$class))), guide = guide_legend(reverse=T)) + 
  theme(axis.title=element_text(size=15), axis.text=element_text(size=13), legend.text=element_text(size=13),
        legend.title=element_blank()) 

print(defor_plot)
filename = paste0("plots/post_katelyn/", "defor_plot", ".png")
ggsave(filename, plot=defor_plot, device="png")




################### EXTRA SECTION TO AID IN SAMPLE REVISION

# Fnc to extract and put together the same column name from multiples df on a list
get_dflist_columns = function(df, colname){
  subcols = as.data.frame(do.call(rbind, lapply(df, function(x) x[,colname])))
  return(t(subcols))
}

ref1 = get_dflist_columns(cm_list, '1')
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

tdf = data.frame()
for(i in 1:length(periods)){
  print(get_condition_rows(i,1,4))
  #tdf = rbind(tdf, get_condition_rows(i, 4, 5)[,c("ID", "PTRW", "period")])
}

arrange(tdf, PTRW)


get_condition_rows(2,1,2)


# MODIFY ROWS, JUST FOR TESTING
# mc = 6
# rc = 3
# 
# for(i in 1:length(periods)){
#   match_rows = which(shp_list_ref[[i]]@data[,"STRATUM"] == mc & shp_list_ref[[i]]@data[,"ref_strata"] == rc)
#   shp_list_ref[[i]]@data[match_rows, "ref_strata"] = mc
# }



############## CREATE TABLES FOR PAPER

# Joint table of areas and stratum percentages. May not be necessary!
# Do calculations first, then assemble table.
melted_pixcount = melt(pixcount_list, id.vars = c('stratum', 'pixels'), value.name = 'value')
melted_pixcount$area_ha = melted_pixcount$pixels * 30^2 / 100^2
melted_pixcount$stratum_percentages=round(melted_pixcount$pixels / tot_area_pix * 100, digits=3) 
pixcount_df = spread(dplyr::select(melted_pixcount, -c(area_ha, stratum_percentages)), L1, pixels)
pixcount_area_ha = spread(dplyr::select(melted_pixcount, -c(pixels, stratum_percentages)), L1, area_ha)
pixcount_stratum_percentages = spread(dplyr::select(melted_pixcount, -c(area_ha, pixels)), L1, stratum_percentages)

fulldf = as.data.frame(matrix(nrow=nrow(pixcount_df), ncol=0))
for (i in 2:ncol(pixcount_df)){
  fulldf = cbind(fulldf, pixcount_area_ha[,i], pixcount_stratum_percentages[,i])
}

rownames(fulldf) = orig_strata_names 
# Need to escape special characters, including backslash itself (e.g. $\\alpha$)
colnames(fulldf) = rep(c("Area [ha]", "Area proportion [%]"), 7) 
# Create table in Latex instead, and produce the pdf there, much easier than grid.table
print(xtable(fulldf, digits=2,type = "latex",sanitize.text.function=function(x){x}))

## Table of producers accuracy
rownames(prod_acc) = periods_long
colnames(prod_acc) = strata_names
print(xtable(t(prod_acc), digits=2,type = "latex",sanitize.text.function=function(x){x}))
