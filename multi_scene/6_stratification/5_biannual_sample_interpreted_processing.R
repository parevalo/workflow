# SCRIPT TO CALCULATE AREAS AND ACCURACIES PER BIANNUAL PERIOD
# USING INDIVIDUALLY COLLECTED SAMPLES.

require(ggplot2)
require(dplyr)
require(rgdal)
require(raster)
require(rgeos)
require(grid)
require(gridExtra)
require(reshape2)
require(xtable)
require(extrafont)

loadfonts() #(device="postscript")

# Set working directories and vars
auxpath = "/media/paulo/785044BD504483BA/test/"
setwd("/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/biannual_sampling/")
strata_config = "/home/paulo/workflow/multi_scene/7_poststratification/strata_calc_config/"
funcs_path = "/home/paulo/workflow/multi_scene/R_functions/"

source(paste0(funcs_path, "area_estimation_fncs.R"))
source(paste0(funcs_path, "plotting_fncs.R"))
source(paste0(strata_config, "input_variables_original_buffer3B.R"))

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
  shp_list[[i]] = readOGR(paste0(auxpath, "biannual_samples/", samples_names[i], ".shp"), samples_names[i])

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
  errorcount = which(samples_uniqueids$count != 1050)
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

# Read basic LUT to compare strata and reference labels and run the comparison
strata_ref_lut = read.csv("/home/paulo/workflow/multi_scene/data/lut_strata_ref_match.csv",
                          as.is = c(3))

shp_list_centroid = list()
for(i in 1:length(shp_list_ref)){
  match_comparison= calc_strata(shp_list_ref[[i]]$STRATUM, 
                                 shp_list_ref[[i]]$ref_strata, 
                                 strata_ref_lut)
  # Add condition for reference = 13 and forest vs any other class except 1,8,9
  match_comparison[shp_list_ref[[i]]$STRATUM == 13] = "to_unclass"
  temp_cond = shp_list_ref[[i]]$STRATUM == 1 & !(shp_list_ref[[i]]$ref_strata %in% c(1,8,9)) 
  match_comparison[temp_cond] = "forest_misclass"
  match_comparison[match_comparison == ""] = "other_missclass"
  shp_list_ref[[i]]@data$match_comp = match_comparison
  
  # Create point version of shapefiles
  temp_centroid = gCentroid(shp_list_ref[[i]], byid=T)
  shp_list_centroid[[i]] = SpatialPointsDataFrame(temp_centroid, shp_list_ref[[i]]@data)
}

# Save updated poly and point shapefiles for analysis
# for(n in 1:length(samples_names)){
#   outname1 = paste0(samples_names[n], "_labels")
#   outname2 = paste0(samples_names[n], "_labels_pts")
#   writeOGR(shp_list_ref[[n]], paste0("shp/",outname1), outname1,
#            driver="ESRI Shapefile", overwrite_layer = T)
#   writeOGR(shp_list_centroid[[n]], paste0("shp/",outname2), outname2,
#            driver="ESRI Shapefile", overwrite_layer = T)
# }

## TEMPORARY DEFORMODE. collapse forest to pasture and forest to secondary forest
## as "forest to pasture" in ref and map labels, pixcount, mapped_areas and strata names
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
# shp_list_ref = lapply(shp_list_ref, apply_deformod)
# strata_names = strata_names[-which(strata_names == "Forest to secondary forest")]
# strata_names[which(strata_names == "Forest to pasture")] = "Deforestation"
# 
# # Recalculate mapped areas
# pixcount_list = lapply(pixcount_list, apply_deformode2)
# mapped_areas = list()
# for (i in 1:(length(periods))){
#   pixcount_list[[i]] = pixcount_list[[i]][!(pixcount_list[[i]]$stratum %in% cr),]
#   mapped_areas[[i]] = pixcount_list[[i]][!(pixcount_list[[i]]$stratum %in% cr_extra),]
# }
# 
# mapped_areas = as.data.frame(do.call(rbind, lapply(mapped_areas, '[[', 2))) * 30^2 / 100^2
  

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

# Calculate optimal sample allocation, had we used one of the confusion matrices to
# minimize the uncertainty in accuracies and areas of one of the change classes
# e.g. forest to pastures. Given just as a reference.
cm_prop_square = as.data.frame.matrix(cm_list[[1]] * as.vector(t(strata_weights[1,]/100)) / strata_pixels$x)
opt_alloc = calc_optimal_sample_alloc(cm_prop_square, 8, 1050)


# Calculate areas and accuracies, cant vectorize it the way the fncs are written
prop_out = list()
se_prop = list()
areas_out = list()
accuracies_out = list()

###### RUN WITH BUFFER
# ref_codes_all used to make sure output tables have same dimensions
for (p in 1:length(periods)){
  prop_out[[p]] = calc_props_and_vars(shp_list_ref[[p]]$STRATUM, shp_list_ref[[p]]$ref_strata, shp_list_ref[[p]]$STRATUM, 
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

area_prop = as.data.frame(do.call(rbind, lapply(prop_out, '[[', 11)))
se_simple_rnd_sampling_area_ha = sqrt(area_prop*(1-area_prop)/1050) * tot_area_ha
se_area_ha = as.data.frame(do.call(rbind, se_prop))* tot_area_ha

# Areas in kha, for article tables and figures
area_kha = area_ha / 1000
ci_kha = ci_ha /1000
area_upper_kha = area_upper / 1000
area_lower_kha = area_lower / 1000
se_area_kha = se_area_ha / 1000
se_simple_rnd_sampling_kha = se_simple_rnd_sampling_area_ha / 1000

# Create table with results for buffer class
buffer_table = as.data.frame(do.call(rbind, lapply(cm_list, function(x) x[13,])))

# Add column and row names
df_list = list(usr_acc, usr_acc_lower, usr_acc_upper,
            prod_acc, prod_acc_lower, prod_acc_upper,
            area_ha, area_upper, area_lower, 
            ci_ha, margin_error, se_area_ha, se_simple_rnd_sampling_area_ha,
            area_kha, ci_kha, area_upper_kha, area_lower_kha, 
            se_area_kha, se_simple_rnd_sampling_kha)

df_names = c("usr_acc", "usr_acc_lower", "usr_acc_upper",
             "prod_acc", "prod_acc_lower", "prod_acc_upper",
             "area_ha", "area_upper", "area_lower", 
             "ci_ha", "margin_error", "se_area_ha", "se_simple_rnd_sampling_area_ha",
             "area_kha", "ci_kha", "area_upper_kha", "area_lower_kha", 
             "se_area_kha", "se_simple_rnd_sampling_kha")

add_names = function(df, cnames, rnames){
  colnames(df) = cnames
  rownames(df) = rnames
  return(df)
}

named_df = lapply(df_list, add_names, cnames=strata_names, rnames=periods_long)
names(named_df) = df_names

# Save tables
# save_tables = function(df, savepath, f_suffix, names){
#   write.csv(df, file=paste0(savepath, names, f_suffix))
# }
# 
# suffix = paste0("_", lut_name,  "_buffered_3B.csv")
# table_savepath = "results/post_katelyn/tables/"
# 
# mapply(save_tables, named_df, savepath=table_savepath, f_suffix=suffix, names=df_names)

 
####### RUN WITHOUT BUFFER

prop_out_nb = list()
se_prop_nb = list()
areas_out_nb = list()
accuracies_out_nb = list()


# Create copies of input data, changing buffer stratum to forest and merge buffer area with forest
shp_list_ref_nb = shp_list_ref
pixcount_list_nb = pixcount_list
strata_pixels_nb = strata_pixels

for (i in 1:7){
  shp_list_ref_nb[[i]]@data[shp_list_ref_nb[[i]]$STRATUM == 16, 'STRATUM'] = 1
  buf_ind = which(pixcount_list_nb[[i]]$stratum == 16)
  for_ind = which(pixcount_list_nb[[i]]$stratum == 1)
  pixcount_list_nb[[i]][for_ind, 'pixels'] = pixcount_list_nb[[i]][for_ind, 'pixels'] + pixcount_list_nb[[i]][buf_ind, 'pixels']
  pixcount_list_nb[[i]] = pixcount_list_nb[[i]][-buf_ind,]
  
}

strata_pixels_nb[strata_pixels_nb$Group.1 == 1, 'x'] =
  strata_pixels_nb[strata_pixels_nb$Group.1 == 1, 'x'] + strata_pixels_nb[buf_ind, 'x']

strata_pixels_nb = strata_pixels_nb[-buf_ind, ]

# Run!
for (p in 1:length(periods)){
  prop_out_nb[[p]] = calc_props_and_vars(shp_list_ref_nb[[p]]$STRATUM, shp_list_ref_nb[[p]]$ref_strata, 
                                            shp_list_ref_nb[[p]]$STRATUM, 
                                            pixcount_list_nb[[p]], strata_pixels_nb, ref_codes_all)
  
  se_prop_nb[[p]] = calc_se_prop(pixcount_list_nb[[p]], strata_pixels_nb, prop_out_nb[[p]][[2]], ref_codes_all, tot_area_pix)
  areas_out_nb[[p]] = calc_unbiased_area(tot_area_pix, prop_out_nb[[p]][[11]], se_prop_nb[[p]]) 
  accuracies_out_nb[[p]] = calc_accuracies(pixcount_list_nb[[p]], strata_pixels_nb, ref_codes_all, tot_area_pix,
                                        prop_out_nb[[p]][[1]], prop_out_nb[[p]][[2]], prop_out_nb[[p]][[3]], prop_out_nb[[p]][[4]],
                                        prop_out_nb[[p]][[5]], prop_out_nb[[p]][[6]], 
                                        prop_out_nb[[p]][[7]], prop_out_nb[[p]][[8]],
                                        prop_out_nb[[p]][[9]], prop_out_nb[[p]][[10]])
  
}

# Get some key results in a readable format

overall_acc_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 1)))
overall_acc_lower_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 2)))
overall_acc_upper_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 3)))

usr_acc_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 4)))
usr_acc_lower_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 5)))
usr_acc_upper_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 6)))

prod_acc_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 7)))
prod_acc_lower_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 8)))
prod_acc_upper_nb = as.data.frame(do.call(rbind, lapply(accuracies_out_nb, '[[', 9)))

area_ha_nb = as.data.frame(do.call(rbind, lapply(areas_out_nb, '[[', 1)))
ci_ha_nb = as.data.frame(do.call(rbind, lapply(areas_out_nb, '[[', 2)))
area_upper_nb = as.data.frame(do.call(rbind, lapply(areas_out_nb, '[[', 3)))
area_lower_nb = as.data.frame(do.call(rbind, lapply(areas_out_nb, '[[', 4)))
margin_error_nb = as.data.frame(do.call(rbind, lapply(areas_out_nb, '[[', 5)))

area_prop_nb = as.data.frame(do.call(rbind, lapply(prop_out_nb, '[[', 11)))
se_area_ha_nb = as.data.frame(do.call(rbind, se_prop_nb))* tot_area_ha
se_simple_rnd_sampling_area_ha_nb = sqrt(area_prop_nb*(1-area_prop_nb)/1050) * tot_area_ha

# Areas in kha, for article tables and figures
area_kha_nb = area_ha_nb / 1000
ci_kha_nb = ci_ha_nb /1000
area_upper_kha_nb = area_upper_nb / 1000
area_lower_kha_nb = area_lower_nb / 1000
se_area_kha_nb = se_area_ha_nb / 1000
se_simple_rnd_sampling_kha_nb = se_simple_rnd_sampling_area_ha_nb / 1000

# Add column and row names
df_list_nb = list(usr_acc_nb, usr_acc_lower_nb, usr_acc_upper_nb,
               prod_acc_nb, prod_acc_lower_nb, prod_acc_upper_nb,
               area_ha_nb, area_upper_nb, area_lower_nb, 
               ci_ha_nb, margin_error_nb, se_area_ha_nb, se_simple_rnd_sampling_area_ha_nb,
               area_kha_nb, ci_kha_nb, area_upper_kha_nb, area_lower_kha_nb, 
               se_area_kha_nb, se_simple_rnd_sampling_kha_nb)

df_names_nb = c("usr_acc_nb", "usr_acc_lower_nb", "usr_acc_upper_nb",
                "prod_acc_nb", "prod_acc_lower_nb", "prod_acc_upper_nb",
                "area_ha_nb", "area_upper_nb", "area_lower_nb", 
                "ci_ha_nb", "margin_error_nb", "se_area_ha_nb", "se_simple_rnd_sampling_area_ha_nb",
                "area_kha_nb", "ci_kha_nb", "area_upper_kha_nb", "area_lower_kha_nb", 
                "se_area_kha_nb", "se_simple_rnd_sampling_kha_nb")

named_df_nb = lapply(df_list_nb, add_names, cnames=strata_names, rnames=periods_long)
names(named_df_nb) = df_names_nb

#mapply(save_tables, named_df_nb, savepath=table_savepath, f_suffix=suffix, names=df_names_nb)
# 
# write.csv(cbind(overall_acc_nb, overall_acc_lower_nb, overall_acc_upper_nb), 
#           file=paste0(table_savepath, "overall_accuracies_minmax_nb", suffix))

# Compare buffer vs no buffer
ci_compare = (ci_ha - ci_ha_nb) / ci_ha_nb * 100
colnames(ci_compare) = strata_names
rownames(ci_compare) = periods_long
  

############################# PLOT AREAS
########## AREA PLOTS USING RESULTS WITH BUFFER

# Vector of max and min y axis values for pontus modes
# Selected to guarantee that one of the breaks (6 total) is zero
maxy_vect1 = c(10, 44000000, 4500000, 90000, 4500000, 1000000, 1000000, 
               300000, 90000, 60000, 60000)
maxy_vect2 = c(12000, 45000000, 4500000, 90000, 4500000, 1000000, 1000000, 
               400000, 400000, 200000, 200000)
miny_vect1 = c(0,39000000,0,0,0,0,0,0,0,0,0)
miny_vect2 = c(-3000, 38000000, 0, -100000, 0, 0, 0, -100000, -100000, -50000, -50000)

# Limits in kha
maxy_vect1 = maxy_vect1 / 1000
miny_vect1 = miny_vect1 / 1000
maxy_vect2 = maxy_vect2 / 1000
miny_vect2 = miny_vect2 / 1000

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
#plot_periods = seq(2002,2014,2)
plot_periods = seq(2,14,2)
# Need a slightly difference letter order to create the full figure in order. 

label_letters = c("a", "b", "d", "c", "e", "f", "g", "h", "i", "j", "k")
# Other to other was removed from paper figures
label_letters2 = c("a", "a", "c", "b", "d", "e", "f", "g", "h", "i", "j")
plot_labels = mapply(paste0, label_letters, ") ", strata_names)
plot_labels2 = mapply(paste0, label_letters2, ") ", strata_names)

tot_area_kha = tot_area_ha / 1000
mapped_areas_kha = mapped_areas / 1000

# Get AREA PLOTS in the original order, for both plot modes plus regular
for(i in 1:length(strata_names)){
  plot_list1[[i]] = plot_areas(tot_area_kha, plot_periods, area_kha[,i],
                               area_lower_kha[,i], area_upper_kha[,i], mapped_areas_kha[,i],
                               margin_error[,i], miny_vect1[i], maxy_vect1[i], plot_labels2[i], plotmode=1)
  plot_list2[[i]] = plot_areas(tot_area_kha, plot_periods, area_kha[,i], 
                               area_lower_kha[,i], area_upper_kha[,i], mapped_areas_kha[,i],
                               margin_error[,i], miny_vect2[i], maxy_vect2[i], plot_labels2[i], plotmode=2)  
  plot_list3[[i]] = plot_areas(tot_area_kha, plot_periods, area_ha[,i], 
                               area_lower[,i], area_upper[,i], mapped_areas_kha[,i],
                               margin_error[,i], 0, maxy_vect2[i], strata_names[i], plotmode=3)  
  
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

gpar_settings = gpar(fontsize=10, 
                     fontfamily="Times New Roman", 
                     fontface="bold")
left_axlabel = textGrob("Area [kha]", gp=gpar_settings, rot=90)
right_axlabel = textGrob("Percentage of total area", gp=gpar_settings, rot=-90)
bottom_axlabel = textGrob("Year", gp=gpar_settings)

# Arrange AREA PLOTS in the NEW grouping order and save multiplots
pontus_multiplot1 = grid.arrange(gpl1[[2]], gpl1[[4]], 
                         gpl1[[3]], gpl1[[5]], gpl1[[6]], gpl1[[7]],
                         gpl1[[8]], gpl1[[9]], gpl1[[10]], gpl1[[11]],ncol=2, 
                         left=left_axlabel, right=right_axlabel, bottom=bottom_axlabel)

# Formatted for article
outfile=paste0("results/post_katelyn/figures/", "ALL_Pontus1_kha_larger_", lut_name, ".pdf")
ggsave(outfile, plot=pontus_multiplot1,  width = 150, height = 200, units='mm') 
embed_fonts(outfile)

pontus_multiplot2 = grid.arrange(gpl2[[2]], gpl2[[4]], 
                                 gpl2[[3]], gpl2[[5]], gpl2[[6]], gpl2[[7]],
                                 gpl2[[8]], gpl2[[9]], gpl2[[10]], gpl2[[11]], ncol=2, 
                                 left=left_axlabel, right=right_axlabel, bottom=bottom_axlabel)

outfile=paste0("results/post_katelyn/figures/", "ALL_Pontus2_kha_larger_", lut_name, ".pdf")
ggsave(outfile, plot=pontus_multiplot2,  width = 150, height = 200, units='mm') 
embed_fonts(outfile)

# Save manually created deforestation plot
deg = grid.arrange(plot_list1[[8]][[1]] + ggtitle("Deforestation"), 
                   left=left_axlabel, right=right_axlabel, bottom=bottom_axlabel)
outfile = paste0("results/post_katelyn/figures/", "Deforestation_kha_", lut_name, ".pdf")
ggsave(outfile, plot=deg,  width = 140, height = 120, units='mm') 
embed_fonts(outfile)


# Arrange MARGIN OF ERROR PLOTS in the NEW grouping order and save multiplots
left_axlabel_me = textGrob("Margin of error [%]", gp=gpar(fontsize=12, fontface="bold"), rot=90)
pontus_multiplotme1 = grid.arrange(textGrob(""), mep1[[1]], mep1[[2]], mep1[[4]], 
                                 mep1[[3]], mep1[[5]], mep1[[6]], mep1[[7]],
                                 mep1[[8]], mep1[[9]], mep1[[10]], mep1[[11]],ncol=4,
                                 left=left_axlabel_me,  bottom=bottom_axlabel)

ggsave(paste0("results/post_katelyn/figures/", "ALL_Pontus1me_kha_", lut_name, ".png"), 
       plot=pontus_multiplotme1,  width = 20, height = 10) 

pontus_multiplotme2 = grid.arrange(textGrob(""), mep2[[1]], mep2[[2]], mep2[[4]], 
                                 mep2[[3]], mep2[[5]], mep2[[6]], mep2[[7]],
                                 mep2[[8]], mep2[[9]], mep2[[10]], mep2[[11]],ncol=4,
                                 left=left_axlabel_me,  bottom=bottom_axlabel)

ggsave(paste0("results/post_katelyn/figures/", "ALL_Pontus2me_kha_", lut_name, ".png"), 
       plot=pontus_multiplotme2,  width = 20, height = 10) 

# Individual regular sized figures for separate saving with margin of error
ap = list()
mep = list()

for(i in 1:length(strata_names)){
  ap[[i]] = ggplotGrob(plot_list3[[i]][[1]])
  mep[[i]] = ggplotGrob(plot_list3[[i]][[2]])
  g = rbind(ap[[i]], mep[[i]], size="first") 
  g$widths = unit.pmax(ap[[i]]$widths, mep[[i]]$widths)
  
  filename = paste0("results/post_katelyn/figures/", strata_names[[i]], "_areas_me_", lut_name, ".png")
  ggsave(filename, plot=g, width=12, height = 15, units = "in")
}


########## AREA PLOTS USING RESULTS WITHOUT BUFFER

# Vector of max and min y axis values for pontus modes
# Selected to guarantee that one of the breaks (6 total) is zero
maxy_vect1_nb = c(300000, 45000000, 4500000, 400000, 5500000, 4500000, 4500000, 
                  2500000, 2500000, 2500000, 2500000)
maxy_vect2_nb = c(300000, 45000000, 4500000, 400000, 5500000, 4500000, 4500000, 
                  2500000, 2500000, 2500000, 2500000)
miny_vect2_nb = c(-85000, 0, 0, -100000, 0, 0, 0, 0, -100000, -100000, -100000)

# Limits in kha
maxy_vect1_nb = maxy_vect1_nb / 1000 
maxy_vect2_nb = maxy_vect2_nb / 1000
miny_vect2_nb = miny_vect2_nb / 1000


# Create each plot in the original order
plot_list1_nb = list()
plot_list2_nb = list()
plot_list3_nb = list()
gpl1_nb = list()
gpl2_nb = list()
gpl3_nb = list()
mep1_nb = list()
mep2_nb = list()
mep3_nb = list()
widths1_nb = list()
widths2_nb = list()
widths3_nb = list()
widths1me_nb = list()
widths2me_nb = list()
widths3me_nb = list()
plot_periods = seq(2002,2014,2)
plot_labels = mapply(paste0, letters[seq(1,11)], ") ", strata_names)

# Get AREA PLOTS in the original order, for both plot modes plus regular
for(i in 1:length(strata_names)){
  plot_list1_nb[[i]] = plot_areas(tot_area_kha, plot_periods, area_kha_nb[,i], 
                                  area_lower_kha_nb[,i], area_upper_kha_nb[,i], mapped_areas_kha[,i],
                                  margin_error_nb[,i], 0, maxy_vect1_nb[i], plot_labels[i], plotmode=1)  
  plot_list2_nb[[i]] = plot_areas(tot_area_kha, plot_periods, area_kha_nb[,i], 
                                  area_lower_kha_nb[,i], area_upper_kha_nb[,i], mapped_areas[,i],
                                  margin_error_nb[,i], miny_vect2[i], maxy_vect2_nb[i], strata_names[i], plotmode=2)  
  plot_list3_nb[[i]] = plot_areas(tot_area_kha, plot_periods, area_kha_nb[,i], 
                                  area_lower_kha_nb[,i], area_upper_kha_nb[,i], mapped_areas[,i],
                                  margin_error_nb[,i], 0, maxy_vect2_nb[i], strata_names[i], plotmode=3)  
  
  gpl1_nb[[i]] = ggplotGrob(plot_list1_nb[[i]][[1]])
  gpl2_nb[[i]] = ggplotGrob(plot_list2_nb[[i]][[1]])
  gpl3_nb[[i]] = ggplotGrob(plot_list3_nb[[i]][[1]])
  mep1_nb[[i]] = ggplotGrob(plot_list1_nb[[i]][[2]])
  mep2_nb[[i]] = ggplotGrob(plot_list2_nb[[i]][[2]])
  mep3_nb[[i]] = ggplotGrob(plot_list3_nb[[i]][[2]])
  widths1_nb[[i]] = gpl1_nb[[i]]$widths[2:5]
  widths2_nb[[i]] = gpl2_nb[[i]]$widths[2:5]
  widths3_nb[[i]] = gpl3_nb[[i]]$widths[2:5]
  widths1me_nb[[i]] = mep1_nb[[i]]$widths[2:5]
  widths2me_nb[[i]] = mep2_nb[[i]]$widths[2:5]
  widths3me_nb[[i]] = mep3_nb[[i]]$widths[2:5]
}

# Calculate max width among all the grobs for each case and use that value for all of them
# This ensures the plotted areas match despite different y axis widths.
maxwidth1_nb = do.call(grid::unit.pmax, widths1_nb)
maxwidth2_nb = do.call(grid::unit.pmax, widths2_nb)
maxwidth3_nb = do.call(grid::unit.pmax, widths3_nb)
maxwidth1me_nb = do.call(grid::unit.pmax, widths1me_nb)
maxwidth2me_nb = do.call(grid::unit.pmax, widths2me_nb)
maxwidth3me_nb = do.call(grid::unit.pmax, widths3me_nb)

for (i in 1:length(gpl1_nb)){
  gpl1_nb[[i]]$widths[2:5] = as.list(maxwidth1_nb)
  gpl2_nb[[i]]$widths[2:5] = as.list(maxwidth2_nb)
  gpl3_nb[[i]]$widths[2:5] = as.list(maxwidth3_nb)
  mep1_nb[[i]]$widths[2:5] = as.list(maxwidth1me_nb)
  mep2_nb[[i]]$widths[2:5] = as.list(maxwidth2me_nb)
  mep3_nb[[i]]$widths[2:5] = as.list(maxwidth3me_nb)
}

left_axlabel = textGrob("Area [kha]", gp=gpar(fontsize=12, fontface="bold"), rot=90)
right_axlabel = textGrob("Percentage of total area", gp=gpar(fontsize=12, fontface="bold"), rot=-90)
bottom_axlabel = textGrob("Year", gp=gpar(fontsize=12, fontface="bold"))

# Arrange AREA PLOTS in the NEW grouping order and save multiplots
pontus_multiplot1_nb = grid.arrange(textGrob(""), gpl1_nb[[1]], gpl1_nb[[2]], gpl1_nb[[4]], 
                                 gpl1_nb[[3]], gpl1_nb[[5]], gpl1_nb[[6]], gpl1_nb[[7]],
                                 gpl1_nb[[8]], gpl1_nb[[9]], gpl1_nb[[10]], gpl1_nb[[11]],ncol=4, 
                                 left=left_axlabel, right=right_axlabel, bottom=bottom_axlabel)

ggsave(paste0("results/post_katelyn/figures/", "ALL_Pontus1_kha_nb_", lut_name, ".png"), 
       plot=pontus_multiplot1_nb,  width = 20, height = 10, units='in') 

pontus_multiplot2_nb = grid.arrange(textGrob(""), gpl2_nb[[1]], gpl2_nb[[2]], gpl2_nb[[4]], 
                                 gpl2_nb[[3]], gpl2_nb[[5]], gpl2_nb[[6]], gpl2_nb[[7]],
                                 gpl2_nb[[8]], gpl2_nb[[9]], gpl2_nb[[10]], gpl2_nb[[11]],ncol=4, 
                                 left=left_axlabel, right=right_axlabel, bottom=bottom_axlabel)

ggsave(paste0("results/post_katelyn/figures/", "ALL_Pontus2_kha_nb_", lut_name, ".png"), 
       plot=pontus_multiplot2_nb,  width = 20, height = 10) 


# Arrange MARGIN OF ERROR PLOTS in the NEW grouping order and save multiplots
pontus_multiplotme1_nb = grid.arrange(textGrob(""), mep1_nb[[1]], mep1_nb[[2]], mep1_nb[[4]], 
                                   mep1_nb[[3]], mep1_nb[[5]], mep1_nb[[6]], mep1_nb[[7]],
                                   mep1_nb[[8]], mep1_nb[[9]], mep1_nb[[10]], mep1_nb[[11]],ncol=4)

ggsave(paste0("results/post_katelyn/figures/", "ALL_Pontus1me_kha_nb_", lut_name, ".png"), 
       plot=pontus_multiplotme1_nb,  width = 20, height = 10) 

pontus_multiplotme2_nb = grid.arrange(textGrob(""), mep2_nb[[1]], mep2_nb[[2]], mep2_nb[[4]], 
                                   mep2_nb[[3]], mep2_nb[[5]], mep2_nb[[6]], mep2_nb[[7]],
                                   mep2_nb[[8]], mep2_nb[[9]], mep2_nb[[10]], mep2_nb[[11]],ncol=4)

ggsave(paste0("results/post_katelyn/figures/", "ALL_Pontus2me_kha_nb_", lut_name, ".png"), 
       plot=pontus_multiplotme2_nb,  width = 20, height = 10) 


# Individual regular sized figures for separate saving with margin of error
ap_nb = list()
mep_nb = list()

for(i in 1:length(strata_names)){
  ap_nb[[i]] = ggplotGrob(plot_list3_nb[[i]][[1]])
  mep_nb[[i]] = ggplotGrob(plot_list3_nb[[i]][[2]])
  g = rbind(ap_nb[[i]], mep_nb[[i]], size="first") 
  g$widths = unit.pmax(ap_nb[[i]]$widths, mep_nb[[i]]$widths)
  
  filename = paste0("results/post_katelyn/figures/", strata_names[[i]], "_areas_me_nb_", lut_name, ".png")
  ggsave(filename, plot=g, width=12, height = 15, units = "in")
}

############## CREATE TABLES FOR PAPER AND PRESENTATIONS

# TABLE OF STRATA AREAS AND PROPORTIONS
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

## TABLES OF AREAS, STANDARD ERRORS, CI COMPARISON, MARGINS OF ERROR
# With buffer, biannual samples
print(xtable(t(named_df$area_kha), digits=1,type = "latex",sanitize.text.function=function(x){x}))
print(xtable(t(named_df$se_area_kha), digits=1,type = "latex",sanitize.text.function=function(x){x}))

# SE With buffer, simulating simple random sampling
print(xtable(t(named_df$se_simple_rnd_sampling_kha), digits=1,type = "latex",sanitize.text.function=function(x){x}))

# Without buffer, biannual samples
print(xtable(t(named_df_nb$area_kha), digits=1,type = "latex",sanitize.text.function=function(x){x}))
print(xtable(t(named_df_nb$se_area_kha), digits=1,type = "latex",sanitize.text.function=function(x){x}))

# Comparison of ci, buffer vs no buffer
print(xtable(t(ci_compare), digits=1,type = "latex",sanitize.text.function=function(x){x}))

# Comparison of standard errors, buffer vs no buffer
compare_se_list = list()
for (i in (1:4)){
  temp = cbind(named_df$se_area_kha[,i+7], named_df_nb$se_area_kha[,i+7])
  compare_se_list[[i]] = temp
}

compare_se = as.data.frame(do.call(cbind, lapply(compare_se_list, '[')))
rownames(compare_se) = periods_long
colnames(compare_se) = rep(c("Buffer", "No buffer"), 4)
print(xtable(format(compare_se, digits=1)))

# Margins of error
print(xtable(t(named_df$margin_error*100), digits=1,type = "latex",sanitize.text.function=function(x){x}))

## TABLE OF USERS AND PRODUCERS ACCURACY
print(xtable(t(named_df$usr_acc), digits=0,type = "latex",sanitize.text.function=function(x){x}))
print(xtable(t(named_df$prod_acc), digits=0,type = "latex",sanitize.text.function=function(x){x}))

## Test printing accuracies with confidence intervals in parenthesis, looks weird...
prod_acc_ci = named_df$prod_acc_upper - named_df$prod_acc_lower
prod_acc_ci_out = format(prod_acc_ci, digits=0)
prod_acc_out = format(named_df$prod_acc, digits=0)
prod_acc_table = mapply(paste0, prod_acc_out, " (", prod_acc_ci_out, ")")
rownames(prod_acc_table) = periods_long
colnames(prod_acc_table) = strata_names
print(xtable(t(prod_acc_table), digits=0, type = "latex",sanitize.text.function=function(x){x}))

## TABLE OF STRATA DESCRIPTION

strata_descript = c("Other transitions that are not relevant.", 
                    "Stable forest.", 
                    "Stable natural grassland.",
                    "Areas that show stable urban cover, as well as other bright surfaces like exposed rock and sand.",
                    "Stable human introduced pasturelands and croplands.", 
                    "Areas that show sustained vegetation regrowth over the course of two years or more.",
                    "Stable water bodies.", 
                    "Areas that experienced conversion from forest to pastures or croplands.",
                    "Areas that experienced a brief conversion to pastures or croplands that were abandoned shortly 
                    thereafter and display a regrowing trend.", 
                    "Areas that experienced a conversion from pastures, grasslands, urban, 
                    water and other to secondary forest.", 
                    "Areas of classes other than forest and secondary forest that experienced a disturbance 
                    but have no class label afterwards.",
                    "Areas of secondary forest that were converted to any other class (except to forest).",
                    "Areas of stable forest that were assigned to a 'buffer' stratum around change areas.")

strata_description_table = cbind(orig_strata_names[-c(11,13)], strata_descript[-c(11,13)]) 
colnames(strata_description_table) = c("Stratum name", "Description")
strata_description_table = rbind(strata_description_table[2:11,], strata_description_table[1,])
print(xtable(strata_description_table,type = "latex",sanitize.text.function=function(x){x}))

## TABLE OF STRATA DESCRIPTION INCLUDING CLASS 13 AND 16, AS WELL AS WEIGHTS AND SAMPLE ALLOCATION
# Load original stratification pixel count
ss=read.csv(paste0(auxpath, pixcount_strata), header=TRUE, col.names=c("stratum", "pixels"))
# Remove classes we don't need in the calculation (i.e. NoData)
ss = ss[!(ss$stratum %in% cr),]
tot_area_pix = sum(ss$pixels) 
str_weight = ss$pixels / tot_area_pix

strata_description_table = cbind(orig_strata_names, strata_descript, format(str_weight*100, digits=1), strata_pixels$x)
colnames(strata_description_table) = c("Stratum name", "Description", "W_h", "n_h")
strata_description_table = rbind(strata_description_table[2:13,], strata_description_table[1,])
print(xtable(as.data.frame(strata_description_table),type = "latex",sanitize.text.function=function(x){x},
             sanitize.colnames.function=bold), include.rownames=FALSE)


## INDIVIDUAL CONFUSION MATRICES FOR APPENDIX. 
orig_strata_names_short = c("Oth. to Oth.", "For.", "Grass.", "Urban", 
                            "Past.", "Sec. For.", "Wat", "For. to Past.", 
                            "For. to Sec. For", "Sec. For. Gain", "To Uncl.", "Sec. For. Loss", "Buff")

strata_names_short = c("Oth. to Oth.", "For.", "Grass.", "Urban", 
                       "Past.", "Sec. For.", "Wat", "For. to Past.", 
                       "For. to Sec. For", "Sec. For. Gain", "Sec. For. Loss")


# I'm removing TO UNCLASS and BUFFER from ref labels bc they only exist for the
# strata, and that also saves space in the table.

cm_out_colnames = c(orig_strata_names_short[-c(11,13)], "Samp, size ($n_h$)", "Strat. weight ($W_h$)")
cm_list_out = lapply(cm_list, cbind, strata_pixels$x)

cm_list_out = lapply(cm_list_out, function(x){x[,-c(11,13)]})

cm_prop_list_out = cm_list_out
col_digits1 = c(rep(0, 13), 2)
col_digits2 = c(0, rep(4, 11), 0,2)

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

for (i in 1:7){ # Couldn't get this to work with mapply!
  cm_list_out[[i]] = cbind(cm_list_out[[i]], t(strata_weights)[,i])
  colnames(cm_list_out[[i]]) = cm_out_colnames
  rownames(cm_list_out[[i]]) = orig_strata_names_short
  print(xtable(cm_list_out[[i]], digits = col_digits1, 
               caption=paste0("Confusion matrix in sample counts for period ", periods_long[i])), 
        type='latex', sanitize.text.function=function(x){x}, sanitize.colnames.function=bold,
        caption.placement="top")
  
  # Calculate proportions, ugly way
  cm_prop_list_out[[i]][,1:11] = (cm_list_out[[i]][,1:11] * cm_list_out[[i]][,13]) / cm_list_out[[i]][,12]
  cm_prop_list_out[[i]] = cbind(cm_prop_list_out[[i]], t(strata_weights)[,i])
  colnames(cm_prop_list_out[[i]]) = cm_out_colnames
  rownames(cm_prop_list_out[[i]]) = orig_strata_names_short
  print(xtable(cm_prop_list_out[[i]], digits = col_digits2,
               caption=paste0("Confusion matrix in area proportions for period ", periods_long[i])), 
        type='latex', sanitize.text.function=function(x){x}, sanitize.colnames.function=bold,
        caption.placement="top")
}

### TABLE OF STRATIFICATION 2001-2016, AREA WEIGHTS AND SAMPLES FOR REFERENCE (SLIDES)

# Calculate original strata weights and area proportions 
orig_strata_weight = (ss$pixels / tot_area_pix) *100
orig_strata_area = ss$pixels * 30^2 / 100^2
orig_strata_table = data.frame(orig_strata_names, orig_strata_area, orig_strata_weight, t(map_sample_count[1,]))
colnames(orig_strata_table) =  c("Strata names", "Area [ha]", "Area $W_h$ [\\%]", "Sample size ($n_h$)") 
orig_strata_table_out = xtable(orig_strata_table, digits=c(0,0,2,4,0), display=c("d", "s", "f", "f", "d"))
align(orig_strata_table_out) = "llrcc"
print(orig_strata_table_out,type = "latex", sanitize.text.function=function(x){x},
      sanitize.colnames.function=bold, format.args=list(big.mark = "'"), include.rownames=F)


######################## OTHER PLOTS

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
filename = paste0("results/post_katelyn/figures/", "secondary_forest_dynamics", ".png")
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
filename = paste0("results/post_katelyn/figures/", "defor_plot", ".png")
ggsave(filename, plot=defor_plot, device="png")


# RATIO of 'Forest to secondary' to 'forest to pasture'. REDO AND CHECK THE VALUES ARE CORRECT

fp = filter(defor_area, class == "Forest to pasture")
fsf = filter(defor_area, class == "Forest to secondary forest")

ratio = fsf$area_ha / fp$area_ha
ratio_table = as.data.frame(cbind(fp$year, ratio, fp$zero + fsf$zero))
colnames(ratio_table) = c("Year", "Ratio", "zero")
ggplot(ratio_table) + geom_line(aes(x=Year, y=Ratio)) + 
  geom_point(data=ratio_table[ratio_table$zero == 0,], aes(x=Year, y=Ratio), shape=3, size=3, stroke=1) 

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



