# Script to join the tables of interpreted samples to the original shapefiles

require(tidyverse)
require(rgdal)
require(raster)
require(rgeos)

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

for (y in 1:(length(years)-1)){
  fname = paste0("sample_", short_years[y], "_", short_years[y+1])

  # Load shapes
  samples_names[y] = fname
  shp_list[[y]] = readOGR(paste0(auxpath, "biannual_samples/"), samples_names[y])

  #Load csvs
  csv_list[[y]] <- read.csv(paste0(samples_names[y], ".csv"))

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

# Crosstab map and ref labels for each period
create_cm = function(shape){
  ref_codes = sort(unique(shape@data$ref_strata))
  map_codes = sort(unique(shape@data$STRATUM))
  class_codes = sort(union(ref_codes, map_codes))
  cm = calc_ct(shape@data$STRATUM, shape@data$ref_strata, class_codes)
  return(cm)
}

cm_list = lapply(shp_list_ref, create_cm)





