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

# Load original sample shapefiles and csvs with interpretations and make the joins
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

# Function to find duplicate ID's and return the whole row
find_duplicates = function(sample_table){
  duplicate_id = unique(sample_table$ID[duplicated(sample_table$ID)])
  duplicate_ind = which(sample_table$ID %in% duplicate_id)
  duplicate_rows = sample_table[duplicate_ind, ]

  return(duplicate_rows)
}

# Find duplicates before we join to the original table in order to avoid further
# duplication in the new table
duplicates = lapply(csv_list, find_duplicates) # 01-03 has 2 duplicates, 13-15 has 1

# Temporarily remove the duplicates in order to find the truly missing samples
remove_ind1 = as.numeric(rownames(duplicates[[1]])[c(2,4)])
csv_list[[1]] = csv_list[[1]][-remove_ind1, ]

remove_ind7 = as.numeric(rownames(duplicates[[7]])[2])
csv_list[[7]] = csv_list[[7]][-remove_ind7, ]

# Get row count per sample file
samples_nrows = as.data.frame(do.call(rbind, lapply(csv_list, nrow)))
colnames(samples_nrows) = "count"
samples_nrows$period = samples_names

# Join csv's to shapefile table in order to find MISSING SAMPLES

missing_samples = list()
for (i in 1:length(csv_list)){
  missing_samples[[i]] = anti_join(shp_list[[i]]@data, csv_list[[i]], by= "ID")
  #shp_list[[i]]@data = full_join(shp_list[[i]]@data, csv_list[[i]], by= "ID")
}


# Temporary code to check if first interpretations match between interpreters

# Read sample
testshp <- readOGR(paste0(auxpath, "biannual_samples"), "sample_01_03")
testshp_copy = testshp

# If we wanted to conver to points in order to read any raster data we needed.
# We don't need to read the strata raster bc that info is already part of the shapefile
#testpts = SpatialPointsDataFrame(gCentroid(testshp, byid=TRUE),testshp@data, match.ID=FALSE)
#writeOGR(testpts, "testpts", "test_pts", driver="ESRI Shapefile", overwrite_layer = T)

# Read reference labels
table1 = read.csv("sample_01_03_chongyang.csv")
table2 = read.csv("sample_01_03_katelyn.csv")
table3 = read.csv("sample_01_03_yihao.csv")

# Find duplicates in each table.Drop one row of each for now to make the tables 
# match the shapefile bc we need to save a new shapefile for the visual analysis.
# These duplicates won't exist in the final interpretations because I will take 
# have to check them and delete them manually first.
duplicates = unique(table1$ID[duplicated(table1$ID)])
dup_ind1 = which(table3$ID %in% duplicates)[c(2,4)]
table1 = table1[-dup_ind1,]
dup_ind2 = which(table2$ID %in% duplicates)[c(2,4)]
table2 = table2[-dup_ind2,]
dup_ind3 = which(table3$ID %in% duplicates)[c(2,4)]
table3 = table3[-dup_ind3,]

# Joint the three tables to simplify comparison
joint_table = full_join(table1, table2, by="ID")
joint_table = full_join(joint_table, table3, by="ID")

# Join codes and display disagreement between interpreters 
# Better calculate disagreement in strata!?
disag_class1 = (joint_table$CODE1 != joint_table$CODE1.x) | (joint_table$CODE1 != joint_table$CODE1.y) | (joint_table$CODE1.x != joint_table$CODE1.y)
disag_class1_ind = which(disag_class1)
imp_cols = c("ID", "PTRW.x", "PTRW.y", "PTRW","CODE1.x", "CODE1.y", "CODE1", "CODE2.x", "CODE2.y", "CODE2")
#joint_table[disag_class1_ind, imp_cols]

# Calculate strata for each interpreter
calc_strata_aux = function(code1, code2, lut){
  ref_strata = vector(length = length(code1))
  ind_na = is.na(code2)
  ref_strata[ind_na] = calc_strata(code1, code1, lut)[ind_na]
  ref_strata[!ind_na] = calc_strata(code1, code2, lut)[!ind_na]
  return(ref_strata)
}

ref_strata1 = calc_strata_aux(table1$CODE1, table1$CODE2, lut)
ref_strata2 = calc_strata_aux(table2$CODE1, table2$CODE2, lut)
ref_strata3 = calc_strata_aux(table3$CODE1, table3$CODE2, lut)

# Calculate confusion matrices
ref_codes = sort(unique(ref_strata1))
class_codes = sort(union(ref_strata1, testshp$STRATUM))

# Crosstab
table1 = inner_join(table1, testshp_copy@data, by="ID")
ct1 = calc_ct(table1$STRATUM, ref_strata1, class_codes)

table2 = inner_join(table2, testshp_copy@data, by="ID")
ct2 = calc_ct(table2$STRATUM, ref_strata2, class_codes)

table3 = inner_join(table3, testshp_copy@data, by="ID")
ct3 = calc_ct(table3$STRATUM, ref_strata3, class_codes)

# Append interpreters codes to shp
testshp@data = full_join(testshp@data, joint_table[,imp_cols], by="ID")
testshp@data = cbind(testshp@data, ref_strata1, ref_strata2, ref_strata3)

# Find pixels where there is an error according to the confusion matrix and add a column to make them easier to filter
# x=1, y=2, none=3
testshp@data[testshp$STRATUM == 1 & (!is.na(testshp$CODE2.x) & testshp$CODE2.x != 1),]

# Write shape to disk to check those samples manually   
writeOGR(testshp, "ref_01_03_review", "ref_01_03_review", driver="ESRI Shapefile", overwrite_layer = T)


# Need to create strata pix count for the buffered biannual samples
ss = read.csv(paste0(auxpath,"strata_buffered_01_03_pixcount.csv"), header=TRUE, col.names=c("stratum", "pixels"))
ss = ss[!(ss$stratum %in% cr),] 
tot_area_pix = sum(ss$pixels)
strata_pixels = aggregate(testshp$STRATUM, by=list(testshp$STRATUM), length)

# Accuracies
prop_out1 = calc_props_and_vars(table1$STRATUM, ref_strata1, table1$STRATUM, 
                               ss, strata_pixels, ref_codes) 

se_prop1 = calc_se_prop(ss, strata_pixels, prop_out1[[2]], ref_codes, tot_area_pix)
areas_out1 = calc_unbiased_area(tot_area_pix, prop_out1[[11]], se_prop1) 
accuracies_out1 = calc_accuracies(ss, strata_pixels, ref_codes, tot_area_pix,
                                 prop_out1[[1]], prop_out1[[2]], prop_out1[[3]], prop_out1[[4]],
                                 prop_out1[[5]], prop_out1[[6]], 
                                 prop_out1[[7]], prop_out1[[8]],
                                 prop_out1[[9]], prop_out1[[10]])

prop_out2 = calc_props_and_vars(table2$STRATUM, ref_strata2, table2$STRATUM, 
                                ss, strata_pixels, ref_codes) 

se_prop2 = calc_se_prop(ss, strata_pixels, prop_out2[[2]], ref_codes, tot_area_pix)
areas_out2 = calc_unbiased_area(tot_area_pix, prop_out2[[11]], se_prop2) 
accuracies_out2 = calc_accuracies(ss, strata_pixels, ref_codes, tot_area_pix,
                                  prop_out2[[1]], prop_out2[[2]], prop_out2[[3]], prop_out2[[4]],
                                  prop_out2[[5]], prop_out2[[6]], 
                                  prop_out2[[7]], prop_out2[[8]],
                                  prop_out2[[9]], prop_out2[[10]])

prop_out3 = calc_props_and_vars(table3$STRATUM, ref_strata3, table3$STRATUM, 
                                ss, strata_pixels, ref_codes) 

se_prop3 = calc_se_prop(ss, strata_pixels, prop_out3[[2]], ref_codes, tot_area_pix)
areas_out3 = calc_unbiased_area(tot_area_pix, prop_out3[[11]], se_prop3) 
accuracies_out3 = calc_accuracies(ss, strata_pixels, ref_codes, tot_area_pix,
                                  prop_out3[[1]], prop_out3[[2]], prop_out3[[3]], prop_out3[[4]],
                                  prop_out3[[5]], prop_out3[[6]], 
                                  prop_out3[[7]], prop_out3[[8]],
                                  prop_out3[[9]], prop_out3[[10]])

