require(rgdal)
require(raster)
library(sp)
library(gstat)
library(caret)
library(pROC)
library(tidyverse)


setwd("/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/biannual_sampling/")

start = 2001
end = 2016
step = 2
years = seq(start, end, step) 
short_years = substr(years, 3,4) # Get years in two digit format
periods = paste0(short_years[-8], "_", short_years[-1])
periods_long = paste0(years[-8], "-", years[-1])
samples_names = character()
shp_list = list()
ref_class_list = list()

# Read shapefiles with reference and map labels
for (i in 1:(length(periods))){
  samples_names[i] = paste0("sample_", periods[i], "_labels")
  # Load shapes with map strata
  shp_list[[i]] = readOGR(paste0("shp/",  samples_names[i], "/", samples_names[i], ".shp"), samples_names[i])
}

# Create new column indicating perfect label match or not and get classes per shp
# Remove class 0 bc it creates problems for the ROC calculation
for (i in 1:length(shp_list)){
  shp_list[[i]]@data$map_strata = shp_list[[i]]@data$STRATUM 
  shp_list[[i]]@data$map_strata[shp_list[[i]]@data$map_strata == 16] = 1
  shp_list[[i]]@data$binmatch = (shp_list[[i]]@data$map_strata == shp_list[[i]]@data$ref_strata)*1
  cl_temp = sort(unique(shp_list[[i]]$ref_strata))
  ref_class_list[[i]] = cl_temp[!(cl_temp %in% 0)]
}

#' Fnc to calculate the AUC-ROC of multiple folds used to run IDW interpolation.
#' Done manually bc krig.cv throws an error
#' @param fold_list: List with n elements, each with a sequence of indexes 
#' to select the rows that belong to that fold.
#' @input_data: SpatialPoints or SpatialPolygons dataframe
#' @nn: Number of neares neighbords to run the interpolation over
#' @idp_val = Inverse Distance Power value

run_cv <- function(fold_list, input_data, nn, idp_val){
  fold_pred = list()
  roc_list =list()
  auc_list = list()
  
  for (i in 1:length(fold_list)) {
    sel = fold_list[[i]]
    train = input_data[-sel,]
    test = input_data[sel,]
    gstat_obj = gstat(id = "binmatch", formula = binmatch~1, 
                      data = train, nmax = nn, 
                      set = list(idp = idp_val))
    fold_pred[[i]] = predict(gstat_obj, test)
    roc_list[[i]] = roc(test$binmatch, fold_pred[[i]]$binmatch.pred)
    auc_list[[i]] = auc(roc_list[[i]])
  }
  
  return(list(fold_pred, auc=unlist(auc_list)))
}

#' Fnc to create stratified random "folds" with proportional sampling
#' A proper stratified kfold cannot be done when there is only a few cases of one class. 
#' Therefore something similar to a stratified bootstrap is required (but without
#' resampling because the IDW doesn't take repeated values. Here it is done
#' doing proportional sampling
create_partitions = function(group_vect, id_vect, npart,fraction){
  class_folds = list()
  for (k in 1:npart){
    df = as.data.frame(cbind(group_vect, id_vect))
    sfrac = df %>% group_by(group_vect) %>% sample_frac(., size=fraction)
    class_folds[[k]] = which(df$id_vect %in% sfrac$id_vect)
  }
  return(class_folds)
}

#' Iterate over ALL shapefiles, classes and number of neighbors
#' If there is only one ocurrence of a value, we need the sample_frac to be >=0.6 
#' so that the single occurrence gets included in every sample we take. 

nsamples = 10
nb = seq(from = 5, to = 6 , by = 1)
sample_frac = 0.6
cv_results = list()
class_results = list()
shp_results = list()
class_folds = list()
set.seed(1512)

for (i in 1:length(shp_list)){
  for (j in 1:length(ref_class_list[[i]])){
    class_ds = shp_list[[i]][shp_list[[i]]$map_strata == ref_class_list[[i]][[j]],]
    
    # Create stratified "folds". See above for explanation.
    class_folds = create_partitions(class_ds@data$binmatch, class_ds@data$ID, nsamples, sample_frac)
    
    # Iterate over multiple neighbors to find the one with the smallest AUC
    for (k in 1:length(nb)){
      cv_results[[k]] = run_cv(class_folds, class_ds, nb[[k]])
    }
    
    names(cv_results) = sprintf("%02d",nb)
    class_results[[j]] = cv_results
  }
  shp_results[[i]] = class_results
  names(shp_results[[i]]) = ref_class_list[[i]]
}


# Or run for a single period-class instead: (COMPLETE!)
subset_ds = shp_list[[1]][shp_list[[1]]$map_strata == ref_class_list[[i]][[j]],]
testpart = create_partitions(subset_ds@data$binmatch, subset_ds@data$ID, 10, 0.6)

# Format and reshape to make it easier to plot
auc_list = as.data.frame(do.call(cbind, lapply(shp_results[[1]]$`8`, function(x){x[[2]]})))
auc_tidy = gather(auc_list, neigh, auc)
auc_boxp <- ggplot(auc_tidy, aes(x = neigh, y = auc)) + geom_boxplot()
auc_boxp

# Select highest mean AUC
mean_auc = auc_tidy %>% group_by(neigh) %>% summarise(avg=mean(auc))
max_nn = as.numeric(mean_auc$neigh[mean_auc$avg == max(mean_auc$avg)])

# Load rasters where interpolation will be performed
testraster = raster("/media/paulo/785044BD504483BA/test/final_strata_annual_09_11_UTM18N.tif")
# Crop to be able to test locally
cropped_raster = crop(testraster, y=extent(testraster, 10000, 14000, 10000, 14000))
cropped_raster[cropped_raster == 15] = NA # To set actual NA's
cropped_raster[cropped_raster < 1 | cropped_raster > 1 ] = NA # Or to interpolate for a single class only

# Convert to spatialpixels, required for idw
crs_string ="+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
area_sp = SpatialPixels(SpatialPoints(coordinates(cropped_raster)[!is.na(values(cropped_raster)),]),
                        proj4string = CRS(crs_string))

# Convert to spatialpixels and run the interpolation

testidw = idw(shp_list[[5]]@data$binmatch~1, locations=shp_list[[5]], 
              newdata=area_sp, idp=1, nmax=max_nn)
testplot = raster(testidw["var1.pred"])
writeRaster(testplot, "testraster_forest_09-11.tif", format="GTiff", overwrite=TRUE)







# ###################### TEST DATA
# data(meuse)
# data("meuse.grid")
# 
# coordinates(meuse) = ~x+y
# coordinates(meuse.grid) = ~x+y
# 
# # Convert from Spatialpoints to Spatialpixels
# gridded(meuse.grid) = TRUE
# 
# # Way 1
# idw2 = idw(formula=zinc~1, locations=meuse, newdata=meuse.grid, idp=1, nmax=7)
# spplot(idw2["var1.pred"])
# 
# 
# # Way 2
# meuse.gstat <- gstat(id = "lime", formula = lime~1,
#                      data = meuse, nmax = 7,
#                      set = list(idp = 1))
# z <- predict(meuse.gstat, meuse.grid)
# spplot(z["lime.pred"])
# 
# # Cross-validation to calculate ROCAUC
# # Need to iterate to select multiple NN values
# 
# set.seed(1010)
# meuse.gstat.cv = gstat.cv(meuse.gstat, nfold=10)
# #Or from scratch without referring to an object
# meuse.idw_cv <- krige.cv(zinc~1, meuse, nmax = 7, nfold=5, set = list(idp = 1))
# 
# # Using variable numbers of NN
# # Nmax has to be different for each class
# nb = seq(from = 3, to = 30 , by = 1)
# results_cv=list()
# 
# for(n in nb){
#   # Saving only the CV results
#   results_cv[[n]] = krige.cv(lime~1, meuse, nmax = n, set = list(idp = 1))
# }
# 
# 
# #ROCAUC, values need to be binary
# category <- c(1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
# prediction <- rev(seq_along(category))
# prediction[9:10] <- mean(prediction[9:10])
# roc_obj <- roc(category, prediction)
# auc_test = auc(roc_obj)
# 
# 
# #####################


