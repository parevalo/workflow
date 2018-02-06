require(rgdal)
require(raster)
library(sp)
library(gstat)
library(caret)
library(pROC)
library(tidyverse)
library(parallel)


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
    print(paste0("Iterating with fold ", i, " out of ", length(fold_list)))
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
  df = as.data.frame(cbind(group_vect, id_vect))
  for (k in 1:npart){
    sfrac = df %>% group_by(group_vect) %>% sample_frac(., size=fraction)
    class_folds[[k]] = which(df$id_vect %in% sfrac$id_vect)
  }
  return(class_folds)
}

#' Iterate over ALL shapefiles, classes and number of neighbors
#' If there is only one ocurrence of a value, we need the sample_frac to be >=0.6 
#' so that the single occurrence gets included in every sample we take. 

nsamples = 10
nb = seq(from = 5, to = 30 , by = 1)
sample_frac = 0.6
idp = 0.5
cv_results = list()
class_results = list()
shp_results = list()
set.seed(1512)

for (i in 1:length(shp_list)){
  for (j in 1:length(ref_class_list[[i]])){
    class_ds = shp_list[[i]][shp_list[[i]]$map_strata == ref_class_list[[i]][[j]],]
    
    # Create stratified "folds". See above for explanation.
    class_folds = create_partitions(class_ds@data$binmatch, class_ds@data$ID, nsamples, sample_frac)
    
    # Iterate over multiple neighbors to find the one with the smallest AUC
    for (k in 1:length(nb)){
      cv_results[[k]] = run_cv(class_folds, class_ds, nb[[k]], idp)
    }
    
    names(cv_results) = sprintf("%02d",nb)
    class_results[[j]] = cv_results
  }
  shp_results[[i]] = class_results
  names(shp_results[[i]]) = ref_class_list[[i]]
}


# Or run for a single period-class instead: (e.g. 09-11)

class_results = list()
for (i in 1:length(ref_class_list[[5]])){
  print(paste0("Iterating with map class ", ref_class_list[[5]][[i]]))
  subset_ds = shp_list[[5]][shp_list[[5]]$map_strata == ref_class_list[[5]][[i]],]
  test_part = create_partitions(subset_ds@data$binmatch, subset_ds@data$ID, 10, 0.6)
  
  # Iterate over multiple neighbors to find the one with the smallest AUC
  # for (k in 1:length(nb)){
  #   print(paste0("Iterating with map class ", ref_class_list[[5]][[i]], ", neighbor ", nb[[k]]))
  #   cv_results[[k]] = run_cv(test_part, subset_ds, nb[[k]], idp)
  # }
  
  # Or better use parallel apply to speed up!
  cv_results = mclapply(nb, run_cv, fold_list=test_part, input_data=subset_ds, idp_val=idp)
  
  names(cv_results) = sprintf("%02d",nb)
  class_results[[i]] = cv_results
}

# Store AUC results per class/partition/neighbor in lists and create boxplots with them.
# Also select highest AUC automatically to input into the interpolation
# Done here for period 09-11

class_auc = list()
auc_boxplot = list()
mean_auc = list()
max_nn = list()
for (i in 1:length(ref_class_list[[5]])){
  
  # Extract AUCs per partition per neighbor
  auc_list_nb = as.data.frame(do.call(cbind, lapply(class_results[[i]], function(x){x[[2]]})))
  class_auc[[i]] = gather(auc_list_nb, neigh, auc)
  auc_boxplot[[i]] <- ggplot(class_auc[[i]], aes(x = neigh, y = auc)) + geom_boxplot()
  
  # Select highest mean AUC
  mean_auc[[i]] = class_auc[[i]] %>% group_by(neigh) %>% summarise(avg=mean(auc))
  max_nn[[i]] = as.numeric(mean_auc[[i]]$neigh[mean_auc[[i]]$avg == max(mean_auc[[i]]$avg)])
}
  

# Load rasters where interpolation will be performed
testraster = raster("/media/paulo/785044BD504483BA/test/final_strata_annual_09_11_UTM18N.tif")
# Crop to be able to test locally
cropped_raster = crop(testraster, y=extent(testraster, 10000, 14000, 10000, 14000))
cropped_raster[cropped_raster == 15] = NA # To set actual NA's

crs_string ="+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

for (i in 1:length(ref_class_list[[5]])){
  iraster = cropped_raster
  # Set other classes as NA to interpolate on that class only
  iraster[iraster < ref_class_list[[5]][[i]] | iraster > ref_class_list[[5]][[i]]] = NA
  
  # Convert to spatialpixels, required for IDW
  area_sp = SpatialPixels(SpatialPoints(coordinates(iraster)[!is.na(values(iraster)),]),
                          proj4string = CRS(crs_string))
  
  # Convert to spatialpixels and run the interpolation
  print(paste0("Running interpolation for class ", ref_class_list[[5]][[i]]))
  class_idw = idw(shp_list[[5]]@data$binmatch~1, locations=shp_list[[5]], 
                newdata=area_sp, idp=idp, nmax=max_nn[[i]])
  class_interp = raster(class_idw["var1.pred"])
  rast_name = paste0("IDW_class_", ref_class_list[[5]][[i]], "_nn", max_nn[[i]], 
                     "_idp", idp, "_09-11.tif")
  writeRaster(class_interp, rast_name, format="GTiff", overwrite=TRUE)
}


