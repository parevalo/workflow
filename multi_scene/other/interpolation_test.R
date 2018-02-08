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
    train = input_data[sel,]
    test = input_data[!sel,] 
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
#' Returns: list of partitions, each with a boolean vector to use to subset the 
#' input data
create_partitions = function(group_vect, id_vect, npart,fraction){
  class_folds = list()
  df = as.data.frame(cbind(group_vect, id_vect))
  for (k in 1:npart){
    sfrac = df %>% group_by(group_vect) %>% sample_frac(., size=fraction)
    class_folds[[k]] = df$id_vect %in% sfrac$id_vect
  }
  return(class_folds)
}

#' Iterate over classes and number of neighbors in a pseudo CV fashion to obtain
#' the NN we will use.
#' If there is only one ocurrence of a value, we need the sample_frac to be >=0.6 
#' so that the single occurrence gets included in every sample we take. 

nsamples = 10
nb = seq(from = 5, to = 30 , by = 1)
sample_frac = 0.6
idp = 1.5
cv_results = list()
class_results = list()
shp_results = list()
crs_obj = CRS("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
set.seed(1512)

# Run for a single period-class: (e.g. 09-11, #5)

# class_results = list()
# for (i in 1:length(ref_class_list[[5]])){
#   print(paste0("Iterating with map class ", ref_class_list[[5]][[i]]))
#   subset_ds = shp_list[[5]][shp_list[[5]]$map_strata == ref_class_list[[5]][[i]],]
#   test_part = create_partitions(subset_ds@data$binmatch, subset_ds@data$ID, 10, 0.6)
#   
#   # Run with multiple nearest neighbors values
#   cv_results = mclapply(nb, run_cv, fold_list=test_part, input_data=subset_ds, idp_val=idp)
# 
#   names(cv_results) = sprintf("%02d",nb)
#   class_results[[i]] = cv_results
# }

# Fnc to get CV results per single period, all classes
get_class_cv_results = function(class_list, full_shp, npart, fraction, nb_seq, idp_val){
  out_class_results = list()
  for (i in 1:length(class_list)){
    print(paste0("Iterating with map class ", class_list[[i]]))
    subset_ds = full_shp[full_shp$map_strata == class_list[[i]],]
    test_part = create_partitions(subset_ds@data$binmatch, subset_ds@data$ID, npart, fraction)
    
    # Run with multiple nearest neighbors values
    cv_results = mclapply(nb_seq, run_cv, fold_list=test_part, input_data=subset_ds, idp_val=idp_val)
    
    names(cv_results) = sprintf("%02d",nb_seq)
    out_class_results[[i]] = cv_results
  }
  return(out_class_results)
}

class_results1 = get_class_cv_results(ref_class_list[[5]][[4]], shp_list[[5]], nsamples, sample_frac, nb, idp)

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
  auc_boxplot[[i]] <- ggplot(class_auc[[i]], aes(x = neigh, y = auc)) + geom_boxplot() + 
    ylab("AUCROC") + xlab("# Neighbors") + ggtitle(paste0("idp=", idp))
  ggsave(paste0("interpolation/auc_vs_nn_class", ref_class_list[[5]][[i]], "_idp", idp, ".png"),
         plot=auc_boxplot[[i]], width = 10, height = 5, units='in') 
  
  # Select highest mean AUC
  mean_auc[[i]] = class_auc[[i]] %>% group_by(neigh) %>% summarise(avg=mean(auc))
  # Rewrite to make it easier to read. But the idea is that if two numbers meet 
  # the criteria, take the smaller neighbor
  max_nn[[i]] = as.numeric(min(mean_auc[[i]]$neigh[mean_auc[[i]]$avg == max(mean_auc[[i]]$avg)]))
}

# Fnc to calculate the # of neighbors per class required to do the interpolation
calc_max_nn = function(class_list, class_results, idp_val){
  class_auc = list()
  auc_boxplot = list()
  mean_auc = list()
  max_nn = list()
  for (i in 1:length(class_list)){
    
    # Extract AUCs per partition per neighbor
    auc_list_nb = as.data.frame(do.call(cbind, lapply(class_results[[i]], function(x){x[[2]]})))
    class_auc[[i]] = gather(auc_list_nb, neigh, auc)
    auc_boxplot[[i]] <- ggplot(class_auc[[i]], aes(x = neigh, y = auc)) + geom_boxplot() + 
      ylab("AUCROC") + xlab("# Neighbors") + ggtitle(paste0("idp=", idp_val))
    ggsave(paste0("interpolation/auc_vs_nn_class", class_list[[i]], "_idp", idp_val, ".png"),
           plot=auc_boxplot[[i]], width = 10, height = 5, units='in') 
    
    # Select highest mean AUC
    mean_auc[[i]] = class_auc[[i]] %>% group_by(neigh) %>% summarise(avg=mean(auc))
    # Rewrite to make it easier to read. But the idea is that if two numbers meet 
    # the criteria, take the smaller neighbor
    max_nn[[i]] = as.numeric(min(mean_auc[[i]]$neigh[mean_auc[[i]]$avg == max(mean_auc[[i]]$avg)]))
  }
  return(max_nn)
}

max_nn_idp15 = calc_max_nn(ref_class_list[[5]], class_results1, idp)

# Load rasters where interpolation will be performed
testraster = raster("/media/paulo/785044BD504483BA/test/final_strata_annual_09_11_UTM18N.tif")
# Crop to be able to test locally
cropped_raster = crop(testraster, y=extent(testraster, 10000, 14000, 10000, 14000))
NAvalue(cropped_raster) = 15

#' Function to run interpolation over a given raster for a list of classes, 
#' given a training ds and idp.
interpolate_fnc = function(cat_raster, train_shp, class_list, nn, idp_val, nsuffix){
  # Needed to convert to SpatialPixels
  crs_string ="+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  for (i in 1:length(class_list)){
    iraster = cat_raster
    # Set other classes as NA to interpolate on that class only
    iraster[iraster < class_list[[i]] | iraster > class_list[[i]]] = NA
    # Use points inside the class we need
    subset_shp = train_shp[train_shp$map_strata == class_list[[i]],]
    
    # Convert raster to spatialpixels, required for IDW
    area_sp = SpatialPixels(SpatialPoints(coordinates(iraster)[!is.na(values(iraster)),]),
                            proj4string = CRS(crs_string))
    
    # Convert to spatialpixels and run the interpolation
    print(paste0("Running interpolation for class ", class_list[[i]]))
    class_idw = idw(subset_shp@data$binmatch~1, locations=subset_shp, 
                    newdata=area_sp, idp=idp_val, nmax=nn[[i]])
    class_interp = raster(class_idw["var1.pred"])
    
    rast_name = paste0("interpolation/IDW_class_", class_list[[i]], "_nn", max_nn[[i]], 
                       "_idp", idp_val, nsuffix)
    
    writeRaster(class_interp, rast_name, format="GTiff", overwrite=TRUE)
  }
}

# In order to test the performace of the definitive interpolation, set aside
# training and testing datasets. 
train_part = list()
train_ds = list()
test_ds = list()
id_list = list()

# Get fractional samples per class
for (i in 1:length(ref_class_list[[5]])){
  classds = shp_list[[5]][shp_list[[5]]$map_strata == ref_class_list[[5]][[i]],]
  # Ensure we can get the same samples every time
  set.seed(1515)
  train_part[i] = create_partitions(classds@data$binmatch, classds@data$ID, 1, 0.6)
  # get ID of the rows and then subset them all at once 
  id_list[[i]] = classds$ID[train_part[[i]]]
}

# Unlist and find boolean indices for train and test datasets. 
train_ids = unlist(id_list)
train_sel = shp_list[[5]]$ID %in% train_ids
train_ds = shp_list[[5]][train_sel,] 
test_ds = shp_list[[5]][!train_sel,] 


# Run with ALL points per class
interpolate_fnc(cropped_raster, shp_list[[5]], ref_class_list[[5]], max_nn, idp, "_ALL_09-11.tif")

# Run with TRAIN points per class only
interpolate_fnc(cropped_raster, train_ds, ref_class_list[[5]], max_nn, idp, "_train_09-11.tif")

# Load interpolated rasters to calculate AUCROC
testlist = list()
int_rasts = list()

for (i in 1:length(ref_class_list[[5]])){
  rast_name = paste0("interpolation/IDW_class_", ref_class_list[[5]][[i]], "_nn", max_nn[[i]], 
                     "_idp", idp, "_train_09-11.tif")
  int_rasts[[i]] = raster(rast_name)
}

# Merge all rasters
interp_raster = do.call(raster::merge, int_rasts)

# Extract values
# Need to convert to spatialpoints for extract to work correctly
sp_temp = SpatialPointsDataFrame(test_ds, proj4string =crs_obj, data=test_ds@data)
sp_temp = raster::extract(interp_raster, sp_temp, sp = TRUE) 

# Calculate ROC for final maps using the TEST dataset
interp_roc = roc(sp_temp$binmatch, sp_temp$layer)
interp_roc$auc
