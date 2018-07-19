#' THIS SCRIPT CONTAINS ALL THE FUNCTIONS REQUIRED FOR THE CALCULATION OF
#' THE UNBIASED AREAS WHEN THE SAMPLES WERE COLLECTED FROM A DIFFERENT
#'  MAP THAN THE ONE BEING ASSESSED:

#' 1) Function to reclassify values of a vector using a lookup table
#' 2) Function to calculate strata for any given pair of years
#' 3) Calculate unbiased area proportions and variance of reference classes
#' 4) Calculate standard error of unbiased area proportion
#' 5) Calculate unbiased areas, their confidence interval and margin of error

#' Function to reclassify any given vector using a lookup table. 
#' 
#' @param vector Vector with integer values to be reclassified.
#' @param lut matrix or dataframe with two columns: the first holding the values
#' to search from, and the second with the replacement values
#' @return Vector of equal lenght than the original with the replaced values. If
#' values in the vector are not found in the codes to be searched for, NA is
#' returned.
reclass_codes <- function(vector, lut){
  temp = vector(mode = "integer", length = length(vector))
  temp[!(vector %in% lut[,1])] = NA
  for(i in 1:dim(lut)[1]){
    temp[vector == lut[i,1]] = lut[i,2] 
  }
  return(temp)
}


#' Calculate new code from element-wise comparison of two vectors 
#' using a lookup table
#'
#' @param v1 First vector of values
#' @param v2 Second vector of values
#' @param lut matrix or dataframe with three columns: the first two holding the 
#' first and second values that v1 and v2 will be compared
#' against, respectively, and the third with the new code to be assigned.
#' @return Vector of equal lenght than the original with the new values. If
#' values in the vector are not found in the codes to be searched for, NA is
#' returned.
#' 

calc_strata <- function(v1, v2, lut){
  if (length(v1) != length(v2)){
    stop("Lengths of the two years are different")
  } else {
    # Get data type from first element of the new values to be assigned
    type_of = typeof(lut[1,3])
    temp = vector(mode=type_of, length = length(v1))
    temp[!(v1 %in% lut[,1]) | !(v2 %in% lut[,2])] = NA
    for(i in 1:dim(lut)[1]){
      temp[(v1 == lut[i,1]) & (v2 == lut[i,2])] = lut[i,3] 
    }
  }
  return(temp)
}


#' Function to calculate the strata for any given pair of years (integers). 
#' If an invalid class code is provided, NA is returned

calculate_strata_old <- function(year1, year2){
  strata = 0
  class_list = c(0,1,2,3,4,5,6,7)
  if (!(year1 %in% class_list) | !(year2 %in% class_list)) {
    strata=NA
  }
  
  else{
    if (year1 == year2) {
      strata = year1 }
    if ((year1 == 5) & (year2 == 1)) {
      strata = 5 }
    if ((year1 == 1) & (year2 == 4)) {
      strata = 8 }
    if ((year1 == 1) & (year2 == 5)) {
      strata = 9 }
    if ((year1 == 1) & !(year2 %in% c(1, 4, 5))) {
      strata = 0 }
    if ((year1 == 5) & !(year2 %in% c(1, 5))) {
      strata = 14 }
    if ((year2 == 5) & !(year1 %in% c(1, 5))) {
      strata = 11 }
    if ((year1 == 7) & (year2 == 7)){
      strata = 3 }
  }  
  return(strata)
}


#' Calculate proportions and variances for area and accuracies calculation per 
#' original strata for a given year/map.
#' 
#' @param orig_strata Vector with numeric codes representing the original 
#' stratification of each sample
#' @param ref_label Vector with numeric codes representing the reference label 
#' for that year/map
#' @param map_label Vector with numeric codes representing the map labels 
#' @param strata_totals Dataframe with two columns and number of rows equal to 
#' the total number of classes in the original strata. The first column must 
#' have the same codes found in the original stratification and the second must 
#' have the total number of PIXELS of each class in that original strata map.
#' @param sample_totals Dataframe with two columns and number of rows equal to 
#' the total number of classes in the original strata. The first column must 
#' have the same codes found in the original stratification, and the second must
#' have the total number of SAMPLE UNITS of each class collected from that 
#' original strata map. 
#' @param rfcodes Vector with numeric values representing the reference codes
#' present in ALL of the periods. 
#' @return List with: dataframe with proportion (mean) of reference labels 
#' present on each sample strata class (ref_prop) and dataframe with its 
#' variance (ref_var), dataframe with proportion (mean) of map labels present 
#' on each sample strata class (map_prop), and dataframe with its variance 
#' (map_var), dataframe with proportion (mean) of matching reference and map 
#' labels present on each sample strata class (map_and_ref_prop) and dataframe 
#' with its variance (map_and_ref_var), dataframe with proportion of all map 
#' labels equal to reference labels present in stratum (overall_acc_prop) and 
#' dataframe with its variance (overall_acc_var) and dataframes associated with 
#  users's and producer's covariance (users_cov, producers_cov), vector of area 
#' proportions per class (class_prop).
#' @export

calc_props_and_vars = function(orig_strata, ref_label, map_label, 
                               strata_totals, sample_totals, rfcodes){
  
  # Obtain unique values in the orig_strata field and get a sequence. We will 
  # use ALL codes even if they are not present for a year, in which case we 
  # would obtain values of 0. This facilitates all the other calculations
  str_codes = sort(unique(orig_strata))
  str_seq = seq_along(str_codes)
  
  # Get a sequence for the reference codes
  ref_seq = seq_along(rfcodes)
  
  # Initialize empty df for proportions per orig_strata, and for sample variance
  # per orig_strata
  ref_prop = data.frame()
  ref_var = data.frame()
  map_prop = data.frame()
  map_var = data.frame()
  map_and_ref_prop = data.frame()
  map_and_ref_var = data.frame()
  overall_acc_prop = vector()
  overall_acc_var = vector()
  users_cov = data.frame()
  producers_cov = data.frame()
  
  # Compare the fields, iterate over "orig_strata" and "ref_label" classes
  for (s in str_seq){
    
    # Get location of rows for the current orig_strata
    str_bool = orig_strata == str_codes[s]
    str_ind = which(str_bool)
    
    for (r in ref_seq){
      # Compare reference vs stratum and get TRUE or FALSE on each row
      ref_vs_str_bool = str_bool & ref_label == rfcodes[r]
      # Get row numbers that meet that condition
      ref_vs_str_ind = which(ref_vs_str_bool)
      
      # Compare map vs stratum and get TRUE or FALSE on each row
      map_vs_str_bool = str_bool & map_label == rfcodes[r]
      # Get row numbers that meet that condition
      map_vs_str_ind = which(map_vs_str_bool)
      
      # Get places where the current reference and map are the same in the 
      # current stratum class
      ref_vs_map_bool = ref_vs_str_bool & map_vs_str_bool
      # Get row numbers that meet that condition
      ref_vs_map_ind = which(ref_vs_map_bool)
      
      # Get place where ANY reference and map labels are the same in the current
      #  stratum class 
      ref_vs_map_all_bool = str_bool & (ref_label == map_label)
      ref_vs_map_all_ind = which(ref_vs_map_all_bool)

      # Calculate proportion (mean) of reference present in stratum
      ref_prop[s, r] = length(ref_vs_str_ind)/sample_totals[,2][sample_totals[,1] == str_codes[s]]
      # Calculate SAMPLE variance of reference in stratum.
      ref_var[s, r] = var(ref_vs_str_bool[str_ind])
      
      # Calculate proportion (mean) of map present in stratum
      map_prop[s, r] = length(map_vs_str_ind)/sample_totals[,2][sample_totals[,1] == str_codes[s]]
      # Calculate SAMPLE variance of map in stratum.
      map_var[s, r] = var(map_vs_str_bool[str_ind])
      
      # Calculate proportion (mean) of map == reference present in stratum for 
      # one particular label
      map_and_ref_prop[s, r] = length(ref_vs_map_ind)/sample_totals[,2][sample_totals[,1] == str_codes[s]]
      # Calculate SAMPLE variance of map == reference.
      map_and_ref_var[s, r] = var(ref_vs_map_bool[str_ind])
      
      # Calculate proportion of map == ref present in stratum for ALL labels
      overall_acc_prop[s] = length(ref_vs_map_all_ind)/sample_totals[,2][sample_totals[,1] == str_codes[s]]
      overall_acc_var[s] = var(ref_vs_map_all_bool[str_ind])
     
      # Calculate covariances associated with users and producers accuracy
      users_cov[s, r] = cov(ref_vs_map_bool[str_ind], map_vs_str_bool[str_ind])
      producers_cov[s, r] = cov(ref_vs_map_bool[str_ind], ref_vs_str_bool[str_ind])
      
    }
  }
  
  # Helper function to assign row and colnames
  assign_names = function(df, rnames, cnames){
    rownames(df) = rnames
    colnames(df) = cnames
    return(df)
  }
  
  # Assign column and row names
  strata_rownames = paste0("strat_",str_codes)
  ref_prop = assign_names(ref_prop, strata_rownames, paste0("ref_", rfcodes))
  ref_var =  assign_names(ref_var, strata_rownames, paste0("ref_", rfcodes))
  map_prop = assign_names(map_prop, strata_rownames, paste0("map_", rfcodes))
  map_var =  assign_names(map_var, strata_rownames, paste0("map_", rfcodes))
  map_and_ref_prop = assign_names(map_and_ref_prop, strata_rownames, 
                                    paste0("mapandref_", rfcodes))
  map_and_ref_var =  assign_names(map_and_ref_var, strata_rownames, 
                                    paste0("mapandref_", rfcodes))
  users_cov = assign_names(users_cov, strata_rownames, paste0("ucov_", rfcodes))
  producers_cov = assign_names(producers_cov, strata_rownames, 
                                paste0("pcov_", rfcodes))
  
  # Calculate total number of pixels in original strata map. Also calculated
  # in the main script, so it's not being returned here.
  totalarea_pix = sum(strata_totals[,2])
  class_prop = vector()
  
  # Calculate ref_label class proportions (i.e. by columns) using total, 
  # original orig_strata areas.
  for (r in 1:ncol(ref_prop)){
    # totalarea_pix is REQUIRED here even if there are no reference counts for 
    # a given stratum
    class_prop[r] = sum(strata_totals[,2] * ref_prop[,r])/totalarea_pix
  }
  
  return(list(ref_prop, ref_var, map_prop, map_var, 
              map_and_ref_prop, map_and_ref_var,
              overall_acc_prop, overall_acc_var,
              users_cov, producers_cov,
              class_prop))
}


#' Function to calculate std error of unbiased area proportions of reference
#' classes for a given year/map
#' 
#' @param strata_totals Dataframe with two columns and number of rows equal to 
#' the total number of classes in the original strata. The first column must 
#' have the same codes found in the original stratification and the second must 
#' have the total number of PIXELS of each class in that original strata map.
#' @param sample_totals Dataframe with two columns and number of rows equal to 
#' the total number of classes in the original strata. The first column must 
#' have the same codes found in the original stratification, and the second must
#' have the total number of SAMPLES of each class collected from that original 
#' strata map.
#' @param ref_var Dataframe with reference class variance for (column) per 
#' original strata class (row).
#' @param rfcodes Vector with all the unique numerical classes present in the 
#' REFERENCE data . This is required to facilitate the calculations for multiple
#' maps/years when not all the reference classes are present in every map.
#' @param totalarea_pix Integer with the total number of pixels present in the 
#' original stratification map
#' @return Vector with standard error of unbiased area proportions per reference
#' class.
#' @export

calc_se_prop = function(strata_totals, sample_totals, ref_var, 
                        rfcodes, totalarea_pix){
  
  # Initialize vector to store results
  se = vector(mode="numeric", length=length(rfcodes))
  
  # Iterate over reference classes
  for (c in 1:length(rfcodes)){
    v = 1/totalarea_pix^2 * (sum(strata_totals[,2]^2 * 
                             (1 - sample_totals[,2]/strata_totals[,2]) * 
                             (ref_var[,c] / sample_totals[,2])))
    se[c] = sqrt(v)
  }
  return(se)
}


#' Function to calculate unbiased area, confidence interval and margin of error
#' 
#' This function takes the area proportions obtained from the function 
#' calc_area_prop and calculates the areas (in ha) as well as the outputs 
#' described below.
#' @param totalarea_pix Integer with the total number of pixels present in the 
#' original stratification map. Assumed to be Landsat (i.e 30 m x 30 m)
#' @param class_prop Vector of area proportions per reference class
#' @param se Vector with standard error of unbiased area proportions per 
#' reference class.
#' @return List with vector of areas in ha (area), vector of HALF the width 
#' confidence interval (ci),
#' vector of higher and lower confidence interval limits (upper_ci, lower_ci) 
#' and margin of error (me)

calc_unbiased_area = function(totarea_pix, class_prop, se){
  # Total area in ha
  N_ha = totarea_pix * 30^2 / 100^2
  # Calculate area in ha from area proportions
  area = class_prop * N_ha
  # Calculate confidence interval in ha
  ci = se * 1.96 * N_ha
  #Upper and lower CI
  upper_ci = area + ci
  lower_ci = area - ci
  me = ci / area 
  return(list(area, ci, upper_ci, lower_ci, me))
}

#' Function to calculate accuracies and their 95% confidence intervals.
#'
#' @param strata_totals Dataframe with two columns and number of rows equal to 
#' the total number of classes in the original strata. The first column must 
#' have the same codes found in the original stratification and the second must 
#' have the total number of PIXELS of each class in that original strata map.
#' @param sample_totals Dataframe with two columns and number of rows equal to 
#' the total number of classes in the original strata. The first column must 
#' have the same codes found in the original stratification, and the second must
#' have the total number of SAMPLE UNITS of each class collected from that 
#' original strata map. 
#' @param rfcodes Vector with numeric values representing the reference codes
#' present in ALL of the periods. 
#' @param totalarea_pix Integer with the total number of pixels present in the 
#' original stratification map. Assumed to be Landsat (i.e 30 m x 30 m)
#' @param ref_prop Dataframe with proportions (mean) of reference labels 
#' (columns) present on each sample strata class (rows)
#' @param ref_var Dataframe of variances of reference labels (columns)
#' present on each sample strata class (rows)
#' @param map_prop Dataframe of proportions (mean) of map labels (columns)
#' present on each sample strata class (rows)
#' @param map_var Dataframe of variances of map labels (columns) present on each
#' sample strata class (rows)
#' @param map_and_ref_prop Dataframe of proportions (mean) of matching reference
#' and map labels (columns) present on each sample strata class (rows)
#' @param map_and_ref_var Dataframe of variances of matching reference and 
#' map labels (columns) present on each sample strata class (rows) 
#' @param overall_acc_prop Vector of proportions of all map labels equal to 
#' reference labels present in stratum 
#' @param overall_acc_var Vector of variances of all map labels equal to 
#' reference labels present in stratum. 
#' @param users_cov Dataframe of user's accuracy covariances per strata
#' @param producers_cov Dataframe of producer's accuracy covariances per strata


#' with its variance (map_and_ref_var), dataframe with proportion of all map 
#' labels equal to reference labels present in stratum (overall_acc_prop) and 
#' dataframe with its variance (overall_acc_var) and dataframes associated with 
#  users's and producer's covariance (users_cov, producers_cov), vector of area 
#' proportions per class (class_prop).
#' @export

calc_accuracies = function(strata_totals, sample_totals, rfcodes, totalarea_pix,
                           ref_prop, ref_var, map_prop, map_var, 
                           map_and_ref_prop, map_and_ref_var,
                           overall_acc_prop, overall_acc_var,
                           users_cov, producers_cov){
  
  # Initialize vector to store results
  se_overall = vector(mode="numeric", length=length(rfcodes))
  se_usr = vector(mode="numeric", length=length(rfcodes))
  se_prod = vector(mode="numeric", length=length(rfcodes))
  
  # Overall accuracy
  overall_acc = sum(strata_totals[,2] * overall_acc_prop) / totalarea_pix
  
  # User's and producers accuracies, common parameter
  param1 = colSums(strata_totals[,2] * map_and_ref_prop)
  
  # Users's accuracy
  uparam = colSums(strata_totals[,2] * map_prop)
  users_acc = param1 / uparam
  names(users_acc) = rfcodes

  # Producer's accuracy
  pparam = colSums(strata_totals[,2] * ref_prop)
  producers_acc = param1 / pparam
  names(producers_acc) = rfcodes

  # Finite population correction term
  corr_term = (1 - sample_totals[,2]/strata_totals[,2])
  
  # Standard error of accuracies
  # Overall accuracy
  vro = 1/totalarea_pix^2 * (sum(strata_totals[,2]^2 * corr_term  * 
                                 overall_acc_var / sample_totals[,2])) 
  se_overall = sqrt(vro)
  
  for (c in 1:length(rfcodes)){
    
    # User's accuracy
    vru = 1/uparam[c]^2 * (sum(strata_totals[,2]^2 * corr_term  * 
                                (map_and_ref_var[,c] + (users_acc[c]^2) * 
                                (map_var[,c]) - 2*users_acc[c]*users_cov[,c]) / 
                                sample_totals[,2]))
    se_usr[c] = sqrt(vru)
    
    # Producer's accuracy, higher than Stehman's paper by 0.002
    vrp = 1/pparam[c]^2 * (sum(strata_totals[,2]^2 * corr_term * 
                               (map_and_ref_var[,c] + (producers_acc[c]^2) * 
                               (ref_var[,c]) - 2*producers_acc[c] * producers_cov[,c]) / 
                               sample_totals[,2]))
    se_prod[c] = sqrt(vrp)
  }
  
  # Calculate confidence intervals
  
  overall_acc_min = (overall_acc - (overall_acc * 1.96 * se_overall))*100
  overall_acc_max = (overall_acc + (overall_acc * 1.96 * se_overall))*100
  users_acc_min =  (users_acc - (users_acc * 1.96 * se_usr))*100
  users_acc_max =  (users_acc + (users_acc * 1.96 * se_usr))*100
  producers_acc_min =  (producers_acc - (producers_acc * 1.96 * se_prod))*100
  producers_acc_max =  (producers_acc + (producers_acc * 1.96 * se_prod))*100
  
  
  return(list(overall_acc*100, overall_acc_min, overall_acc_max,
              users_acc*100, users_acc_min, users_acc_max,
              producers_acc*100, producers_acc_min, producers_acc_max))
  
}


#' Function to obtain a confusion matrix based on two vectors of equal lenght
#' Output matrix is square and has labels resulting from the union of the
#' unique values of the two vectors. First vector is rows, second is columns.
#' @param v1 First vector of integer values
#' @param v2 Second vector of integer values
#' @param code_levels List of integers with all the factors to be used as levels
#' @return ctab Crosstabulation of the two vectors

calc_ct = function(v1, v2, code_levels){
  if(missing(code_levels)){
    l1 = unique(v1)  
    l2 = unique(v2)
    code_levels = sort(union(l1,l2))
    f1 = factor(v1, levels=code_levels)
    f2 = factor(v2, levels=code_levels)
    ct = table(f1, f2)
  } else { 
    f1 = factor(v1, levels=code_levels)
    f2 = factor(v2, levels=code_levels)
    ct = table(f1, f2)
  }
  return(ct)
}


#' Function to calculate optimal sample size for an expected confusion matrix
#' Taken from the Excel spreadsheet provided in Wagner and Stehman 2015. Test
#' data taken from the 4x4 spreadsheet
#' @param cm: Confusion matrix as a dataframe, in proportions. Must be square. 
#' Classes in rows and cols must be in the same order.
#' @param i: Column/row number to be used as the target for the estimation
#' @param n: Total sample size

calc_optimal_sample_alloc = function(cm, i, n){
  if(1/sum(cm) < 0.99){stop("Confusion matrix does not add to 1")}
     
  # Calc variances and sort them according to position of
  # class of interest (i.e. 'i')
  csum = sum(cm[,i])
  var_user = (cm[i,i] * (sum(cm[i,-i]))) / (sum(cm[i,]))^2
  
  var_prod_class = (sum(cm[-i,i])^2 * cm[i,i] * (sum(cm[i, -i]))) / 
        csum^4
  var_prod_others = (cm[i,i]^2 * cm[-i,i] * rowSums(cm[-i,-i])) / 
        csum^4
  
  var_prod = vector(length = length(cm))
  var_prod[i] = var_prod_class
  var_prod[-i] = var_prod_others
  
  var_parea =  cm[,i] * (rowSums(cm[, -i]))
  
  # Calculate K's and min K =! 0
  k = vector(length = length(cm))
  k[-i] = var_prod[-i] + var_parea[-i]
  k[i] = var_user + var_prod[i] + var_parea[i]
  mink = min(k[k>0])
  
  # Calc theoretical smallest sample size
  nhat = (k / mink)^0.5
  
  # Calc optimal sample size
  nh = (nhat / sum(nhat)) * n
  
  # Calc standard errors
  if (nh[i] == 0){
    se_user = 0
  } else {
    se_user = sqrt(var_user / nh[i])
  }
  
  se_prod = sum(var_prod[nh > 0] / nh[nh > 0])^0.5
  se_parea = sum(var_parea[nh > 0] / nh[nh > 0])^0.5

  out = list(sample_alloc=nh, 
             se_user=se_user,
             se_prod=se_prod,
             se_parea=se_parea)

  return(out)
}

# Test data for the previous function. 
testcm = data.frame(c(0.01400, 0.00100, 0.00200, 0.00400), 
                    c(0.00000,0.00900, 0.00000, 0.00200),
                    c(0.00300, 0.00250, 0.28800, 0.02500),
                    c(0.00300, 0.00250, 0.03000, 0.61400))
colnames(testcm) = c("1", "2", "3", "4")


