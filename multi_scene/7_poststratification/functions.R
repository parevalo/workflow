#' THIS SCRIPT CONTAINS ALL THE FUNCTIONS REQUIRED FOR THE CALCULATION OF
#' THE UNBIASED AREAS IN THE CALCULATE_STRATA_PER_YEAR.R SCRIPT:

#' 1) Function to calculate strata for any given pair of years
#' 2) Calculate unbiased area proportions and variance of reference classes
#' 3) Calculate standard error of unbiased area proportion
#' 4) Calculate unbiased areas, their confidence interval and margin of error

#' Function to calculate the strata for any given pair of years (integers). 
#' If an invalid class code is provided, NA is returned

calculate_strata <- function(year1, year2){
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


#' Calculate area proportions of reference classes and variance of reference samples per original strata
#' for a given year/map.
#'
#' @param samp_strata Vector with numeric codes representing the original stratification of each sample
#' @param samp_reference Vector with numeric codes representing the reference stratification for that year/map
#' @param strata_totals Dataframe with two columns and number of rows equal to the total number of classes 
#' in the original strata. The first column must have the same codes found in the original stratification 
#' and the second must have the total number of PIXELS of each class in that original strata map.
#' @param sample_totals Dataframe with two columns and number of rows equal to the total number of classes 
#' in the original strata. The first column must have the same codes found in the original stratification, 
#' and the second must have the total number of SAMPLES of each class collected from that original strata map. 
#' @param rfcodes Vector with all the unique numerical classes present in the REFERENCE data . This is required
#' to facilitate the calculations for multiple maps/years when not all the reference classes are present in every map.
#' @return List with: Vector of area proportions per class (class_prop), dataframe with sample variance per sample strata (ref_var), 
#' dataframe with original strata codes and total area in pixels present in the map/year being evaluated (fss), 
#' vector with proportion of sample reference class present on each sample strata class (ref_prop)
#' @export

calc_area_prop = function(samp_strata, samp_reference, strata_totals, sample_totals, rfcodes){
  #TODO: Test if class codes match when they should?
  
  # Obtain unique values in the samp_strata field
  str_codes = sort(unique(samp_strata))
  str_ind = seq_along(str_codes)
  
  # Get a sequence for the reference codes
  ref_ind = seq_along(rfcodes)
  
  # Initialize empty df for proportions per samp_strata, and for sample variance per samp_strata
  ref_prop = data.frame()
  ref_var = data.frame()
  
  # Compare the fields, iterate over "samp_strata" and "samp_reference" classes
  for (s in str_ind){
    for (r in ref_ind){
      # Compared fields and get TRUE or FALSE on each row
      cond_bool = samp_strata == str_codes[s] & samp_reference == rfcodes[r]
      # Get row numbers that meet that condition
      ind = which(cond_bool)
      # Get location of rows for the current samp_strata
      str_ref = which(samp_strata == str_codes[s])
      # Get proportion of samp_reference class present on that samp_strata
      ref_prop[s, r] = length(ind)/sample_totals[,2][sample_totals[,1] == str_codes[s]]
      # Calculate SAMPLE variance of current samp_reference code in current samp_strata (needed later)
      # which is the same formula specified in the paper
      ref_var[s, r] = var(cond_bool[str_ref]) 
      
    }
  }
  
  # Assign column and row names
  rownames(ref_prop) = paste0("strat_",str_codes)
  colnames(ref_prop) = paste0("ref_", rfcodes)
  rownames(ref_var) = paste0("strat_",str_codes)
  colnames(ref_var) = paste0("ref_", rfcodes)
  
  # Calculate total number of pixels in original strata map
  totalarea_pix = sum(strata_totals[,2])
  class_prop = vector()
  # Filter only total sample sizes that are present in the samp_strata for that year
  fss = strata_totals[strata_totals[,1] %in% str_codes,]
  # Calculate samp_reference class proportions (i.e. by columns) using total, original samp_strata areas.
  for (r in 1:ncol(ref_prop)){
    # LEAVE THE SUM OF THE ENTIRE samp_strata. #CHECK!!!!!!
    class_prop[r] = sum(fss[,2] * ref_prop[,r])/totalarea_pix
  }
  return(list(class_prop, ref_var, fss, ref_prop, totalarea_pix)) 
}


#' Function to calculate standard error of unbiased area proportions of reference classes for a given year/map
#' 
#' @param strata_totals Dataframe with two columns and number of rows equal to the total number of classes 
#' in the original strata. The first column must have the same codes found in the original stratification 
#' and the second must have the total number of PIXELS of each class in that original strata map.
#' @param sample_totals Dataframe with two columns and number of rows equal to the total number of classes 
#' in the original strata. The first column must have the same codes found in the original stratification, 
#' and the second must have the total number of SAMPLES of each class collected from that original strata map.
#' @param ref_var Dataframe with reference class variance for (column) per original strata class (row).
#' @param rfcodes Vector with all the unique numerical classes present in the REFERENCE data . This is required
#' to facilitate the calculations for multiple maps/years when not all the reference classes are present in every map.
#' @param totalarea_pix Integer with the total number of pixels present in the original stratification map
#' @return Vector with standard error of unbiased area proportions per reference class.
#' @export

calc_se_prop = function(strata_totals, sample_totals, ref_var, rfcodes, totalarea_pix){
  
  # Initialize vector to store results
  se = vector(mode="numeric", length=length(rfcodes))
  
  # Iterate over reference classes
  for (c in 1:length(rfcodes)){
    v = 1/totalarea_pix^2 * (sum(strata_totals[,2]^2 * (1 - sample_totals[,2]/strata_totals[,2]) * (ref_var[,c] / sample_totals[,2])))
    se[c] = sqrt(v)
  }
  return(se)
}


#' Function to calculate unbiased area, confidence interval and margin of error
#' 
#' This function takes the area proportions obtained from the function calc_area_prop and
#' calculates the areas (in ha) as well as the outputs described below.
#' @param totalarea_pix Integer with the total number of pixels present in the original stratification map.
#' Assumed to be Landsat (i.e 30 m x 30 m)
#' @param class_prop Vector of area proportions per reference class
#' @param se Vector with standard error of unbiased area proportions per reference class.
#' @return List with vector of areas in ha (area), vector of HALF the width confidence interval (ci),
#' vector of higher and lower confidence interval limits (upper_ci, lower_ci) and margin of error (me)

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

