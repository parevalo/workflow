# Function to calculate the strata for any given pair of years. If an invalid
# class code is provided, NA is returned

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

