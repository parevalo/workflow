# Script to test the script that calculates areas, accuracies and confidence
# interval to make sure they match the results in Stehman 2014.

setwd("/home/paulo/workflow/multi_scene/")
source("7_poststratification/functions.R")
test_data = read.csv("other/stehman_test_data.csv")
test_data$stratum = as.factor(test_data$stratum)

# Stratum sizes
stratum_sizes = cbind(c(1,2,3,4), c(40000, 30000, 20000, 10000))
totalsize = sum(stratum_sizes)
sample_sizes = cbind(c(1,2,3,4), c(10, 10, 10, 10))
codes = c("A", "B", "C", "D")

# Area proportions
test_aprop = calc_props_and_vars(test_data$stratum, test_data$ref_class, test_data$map_class, stratum_sizes, sample_sizes, codes)
#test_aprop

# Standard errors
test_se_prop = calc_se_prop(stratum_sizes, sample_sizes, test_aprop[[2]], codes, totalsize)
#test_se_prop

# Areas
test_areas = calc_unbiased_area(totalsize, test_aprop[[9]], test_se_prop)
#test_areas

# Accuracies and SE, update with new outputs for overall acc
test_accuracies = calc_accuracies(stratum_sizes, sample_sizes, codes, totalsize, 
                                  test_aprop[[1]], test_aprop[[2]], test_aprop[[3]], test_aprop[[4]],
                                  test_aprop[[5]], test_aprop[[6]], test_aprop[[7]],test_aprop[[8]],
                                  test_aprop[[9]],test_aprop[[10]])
test_accuracies

