# Script to test the script that calculates areas, accuracies and confidence
# interval to make sure they match the results in Stehman 2014.

setwd("/home/paulo/workflow/multi_scene/")
source("7_poststratification/functions.R")
test_data = read.csv("other/stehman_test_data.csv")
test_data$stratum = as.factor(test_data$stratum)

# Stratum sizes
stratum_sizes = cbind(c(1,2,3,4), c(40000, 30000, 20000, 10000))
sample_sizes = cbind(c(1,2,3,4), c(10, 10, 10, 10))
codes = c("A", "B", "C", "D")

test_output = calc_area_prop(test_data$stratum, test_data$ref_class, test_data$map_class, stratum_sizes, sample_sizes, codes)
test_output

