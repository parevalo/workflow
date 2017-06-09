auxpath = "/media/paulo/785044BD504483BA/test/"

setwd(auxpath)
mapped_areas_list = list()
years=c("01_03", "03_05", "05_07", "07_09", "09_11", "11_13", "13_15")
strata_03_05_pixcount_original_buff3
cr = c(7, 10, 12, 15)

for(i in 1:7){
  fname = paste0("strata_", years[i], "_pixcount_original_buff3.csv")
  mapped_areas_list[[i]] = read.csv(paste0(auxpath,fname), header=TRUE, col.names=c("stratum", "pixels"))
  # Get rid of rows we don't need, eg class 15, also done below for ss
  mapped_areas_list[[i]] = mapped_areas_list[[i]][!(mapped_areas_list[[i]]$stratum %in% cr),]
}

# Create empty matrix to store mapped values and fill
mapped_areas = matrix(0, nrow=length(years), ncol=nrow(mapped_areas_list[[1]]), byrow=T, dimnames = list(years, mapped_areas_list[[1]][,1]))

for (i in 1:length(mapped_areas_list)){
  mapped_areas[i,] = mapped_areas_list[[i]][,2]  
}

# Calculate strata weights for each of them
mapped_weights = mapped_areas / rowSums(mapped_areas)

# Convert areas to ha
mapped_areas = mapped_areas * 30^2 / 100^2

# Save mapped areas and strata weights to file.
write.csv(mapped_weights, "mapped_weights_biannual_stratification.csv")
