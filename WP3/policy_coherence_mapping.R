# policy mapping
# takes a land cover map and the policy coherence data

library(terra)
library(sf)
library(dplyr)
## map policy according to influence on landscape type
proj_id<-"FRL04"
main_dir<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",proj_id,"/raw_data_backup/")

## CORINE land use land cover with 24 classes
lulc<-rast(paste0(main_dir,"/2_env_var/lulc.tif"))
lulc[lulc == 0] <- NA

#import lulc - policy rating
lulc_rating<-read.csv("P:/312204_pareus/WP3/lulc_policy_FRL04.csv", header = T, sep=";",dec = ",",row.names = "X")
colnames(lulc_rating)<-lulc_rating$X

stack_imp_pol_lulc <- rast()
lulc_rating <- as.matrix(lulc_rating)
for (i in 1:nrow(lulc_rating)) {
  tmp_policy_name<-rownames(lulc_rating)[i]

  rcl <- matrix(c(
    1,  139,  lulc_rating[i,1],   # artificial
    140, 199, lulc_rating[i,1],  # urban green
    200, 219,  lulc_rating[i,2],   # aggriculture 
    220, 229,  lulc_rating[i,2],   # aggriculture
    230, 239,  lulc_rating[i,2],   # aggriculture
    240, 299, lulc_rating[i,2],    # aggriculture
    300, 319, lulc_rating[i,3], ## Forest
    320, 329, lulc_rating[i,3], ## scrub
    330, 399, lulc_rating[i,3], # open spaces
    400, 499,  lulc_rating[i,4],    # wetlands
    500, 599, lulc_rating[i,5] #water
  ), ncol = 3, byrow = TRUE)
  
  # Apply reclassification
  lulc_policy <- classify(lulc, rcl)
  lulc_policy_scaled<-lulc_policy/3
  plot(lulc_policy_scaled)

  names(lulc_policy_scaled)<-tmp_policy_name
  stack_imp_pol_lulc<-c(stack_imp_pol_lulc,lulc_policy_scaled)
}


# Loop over layers and export
for (i in seq_along(names(stack_imp_pol_lulc))) {
  writeRaster(
    stack_imp_pol_lulc[[i]],
    filename = file.path("P:/312204_pareus/WP3/pol_impact_lulc",proj_id,"/", paste0(names(stack_imp_pol_lulc)[i], ".tif")),
    overwrite = TRUE
  )
}



