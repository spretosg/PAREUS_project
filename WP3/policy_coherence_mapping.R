# policy mapping
# takes a land cover map and the policy coherence data

library(terra)
library(sf)
library(dplyr)
## map policy according to influence on landscape type
proj_id<-"SK021"
main_dir<-paste0("P:/312204_pareus/WP2/PGIS_ES_mapping/raw_data_backup/",proj_id)

## CORINE land use land cover with 24 classes
lulc<-rast(paste0(main_dir,"/2_env_var/lulc.tif"))

#import lulc - policy rating
lulc_rating<-read.csv("P:/312204_pareus/WP3/lulc_policy.csv")
stack_imp_pol_lulc<-c()

for (i in 1:nrow(lulc_rating)) {
  tmp_policy_name<-lulc_rating[i,1]
  
  # reclassified policy lulc with 5 classes
  rcl <- matrix(c(
    1,  139,  lulc_rating[i,2],   # urban sealed
    140, 199, lulc_rating[i,3],  # urban green
    200, 219,  lulc_rating[i,4],   # aggriculture arable land
    220, 229,  lulc_rating[i,5],   # permanent crops
    230, 239,  lulc_rating[i,6],   # pastures
    240, 299, lulc_rating[i,7],    # heterogenous agri
    300, 319, lulc_rating[i,8], ## Forest
    320, 329, lulc_rating[i,9], ## scrub
    330, 399, lulc_rating[i,10], # open spaces
    400, 499,  lulc_rating[i,11],    # wetlands
    500, 599, lulc_rating[i,12] #water
  ), ncol = 3, byrow = TRUE)
  
  # Apply reclassification
  lulc_policy <- classify(lulc, rcl)
  names(lulc_policy)<-tmp_policy_name
  stack_imp_pol_lulc<-c(stack_imp_pol_lulc,lulc_policy)
}

layer_names <- c("eu_bio_div","ELC","terr_system_ecol_stab","env_strat_2023","svk_agri_plan","act_nat_land_prot")

# Loop over layers and export
for (i in seq_along(layer_names)) {
  writeRaster(
    stack_imp_pol_lulc[[i]],
    filename = file.path("P:/312204_pareus/WP3", paste0(layer_names[i], ".tif")),
    overwrite = TRUE
  )
}



