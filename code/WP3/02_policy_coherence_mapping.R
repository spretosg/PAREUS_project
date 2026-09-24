############################################################################
# Pareus WP3 - mapping policy coherence
# Date: 06.08.2026
# Author: Reto Spielhofer (Norwegian Institute for Nature Research)
############################################################################
# takes a land cover map and the policy coherence data

library(terra)
library(sf)
library(dplyr)
## map policy according to influence on landscape type
proj_id<-"TRD"
main_dir<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",proj_id,"/raw_data_backup/")

## CORINE land use land cover with 24 classes
lulc<-rast(paste0("data/shared/",proj_id,"_lulc.tif"))
lulc[lulc == 0] <- NA

#import lulc - policy rating
lulc_resistance<-read.csv("P:/312204_pareus/pareus_repository/WP3/policy_resistance_mapping.csv",sep=";")
lulc_resistance<-lulc_resistance%>%filter(studyID == proj_id)
lulc_resistance$policy_resistance <- as.numeric(gsub(",", ".", lulc_resistance$policy_resistance))

#recl lulc to 5 classes
lulc <- lulc %>%
  (\(x) floor(x / 100))() %>%
  disagg(fact = 5, method = "near")
rcl <- as.matrix(lulc_resistance[, c("CLC", "policy_resistance")])

# Map LULC classes to lookup values
pol_res <- classify(
  lulc,
  rcl,
  others = NA
)
plot(pol_res)
writeRaster(pol_res,paste0("data/WP3/pol_resistance_",proj_id,".tif"))
