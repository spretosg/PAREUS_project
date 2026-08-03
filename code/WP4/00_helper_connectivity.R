library(circuitscaper)
#cs_install_julia()
Sys.which("julia")
system("julia --version")
#cs_setup(threads = 4)

# library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
# library(gdistance)

source("WP4/1_code/wp4_functions_utils.R")
main_dir<-"P:/312204_pareus/"
siteID<-"FRL04"
####---- user parameter ----####
##factors to multiply the ecosystem condition resistance for movement based on LULC classes 
## higher values == easier to move through area general assumptions
factors <- data.frame(
  class = c(1, 2, 3, 4, 5),
  factor = c(0.1, 0.5, 1, 0.7,0.7)
)

# raster aggregation factor to speed up connectivity calculation
agg_factor = 1

####---- Input and processing ----####
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRL04")
# ecosystem condition and landcover
es_cond<-rast(paste0(main_dir,"WP4/features/",siteID,"_ec.tif"))
lulc<-rast(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/2_env_var/lulc.tif"))

lulc[lulc == 0] <- NA
lulc <- trunc(lulc / 100)

#### ---- resistance grid and connectivity ----
lulc <- lulc %>%
  project(y = crs(es_cond)) %>%
  resample(y = es_cond, method = "near")


factor_raster <- classify(
  lulc,
  rcl = as.matrix(factors)
)

# --- Multiply ec raster by factors ---
r_adjusted <- es_cond * factor_raster

resistance<-log(1/r_adjusted)

r_coarse <- aggregate(resistance, fact = agg_factor)
r_coarse[r_coarse <= 0] <- 1e-6

# calculate moving window connectivity based on ecosystem condition in landscape
start<-Sys.time()
mw_result <- os_run(r_coarse, radius = 20, block_size = 10)
plot(mw_result$normalized_current)
print(Sys.time()-start)
writeRaster(mw_result,paste0("WP4/2_output/02_optim/mw_connectivity_",siteID,".tif"))
