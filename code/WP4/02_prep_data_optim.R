## test data read tif
## make pu of stud area
## sample cost, es_cond and es_capacity per PU

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("WP4/1_code/wp4_functions_utils.R")
main_dir<-"P:/312204_pareus/"
siteID<-"FRL04"


####---- Input and processing ----####
## read PU from 01_pa_status.R
grid<-st_read(paste0("WP4/2_output/01_PA_analysis/",siteID,"_input_grid.json"))

### Features
## Mean ecosystem service (cookbook nr.2)
grand_mean_es<-terra::rast(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/es_weight_mean.tif"))
# ES groups
es_raster <- terra::rast(list.files(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/4_mean_R1"),pattern = "\\.tif$", full.names = TRUE))
#prov es mean
mean_prov <- mean(es_raster[[c(3, 6, 10, 11)]])
#cult es mean
mean_cult <- mean(es_raster[[c(1, 7, 9)]])
#reg es mean
mean_reg <- mean(es_raster[[c(2, 4, 5)]])

## Ecosystem condition (cookbook nr.2)
es_cond<-terra::rast(paste0(main_dir,"WP4/features/",siteID,"_ec.tif"))

### costs
##based on ecosystem service exclusiveness
cost_es<-terra::rast(paste0(main_dir,"WP4/cost_raster_es/",siteID,"_cost_raster_es.tif"))
## based on policy coherence (cookbook nr. 3)
cost_policy<-terra::rast(paste0(main_dir,"WP4/cost_raster_policy/",siteID,"_cost_raster_pol.tif"))

### Further features of potential interest
##pop (world pop cover 2020)
pop25<-terra::rast(paste0(main_dir,"WP4/lock_out/pop_2025_",siteID,".tif"))
pop30<-terra::rast(paste0(main_dir,"WP4/lock_out/pop_2030_",siteID,".tif"))

target_crs<-st_crs(grid)$wkt

rasters <- list(
  cost_es       = cost_es,
  cost_policy   = cost_policy,
  grand_mean_es = grand_mean_es,
  es_cult       = mean_cult,
  es_reg        = mean_reg,
  es_prov       = mean_prov,
  pop25         = pop25,
  pop30         = pop30,
  eco_cond = es_cond
)

rasters <- lapply(rasters, project, y = target_crs)

list2env(rasters, envir = environment())

grid[paste0("sampled_", names(rasters))] <-
  lapply(rasters, \(r)
         terra::extract(r, grid, fun = mean, na.rm = TRUE)[, 2]
  )

#calculate distance from cells outside core PA to core PA as defined in step 01_pa_status.R
out   <- grid%>%filter(exisiting_corePA == F)
pa_core <- grid%>%filter(exisiting_corePA == T)

out$min_distance <- apply(st_distance(out, pa_core), 1, min)
pa_core$min_distance<-0
grid<-rbind(out,pa_core)

## remove na grid cells based on features and costs (not)
grid_clean <- grid %>%
  filter(!if_any(all_of(c("sampled_grand_mean_es","sampled_es_cult","sampled_es_prov","sampled_es_reg","sampled_cost_es","sampled_cost_policy","sampled_eco_cond")), ~ is.na(.) | is.nan(.)))

grid_clean<- zero_one_scale(
  grid_clean,
  cols = c("sampled_pop25", "sampled_pop30", "sampled_grand_mean_es","sampled_es_cult","sampled_es_prov","sampled_es_reg","sampled_cost_es","sampled_eco_cond","min_distance")
)


st_write(grid_clean, paste0("WP4/2_output/02_optim/",siteID,"_input_final_grid.json"), driver = "GeoJSON", overwrite = T)
