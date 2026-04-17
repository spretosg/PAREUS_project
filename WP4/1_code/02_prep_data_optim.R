## test data read tif
## make pu of stud area
## sample cost, es_cond and es_capacity per PU

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("WP4/1_code/wp4_functions_utils.R")


main_dir<-"P:/312204_pareus/"
siteID<-"FRA_BAR2"

## read grid from 01_pa_status.R
grid<-st_read(paste0("WP4/2_output/02_optim/",siteID,"_input_grid.json"))
target_crs<-st_crs(grid)$wkt

## read cost based on es
cost<-terra::rast(paste0(main_dir,"WP4/cost_raster_es/",siteID,"_cost_raster_es.tif"))

#read all es_make mean es as features 1
es_raster_files <- list.files(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/4_mean_R1"),pattern = "\\.tif$", full.names = TRUE)
# 2. Read all rasters
es_raster <- terra::rast(es_raster_files)

grand_mean_es<-mean(es_raster)
#prov es mean
mean_prov <- mean(es_raster[[c(3, 6, 10, 11)]])
#cult es mean
mean_cult <- mean(es_raster[[c(1, 7, 9)]])
#reg es mean
mean_reg <- mean(es_raster[[c(2, 4, 5)]])

## es condition as features 2
es_cond<-paste0(main_dir,"WP4/features/",siteID,"_ec.tif")
es_cond<-terra::rast(es_cond)

##habitat (from lulc)
habitat<-paste0(main_dir,"WP4/habitat/",siteID,"_lulc.tif")
habitat<-terra::rast(habitat)

##pop (world pop cover 2020)
pop25<-paste0(main_dir,"WP4/lock_out/pop_2025_",siteID,".tif")
pop25<-terra::rast(pop25)

pop30<-paste0(main_dir,"WP4/lock_out/pop_2030_",siteID,".tif")
pop30<-terra::rast(pop30)



## transform and resample rasters
cost <- project(cost, target_crs)
grand_mean_es <- project(grand_mean_es, target_crs)
mean_cult <- project(mean_cult, target_crs)
mean_reg <- project(mean_reg, target_crs)
mean_prov <- project(mean_prov, target_crs)

# habitat <-project(habitat, target_crs)

pop25 <-project(pop25, target_crs)
pop30 <-project(pop30, target_crs)

## existing PA



# 5. Sample mean value from each raster layer into the grid cells
grid$sampled_es <- terra::extract(grand_mean_es, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_cult <- terra::extract(mean_cult, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_prov <- terra::extract(mean_prov, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_reg <- terra::extract(mean_reg, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_pop25 <- terra::extract(pop25, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_pop30 <- terra::extract(pop30, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_cost <- terra::extract(cost, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_condition<- terra::extract(es_cond, grid, fun = mean, na.rm = TRUE)[,2]




#calculate distance from cells outside core PA to core PA (IUCN Ia II)
out   <- grid%>%filter(is.na(class)|class<6)
pa_core <- grid%>%filter(class>5)
dist_matrix <- st_distance(out, pa_core)

out$min_distance <- apply(dist_matrix, 1, min)
pa_core$min_distance<-0
grid<-rbind(out,pa_core)


grid_clean <- grid %>%
  filter(!if_any(all_of(c("sampled_es","sampled_cult","sampled_prov","sampled_reg","sampled_cost","sampled_condition","sampled_habitat")), ~ is.na(.) | is.nan(.)))



grid_clean<- zero_one_scale(
  grid_clean,
  cols = c("sampled_pop25", "sampled_pop30", "sampled_es","sampled_cult","sampled_prov","sampled_reg","sampled_cost","sampled_condition","min_distance")
)


st_write(grid_clean, paste0("WP4/2_output/02_optim/",siteID,"_input_final_grid.json"), driver = "GeoJSON", overwrite = T)
