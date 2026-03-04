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

## stud_area 
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/stud_site.gpkg"))
stud_area<-st_transform(stud_area,2154)%>%filter(siteID =="FRA_BAR2")
st_area(stud_area)/ 1e6
### make pu
target_crs <- st_crs(stud_area)$wkt  # Use WKT for terra

## transform and resample rasters
cost <- project(cost, target_crs)
grand_mean_es <- project(grand_mean_es, target_crs)
mean_cult <- project(mean_cult, target_crs)
mean_reg <- project(mean_reg, target_crs)
mean_prov <- project(mean_prov, target_crs)

habitat <-project(habitat, target_crs)
habitat_main<-floor(habitat / 100)
pop25 <-project(pop25, target_crs)
pop30 <-project(pop30, target_crs)

## existing PA
PA<-st_read(paste0(main_dir,"WP4/pa_existing/WDPA_FRL04.shp"))
PA<-st_transform(PA,st_crs(target_crs))
PA <- st_intersection(PA, stud_area)
PA<-st_make_valid(PA)


PA <- PA %>%
  mutate(class = case_when(
    IUCN_CAT == "Not Applicable" ~ 1,
    IUCN_CAT == "Not Assigned" ~ 1,
    IUCN_CAT == "Not Reported" ~ 1,
    IUCN_CAT == "IV" ~ 2,
    IUCN_CAT == "V" ~ 3,
    IUCN_CAT == "IV" ~ 4,
    IUCN_CAT == "III" ~ 5,
    IUCN_CAT == "II" ~ 6,
    IUCN_CAT == "Ia" ~ 7,
    TRUE ~ 0
  ))


# 3. Create a regular 1km x 1km grid over the polygon
#grid <- st_make_grid(stud_area, cellsize = 0.02, square = TRUE)

grid <- st_make_grid(stud_area, cellsize = 1000, square = TRUE)
grid <- st_sf(geometry = grid)
grid$area<-st_area(grid)

# 5. Sample mean value from each raster layer into the grid cells
grid$sampled_es <- terra::extract(grand_mean_es, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_cult <- terra::extract(mean_cult, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_prov <- terra::extract(mean_prov, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_reg <- terra::extract(mean_reg, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_pop25 <- terra::extract(pop25, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_pop30 <- terra::extract(pop30, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_cost <- terra::extract(cost, grid, fun = mean, na.rm = TRUE)[,2]
grid$sampled_condition<- terra::extract(es_cond, grid, fun = mean, na.rm = TRUE)[,2]

grid$sampled_habitat <- terra::extract(
  habitat_main,
  vect(st_centroid(grid))
)[,2]


#sampling number of PA in one cell
# get intersection index list
ints <- st_intersects(grid, PA)

# count unique IUCN categories per grid cell
grid$n_pa <- sapply(ints, function(i) {
  length(unique(PA$IUCN_CAT[i]))
})

#sample highest status of protection
grid <- st_join(grid, PA["class"], join = st_intersects)

grid <- grid %>%
  group_by(across(-class)) %>%
  summarise(
    max_IUCN_class = if (all(is.na(class))) NA_real_ else max(class, na.rm = TRUE)
  ) %>%
  ungroup()


grid$area<-st_area(grid)
#calculate distance from cells outside core PA to core PA (IUCN Ia II)
out   <- grid%>%filter(is.na(max_IUCN_class)|max_IUCN_class<6)
pa_core <- grid%>%filter(max_IUCN_class>5)
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


st_write(grid_clean, paste0("WP4/2_output/02_optim/",siteID,"_input_grid.json"), driver = "GeoJSON", overwrite = T)
