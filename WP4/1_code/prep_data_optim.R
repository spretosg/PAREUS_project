## test data read tif
## make pu of stud area
## sample cost, es_cond and es_capacity per PU
## save as json upload gcs bucket for python
library(sf)
library(terra)
library(dplyr)


main_dir<-"P:/312204_pareus/"
siteID<-"SK021"

## read cost
cost<-terra::rast(paste0(main_dir,"WP4/cost_raster_es/",siteID,"_cost_raster_es.tif"))

#read all es_make mean es as features 1
es_raster_files <- list.files(paste0(main_dir,"WP2/PGIS_ES_mapping/raw_data_backup/",siteID,"/6_mean_R2"), full.names = TRUE)
# 2. Read all rasters
es_raster <- terra::rast(es_raster_files)

mean_es<-mean(es_raster)


## es condition as features 2
es_cond<-paste0(main_dir,"WP2/es_condition/",siteID,"_es_cond.tif")
es_cond<-terra::rast(es_cond)

##habitat
habitat<-paste0(main_dir,"WP4/habitat/",siteID,"_habitat.tif")
habitat<-terra::rast(habitat)


## transform and resample
cost <- project(cost, es_cond)
mean_es <- project(mean_es, es_cond)
habitat <-project(habitat, es_cond)

cost <- resample(cost, es_cond, method = "bilinear") 
mean_es <- resample(mean_es, es_cond, method = "bilinear") 
habitat <- resample(habitat, es_cond, method = "min") 
### combine feat 1 / 2 into i
feat <- c(es_cond, mean_es, habitat)
names(feat) <- c("es_conditon", "mean_es", "habitat")

terra::writeRaster(feat,paste0(main_dir,"WP4/features/",siteID,"_optim_features.tif"), overwrite=T)
terra::writeRaster(cost,paste0(main_dir,"WP4/cost_raster_es/",siteID,"_optim_cost.tif"))

## lock in (existing PA)
# PA1<-st_read(paste0(main_dir,"WP4/pa_existing/PA_",siteID,"/PA1.shp"))
# PA2<-st_read(paste0(main_dir,"WP4/pa_existing/PA_",siteID,"/PA2.shp"))
# PA3<-st_read(paste0(main_dir,"WP4/pa_existing/PA_",siteID,"/PA3.shp"))
# 
# PA<-rbind(PA1,PA2,PA3)
# PA<-st_transform(PA,crs(es_cond))
# # Crop and mask with raster
# PA <- vect(PA)
# PA <- crop(PA, es_cond)
# pa_raster <- rasterize(PA, es_cond, field = 1, background = 0)
# terra::writeRaster(pa_raster,paste0(main_dir,"WP4/pa_existing/",siteID,"_pa.tif"))

## stud_area svk
svk<-read_sf(paste0(main_dir,"WP2/PGIS_ES_mapping/raw_data_backup/",siteID,"/study_site.gpkg"))%>%filter(cntrID == "SVK")

### make pu
target_crs <- st_crs(svk)$wkt  # Use WKT for terra

# 3. Reproject rasters to the target CRS
cost <- project(cost, target_crs)
mean_es <- project(mean_es, target_crs)
es_cond <- project(es_cond, target_crs)
habitat <- project(habitat, target_crs)

# Resample all to match the first projected raster
mean_es <- terra::resample(mean_es, cost)
es_cond <- terra::resample(es_cond,cost)
habitat3 <-terra::resample(habitat,cost, method = "min")

min_max_normalize <- function(r) {
  # Normalize each layer individually
  norm_r <- lapp(r, function(x) {
    r_min <- min(x, na.rm = TRUE)
    r_max <- max(x, na.rm = TRUE)
    (x - r_min) / (r_max - r_min)
  })
  return(norm_r)
}
# 2. Transform CRS to match raster (optional but recommended)
# Assuming your raster is already loaded as `r_stack`
mean_es<-min_max_normalize(mean_es)
cost<-min_max_normalize(cost)
es_cond<-min_max_normalize(es_cond)

# 3. Create a regular 1km x 1km grid over the polygon
grid <- st_make_grid(svk, cellsize = 0.02, square = TRUE)
grid <- st_sf(geometry = grid)

# Clip the grid to polygon boundary
grid_clipped <- st_intersection(grid, svk)%>%dplyr::select()

# 4. Convert grid to SpatVector for terra compatibility
grid_vect <- vect(grid_clipped)

# 5. Sample mean value from each raster layer into the grid cells
sampled_es <- terra::extract(mean_es, grid_vect, fun = mean, na.rm = TRUE)
sampled_cost <- terra::extract(cost, grid_vect, fun = mean, na.rm = TRUE)
sampled_condition <- terra::extract(es_cond, grid_vect, fun = mean, na.rm = TRUE)
sampled_habitat <- terra::extract(habitat3, grid_vect, fun = max, na.rm = TRUE)
sampled_habitat$habitat<-as.integer(sampled_habitat$habitat)

# 6. Combine sampled data back with grid for NSGA py optimization
grid_data <- cbind(grid_clipped, sampled_es[, -1], sampled_cost[, -1],sampled_condition[, -1],sampled_habitat[, -1])  # remove ID column
colnames(grid_data)<-c("es_service","cost","es_condition","habitat","geometry")
st_write(grid_data, "SVK_test4.json", driver = "GeoJSON")
