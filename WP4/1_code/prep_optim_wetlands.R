## test data read tif
## make pu of stud area
## sample cost, es_cond and es_capacity per PU
## save as json upload gcs bucket for python

library(sf)
library(terra)
library(dplyr)
source("wp4_functions_utils.R")

########################################
############ GENERAL ###################
########################################

main_dir<-"P:/312204_pareus/"
sid<-"SK021"
stud_area<-read_sf(paste0(main_dir,"WP2/PGIS_ES_mapping/raw_data_backup/",sid,"/study_site.gpkg"))%>%filter(siteID %in% sid)
target_crs <- st_crs(stud_area)$wkt

PA<-st_read(paste0(main_dir,"WP4/pa_existing/PA_",sid,"/PA_WGS84.shp"))
PA<-st_transform(PA,target_crs)
PA<-st_union(PA)
PA <- st_sf(geometry = PA)

# create a grid for defining the planning units PU for optimization
grid <- st_make_grid(stud_area, cellsize = 0.01, square = TRUE)
grid <- st_sf(geometry = grid)
# Clip the grid to polygon boundary
grid_clipped <- st_intersection(grid, stud_area)%>%dplyr::select()
grid_clipped<-st_make_valid(grid_clipped)

#### sample if cell is in PA (TRUE) or not (FALSE)
pa_intersects_cell<-st_intersects(grid_clipped,PA, sparse = FALSE)
pa_intersects_cell<-pa_intersects_cell[,1]

#Convert grid to SpatVector for terra compatibility
grid_vect <- vect(grid_clipped)

######################
#### OPTIM VARS #######
#######################

#### es condition from INRAE (maximize)
es_cond<-paste0(main_dir,"WP2/es_condition/",sid,"_es_cond.tif")
es_cond<-terra::rast(es_cond)
es_cond <- project(es_cond, target_crs)

#Mean ES from PGIS (maximize)
es_raster_files <- list.files(paste0(main_dir,"WP2/PGIS_ES_mapping/raw_data_backup/",sid,"/6_mean_R2"), full.names = TRUE)
#filtered_files <- es_raster_files[!basename(es_raster_files) %in% "farm_mean.tif"]
# 2. Read all rasters
es_raster <- terra::rast(es_raster_files)

mean_es<-mean(es_raster)
mean_es <- project(mean_es, es_cond)
mean_es <- resample(mean_es, es_cond, method = "bilinear")


#Cost in form of ES exclusivity (minimize)
cost_syn<-terra::rast(paste0(main_dir,"WP4/cost_raster_es/",sid,"_cost_raster_es.tif"))
cost_syn <- project(cost_syn, es_cond)
cost_syn <- resample(cost_syn, es_cond, method = "bilinear")

##habitat for classification (from lulc)
habitat<-paste0(main_dir,"WP4/habitat/",sid,"_lulc.tif")
habitat<-terra::rast(habitat)

## assure that the spatial data determining optimization goals are normalized
mean_es<-min_max_normalize(mean_es)
cost_syn<-min_max_normalize(cost_syn)
es_cond<-min_max_normalize(es_cond)


cost_policy<-terra::rast(paste0(main_dir,"WP4/cost_raster_policy/",sid,"_cost_pol_sum.tif"))
cost_policy <- project(cost_policy, es_cond)
cost_policy <- resample(cost_policy, es_cond, method = "bilinear")
cost_policy[cost_policy == 0] <- NA
cost_policy<-min_max_normalize(cost_policy)# to minimize costs later
sampled_cost_pol <- terra::extract(cost_policy, grid_vect, fun = mean, na.rm = TRUE)

###############################################
############ SAMPLE GRID and OPTIM VARS #######
###############################################
sampled_es <- terra::extract(mean_es, grid_vect, fun = mean, na.rm = TRUE)
sampled_cost_syn <- terra::extract(cost_syn, grid_vect, fun = mean, na.rm = TRUE)
sampled_condition <- terra::extract(es_cond, grid_vect, fun = mean, na.rm = TRUE)

sampled_habitat <- terra::extract(habitat, grid_vect, fun = max, na.rm = TRUE)


# 6. Combine sampled data back with grid for NSGA py optimization
grid_data <- cbind(grid_clipped, sampled_es[, -1], 1/sampled_cost_pol[, -1], sampled_cost_syn[, -1], sampled_condition[, -1], sampled_habitat[, -1], pa_intersects_cell)  # remove ID column
colnames(grid_data)<-c("es_service","cost_policy","cost_es_synergy","es_condition","habitat","in_pa","geometry")

##for each cell calculate how many different lulc types they have borders with (maximize in optim)

# Step 1: Create a spatial index for fast lookup
grid_data <- grid_data %>% mutate(id = row_number())  # Add unique ID

# Step 2: Find touching neighbors (queen = FALSE => rook adjacency: cardinal directions only)
nb <- st_relate(grid_data, pattern = "F***1****")  # rook-style adjacency

# Step 3: Count neighbors with different habitat
diff_neighbors <- sapply(seq_along(nb), function(i) {
  this_habitat <- grid_data$habitat[i]
  neighbor_ids <- nb[[i]]
  sum(grid_data$habitat[neighbor_ids] != this_habitat)
})

# Step 4: Add result as a new column
grid_data$habitat_border_count <- diff_neighbors


##############################
######### subset to wetlands #####
##############################
target_lulc<-c(300:399)
out_dat<-grid_data%>%filter(habitat %in% target_lulc)
st_write(out_dat, "SK021_forest.json", driver = "GeoJSON", overwrite = T)


