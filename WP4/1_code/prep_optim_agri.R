## test data read tif
## make pu of stud area
## sample cost, es_cond and es_capacity per PU
## save as json upload gcs bucket for python

library(sf)
library(terra)
library(dplyr)
source("WP4/1_code/wp4_functions_utils.R")

########################################
############ GENERAL ###################
########################################

main_dir<-"P:/312204_pareus/"
sid<-"SK021"
stud_area<-read_sf(paste0(main_dir,"WP2/PGIS_ES_mapping/raw_data_backup/",sid,"/study_site.gpkg"))%>%filter(siteID %in% sid)
target_crs <- st_crs(stud_area)$wkt

PA<-st_read(paste0(main_dir,"WP4/pa_existing/PA_",sid,"/PA_WGS84.shp"))
PA<-st_transform(PA,target_crs)
#PA<-PA[1:3,]


cent_PA<-sf::st_centroid(PA)



PA<-st_union(PA)
PA <- st_sf(geometry = PA)

# create a grid for defining the planning units PU for optimization
grid <- st_make_grid(stud_area, cellsize = 0.1, square = TRUE)
grid <- st_sf(geometry = grid)
# Clip the grid to polygon boundary
grid_clipped <- st_intersection(grid, stud_area)%>%dplyr::select()
grid_clipped<-st_make_valid(grid_clipped)

#### sample if cell is in PA (TRUE) or not (FALSE)
pa_intersects_cell<-st_intersects(grid_clipped,PA, sparse = FALSE)
pa_intersects_cell<-pa_intersects_cell[,1]

#### and PA centroids for faster path calculation
# pa_intersects_cent<-st_within(grid_clipped,cent_PA, sparse = FALSE)
# pa_intersects_cent<-pa_intersects_cent[,1]

# Spatial join to find which grid cells contain any points



#Convert grid to SpatVector for terra compatibility
grid_vect <- vect(grid_clipped)

#### POPULATION grid from projects/sat-io/open-datasets/GHS/GHS_POP/GHS_POP_E2025 for least cost path lock out (not establish OECM in built-up areas) # number of people per cell (0.01km2)
lock_out_pop<-terra::rast(paste0(main_dir,"WP4/lock_out/pop_",sid,".tif"))
lock_out_pop <- project(lock_out_pop, target_crs)
## and a reclassified raster for the barriers in the least cost path
n_target = 6 #4 people per 0.01km2
m <- c(0, n_target, FALSE,
       n_target, max(values(lock_out_pop,na.rm=T)), TRUE,NaN,NaN,FALSE)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
lock_out<- terra::classify(lock_out_pop,rclmat,include.lowest = T)
plot(lock_out)



######################
#### OPTIM VARS #######
#######################

#### es condition from INRAE (maximize)
es_cond<-paste0(main_dir,"WP2/es_condition/",sid,"_es_cond.tif")
es_cond<-terra::rast(es_cond)
es_cond <- project(es_cond, target_crs)

#rasterize
r_rasterized <- stars::st_rasterize(cent_PA)
r_rasterized<-terra::rast(r_rasterized)
pa_cent<-r_rasterized$WDPAID_WDPAID
m <- c(1, 150000, TRUE,
       0, 1, FALSE,NaN,NaN,FALSE)
rclmat1 <- matrix(m, ncol=3, byrow=TRUE)
pa_cent<- terra::classify(pa_cent,rclmat1,include.lowest = T)


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

####################################
########## AGRI (LULC 200-299) #####
####################################

cost_es<-terra::rast(es_raster_files[basename(es_raster_files) %in% "farm_mean.tif"])
cost_es <- project(cost_es, es_cond)
cost_es <- resample(cost_es, es_cond, method = "bilinear")
cost_es<-min_max_normalize(cost_es)
sampled_cost_es <- terra::extract(cost_es, grid_vect, fun = mean, na.rm = TRUE)

###############################################
############ SAMPLE GRID and OPTIM VARS #######
###############################################
sampled_es <- terra::extract(mean_es, grid_vect, fun = mean, na.rm = TRUE)
sampled_cost_syn <- terra::extract(cost_syn, grid_vect, fun = mean, na.rm = TRUE)
sampled_condition <- terra::extract(es_cond, grid_vect, fun = mean, na.rm = TRUE)
#take min to account for small patches as well
sampled_habitat <- terra::extract(habitat, grid_vect, fun = min, na.rm = TRUE)


###############################################
############ SAMPLE GRID LOCK OUT VARS ########
###############################################
sampled_pop <- terra::extract(lock_out_pop, grid_vect, fun = min, na.rm = TRUE)
sampled_lock_out <- terra::extract(lock_out, grid_vect, fun = max, na.rm = TRUE)

sampled_pa_cent<-terra::extract(pa_cent, grid_vect, fun = max, na.rm = TRUE)

# 6. Combine sampled data back with grid for NSGA py optimization
grid_data <- cbind(grid_clipped, sampled_es[, -1], sampled_cost_es[, -1], sampled_cost_syn[, -1], sampled_condition[, -1], sampled_habitat[, -1], pa_intersects_cell,sampled_pa_cent[,-1],sampled_pop[,-1], sampled_lock_out[,-1])  # remove ID column
colnames(grid_data)<-c("es_service","cost_agri_es","cost_es_synergy","es_condition","habitat","in_pa","sampled_pa_cent","population","lock_out_pop" ,"geometry")

##for each cell calculate how many different lulc types they have borders with (maximize in optim)

# Step 1: Create a spatial index for fast lookup
grid_data <- grid_data %>% mutate(id = row_number())  # Add unique ID

# 
nb <- st_relate(grid_data, pattern = "F***T****")  # queen adjacency

# Step 3: Count neighbors with different habitat
diff_neighbors <- sapply(seq_along(nb), function(i) {
  this_habitat <- grid_data$habitat[i]
  neighbor_ids <- nb[[i]]
  sum(grid_data$habitat[neighbor_ids] != this_habitat)
})

# Step 4: Add result as a new column
grid_data$habitat_border_count <- diff_neighbors


# ### add queen adjecency to PA
# diff_neighbors_pa <- sapply(seq_along(nb), function(i) {
#   if (!grid_data$in_pa[i]) {
#     return(0)
#   }
#   pa <- grid_data$in_pa[i]
#   neighbor_ids <- nb[[i]]
#   sum(grid_data$in_pa[neighbor_ids] != pa)
# })
# 
# grid_data$pa_border <- diff_neighbors_pa

st_write(grid_data, "SK021_all_low2.json", driver = "GeoJSON", overwrite = T)

##############################
######### subset to agri #####
##############################
target_lulc<-c(200:299)
out_dat<-grid_data%>%filter(habitat %in% target_lulc)
st_write(out_dat, "SK021_agri.json", driver = "GeoJSON", overwrite = T)


