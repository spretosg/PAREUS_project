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

area_shp<-st_area(stud_area)/10^6

PA<-st_read(paste0(main_dir,"WP4/pa_existing/PA_",sid,"/PA_WGS84.shp"))
PA<-st_transform(PA,target_crs)

#the hard Ia or II protected areas according to IUCN (unified)
PA_hard<-PA%>%filter(IUCN_CAT %in% c("Ia","II"))%>%st_union()%>%st_sf()

#centroids of all PAs
cent_PA<-sf::st_centroid(PA)

PA_union<-st_union(PA)
PA_union <- st_sf(geometry = PA_union)

# create a grid for defining the planning units PU for optimization
grid <- st_make_grid(stud_area, cellsize = 0.01, square = TRUE)
grid <- st_sf(geometry = grid)

# Clip the grid to study area
grid_clipped <- st_intersection(grid, stud_area)%>%dplyr::select()
grid_clipped<-st_make_valid(grid_clipped)



## sample all IUCN categories if itersection
grid_clipped <- st_join(grid_clipped, PA[, c("IUCN_CAT")], left = TRUE)

## attach area
grid_clipped$area<-as.numeric(st_area(grid_clipped)/10^6)



#### sample if cell is in PA (TRUE) or not (FALSE)
pa_intersects_cell<-st_intersects(grid_clipped,PA_union, sparse = FALSE)
pa_intersects_cell<-pa_intersects_cell[,1]

#### sample if cell is in hard_PA (TRUE) or not (FALSE)
pa_intersects_cell_hard<-st_intersects(grid_clipped,PA_hard, sparse = FALSE)
pa_intersects_cell_hard<-pa_intersects_cell_hard[,1]

#### and PA centroids for faster path calculation
pa_intersects_cent<-st_within(grid_clipped,cent_PA, sparse = FALSE)
pa_intersects_cent<-pa_intersects_cent[,1]

# For each cell calc closest dist do PA hard
grid_clipped$dist_pa_hard<-st_distance(grid_clipped,PA_hard)
# and distance to all kinds of PA
grid_clipped$dist_pa<-st_distance(grid_clipped,PA_union)


#Convert grid to SpatVector for terra compatibility
grid_vect <- vect(grid_clipped)

#### POPULATION grid from projects/sat-io/open-datasets/GHS/GHS_POP/GHS_POP_E2025 for least cost path lock out (not establish OECM in built-up areas) # number of people per cell (0.01km2)
pop<-terra::rast(paste0(main_dir,"WP4/lock_out/pop_",sid,".tif"))
pop <- project(pop, target_crs)

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

# cost regarding agricultural es benefits
cost_agri<-terra::rast(es_raster_files[basename(es_raster_files) %in% "farm_mean.tif"])
cost_agri <- project(cost_agri, es_cond)
cost_agri <- resample(cost_agri, es_cond, method = "bilinear")
cost_agri<-min_max_normalize(cost_agri)
sampled_cost_agri <- terra::extract(cost_agri, grid_vect, fun = mean, na.rm = TRUE)

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
sampled_pop <- terra::extract(pop, grid_vect, fun = min, na.rm = TRUE)


sampled_pa_cent<-terra::extract(pa_cent, grid_vect, fun = max, na.rm = TRUE)

# 6. Combine sampled data back with grid for NSGA py optimization
grid_data <- cbind(grid_clipped, sampled_es[, -1], sampled_cost_agri[, -1], sampled_cost_syn[, -1], sampled_condition[, -1], sampled_habitat[, -1], pa_intersects_cell,sampled_pa_cent[,-1],sampled_pop[,-1])  # remove ID column
colnames(grid_data)<-c("iucn_cat","area","dist_pa_hard","dist_pa","es_mean","cost_agri_es","cost_es_synergy","ec_mean","habitat_class","in_pa","sampled_pa_cent","population","geometry")

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
# main habitat classes
grid_data$main_habitat<-grid_data$habitat_class %/%  100
sum(grid_data$area)
st_write(grid_data, paste0(sid,"_all9.json"), driver = "GeoJSON", overwrite = T)


### gap analysis
gap_stats <- grid_data %>%st_drop_geometry()%>%
  group_by(main_habitat) %>%
  summarise(
    km2_tot = sum(area, na.rm = TRUE),
    km2_core_prot = sum(area[iucn_cat %in% c("Ia", "II")], na.rm = TRUE),
    km2_target_core_prot = 0.1*km2_tot,
    km2_gap_core_pa = km2_target_core_prot - km2_core_prot,
    km2_other_pa = sum(area[iucn_cat %in% c("III", "IV","V")], na.rm = TRUE),
    km2_target_other = 0.2*km2_tot,
    km2_gap_other_pa = km2_target_other - km2_other_pa
  )%>%filter(main_habitat %in% c(3,4,5))

# main LULC typers
required_ids <- c(3, 4, 5)

# Find which required IDs are missing
missing_ids <- setdiff(required_ids, gap_stats$main_habitat)

# Create rows with missing IDs and 0s
if (length(missing_ids) > 0) {
  new_rows <- data.frame(
    main_habitat = missing_ids,
    km2_tot = 0,
    km2_core_prot = 0,
    km2_target_core_prot = 0,
    km2_gap_core_pa = 0,
    km2_other_pa = 0,
    km2_target_other =0,
    km2_gap_other_pa = 0
  )
  
  # Append to original data frame
  gap_stats <- bind_rows(gap_stats, new_rows)
}

write.csv(gap_stats,paste0(sid,"_gap_stats.csv"))
