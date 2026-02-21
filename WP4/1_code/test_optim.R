# test prioritizr
library(prioritizr)
library(sf)
library(terra)
library(dplyr)
library(gdistance)

main_dir<-"P:/312204_pareus/"
siteID<-"SK021"
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/stud_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="SK021")

pu<-st_read("SK021_optim_grid_210226a.json")
pu$ID<-c(1:nrow(pu))
pu<-pu%>%mutate(lock_in = case_when(
  is.na(max_IUCN_class) ~  FALSE,
  max_IUCN_class>5 ~ TRUE,
  max_IUCN_class<6 ~ FALSE
))
pu_for<-pu%>%filter(sampled_habitat==3)
pu_wetl<-pu%>%filter(sampled_habitat==4)
pu_wat<-pu%>%filter(sampled_habitat==5)

# plot map of planning unit coverage by protected areas
plot(st_as_sf(pu[, "sampled_reg"]), main = "Ecosytem services")
plot(st_as_sf(pu[, "sampled_condition"]), main = "Ecosytem condition")
plot(st_as_sf(pu[, "sampled_cost"]), main = "costs")
#plot(st_as_sf(pu_sel[, "habitat"]), main = "habitat")
plot(st_as_sf(pu[, "lock_in"]), main = "lock_in")

# build problem
p_for <-
  problem(pu_for, c("sampled_reg","sampled_condition"), cost_column = "sampled_cost") %>%
  add_min_set_objective() %>%
  #add_boundary_penalties(penalty = 0.005) %>%
  add_neighbor_constraints(k = 3) %>%
  add_relative_targets(0.1) %>%
  add_locked_out_constraints("lock_in") %>%
  add_binary_decisions()

# solve problem
forest <- solve(p_for)
# plot map of prioritization
plot(
  st_as_sf(forest[, "solution_1"]), main = "Prioritization",
  pal = c("grey90", "darkgreen")
)

#additional area
add_core_pa_for<-forest%>%filter(solution_1==1)

p_wat <-
  problem(pu_wat, c("sampled_reg","sampled_condition"), cost_column = "sampled_cost") %>%
  add_min_set_objective() %>%
  #add_boundary_penalties(penalty = 0.005) %>%
  add_neighbor_constraints(k = 3) %>%
  add_relative_targets(0.1) %>%
  add_binary_decisions()

wat <- solve(p_wat)
plot(
  st_as_sf(wat[, "solution_1"]), main = "Prioritization",
  pal = c("grey90", "darkgreen")
)
add_core_pa_wat<-wat%>%filter(solution_1==1)

p_wetl <-
  problem(pu_wetl, c("sampled_reg","sampled_condition"), cost_column = "sampled_cost") %>%
  add_min_set_objective() %>%
  #add_boundary_penalties(penalty = 0.005) %>%
  add_neighbor_constraints(k = 3) %>%
  add_relative_targets(0.1) %>%
  #add_locked_out_constraints("lock_in") %>%
  add_binary_decisions()
wetl <- solve(p_wetl)

add_core_pa_wetl<-wetl%>%filter(solution_1==1)

PA_new<-rbind(forest,wetl,wat)

pu <- pu %>%
  left_join(PA_new %>%st_drop_geometry()%>% dplyr::select(ID, solution_1),
            by = "ID")

## make cost raster
es_cond<-paste0(main_dir,"WP4/features/",siteID,"_ec.tif")
es_cond<-rast(es_cond)

es_cond<-project(es_cond,st_crs(pu)$wkt) 
es_cond<-raster(es_cond)
resistance<-1/es_cond

## transition raster
tr <- transition(resistance,
                 function(x) 1/mean(x),
                 directions = 8)

tr <- geoCorrection(tr, type = "c")

# make core pa poly

# 1. Keep only cells with pu > 0
cells <- pu %>% 
  filter(solution_1>0)

# 2. Union touching cells
core_pa <- cells %>%
  summarise(geometry = st_union(geometry)) %>%
  st_cast("POLYGON")%>%st_transform(st_crs(resistance))


## make core pa centroids
core_cent <- st_coordinates(st_centroid(core_pa))

#cost distance
cost_matrix <- costDistance(tr, core_cent)
cost_matrix <- as.matrix(cost_matrix)
cost_matrix[is.infinite(cost_matrix)]<-NA
# nearest n
# nearest_index <- apply(cost_matrix, 1, function(x) {
#   x[x == 0] <- NA   # ignore self-distance
#   which.min(x)
# })


# lcp
paths <- list()

# for(i in seq_along(nearest_index)) {
#   
#   j <- nearest_index[i]
#   
#   paths[[i]] <- shortestPath(tr,
#                              origin = core_cent[i,],
#                              goal   = core_cent[j,],
#                              output = "SpatialLines")
# }

#all pairs
pairs <- combn(nrow(core_cent), 2)


for(k in 1:ncol(pairs)) {
  
  i <- pairs[1, k]
  j <- pairs[2, k]
  
  # skip if unreachable
  if(is.na(cost_matrix[i, j])) next
  
  p <- shortestPath(tr,
                    origin = core_cent[i,],
                    goal   = core_cent[j,],
                    output = "SpatialLines")
  
  paths[[k]] <- p
}

paths <- paths[!sapply(paths, is.null)]

paths_sf <- st_as_sf(do.call(rbind, paths))
paths_sf <- st_set_crs(paths_sf, st_crs(pu)$wkt)

#path area

#path area for LULC
# pu_path <- pu %>%
#   st_filter(paths_sf, .predicate = st_intersects)
#number of paths
# intersection list
ints <- st_intersects(pu, paths_sf)

# count how many paths intersect each PU
pu$n_paths <- lengths(ints)
pu$n_paths_rel <- pu$n_paths/nrow(paths_sf)
q<-quantile(pu$n_paths_rel, probs = seq(0, 1, 0.25))
pu<-pu%>%mutate(imp_corridor = case_when(
  n_paths_rel>=q[4]~ T,
  n_paths_rel<q[4]~ F,
))



#stats on important corridors
stats_path <- pu %>%filter(imp_corridor==T & max_IUCN_class<6 & solution_1<1)%>%
  group_by(sampled_habitat) %>%
  summarise(area_km2 = n())%>%st_drop_geometry()

plot(st_as_sf(pu[, "n_paths_rel"]), main = "costs")


#save PU for OECM analysis
st_write(pu,paste0("WP4/2_output/1_PA_corridor_optim_",siteID,".geojson"))


