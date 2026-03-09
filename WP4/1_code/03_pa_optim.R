# test prioritizr
library(prioritizr)
library(sf)
library(terra)
library(dplyr)
library(gdistance)

main_dir<-"P:/312204_pareus/"
siteID<-"FRA_BAR2"
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/stud_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRA_BAR2")

pu<-st_read(paste0("WP4/2_output/02_optim/",siteID,"_input_grid.json"))
pu$ID<-c(1:nrow(pu))
pu<-pu%>%mutate(lock_in = case_when(
  is.na(max_IUCN_class) ~  FALSE,
  max_IUCN_class>5 ~ TRUE,
  max_IUCN_class<6 ~ FALSE
))


pu<-pu%>%mutate(inv_dist = case_when(
  min_distance_scaled == 0 ~  0,
  min_distance_scaled>0 ~ 1/min_distance_scaled
))

# make subsets for forest, wetlands and water
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
  problem(pu_for, c("sampled_reg_scaled","sampled_condition_scaled","inv_dist"), cost_column = c("sampled_cost_scaled")) %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.0005) %>%
  add_neighbor_constraints(k = 7) %>%
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
  problem(pu_wat, c("sampled_reg","sampled_condition","inv_dist"), cost_column = "sampled_cost") %>%
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
  problem(pu_wetl, c("sampled_reg","sampled_condition","inv_dist"), cost_column = "sampled_cost") %>%
  add_min_set_objective() %>%
  #add_boundary_penalties(penalty = 0.005) %>%
  add_neighbor_constraints(k = 3) %>%
  add_relative_targets(0.1) %>%
  #add_locked_out_constraints("lock_in") %>%
  add_binary_decisions()
wetl <- solve(p_wetl)
plot(
  st_as_sf(wetl[, "solution_1"]), main = "Prioritization",
  pal = c("grey90", "darkgreen")
)

add_core_pa_wetl<-wetl%>%filter(solution_1==1)

PA_new<-rbind(forest,wetl,wat)

pu <- pu %>%
  left_join(PA_new %>%st_drop_geometry()%>% dplyr::select(ID, solution_1),
            by = "ID")

pu <- pu %>%
  mutate(fill_col = case_when(
    lock_in == TRUE            ~ "existing core pa",
    solution_1 == 1            ~ "new core pa",
    TRUE                       ~ "other"
  ))

p<-ggplot(pu) +
  geom_sf(aes(fill = fill_col), color = NA) +
  scale_fill_manual(
    values = c(
      "existing core pa" = "orange",
      "new core pa" = "green",
      "other" = "white"
    ),
    name = NULL
  ) +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()
ggsave(paste0("WP4/2_output/02_optim/",siteID,"_optim_corPA.png"), plot = p, width = 8, height = 6, dpi = 300)


stats_new_pa<-pu%>%st_drop_geometry()%>%dplyr::filter(is.na(max_IUCN_class) & solution_1>=1)%>%group_by(sampled_habitat)%>%
  summarise(area = n())

####----2. Connect PA ####
#cost raster
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



# lcp
paths <- list()

## closest pair
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

pu<-pu%>%mutate(prot_class = case_when(
  is.na(max_IUCN_class) ~ "no_prot",
  max_IUCN_class<6 | solution_1<1 ~ "weak",
  max_IUCN_class > 5 | solution_1>1 ~ "core pa"
))

#stats on important corridors
stats_path <- pu %>%filter(imp_corridor==T)%>%
  group_by(sampled_habitat, prot_class) %>%
  summarise(area_km2 = n())%>%st_drop_geometry()

p<-ggplot(pu %>% filter(n_paths_rel > 0.02)%>%dplyr::select(n_paths_rel)) +
  geom_sf(aes(fill = n_paths_rel), color = NA) +
  labs(title = "Important corridors",
       fill = "n_paths_rel") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()
ggsave(paste0("WP4/2_output/02_optim/",siteID,"_imp_paths.png"), plot = p, width = 8, height = 6, dpi = 300)


#save PU for OECM analysis
st_write(pu,paste0("WP4/2_output/02_optim/PA_corridor_optim_",siteID,".geojson"))


