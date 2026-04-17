# test prioritizr
library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(gdistance)

main_dir<-"P:/312204_pareus/"
siteID<-"SK021"
gap<-read.csv(paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap_analysis.csv"))

stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/stud_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="SK021")

pu<-st_read(paste0("WP4/2_output/02_optim/",siteID,"_input_final_grid.json"))
pu$ID<-c(1:nrow(pu))
pu$area<-as.numeric(st_area(pu))

pu<-pu%>%mutate(lock_in = case_when(
  is.na(max_IUCN_class) ~  FALSE,
  max_IUCN_class>5 ~ TRUE,
  max_IUCN_class<6 ~ FALSE
))


pu<-pu%>%mutate(inv_dist = case_when(
  min_distance_scaled == 0 ~  0,
  min_distance_scaled>0 ~ 1/min_distance_scaled
))

gap_all<-0.1-sum(gap$area_km2)/sum(gap$tot_habitat_area)

# make subsets for forest, wetlands and water
pu_for<-pu%>%filter(sampled_habitat==3)
pu_wetl<-pu%>%filter(sampled_habitat==4)
pu_wat<-pu%>%filter(sampled_habitat==5)

# ---- single LULC optim ----
p_for <-
  prioritizr::problem(pu_for, c("sampled_reg_scaled","sampled_condition_scaled"), cost_column = c("sampled_cost_scaled")) %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.0005) %>%
  #add_neighbor_constraints(k = 5) %>%
  add_relative_targets(gap[gap$sampled_habitat == 3, ]$rel_gap) %>% # target = existing plus new = 10%
  #add_absolute_targets(82)%>%
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

#### WATER ####

p_wat <-
  problem(pu_wat, c("sampled_reg","sampled_condition","inv_dist"), cost_column = "sampled_cost") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.0005) %>% # spatially clump planning units togethe
  #add_neighbor_constraints(k = 3) %>%
  add_relative_targets(gap[gap$sampled_habitat == 5, ]$rel_gap) %>%
  #add_relative_targets(0.1) %>%
  #add_absolute_targets(40)%>%
  #add_locked_out_constraints("lock_in") %>%
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
  add_boundary_penalties(penalty = 0.0005) %>%
  #add_neighbor_constraints(k = 3) %>%
  add_relative_targets(gap[gap$sampled_habitat == 4, ]$rel_gap) %>%
  #add_absolute_targets(40)%>%
  
  #add_locked_out_constraints("lock_in") %>%
  add_binary_decisions()
wetl <- solve(p_wetl)
plot(
  st_as_sf(wetl[, "solution_1"]), main = "Prioritization",
  pal = c("grey90", "darkgreen")
)

add_core_pa_wetl<-wetl%>%filter(solution_1==1)

PA_single_lulc<-rbind(forest,wetl,wat)

pu <- pu %>%
  left_join(PA_single_lulc %>%st_drop_geometry()%>% dplyr::select(ID, solution_1),
            by = "ID")
names(pu)[names(pu) == "solution_1"] <- "optim_PA_single_lulc"


pu <- pu %>%
  mutate(core_pa_lulc = case_when(
    lock_in == TRUE            ~ "existing core pa",
    optim_PA_single_lulc == 1 & n_pa == 0 ~ "new core PA",
    optim_PA_single_lulc == 1 & n_pa >0 ~ "proposed upgrade existing PA",
    TRUE                       ~ "other"
  ))

pu$core_pa_lulc <- factor(pu$core_pa_lulc,
                         levels = c(
                           "existing core pa",
                           "new core PA",
                           "proposed upgrade existing PA",
                           "other"
                         )
)

pu<-pu%>%mutate(is_core_prot_lulc = case_when(
  core_pa_lulc == "other" ~F,
  core_pa_lulc != "other" ~T
))

cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#00FF00",
  "proposed upgrade existing PA" = "#B8FFB8",
  "other" = NA
)

map_lulc <- ggplot(pu) +
  geom_sf(aes(fill = core_pa_lulc), color = NA) +
  scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))

ggsave(paste0("WP4/2_output/02_optim/",siteID,"_pa_optim.png"), plot = map_lulc, width = 18, height = 10, dpi = 300)


stats_lulc_optim<-pu%>%filter(core_pa_lulc!="other" & core_pa_lulc!="existing core pa")%>%group_by(sampled_habitat,core_pa_lulc)%>%
  summarise(ec_mean = mean(sampled_condition_scaled),
            ec_sd = sd(sampled_condition_scaled),
            km2 = sum(area)/10^6)%>%st_drop_geometry()




p_box <- ggplot(pu, aes(x = core_pa_lulc, y = LULC_class, fill = core_pa_lulc)) +
  geom_boxplot() +
  scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "none",text = element_text(size = 20))+
  labs(
    x = "",
    y = "Ecosystem condition",
    fill = ""
  )

ggsave(paste0("WP4/2_output/02_optim/",siteID,"_ec_pa_optim.png"), plot = p_box, width = 18, height = 10, dpi = 300)

kruskal.test(sampled_condition ~ core_prot, data = pu_stats)

pairwise.wilcox.test(pu_stats$sampled_condition, pu_stats$core_prot, p.adjust.method = "BH")

#save PU for OECM analysis
st_write(pu,paste0("WP4/2_output/02_optim/PA_optim_",siteID,".geojson"))
