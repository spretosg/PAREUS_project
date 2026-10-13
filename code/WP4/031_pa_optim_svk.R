### PA_optim SVK

library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
source("code/WP4/wp4_functions_utils.R")
target_site<-"SK021"

scenario<-"nature"# "nature" "multifunctional"

#the low nature scenario: maximizes condition and regulation and simultaneously maximizes for ecosystem services
#intermediate just maximizes condition and regulation minimizes policy
#conflicting: maximizes condition and regulation, minimizes ecosystem services

#select optimization features according to the scenario
if(scenario == "multifunctional"){
  features<-c("sampled_eco_cond_scaled","sampled_es_reg_scaled","sampled_cost_policy_scaled","sampled_es_prov_scaled","sampled_es_cult_scaled","inv_dist","inv_ghm")
}else{
  features<-c("sampled_eco_cond_scaled","sampled_es_reg_scaled","sampled_cost_policy_scaled","inv_prov","inv_cult","inv_dist","inv_ghm")
  
}


gap<-read.csv(paste0("outputs/WP4/01_PA_analysis/",target_site,"_gap_analysis.csv"))

stud_area<-read_sf(paste0("data/shared/",target_site,".gpkg"))
# planning units
pu<-st_read(paste0("outputs/WP4/02_optim/",target_site,"_input_final_grid.json"))

## as costs policy reflects policy coherence 1-coherence = policy friction, fiction needs to be minimizede thus
#pu$sampled_cost_policy

# the others
pu$inv_prov<-1-pu$sampled_es_prov_scaled
pu$inv_cult<-1-pu$sampled_es_cult_scaled
pu$inv_ghm<-1-pu$sampled_ghm_scaled
#scale distance to pa
pu<-zero_one_scale(pu,cols = "min_distance")
pu$inv_dist<-1-pu$min_distance_scaled

pu$area_km2<-pu$area/10^6
pu<-pu%>%filter(area_km2 > 0.28)

##neighbor cells to core_pa
# Identify core cells
core <- pu[pu$existing_corePA == TRUE, ]

# Buffer core cells by 5 km (use a projected CRS in metres)
core_buffer <- st_buffer(core, dist = 5000) |>
  st_union()

# Select cells within 5 km of core cells, excluding core cells
buffer <- pu[
  lengths(st_intersects(pu, core_buffer)) > 0 &
    pu$existing_corePA != TRUE,
]





### SCEN 1 nature global 10%
core_pa_global<-pa_optim(pu=pu,lockout = c("built-up"),lockin_col="existing_corePA", features =features, 
                         area_budget = gap[gap$lulc_name == "global_0.1_prot", ]$target_area/10^6, boundary_penalty = 5e-4)

core_pa_global_agri<-pa_optim(pu=pu,lockout = c("built-up","agricultural"),lockin_col="existing_corePA", features =features, 
                         area_budget = gap[gap$lulc_name == "global_0.1_prot", ]$target_area/10^6, boundary_penalty = 5e-4)



pu <- pu %>%
  left_join(core_pa_global %>%st_drop_geometry()%>% dplyr::select(id, solution_1),
            by = "id")
names(pu)[names(pu) == "solution_1"] <- "optim_PA_global"

pu <- pu %>%
  left_join(core_pa_global_agri %>%st_drop_geometry()%>% dplyr::select(id, solution_1),
            by = "id")
names(pu)[names(pu) == "solution_1"] <- "optim_PA_global_agri"

## classify the different areas old and new core PAs
pa_levels <- c("existing core pa", "new core PA", "proposed upgrade existing PA",
               "other protected areas", "not protected")


pu$optim_PA_global      <- pu$id %in% core_pa_global$id
pu$optim_PA_global_agri <- pu$id %in% core_pa_global_agri$id   # rename to ..._noagri

pu$core_pa_global <- classify_pa(pu, core_pa_global$id)
pu$core_pa_agri   <- classify_pa(pu, core_pa_global_agri$id)    # rename to core_pa_noagri

pu$is_core_prot_global <- pu$core_pa_global %in% pa_levels[1:3]
pu$is_core_prot_agri   <- pu$core_pa_agri   %in% pa_levels[1:3]

# Sanity check per scenario
sapply(list(global = core_pa_global, noagri = core_pa_global_agri),
       function(x) if (is.null(x)) 0 else sum(x$area_km2))


st_write(pu,paste0("outputs/WP4/02_optim/PA_optim_grid_",target_site,"_",scenario,".geojson"))

cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#bf8bff",
  "proposed upgrade existing PA" = "#e5d0ff",
  "other" = NA
)

map_pa <- ggplot(pu) +
  geom_sf(aes(fill = core_pa_global), color = NA) +
  scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  geom_sf(data=pu%>%filter(lulc_name == "built-up"),aes(fill = lulc_name),
          color = "grey", fill="grey" )+
  theme_minimal()+
  theme(
    legend.position = "top",
    text = element_text(size = 15)
  )+ggtitle("Including agricultural areas")




stats_agri_optim<-pu%>%mutate(core_pa = core_pa_agri)%>%group_by(lulc_name,core_pa)%>%
  summarise(ec_mean = mean(sampled_eco_cond_scaled),
            ec_sd = sd(sampled_eco_cond_scaled),
            km2 = sum(area_km2),
            reg_mean = mean(sampled_es_reg_scaled),
            reg_sd = sd(sampled_es_reg_scaled),
            pol_cost_mean = mean(sampled_cost_policy),
            cost_sd = sd(sampled_cost_policy),
            dist_mean = mean(min_distance),
            dist_sd = sd(min_distance))%>%st_drop_geometry() %>%mutate(scenario = "AGRI")

stats_global_optim<-pu%>%mutate(core_pa = core_pa_global)%>%group_by(lulc_name,core_pa)%>%
  summarise(ec_mean = mean(sampled_eco_cond_scaled),
            ec_sd = sd(sampled_eco_cond_scaled),
            km2 = sum(area_km2),
            reg_mean = mean(sampled_es_reg_scaled),
            reg_sd = sd(sampled_es_reg_scaled),
            pol_cost_mean = mean(sampled_cost_policy),
            cost_sd = sd(sampled_cost_policy),
            dist_mean = mean(min_distance),
            dist_sd = sd(min_distance))%>%st_drop_geometry()%>%mutate(scenario = "GLOBAL")


stats<-rbind(stats_agri_optim,stats_global_optim)

#area statistics
cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#bf8bff",
  "proposed upgrade existing PA" = "#e5d0ff",
  "other protected areas"        = "#ADD8E6",
  "not protected"                = "#1A1A1A"
)
area_plot<-ggplot(stats%>%filter(core_pa != "other protected areas" & !is.na(lulc_name) & !lulc_name %in% c("built-up") & core_pa != "not protected" ),
                  aes(x = scenario,
                      y = km2,
                      fill = core_pa)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(~ lulc_name, scales = "free_y") +
  scale_fill_manual(values = cols)+
  labs(
    x = NULL,
    y = expression("Area (km"^2*")"),
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
ggsave(paste0("outputs/WP4/02_optim/",target_site,"_",scenario,"_pa_area_statsProt.png"), plot = area_plot, width = 13, height = 10, dpi = 300)







### SCEN 2 nature 10% within buffer zone

##check available land first
# cand <- buffer[!(buffer$lulc_name %in% c("built-up")) &
#                      !buffer$existing_corePA, ]
# sum(cand$area_km2)
# gap[gap$lulc_name == "global_0.1_prot", ]$target_area / 1e6
# 
# c(available = avail_km2, budget = budget_km2)
# 
# 
# core_pa_neigh<-pa_optim(pu=buffer,lockout = c("built-up"),lockin_col="existing_corePA", features =features, 
#                          area_budget = gap[gap$lulc_name == "global_0.1_prot", ]$target_area/10^6, boundary_penalty = 5e-4)


