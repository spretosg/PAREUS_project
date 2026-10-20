##FRL04 and FRA_BAR are the smaller areas
library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(tidyverse)
library(ggnewscale)
## FRA is the department

## calculate PA and OECM globally (not by LULC) for all three. Crop the FRA with BAR and FRL and count area stats (LULC by protection status)
## compare the regional protection with the two stats at local level

#FRA
stud_site<-"FRA_BAR2"
shared_dir = "data/shared"
pu_size = 800
stud_area <- sf::read_sf(file.path(shared_dir, paste0(stud_site, ".gpkg")))

source("code/WP4/01_pa_analysis_funtions.R")
## ---- Defaults ---------------------------------------------------------------
default_iucn_rank <- c(
  "Not Applicable" = 1, "Not Assigned" = 1, "Not Reported" = 1,
  "VI" = 2, "V" = 3, "IV" = 4, "III" = 5, "II" = 6, "Ib" = 7, "Ia" = 8
)
default_lulc_names <- c("1" = "built-up", "2" = "agricultural", "3" = "forest",
                        "4" = "wetland", "5" = "water")
default_iucn_colours <- c(
  "Ia" = "#00A300", "Ib" = "#007A00", "II" = "#00FF00",
  "III" = "#804FB3", "IV" = "#9969C7", "V" = "#B589D6", "VI" = "#CDB4E6",
  "Not Applicable" = "#F1EE8D", "Not Assigned" = "#F1EE8D", "Not Reported" = "#F1EE8D"
)

out <- "outputs/WP4/01_PA_analysis"

x2 <- build_pa_grid(stud_site, lulc_method = "majority", pu_size =  pu_size, out_dir = out)
res <- analyse_pa_gap(x2, global_target = 0.1,
                      lulc_targets = c(forest = 0.1, water = 0.1, wetland = 0.1, agricultural = 0),
                      out_dir = out)


#### 2. prep optimization
source("code/WP4/02_prep_optim_funtions.R")
op <- prepare_optim_grid(stud_site, na_action = "neighbour_mean", skip_missing = TRUE,     layer_files  = c(
  eco_cond    = "data/WP2/{site}/int.tif",
  cost_policy = "data/WP3/pol_resistance_{site}.tif"
))


### 3. optim PA
source("code/WP4/03_pa_optim_functions.R")
### NOR case LULC vs global

feat_cols <- c("sampled_eco_cond_scaled", "sampled_es_reg_scaled", "sampled_cost_policy_scaled",
               "sampled_es_prov_scaled", "sampled_es_cult_scaled",
               "min_distance_scaled")

FRA_feature_sets <- list(
  nature = data.frame(column = feat_cols, weight = 1,
                      direction = c("max", "max", "max", "min", "min", "min")),
  multifunctional = data.frame(column = feat_cols, weight = 1,
                               direction = c("max", "max", "max", "max", "max", "min"))
)


res_FR <- run_pa_optimization(
  stud_site,
  scenario         = "nature",
  feature_sets     = FRA_feature_sets,
  pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
  gap_dir          = "outputs/WP4/01_PA_analysis",
  pu_dir           = "outputs/WP4/02_optim",
  shared_dir       = "data/shared",
  out_dir          = "outputs/WP4/02_optim",
  budget_col       = "gap_area",           # or "gap_area"
  dist_transform   = "log",
  min_area_km2     = 0.001, keep_core = TRUE,
  area_weighted    = T,
  lulcs            = c("water","forest"), run_global = TRUE,
  global_lockout   = c("built-up", "agricultural"),
  boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
  mip_gap = 0.01, time_limit = 600,
  write            = TRUE
)
plot(res_FR$plots$area)

plot(res_FR$plots$map_global)


source("code/WP4/04_OECM_functions.R")
scenario<-"nature"# "nature" #"intermediate" "multifunctional"
####---- User parameter ----####
#MCDA weighting 
human_c<-list(
  w_cult = 0.4,
  w_prov =0.4,
  w_reg = 0.1,
  w_connect = 0.1
)

coverage <- c(
  forest = 0.2,
  water = 0.2,
  wetland = 0.2,
  agricultural = 0.2
)

# PUs containing the optimized core PA
pu<-st_read(paste0("outputs/WP4/02_optim/PA_optim_grid_",stud_site,"_",scenario,".geojson"))

#connectivity raster (run 00_helper_connectivity.R first)
connectivity<-rast(paste0("outputs/WP4/02_optim/conn_",stud_site,".tif"))

# sample connetivity values to pu
pu$connectivity<- terra::extract(connectivity$normalized_current, pu, fun = mean, na.rm = TRUE)[,2]
pu<-zero_one_scale(
  pu,
  cols = c("connectivity")
)


#human centric
pu$oecm_suit<-oecm_lin_w(pu$sampled_es_cult_scaled,pu$sampled_es_prov_scaled,
                         pu$sampled_es_reg_scaled,pu$connectivity_scaled,human_c$w_cult,human_c$w_prov,human_c$w_reg,human_c$w_connect)


oecm_global<-select_oecm(pu=pu, 
                         mode= "global", 
                         coverage = 0.2,
                         suitability_col = "oecm_suit",
                         corePA_col = "core_pa_global",
                         search_oecm_in = c("not protected","other protected areas"),
                         w= 0.9, order = 2)

base_glob<-plot_pca_map(pu=pu,
                        corePA_col = "core_pa_global",
                        oecm_df=oecm_global,
                        stud_area,
                        scen = "GLOB_10",
                        save_output = T)


pca_labels <- c(
  "Other PA (IUCN III-VI) high suitability for OECM"     = "Other PA (IUCN III-VI)",
  "other_pa"                                             = "Other PA (IUCN III-VI)",
  "core_PA"                                              = "Optimized IUCN Ia or II",
  "not_protected"                                        = "No protection",
  "High OECM suitability not protected - potential OECM" = "Potential OECM"
)




#### compare the two local cases with the regional case:

BAR<-"FRA_BAR2"
FRL04<-"FRL04"
region<-"FRA"
shared_dir = "data/shared"



pu_reg<-sf::read_sf(paste0("outputs/WP4/03_pca_landscape/PCA_PU_",region,"_GLOB_10.geojson"))
pu_BAR<-sf::read_sf(paste0("outputs/WP4/03_pca_landscape/PCA_PU_",BAR,"_GLOB_10.geojson"))
pu_FRL04<-sf::read_sf(paste0("outputs/WP4/03_pca_landscape/PCA_PU_",FRL04,"_GLOB_10.geojson"))


area_BAR <- sf::read_sf(file.path(shared_dir, paste0(BAR, ".gpkg")))%>%st_transform(st_crs(pu_reg))
area_FRL04 <- sf::read_sf(file.path(shared_dir, paste0(FRL04, ".gpkg")))%>%st_transform(st_crs(pu_reg))

## intersect the PU reg with BAR and FRL for regionally optimized BAR and FRL
pu_reg <- pu_reg %>%
  mutate(
    in_BAR = lengths(st_intersects(geometry, area_BAR)) > 0,
    in_FRL04 = lengths(st_intersects(geometry, area_FRL04)) > 0
  )

pu_reg_BAR<-pu_reg%>%filter(in_BAR == T)
pu_reg_FRL<-pu_reg%>%filter(in_FRL04 == T)

stats_FRL_reg<-pu_reg_FRL%>%group_by(lulc_name, pca)%>%
  st_drop_geometry()%>%summarise(area = sum(area_km2))%>%mutate(optim_scenario = "regional")
stats_FRL_local<-pu_FRL04%>%group_by(lulc_name, pca)%>%
  st_drop_geometry()%>%summarise(area = sum(area_km2))%>%mutate(optim_scenario = "local")
stats_FRL<-rbind(stats_FRL_reg,stats_FRL_local)


stats_BAR_reg<-pu_reg_BAR%>%group_by(lulc_name, pca)%>%
  st_drop_geometry()%>%summarise(area = sum(area_km2))%>%mutate(optim_scenario = "regional")
stats_BAR_local<-pu_BAR%>%group_by(lulc_name, pca)%>%
  st_drop_geometry()%>%summarise(area = sum(area_km2))%>%mutate(optim_scenario = "local")
stats_BAR<-rbind(stats_BAR_reg,stats_BAR_local)

plot_df <- stats_BAR %>%
  group_by(optim_scenario, pca, lulc_name) %>%
  summarise(km2 = sum(area, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    pca = factor(pca),
    lulc_name = factor(lulc_name)
  )

ggplot(plot_df%>%filter(pca %in% c("Other PA (IUCN III-VI) high suitability for OECM",
                                   "new core PA","proposed upgrade existing PA")&lulc_name !="built-up"), aes(
  x = optim_scenario,
  y = km2,
  fill = pca
)) +
  geom_col(width = 0.75) +
  facet_wrap(~lulc_name, ncol = 3,
             labeller = labeller(
               pca = function(x) paste("Land-cover type", x)
             )) +
  scale_fill_manual(
    values = c(    "existing core pa" = "#00A300",
                   "new core PA" = "#bf8bff",
                   "proposed upgrade existing PA" = "#e5d0ff",
                   "Other PA (IUCN III-VI) high suitability for OECM" = "#44a6c6",
                   "High OECM suitability not protected - potential OECM" = "#194553",
                   "Other PA (IUCN III-VI)" = "#ADD8E6"),
    name = "Protection"
  ) +
  scale_y_continuous(
    name = expression("Area (" * km^2 * ")"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Optimisation scenario") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(forcats)
# 
# # Summarise area by scenario and category
# plot_df <- stats_FRL %>%
#   group_by(optim_scenario, lulc_name, pca) %>%
#   summarise(km2 = sum(area, na.rm = TRUE), .groups = "drop") %>%
#   mutate(
#     lulc_name = as.character(lulc_name),
#     pca = factor(pca, levels = 1:5)
#   ) %>%
#   pivot_wider(
#     names_from = optim_scenario,
#     values_from = km2,
#     values_fill = 0
#   )
# 
# plot_df <- bind_rows(lulc, protection) %>%
#   mutate(group = factor(
#     group,
#     levels = c("Land cover", "Protection class")
#   ))
# 
# plot_pairs <- plot_df %>%
#   pivot_wider(names_from = optim_scenario, values_from = km2)
# 
# ggplot(plot_pairs, aes(y = fct_reorder(category, local))) +
#   geom_segment(
#     aes(x = local, xend = regional,
#         yend = category),
#     colour = "grey75",
#     linewidth = 1
#   ) +
#   geom_point(aes(x = local, colour = "Local optim"),
#              size = 3.2) +
#   geom_point(aes(x = regional, colour = "Regional optim"),
#              size = 3.2) +
#   facet_wrap(~group, scales = "free_y", ncol = 2) +
#   scale_colour_manual(
#     values = c("Local optim" = "#2878B5",
#                "Regional optim" = "#E07A3F")
#   ) +
#   scale_x_continuous(
#     name = expression("Area (" * km^2 * ")"),
#     expand = expansion(mult = c(0.03, 0.10))
#   ) +
#   labs(y = NULL, colour = NULL) +
#   theme_minimal(base_size = 12) +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.text = element_text(face = "bold", size = 12),
#     legend.position = "top",
#     axis.text.y = element_text(colour = "grey20"),
#     panel.spacing = grid::unit(1.5, "lines")
#   )
