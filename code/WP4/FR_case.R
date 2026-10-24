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
area_FRA <- sf::read_sf(file.path(shared_dir, paste0(region, ".gpkg")))%>%st_transform(st_crs(pu_reg))


### plot IUCN region with counties
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

a<-ggplot2::ggplot() +
  ggplot2::geom_sf(data = pu_reg, ggplot2::aes(fill = .data$IUCN_CAT), colour = NA) +
  ggplot2::scale_fill_manual(values = default_iucn_colours, na.value = NA, name = "IUCN") +
  ggplot2::geom_sf(data = pu_reg[pu_reg$lulc_name =="built-up", ], fill = "grey", colour = NA) +
  ggplot2::geom_sf(data = area_FRA, fill = NA, colour = "black") +
  ggplot2::geom_sf(data = area_BAR, fill = NA, colour = "red",lty=2,linewidth=1.7) +
  ggplot2::geom_sf(data = area_FRL04, fill = NA, colour = "red", lty= 2, linewidth =1.7) +
  ggplot2::theme_minimal()

out <- "outputs/WP4/01_PA_analysis"
ggplot2::ggsave(file.path(out, "FRA_all_iucn.png"), a,
                width = 8, height = 6, dpi = 300)

## gap for all three areas core PA to reach 10% as tables

## PU optim global with counties
pa_levels <- c("existing core pa", "new core PA", "proposed upgrade existing PA",
               "other protected areas", "not protected")
pa_cols <- c("existing core pa" = "#00A300", "new core PA" = "#bf8bff",
             "proposed upgrade existing PA" = "#e5d0ff",
             "other protected areas" = "#ADD8E6", "not protected" = "#1A1A1A")
show <- c(pa_cols[1:3], "other protected areas" = "transparent", "not protected" = "transparent")
class_col="core_pa_global"
a<-ggplot2::ggplot() +
  ggplot2::geom_sf(data = pu_reg, ggplot2::aes(fill = .data[[class_col]]), colour = NA) +
  ggplot2::geom_sf(data = pu_reg[pu_reg$lulc_name =="built-up", ], fill = "grey", colour = NA) +
  ggplot2::geom_sf(data = area_FRA, fill = NA, colour = "black") +
  ggplot2::scale_fill_manual(values = show, breaks = names(pa_cols)[1:3], name = NULL) +
  ggplot2::ggtitle("Regional optimization") +
  ggplot2::geom_sf(data = area_BAR, fill = NA, colour = "red",lty=2,linewidth=1.7) +
  ggplot2::geom_sf(data = area_FRL04, fill = NA, colour = "red", lty= 2, linewidth =1.7) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "top", text = ggplot2::element_text(size = 15))

out <- "outputs/WP4/02_optim"
ggplot2::ggsave(file.path(out, "FRA_regional_pa_optim.png"), a,
                width = 8, height = 6, dpi = 300)


### OECM stats
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
stats_FRL$study_site<-"FRL04"


stats_BAR_reg<-pu_reg_BAR%>%group_by(lulc_name, pca)%>%
  st_drop_geometry()%>%summarise(area = sum(area_km2))%>%mutate(optim_scenario = "regional")
stats_BAR_local<-pu_BAR%>%group_by(lulc_name, pca)%>%
  st_drop_geometry()%>%summarise(area = sum(area_km2))%>%mutate(optim_scenario = "local")
stats_BAR<-rbind(stats_BAR_reg,stats_BAR_local)
stats_BAR$study_site<-"BAR"

stats_all<-rbind(stats_BAR,stats_FRL)

plot_df <- stats_all %>% filter(lulc_name != "built-up")%>%
  group_by(optim_scenario, pca,study_site) %>%
  summarise(km2 = sum(area, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    pca = factor(pca)
  )

#"Other PA (IUCN III-VI) high suitability for OECM",
ggplot(
  plot_df %>%
    filter(
      pca %in% c("new core PA", "proposed upgrade existing PA") 
    ),
  aes(
    x = optim_scenario,
    y = km2,
    fill = pca
  )
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = round(km2, 1)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  facet_wrap(~study_site, ncol = 3) +
  scale_fill_manual(
    values = c(
      "existing core pa" = "#00A300",
      "new core PA" = "#bf8bff",
      "proposed upgrade existing PA" = "#e5d0ff",
      "Other PA (IUCN III-VI) high suitability for OECM" = "#44a6c6",
      "High OECM suitability not protected - potential OECM" = "#194553",
      "Other PA (IUCN III-VI)" = "#ADD8E6"
    ),
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


### PCA regional with subareas
a<-ggplot() +
  # Layer 0: other pa
  geom_sf(data = pu_reg%>%filter(pca != "not protected"),
          aes(fill = pca),
          color = NA) +
  scale_fill_manual(
    values = c(    "existing core pa" = "#00A300",
                   "new core PA" = "#bf8bff",
                   "proposed upgrade existing PA" = "#e5d0ff",
                   "Other PA (IUCN III-VI) high suitability for OECM" = "#44a6c6",
                   "High OECM suitability not protected - potential OECM" = "#194553",
                   "Other PA (IUCN III-VI)" = "#ADD8E6"),
    name = NULL
  )  +
  geom_sf(
    data = pu_reg%>%filter(lulc_name == "built-up"),
    aes(fill = lulc_name),
    color = "grey", fill="grey"  )+
  
  new_scale_fill() +
  # stud area
  geom_sf(data = area_FRA, fill = NA, color = "black") +
  ggplot2::geom_sf(data = area_BAR, fill = NA, colour = "red",lty=2,linewidth=1.7) +
  ggplot2::geom_sf(data = area_FRL04, fill = NA, colour = "red", lty= 2, linewidth =1.7) +
  theme_minimal()+
  ggplot2::theme(legend.position = "top", text = ggplot2::element_text(size = 15))

out <- "outputs/WP4/03_pca_landscape"
ggplot2::ggsave(file.path(out, "FRA_regional_pca_optim.png"), a,
                width = 8, height = 6, dpi = 300)


## pie plot
fill_values = c(  "existing core pa" = "#00A300",
                  "new core PA" = "#bf8bff",
                  "proposed upgrade existing PA" = "#e5d0ff",
                  "Other PA (IUCN III-VI) high suitability for OECM" = "#44a6c6",
                  "High OECM suitability not protected - potential OECM" = "#194553",
                  "Other PA (IUCN III-VI)" = "#ADD8E6",
                  "not protected" = "#999999"
)
source("code/WP4/04_OECM_functions.R")
ggplot2::ggsave("outputs/WP4/03_pca_landscape/TRD_pie_pca.png",
                plot_area_pies(stats_all, value_col       = "area",
                               scenario_col    = "optim_scenario", slice_col = "pca", label = "pct", keep_slices = NULL, 
                               facet_col= "study_site", fill_values = fill_values),
                width = 11, height = 6, dpi = 300)