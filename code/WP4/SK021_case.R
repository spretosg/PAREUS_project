library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(tidyverse)
stud_site<-"SK021"
shared_dir = "data/shared"
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

x2 <- build_pa_grid(stud_site, lulc_method = "majority", out_dir = out)
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

feat_cols <- c("sampled_eco_cond_scaled", "sampled_es_reg_scaled", "sampled_cost_policy_scaled",
               "sampled_es_prov_scaled", "sampled_es_cult_scaled",
               "min_distance_scaled")

SK_feature_sets <- list(
  nature = data.frame(column = feat_cols, weight = 1,
                      direction = c("max", "max", "max", "min", "min", "min")),
  multifunctional = data.frame(column = feat_cols, weight = 1,
                               direction = c("max", "max", "max", "max", "max", "min"))
)


res_SK_glob <- run_pa_optimization(
  stud_site,
  scenario         = "nature",
  feature_sets     = SK_feature_sets,
  pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
  gap_dir          = "outputs/WP4/01_PA_analysis",
  pu_dir           = "outputs/WP4/02_optim",
  shared_dir       = "data/shared",
  out_dir          = "outputs/WP4/02_optim",
  budget_col       = "gap_area",           # or "gap_area"
  dist_transform   = "log",
  min_area_km2     = 0.075, keep_core = TRUE,
  area_weighted    = T,
  lulcs            = NULL, run_global = TRUE,
  global_lockout   = c("built-up", "agricultural"),
  boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
  mip_gap = 0.01, time_limit = 600,
  write            = T
)

res_SK_agri <- run_pa_optimization(
  stud_site,
  scenario         = "nature",
  feature_sets     = SK_feature_sets,
  pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
  gap_dir          = "outputs/WP4/01_PA_analysis",
  pu_dir           = "outputs/WP4/02_optim",
  shared_dir       = "data/shared",
  out_dir          = "outputs/WP4/02_optim",
  budget_col       = "gap_area",           # or "gap_area"
  dist_transform   = "log",
  min_area_km2     = 0.075, keep_core = TRUE,
  area_weighted    = T,
  lulcs            = NULL, run_global = TRUE,
  global_lockout   = c("built-up"),
  boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
  mip_gap = 0.01, time_limit = 600,
  write            = T
)

plot(res_SK_agri$plots$map_global)
res_SK_NOagri$log

## PA area stats:

stats<-rbind(res_SK_agri$stats%>%filter(scenario == "GLOBAL")%>%mutate(scenario = case_when(scenario == "GLOBAL" ~"AGRI_ALLOWED")),
             res_SK_glob$stats%>%filter(scenario == "GLOBAL"))

exclude_lulc = c("built-up")
d <- stats[stats$core_pa %in% pa_levels[1:3] & !stats$lulc_name %in% exclude_lulc, ]
ggplot2::ggplot(d, ggplot2::aes(x = scenario, y = km2, fill = core_pa)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::facet_wrap(~ lulc_name, scales = "free_y") +
  ggplot2::scale_fill_manual(values = pa_cols, breaks = pa_levels[1:3]) +
  ggplot2::labs(x = NULL, y = expression("Area (km"^2*")"), fill = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top", strip.background = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))




source("code/WP4/04_OECM_functions.R")
scenario<-"nature"# "nature" #"intermediate" "multifunctional"

# PUs containing the optimized core PA
pu_agri<-st_read(paste0("outputs/WP4/02_optim/PA_optim_grid_",stud_site,"_agri.geojson"))
#pu_glob<-st_read(paste0("outputs/WP4/02_optim/PA_optim_grid_",stud_site,"_global.geojson"))
#agri <- run_pca_corridors("SK021", "agri")
agri<-run_pca_corridors_grid("SK021", scenario = "agri",
                                         levels = c(low = 0.1, medium = 0.5, high = 0.8),
                                         n_rep = 3, seed = 1, k_nearest = 3,
                                         optim_dir = "outputs/WP4/02_optim",
                                         shared_dir = "data/shared",
                                         out_dir = "outputs/WP4/03_pca_landscape",
                                         write_maps = FALSE)

glob<-run_pca_corridors_grid("SK021", scenario = "global",
                             levels = c(low = 0.1, medium = 0.5, high = 0.8),
                             n_rep = 3, seed = 1, k_nearest = 3,
                             optim_dir = "outputs/WP4/02_optim",
                             shared_dir = "data/shared",
                             out_dir = "outputs/WP4/03_pca_landscape",
                             write_maps = FALSE)

plot(glob$map)

agri$corridor_km2; agri$path_km

# pies and comparison table (summarise_area() returns the area as `km2`)
stats <- dplyr::bind_rows(agri$stats, glob$stats)
area_table(stats, slice_col = "lulc_name")
plot_area_pies(stats, slice_col = "lulc_name", value_col = "km2",
               scenario_levels = c(agri$label, glob$label))


