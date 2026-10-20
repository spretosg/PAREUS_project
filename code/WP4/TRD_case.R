#### PCA TRD

library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(tidyverse)
library(ggnewscale)
stud_site<-"TRD"
shared_dir = "data/shared"
pu_size = 350
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

TRD_feature_sets <- list(
  nature = data.frame(column = feat_cols, weight = 1,
                      direction = c("max", "max", "max", "min", "min", "min")),
  multifunctional = data.frame(column = feat_cols, weight = 1,
                               direction = c("max", "max", "max", "max", "max", "min"))
)


res_TRD <- run_pa_optimization(
  stud_site,
  scenario         = "nature",
  feature_sets     = TRD_feature_sets,
  pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
  gap_dir          = "outputs/WP4/01_PA_analysis",
  pu_dir           = "outputs/WP4/02_optim",
  shared_dir       = "data/shared",
  out_dir          = "outputs/WP4/02_optim",
  budget_col       = "gap_area",           # or "gap_area"
  dist_transform   = "log",
  min_area_km2     = 0.001, keep_core = TRUE,
  area_weighted    = T,
  lulcs            = c("water","forest","wetland"), run_global = TRUE,
  global_lockout   = c("built-up", "agricultural"),
  boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
  mip_gap = 0.01, time_limit = 600,
  write            = TRUE
)
plot(res_TRD$plots$area)


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



oecm_lulc<-select_oecm(pu=pu, 
                       mode= "class", 
                       coverage = coverage,
                       lulc_col = "lulc_name",
                       suitability_col = "oecm_suit",
                       corePA_col = "core_pa_lulc",
                       search_oecm_in = c("not protected","other protected areas"),
                       w= 0.9, order = 2)

oecm_global<-select_oecm(pu=pu, 
                         mode= "global", 
                         coverage = 0.2,
                         suitability_col = "oecm_suit",
                         corePA_col = "core_pa_global",
                         search_oecm_in = c("not protected","other protected areas"),
                         w= 0.9, order = 2)

base_lulc<-plot_pca_map(pu=pu,
                        corePA_col = "core_pa_lulc",
                        oecm_df=oecm_lulc,
                        stud_area,
                        scen = "LULC_balanced",
                        save_output = T)

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

stats_all <- rbind(base_lulc$stats, base_glob$stats) 
stats_all$AREA_SQKM<-stats_all$area/10^6
  # dplyr::mutate(pca = dplyr::recode(pca, !!!pca_labels)) |>
  # dplyr::group_by(lulc_name, scenario, pca) |>
  # dplyr::summarise(km2 = sum(area) / 1e6, .groups = "drop")  
# 
#net <- c("Optimized IUCN Ia or II", "Potential OECM", "No protection", "Other PA (IUCN III-VI)")   # the PA/OECM network
fill_values = c(  "existing core pa" = "#00A300",
               "new core PA" = "#bf8bff",
               "proposed upgrade existing PA" = "#e5d0ff",
               "Other PA (IUCN III-VI) high suitability for OECM" = "#44a6c6",
               "High OECM suitability not protected - potential OECM" = "#194553",
               "Other PA (IUCN III-VI)" = "#ADD8E6",
               "not protected" = "#999999"
               )


# A) core PA + OECM per scenario; pie area = total network area
# plot_area_pies(stats_all, slice_col = "pca", keep_slices = net, slice_levels = net,
#                fill_values = pca_cols, label = "both")
# 
# # B) the same, one row per LULC class
# plot_area_pies(stats_all, slice_col = "pca", facet_col = "lulc_name",
#                keep_slices = net, slice_levels = net, fill_values = pca_cols)
# 
# # C) all four classes (including other PA and no protection)
# plot_area_pies(stats_all, slice_col = "pca", fill_values = pca_cols,
#                slice_levels = names(pca_cols))
# 
# # D) slices by LULC instead of by class (as before)



# Save
ggplot2::ggsave("outputs/WP4/03_pca_landscape/TRD_pie_pca.png",
                plot_area_pies(stats_all, slice_col = "pca", label = "pct", keep_slices = NULL, 
                               facet_col= "lulc_name", fill_values = fill_values),
                width = 11, height = 6, dpi = 300)
