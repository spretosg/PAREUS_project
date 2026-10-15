library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
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

out <- "outputs/WP4/01_PA_analysis_new"

# Original settings
x   <- build_pa_grid("SK021", out_dir = out)
res <- analyse_pa_gap(x, global_target = 0.1,
                      lulc_targets = c(forest = 0.1, water = 0.1, wetland = 0.1, agricultural = 0),
                      out_dir = out)
res$gap                      # table (no longer overwritten by the plot)
res$plots$gap                # plot
x$grid                       # planning units, same columns as before + area_km2, n_pa_poly

# Another site, different core PA definition, majority LULC, coarser cells
x2 <- build_pa_grid("TRD", lulc_method = "majority", out_dir = out)
res <- analyse_pa_gap(x2, global_target = 0.1,
                      lulc_targets = c(forest = 0.1, water = 0.1, wetland = 0.1, agricultural = 0),
                      out_dir = out)

# Target sensitivity: the grid is built once, only the fast part is repeated
#gaps <- lapply(c(0.1, 0.2, 0.3), function(t) analyse_pa_gap(x, global_target = t)$gap)

#### 2. prep optimization
source("code/WP4/02_prep_optim_funtions.R")

# Chained to step 1

#op <- prepare_optim_grid("SK021", grid = x$grid)

# Or from the written file, as in the original script
#op <- prepare_optim_grid("SK021")

op$na_summary                 # NAs per feature before cleaning
length(op$dropped_ids)        # cells removed

# Fill gaps from neighbours instead of dropping, skip layers that don't exist for a site
op2 <- prepare_optim_grid("SK021", na_action = "neighbour_mean", skip_missing = TRUE)

# NA-tolerant ES means
#op3 <- prepare_optim_grid("SK021", es_na_rm = TRUE)

### 3. optim PA
source("code/WP4/03_pa_optim_functions.R")
### NOR case LULC vs global
target_site<-"TRD"
feat_cols <- c("sampled_eco_cond_scaled", "sampled_es_reg_scaled", "sampled_cost_policy_scaled",
               "sampled_es_prov_scaled", "sampled_es_cult_scaled",
               "min_distance_scaled", "sampled_ghm_scaled")

TRD_feature_sets <- list(
  nature = data.frame(column = feat_cols, weight = 1,
                      direction = c("max", "max", "max", "min", "min", "min", "min")),
  multifunctional = data.frame(column = feat_cols, weight = 1,
                               direction = c("max", "max", "max", "max", "max", "min", "min"))
)


res4 <- run_pa_optimization(
    target_site,
    scenario         = "nature",
    feature_sets     = TRD_feature_sets,
    pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
    gap_dir          = "outputs/WP4/01_PA_analysis",
    pu_dir           = "outputs/WP4/02_optim",
    shared_dir       = "data/shared",
    out_dir          = "outputs/WP4/02_optim",
    budget_col       = "target_area",           # or "gap_area"
    dist_transform   = "log",
    min_area_km2     = 0.075, keep_core = TRUE,
    area_weighted    = T,
    lulcs            = c("water","forest","wetland"), run_global = TRUE,
    global_lockout   = c("built-up", "agricultural"),
    boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
    mip_gap = 0.01, time_limit = 600,
    write            = TRUE
)

### SK021 case agri - no agri (close to existing)

##FRA
