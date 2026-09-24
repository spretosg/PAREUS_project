############################################################################
# Pareus WP4 - Prepare data for PA optimization
# Date: 06.08.2026
# Author: Reto Spielhofer (Norwegian Institute for Nature Research)
############################################################################

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("code/WP4/wp4_functions_utils.R")
target_site<-"FRA"


####---- Input and processing ----####
## read PU from 01_pa_status.R
grid<-st_read(paste0("outputs/WP4/01_PA_analysis/",target_site,"_input_grid.json"))

### Features
## Mean ecosystem service (cookbook nr.2)
es_path<-paste0("data/WP2/",target_site,"/es_all/")
# ES groups
files <- list.files(
  es_path,
  pattern = "\\.tif$",
  full.names = TRUE
)

r <- rast(files)

# Check layer names
names(r) <- sub(
  paste0("^", target_site, "_(.*)_pred$"),
  "\\1",
  names(r)
)

es_groups <- list(
  provisioning = c("wild_plant", "wild_hunt", "farm", "mat"),
  regulating   = c("erosion", "flood", "habitat"),
  cultural     = c("aest", "recr","sense")
)

#prov es mean
mean_prov <- min_max_normalize(mean(r[[es_groups$provisioning]]))
mean_reg  <- min_max_normalize(mean(r[[es_groups$regulating]]))
mean_cult <- min_max_normalize(mean(r[[es_groups$cultural]]))
mean_prov
mean_cult
mean_reg

# writeRaster(mean_prov,paste0("data/WP4/prov_",target_site,".tif"))
# writeRaster(mean_reg,paste0("data/WP4/reg_",target_site,".tif"))
# writeRaster(mean_cult,paste0("data/WP4/cult_",target_site,".tif"))

## Ecosystem condition (cookbook nr.2)
es_cond<-terra::rast(paste0("data/WP2/",target_site,"/int.tif"))

### costs
## based on policy coherence (cookbook nr. 3)
cost_policy_file <- file.path(
  "data/WP3/",
  paste0("pol_resistance_",target_site, ".tif")
)
cost_policy<-rast(cost_policy_file)

target_crs<-st_crs(grid)$wkt

rasters <- list(
  cost_policy   = cost_policy,
  es_cult       = mean_cult,
  es_reg        = mean_reg,
  es_prov       = mean_prov,
  eco_cond = es_cond
)

rasters <- Filter(Negate(is.null), rasters)

list2env(rasters, envir = environment())

grid[paste0("sampled_", names(rasters))] <-
  lapply(rasters, \(r)
         terra::extract(r, grid, fun = mean, na.rm = TRUE)[, 2]
  )

#calculate distance from cells outside core PA to core PA as defined in step 01_pa_status.R
out   <- grid%>%filter(existing_corePA == F)
pa_core <- grid%>%filter(existing_corePA == T)

out$min_distance <- apply(st_distance(out, pa_core), 1, min)
pa_core$min_distance<-0
grid<-rbind(out,pa_core)

## remove na grid cells based on features and costs (not)
cols_to_check <- c(
  "sampled_es_cult",
  "sampled_es_prov",
  "sampled_es_reg",
  "sampled_eco_cond",
  "sampled_cost_policy"
)

if ("sampled_cost_policy" %in% names(grid)) {
  cols_to_check <- c(cols_to_check, "sampled_cost_policy")
}

grid_clean <- grid %>%
  filter(!if_any(all_of(cols_to_check), ~ is.na(.) | is.nan(.)))

## scale
if ("sampled_cost_policy" %in% names(grid_clean)) {
  cols_to_check <- c(cols_to_check, "sampled_cost_policy")
}

grid_clean <- zero_one_scale(
  grid_clean,
  cols = cols_to_check
)


st_write(grid_clean, paste0("outputs/WP4/02_optim/",target_site,"_input_final_grid.json"), driver = "GeoJSON", overwrite = T)
